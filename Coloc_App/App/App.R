library(shiny)
library(rgl)
library(geometry)
library(parallel)
library(FNN)

# Load custom functions
source("Core_Function.R")

# Define UI
ui <- fluidPage(
  titlePanel("Colocalization Runner"),
  sidebarLayout(
    sidebarPanel(
      fileInput("obj_file", "Upload OBJ File", accept = c(".obj")),
      fileInput("txt_file_1a", "Upload 1a XYZ TXT File", accept = c(".txt")),
      fileInput("txt_file_1b", "Upload 1b XYZ TXT File", accept = c(".txt")),
      numericInput("num_trials", "Number of Trials", value = 200, min = 1),
      numericInput("threshold", "Threshold", value = 1, min = 0),
      actionButton("run", "Run")
    ),
    mainPanel(
      verbatimTextOutput("results")
    )
  )
)

# Define server logic
server <- function(input, output) {
  observeEvent(input$run, {
    req(input$obj_file, input$txt_file_1a, input$txt_file_1b)
    
    num_trials <- input$num_trials
    threshold <- input$threshold
    
    start_time <- Sys.time()
    
    obj_file <- input$obj_file$datapath
    txt_file_1a <- input$txt_file_1a$datapath
    txt_file_1b <- input$txt_file_1b$datapath
    
    vertices <- extract_vertices(obj_file)
    rna_spots_1a <- extract_rna_spots(txt_file_1a)
    rna_spots_1b <- extract_rna_spots(txt_file_1b)
    
    # Debugging statements
    print("Vertices:")
    print(dim(vertices))
    print(head(vertices))
    
    print("RNA Spots 1a:")
    print(dim(rna_spots_1a))
    print(head(rna_spots_1a))
    
    print("RNA Spots 1b:")
    print(dim(rna_spots_1b))
    print(head(rna_spots_1b))
    
    colocog <- colocper(rna_spots_1a, rna_spots_1b, threshold = threshold)
    
    num_cores <- max(1, detectCores() - 1)
    cl <- makeCluster(num_cores, type = "PSOCK")
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterEvalQ(cl, {
      library(geometry)
      library(FNN)
    })
    
    withProgress(message = 'Running analysis', value = 0, {
      run_trial_randomly_random <- function(i) {
        random_rna_spots_1a <- generate_random_rna_spots(rna_spots_1a, vertices)
        random_rna_spots_1b <- generate_random_rna_spots(rna_spots_1b, vertices)
        coloc_values <- colocper(random_rna_spots_1a, random_rna_spots_1b, threshold = threshold)
        c(coloc_1a_with_1b = as.numeric(coloc_values[1]), coloc_1b_with_1a = as.numeric(coloc_values[2]))
      }
      
      clusterExport(cl, varlist = c("run_trial_randomly_random", "generate_random_rna_spots", "generate_randomized_spots", "generate_bin_random_rna_spots", "rna_spots_1a", "rna_spots_1b", "vertices", "threshold", "colocper", "inhulln", "convhulln", "extract_vertices", "extract_rna_spots", "plot_rna_spots"), envir = environment())
      
      results_randomly_random <- parLapply(cl, 1:num_trials, run_trial_randomly_random)
      coloc_results_randomly_random <- as.data.frame(do.call(rbind, results_randomly_random))
      coloc_results_randomly_random$Trial <- 1:num_trials
      randomly_random <- colMeans(coloc_results_randomly_random[, 1:2], na.rm = TRUE)
      incProgress(1 / 3, detail = "Completed Randomly Random Analysis")
      
      run_trial_equal_bins <- function(trial_index) {
        bin_random_rna_spots_1a <- generate_bin_random_rna_spots(rna_spots_1a, vertices)
        bin_random_rna_spots_1b <- generate_bin_random_rna_spots(rna_spots_1b, vertices)
        coloc_values <- colocper(bin_random_rna_spots_1a, bin_random_rna_spots_1b, threshold = threshold)
        return(c(coloc_values$Dataset1_to_Dataset2, coloc_values$Dataset2_to_Dataset1))
      }
      
      results_equal_bins <- parLapply(cl, 1:num_trials, run_trial_equal_bins)
      coloc_results_equal_bins <- as.data.frame(do.call(rbind, results_equal_bins))
      coloc_results_equal_bins$Trial <- 1:num_trials
      equal_bins <- colMeans(coloc_results_equal_bins[, 1:2], na.rm = TRUE)
      incProgress(1 / 3, detail = "Completed Equal Bins Analysis")
      
      run_trial_density_bins <- function(trial_index) {
        knn_randomized_rna_spots_1a <- generate_randomized_spots(vertices, as.matrix(rna_spots_1a))
        knn_randomized_rna_spots_1b <- generate_randomized_spots(vertices, as.matrix(rna_spots_1b))
        knn_randomized_rna_spots_1a <- knn_randomized_rna_spots_1a[, c("x", "y", "z")]
        knn_randomized_rna_spots_1b <- knn_randomized_rna_spots_1b[, c("x", "y", "z")]
        coloc_values <- colocper(knn_randomized_rna_spots_1a, knn_randomized_rna_spots_1b, threshold = threshold)
        return(c(coloc_values$Dataset1_to_Dataset2, coloc_values$Dataset2_to_Dataset1))
      }
      
      results_density_bins <- parLapply(cl, 1:num_trials, run_trial_density_bins)
      coloc_results_density_bins <- as.data.frame(do.call(rbind, results_density_bins))
      coloc_results_density_bins$Trial <- 1:num_trials
      density_bins <- colMeans(coloc_results_density_bins[, 1:2], na.rm = TRUE)
      incProgress(1 / 3, detail = "Completed Density Bins Analysis")
    })
    
    end_time <- Sys.time()
    execution_time <- end_time - start_time
    
    output$results <- renderPrint({
      cat("Total Execution Time:", execution_time, "seconds\n")
      cat("Results:\n")
      print(colocog)
      print(randomly_random)
      print(equal_bins)
      print(density_bins)
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)