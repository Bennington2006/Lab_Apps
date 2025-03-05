############################################################################################################
# Function to extract vertices from OBJ file
extract_vertices <- function(file) {
  lines <- readLines(file)
  vertices <- do.call(rbind, lapply(lines[startsWith(lines, "v ")], function(line) {
    as.numeric(unlist(strsplit(line, " "))[2:4])
  }))
  vertices <- round(vertices, 4)
  return(vertices)
}

############################################################################################################

# Function to extract RNA spots with 4 decimal places from a .txt file
extract_rna_spots <- function(file) {
  rna_data <- read.table(file, header = TRUE, sep = "\t")
  rna_spots <- rna_data[, c("Position.X", "Position.Y", "Position.Z")]
  colnames(rna_spots) <- c("x", "y", "z")
  rna_spots <- round(rna_spots, 4)
  return(rna_spots)
}

############################################################################################################

# Function to generate random spots within the cell
generate_random_rna_spots <- function(reference_spots, vertices) {
  # Calculate x, y, z ranges directly from vertices
  x_range <- range(vertices[, 1])
  y_range <- range(vertices[, 2])
  z_range <- range(vertices[, 3])
  n <- nrow(reference_spots)  # Match the number of reference spots
  
  # Generate random points within the bounding box
  random_rna_spots <- data.frame(
    x = runif(n, min = x_range[1], max = x_range[2]),
    y = runif(n, min = y_range[1], max = y_range[2]),
    z = runif(n, min = z_range[1], max = z_range[2])
  )
  
  # Filter points to include only those inside the convex hull
  random_rna_spots <- random_rna_spots[inhulln(convhulln(vertices), as.matrix(random_rna_spots)), ]
  
  # If not enough points, regenerate to match the required number
  while (nrow(random_rna_spots) < n) {
    additional_points <- data.frame(
      x = runif(n - nrow(random_rna_spots), min = x_range[1], max = x_range[2]),
      y = runif(n - nrow(random_rna_spots), min = y_range[1], max = y_range[2]),
      z = runif(n - nrow(random_rna_spots), min = z_range[1], max = z_range[2])
    )
    additional_points <- additional_points[inhulln(convhulln(vertices), as.matrix(additional_points)), ]
    random_rna_spots <- rbind(random_rna_spots, additional_points)
  }
  
  # Return exactly n points
  return(random_rna_spots[1:n, ])
}

############################################################################################################

# Generate randomized RNA spots with 5 fixed bins
generate_bin_random_rna_spots <- function(rna_spots, vertices, num_bins = 5) {
  # Calculate x, y, z ranges directly from vertices
  x_range <- range(vertices[, 1])
  y_range <- range(vertices[, 2])
  z_range <- range(vertices[, 3])
  
  # Define bins and bin edges
  x_bins <- seq(x_range[1], x_range[2], length.out = num_bins + 1)
  y_bins <- seq(y_range[1], y_range[2], length.out = num_bins + 1)
  z_bins <- seq(z_range[1], z_range[2], length.out = num_bins + 1)
  
  # Assign RNA spots to bins
  rna_spots$bin_x <- findInterval(rna_spots$x, x_bins)
  rna_spots$bin_y <- findInterval(rna_spots$y, y_bins)
  rna_spots$bin_z <- findInterval(rna_spots$z, z_bins)
  
  # Count spots in bins and randomize within each bin
  random_rna_spots <- do.call(rbind, lapply(split(rna_spots, interaction(rna_spots$bin_x, rna_spots$bin_y, rna_spots$bin_z)), function(bin) {
    if (nrow(bin) == 0) return(NULL)
    data.frame(
      x = runif(nrow(bin), x_bins[bin$bin_x[1]], x_bins[bin$bin_x[1] + 1]),
      y = runif(nrow(bin), y_bins[bin$bin_y[1]], y_bins[bin$bin_y[1] + 1]),
      z = runif(nrow(bin), z_bins[bin$bin_z[1]], z_bins[bin$bin_z[1] + 1])
    )
  }))
  
  # Filter random spots inside the convex hull
  random_rna_spots <- random_rna_spots[inhulln(convhulln(vertices), as.matrix(random_rna_spots)), ]
  
  return(random_rna_spots)
}

####################################################################################################################

# Density bins
generate_randomized_spots <- function(vertices, rna_spots) {
  library(geometry)
  library(FNN)
  
  # Set k dynamically based on RNA spots
  k <- max(1, round(nrow(rna_spots) / 5))
  
  # Step 1: Define grid and calculate density
  x_bins <- seq(min(vertices[, 1]), max(vertices[, 1]), length.out = 101)
  y_bins <- seq(min(vertices[, 2]), max(vertices[, 2]), length.out = 101)
  z_bins <- seq(min(vertices[, 3]), max(vertices[, 3]), length.out = 11)
  grid_points <- expand.grid(x = x_bins, y = y_bins, z = z_bins)
  
  hull <- convhulln(vertices)
  inside <- inhulln(hull, as.matrix(grid_points))
  grid_points_inside <- grid_points[inside, ]
  
  # Estimate density using k-NN
  grid_density <- knnx.dist(data = as.matrix(rna_spots), query = as.matrix(grid_points_inside), k = k)
  grid_points_inside$density <- 1 / (rowMeans(grid_density) * 10)
  grid_points_inside$normalized_density <- grid_points_inside$density / max(grid_points_inside$density)
  
  # Step 2: Assign density regions based on quantiles 
  density_thresholds <- quantile(grid_points_inside$normalized_density, probs = seq(0.2, 0.8, by = 0.2))
  
  grid_points_inside$region <- cut(
    grid_points_inside$normalized_density,
    breaks = c(-Inf, density_thresholds, Inf),
    labels = paste("Region", 5:1),  # High to Low Density
    include.lowest = TRUE
  )
  
  # Map vertices to regions
  vertex_density <- knnx.index(data = as.matrix(grid_points_inside[, c("x", "y", "z")]), query = as.matrix(vertices), k = 1)
  vertices_density_region <- grid_points_inside$region[vertex_density]
  
  # Step 3: Randomize spots by region (keeping empty regions empty)
  randomize_by_region <- function(rna_spots, outer_hull, inner_hull = NULL, region_name) {
    inside_outer <- inhulln(outer_hull, as.matrix(rna_spots[, c("x", "y", "z")]))
    inside_inner <- if (!is.null(inner_hull)) inhulln(inner_hull, as.matrix(rna_spots[, c("x", "y", "z")])) else rep(FALSE, nrow(rna_spots))
    
    # Keep original region empty if no spots exist
    filtered_spots <- rna_spots[inside_outer & !inside_inner, , drop = FALSE]
    num_spots <- nrow(filtered_spots)
    
    if (num_spots == 0) {
      return(data.frame(x = numeric(0), y = numeric(0), z = numeric(0), region = character(0)))  # Fix
    }
    
    random_points <- data.frame(x = numeric(), y = numeric(), z = numeric())
    
    while (nrow(random_points) < num_spots) {
      candidates <- data.frame(
        x = runif(num_spots, min(vertices[, 1]), max(vertices[, 1])),
        y = runif(num_spots, min(vertices[, 2]), max(vertices[, 2])),
        z = runif(num_spots, min(vertices[, 3]), max(vertices[, 3]))
      )
      inside <- inhulln(outer_hull, as.matrix(candidates)) & !(if (!is.null(inner_hull)) inhulln(inner_hull, as.matrix(candidates)) else FALSE)
      random_points <- rbind(random_points, candidates[inside, ])
    }
    
    random_points <- random_points[1:num_spots, ]
    random_points$region <- region_name
    return(random_points)
  }
  
  # Generate random spots for each density region
  region_labels <- levels(grid_points_inside$region)
  random_spots_list <- list()
  
  for (i in seq_along(region_labels)) {
    outer_hull <- convhulln(vertices[vertices_density_region == region_labels[i], ])
    inner_hull <- if (i < length(region_labels)) convhulln(vertices[vertices_density_region == region_labels[i + 1], ]) else NULL
    random_spots_list[[i]] <- randomize_by_region(rna_spots, outer_hull, inner_hull, region_labels[i])
  }
  
  # Combine all randomized spots
  randomized_spots <- do.call(rbind, random_spots_list)
  colnames(randomized_spots) <- c("x", "y", "z", "region")
  
  return(randomized_spots)
}

############################################################################################################

# Function to calculate and print the percentage of spots within a certain distance
colocper <- function(spots1, spots2, threshold = 1) {
  # Calculate pairwise Euclidean distances
  distances <- as.matrix(dist(rbind(spots1, spots2)))
  distances <- distances[1:nrow(spots1), (nrow(spots1) + 1):nrow(distances)]
  
  # Count the number of spots in spots1 within the threshold distance of spots2
  close_spots1 <- apply(distances, 1, function(row) any(row <= threshold))
  percentage_spots1 <- mean(close_spots1) * 100
  
  # Count the number of spots in spots2 within the threshold distance of spots1
  close_spots2 <- apply(distances, 2, function(col) any(col <= threshold))
  percentage_spots2 <- mean(close_spots2) * 100
  # Create a result table
  result_table <- data.frame(
    Dataset1_to_Dataset2 = percentage_spots1,
    Dataset2_to_Dataset1 = percentage_spots2
  )
  
  # Return the result table
  return(result_table)
}

############################################################################################################
# Function to plot RNA spots and cell boundary
plot_rna_spots <- function(vertices, spots_1a, spots_1b, title, legend_labels = c("1a Spots", "1b Spots", "Cell Boundary")) {
  # Open a new 3D plot
  open3d()
  # Plot the cell boundary
  plot3d(vertices, col = "lightblue", size = 1, alpha = 0.3, xlab = "X", ylab = "Y", zlab = "Z")
  aspect3d(diff(range(vertices[, 1])), diff(range(vertices[, 2])), diff(range(vertices[, 3])))
  # Plot RNA spots for 1a and 1b
  points3d(spots_1a$x, spots_1a$y, spots_1a$z, col = "green", size = 3)
  points3d(spots_1b$x, spots_1b$y, spots_1b$z, col = "red", size = 3)
  # Add a title
  title3d(title, xlab = "X", ylab = "Y", zlab = "Z")
  # Add a legend
  legend3d("topright", legend = legend_labels, col = c("green", "red", "lightblue"), pch = 16, cex = 1.2)
}