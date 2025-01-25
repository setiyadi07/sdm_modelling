# Load Required Packages
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
if (!requireNamespace("rpart", quietly = TRUE)) install.packages("rpart")
if (!requireNamespace("sp", quietly = TRUE)) install.packages("sp")
if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("rspat", quietly = TRUE)) remotes::install_github("rspatial/rspat")
if (!requireNamespace("rgbif", quietly = TRUE)) install.packages("rgbif")
if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata")
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("predicts", quietly = TRUE)) install.packages("predicts")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("raster", quietly = TRUE)) install.packages("raster")

# Load libraries
library(caret)
library(randomForest)
library(rpart)
library(sp)
library(terra)
library(rspat)
library(rgbif)
library(geodata)
library(sf)
library(predicts)
library(pROC)
library(raster)

# Define folder location - set based on your desire location
folder_location <- "D:/Office/SDM_Biodiv/Data/"

# Step 1: Fetch Occurrence Data for Pongo pygmaeus (Orangutan)
occurrences <- occ_search(scientificName = "Pongo pygmaeus", limit = 10000)
occ_data <- occurrences$data

# Filter for Valid Coordinates
coords <- occ_data[!is.na(occ_data$decimalLongitude) & !is.na(occ_data$decimalLatitude),
                   c("decimalLongitude", "decimalLatitude")]

# Step 2: Define Borneo Bounding Box
borneo_bounds <- list(lon_min = 108, lon_max = 119, lat_min = 0, lat_max = 7)

# Step 3: Filter Points within Borneo Bounds
coords_filtered <- coords[
  coords$decimalLongitude >= borneo_bounds$lon_min &
    coords$decimalLongitude <= borneo_bounds$lon_max &
    coords$decimalLatitude >= borneo_bounds$lat_min &
    coords$decimalLatitude <= borneo_bounds$lat_max, ]

# Remove Duplicate Points
coords_filtered <- coords_filtered[!duplicated(coords_filtered), ]

# Step 4: Load and Process Map Data
world <- geodata::world(path = ".")
borneo_countries <- world[world$NAME_0 %in% c("Indonesia", "Malaysia", "Brunei Darussalam"), ]
borneo_sf <- st_as_sf(borneo_countries)

# Ensure CRS Consistency
coords_sf <- st_as_sf(coords_filtered, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
coords_sf <- st_transform(coords_sf, st_crs(borneo_sf))

# Step 5: Visualize Occurrences on Borneo Map (with zoom on Borneo)
plot(st_geometry(borneo_sf), col = "white", main = "Pongo pygmaeus in Borneo",xlab = "Longitude", ylab = "Latitude",
     xlim = c(108, 119), ylim = c(0, 7), axes = FALSE)

# Add Points
points(coords_filtered$decimalLongitude, coords_filtered$decimalLatitude, 
       col = "red", pch = 20, cex = 1)

# Add Box and Axes
box(lwd = 2)
axis(1, at = seq(108, 119, by = 1), labels = seq(108, 119, by = 1)) # Longitude axis
axis(2, at = seq(0, 7, by = 1), labels = seq(0, 7, by = 1), las = 1) # Latitude axis


# Step 6: Reduce Sampling Bias with SpatRaster Grid
acv <- vect(coords_sf)
r <- rast(acv)
res(r) <- 1
r <- extend(r, ext(r) + 1)
set.seed(13)
acsel <- spatSample(acv, size = 1, method = "random", strata = r)

# Step 7: Visualize Sampling Bias Reduction
p <- as.polygons(r)
plot(st_geometry(borneo_sf), col = "white", main = "Sampling Bias Reduction for Pongo pygmaeus",
     xlab = "Longitude", ylab = "Latitude", axes = TRUE, xlim = c(108, 119), ylim = c(0, 7))
plot(p, add = TRUE, border = "gray")
points(acv, col = "blue", pch = 19, cex = 0.5)
points(acsel, col = "red", pch = "x", cex = 1)

# Step 8: Save Cleaned Dataset
acsel_df <- as.data.frame(acsel)
write.csv(acsel_df, "cleaned_occurrence_data.csv", row.names = FALSE)
saveRDS(acsel, "cleaned_occurrence_data.rds")

# Step 9: Fetch and Process WorldClim Data
unlink(list.files(".", pattern = "worldclim", full.names = TRUE))
wc <- geodata::worldclim_global("bio", res = 2.5, path = ".")
borneo_extent <- ext(c(109, 119, 0, 7))
wc_borneo <- crop(wc, borneo_extent)

# Step 10: Visualize Selected Bioclimatic Variables
selected_indices <- c(1, 5, 6, 7, 8, 9, 12, 16, 17)
wc_selected <- wc_borneo[[selected_indices]]
plot(wc_selected)

# Step 11: Generate Random Samples from Bioclimatic Data
set.seed(0)
bg <- spatSample(wc_selected, 800, "random", na.rm = TRUE, xy = TRUE)

# Check the first few rows of the sampled background data
head(bg)

# Step 12: Plot Sampled Points on Bioclimatic Data
plot(wc_selected[[1]], main = "Bioclimatic Variable with Background Sampled Points")
points(bg[, c("x", "y")], col = "black", pch = 20, cex = 0.5)

# Step 13: Analyze Relationships between Variables
# Ensure the required variables exist
required_vars <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_12")
if (all(required_vars %in% names(bg))) {
  plot(main = "Relationships Annual (Mean Temperature - Precipitation) Variables",
       bg[, required_vars[1]],  # Annual Mean Temperature
       bg[, required_vars[2]],  # Annual Precipitation
       xlab = "Annual Mean Temperature (°C)",
       ylab = "Annual Precipitation (mm)",
       cex = 0.8, col = "blue", pch = 20)
} else {
  cat("Required variables not found in bg.\n")
}

# Step 14: Create a 50 km Buffer around Each Sampled Point
# Convert the random sample points to a SpatialPoints object
ac <- SpatialPoints(cbind(bg$x, bg$y))

# Define the CRS as WGS 84 (latitude/longitude)
proj4string(ac) <- CRS("+proj=longlat +datum=WGS84")

# Create a 50 km buffer around each sampled point
buffer_area <- buffer(ac, 5000)  # Buffer with 50 km radius

# Step 15: Randomly Sample Points within Each Buffer
set.seed(999)
# Sample random points from the buffer areas (ensure one sample per buffer)
samp1 <- spsample(buffer_area, n = 800, type = "random")

# Step 16: Plot the Results
# Plot the 50 km Buffer with Background wc_selected[[1]] and Title
plot(wc_selected[[1]], main = "Buffer Point for Background Sample Extracted", axes = TRUE)  # Plot the bioclimatic background
plot(buffer_area, add = TRUE, border = "darkred")  # Add the buffer areas on top of the background
points(samp1, pch = "+", cex = 0.5)  # Plot random sampled points within buffers


# Step 18: Extract environmental data at presence and background points
# Convert buffer_area to a 'terra' object
buffer_area_terra <- vect(buffer_area)
# Extract environmental data at the presence points
presvals <- extract(wc_selected, coords_filtered)
head(presvals)

# Remove the ID column and NA values in the bioclimatic value columns
presvals_clean <- na.omit(presvals[, -1])

# Desired sample size
desired_sample_size <- 800

# Adjust `presvals_clean` to the desired size
set.seed(123)  # Set seed for reproducibility

if (nrow(presvals_clean) > desired_sample_size) {
  # If `presvals_clean` has more than 800 rows, randomly sample 800
  presvals_clean <- presvals_clean[sample(nrow(presvals_clean), desired_sample_size), ]
} else if (nrow(presvals_clean) < desired_sample_size) {
  # If `presvals_clean` has fewer than 800 rows, replicate rows to reach 800
  presvals_clean <- presvals_clean[rep(1:nrow(presvals_clean), length.out = desired_sample_size), ]
}

# Confirm the number of rows
nrow(presvals_clean)


# Extract environmental data at the background points
# Extract only the x and y coordinates from bg
bg_coords <- bg[, 1:2]

# Extract values from the raster using the coordinates
bgvals <- extract(wc_selected, bg_coords)

# Remove the ID column and NA values in the bioclimatic value columns
bgvals_clean <- na.omit(bgvals[, -1])

# Check the extracted values for background data
head(bgvals_clean)

str(bgvals_clean)  # Check structure
colnames(bgvals_clean)  # Check column names


# Create a plot comparing background and presence points

# Plot background data (using Annual Mean Temperature vs. Annual Precipitation)
# Corrected plot using the actual columns for Annual Mean Temperature and Annual Precipitation
plot(bgvals_clean$wc2.1_2.5m_bio_1,  # Annual Mean Temperature
     bgvals_clean$wc2.1_2.5m_bio_12, # Annual Precipitation
     xlab = "Annual Mean Temperature (°C)",
     ylab = "Annual Precipitation (mm)", 
     cex = 0.8, 
     main = "Background vs Presence Data")

# Add presence points as red circles
points(presvals_clean$wc2.1_2.5m_bio_1,  # Annual Mean Temperature
       presvals_clean$wc2.1_2.5m_bio_12, # Annual Precipitation
       col = "red", 
       cex = 0.7, 
       pch = 1)

# Add a legend
legend("topleft", 
       legend = c("Presence", "Background"), 
       col = c("red", "black"), 
       pch = c(1, 1), 
       pt.cex = c(0.7, 0.8))


# Step 19: Combine presence and background data with a pb column

# Create the presence/background indicator column (pb)
pb <- c(rep(1, nrow(presvals_clean)),  # 1 for presence points
        rep(0, nrow(bgvals_clean)))   # 0 for background points

# Combine presence and background data with the pb column
sdmdata <- data.frame(pb = pb, 
                      rbind(presvals_clean, bgvals_clean))

# Inspect the combined dataset
head(sdmdata)   # View the first few rows
tail(sdmdata)   # View the last few rows
summary(sdmdata)  # Summary statistics for all variables

# Visualize relationships between selected environmental variables
pairs(sdmdata[, 2:5], 
      cex = 0.1, 
      main = "Scatterplot Matrix of Environmental Variables")

# Save the combined data and presence values for future use
saveRDS(sdmdata, "sdm.Rds")   # Save the combined data
saveRDS(presvals, "pvals.Rds") # Save the raw presence values

# Step 20: Train a Random Forest Model

# Set seed for reproducibility
set.seed(123)

# Split the data into training and testing sets using createDataPartition (caret)
trainIndex <- createDataPartition(sdmdata$pb, p = .7, list = FALSE)
trainData <- sdmdata[trainIndex, ]
testData <- sdmdata[-trainIndex, ]

# Convert response variable to factor for classification
trainData$pb <- as.factor(trainData$pb)
testData$pb <- as.factor(testData$pb)

# Train Random Forest Model with 500 trees and importance calculation
rf_model <- randomForest(pb ~ ., data = trainData, ntree = 500, importance = TRUE)

# Display Random Forest Model Summary
cat("Random Forest Model Summary:\n")
print(rf_model)

# Evaluate Model Performance on the Entire Dataset (Accuracy, RMSE)
# 1. Accuracy on training data
predictions_rf <- predict(rf_model, sdmdata, type = "response")
accuracy_rf <- sum(predictions_rf == sdmdata$pb) / length(sdmdata$pb)
cat("Random Forest Accuracy: ", round(accuracy_rf, 4), "\n")
# Confusion Matrix with Kappa
confusion_rf <- confusionMatrix(predictions_rf_test, testData$pb)
cat("Random Forest Kappa: ", round(confusion_rf$overall['Kappa'], 4), "\n")

# Print the confusion matrix
cat("Confusion Matrix:\n")
print(confusion_rf$table)


# 2. RMSE on predicted probabilities
pred_probs_rf <- predict(rf_model, sdmdata, type = "prob")[, 2]  # Probability for the '1' class
rmse_rf <- sqrt(mean((pred_probs_rf - sdmdata$pb)^2))
cat("Random Forest RMSE: ", round(rmse_rf, 4), "\n")

# Evaluate Test Data
# 1. Confusion Matrix (Test Data)
predictions_rf_test <- predict(rf_model, testData, type = "response")
confusion_rf <- confusionMatrix(predictions_rf_test, testData$pb)
cat("Random Forest Confusion Matrix (Test Data):\n")
print(confusion_matrix_rf)

# 2. Accuracy on Test Data
accuracy_rf_test <- sum(predictions_rf_test == testData$pb) / length(testData$pb)
cat("Random Forest Test Accuracy: ", round(accuracy_rf_test, 4), "\n")

# 3. ROC Curve and AUC for Test Data
pred_probs_rf_test <- predict(rf_model, testData, type = "prob")[, 2]  # Probability for the '1' class
roc_curve_rf <- roc(testData$pb, pred_probs_rf_test)
# Plot the ROC curve
plot(roc_curve_rf, main = "ROC Curve for Random Forest Model")

# Add the AUC value at the center-left of the plot
auc_value <- auc(roc_curve_rf)
text(x = 0.1, y = 0.5, labels = paste("AUC = ", round(auc_value, 2)), col = "red", cex = 1.5, pos = 4)

# Variable Importance
importance_values <- importance(rf_model)
importance_table <- data.frame(
  Variable = rownames(importance_values),
  X.incMSE = as.numeric(importance_values[, 1]),  # X-incMSE for Mean Squared Error
  IncNodePurity = as.numeric(importance_values[, 2])  # IncNodePurity for Node Purity
)
importance_table <- importance_table[order(-importance_table$X.incMSE), ]
top_variables <- importance_table$Variable[1:5]  # Top 5 variables by X.incMSE

cat("Top Variables by RF Variable Importance:\n")
print(top_variables)

# Plot Variable Importance
cat("Plotting Variable Importance:\n")
varImpPlot(rf_model, main = "Variable Importance in Random Forest Model")

# Predict Habitat Suitability for wc_borneo (or wc_selected, based on your dataset)
predicted_suitability <- predict(wc_borneo, rf_model, type = "response")  # Change 'wc_borneo' to 'wc_selected' if needed
predicted_suitability_prob <- predict(wc_borneo, rf_model, type = "prob")  # Change 'wc_borneo' to 'wc_selected' if needed

par(mfrow = c(1, 2))  # Adjust the layout to 2x2 grid
par(mar = c(5, 5, 2, 2))  # Adjust margins

# Visualize Predicted Habitat Suitability
plot(predicted_suitability, main = "Predicted Habitat Suitability for Pongo pygmaeus (Current)")
plot(predicted_suitability_prob, main = "Probable Habitat Suitability for Pongo pygmaeus (Current)",2)
par(mfrow = c(1, 1))


# Write raster files to the defined location
writeRaster(predicted_suitability, paste0(folder_location, "predicted_habitat_suitability_orangutan.tif"), overwrite = TRUE)
writeRaster(predicted_suitability_prob, paste0(folder_location, "probable_habitat_suitability_orangutan.tif"), overwrite = TRUE)


# Save the Random Forest Model for later use
saveRDS(rf_model, "rf_model_updated.Rds")


### Step 22: Future Prediction 2061-2080 (SSP245 & SSP585)

# Download predicted climate data
forecast_data_ssp245_2041_2060 <- cmip6_world(
  model = "MPI-ESM1-2-HR",
  ssp = "245",
  time = "2041-2060",
  var = "bioc",
  res = 2.5,
  path = folder_location
)

forecast_data_ssp245_2061_2080 <- cmip6_world(
  model = "MPI-ESM1-2-HR",
  ssp = "245",
  time = "2061-2080",
  var = "bioc",
  res = 2.5,
  path = folder_location
)

# Check if both SSP245 files exist
ssp245_file_path_2041 <- file.path(folder_location, "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060.tif")
ssp245_file_path_2061 <- file.path(folder_location, "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2061-2080.tif")

if (!file.exists(ssp245_file_path_2041) || !file.exists(ssp245_file_path_2061)) {
  stop("One or both SSP245 climate data files are missing. Please check the file paths.")
}

# Load SSP245 datasets
bioclim_data_ssp245_2041 <- brick(ssp245_file_path_2041)  # Load SSP245 (2041-2060) as RasterBrick
bioclim_data_ssp245_2061 <- brick(ssp245_file_path_2061)  # Load SSP245 (2061-2080) as RasterBrick

# Select bands to extract
selected_bands <- c(1, 5, 6, 7, 8, 9, 12, 16, 17)

# Extract and save selected bands for SSP245 (2041-2060)
selected_layers_ssp245_2041 <- bioclim_data_ssp245_2041[[selected_bands]]
ssp245_2041_output_file <- file.path(folder_location, "ssp245_2041_selected_bands.tif")
writeRaster(selected_layers_ssp245_2041, ssp245_2041_output_file, format = "GTiff", overwrite = TRUE)

# Extract and save selected bands for SSP245 (2061-2080)
selected_layers_ssp245_2061 <- bioclim_data_ssp245_2061[[selected_bands]]
ssp245_2061_output_file <- file.path(folder_location, "ssp245_2061_selected_bands.tif")
writeRaster(selected_layers_ssp245_2061, ssp245_2061_output_file, format = "GTiff", overwrite = TRUE)

cat("SSP245 (2041-2060) output saved at:", ssp245_2041_output_file, "\n")
cat("SSP245 (2061-2080) output saved at:", ssp245_2061_output_file, "\n")

# Rename layers to match trainData
names(selected_layers_ssp245_2041) <- names(trainData)[2:10]
names(selected_layers_ssp245_2061) <- names(trainData)[2:10]

# Crop data to desired extent for both datasets
selected_layers_ssp245_2041 <- terra::rast(selected_layers_ssp245_2041)  # Convert to SpatRaster
forecast_data_ssp245_2041 <- terra::crop(selected_layers_ssp245_2041, borneo_extent)

selected_layers_ssp245_2061 <- terra::rast(selected_layers_ssp245_2061)  # Convert to SpatRaster
forecast_data_ssp245_2061 <- terra::crop(selected_layers_ssp245_2061, borneo_extent)

# Summarize cropped data
cat("Summary of SSP245 forecast data (2041-2060):\n")
summary(forecast_data_ssp245_2041)

cat("Summary of SSP245 forecast data (2061-2080):\n")
summary(forecast_data_ssp245_2061)

# Predict presence using trained model
forecast_presence_ssp245_2041 <- predict(forecast_data_ssp245_2041, rf_model, type = "response")
forecast_presence_ssp245_2061 <- predict(forecast_data_ssp245_2061, rf_model, type = "response")

forecast_prob_ssp245_2041 <- predict(forecast_data_ssp245_2041, rf_model, type = "prob")
forecast_prob_ssp245_2061 <- predict(forecast_data_ssp245_2061, rf_model, type = "prob")

# Plot predictions
par(mfrow = c(2, 2))  # Adjust layout for 2x2 grid
par(mar = c(5, 5, 2, 2))  # Adjust margins

# Plot SSP245 (2041-2060) results
plot(forecast_presence_ssp245_2041, 
     main = "Predicted Habitat Suitability 2041-2060 - SSP245")
plot(forecast_prob_ssp245_2041, 
     main = "Probable Habitat Suitability 2041-2060 - SSP245", 2)

# Plot SSP245 (2061-2080) results
plot(forecast_presence_ssp245_2061, 
     main = "Predicted Habitat Suitability 2061-2080 - SSP245")
plot(forecast_prob_ssp245_2061, 
     main = "Probable Habitat Suitability 2061-2080 - SSP245", 2)
