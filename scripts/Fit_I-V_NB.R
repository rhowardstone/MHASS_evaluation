# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if output directory is provided
if (length(args) < 1) {
  cat("Usage: Rscript script.R <output_directory>\n")
  cat("Example: Rscript script.R ./results\n")
  stop("Please provide an output directory as an argument.")
}

# Get output directory from arguments
output_dir <- args[1]

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load required libraries
library(metaSPARSim)

# Load R1 data
data(R1)

# Combine all groups into one dataset
all_data <- data.frame()
for (i in 1:length(R1)) {
  group_data <- R1[[i]]
  temp_df <- data.frame(
    group = names(R1)[i],
    intensity = group_data$intensity,
    variability = as.numeric(group_data$variability),
    stringsAsFactors = FALSE
  )
  all_data <- rbind(all_data, temp_df)
}

# Filter data: remove NA variability and zero intensity
cat("Original data points:", nrow(all_data), "\n")
filtered_data <- all_data[!is.na(all_data$variability) & 
                         all_data$intensity > 0 & 
                         all_data$variability > 0, ]
cat("Filtered data points:", nrow(filtered_data), "\n")

I <- filtered_data$intensity
V <- filtered_data$variability

# Fit Negative Binomial CV model: CV² = 1/μ + φ
# Since V is already CV² in metaSPARSim: V = β₁*(1/I) + β₀
# Where β₀ = φ (dispersion) and β₁ = 1 (in theory)
x_nb <- 1/I
fit_nb <- lm(V ~ x_nb)
beta_nb <- coef(fit_nb)

# Extract NB parameters
phi_estimate <- beta_nb[1]  # Intercept = dispersion parameter
inverse_mean_coef <- beta_nb[2]  # Should be ~1 for true NB

# Print model parameters in NB context
cat("\n=== Negative Binomial Model ===\n")
cat("Model: CV² = ", inverse_mean_coef, " * (1/μ) + ", phi_estimate, "\n", sep="")
cat("\nNB Parameter Estimates:\n")
cat("Dispersion (φ):", phi_estimate, "\n")
cat("1/μ coefficient:", inverse_mean_coef, " (expect 1.0 for true NB)\n")

# Check if this follows NB assumptions
if (inverse_mean_coef < 0) {
  cat("\nWARNING: Negative coefficient for 1/μ term (", inverse_mean_coef, 
      ") violates NB assumptions!\n", sep="")
  cat("This suggests the data does NOT follow a negative binomial distribution.\n")
}

# Calculate R-squared
r2 <- summary(fit_nb)$r.squared
cat("\nModel fit: R² =", r2, "\n")

# Convert to variance parameterization
# From CV² = inverse_mean_coef/μ + phi_estimate
# We get: variance = σ² = μ² × CV² = inverse_mean_coef*μ + phi_estimate*μ²
linear_coef <- inverse_mean_coef      # Coefficient of μ
quadratic_coef <- phi_estimate        # Coefficient of μ²

cat("\n=== Variance Model Parameters ===\n")
cat("Variance = ", linear_coef, " * μ + ", quadratic_coef, " * μ²\n", sep="")
cat("Standard deviation = sqrt(variance)\n")

# Prepare prediction grid
I_grid <- 10^seq(log10(min(I)), log10(max(I)), length.out = 300)
V_nb <- pmax(0, phi_estimate + inverse_mean_coef * (1/I_grid))

# Create single log-log plot
png_file <- file.path(output_dir, "nb_model_fit.png")
png(png_file, width = 8, height = 6, units = "in", res = 150)
par(mar = c(4, 4, 3, 1))

# Log-log plot
plot(I, V, log = "xy", pch = 16, cex = 0.3, col = rgb(0.5, 0.5, 0.5, 0.5),
     xlab = "Intensity (μ)", ylab = "Variability (CV²)", 
     main = "Intensity vs. Variability: Negative Binomial Model",
     xlim = range(I), ylim = range(V[V>0]))
lines(I_grid, V_nb, col = "red", lwd = 3)
legend("topleft", 
       legend = c("Data", 
                  paste0("CV² = ", round(inverse_mean_coef, 2), "/μ + ", 
                         round(phi_estimate, 3))),
       col = c("gray", "red"),
       lty = c(NA, 1),
       pch = c(16, NA),
       lwd = c(NA, 3),
       bty = "n")

dev.off()

# Save the model summary
cat("\n=== Full Model Summary ===\n")
print(summary(fit_nb))

# Show some predictions at key intensity values
cat("\n=== Predictions at Key Mean Values ===\n")
cat("μ (mean) | CV² (predicted) | CV | Variance (σ²)\n")
test_means <- c(0.001, 0.01, 0.1, 1, 10, 100)
for (mu in test_means) {
  cv2_pred <- phi_estimate + inverse_mean_coef / mu
  cv_pred <- sqrt(max(0, cv2_pred))
  # Using our fitted model
  var_pred <- linear_coef * mu + quadratic_coef * mu^2
  cat(sprintf("%8.3f | %14.4f | %6.3f | %16.4f\n", 
              mu, cv2_pred, cv_pred, var_pred))
}

# Save the original NB-style parameters
nb_params_file <- file.path(output_dir, "nb_model_parameters.csv")
write.csv(data.frame(
  parameter = c("dispersion_phi", "inverse_mean_coefficient"), 
  value = c(phi_estimate, inverse_mean_coef),
  expected_for_NB = c("positive", "1.0"),
  interpretation = c("NB dispersion parameter", "Coefficient of 1/μ term")
), nb_params_file, row.names = FALSE)

# Save the variance model parameters (more intuitive)
var_params_file <- file.path(output_dir, "variance_model_parameters.csv")
write.csv(data.frame(
  parameter = c("linear_coefficient", "quadratic_coefficient"),
  value = c(linear_coef, quadratic_coef),
  formula_term = c("μ", "μ²"),
  interpretation = c(
    "Coefficient of μ in variance formula",
    "Coefficient of μ² in variance formula"
  ),
  formula = c(
    "variance = linear_coef * μ + quadratic_coef * μ²",
    "CV² = variance / μ²"
  )
), var_params_file, row.names = FALSE)

# Save filtered data
filtered_data_file <- file.path(output_dir, "R1_filtered_data.csv")
write.csv(filtered_data, filtered_data_file, row.names = FALSE)

cat("\nFiles saved in", output_dir, ":\n")
cat("- nb_model_fit.png (plot)\n")
cat("- R1_filtered_data.csv (filtered data)\n")
cat("- nb_model_parameters.csv (original NB-style parameters)\n")
cat("- variance_model_parameters.csv (variance parameterization)\n")
