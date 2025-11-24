# DamFac IFE with interFE() r=70 with lines through dots

DamFuc <- readRDS("/Users/CanKaraarslan/Library/Mobile Documents/com~apple~CloudDocs/Benjamin Blanz/Benjamin Data Set 2/1756737438341_asRegDF.RDS")
# DamFuc <- readRDS("1756737438341_asRegDF.RDS")

# Packages ##########################################################################################

# for (p in c("plm","dplyr","sandwich","lmtest","flextable","stargazer","caret","future.apply","doParallel","foreach","scales","panelView")) {if (!requireNamespace(p, quietly = TRUE)) install.packages(p)}
# if (!requireNamespace("gsynth", quietly = TRUE)) install.packages("gsynth")
library(scales);library(gsynth); library(plm); library(dplyr); library(sandwich); library(lmtest); library(caret); library(future.apply);
library(doParallel); library(foreach); library(panelView)

# Filter for Collinear Variables (one of two is randomly dropped) ###################################

num_ok  <- sapply(DamFuc, is.numeric) & sapply(DamFuc, sd, na.rm = TRUE) > 0 # We keep numeric vectors and sd greater zero
R       <- cor(DamFuc[, num_ok, drop = FALSE], use = "pairwise.complete.obs")
drops   <- setdiff(findCorrelation(R, cutoff = 0.10, names = TRUE), # (We can increase the cutoff value to 0.95 (cutoff=0.99) if we want more variables in the data frame.)
                   c("gdp","l1gdp","sta","id","year")) # If we want to keep more variables for sure, we can add those here.
DamFuc  <- DamFuc[, setdiff(names(DamFuc), drops), drop = FALSE]

DamFuc$sta2<-DamFuc$sta^2
DamFuc$sta3<-DamFuc$sta^3
DamFuc$l1gdpXsta <- DamFuc$l1gdp * DamFuc$sta
DamFuc$l1gdp2Xsta2 <- DamFuc$l1gdp^2 * DamFuc$sta^2
DamFuc$l1gdp3Xsta3 <- DamFuc$l1gdp^3 * DamFuc$sta^3

# Dependent Variable and Vector of Covariates 

y <- "gdp"
c <- setdiff(names(DamFuc), c("id","year","ccs_storing_co2","emissions_so2_emissions_from_food_and_land_use",
                              "sea_level_rise_costs_and_impacts_total_slr_protection_costs",
                              "terrestrial_carbon_balance_annual_carbon_uptake_in_peatlands", y))

# Function for easy-use in all regressions. We create right-hand-side (rhs) and formula #############

mk_rhs     <- function(vars) if (length(vars)==0) "1" else paste(vars, collapse = " + ")
mk_formula <- function(y, vars) as.formula(paste(y, "~", mk_rhs(vars)))
fml_base   <- mk_formula(y, c)

### Ich habe die Modelwahl via BIC erstmal ausgelassen ------------------------------------

# Interactive Fixed Effects (IFE) with interFE (not gsynth) ################################################################################

IFE_fit <- interFE(
  formula = fml_base,
  data    = as.data.frame(DamFuc),
  index   = c("id", "year"),
  r       = 70,           # We take 70 interactive fixed effects but there is not much difference.
  force   = "two-way",   
  se      = FALSE,        # bootstrap standard errors put off for faster computation
  nboots  = 100,         
  normalize = FALSE      # We keep scale as in in our data
)

# Plot ####
## Counterfactual data ####

cfdat <- DamFuc
sta_cols <- grepl("sta", names(cfdat), ignore.case = TRUE) # We detect all sta s in the data
if (any(sta_cols)) {cfdat[, sta_cols] <- 0}                # We set sta s equal to zero

## Prediction ####
# We define a predictor function because we cannot use predict() for a gsynth object
## Predictor Function ####

predict_IFE <- function(model, data, formula, id_col = "id", time_col = "year") { # We assign the prediction function to "predict_IFE"
  mf <- model.frame(formula, data = data, na.action = na.omit)                    # We build the model frame, without missing values. Out fml_base will be used for any data we provide.
  rows_used <- as.integer(rownames(mf))                                           # We convert the character vector of row names into integers
  X <- model.matrix(formula, data = mf)                                           # We build the model matrix
  b <- as.numeric(model$beta)                                                     # We make sure, our coefficient estimates, from gsynth and IFE are numeric, ready to be multiplied.
  stopifnot(length(b) > 0)                                                        # We make sure, coefficients exist - at least one.
  int_col <- which(colnames(X) == "(Intercept)")                                  # Finds the position of the intercept term.
  if (length(int_col)) X <- X[, -int_col, drop = FALSE]                           # We drop our intercept, if we have one, such that matrix multiplication with IFE object works,
  p <- min(ncol(X), length(b))                                                    # We ensure that the number of our columns and coefficients is equal.
  X <- X[, seq_len(p), drop = FALSE]                                              # We keep the first p columns of the matrix and ensure it stays a matrix. (X and b are compatible)
  b <- b[seq_len(p)]                                                              # We keep the first p coefficients from b, thus, X and b are compatible.
  pred <- as.numeric(X %*% b)                                                     # We multiply (n x p) matrix X with (p x 1) vector b and store the (n x 1) result as numeric.
  data.frame(                                                                     # Our data frame with the predicted values is created.
    id   = data[[id_col]][rows_used],
    year = data[[time_col]][rows_used],
    pred = pred,
    stringsAsFactors = FALSE)} 

## Prediction of Observed and Counterfactual ####

pred_obs_df <- predict_IFE(IFE_fit, DamFuc, fml_base) # Prediction of Observed

pred_cf_df  <- predict_IFE(IFE_fit, cfdat,  fml_base) # Prediction of Counterfactual

## Difference relative to Counterfactual, Relative Damage ####
key_obs <- paste(pred_obs_df$id, pred_obs_df$year, sep = "_") # We combine the vectors by elements and create id_year for the observations.
key_cf <- paste(pred_cf_df$id, pred_cf_df$year, sep = "_")    # Same for the counterfactuals.
key_dam <- paste(DamFuc$id, DamFuc$year, sep = "_")           

idx_cf <- match(key_obs, key_cf)      # We match the indices of the observed units with the counterfactual.
idx_sta <- match(key_obs, key_dam)    # Find for every i in key_obs the position of it in key_dam.
sta_vec <- DamFuc$sta[idx_sta]        # We take the sta values from DamFac at the positions indicated in idx_sta, thus, where key_obs and key_dam intersect. 
                                      # So, we reorder the sta values such that these fit the order of the predicted values.
# Damage Factor
damage <- (pred_cf_df$pred[idx_cf] - pred_obs_df$pred) / pred_cf_df$pred[idx_cf] 

ok <- is.finite(sta_vec) & is.finite(damage) & (pred_cf_df$pred[idx_cf] != 0)

### Final Scatter Plot ####
 
par(family = "Times New Roman", mgp = c(2, 0.8, 0))
plot(pred_merged$sta[ok], damage[ok], pch = 20,
     xlab = "Surface Temperature Anomaly",
     ylab = "Damage = (Cf - Obs)/CF",
     main = "Damage Factor versus Tempertaure (IFE), r=70")

# Lines through dots

df_plot <- data.frame(
  id     = pred_merged$id[ok],
  sta    = pred_merged$sta[ok],
  damage = damage[ok]
)
scenario_ids <- unique(df_plot$id)              # We take all scenario id s that appear in our data.

for (scenario_id in scenario_ids) {             # We take each scenario and order it by temperature.
  tmp <- df_plot[df_plot$id == scenario_id, ]   # We take a single scenario.
  tmp <- tmp[order(tmp$sta), ]                  # We order by temperature.
  lines(tmp$sta, tmp$damage,
        col = "red", lwd = 0.6)                 # We draw red lines through the dots.
}
