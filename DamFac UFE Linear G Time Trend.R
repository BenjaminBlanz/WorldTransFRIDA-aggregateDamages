# DamFac Diagnostics One-way (Unit) Fixed Effects 021225, plm()
# with Linear Global Time Trend

# Packages ##########################################################################################

# for (p in c("plm","dplyr","sandwich","lmtest","flextable","stargazer","caret","future.apply","doParallel","foreach","scales","panelView","tidyverse")) {if (!requireNamespace(p, quietly = TRUE)) install.packages(p)}
# if (!requireNamespace("gsynth", quietly = TRUE)) install.packages("gsynth")
library(scales);library(gsynth); library(plm); library(dplyr); library(sandwich); library(lmtest); library(caret); library(future.apply);
library(doParallel); library(foreach); library(panelView); library(tidyverse)

# Merge Data ####

path <- "/Users/CanKaraarslan/Library/Mobile Documents/com~apple~CloudDocs/Benjamin Blanz/Benjamin Data 3/temperaturID318"

## Single Files ####

temp_raw <- readRDS(file.path(path,"energy_balance_model_surface_temperature_anomaly.RDS"))
gdp_raw <- readRDS(file.path(path,"gdp_real_gdp_in_2021c.RDS"))
infl_raw <- readRDS(file.path(path,"inflation_inflation_rate.RDS"))

## Format from Wide to Long Panel ####

temp_long <- temp_raw %>% pivot_longer(cols = -id, names_to  = "year",values_to = "sta")
gdp_long <- gdp_raw %>% pivot_longer(cols = -id, names_to  = "year", values_to = "gdp")
infl_long <- infl_raw %>% pivot_longer(cols = -id, names_to  = "year", values_to = "inflation")

temp_long <- temp_long %>% mutate(year = as.integer(year)) # We make year numeric.
gdp_long <- gdp_long  %>% mutate(year = as.integer(year))
infl_long <- infl_long %>% mutate(year = as.integer(year))

## We merge everything into DamFuc and order by id and year

DamFuc <- temp_long %>% left_join(gdp_long, by = c("id","year")) %>% left_join(infl_long, by = c("id","year")) %>% 
  arrange(id, year)

# interFE() doesn't work with missing data, thus, I take a random draw without missing values

set.seed(123)
DamFuc <- DamFuc %>% drop_na()
ids_to_keep <- DamFuc %>% 
  distinct(id) %>%            # We take distinct id s otherwise rows are selected
  slice_sample(n = 100) %>%   # We take 1000 scenarios randomly
  pull(id)                    # We convert the tibble into a vector for filtering in the next step

DamFuc <- DamFuc %>% filter(id %in% ids_to_keep) # We take only the randomly drawn ids. 

# Filter for Collinear Variables (one of two is randomly dropped) ###################################

num_ok  <- sapply(DamFuc, is.numeric) & sapply(DamFuc, sd, na.rm = TRUE) > 0 # We keep numeric vectors and sd greater zero
R       <- cor(DamFuc[, num_ok, drop = FALSE], use = "pairwise.complete.obs")
drops   <- setdiff(findCorrelation(R, cutoff = 0.90, names = TRUE), # (We can increase the cutoff value to 0.95 (cutoff=0.99) if we want more variables in the data frame.)
                   c("gdp","l1gdp","sta","id","year")) # If we want to keep more variables for sure, we can add those here.
DamFuc  <- DamFuc[, setdiff(names(DamFuc), drops), drop = FALSE]

# New Variables 

DamFuc <- DamFuc %>% group_by(id) %>% mutate(l1gdp = lag(gdp, 1),l2gdp = lag(gdp, 2)) %>% ungroup()
DamFuc$sta2 <- DamFuc$sta^2
DamFuc$sta3 <- DamFuc$sta^3
DamFuc$l1gdpXsta <- DamFuc$l1gdp * DamFuc$sta
DamFuc$l1gdp2Xsta2 <- DamFuc$l1gdp^2 * DamFuc$sta^2
DamFuc$l1gdp3Xsta3 <- DamFuc$l1gdp^3 * DamFuc$sta^3

# Time Trends ####
# For Linear Global Time Trend

DamFuc$t <- DamFuc$year

# For Quadratic Global Time Trend

# DamFuc$t2 <- DamFuc$year^2

# Dependent Variable and Vector of Covariates 

y <- "gdp"
c <- setdiff(names(DamFuc), c("id","year", y))

# Function for easy-use in all regressions. We create right-hand-side (rhs) and formula #############

mk_rhs     <- function(vars) if (length(vars)==0) "1" else paste(vars, collapse = " + ")
mk_formula <- function(y, vars) as.formula(paste(y, "~", mk_rhs(vars)))
fml_base   <- mk_formula(y, c)

### Ich habe die Modelwahl via BIC erstmal ausgelassen ------------------------------------

# Unit Fixed Effects (Name IFE kept only at the moment)

UFE_fit<- plm(
  formula = fml_base,
  data    = as.data.frame(DamFuc),
  index   = c("id", "year"),
  model   = "within",     # within estimation
  effect  = "individual"  # only unit fixed effects
)

# Plot ####
## Counterfactual data ####

cfdat <- DamFuc
sta_cols <- grepl("sta", names(cfdat), ignore.case = TRUE) # We detect all sta s in the data
if (any(sta_cols)) {cfdat[, sta_cols] <- 0}                # We set sta s equal to zero

## Prediction ####
# We define a predictor function because we cannot use predict() for a gsynth object
## Predictor Function ####

predict_UFE <- function(model, data, formula, id_col = "id", time_col = "year") { # We assign the prediction function to "predict_IFE"
  mf <- model.frame(formula, data = data, na.action = na.omit)                    # We build the model frame, without missing values. Out fml_base will be used for any data we provide.
  rows_used <- as.integer(rownames(mf))                                           # We convert the character vector of row names into integers
  X <- model.matrix(formula, data = mf)                                           # We build the model matrix
  b <- as.numeric(coef(model))                                                    # We extract the model parameters from plm
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

pred_obs_df <- predict_UFE(UFE_fit, DamFuc, fml_base) # Prediction of Observed

pred_cf_df  <- predict_UFE(UFE_fit, cfdat,  fml_base) # Prediction of Counterfactual

## Difference relative to Counterfactual, Relative Damage ####

key_obs <- paste(pred_obs_df$id, pred_obs_df$year, sep = "_") # We combine the vectors by elements and create id_year for the observations.
key_cf <- paste(pred_cf_df$id, pred_cf_df$year, sep = "_")    # Same for the counterfactuals.
key_dam <- paste(DamFuc$id, DamFuc$year, sep = "_")           # 

idx_cf <- match(key_obs, key_cf)      # We match the indices of the observed units with the counterfactual.
idx_sta <- match(key_obs, key_dam)    # Find for every i in key_obs the position of it in key_dam.
sta_vec <- DamFuc$sta[idx_sta]        # We take the sta values from DamFac at the positions indicated in idx_sta, thus, where key_obs and key_dam intersect. 
# So, we reorder the sta values such that these fit the order of the predicted values.

# Damage Factor

damage <- (pred_cf_df$pred[idx_cf] - pred_obs_df$pred) / pred_cf_df$pred[idx_cf] 

ok <- is.finite(sta_vec) & is.finite(damage) & (pred_cf_df$pred[idx_cf] != 0)

### Final Scatter Plot ####

par(family = "Times New Roman", mgp = c(2, 0.8, 0))
plot(sta_vec[ok], damage[ok], pch = 20,
     xlab = "Surface Temperature Anomaly",
     ylab = "Damage = (Cf - Obs)/CF",
     main = "Damage Factor versus Tempertaure (UFE) with Linear Global Time Trend")

# Lines through dots

df_plot <- data.frame(
  id     = pred_obs_df$id[ok],
  sta    = sta_vec[ok],
  damage = damage[ok]
)
scenario_ids <- unique(df_plot$id)              # We take all scenario id s that appear in our data.

for (scenario_id in scenario_ids) {             # We take each scenario and order it by temperature.
  tmp <- df_plot[df_plot$id == scenario_id, ]   # We take a single scenario.
  tmp <- tmp[order(tmp$sta), ]                  # We order by temperature.
  lines(tmp$sta, tmp$damage,
        col = "red", lwd = 0.6)                 # We draw red lines through the dots.
}

# Damage Factor Time Trend

mean_damage <- sapply(sort(unique(sta_vec[ok])), function(t) mean(damage[ok][sta_vec[ok] == t], na.rm = TRUE))
lines(sort(unique(sta_vec[ok])), mean_damage, col = "blue", lwd = 4)

# Diagnostics 1 UFE GDP comparison via Spaghetti Plot

obs_clean <- pred_obs_df[is.finite(pred_obs_df$year) & is.finite(pred_obs_df$pred), ] # We keep only the rows where year and predicted gdp is present.
cf_clean  <- pred_cf_df[is.finite(pred_cf_df$year) & is.finite(pred_cf_df$pred), ]    # Also for the counterfactual, we keep only the rows with existent values.

xlim <- range(c(obs_clean$year, cf_clean$year)) # We define the plot limits by taking from both vectors the range, defined between the lowest and the highest values.
ylim <- range(c(obs_clean$pred, cf_clean$pred)) # Plot limits defined as c(min_year, max_year) from observed and counterfactual vectors of predicted gdp.

# Empty Plot
par(family = "Times New Roman", mgp = c(2, 0.8, 0))      # Times New Roman (Font), Labels of Axis a bit closer
plot(xlim, ylim, type = "n",                             # We draw an empty plot, Cartesian plane, first. 
     xlab = "Year",
     ylab = "GDP",
     main = "\nDiagnostics 3: GDP Comparison UFE with Linear Global Time Trend",
     cex.main = 1.8)

ids <- intersect(unique(obs_clean$id), unique(cf_clean$id)) # We take the id s from both columns once and only those that appear in both data frames. (Set Intersection)

## We draw two lines for each scenario id
for (scenario_id in ids) {                               # From the vector of ids we pick recursively one scenario at a time.
  obs_i <- obs_clean[obs_clean$id == scenario_id, ]      
  cf_i  <- cf_clean[cf_clean$id == scenario_id, ]
  obs_i <- obs_i[order(obs_i$year), ]                    # We order the chosen ids by year. 
  cf_i  <- cf_i[order(cf_i$year), ]
  lines(cf_i$year, cf_i$pred, col = "black", lwd = 1)    # Counterfactual gets the black line
  lines(obs_i$year, obs_i$pred, col = "red", lwd = 1)    # Observed gets the red line
}

## Diagnostics Time trends ####

lines(years <- sort(unique(obs_clean$year)), 
      sapply(years, \(yy) mean(obs_clean$pred[obs_clean$year == yy], na.rm = TRUE)),
      col = "red",   lwd = 4)

lines(years, 
      sapply(years, \(yy) mean(cf_clean$pred[cf_clean$year == yy], na.rm = TRUE)),
      col = "black", lwd = 4)

## Legend ####

legend("topleft", legend = c("Observed", "Counterfactual"), # Label
       col = c("red", "black"),
       lwd = 2,
       bty = "n")
