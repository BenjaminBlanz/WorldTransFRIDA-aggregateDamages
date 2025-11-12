# Dam Fuc Gsynth small Version 

DamFuc <- readRDS("data/asRegDF.RDS")

# Packages ##########################################################################################

# for (p in c("plm","dplyr","sandwich","lmtest","flextable","stargazer","caret","future.apply","doParallel","foreach","scales","panelView")) {if (!requireNamespace(p, quietly = TRUE)) install.packages(p)}
# if (!requireNamespace("gsynth", quietly = TRUE)) install.packages("gsynth")
library(scales);library(gsynth); library(plm); library(dplyr); library(sandwich); library(lmtest); library(caret); library(future.apply);
library(doParallel); library(foreach); library(panelView)

# Filter for Collinear Variables (one of two is randomly dropped) ###################################

num_ok  <- sapply(DamFuc, is.numeric) & sapply(DamFuc, sd, na.rm = TRUE) > 0 # We keep numeric vectors and sd greater zero
R       <- cor(DamFuc[, num_ok, drop = FALSE], use = "pairwise.complete.obs")
drops   <- setdiff(findCorrelation(R, cutoff = 0.90, names = TRUE), # (We can increase the cutoff value to 0.95 (cutoff=0.99) if we want more variables in the data frame.)
                   c("gdp","l1gdp","sta","id","year")) # If we want to keep more variables for sure, we can add those here.
DamFuc  <- DamFuc[, setdiff(names(DamFuc), drops), drop = FALSE]

DamFuc$sta2<-DamFuc$sta^2
DamFuc$sta3<-DamFuc$sta^3
DamFuc$l1gdpXsta <- DamFuc$l1gdp * DamFuc$sta

# We create a dummy for an arbitrary scenario, i.e. scenario 1

DamFuc$D<-ifelse(DamFuc$year>=2100& DamFuc$id=='1',1,0)

# We sort DamFuc. D has to be in the front.

DamFuc <- DamFuc[, c("D", setdiff(names(DamFuc), "D"))]

# Panel #############################################################################################

stopifnot(all(c("id","year") %in% names(DamFuc)))     # We make sure, indices exist.
pdata   <- pdata.frame(DamFuc, index = c("id","year"))

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

# Interactive Fixed Effects (IFE) with gsynth ################################################################################

IFE_fit <- gsynth(
  formula = fml_base,            
  data    = as.data.frame(DamFuc),                   # I plugged in DamFuc not pdata because pdata is already a panel
  index   = c("id","year"),
  force   = "two-way",                 # this is for the interactive part
  r       = c(0, 50),                   # I set max to 5 but in Levante we can chose 10 such that it is searched for more interactive fixed effects
  se      = TRUE, 
  inference = "jackknife",             # jackknife is the most robust option among other possible methods
  CV      = TRUE,                      # We use cross validation to find the best r 
  criterion = "mspe",                 # We could also use the Bai Ng Information criterion 2002 IC1 or IC3 which are least or most punishing.
  nboots = 200,                        # I chose the number of bootstraps 100 but we can increase it in Levante
  parallel = TRUE,
  estimator = "ife",
  cores    = max(1L, parallel::detectCores() - 1L)
)

# Plot ####

