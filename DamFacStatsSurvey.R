# Dam Fuc 121125 b POLS only

# Load Data ####

DamFuc <- readRDS("/Users/CanKaraarslan/Library/Mobile Documents/com~apple~CloudDocs/Benjamin Blanz/Benjamin Data Set 2/1756737438341_asRegDF.RDS")

# Packages ##########################################################################################

for (p in c("plm","dplyr","sandwich","lmtest","flextable","stargazer","caret","future.apply","doParallel","foreach")) {if (!requireNamespace(p, quietly = TRUE)) install.packages(p)}
if (!requireNamespace("gsynth", quietly = TRUE)) install.packages("gsynth")
library(gsynth); library(plm); library(dplyr); library(sandwich); library(lmtest); library(flextable); library(stargazer); library(caret); library(future.apply)
library(doParallel); library(foreach)

# Filter for Collinear Variables (one of two is randomly dropped) ###################################

num_ok  <- sapply(DamFuc, is.numeric) & sapply(DamFuc, sd, na.rm = TRUE) > 0 # We keep numeric vectors and sd greater zero
R       <- cor(DamFuc[, num_ok, drop = FALSE], use = "pairwise.complete.obs")
drops   <- setdiff(findCorrelation(R, cutoff = 0.90, names = TRUE), # (We can increase the cutoff value to 0.95 (cutoff=0.99) if we want more variables in the data frame.)
                   c("gdp","l1gdp","sta","id","year")) # If we want to keep more variables for sure, we can add those here.
DamFuc  <- DamFuc[, setdiff(names(DamFuc), drops), drop = FALSE]

# Panel #############################################################################################

stopifnot(all(c("id","year") %in% names(DamFuc)))     # We make sure, indices exist.
pdata   <- pdata.frame(DamFuc, index = c("id","year"))

# Dependent Variable and Vector of Covariates 

y <- "gdp"
c <- setdiff(names(DamFuc), c("id","year", y))

# Function for easy-use in all regressions. We create right-hand-side (rhs) and formula #############

mk_rhs     <- function(vars) if (length(vars)==0) "1" else paste(vars, collapse = " + ")
mk_formula <- function(y, vars) as.formula(paste(y, "~", mk_rhs(vars)))
fml_base   <- mk_formula(y, c)

# Model Selection via BIC #########################################################################

anchor_vars <- intersect(c, c("sta","l1gdp")) # These variables are always part of the model. We can include those we want to have for sure.
pool        <- setdiff(c, anchor_vars)
max_k <- min(3L, length(pool))  # I set it to 3, but in Levante we can set it higher. We have to be careful, because of the exponential combinatorics.

rhs_list <- unique(c(if (length(anchor_vars)) paste(anchor_vars, collapse = " + ") else NULL, unlist(lapply(1:max_k, function(k)
    apply(combn(pool, k), 2, function(cols) paste(c(anchor_vars, cols), collapse = " + "))), use.names = FALSE)))

# Cluster Start and Stop

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)

# We fit all candidates models in parallel and compute BIC

fits <- foreach(rhs = rhs_list, .packages = c("plm")) %dopar% {
  m <- tryCatch(
    plm(as.formula(paste(y, "~", rhs)), data = pdata, model = "pooling"),  # We can change the model to FE or RE etc. and adjust it to the other methods
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  r <- residuals(m); n <- sum(!is.na(r)); rss <- sum(r^2, na.rm = TRUE); k <- length(coef(m))
  list(rhs = rhs, model = m, bic = n*log(rss/n) + k*log(n))
}

parallel::stopCluster(cl)

# We keep the successful fits 
fits <- Filter(Negate(is.null), fits)
stopifnot(length(fits) > 0L)

# The "best of all" models is selected
bics <- vapply(fits, `[[`, numeric(1), "bic")
best <- fits[[ which.min(bics) ]]
result.lst <- setNames(list(best$model), best$rhs)
cat("Best RHS:", names(result.lst), "  BIC:", min(bics), "\n")

# POLS (pooled OLS) ################################################################################

POLS <- plm(fml_base, data = pdata, model = "pooling")

# Plot ################

for(i in 1:length(result.lst)){
  model <- result.lst[[i]]
  pred <- predict(model)
  cfPred <- predict(model,newdata = cfdat2)
  moddiff <- cfPred-pred
  damfac <- moddiff/cfPred
  damfac.lst[[names(result.lst)[i]]] <- damfac
  plot(c(pdat2$sta),c(damfac),pch=20,
       main=names(result.lst)[i])
}







