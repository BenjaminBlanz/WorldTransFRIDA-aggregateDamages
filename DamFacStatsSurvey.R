# Dam Fuc 101125

# Plan:
# 1.We optimize the model under each estimation technique (estimator) using Hannan-Quin (HQ) and BIC. 
#   We may use AIC only as a base, because the punishment factor does not change with the size of the data set n. 
# 2.We use the different estimates for robustness checks.
# 3.The "best" estimation technique can only be derived by argument and convincing robustness checks.

DamFuc <- readRDS("data/asRegDF.RDS")

# Packages
for (p in c("plm","dplyr","sandwich","lmtest","flextable","stargazer")) {if (!requireNamespace(p, quietly = TRUE)) install.packages(p)}
if (!requireNamespace("gsynth", quietly = TRUE)) install.packages("gsynth")
library(gsynth); library(plm); library(dplyr); library(sandwich); library(lmtest); library(flextable) ; library(stargazer)

# Panel
pdat <- pdata.frame(DamFuc, index = c("id","year"))

# Small data set
pdat2 <- subset(pdat, !is.na(gdp) & !is.na(l1gdp) & !is.na(sta))

# counter factrual data set
cfdat2 <- pdat2
cfdat2[,grep('sta$',colnames(cfdat2))] <- 0
cfdat <- pdat
cfdat[,grep('sta$',colnames(cfdat))] <- 0
years <- 2020:2149

# Vector of Covariates which can easily be modified.
y <- "gdp"                                  
c <- c("sta","l1gdp")                       # for easy editing. Just add the covariates here <------------------------------

# Function for easy-use in all regressions. We create right-hand-side (rhs) and formula
mk_rhs     <- function(vars) if (length(vars)==0) "1" else paste(vars, collapse = " + ")
mk_formula <- function(y, vars) as.formula(paste(y, "~", mk_rhs(vars)))
fml_base   <- mk_formula(y, c)

result.lst <- list()
damfac.lst <- list()

# POLS (pooled OLS)
result.lst$POLS <- POLS <- plm(fml_base, data = pdat2, model = "pooling")

# POLS.pred <- predict(POLS)
# plot(c(pdat2$year),c(POLS.pred),pch=20)
# points(c(pdat2$year),c(POLS.pred),pch=20)
# points(c(pdat2$year),c(POLS.cfPred),pch=20,col='red')
# plot(c(pdat2$year),c(pdat2$sta),pch=20)
# points(c(pdat2$year),c(pdat2$sta),pch=20)
# POLS.cfPred <- predict(POLS,newdata = cfdat2)
# points(rep(years,each=nrow(POLS.cfpred)),c(POLS.cfpred),pch=20,col='red')
# 
# POLS.diff <- POLS.cfPred-POLS.pred
# plot(c(pdat2$year),c(POLS.diff),pch=20)
# plot(c(pdat2$sta),c(POLS.diff),pch=20)
# 
# POLS.damfac <- POLS.diff/POLS.cfPred
# plot(c(pdat2$year),c(POLS.damfac),pch=20)
# points(c(pdat2$year),c(POLS.damfac),pch=20)
# plot(c(pdat2$sta),c(POLS.damfac),pch=20)



# cluster Robust POLS
RPOLS <- coeftest(POLS, vcov = vcovHC(POLS, method = "arellano", type = "HC1", cluster = "group"))

# Fixed Effects (FE)
result.lst$FE <- FE <- plm(fml_base, data = pdat2, model = "within")

# Random Effects (RE)
result.lst$RE <- RE <- plm(fml_base, data = pdat2, model = "random")

# First Differencing (FD)
result.lst$FD <- FD <- plm(fml_base, data = pdat2, model = "fd")

# first set of plots ####
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

# Twice First Differenced (TFD) ####

## 1) We create twice-differenced columns for y and all covariates in c ####
vars_tfd <- c(y, c)
for (v in vars_tfd) {pdat2[[paste0(v, "_tfd")]] <- diff(pdat2[[v]], differences = 2)}

## 2) We build the TFD vectors ####
y_tfd  <- paste0(y, "_tfd")
ctfd   <- paste0(c, "_tfd")

## 3) We drop rows with NAs, created because of differencing ####
keep <- complete.cases(pdat2[, c(y_tfd, ctfd)])
pdat_tfd <- pdat2[keep, ]

## 4) TFD regression ####
TFD <- plm(as.formula(paste(y_tfd, "~", paste(ctfd, collapse = " + "))),
           data = pdat_tfd, model = "pooling")

# Two-Way Fixed Effects (TFE) ####
TFE <- plm(fml_base, data = pdat2, model = "within", effect = "twoways")

# Interactive Fixed Effects (IFE), Two-Factor Interactive Fixed Effects (IFE2), Three-Factor Interactive Fixed Effects (IFE3)
.pick <- function(x, nm) if (nm %in% names(x)) unname(x[[nm]]) else NA_real_

run_ife <- function(r) {
	fit <- interFE(fml_base, data = as.data.frame(pdat2), index = c("id","year"), r = r,se = TRUE, nboots = 200, 
                                       seed= 123,normalize = FALSE)
	ET  <- fit$est.table
	cf  <- setNames(ET[, "Coef"], rownames(ET))
	se  <- setNames(ET[, "S.E."], rownames(ET))
	  
	# vectors in the order of intercept then covariates
	cf_vec <- c(.pick(cf, "_const"), sapply(c, function(v) .pick(cf, v)))
	se_vec <- c(.pick(se, "_const"), sapply(c, function(v) .pick(se, v)))
  
	list(model = fit, coef = cf, se = se, coef_vec = cf_vec, se_vec = se_vec)
}

# We fit the IFEs
IFE1_out <- run_ife(1)
IFE2_out <- run_ife(2)
IFE3_out <- run_ife(3)

IFE  <- IFE1_out$model
IFE2 <- IFE2_out$model
IFE3 <- IFE3_out$model

# Coefficient and SEs for Stargazer
ife_coef_vec  <- IFE1_out$coef_vec
ife_se_vec    <- IFE1_out$se_vec

ife2_coef_vec <- IFE2_out$coef_vec
ife2_se_vec   <- IFE2_out$se_vec

ife3_coef_vec <- IFE3_out$coef_vec
ife3_se_vec   <- IFE3_out$se_vec

# Rolling Fixed Effects (RFE)

RFE <- rfe_mde(formula = fml_base, data = pdat2, index = c("id","year"), j = 5,vcov_type = "HC1", cluster = "group")

# helper to safe getter
g <- function(x, nm) if (!is.null(x) && nm %in% names(x)) unname(x[[nm]]) else NA_real_

# pulls fields in a way that works for both return shapes:
get_first <- function(...) { L <- list(...); for (z in L) if (!is.null(z)) return(z); NULL }

# rfe output 
beta_raw <- get_first(RFE$beta_rfe, RFE$beta)
se_raw   <- get_first(RFE$se_rfe,   RFE$se)

# coefficient names to align on
coef_names <- if (!is.null(RFE$coef_names)) RFE$coef_names else names(beta_raw)

# We name the vectors 
if (!is.null(beta_raw) && !is.null(coef_names)) {rfe_coef_named <- setNames(beta_raw, coef_names)
} else {rfe_coef_named <- beta_raw}
if (!is.null(se_raw) && !is.null(coef_names)) {rfe_se_named <- setNames(se_raw, coef_names)
} else {rfe_se_named <- se_raw}

# We build vectors for stargazer
rfe_coef_vec <- c(NA_real_, sapply(c, function(v) g(rfe_coef_named, v)))
rfe_se_vec   <- c(NA_real_, sapply(c, function(v) g(rfe_se_named,   v)))

rfe_coef_vec <- unname(rfe_coef_vec)
rfe_se_vec   <- unname(rfe_se_vec)

# Table 
stargazer(POLS, RPOLS, FE, RE, FD, TFD, TFE, POLS, POLS, POLS, POLS,
          type = "text",
          title = "Vergleich von Panel-Verfahren",
          dep.var.labels.include = FALSE,
          column.labels = c("POLS","RPOLS","FE","RE","FD","TFD","TFE","IFE1","IFE2","IFE3","RFE"),
          coef = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL, 
                      ife_coef_vec, ife2_coef_vec, ife3_coef_vec, rfe_coef_vec),
          se   = list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                      ife_se_vec,   ife2_se_vec,   ife3_se_vec, rfe_se_vec),
          keep.stat = c("n"),
          align = TRUE, no.space = TRUE, header = FALSE)

# Rolling Two-Way Fixed Effects (RTFE)

