# Dam Fuc 121125 b POLS only

# Load Data ####

# Packages ##########################################################################################

for (p in c("plm","dplyr","sandwich","lmtest","flextable","stargazer","caret","future.apply","doParallel","foreach")) {if (!requireNamespace(p, quietly = TRUE)) install.packages(p)}
if (!requireNamespace("gsynth", quietly = TRUE)) install.packages("gsynth")
library(gsynth); library(flextable)
library(plm); library(dplyr); library(sandwich); library(lmtest); library(stargazer); library(caret); library(future.apply)
library(doParallel); library(foreach)

# Filter for Collinear Variables (one of two is randomly dropped) ###################################

# (We can increase the cutoff value to 0.95 (cutoff=0.99) if we want more variables in the data frame.)
cutoffPar <- 0.2

DamFuc <- readRDS("data/asRegDF.RDS")
DamFuc  <- DamFuc[,-which(colnames(DamFuc)%in%c('gdpGro','gdp_future_recession_counter','gdppcGroRt','gdp_nominal_gdp_growth_rate','ccs_storing_co2'))]

num_ok  <- sapply(DamFuc, is.numeric) & sapply(DamFuc, sd, na.rm = TRUE) > 0 # We keep numeric vectors and sd greater zero
R       <- cor(DamFuc[, num_ok, drop = FALSE], use = "pairwise.complete.obs")
drops   <- setdiff(findCorrelation(R, cutoff = cutoffPar, names = TRUE), 
                   c("gdp","l1gdp","sta","id","year","l1gdp","l1sta",
                   	"l2gdp","l2sta","l3gdp","l3sta","l10gdp","l10sta",
                   	"l20gdp","l20sta","l30gdp","l30sta","l40gdp","l40sta")) # If we want to keep more variables for sure, we can add those here.
DamFuc  <- DamFuc[, setdiff(names(DamFuc), drops), drop = FALSE]

DamFuc$l1gdpXsta <- DamFuc$l1gdp * DamFuc$sta 
DamFuc$l2gdpXsta <- DamFuc$l2gdp * DamFuc$sta 
DamFuc$l3gdpXsta <- DamFuc$l3gdp * DamFuc$sta 
DamFuc$l1gdpXl1sta <- DamFuc$l1gdp * DamFuc$l1sta 
DamFuc$l1gdpXl2sta <- DamFuc$l1gdp * DamFuc$l2sta 
DamFuc$l1gdpXl3sta <- DamFuc$l1gdp * DamFuc$l3sta
DamFuc$l1gdp2 <- DamFuc$l1gdp^2
DamFuc$l1gdp3 <- DamFuc$l1gdp^3
DamFuc$sta2 <- DamFuc$sta^2
DamFuc$sta3 <- DamFuc$sta^3
DamFuc$l1sta2 <- DamFuc$l1sta^2
DamFuc$l1sta3 <- DamFuc$l1sta^3
DamFuc$l2sta2 <- DamFuc$l2sta^2
DamFuc$l2sta3 <- DamFuc$l2sta^3
DamFuc$l3sta2 <- DamFuc$l3sta^2
DamFuc$l3sta3 <- DamFuc$l3sta^3
DamFuc$l1gdpXl1sta2 <- DamFuc$l1gdp * DamFuc$l1sta2
DamFuc$l1gdpXl2sta2 <- DamFuc$l1gdp * DamFuc$l2sta2 
DamFuc$l1gdpXl3sta2 <- DamFuc$l1gdp * DamFuc$l3sta2

cat(sprintf('After kicking out variables with multicolinearity of %0.1f%% or below we are left with %i variables.\n',
						cutoffPar*100,ncol(DamFuc)-3))

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

anchor_vars <- intersect(c, c("sta","sta2","l1gdp")) # These variables are always part of the model. We can include those we want to have for sure.
pool        <- setdiff(c, anchor_vars)
max_k <- min(6L, length(pool))  # I set it to 3, but in Levante we can set it higher. We have to be careful, because of the exponential combinatorics.

rhs_list <- unique(c(if (length(anchor_vars)) paste(anchor_vars, collapse = " + ") else NULL, unlist(lapply(1:max_k, function(k)
    apply(combn(pool, k), 2, function(cols) paste(c(anchor_vars, cols), collapse = " + "))), use.names = FALSE)))

cat(sprintf('Allowing models with up to %i covariates, we have a list of %i candidate models\n',
						max_k, length(rhs_list)))

# Cluster Start and Stop
library(parallel)
n_workers <- detectCores()
chunkSize <- 100

# We fit all candidates models in parallel and compute BIC

parfun <- function(i){
	rhs <- rhs_list[[i]]
	sucess <- F
	tryCatch({
		m <- plm(as.formula(paste(y, "~", rhs)), data = pdata, model = "pooling") # We can change the model to FE or RE etc. and adjust it to the other methods
		r <- residuals(m)
		n <- sum(!is.na(r))
		rss <- sum(r^2, na.rm = TRUE)
		tss <- sum((pdata[[y]] - mean(pdata[[y]]))^2)
		k <- length(coef(m))
		# list(rhs = rhs, model = m, bic = n*log(rss/n) + k*log(n))
		bic <- -n*log(rss/n)+k*log(n)
		r2 <- 1-rss/tss
		adjr2 <- 1-((1-r2)*(n-1))/(n-k-1)
		sucess <- T
		},  
		error = function(e){})
	if (sucess) {
		return(list(i=i,bic=bic,r2=r2,adjr2=adjr2))
	} else {
		return(NULL)
	}
}

cl <- makeForkCluster(min(detectCores(),length(rhs_list)/chunkSize))
fits <- parLapplyLB(cl,1:length(rhs_list),parfun,chunk.size = 100)
stopCluster(cl)

# doParallel::registerDoParallel(cl)
# fits <- foreach(rhs = rhs_list, .packages = c("plm")) %dopar% {
#   m <- tryCatch(
#     plm(as.formula(paste(y, "~", rhs)), data = pdata, model = "pooling"),  # We can change the model to FE or RE etc. and adjust it to the other methods
#     error = function(e) NULL)
#   if (is.null(m)) return(NULL)
#   r <- residuals(m); n <- sum(!is.na(r)); rss <- sum(r^2, na.rm = TRUE); k <- length(coef(m))
#   list(rhs = rhs, model = m, bic = n*log(rss/n) + k*log(n))
# }
# parallel::stopCluster(cl)


# We keep the successful fits 
fits <- Filter(Negate(is.null), fits)
stopifnot(length(fits) > 0L)

# The "best of all" models is selected
bics <- vapply(fits, `[[`, numeric(1), "bic")
r2s <- vapply(fits, `[[`, numeric(1), "r2")
adjr2s <- vapply(fits, `[[`, numeric(1), "adjr2")
best.idx.bic <- which.min(bics)
best.idx.bic <- fits[[best.idx.bic]]$i
best.idx.r2 <- which.max(r2s)
best.idx.r2 <- fits[[best.idx.r2]]$i
best.idx.adjr2 <- which.max(adjr2s)
best.idx.adjr2 <- fits[[best.idx.adjr2]]$i
bestModels <- list()
bestModels$r2 <- plm(as.formula(paste(y, "~", rhs_list[[best.idx.r2]])), data = pdata, model = "pooling")
bestModels$adjr2 <- plm(as.formula(paste(y, "~", rhs_list[[best.idx.adjr2]])), data = pdata, model = "pooling")
bestModels$bic <- plm(as.formula(paste(y, "~", rhs_list[[best.idx.bic]])), data = pdata, model = "pooling")

for(i in 1:length(bestModels)){
	cat(sprintf("Best RHS by %s: %s\n", names(bestModels)[i], as.character(bestModels[[i]]$formula)[3]))
}
# Plot ################
figFolder<-file.path('figures','damFacStatsSurvey',paste0('upTo',max_k,'Covariates'))
dir.create(figFolder,showWarnings = F,recursive = T)
cfdata <- pdata
cfdata[,grep('sta$',colnames(cfdata))] <- 0
cfdata[,grep('sta2$',colnames(cfdata))] <- 0
cfdata[,grep('sta3$',colnames(cfdata))] <- 0
damfacs <- list()
for(i in 1:length(bestModels)){
	name <- names(bestModels)[i]
	years <- 2020:2149
	pred <- predict(bestModels[[i]])
	cfPred <- predict(bestModels[[i]],newdata = cfdata)
	moddiff <- cfPred-pred
	damfacs[[name]] <- moddiff/cfPred
	png(file.path(figFolder,paste0('modelPOLS-',name,'.png')),400,400)
	plot(c(pdata$sta),c(damfacs[[name]]),pch=20,
	     main=name)
	dev.off()
}
