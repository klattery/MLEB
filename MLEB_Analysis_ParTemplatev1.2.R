#################################################################################
#   Template for Machine Learning Emp Bayes                                     #
#   Built using MLEB_Functions_Finalv1.4                                        #
#                                                                               #
#   Steps 1, 2, & 3 must be specified                                           #
#     1) Directories and Read in Data                                           #
#     2) Coding                                                                 #
#     3) Constraints, Data with Id_Task_Dep, rowfilter                          #
#        4) Initial x0 if unconstrained aggregate solution violates constraints #
#################################################################################
dir_mleb <- "C:/Users/k.lattery/OneDrive - SKIM/Development/MLEB/R_Parallel/"
source(paste0(dir_mleb, "MLEB_Functions_FinalPar_v1.2.R")) # Install Functions

# install.packages
library(readxl) # Optional to read Excel files
library(doParallel) # Also

# Setup multiple cores for MLEB
ncores <- max(detectCores() -1,1) # CHECK ncores
if (file.exists("k_multi_core")){
  stopCluster(k_multi_core)
  remove(k_multi_core)
}
k_multi_core <- makeCluster(ncores)
registerDoParallel(k_multi_core)

####################################################################
##       1.  SPECIFY DIRECTORIES & READ IN DATA                   ##
####################################################################
dir_data <- "C:/Users/k.lattery/OneDrive - SKIM/Development/Training/MLEB/Data/"
dir_out <- "C:/Users/k.lattery/OneDrive - SKIM/Development/MLEB/SawtoothPrize/"

data_all <- read.csv(paste0(dir_data, file = "Data_Full.csv"), colClasses = c("numeric"))
train <- data_all$Train == 1 # Filter for training data

# Install MLEB Functions
# Below uses sync to Resources >> Tools and Templates >> Tools >> 03 Analysis >> MLEB 
###################################################################
#     END STEP 1                                                  #
###################################################################
####################################################################
##       2.  SPECIFY CODING                                       ##
####################################################################
# For catcode:
# codetype = (1 indicator, 2 dummy, 3 effects default)
# reflev = value to use
indcode_list <- list()
indcode_list[[1]] <- catcode(data_all, 5) 
indcode_list[[2]] <- catcode(data_all, 6) 
indcode_list[[3]] <- catcode(data_all, 7) 
indcode_list[[4]] <- catcode(data_all, 8) 
indcode_list[[5]] <- catcode(data_all, 9) 
indcode_list[[6]] <- ordcode(data_all, 10, cut_pts = c(1,2,3,4,5), varout = "ppd") # Ordinal Code

###################################################################
#     END STEP 2                                                  #
###################################################################

# Run code to compile your list above
make_codefiles(indcode_list) # Create: code_master, indcode, indprior 
write.table(cbind(rownames(code_master), code_master), file = paste0(dir_out, "code_master.csv"), sep = ",", na = ".", row.names = FALSE)
# STOP and Review code_master in dir_out

###################################################################
##    3.  SPECIFY CONSTRAINTS, IDTASKDEP, FILTER                 ##
####################################################################

# Must specify columns for id, task, dep to prep file
###############################################
idtaskdep <- data_all[,c(2,3,11)]  # Specify data with id, task, dep in that order
constrain <- make_con(code_master,
                      col_pos = NULL,
                      col_neg = 18:21,
                      row_rel = NULL)
rowfilter <- data_all$Train == 1 # Logical vector, rowfilter <- TRUE to keep everything

###############################################


###################################################################
#     END STEP 3                                                  #
# No other necessary steps, except possibly Step4                 #
###################################################################
rm(indcode_list) # Free RAM
data_conjoint <- prep_file(idtaskdep[rowfilter,],  indcode[rowfilter,]) # 
rm(idtaskdep, indcode) # Free RAM
nresp <- length(data_conjoint$resp_id) # number of respondents used for control
npar <- ncol(code_master) # number of parameters

#alpha_old <- as.vector(as.matrix(read_xlsx(paste0(dir_data,"Recoding_v5.xlsm"), sheet = "alpha")))
# Previous solution to start from
x0 <- rep(-.02, npar)
check <- check_beta(x0, constrain)

# Base model for Agg Model Unconstrained
model_base <- list(
  func = list(
    pred = PredMNL,
    min = LL_Neg,
    gr = grad_MNL
  ),
  con = as.matrix(constrain), # must be matrix
  #x0 = as.vector(as.matrix(xstart)) # must be vector in interior
  #con = con_trivial(npar), # no constraints
  x0 =  x0 # start at 0
  #x0 = alpha_old
)

agg_result <- agg_solve(data_conjoint, model_base) # Aggregate solution

write.table(agg_result$par, file = paste0(dir_out, "agg_result.csv"), sep = ",", na = ".", row.names = TRUE)
aggbeta_r <- t(agg_result$par %*% t(code_master)) # recoded to all levels
write.table(aggbeta_r, file = paste0(dir_out, "aggbeta_r.csv"), sep = ",", na = ".", row.names = TRUE)

con_check <- check_beta(agg_result$par, constrain) # check if within constraint 
consistent <- con_check$consistent # Is unconstrained model consistent with constraints
# STOP If consistent = FALSE consider recoding or new constraints 
# Check aggbeta_r for conceptual sense

if (!consistent){ # If not consistent run again with constraints
  model_base_con <- model_base
  model_base_con$con <- as.matrix(constrain) # With constraints now
  
  ####################################################################    
  ##    Step 4. ONLY NEEDED IF consistent = FALSE                   ##
  model_base_con$x0 <- rep(-0.02, npar)     ##
  ####################################################################  
  
  agg_result <- agg_solve(data_conjoint, model_base_con)
  aggbeta_r <- t(agg_result$par %*% t(code_master)) # recoded to all levels
  con_check <- check_beta(agg_result$par, constrain) # check if within constraint 
}

# Get Prior Covariance: cov_prior
# For big data set this may take like 30 minutes
JeffPriors <- JeffPrior(agg_result$par, data_conjoint, model_base,
                        ObsFish = consistent, ExpFish = !consistent, score_n = 500)
JeffPriors[JeffPriors < .1] <- .1
if (consistent){ 
  cov_prior <- diag(sqrt(1/JeffPriors$ObsFish)) %*% cov2cor(indprior) %*% diag(sqrt(1/JeffPriors$ObsFish))
} else cov_prior <- diag(sqrt(1/JeffPriors$ExpFish)) %*% cov2cor(indprior) %*% diag(sqrt(1/JeffPriors$ExpFish))


# Specify Standard MLEB Model.
# x0 and alpha are normally the same and both are within your constraints
model_mleb <- list(
  func = list(
    pred = PredMNL,
    min = LL_wPriorPDF,
    gr = grad_MNLwMVN,
    logpdf = logpdf_mvnorm
  ),
  prior = list(
    alpha = agg_result$par,
    cov_inv = solve(cov_prior), # Cov Inverse 
    upper_model = rep(TRUE, length(agg_result$par)),
    scale = 1 # scale factor that we will solve for
  ),
  con = as.matrix(constrain), # must be matrix, 
  x0 = agg_result$par # should be same as prior$alpha
)
save(model_mleb, file = paste0(dir_out, "model_mleb.RData")) # Save result to file

# MLEB Control
mleb_control <- list(cal_resp = min(300,nresp), cal_tasks = 1,
                     hold_resp = min(500, nresp), hold_tasks = 2,
                     maxiter = 20, conv_n = 3,
                     tolerance = .1, hold_poss = TRUE,
                     dir_pdf = dir_out, solveonly = FALSE)

mleb(data_list = data_conjoint, model_list = model_mleb, mleb_control = mleb_control)
load(paste0(dir_out, "mleb_result.RData"))
#mleb_result <- model_env$mleb_result
if (!is.null(dev.list())) dev.off() # Reset plot area
kscale_hist <- mleb_result$kscale_hist
kiter <- 1:sum(!is.na(kscale_hist[,1]))
plot(kiter, kscale_hist[kiter, 3])
best_fit <- max(kscale_hist[kiter, 3])
best_iter <- match(best_fit, kscale_hist[, 3]) # Picks iteration with best rlh
# best_iter <- You Can Override with Manual Choice

# Get betas, back code them, and export
best_betas <- mleb_result$eb_betas_hist[[best_iter]] #coded betas
betas_final_r <- cbind(best_betas[,1:2], best_betas[,-1:-2] %*% t(code_master)) #back coded
write.table(betas_final_r, file = paste0(dir_out, "betas_final_r.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(mleb_result$prior_cov_hist[[best_iter]], file = paste0(dir_out, "best_cov.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(mleb_result$prior_alpha_hist[[best_iter]], file = paste0(dir_out, "best_alpha.csv"), sep = ",", na = ".", row.names = FALSE)

# Save mleb_result
mleb_result$best_iter <- best_iter
mleb_result$code_master <- code_master
save(mleb_result, file = paste0(dir_out, "mleb_result.RData")) # Save result to file

