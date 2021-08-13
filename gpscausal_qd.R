library(devtools)
require(doParallel)
library(data.table)
library("parallel")
require(xgboost)
require(dplyr)
require(tidyr)
#try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
#install_github("fasrc/CausalGPS", ref="master", force = TRUE)
require(CausalGPS)
library(fst)
library(data.table)
library("mgcv")
library("gnm")
require(dplyr)
options(stringsAsFactors = FALSE)
require(parallel)
require(KernSmooth)
library(fst)
library("parallel")
require(ggplot2)
require(cowplot)
require(ggExtra)
set.seed(1)
matching_smooth<-function(pseudo.out=dose.response.mean,
                          a,
                          bw.seq=seq(1,1,length.out=10),
                          a.vals){
  kern <- function(t){ dnorm(t) }
  w.fn <- function(bw){ w.avals <- NULL; for (a.val in a.vals){
    a.std <- (a-a.val)/bw; kern.std <- kern(a.std)/bw
    w.avals <- c(w.avals, mean(a.std^2*kern.std)*(kern(0)/bw) /
                   (mean(kern.std)*mean(a.std^2*kern.std)-mean(a.std*kern.std)^2))
  }; return(w.avals/length(a)) }
  hatvals <- function(bw){ approx(a.vals,w.fn(bw),xout=a,rule=2)$y }
  cts.eff.fn <- function(out,bw){
    approx(locpoly(a,out,bandwidth=bw, gridsize=1000),xout=a,rule=2)$y }
  # note: choice of bandwidth range depends on specific problem,
  # make sure to inspect plot of risk as function of bandwidth
  risk.fn <- function(h){ hats <- hatvals(h); mean( ((pseudo.out - cts.eff.fn(pseudo.out,bw=h))/(1-hats))^2) }
  risk.est <- sapply(bw.seq,risk.fn);
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw=bw.seq, risk=risk.est)
 
  est <- approx(locpoly(a,pseudo.out,bandwidth=h.opt),xout=a.vals)$y
  return(est)
}

dir_data = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'


#All
load("/nfs/home/P/prd789/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_qd.RData")
covariates_qd$year<-as.factor(covariates_qd$year)
covariates_qd$region<-as.factor(covariates_qd$region)
a.vals <- seq(min(covariates_qd$pm25_ensemble), max(covariates_qd$pm25_ensemble), length.out = 100)
delta_n <- (a.vals[2] - a.vals[1])

load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd$year<-as.factor(aggregate_data_qd$year)
aggregate_data_qd$region<-as.factor(aggregate_data_qd$region)

dead_personyear<-aggregate(cbind(aggregate_data_qd$dead,
                                 aggregate_data_qd$time_count),
                           by=list( aggregate_data_qd$zip,
                                    aggregate_data_qd$year),
                           FUN=sum)
colnames(dead_personyear)[1:4]<-c("zip","year","dead","time_count")
dead_personyear[,"mortality"] <- dead_personyear[,"dead"]/dead_personyear[,"time_count"]

prematch_data <- merge(dead_personyear, covariates_qd,
                       by=c("zip", "year"))

match_pop_all1 <- generate_pseudo_pop(Y=prematch_data$mortality,
                                     w=prematch_data$pm25_ensemble,
                                     c=prematch_data[, c(7:22)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = FALSE,
                                     #transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list(xgb_nrounds=c(50)),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.01,0.99), # trimed, you can change
                                     optimized_compile = FALSE, #created a column counter for how many times matched,
                                     max_attempt = 1,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0,
                                     set_logger("DEBUG"))

match_pop_data1<-match_pop_all1$pseudo_pop
a.vals <- seq(min(match_pop_data1$w), max(match_pop_data1$w), length.out = 100)
erf1<-estimate_npmetric_erf(matched_Y=match_pop_data1$Y,
                           matched_w=match_pop_data1$w,
                           bw_seq= seq(8*delta_n,40*delta_n,2*delta_n),
                           w_vals= a.vals, nthread=10)
plot(erf1)
plot(a.vals, erf1$erf/erf1$erf[1])
save(match_pop_data1, erf1, file=paste0(dir_data, "erf_replicate.RData"))
