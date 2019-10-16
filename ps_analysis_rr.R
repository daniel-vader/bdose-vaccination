# Description: Propensity score analysis using boosted CART of the association 
#              between hepatitis B birth dose and up-to-date vaccination 
#              status at 18 months. RR analysis per reviewer request. (Functions
#              adapted from ps_analysis.R file)
# Author:      Daniel Vader
# Date:        7-29-2019

################################################################################
# Load libraries, outside code, and datasets
################################################################################
library(dplyr)
library(twang)
library(survey)
library(tidyr)

source("formats_nis2017.R") # load format function

# load data and apply formats
NIS2017 <- readRDS("NIS2017finalps.rds") %>% apply_nis_formats() #apply formats

### Check balance
balance.vars <- c("AGEGRP", "CHILDNM", "SEX", "CWIC_01", "EDUC1", "FRSTBRN",
                  "INCPORAR_I", "M_AGEGRP2", "MARITAL2", "MOBIL_I", "RACEETHK",
                  "RENT_OWN", "N_PRVR", "PROV_FAC", "REGISTRY", 
                  "VFC_ORDER", "INS_BREAK_I", "INS_STAT2_I")


################################################################################
# Create outcome tables
################################################################################

#NIS2017 <- readRDS("NIS2017ps.rds") %>% apply_nis_formats() #apply formats

# Initialize survey design with ps composite weights
dsgn.cart.ATE <- svydesign(ids=~SEQNUMC, strata=~STRATUM, weights=~cart.ATE, data=NIS2017)

# Set basic survey design
dsgn.svy <- svydesign(ids=~SEQNUMC, strata=~STRATUM, weights=~PROVWT_D, data=NIS2017)


# Setup basic model functions. We have two model structures: simple - which relies
# on propensity score balance - and robust - which addidtionally adjusts for
# covariables in a multivariable glm.

modelrun.rr <- function(x, dsgn, simp=T, verbose=F){
  if(simp){
    f <- as.formula(paste(x, "~ birthhep"))
  } else{
    f <- as.formula(paste(x, "~ birthhep + AGEGRP + SEX + CHILDNM + CWIC_01 + EDUC1 + 
                          FRSTBRN + INCPORAR_I + M_AGEGRP2 + MARITAL2 + 
                          MOBIL_I + RACEETHK + RENT_OWN + N_PRVR + 
                          PROV_FAC + REGISTRY + VFC_ORDER +
                          INS_STAT2_I"))
  }
  m <- svyglm(f, design=dsgn, family=quasipoisson(link="log"))
  ci <- exp(confint(m, level=0.95))
  coef <- exp(coef(m))
  
  if(verbose){
    print(summary(m))
    print(m)
    print(ci)
  }
  
  return(c(x, coef[2], ci[2,1], ci[2,2])) #Save the results for the variable of interest (birth dose)
  }

# Function takes a vector of strings pointing to column names and returns a 
# combined table of results
build.effecttab.rr <- function(x, dsgn, simp=T){
  # Initialize empty data frame
  t <- data.frame(matrix(vector(), length(x), 4,
                         dimnames=list(c(), c("Outcome", "Coef", "L.CI", "U.CI"))),
                  stringsAsFactors=F)
  for(i in 1:length(x)){
    m <- modelrun.rr(x[i], dsgn, simp)
    t[i,] <- m
  }
  t$Coef <- as.numeric(t$Coef)
  t$L.CI <- as.numeric(t$L.CI)
  t$U.CI <- as.numeric(t$U.CI)
  return(t)
}

# Set outcome variables to evaluate for series with and without the hepatitis B series.
outcomevars <- c("utd.dt", "utd.pol", "utd.mmr", "utd.hib", "utd.hepb", 
                 "utd.var", "utd.pcv", "utd.comb7")
outcomevars.nohep <- c("utd.comb7.nohep")

# Unadjusted survey-weighted estimates
tab.svy.hep.rr <- build.effecttab.rr(outcomevars, dsgn.svy, simp=T)
tab.svy.nohep.rr <- build.effecttab.rr(outcomevars.nohep, dsgn.svy, simp=T)

# Simple PS estimates
tab.simple.hep.rr <- build.effecttab.rr(outcomevars, dsgn.cart.ATE, simp=T)
tab.simple.nohep.rr <- build.effecttab.rr(outcomevars.nohep, dsgn.cart.ATE, simp=T)

# Robust estimates
tab.robust.hep.rr <- build.effecttab.rr(outcomevars, dsgn.cart.ATE, simp=F)
tab.robust.nohep.rr <- build.effecttab.rr(outcomevars.nohep, dsgn.cart.ATE, simp=F)


library(writexl)
sheets <- list("survey" = tab.svy.hep.rr, "survey nohep" = tab.svy.nohep.rr,
               "simple" = tab.simple.hep.rr, "simple nohep" = tab.simple.nohep.rr,
               "robust" = tab.robust.hep.rr, "robust nohep" = tab.robust.nohep.rr)
write_xlsx(sheets, "effect_estimates_rr.xlsx")

# Save final estimates
est.final <- list(tab.simple.hep, tab.simple.nohep, tab.robust.hep, tab.robust.nohep)
saveRDS(est.final, "final_estimate_tables_rr.rds")