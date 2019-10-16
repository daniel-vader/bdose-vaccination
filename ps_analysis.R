# Description: Propensity score analysis using boosted CART of the association 
#              between hepatitis B birth dose and up-to-date vaccination 
#              status at 18 months. 
# Author:      Daniel Vader
# Date:        1-10-2019
# Requires:    formats_nis2017.R, data set created in data_step.R

################################################################################
# Load libraries, outside code, and datasets
################################################################################
library(dplyr)
library(twang)
library(survey)
library(tidyr)

source("formats_nis2017.R") # load format function

# load data and apply formats
NIS2017 <- readRDS("NIS2017final.rds") %>% apply_nis_formats() #apply formats

################################################################################
# Check data
################################################################################

varlist <-  names(NIS2017) #create a vector of variable names
for(n in 1:length(varlist)){
  x <- select(NIS2017, varlist[n]) #select the nth variable
  if(is.factor(x)){ #if variable is a factor, output a table
    table(x, useNA = "ifany") %>% print()
  }else{summary(x) %>% print()} # if variable is not a factor, output a summary
  cat("\n")
}
rm(x) #remove x from the workspace


################################################################################
# Generate Propensity scores using boosted cart with survey weights (no strata)
################################################################################

### Generate boosted cart propensity scores ###
# NOTE: ps() fails silently on unused argument errors. eg. if you use an argument 
# that is not listed in the documentation, the function will still run. 

# Choose seed
# sample(1:100000, 1) # Seed generated = 12900
set.seed(12900)

# Average Treatment Effect (ATE)
# Note: can use ps() perm.test.iters = x option to estimate ks p-values using 
#       monte carlo methods (but will be SLOW)
pscore.cart <- ps(birthhep~ AGEGRP + CHILDNM + SEX + CWIC_01 + EDUC1 + FRSTBRN +
                    INCPORAR_I + M_AGEGRP2 + MARITAL2 + MOBIL_I + RACEETHK +
                    RENT_OWN + N_PRVR + PROV_FAC + REGISTRY + VFC_ORDER +
                    INS_BREAK_I + INS_STAT2_I
                  , 
                  sampw = NIS2017$PROVWT_D, shrinkage=0.0005,
                  n.trees = 20000, stop.method = "ks.mean",
                  bag.fraction = .5, 
                  data=NIS2017) 
summary(pscore.cart)
NIS2017$pscore.cart <- pscore.cart$ps[,1] ; #save prob to data set
saveRDS(pscore.cart, "pscoreCART.rds")

### Check relative influence
summary(pscore.cart$gbm.obj)

### Check balance
balance.vars <- c("AGEGRP", "CHILDNM", "FRSTBRN", "SEX", "RACEETHK", "CWIC_01", 
                  "INS_BREAK_I", "INS_STAT2_I", "EDUC1", 
                  "INCPORAR_I", "M_AGEGRP2", "MARITAL2", "MOBIL_I", 
                  "RENT_OWN", "N_PRVR", "PROV_FAC", "REGISTRY", 
                  "VFC_ORDER")

bal.cart.ATE <- dx.wts(pscore.cart, data=NIS2017, estimand = "ATE", treat.var="birthhep",
                       vars= balance.vars)
bal.table(bal.cart.ATE)


NIS2017$cart.ATE <- ifelse(NIS2017$birthhep == 1, 
                   1/NIS2017$pscore.cart * NIS2017$PROVWT_D, 
                   1/(1- NIS2017$pscore.cart) * NIS2017$PROVWT_D)

bal.cart.ATEcomp <- dx.wts(NIS2017$cart.ATE, data=NIS2017, estimand = "ATE", treat.var="birthhep",
                       vars= balance.vars, sampw = NIS2017$PROVWT_D)
bal.table(bal.cart.ATEcomp)

# Output results to data file so we don't have to run PS again ###
saveRDS(NIS2017, file="NIS2017ps.rds")



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
modelrun <- function(x, dsgn, simp=T, verbose=F){
  if(simp){
    f <- as.formula(paste(x, "~ birthhep"))
  } else{
    f <- as.formula(paste(x, "~ birthhep + AGEGRP + SEX + CHILDNM + CWIC_01 + EDUC1 + 
                              FRSTBRN + INCPORAR_I + M_AGEGRP2 + MARITAL2 + 
                              MOBIL_I + RACEETHK + RENT_OWN + N_PRVR + 
                              PROV_FAC + REGISTRY + VFC_ORDER +
                              INS_STAT2_I"))
  }
  m <- svyglm(f, design=dsgn, family=quasibinomial(link=logit))
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
build.effecttab <- function(x, dsgn, simp=T){
  # Initialize empty data frame
  t <- data.frame(matrix(vector(), length(x), 4,
                                     dimnames=list(c(), c("Outcome", "Coef", "L.CI", "U.CI"))),
                              stringsAsFactors=F)
  for(i in 1:length(x)){
    m <- modelrun(x[i], dsgn, simp)
    t[i,] <- m
  }
  t$Coef <- as.numeric(t$Coef)
  t$L.CI <- as.numeric(t$L.CI)
  t$U.CI <- as.numeric(t$U.CI)
  return(t)
}

# Set outcome variables to evaluate for series with and without the hepatitis B series.
outcomevars <- c("utd.dt", "utd.pol", "utd.mmr", "utd.hib", "utd.hepb", 
                 "utd.var", "utd.pcv", "utd.comb3", "utd.comb5", "utd.comb7")
outcomevars.nohep <- c("utd.comb3", "utd.comb5.nohep", "utd.comb7.nohep")

# Unadjusted survey-weighted estimates
tab.svy.hep <- build.effecttab(outcomevars, dsgn.svy, simp=T)
tab.svy.nohep <- build.effecttab(outcomevars.nohep, dsgn.svy, simp=T)

# Simple PS estimates
tab.simple.hep <- build.effecttab(outcomevars, dsgn.cart.ATE, simp=T)
tab.simple.nohep <- build.effecttab(outcomevars.nohep, dsgn.cart.ATE, simp=T)

# Robust estimates
tab.robust.hep <- build.effecttab(outcomevars, dsgn.cart.ATE, simp=F)
tab.robust.nohep <- build.effecttab(outcomevars.nohep, dsgn.cart.ATE, simp=F)

library(writexl)
sheets <- list("survey" = tab.svy.hep, "survey nohep" = tab.svy.nohep,
               "simple" = tab.simple.hep, "simple nohep" = tab.simple.nohep,
               "robust" = tab.robust.hep, "robust nohep" = tab.robust.nohep)
write_xlsx(sheets, "effect_estimates.xlsx")

# Save final estimates
est.final <- list(tab.simple.hep, tab.simple.nohep, tab.robust.hep, tab.robust.nohep)
saveRDS(est.final, "final_estimate_tables.rds")


################################################################################
## Generate summary tables showing balance of covariables across birth dose of 
## survey weighted and ps score weighted data
################################################################################

# Generate survey table based on user-specified design
evaltab <- function(x, dsgn){
  summary(svytable(as.formula(paste("~",x,"+birthhep")), dsgn))
}

# Generate aggregate table summarizing distribution of multiple characteristics 
# by a categorical varible as column percents. 
balancetab <- function(dsgn){
  # Initialize data frame
  bigtab <- data.frame(var=character(),
                       group=character(), 
                       bdose.no=numeric(), 
                       bdose.yes=numeric(),
                       pval=character())
  
  # Build survey-weighted table for each covariable
  for(i in 1:length(balance.vars)){
    if(is.factor(NIS2017[,balance.vars[i]])){
      tsum <- evaltab(balance.vars[i], dsgn) #generate summary of srvytable
      
      t <- as.data.frame(prop.table(tsum$table, margin=2)) # extract table
      t <- spread(t, birthhep, Freq) # convert to wide format
      
      pval <- c(tsum$statistic$p.value) # extract p-val (Rao-Scott)
      length(pval) <- nrow(t) # Match pvalue vector length to number of table rows
      t <- cbind(t, pval) # append p-value to table
      
      varname <- c(balance.vars[i]) # extract name
      length(varname) <- nrow(t) # Match name vector length to number of table rows
      t <- cbind(varname, t)
      
      # Set column names to match main table and remove row names
      colnames(t) <- c("var", "group", "bdose.no", "bdose.yes", "pval")
      rownames(t) <- NULL
      
      # Add covariable table i to main table
      if(i == 1){
        bigtab <- t 
      }else{
        bigtab <- rbind(bigtab, t)
      }
    }
  }
  
  # create empty cells for NAs
  bigtab$var <- as.character(bigtab$var)
  bigtab$var <- ifelse(is.na(bigtab$var), "", bigtab$var)
  bigtab$pval <- ifelse(is.na(bigtab$pval), "", as.character(bigtab$pval)) 
  
  # return main table
  return(bigtab)
}

write.csv(balancetab(dsgn.cart.ATE), "cart-balance.csv")
write.csv(balancetab(dsgn.svy), "svy-balance.csv")


# Compare means of income-poverty ratio across levels of
# hepatitis B birth dose
svyby(~INCPORAR_I, ~birthhep, dsgn.svy, svymean)
svyttest(INCPORAR_I ~ birthhep, dsgn.svy)

svyby(~INCPORAR_I, ~birthhep, dsgn.cart.ATE, svymean)
svyttest(INCPORAR_I ~ birthhep, dsgn.cart.ATE)

