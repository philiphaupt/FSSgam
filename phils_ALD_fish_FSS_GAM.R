rm(list=ls())
#  run FSS GAM on Phil s data

# librarys----
# detach("package:plyr", unload=TRUE) # this is to ensure that plyr is not loaded prior to dplyr which will cause errors in this script - ignore if youget an error message
library(tidyverse) # probably makes sense to load tidyverse
options(dplyr.width = Inf) #enables head() to display all columns
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
library(doParallel) #this can removed?
library(doSNOW)
library(gamm4)
library(RCurl) #needed to download data from GitHub

# install package----
# devtools::install_github("beckyfisher/FSSgam_package") #run once - Phil: done
library(FSSgam)

#  REad and arrange my data
source("./read_and_arrange_phils_data.R")

# Define the preditcor variable names.
pred.vars <- c("treatment", "depth", "rugosity", "sample_area", "temperature", "we2_mean", "f_effort", "perc_water_column", "has_score", "treatment.partner")

# Check for correlation between numerical variables
round(cor(dat[,pred.vars[2:8]]),2)
# nothing of concern

# Plot of likely transformations - thanks to Anna Cresswell for this loop!#Phil comment: I s this correct when data is arranged long? It will be using multiples for every site.
par(mfrow=c(3,2))
for (i in pred.vars[2:8]) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

# Review of individual predictors - we have to make sure they have an even distribution---
# If the data are squewed to low numbers try sqrt>log or if squewed to high numbers try ^2 of ^3
# Decided that ... needed a sqrt transformation
# Decided ... were not informative variables, and drop.

# I am using these as test case so just continuing for now...

# Check to make sure Response vector has not more than 80% zeros----
# Phil's comment: this will get rid of rarer species like sharks...not sure that this analysis is correct.
unique.vars=unique(as.character(dat$Taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$Taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<0.8){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}

# Run the full subset model selection----
dir.create("./output_phil")
setwd("./output_phil") #Set wd for example outputs - will differ on your computer
resp.vars=unique.vars.use
use.dat=dat #Why is this necessary, as it is defined in the loop below?
factor.vars=c("treatment")# Status as a Factor with two levels
out.all=list()
var.imp=list()

# Loop through the FSS function for each Taxa----
for(i in 1:length(resp.vars)){
  use.dat=dat[which(dat$Taxa==resp.vars[i]),]
  
  Model1=mgcv::gam(response~s(depth,k=3,bs='cr')+ s(treatment, bs="re"), # Phil Haupt:Fails here
                   family=tw(),  data=use.dat)
  
  model.set=generate.model.set(use.dat=use.dat,
                               test.fit=Model1,
                               pred.vars.cont=pred.vars,
                               pred.vars.fact=factor.vars,
                               # linear.vars="",
                               k=3,
                               null.terms="s(treament,bs='re')")
  out.list=fit.model.set(model.set,
                         max.models=600,
                         parallel=T)
  names(out.list)
  
  out.list$failed.models # examine the list of failed models
  mod.table=out.list$mod.data.out  # look at the model selection table
  mod.table=mod.table[order(mod.table$AICc),]
  mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
  out.i=mod.table[which(mod.table$delta.AICc<=3),]
  out.all=c(out.all,list(out.i))
  # var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2
  
  # plot the best models
  for(m in 1:nrow(out.i)){
    best.model.name=as.character(out.i$modname[m])
    
    png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
    if(best.model.name!="null"){
      par(mfrow=c(3,1),mar=c(9,4,3,1))
      best.model=out.list$success.models[[best.model.name]]
      plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
      mtext(side=2,text=resp.vars[i],outer=F)}  
    dev.off()
  }
}
