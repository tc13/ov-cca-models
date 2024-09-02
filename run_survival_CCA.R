#Survival analysis script - inferring various things from CCA mortality data

#Load packages
library(cmdstanr)
library(posterior)
library(scales)
library(Hmisc)
library(rethinking)

#Round function
round2 = function(x, n=0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

## Counts of CCA cases from Khon Kaen Province, published in Kamsa-Ard, S. et al. 2011. 
## Trends in liver cancer incidence between 1985 and 2009, Khon Kaen, Thailand: cholangiocarcinoma. Asian Pac J Cancer Prev, 12(9), pp.2209-2213
kk_cases <- data.frame(
  AgeGroup = c("0-14","15-19", "20-24", "25-29", "30-34", "35-39", "40-44", 
               "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+"),
  AgeMid = c(7, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 80),
  CCA_1985_1997 = c(0, 0, 6, 11, 47, 117, 218, 414, 582, 729, 749, 556, 378, 320),
  CCA_1998_2009 = c(0, 3, 3, 9, 41, 100, 238, 395, 748, 948, 1140, 1161, 900, 924)
  )

#Denominators, Khon Kaen Province
#From census data, https://www.citypopulation.de/php/thailand-prov-admin.php?adm2id=40
pop_2000 <- 1733434
#Estimated from the above, assuming a 0.5% annual population increase http://web.nso.go.th/census/poph/popmap/nemap.pdf
pop_1991 <- 1656971

#Age distribution, UN data from Thai census http://data.un.org/
aged <- read.csv("data/who/Thai_demographics_UN.csv", header=T)
aged$Age <- as.character(aged$Age)

#Thai population by age group 1990
age_1990 <- aged[aged$Year=="1990"&aged$Sex=="Both Sexes"&aged$Area=="Total",]
ThaiPop1990 <- c(sum(age_1990$Value[2:4]), age_1990$Value[5:16], sum(age_1990$Value[17:18]))
#Total Thai population 1990
ThaiTotalPop1990 <- age_1990$Value[1]
#Age group as a proportion of total
ThaiPopProp1990 <- ThaiPop1990/ThaiTotalPop1990
#Add 1990 population values to cases data.frame
kk_cases$KKP_TotalPop_1991 <- pop_1991
kk_cases$KKP_PropPop_1991 <- ThaiPopProp1990
kk_cases$KKP_Pop_1991 <- round2(kk_cases$KKP_TotalPop_1991*kk_cases$KKP_PropPop_1991)

#Thai population by age group 2000
age_2000 <- aged[aged$Year=="2000"&aged$Sex=="Both Sexes"&aged$Area=="Total",]
age_2000 <- age_2000[grepl("-",age_2000$Age),]
ThaiPop2000 <- c(sum(age_2000$Value[1:3]), age_2000$Value[4:15], sum(age_2000$Value[16:19]))
#Total Thai population 1990
ThaiTotalPop2000 <- sum(age_2000$Value)
#Age group as a proportion of total
ThaiPopProp2000 <- ThaiPop2000/ThaiTotalPop2000
#Add 1990 population values to cases data.frame
kk_cases$KKP_TotalPop_2002 <- pop_2000
kk_cases$KKP_PropPop_2002 <- ThaiPopProp2000
kk_cases$KKP_Pop_2002 <- round2(kk_cases$KKP_TotalPop_2002*kk_cases$KKP_PropPop_2002)
kk_cases$AgeLower = c(1, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75)
kk_cases$AgeUpper = c(14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 85)

#synthetic data for 11 years, given grouping
set.seed(1234)
kk_case_mat <- sapply(kk_cases$CCA_1998_2009, FUN= function(x) rmultinom(1, x, prob=rep(1/11,11)))
kk_case_vec <- as.vector(t(kk_case_mat))

#Baseline data from Malaysia (code=3236) - WHO mortality data 
#https://www.who.int/data/data-collection-tools/who-mortality-database

ccs <- read.csv("data/who/country_codes", header=T)
pop <- read.csv("data/who/pop")

# ICD code (V10) 
# C22 = Malignant neoplasm of liver and intrahepatic bile ducts
# C24 = Malignant neoplasm of other and unspecified parts of biliary tract
# https://icd.who.int/browse10/2019/en#C22.1

mort1 <- read.csv("data/who/Morticd10_part1", header=T) #up to 2002
mort2 <- read.csv("data/who/Morticd10_part2", header=T) #2003-2007
mort3 <- read.csv("data/who/Morticd10_part3", header=T) #2008-2012

#2007
m_07_ICD <- mort2[(mort2$Cause %in% c("C22", "C24") & mort2$Year==2007 & mort2$Country==3236),]
m_07_pop <- pop[(pop$Country==3236 & pop$Year==2007) , c(8:29)] #population data
#2008
m_08_ICD <- mort3[(mort3$Cause %in% c("C22", "C24") & mort3$Year==2008 & mort3$Country==3236),]
m_08_pop <- pop[(pop$Country==3236 & pop$Year==2008) , c(8:29)] #population data 
#2009
m_09_ICD <- mort3[(mort3$Cause %in% c("C22", "C24") & mort3$Year==2009 & mort3$Country==3236),]
m_09_pop <- pop[(pop$Country==3236 & pop$Year==2009) , c(8:29)] #population data 

#Aggregate ICD codes and sexes
m_cases_07 <- apply(m_07_ICD[c(11:32)], 2, sum) 
m_cases_08 <- apply(m_08_ICD[c(11:32)], 2, sum) 
m_cases_09 <- apply(m_09_ICD[c(11:32)], 2, sum) 

m_pop_07 <- apply(m_07_pop, 2, sum)
m_pop_08 <- apply(m_08_pop, 2, sum)
m_pop_09 <- apply(m_09_pop, 2, sum)

#Set up data frames for Malaysia by year (2007-2009)
malaysia07 <- data.frame(
  AgeGroup = c("0-14","15-19", "20-24", "25-29", "30-34", "35-39", "40-44", 
               "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+"),
  AgeMid = c(7, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 80),
  AgeLower = c(1, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75),
  AgeUpper = c(14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 85),
  cases = c(sum(m_cases_07[1:7]), m_cases_07[8], m_cases_07[9], m_cases_07[10],  m_cases_07[11], 
            m_cases_07[12], m_cases_07[13], m_cases_07[14], m_cases_07[15], m_cases_07[16],
            m_cases_07[17], m_cases_07[18], m_cases_07[19], sum(m_cases_07[20:22])),
  pop = c(sum(m_pop_07[1:7]), m_pop_07[8], m_pop_07[9], m_pop_07[10],  m_pop_07[11], 
          m_pop_07[12], m_pop_07[13], m_pop_07[14], m_pop_07[15], m_pop_07[16],
          m_pop_07[17], m_pop_07[18], m_pop_07[19], sum(m_pop_07[20:22]))
)

malaysia08 <- data.frame(
  AgeGroup = c("0-14","15-19", "20-24", "25-29", "30-34", "35-39", "40-44", 
               "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+"),
  AgeMid = c(7, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 80),
  AgeLower = c(1, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75),
  AgeUpper = c(14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 85),
  cases = c(sum(m_cases_08[1:7]), m_cases_08[8], m_cases_08[9], m_cases_08[10],  m_cases_08[11], 
            m_cases_08[12], m_cases_08[13], m_cases_08[14], m_cases_08[15], m_cases_08[16],
            m_cases_08[17], m_cases_08[18], m_cases_08[19], sum(m_cases_08[20:22])),
  pop = c(sum(m_pop_08[1:7]), m_pop_08[8], m_pop_08[9], m_pop_08[10],  m_pop_08[11], 
          m_pop_08[12], m_pop_08[13], m_pop_08[14], m_pop_08[15], m_pop_08[16],
          m_pop_08[17], m_pop_08[18], m_pop_08[19], sum(m_pop_08[20:22]))
)

malaysia09 <- data.frame(
  AgeGroup = c("0-14","15-19", "20-24", "25-29", "30-34", "35-39", "40-44", 
               "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+"),
  AgeMid = c(7, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 80),
  AgeLower = c(1, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75),
  AgeUpper = c(14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 85),
  cases = c(sum(m_cases_09[1:7]), m_cases_09[8], m_cases_09[9], m_cases_09[10],  m_cases_09[11], 
            m_cases_09[12], m_cases_09[13], m_cases_09[14], m_cases_09[15], m_cases_09[16],
            m_cases_09[17], m_cases_09[18], m_cases_09[19], sum(m_cases_09[20:22])),
  pop = c(sum(m_pop_09[1:7]), m_pop_09[8], m_pop_09[9], m_pop_09[10],  m_pop_09[11], 
          m_pop_09[12], m_pop_09[13], m_pop_09[14], m_pop_09[15], m_pop_09[16],
          m_pop_09[17], m_pop_09[18], m_pop_09[19], sum(m_pop_09[20:22]))
)

###############################
## Prevalence of liver fluke ##
###############################
#Load in R data from model posteriors
load("worm-mod-out.RData")
age <- c(1:85)

M_pre_mat <- matrix(data=0, ncol=length(post$alpha), nrow=length(age))
M_post_mat <- matrix(data=0, ncol=length(post2$alpha), nrow=length(age))
prev_pre_mat <- matrix(data=0, ncol=length(post$alpha), nrow=length(age))
prev_post_mat <- matrix(data=0, ncol=length(post2$alpha), nrow=length(age))

for(a in 1:max(age)){
  for(i in 1:length(post$alpha)){
    M_pre_mat[a,i] = worm_age(a, alpha=post$alpha[i], beta=post$beta[i], mu=post$mu[i])
    M_post_mat[a,i] = worm_age(a, alpha=post2$alpha[i], beta=post2$beta[i], mu=post$mu[i])
    
    prev_pre_mat[a,i] = 1-dnbinom(0, mu=M_pre_mat[a,i], size=post$k_mu[i])
    prev_post_mat[a,i] = 1-dnbinom(0, mu=M_post_mat[a,i], size=post2$k_mu[i])
  }
}

prev_pre <- apply(prev_pre_mat, 1, median)
prev_post <- apply(prev_post_mat, 1, median)

plot(prev_pre~age)
plot(prev_post~age)

###############################
## Examination of cases only ##
###############################
malaysia_all <- data.frame(
                  cases = malaysia07$cases+malaysia08$cases+malaysia09$cases,
                  AgeMid = malaysia07$AgeMid,
                  AgeLower = malaysia07$AgeLower,
                  AgeUpper = malaysia07$AgeUpper
)

par(mfrow=c(1,3))
barplot(kk_cases$CCA_1985_1997~kk_cases$AgeMid, ylim=c(0,1100))
barplot(kk_cases$CCA_1998_2009~kk_cases$AgeMid, ylim=c(0, 1100))
barplot(malaysia_all$cases~malaysia_all$AgeMid, ylim=c(0, 1100))

###############
## Figure 3A ##
###############

#pdf("CA-Cases-TH.pdf", height=6, width=7.6)
thai_cases <- kk_cases$CCA_1985_1997+kk_cases$CCA_1998_2009
par(mfrow=c(1,1), mar=c(6.5,6,2,1))
barplot(thai_cases~kk_cases$AgeMid, ylim=c(0,2000), col=alpha("red2", alpha=0.5),
        xlab="", ylab="", las=2, cex=1.3,
        names.arg=kk_cases$AgeGroup,
        cex.axis=1.6, cex.lab=1.6, axes=F)
axis(2, las=2, cex.axis=1.4, lwd=1.2)
title(ylab="Cases of cholangiocarcinoma", cex.lab=1.65, line=4.5)
title(xlab="Age group (years)", cex.lab=1.65, line=5.2)
#dev.off()

#################################################
## Stan model for latent and induction periods ##
#################################################

d_cases <- list(
    N = sum(thai_cases),
    age_groups = length(thai_cases),
    cases = thai_cases,
    age_lower = kk_cases$AgeLower,
    age_upper = kk_cases$AgeUpper,
    #median values from 'pre-intervention' model (1980s)
    foi_alpha = 0.961,
    foi_beta = 0.0303,
    k = 0.348,
    max_age = max(kk_cases$AgeUpper)
)
    
stan_file_time <- file.path("stan", "time-to-cca.stan")
mod_time <- cmdstan_model(stan_file_time, compile = FALSE)
mod_time$compile(force_recompile = TRUE, cpp_options = list(stan_threads = TRUE))

fit_time <- mod_time$sample(
  data = d_cases,
  seed = 123,
  chains =  4,
  parallel_chains = 4,
  threads_per_chain = 1,
  refresh = 10, # print update every 50 iters
  iter_warmup = 1200,
  iter_sampling = 500,
  init = 5
)

fit_time$summary(
  variables = c("induction_mean", "induction_shape", 
                "latent_mean", "latent_shape")
)

tpars <- fit_time$draws(
  variables=c("induction_mean", "induction_shape", 
              "latent_mean", "latent_shape"),
                        format="data.frame")

induction_shape <- median(tpars$induction_shape)
induction_rate <- median(tpars$induction_shape/tpars$induction_mean)

latent_shape <- median(tpars$latent_shape)
latent_rate <- median(tpars$latent_shape/tpars$latent_mean)

induction_pred <- dgamma(seq(0,55, 0.1), shape=induction_shape, rate=induction_rate)
latent_pred <- dgamma(seq(0,55, 0.1), shape=latent_shape, rate=latent_rate)

induction_CrI90 <- quantile(tpars$induction_mean, probs=c(0.05, 0.5, 0.95))
latent_CrI90 <- quantile(tpars$latent_mean, probs=c(0.05, 0.5, 0.95))

#########################################
## Figure 3B - Plot both distributions ##
#########################################

#pdf("induction-latent-posterior.pdf", height=7, width=7)

xseq <- seq(0,55,0.1)
par(mfrow=c(2,1), xpd=F, mar=c(6,6,2,1))

plot(induction_pred~xseq, type="l",
      lwd=1.2, axes=F,
      ylab="", yaxs="i", ylim=c(0,0.055),
      xlab="Time (years)", cex.lab=1.4)
#abline(v=median(tpars$induction_mean), lty=2, lwd=1.4)

polygon(x=c(induction_CrI90[1], 
            xseq[which(xseq==round(induction_CrI90[1])):which(xseq==round(induction_CrI90[3]))], induction_CrI90[3]), 
        y=c(0, induction_pred[which(xseq==round(induction_CrI90[1])):which(xseq==round(induction_CrI90[3]))],0),
        col=alpha("orange1",0.7), border=alpha("orange1",0.7))

lines(x=c(median(tpars$induction_mean),median(tpars$induction_mean)), 
      y=c(0, induction_pred[xseq==round(median(tpars$induction_mean))]-0.001), 
      lty=2, lwd=1.8)

axis(1, lwd=1.4, cex.axis=1.3)
axis(2, lwd=1.4, cex.axis=1.3, las=2)

curve(dgamma(x, shape=induction_shape, rate=induction_rate), 
      from=0, to=55, lwd=1.5, add=T)

title(ylab="Posterior prob. density", cex.lab=1.4, line=4)
title(main="Time from parasite exposure to driver mutation", cex.main=1.3)

curve(dgamma(x, shape=latent_shape, rate=latent_rate), 
            from=0, to=55, lwd=1.2, axes=F,
            ylab="", yaxs="i", ylim=c(0,0.055),
            xlab="Time (years)", cex.lab=1.4)

polygon(x=c(latent_CrI90[1], 
            xseq[which(xseq==round(latent_CrI90[1])):which(xseq==round(latent_CrI90[3]))], latent_CrI90[3]), 
        y=c(0, latent_pred[which(xseq==round(latent_CrI90[1])):which(xseq==round(latent_CrI90[3]))],0),
        col=alpha("orange1",0.7), border=alpha("orange1",0.7))

lines(x=c(median(tpars$latent_mean),median(tpars$latent_mean)), 
      y=c(0, latent_pred[xseq==round(median(tpars$latent_mean))]), 
      lty=2, lwd=1.8)

curve(dgamma(x, shape=latent_shape, rate=latent_rate), 
      from=0, to=55, lwd=1.5, add=T)
axis(1, lwd=1.4, cex.axis=1.3)
axis(2, lwd=1.4, cex.axis=1.3, las=2)

title(ylab="Posterior prob. density", cex.lab=1.4, line=4)
title(main="Time from driver mutation to cancer", cex.main=1.3)

#dev.off()

##############################
## Data list for stan model ##
##############################

baseline <- rbind(malaysia07, malaysia08, malaysia09)

cca_stan <- list(
  N_baseline = nrow(baseline),
  N_fluke = length(kk_case_vec),
  N_age_class = length(unique(baseline$AgeGroup)),
  #age_class_b = match(baseline$AgeGroup, unique(baseline$AgeGroup)),
  #age_class_f = match(kk_cases$AgeMid, unique(kk_cases$AgeMid)),
  
  tested_baseline = baseline$pop,
  positive_baseline = round2(baseline$cases/2), #assumption that 50% of liver cancers are CCA
  age_mid_baseline = baseline$AgeMid,
  age_lower_baseline = baseline$AgeLower,
  age_upper_baseline = baseline$AgeUpper,
  
  tested_fluke = rep(kk_cases$KKP_Pop_2002, 11),
  positive_fluke = kk_case_vec,
  age_mid_fluke = rep(kk_cases$AgeMid, 11),
  age_lower_fluke = rep(kk_cases$AgeLower, 11),
  age_upper_fluke = rep(kk_cases$AgeUpper, 11),
  
  shape_mut = 3, #shape parameter for neg_binomial dist (analogous to gamma)
  shape_can = 10, #shape parameter for neg_binomial dist (analogous to gamma)
  
  foi_alpha = 0.0129, #value from model run - posterior median
  foi_beta = 0.0222, #value from model run - posterior median
  k = 0.109, #value from model run - posterior median
  max_age = max(baseline$AgeUpper)
)

dist <- matrix(ncol=cca_stan$max_age, nrow=cca_stan$max_age, data = 0)
for(i in 1:cca_stan$max_age){
  for(j in 1:cca_stan$max_age)
  {dist[i,j] <- abs(i-j)}
}

cca_stan$Dmat <- dist #add distance matrix for correlations

########################
## compile stan model ##
########################

stan_file_cca <- file.path("stan","prob-cca-discrete-corr-matrix.stan")
mod_cca <- cmdstan_model(stan_file_cca, compile = FALSE)
mod_cca$compile(force_recompile = TRUE, cpp_options = list(stan_threads = TRUE))

fit_cca <- mod_cca$sample(
  data = cca_stan,
  seed = 1234,
  chains =  4,
  parallel_chains = 4,
  threads_per_chain = 1,
  refresh = 10, # print update every 50 iters
  iter_warmup = 1700,
  iter_sampling = 300
)

fit_cca$summary(variables = c("beta_fluke", "etasq", "rhosq"),
                "mean", "median","sd", "rhat", "ess_bulk", "ess_tail")
                

# Posterior extraction
lpb <- fit_cca$draws(variables=c("prob_baseline"), format="matrix")
lpf <- fit_cca$draws(variables=c("prob_fluke"), format="matrix")

pb <- apply(lpb, 2, median)
pb_90CrI <- apply(lpb, 2, FUN = function(x) quantile(x, probs=c(0.05, 0.95)))
pf <- apply(lpf, 2, median)
pf_90CrI <- apply(lpf, 2, FUN = function(x) quantile(x, probs=c(0.05, 0.95)))

par(mfrow=c(1,2))
plot(pb)
plot(pf)

pred_p_b <- pred_p_f <- c()

for(i in 1:length(kk_cases$AgeMid)){
  pred_p_b[i] =  1-prod(1-pb[kk_cases$AgeLower[i]:kk_cases$AgeUpper[i]])
  
  pred_p_f[i] =  1-prod(1-pf[kk_cases$AgeLower[i]:kk_cases$AgeUpper[i]])
}

###############
## Figure 3C ##
###############
malaysia_age_offset <- malaysia08$AgeMid-0.5
thai_age_offset <- kk_cases$AgeMid+0.5
kk_cases$CCA_1998_2009_100K <- ((kk_cases$CCA_1998_2009/11)/kk_cases$KKP_Pop_2002)*100000

#pdf("CCA-incidence-TH-MAL.pdf",  height=6, width=7.5)
par(mfrow=c(1,1), mar=c(6,6,1,1))
plot((malaysia08$cases/2)/malaysia08$pop*100000~malaysia_age_offset, 
     pch=16, cex=2.2, xlim=c(0,80), ylim=c(0,300), axes=F, col=alpha("navy",0.55),
     xlab="Age (years)", cex.lab=1.6, ylab="")
points((malaysia08$cases/2)/malaysia08$pop*100000~malaysia_age_offset,
       cex=2.2)
points(kk_cases$CCA_1998_2009_100K~thai_age_offset,
       pch=16, cex=2.2, col=alpha("red2", 0.55))
points(kk_cases$CCA_1998_2009_100K~thai_age_offset,
      cex=2)
axis(1, at = c(0,10,20,30,40,50,60,70,80), cex.axis=1.35)
axis(2, las=2, cex.axis=1.35)
title(ylab="Annual incidence of CCA per 100,000",
      cex.lab=1.6, line=4)

legend(y=c(250,300), x=c(2,22), cex=1.3, y.intersp=1.4,
       #legend=c("N.E. Thailand (1998-2009)", "Malaysia (2007-2009)"),
       legend=c("N.E. Thailand", "Malaysia"),
       pch=16, pt.cex=2.4, col=c(alpha("red2", 0.6), alpha("navy",0.6)), bty="n")

legend(y=c(250,300), x=c(2,22), cex=1.3,
       legend=c("", ""), y.intersp=1.4,
       pch=1, pt.cex=2.4, col=c("black","black"), bty="n")
#dev.off()

points(pred_p_f*100000~kk_cases$AgeMid, col="orange2", pch=16)
points(pred_p_b*100000~kk_cases$AgeMid, col="blue1", pch=16)

#Probability from model output
cpb <- cpf <- c()
cpb_int <- cpf_int <- matrix(nrow=2, ncol=85, data=0)
for(a in 1:85){
  cpb[a] = 1-prod(1-pb[1:a])
  
  cpb_int[1,a] <- 1-prod(1-pb_90CrI[1,1:a])
  cpb_int[2,a] <- 1-prod(1-pb_90CrI[2,1:a])
  
  cpf[a] = 1-prod(1-pf[1:a])
  
  cpf_int[1,a] <- 1-prod(1-pf_90CrI[1,1:a])
  cpf_int[2,a] <- 1-prod(1-pf_90CrI[2,1:a])
}

###############
## Figure 3D ##
###############

#pdf("cumulative-prob-cca.pdf", width=7.5, height=6)
par(mfrow=c(1,1), mar=c(5.5, 5.5, 2,2))
plot(cpf*100~seq(1,85,1), type="l", axes=F, xlab="Age (years)", 
     ylab="", lwd=2, ylim=c(0,1.4),
     cex.lab=1.4)
lines(cpb*100~seq(1,85,1), type="l", lwd=2)
axis(1, lwd=1.4, cex.axis=1.3)
axis(2, lwd=1.4, cex.axis=1.3, las=2)
title(ylab="Cumulative probability of CCA (%)", cex.lab=1.4, line=4)
shade((cpb_int*100), c(1:85), col=alpha("navy", 0.35))
shade((cpf_int*100), c(1:85), col=alpha("red2", 0.35))

legend(y=c(1.25,1.45), x=c(2,22), cex=1.3, y.intersp=1.4,
       legend=c("", ""), lwd=1,
       lty=1, col="black", bty="n")

legend(y=c(1.25,1.45), x=c(2,22), cex=1.3, y.intersp=1.4,
       legend=c("N.E. Thailand", "Malaysia"), lwd=12,
       lty=1, col=c(alpha("red2", 0.5), alpha("navy", 0.5)), bty="n")                  

#dev.off()

## Calculation for lifetime probability given infection
((cpf[75]-cpb[75])*0.9/0.88 + cpb[75]) *100
(((cpf_int[1,75]-cpb_int[1,75])*0.9/0.88) + cpb_int[1,75]) *100
(((cpf_int[2,75]-cpb_int[2,75])*0.9/0.88) + cpb_int[2,75]) *100

((cpf[85]-cpb[85])*0.9/0.88 + cpb[75]) *100
(((cpf_int[1,85]-cpb_int[1,85])*0.9/0.88) + cpb_int[1,85]) *100
(((cpf_int[2,85]-cpb_int[2,85])*0.9/0.88) + cpb_int[2,85]) *100
