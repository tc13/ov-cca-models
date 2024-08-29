## Script to run models for worm burden by age
## Generates figures 2A-D for manuscript
## Thomas Crellen, August 2024
## thomas.crellen@glasgow.ac.uk

#Load packages
require("cmdstanr")
require("posterior")
require("scales") 
require("Hmisc")
require("rethinking")
require("ggplot2")
require("bayesplot")

#Functions relating to worm burden
worm_func <- function(x, lambda, mu) (lambda/mu)*(1-exp(-mu*x))

worm_age <- function(x, alpha=1.34, beta=0.036, mu=0.05){ 
  e1 = alpha*exp(-mu*x)/((beta-mu)^2)
  e2 = (1+exp(x*(mu-beta))*(x*(mu-beta)-1))
  return(e1*e2)
}

foi <- function(x, alpha=1.344, beta=0.0364) alpha*x*exp(-x*beta)

foi_def_integral <- function(x, alpha=1.344, beta=0.0364){
  y = -(alpha*exp(-beta*x)*(beta*x+1))/beta^2
  y0 = -(alpha*exp(-beta*0)*(beta*0+1))/beta^2
  return(y-y0)
}

epg_to_worm <- function(x, L1=78, gamma=0.76) (x/L1)^(1/gamma) 

#Round function
round2 = function(x, n=0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

#Read in synthetic worm burden data for model
worm_data <- read.csv("ov-data/worms-age-data.csv", header=T)
eggs_data <- read.csv("ov-data/eggs-age-data.csv", header=T)

#These data frame contain studies S1-S5 as shown in Table 1 of manuscript
#Get indices for studies

#Study 1 = Autopsy
aut.idx <- which(worm_data$study_type==1)
#Study 2 = Expulsion studies (Ramsay and Elkins) combined
expul.idx <- which(worm_data$study_type==2)
#Study 3 = Upatham egg counts
upatham.idx <- which(eggs_data$study_id==3)
#Study 4 = Rural egg counts (Kurathong)
urban.idx <- which(eggs_data$study_id==4)
#Study 5 = Urban egg counts (Kurathong)
rural.idx <- which(eggs_data$study_id==5)

####################################
## plot underlying data for S1-S5 ##
####################################
par(mfrow=c(1,2),
mar=c(5,5,3,3))
#Worm burdens (autopsy + expulsion)
plot(log(worm_data$worms[aut.idx]+1)~worm_data$age[aut.idx], pch=16, 
     xlim=c(0,80), axes=F, ylim=c(0, log(3000)), cex=1.1, cex.lab=1.4,
     xlab="Host age (years)", ylab="Worm burden", col=alpha("skyblue3", alpha=0.65))
axis(1, lwd=1.2, cex.axis=1.2)
axis(2, las=2, at=log(c(0,1,10, 100, 1000, 3000)+1), 
     labels = c(0,1,10, 100,1000,3000), lwd=1.2, cex.axis=1.2)
title("Liver fluke worm burdens", cex.main=1.2)
points(log(worm_data$worms[expul.idx]+1)~worm_data$age[expul.idx], pch=16, 
       col=alpha("tomato3", alpha=0.65), cex=1.1)
legend("topright", bty="n", title="Study type", legend=c("Autopsy", "Expulsion"),
       pch=16, pt.cex=1.5, col=c("skyblue3", "tomato3"), title.cex = 1.2)

#Egg counts
plot(log(eggs_data$eggs[upatham.idx]+1)~eggs_data$age[upatham.idx], pch=16,
     xlab="Host age (years)", ylab="", axes=F,
     col=alpha("seagreen4", alpha=0.65), cex=1.1, cex.lab=1.4)
axis(1, lwd=1.2, cex.axis=1.2)
axis(2, las=2, at=log(c(0,10,100,1000,10000,50000)+1), 
     labels=c(c(0,10,100,1000,10000,50000)), lwd=1.2, cex.axis=1.2)
title("Liver fluke faecal egg counts", cex.main=1.2)
title(ylab="Egg count", cex.lab=1.4, line=3.8)
points(log(eggs_data$eggs[rural.idx]+1)~eggs_data$age[rural.idx], pch=16,
       col=alpha("seagreen4", alpha=0.65))
points(log(eggs_data$eggs[urban.idx]+1)~eggs_data$age[urban.idx], pch=16,
       col=alpha("seagreen4", alpha=0.65))
#dev.off()

##Study 6 - Sornmani et al. A Pilot Project for Controlling 
## O. viverrini Infection in Nong Wai, Northeast Thailand, 
## by Applying Praziquantel and Other Measures

sorn_pilot <- data.frame(age_mid = c(2, 9.5, 17, 24.5, 34.5, 49.5, 65),
                         tested = c(88, 250, 83, 192, 193, 136, 59),
                         positive = c(3, 118, 57, 105, 106, 91, 45),
                         study_id = rep(6, 7),
                         age_class = 1:7)

sorn_control <- data.frame(age_mid = c(2, 9.5, 17, 24.5, 34.5, 49.5, 65),
                           tested = c(42, 104, 37, 49, 49, 79, 30),
                           positive = c(7, 30, 24, 31, 32, 67, 20),
                           study_id = rep(6, 7),
                           age_class = 1:7)

#Binomial distribution
age_class_mid <- round2(aggregate(c(worm_data$age, eggs_data$age), 
                                  by=list(c(worm_data$age_class, eggs_data$age_class)),
                                  FUN=median))

sornmani <- data.frame(age_mid = age_class_mid$x[1:8],
                       tested = c((88+42), (250+104), (83+37), 285, 98, 136, 79, (59+30)),
                       positive = c((3+7), (118+30), (57+24), 211, 63, 91, 67, (45+20)),
                       study_id = rep(6, 8))

sornmani_binom <- binconf(x=sornmani$positive, n=sornmani$tested)

## Visualise datasets for prevalence
#worms
worm_binary <- ifelse(worm_data$worms==0,0,1)
worm_pos <- aggregate(worm_binary, by=list(worm_data$age_class), sum)
worm_tested <- aggregate(worm_binary, by=list(worm_data$age_class), length)
worm_binom <- binconf(x=worm_pos$x, n=worm_tested$x)

#eggs
egg_binary <- ifelse(eggs_data$eggs==0,0,1)
egg_pos <- aggregate(egg_binary, by=list(eggs_data$age_class), sum)
egg_tested <- aggregate(egg_binary, by=list(eggs_data$age_class), length)
egg_binom <- binconf(x=egg_pos$x, n=egg_tested$x)

#############################
## plot prevalence and CIs ##
#############################
par(mfrow=c(1,3),
    mar=c(5,5,3,3))

plot(worm_binom[,1]~age_class_mid$x, pch=16, 
     ylim=c(0,1),xlim=c(0,80), xlab="Age (years)",
     ylab="Prevalence", cex=1.2, axes=F, cex.lab=1.4)
arrows(x0=age_class_mid$x, x1= age_class_mid$x,
       y1 = worm_binom[,2], y0 = worm_binom[,1],
       angle=90, length = 0.05, lwd=1.2)
arrows(x0=age_class_mid$x, x1= age_class_mid$x,
       y1 = worm_binom[,3], y0 = worm_binom[,1],
       angle=90, length = 0.05, lwd=1.2)
axis(1, cex.axis=1.3, lwd=1.2)
axis(2, las=2, cex.axis=1.3, lwd=1.2)
title("Worm burden studies", cex.main=1.2)

plot(egg_binom[,1]~age_class_mid$x, pch=16, 
     ylim=c(0,1),xlim=c(0,80), xlab="Age (years)",
     ylab="Prevalence", cex=1.2, axes=F, cex.lab=1.4)
arrows(x0=age_class_mid$x, x1= age_class_mid$x,
       y1 = egg_binom[,2], y0 = egg_binom[,1],
       angle=90, length = 0.05, lwd=1.2)
arrows(x0=age_class_mid$x, x1= age_class_mid$x,
       y1 = egg_binom[,3], y0 = egg_binom[,1],
       angle=90, length = 0.05, lwd=1.2)
axis(1, cex.axis=1.3, lwd=1.2)
axis(2, las=2, cex.axis=1.3, lwd=1.2)
title("Egg count studies", cex.main=1.2)

plot(sornmani_binom[,1]~sornmani$age_mid, pch=16, 
     ylim=c(0,1),xlim=c(0,80), xlab="Age (years)",
     ylab="Prevalence", cex=1.2, axes=F, cex.lab=1.4)
axis(1, cex.axis=1.3, lwd=1.2)
axis(2, las=2, cex.axis=1.3, lwd=1.2)
arrows(x0=sornmani$age_mid, x1= sornmani$age_mid,
       y1 = sornmani_binom[,2], y0 = sornmani_binom[,1],
       angle=90, length = 0.05, lwd=1.2)
arrows(x0=sornmani$age_mid, x1= sornmani$age_mid,
       y1 = sornmani_binom[,3], y0 = sornmani_binom[,1],
       angle=90, length = 0.05, lwd=1.2)
title("Sornmani survey", cex.main=1.2)

#Collect data for model into a list
d_1 <- list(
  N_studies = 5,
  N_worm_studies = 2,
  N_egg_studies = 3,
  N_prev_studies = 1,
  N_worms = nrow(worm_data),
  N_eggs = nrow(eggs_data),
  study_id_worm = worm_data$study_id,
  study_id_egg = eggs_data$study_id,
  worms = worm_data$worms,
  eggs = eggs_data$eggs,
  N_classes = nrow(age_class_mid),
  tested = sornmani$tested,
  positive = sornmani$positive,
  age_class_mid = age_class_mid$x,
  age_class_worm = worm_data$age_class,
  age_class_egg = eggs_data$age_class,
  max_worm = 2000
)

## First model - single age-variable FOI across four studies from the 1980s
## Curve fitted to mean worm burden by age categories
## Worm aggregation (k) estimated independently for each age category (with hierarchical prior)
stan_file_1 <- file.path("dynamic-worm-mod.stan")
mod1 <- cmdstan_model(stan_file_1, compile = FALSE)
mod1$compile(force_recompile = TRUE, cpp_options = list(stan_threads = TRUE))

## Fit model using stan, with 4 chains on 4 cores
## Note this takes around 12 hours to run
## Alternatively load("worm-mod-out.RData") to access posteriors
fit1 <- mod1$sample(
  data = d_1,
  seed = 999,
  chains =  4,
  parallel_chains = 4,
  threads_per_chain = 1,
  refresh = 10, # print update every 50 iters
  iter_warmup = 1700,
  iter_sampling = 300
)

# Examine model output
fit1$summary(variables = c("alpha", "beta", "mu", "k_mu", "k_var", "L1", "gamma"),
             "mean", "median","sd", "rhat", "ess_bulk", "ess_tail")


# Posterior extraction
post <- fit1$draws(variables = c("alpha", "beta", "mu", "L1", "gamma", "k_mu"), 
                   format = "df")

loglik <- fit1$draws(variables = c("log_lik"), format = "matrix")
LL <- apply(loglik, 1, mean)

# Parameter correlations
with(post, plot(alpha~beta))
with(post, plot(alpha~(log(2)/mu)))
post$lifespan <- log(2)/post$mu
with(post, plot(lifespan, beta, xlab="Worm life expectancy (years)", 
                xlim=c(0,50), pch=16, col=alpha("grey20", 0.6)))
with(post, plot(L1~gamma))

p2 <- mcmc_scatter(post, pars = c("lifespan", "beta"), 
                              size = 3.5, alpha = 0.25)
p2 + stat_density_2d(color = "black", size = 0.4) + xlim(0,50) + 
  xlab("Worm mean lifespan") + theme_bw()

# Worm life expectancy
summary(log(2)/post$mu)
quantile(log(2)/post$mu, probs=c(0.05,0.5,0.95))
hdpi_90 <- HPDI(log(2)/post$mu, prob=0.90)

###############################################################
## Figure 2D - posterior probability of worm life expectancy ##
###############################################################

#pdf("worm-life-expectancy.pdf", height=6, width=7.5)
par(mar=c(5,6.4,3,2), xpd=0, mfrow=c(1,1))
dens <- density(log(2)/post$mu, adjust=1.8, kernel="gaussian")
plot(dens,
     xlim=c(0, 40),
     main=expression(paste(bold("Lifespan of "), bolditalic("O. viverrini"), bold(" worms in human infections"))), 
     xlab="Parasite average lifespan (years)",
     cex.lab=1.6,
     ylab="", yaxs="i", xaxs="i", axes=F,
     lwd=2, ylim=c(0,0.082),
     cex.main=1.6)

polygon(x=c(hdpi_90[1], dens$x[dens$x>=hdpi_90[1]&dens$x<=hdpi_90[2]], hdpi_90[2]), 
        y=c(0, dens$y[dens$x>=hdpi_90[1]&dens$x<=hdpi_90[2]],0),
        col=alpha("orange2",0.7), border=alpha("orange2",0.7))

lines(x=c(median(log(2)/post$mu),median(log(2)/post$mu)), 
      y=c(0, dens$y[which.min(abs(dens$x - median(log(2)/post$mu)))]+0.0002), 
          lty=2, lwd=2.8)
lines(dens, lwd=2)
axis(1, lwd=2, cex.axis=1.5)
axis(2, las=2, lwd=2, cex.axis=1.5)
title(ylab="Posterior probability density", cex.lab=1.6, line = 4.5)
#dev.off()

# Relationship between M and k
k.post <- fit1$draws(variables = c("k"), format = "matrix")
M.post <- fit1$draws(variables = c("M_age"), format = "matrix")
k.mean <- apply(k.post, 2, mean)
M.mean <- apply(M.post, 2, mean)
plot(k.mean~M.mean, ylim=c(0,max(k.mean)), xlim=c(0,max(M.mean)), pch=16, cex=1.2)

#############################
## Plot data and model fits
#############################
worms_pr_worm <- fit1$draws(variables = c("wb_infered_worms"), format = "matrix")
worms_pr_egg <- fit1$draws(variables = c("wb_infered_eggs"), format = "matrix")
wiw <- apply(worms_pr_worm, 2, median)
wie <- apply(worms_pr_egg, 2, median)

#Plot individual worm burdens by age
aut.idx <- which(worm_data$study_id==1)
ram.idx <- which(worm_data$study_id==2)
elk.idx <- which(worm_data$study_id==3)

#Simulate worm values
set.seed(1234)
x_age <- seq(0,80,0.25)
M_mu <- matrix(data=0, nrow=length(x_age), ncol=length(post$alpha))
M_sim <- matrix(data=0, nrow=length(x_age), ncol=length(post$alpha))
p_sim <- matrix(data=0, nrow=length(x_age), ncol=length(post$alpha))
for(a in 1:length(x_age)){
  for(j in 1:length(post$alpha)){
    M_mu[a,j] = worm_age(x_age[a], alpha=post$alpha[j], beta=post$beta[j], mu=post$mu[j])
    M_sim[a,j] = median(rnbinom(n=2000, mu = M_mu[a,j], size=post$k_mu[j]))
    p_sim[a,j] = 1-dnbinom(0, mu=M_mu[a,j], size=post$k_mu[j])
  }
}

M_mean <- apply(M_mu, 1, mean)
M_median <- apply(M_sim, 1, median)
M_interval <- apply(M_sim, 1, FUN=function(x) quantile(x, probs=c(0.05, 0.95)))
M_mu_interval <- apply(M_mu, 1, FUN=function(x) quantile(x, probs=c(0.05, 0.95)))
p_mean <- apply(p_sim, 1, mean)

####################################
## Figure 2A - worm burden by age ##
####################################

#pdf("pre-intervention-data-curve.pdf", height=7.6, width=9.5)
par(mar=c(5,6.5,1.8,1.8), xpd=0, mfrow=c(1,1))
plot(log(wiw[aut.idx]+1) ~ worm_data$age[aut.idx],
     xlim=c(0,80),
     axes=F, ylim=c(0.01,8.3),pch=16, cex=1.2,cex.lab=1.8,
     xlab="Age (years)", ylab="",col=alpha("skyblue3", 0.75))
points(log(wiw[ram.idx]+1) ~ worm_data$age[ram.idx],
       cex=1.2,cex.lab=1.2, col=alpha("tomato3",0.75), pch=16)
points(log(wiw[elk.idx]+1) ~ worm_data$age[elk.idx],
       cex=1.2,cex.lab=1.2, col=alpha("tomato3",0.75), pch=16)
points(log(wie+1) ~ eggs_data$age,
       cex=1.2,cex.lab=1.2, col=alpha("seagreen4",0.6), pch=16)
title(ylab = "Worm burden (inferred)", line=4, cex.lab=1.8)
axis(1, at=c(0, 10, 20, 30, 40, 50, 60, 70, 80), lwd=1.4,
     labels=c(0, 10, 20, 30, 40, 50, 60, 70, 80), cex.axis=1.65)
axis(2, las=2, at=log(c(0, 1, 10, 100, 1000, 3000)+1),
     labels=c(0, 1, 10, 100, 1000, 3000), cex.axis=1.65, lwd=1.4)

shade(log(M_interval+1), x_age, col=alpha("grey50", 0.4))
lines(log(M_median+1)~x_age, lwd=2.8, col="black")

legend(title="Study type", legend = c("Autopsy", "Worm expulsion", "Faecal egg diagnostics"), 
       x=c(48, 58), y=log(c(5000,6000)), pch=16, cex=1.4, bty="n", pt.cex=2,
       col=c("skyblue3", "tomato3", "seagreen4"))
#dev.off()

######################################
## Second model with later surveys
######################################

stan_file_2 <- file.path("dynamic-worm-age-binomial.stan")
mod2 <- cmdstan_model(stan_file_2, compile = FALSE)
mod2$compile(force_recompile = TRUE, cpp_options = list(stan_threads = TRUE))

## Prev age classes
prev_mid_age <- c(5, 15, 25, 35, 45, 55, 65)

jong_data <- read.csv("ov-data/jong-eggs-age-data.csv", header=T)

# Jongsuksuntigul & Imsomboon 1997 - study conducted 1994
#Egg count data also available
jongsuksuntigul <- data.frame(age_class = 1:7,
                              tested = c(262, 271, 294, 363, 294, 201, 193),
                              positive = c(8, 38, 67, 75, 69, 50, 48),
                              study=1)


#Kaewpitoon et al. 2012 - https://koreascience.kr/article/JAKO201205061576159.pdf
# Nakhon Ratchasima
kaewpitoon <- data.frame(age_class = 1:7,
                         tested = c(62, 47, 51, 155, 451, 218, 184),
                         positive = c(0, 1, 1, 4, 12, 6, 5),
                         study=2)

#Thaewnongiew et al. 2014 - https://koreascience.kr/article/JAKO201429765166986.pdf
#Study conducted 2013
thaewnongiew <- data.frame(age_class = 1:7,
                           tested = c(0, 455, 500, 665, 881, 750, 665),
                           positive = c(0, 65, 94, 150, 231, 241, 196),
                           study=3)

#Kaewpitoon et al. 2016 - https://journal.waocp.org/article_16378_dce94accc908820def85cf4ae60dad6c.pdf
# Surveys in 2016 at Kaeng Sanam Nang district, Nakhon Ratchasima province, 
#Khon Sawan district of Chaiyaphum province, 
#and Waeng Noi district of Khon Kaen province
kaewpitoon_border <- data.frame(age_class = 1:7,
                                tested = c(0, 0, 0, 80, 331, 357, 210),
                                positive = c(0, 0, 0, 2, 4, 7, 4),
                                study=4)

#Laoraksawong et al. 2018 - https://link.springer.com/content/pdf/10.1186/s12889-018-5871-1.pdf
# Study conducted 2017 in Mueang Khon Kaen district
laoraksawong <- data.frame(age_class = 1:7,
                           tested = c(0, 0, 17, 14, 73, 94, 189),
                           positive = c(0, 0, 2, 3, 12, 17, 45),
                           study=5)

#get data as array
tested_array <- t(array(data=c(jongsuksuntigul$tested, 
                               kaewpitoon$tested, 
                               thaewnongiew$tested, 
                               kaewpitoon_border$tested, 
                               laoraksawong$tested), 
                          dim=c(7,5)))

positive_array <- t(array(data=c(jongsuksuntigul$positive, 
                                 kaewpitoon$positive, 
                                 thaewnongiew$positive, 
                                 kaewpitoon_border$positive, 
                                 laoraksawong$positive), 
                          dim=c(7,5)))

# Collect data for binomial model
d_2 <- list(
 N_studies = 5,
 N_classes = 7,
 N_eggs = nrow(jong_data),
 eggs = jong_data$eggs,
 age_class_egg = jong_data$age_class,
 age_class_mid = prev_mid_age,
 tested = apply(tested_array, 2, sum),
 positive = apply(positive_array, 2, sum),
 mu = quantile(post$mu, probs=(0.99)),
 L1 = mean(post$L1),
 gamma = mean(post$gamma),
 max_worm = 2000
)

#Quick plot
prev <- d_2$positive/d_2$tested
plot(prev~prev_mid_age, pch=16, cex=1.2, axes=F, xlim=c(0,70), ylim=c(0,0.25))
axis(1)
axis(2, las=2)

#Fit model 2
fit2 <- mod2$sample(
  data = d_2,
  seed = 999,
  chains =  4,
  parallel_chains = 4,
  threads_per_chain = 1,
  refresh = 10, # print update every 50 iters
  iter_warmup = 1700,
  iter_sampling = 300
)

# Examine model output
fit2$summary(variables = c("alpha" ,"beta", "k_mu", "k_var"),
             "mean", "median","sd", "rhat", "ess_bulk", "ess_tail")

# Posterior extraction
post2 <- fit2$draws(variables = c("alpha", "beta", "k_mu"), format = "df")

k2 <- fit2$draws(variables = "k", format = "matrix")
prev2 <- fit2$draws(variables = "prev", format="matrix")
prev.mean <- apply(prev2, 2, mean)
k2.mean <- apply(k2, 2, mean)

#Worm burden by age
curve(worm_age(x, alpha=mean(post2$alpha), beta=mean(post2$beta), mu=d_2$mu), 
      from=0, to=60, ylab="Worm burden", xlab="Age (years)", 
      lwd=1.2, axes=F, cex.lab=1.2)
axis(1, lwd=1.2, cex.axis=1.2)
axis(2, las=2, cex.axis=1.2)

#######################
## Prevalence by age ##
#######################
M_post_mat <- matrix(data=0, ncol=length(post2$alpha), nrow=length(x_age))
prev_post_mat <- matrix(data=0, ncol=length(post2$alpha), nrow=length(x_age))
  
for(a in 1:length(x_age)){
  for(i in 1:length(post2$alpha)){
    M_post_mat[a,i] = worm_age(x_age[a], alpha=post2$alpha[i], beta=post2$beta[i], mu=post$mu[i])
    prev_post_mat[a,i] = 1-dnbinom(0, mu=M_post_mat[a,i], size=post2$k_mu[i])
  }
}

prev_post <- apply(prev_post_mat, 1, median)

#Mean p
plot(p_mean ~ x_age, type="l", ylim=c(0,1), axes=F, xlab="Age (years)", ylab="Liver fluke prevalence",
     cex.lab=1.6, lwd=2)
lines(prev_post ~ x_age, lwd=2, lty=2)
legend(x=c(1,5), y=c(0.95,1.05), lty=c(1,2), lwd=2.5, bty="n", cex=1.3,
       legend = c("Pre-intervention (born 1960)", "Post-intervention (born 1990)"))

axis(1, cex.axis=1.4, lwd=1.6)
axis(2, cex.axis=1.4, lwd=1.6, las=2)

max(p_mean)

#integrals - cumulative number of worms acquired
more_age_x <- seq(0,30,0.1)
cworm_pre <- matrix(data=0, nrow=length(more_age_x), ncol=length(post$alpha))
cworm_post <- matrix(data=0, nrow=length(more_age_x), ncol=length(post2$alpha))
for(i in 1:length(more_age_x)){
  for(j in 1:length(post$alpha)){
    cworm_pre[i,j] = foi_def_integral(more_age_x[i], alpha=post$alpha[j], beta=post$beta[j])
  }
  for(j in 1:length(post2$alpha)){
    cworm_post[i,j] = foi_def_integral(more_age_x[i], alpha=post2$alpha[j], beta=post2$beta[j])
  }
}

cworm_pre_median <- apply(cworm_pre, 1, median)
cworm_post_median <- apply(cworm_post, 1, median)
idx_25 <- which(more_age_x==25)
cworm_pre_median[idx_25]
cworm_post_median[idx_25]

plot(cworm_pre_median ~ more_age_x, type="l")

#Cumulative burden in the top ten percent
wb_decile <- c()
for(i in 1:length(post2$alpha)){
  wb_decile[i] <- quantile(rnbinom(1000, 
                                   mu=cworm_post_median[which(more_age_x==25)], 
                                   size=post2$k_mu[i]),
                                  probs=0.9)
}
summary(wb_decile)

#Integral of later FoI
max(prev.mean)/2
foi_def_integral(9, alpha=median(post2$alpha), beta=median(post2$beta))
1-dnbinom(0, mu=0.329, size=median(k2[,1]))

#Lower
foi_def_integral(8, alpha=median(post2$alpha), beta=median(post2$beta))
foi_sim_8yr <- foi_sim_7yr <- c()

for(i in 1:length(post2$alpha)){
  #Upper age estimate
  foi_sim_8yr[i] = foi_def_integral(8, alpha=post2$alpha[i], beta=post2$beta[i])
  foi_sim_7yr[i] = foi_def_integral(7, alpha=post2$alpha[i], beta=post2$beta[i])
}

1-dnbinom(0, mu=quantile(foi_sim_8yr, 0.95), size=quantile(k2[,1], 0.95))
1-dnbinom(0, mu=quantile(foi_sim_7yr, 0.95), size=quantile(k2[,1], 0.5))

#############################
## Age of first infection ##
############################

#Simulate uncertainty in FOI from pre-intervention data
foi_pre_mat <- matrix(data = 0, nrow = length(post$alpha), ncol = 121)
prop_pre_mat <- matrix(data = 0, nrow = length(post$alpha), ncol = 121)

for(i in 1:length(post$alpha)){
  for(j in 1:121){
    #From 0 month to 24 months
    foi_pre_mat[i,j] <- foi_def_integral(((j-1)/12), alpha=post$alpha[i], beta=post$beta[i])
    #proportion infected
    prop_pre_mat[i,j] <- 1-dnbinom(0, mu=foi_pre_mat[i,j], size=post$k_mu[i])
  }
}

#Get median and 90% Credible interval
prep_pre_median <- apply(prop_pre_mat, 2, median)/max(p_mean)
prop_interval_90 <- apply(prop_pre_mat/max(p_mean), 2, 
                          FUN=function(x) quantile(x, probs=c(0.05,0.95)))

#Simulate uncertainty in FOI in post-intervention data
foi_post_mat <- matrix(data = 0, nrow = length(post2$alpha), ncol = 22)
prop_post_mat <- matrix(data = 0, nrow = length(post2$alpha), ncol = 22)

for(i in 1:length(post2$alpha)){
  for(j in 1:22){
    #From 0 to 10 years in 6 month intervals
    foi_post_mat[i,j] <- foi_def_integral((j-1)/2, alpha=post2$alpha[i], beta=post2$beta[i])
    #proportion infected
    prop_post_mat[i,j] <- 1-dnbinom(0, mu=foi_post_mat[i,j], size=post2$k_mu[i])
  }
}

prop_interval_post_90 <- apply(prop_post_mat/max(p_mean_post), 2, 
                               FUN=function(x) quantile(x, probs=c(0.05,0.95), na.rm=T))

prep_post_median <- apply(prop_post_mat/max(p_mean_post), 2, median)

###############
## Figure 2B ##
###############
# Plot both time periods together as a relative proportion of max prevalence
# Using area under the foi curve (cumulative exposure)
#pdf("exposed-prop-by-age.pdf",  height=6, width=7.5)
par(mar=c(5,6,3,2), xpd=F)
plot(1, type="n", axes=F, cex.lab=1.6,
     ylim=c(0, 1),
     xlim=c(0,10),
     xlab="Age (years)",
     ylab="", xaxs="i", 
     cex.main=1.6, yaxs="i",
     main="First exposure to liver fluke by age")
rethinking::shade(prop_interval_90, c((0:120)/12), col=alpha("maroon4", 0.35))
rethinking::shade(prop_interval_post_90, c((0:21)/2), col=alpha("slateblue4", 0.4))
lines(prep_pre_median ~ c((0:120)/12), lwd=2.8)
lines(prep_post_median ~ c((0:21)/2), lwd=2.8, col="black")
#axis(1, at=seq(1,13,1), labels=rep("",13), lwd=1.2, tck=-0.012)
axis(1, at=seq(0,10,1), labels=seq(0,10,1), cex.axis=1.45, lwd=1.6)
axis(2, las=2, cex.axis=1.45, lwd=1.6, at=c(0,0.25, 0.5,0.75,1))
title(ylab = expression(paste("Proportion infected with", "">=1, " worm", sep="")), 
      line=3.8, cex.lab=1.6)
abline(h=0.5, lty=2, lwd=2.8)
legend(x=c(0.2,3.2), y=c(0.98,1.035), 
       legend=c("Pre-intervention (born 1960)", "Post-intervention (born 1990)"), cex=1.1, bty="n",
       lty=c(1,1), col=c(alpha("maroon4", 0.8), alpha("slateblue4", 0.8)), lwd=6)
#dev.off()

which.min(abs(prep_pre_median - 0.5)) #index 23 = 22 months
which.min(abs(prop_interval_90[1,] - 0.5)) #index 26 = 25 months (lower)
which.min(abs(prop_interval_90[2,] - 0.5)) #index 20 = 19 months (upper)

which.min(abs(prep_post_median - 0.5)) #index 17 = 8.0 years
which.min(abs(prop_interval_post_90[1,] - 0.5)) #index 20 = 9.5 years
which.min(abs(prop_interval_post_90[2,] - 0.5)) #index 14 = 6.5 years

############################ 
## Cumulative worm burden ##
############################

cwb_1960_mean <- cwb_1960_dist <- cwb_1960_upper <- matrix(data=0, nrow=length(post$alpha), ncol=61) #cumulative mean burden up to 30 years
cwb_1990_mean <- cwb_1990_dist <- cwb_1990_upper <- matrix(data=0, nrow=length(post$alpha), ncol = 61) #cumulative mean burden up to 30 years

set.seed(1234)
for(i in 1:length(post$alpha)){
  for(j in 1:61){
    #From 0 to 30 years in 6 month intervals
    cwb_1960_mean[i,j] <- foi_def_integral((j-1)/2, alpha=post$alpha[i], beta=post$beta[i])
    cwb_1990_mean[i,j] <- foi_def_integral((j-1)/2, alpha=post2$alpha[i], beta=post2$beta[i])
    #proportion infected
    sims_1960 <- rnbinom(2500, mu=cwb_1960_mean[i,j], size=post$k_mu[i])
    sims_1990 <- rnbinom(2500, mu=cwb_1990_mean[i,j], size=post2$k_mu[i])
    #Get distributional quantiles
    cwb_1960_dist[i,j] <- median(sims_1960)
    cwb_1990_dist[i,j] <- median(sims_1990)
    cwb_1960_upper[i,j] <- quantile(sims_1960, probs=c(0.9))
    cwb_1990_upper[i,j] <- quantile(sims_1990, probs=c(0.9))
  }
}

#Get uncertainty intervals and median for plotting
cwb_1960_int <- apply(cwb_1960_dist, 2, FUN=function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
cwb_1960_upint <- apply(cwb_1960_upper, 2, FUN=function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))

cwb_1990_int <- apply(cwb_1990_dist, 2, FUN=function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
cwb_1990_upint <- apply(cwb_1990_upper, 2, FUN=function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))

#Plot
x_cwb <- seq(0,30, 0.5)

###############
## Figure 2C ##
###############
#pdf("cumulative-exposure-by-age.pdf", height=6, width=7.5)
par(mar=c(5,6,3,2), xpd=F, mfrow=c(1,1))
plot(1, type="n", axes=F, cex.lab=1.6,
     ylim=c(0, log(800+1)),
     xlim=c(0,30),
     xlab="Age (years)",
     ylab="", #xaxs="i", 
     cex.main=1.6, #yaxs="i",
     main="Cumulative exposure to liver fluke by age")
axis(1, cex.axis=1.45, lwd=1.5)
axis(2, las=2, cex.axis=1.45, lwd=1.5, at=log(c(0,1,10,100,600)+1),
     labels=c(0,1,10,100,600))
title(ylab = "Cumulative number of worms", 
      line=4, cex.lab=1.6)

rethinking::shade(log(cwb_1960_int[c(1,3),]+1), x_cwb, col=alpha("maroon4", 0.35))
lines(log(cwb_1960_int[2,]+1)  ~ x_cwb, lwd=2.8)
rethinking::shade(log(cwb_1990_upint[c(1,3),]+1), x_cwb, col=alpha("slateblue4", 0.5))
lines(log(cwb_1990_upint[2,]+1)  ~ x_cwb, lwd=2.8, lty=4)
rethinking::shade(log(cwb_1960_upint[c(1,3),]+1), x_cwb, col=alpha("maroon4", 0.35))
lines(log(cwb_1960_upint[2,]+1)  ~ x_cwb, lwd=2.8, lty=4)
legend(x=c(0,2), y=log(c(900,1200)), 
       legend=c("Pre-intervention (born 1960)", 
                "Post-intervention (born 1990)",
                "Median",
                "90% percentile"), 
       cex=1.1, bty="n", lty=c(1,1,1,4), lwd=c(6,6,2.5,2.5),
       col=c(alpha("maroon4", 0.8), alpha("slateblue4", 0.8), "black", "black"))
#dev.off()
