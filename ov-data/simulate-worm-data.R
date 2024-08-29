## Generate data on worm burdens by age for dynamic model

#Round function
round2 = function(x, n=0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

#Generate  values from beta distribution
#https://stackoverflow.com/questions/37773469/r-random-distribution-with-predefined-min-max-mean-and-sd-values
rgbeta <- function(n, mean, var, min = 0, max = 1){
  dmin <- mean - min
  dmax <- max - mean
  
  if (dmin <= 0 || dmax <= 0){
    stop(paste("mean must be between min =", min, "and max =", max)) 
  }
  
  if (var >= dmin * dmax){
    stop(paste("var must be less than (mean - min) * (max - mean) =", dmin * dmax))
  }
  
  # mean and variance of the standard beta distributed variable
  mx <- (mean - min) / (max - min)
  vx <- var / (max - min)^2
  
  # find the corresponding alpha-beta parameterization
  a <- ((1 - mx) / vx - 1 / mx) * mx^2
  b <- a * (1 / mx - 1)
  
  # generate standard beta observations and transform
  x <- rbeta(n, a, b)
  y <- (max - min) * x + min
  
  return(y)
}

se_to_var <- function(se, n) se^2*n

##################################
# Datasets on worm burden by age #
##################################

# Study 1) Ramsay et al. 1989 | Expulsion
ramsay <- read.delim("Google Drive/My Drive/Opisthorchis/worm-output-digitise/extracted-data/Ramsay-expulsion.txt")

# Study 2) Sithithaworn et al 1991 | Autopsy 
# Quantitative post-mortem study of Opisthorchis viverrini in man in north- east Thailand
#set seed
set.seed(12345)

#cases -  ages
age_cases <- list()
age_cases[[1]] <- rgbeta(9, mean = 5, var = 8, min = 2, max = 10.4)
age_cases[[2]] <- rgbeta(31, mean = 15, var = 8, min = 10.5, max = 20.4)
age_cases[[3]] <- rgbeta(48, mean = 25, var = 8, min = 20.5, max = 30.4)
age_cases[[4]] <- rgbeta(25, mean = 35, var = 8, min = 30.5, max = 40.4)
age_cases[[5]] <- rgbeta(21, mean = 45, var = 8, min = 40.5, max = 50.4)
age_cases[[6]] <- rgbeta(25, mean = 56, var = 22, min = 50.5, max = 78)

#cases - worms
worms_cases <- list()
worms_cases[[1]] <- c(rep(0, 6), 
                      round2(rgbeta(3, mean = 4, var = 10, min = 1, max = 50)))

worms_cases[[2]] <- round2(c(rep(0, 6), 
                             rgbeta(16, mean = 20, var = 300, min = 1, max = 50),
                             rgbeta(1, mean = 75, var = 500, min = 51, max = 100),
                             rgbeta(3, mean = 150, var = 1000, min = 101, max = 200),
                             rgbeta(2, mean = 300, var = 3000, min = 201, max = 400),
                             rgbeta(3, mean = 720, var = 6000, min = 401, max = 2500)))

worms_cases[[3]] <- round2(c(rep(0,8), 
                             rgbeta(17, mean = 30, var = 300, min = 1, max = 50),
                             rgbeta(6, mean = 80, var = 500, min = 51, max = 100),
                             rgbeta(9, mean = 180, var = 1000, min = 101, max = 200),
                             rgbeta(3, mean = 350, var = 3000, min = 201, max = 400),
                             rgbeta(5, mean = 2000, var = 6000, min = 401, max = 2800)))

worms_cases[[4]] <- round2(c(rep(0, 3), 
                             rgbeta(8, mean = 25, var = 300, min = 1, max = 50),
                             rgbeta(1, mean = 75, var = 500, min = 51, max = 100),
                             rgbeta(5, mean = 150, var = 1000, min = 101, max = 200),
                             rgbeta(3, mean = 300, var = 3000, min = 201, max = 400),
                             rgbeta(5, mean = 1200, var = 6000, min = 401, max = 2100)))

worms_cases[[5]] <- round2(c(rep(0, 3), 
                             rgbeta(10, mean = 20, var = 300, min = 1, max = 50),
                             rgbeta(3, mean = 75, var = 500, min = 51, max = 100),
                             rgbeta(2, mean = 150, var = 1000, min = 101, max = 200),
                             rgbeta(2, mean = 300, var = 3000, min = 201, max = 400)))

worms_cases[[6]] <- round2(c(rep(0, 7), 
                             rgbeta(9, mean = 25, var = 300, min = 1, max = 50),
                             rgbeta(2, mean = 75, var = 500, min = 51, max = 100),
                             rgbeta(4, mean = 150, var = 1000, min = 101, max = 200),
                             rgbeta(4, mean = 320, var = 2000, min = 201, max = 400)))


# Study 3) Haswell Elkins et al 1991 | expulsion
#Distribution patterns of Opisthorchis viverrini within a human community

#ages
age_cases_E <- list()
age_cases_E[[1]] <- rgbeta(117, mean = 7, var = 3.5, min = 5, max = 9) #5-9
age_cases_E[[2]] <- rgbeta(68, mean = 15, var = 8, min = 10, max = 19) #10-19
age_cases_E[[3]] <- rgbeta(53, mean = 25, var = 8, min = 20, max = 29) #20-29
age_cases_E[[4]] <- rgbeta(43, mean = 35, var = 8, min = 30, max = 39) #30-39
age_cases_E[[5]] <- rgbeta(41, mean = 45, var = 8, min = 40, max = 49) #40-49
age_cases_E[[6]] <- rgbeta(51, mean = 55, var = 8, min = 50, max = 60) #50-60

#worms
worm_burden_E <- list()
worm_burden_E[[1]] <- round2(c(rep(0, 107), 
                        rgbeta(7, mean = 5, var = 4, min = 1, max = 9), 
                        rgbeta(3, mean = 20, var = 8, min = 10, max = 29)))

worm_burden_E[[2]] <- round2(c(rep(0, 38), 
                        rgbeta(19, mean = 5, var = 4, min = 1, max = 9), 
                        rgbeta(6, mean = 20, var = 8, min = 10, max = 29),
                        rgbeta(4, mean = 65, var = 100, min = 30, max = 99),
                        rgbeta(1, mean = 300, var = 6000, min = 100, max = 830)))

worm_burden_E[[3]] <- round2(c(rep(0, 25), 
                        rgbeta(11, mean = 5, var = 4, min = 1, max = 9), 
                        rgbeta(7, mean = 20, var = 8, min = 10, max = 29),
                        rgbeta(7, mean = 65, var = 100, min = 30, max = 99),
                        rgbeta(3, mean = 300, var = 6000, min = 100, max = 830)))

worm_burden_E[[4]] <- round2(c(rep(0, 17), 
                        rgbeta(10, mean = 5, var = 4, min = 1, max = 9), 
                        rgbeta(2, mean = 20, var = 8, min = 10, max = 29),
                        rgbeta(6, mean = 65, var = 100, min = 30, max = 99),
                        rgbeta(8, mean = 300, var = 6000, min = 100, max = 830)))

worm_burden_E[[5]] <- round2(c(rep(0, 17), 
                        rgbeta(5, mean = 5, var = 4, min = 1, max = 9), 
                        rgbeta(6, mean = 20, var = 8, min = 10, max = 29),
                        rgbeta(5, mean = 65, var = 100, min = 30, max = 99),
                        rgbeta(8, mean = 300, var = 6000, min = 100, max = 830)))

worm_burden_E[[6]] <- round2(c(rep(0, 22), 
                        rgbeta(12, mean = 5, var = 4, min = 1, max = 9), 
                        rgbeta(6, mean = 20, var = 8, min = 10, max = 29),
                        rgbeta(4, mean = 65, var = 100, min = 30, max = 99),
                        rgbeta(7, mean = 300, var = 6000, min = 100, max = 830)))

# Study 4) Upatham et al 1984 | egg counts
# Relationship between prevalence and intensity of Opisthorchis viverrini infection
# and clinical symptoms and signs in a rural community in north-east Thailand
age_cases_U <- list()
egg_burden_U <- list()

#Simulate individual ages within categories - Brockelman table
age_cases_U[[1]] <- rep(1, 46) #1
age_cases_U[[2]] <- rep(2, 36) #2
age_cases_U[[3]] <- rep(3, 42) #3
age_cases_U[[4]] <- rep(4, 31) #4
age_cases_U[[5]] <- rep(5, 43) #5
age_cases_U[[6]] <- rep(6, 33) #6
age_cases_U[[7]] <- rep(7, 45) #7
age_cases_U[[8]] <- rep(8, 45) #8
age_cases_U[[9]] <- rep(9, 52) #9
age_cases_U[[10]] <- rgbeta(228, mean = 12, var = 2.5, min = 9.5, max = 14.4) #10-14
age_cases_U[[11]] <- rgbeta(205, mean = 17, var = 2.5, min = 14.5, max = 19.4) #15-19
age_cases_U[[12]] <- rgbeta(115, mean = 22, var = 2.0, min = 19.5, max = 24.4) #20-24
age_cases_U[[13]] <- rgbeta(111, mean = 27, var = 2.5, min = 25.5, max = 29.4) #25-29
age_cases_U[[14]] <- rgbeta(197, mean = 35, var = 10, min = 29.5, max = 39.4) #30-39
age_cases_U[[15]] <- rgbeta(162, mean = 45, var = 10, min = 39.5, max = 49.4) #40-49
age_cases_U[[16]] <- rgbeta(130, mean = 55, var = 10, min = 49.5, max = 59.4) #50-59
age_cases_U[[17]] <- rgbeta(94, mean = 65, var = 10, min = 59.5, max = 69.4) #60-69
age_cases_U[[18]] <- rgbeta(36, mean = 74, var = 6, min = 69.5, max = 78) #70+

#Simulate individual egg counts within categories
egg_burden_U[[1]] <- c(rep(0, 41), round2(rgbeta(5, mean = 2.9, var = se_to_var(1.7, 5), min = 0.006, max = 500)*1000/6)) #0-1
egg_burden_U[[2]] <- c(rep(0, 26), round2(rgbeta(10, mean = 2.9, var = se_to_var(1.7, 10), min = 0.006, max = 500)*1000/6)) #2
egg_burden_U[[3]] <- c(rep(0, 24), round2(rgbeta(18, mean = 2.9, var = se_to_var(1.7, 18), min = 0.006, max = 500)*1000/6)) #3
egg_burden_U[[4]] <- c(rep(0, 14), round2(rgbeta(17, mean = 2.9, var = se_to_var(1.7, 17), min = 0.006, max = 500)*1000/6)) #4
egg_burden_U[[5]] <- c(rep(0, 7), round2(rgbeta(36, mean = 9.3, var = se_to_var(2.0, 36), min = 0.006, max = 500)*1000/6)) #5
egg_burden_U[[6]] <- c(rep(0, 3), round2(rgbeta(30, mean = 9.3, var = se_to_var(2.0, 30), min = 0.006, max = 500)*1000/6)) #6
egg_burden_U[[7]] <- c(rep(0, 5), round2(rgbeta(40, mean = 9.3, var = se_to_var(2.0, 40), min = 0.006, max = 500)*1000/6)) #7
egg_burden_U[[8]] <- c(rep(0, 4), round2(rgbeta(41, mean = 9.3, var = se_to_var(2.0, 41), min = 0.006, max = 500)*1000/6)) #8
egg_burden_U[[9]] <- c(rep(0, 2), round2(rgbeta(50, mean = 9.3, var = se_to_var(2.0, 50), min = 0.006, max = 500)*1000/6)) #9
egg_burden_U[[10]] <- c(rep(0, 10), round2(rgbeta(218, mean = 10.8, var = se_to_var(2.0, 218), min = 0.006, max = 500)*1000/6)) #10-14
egg_burden_U[[11]] <- c(rep(0, 8), round2(rgbeta(197, mean = 17.6, var = se_to_var(1.8, 197), min = 0.006, max = 500)*1000/6)) #15-19
egg_burden_U[[12]] <- c(rep(0, 6), round2(rgbeta(109, mean = 22.9, var = se_to_var(3.0, 109), min = 0.006, max = 500)*1000/6)) #20-24
egg_burden_U[[13]] <- c(rep(0, 2), round2(rgbeta(109, mean = 20.5, var = se_to_var(2.9, 109), min = 0.006, max = 500)*1000/6)) #25-29
egg_burden_U[[14]] <- c(rep(0, 7), round2(rgbeta(190, mean = 22.7, var = se_to_var(2.3, 190), min = 0.006, max = 500)*1000/6)) #30-39
egg_burden_U[[15]] <- c(rep(0, 5), round2(rgbeta(157, mean = 31.7, var = se_to_var(3.8, 157), min = 0.006, max = 500)*1000/6)) #40-49
egg_burden_U[[16]] <- c(rep(0, 4), round2(rgbeta(126, mean = 29.1, var = se_to_var(3.8, 126), min = 0.006, max = 500)*1000/6)) #50-59
egg_burden_U[[17]] <- c(rep(0, 3), round2(rgbeta(91, mean = 26.9, var = se_to_var(4.7, 91), min = 0.006, max = 500)*1000/6)) #60-69
egg_burden_U[[18]] <- c(rep(0, 2), round2(rgbeta(34, mean = 26.9, var = se_to_var(4.7, 34), min = 0.006, max = 500)*1000/6)) #70+

length(unlist(age_cases_U))
length(unlist(egg_burden_U))
sapply(age_cases_U, length)
sapply(egg_burden_U, length)

#Geometric mean corrections
am <- c(2.9, 9.3, 10.8, 17.6, 22.9, 20.5, 22.7, 31.7, 29.1, 26.9)
gm <- c(0.5, 2.9, 4.9, 8.2, 8.5, 8.7, 9.4, 12.2, 13.2, 11.9)
#par(mar=c(5,5,3,3))
#plot(am~gm, ylim=c(0,35), xlim=c(0,14), pch=16, cex=1.4, 
#     ylab="Arithmetic Mean", xlab="Geometic Mean")
m1 <- lm(am~gm-1)
#abline(a=0, b=m1$coefficients[1])
gmx <- m1$coefficients[1]
#multiply gm by 2.38 to get arithmetic mean

# Study 5A) Kurathong et al 1987 | egg counts
# Opisthorchis viverfini infection in
#rural and urban communities in northeast Thailand

## Note the transformation applied to the geometric mean

## Rural population (n = 433)
age_cases_rural <- list()
egg_burden_rural <- list()

age_cases_rural[[1]] <- rgbeta(17, mean = 2, var = 1.5, min = 1, max = 4.9) #0-4
age_cases_rural[[2]] <- rgbeta(18, mean = 7, var = 2.5, min = 5, max = 9.9) #5-9
age_cases_rural[[3]] <- rgbeta(55, mean = 15, var = 7.5, min = 10, max = 19.9) #10-19
age_cases_rural[[4]] <- rgbeta(73, mean = 25, var = 10, min = 20, max = 29.9) #20-29
age_cases_rural[[5]] <- rgbeta(77, mean = 35, var = 10, min = 30, max = 39.9) #30-39
age_cases_rural[[6]] <- rgbeta(72, mean = 45, var = 10, min = 40, max = 49.9) #40-49
age_cases_rural[[7]] <- rgbeta(74, mean = 55, var = 10, min = 50, max = 59.9) #50-59
age_cases_rural[[8]] <- rgbeta(47, mean = 64, var = 7.5, min = 60, max = 70) #60+

egg_burden_rural[[1]] <- c(rep(0, 15), 
                           round2(rgbeta(2, mean = 0.2*gmx, var = se_to_var(0.1, 2), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[2]] <- c(rep(0, 9), 
                           round2(rgbeta(9, mean = 11.6*gmx, var = se_to_var(9.8, 9), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[3]] <- c(rep(0, 11), 
                           round2(rgbeta(44, mean = 7.5*gmx, var = se_to_var(2.1, 44), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[4]] <- c(rep(0, 17), 
                           round2(rgbeta(56, mean = 15.2*gmx, var = se_to_var(4.4, 56), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[5]] <- c(rep(0, 9), 
                           round2(rgbeta(68, mean = 17.4*gmx, var = se_to_var(4.7, 68), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[6]] <- c(rep(0, 10), 
                           round2(rgbeta(62, mean = 30.7*gmx, var = se_to_var(7.9, 62), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[7]] <- c(rep(0, 10), 
                           round2(rgbeta(64, mean = 21.0*gmx, var = se_to_var(6.7, 64), min = 0.006, max = 500)*1000/6))

egg_burden_rural[[8]] <- c(rep(0, 8), 
                           round2(rgbeta(39, mean = 35.8*gmx, var = se_to_var(3, 39), min = 0.006, max = 500)*1000/6))

length(unlist(age_cases_rural))
length(unlist(egg_burden_rural))
sapply(age_cases_rural, length)
sapply(egg_burden_rural, length)

# Study 5B) Kurathong et al 1987 | egg counts
# Opisthorchis viverfini infection in
#rural and urban communities in northeast Thailand

## Note the transformation applied to the geometric mean

## Urban population (n = 126)

age_cases_urban <- list()
egg_burden_urban <- list()

age_cases_urban[[1]] <- rgbeta(6, mean = 2, var = 1.5, min = 1, max = 4.9) #0-5
age_cases_urban[[2]] <- rgbeta(1, mean = 7, var = 2.5, min = 5, max = 9.9) #5-9
age_cases_urban[[3]] <- rgbeta(21, mean = 15, var = 7.5, min = 10, max = 19.9) #10-19
age_cases_urban[[4]] <- rgbeta(23, mean = 25, var = 10, min = 20, max = 29.9) #20-29
age_cases_urban[[5]] <- rgbeta(27, mean = 35, var = 10, min = 30, max = 39.9) #30-39
age_cases_urban[[6]] <- rgbeta(13, mean = 45, var = 10, min = 40, max = 49.9) #40-49
age_cases_urban[[7]] <- rgbeta(26, mean = 55, var = 10, min = 50, max = 59.9) #50-59
age_cases_urban[[8]] <- rgbeta(9, mean = 64, var = 7.5, min = 60, max = 70) #60+

egg_burden_urban[[1]] <- rep(0, 6) #0-4
egg_burden_urban[[2]] <- rep(0, 1) #5-9
egg_burden_urban[[3]] <- c(rep(0, 10), 
                           round2(rgbeta(11, mean = 5.1*gmx, var = se_to_var(3.1, 11), min = 0.006, max = 100)*1000/6)) #10-19
egg_burden_urban[[4]] <- c(rep(0, 6), 
                           round2(rgbeta(17, mean = 3.1*gmx, var = se_to_var(1.2, 17), min = 0.006, max = 100)*1000/6)) #20-29
egg_burden_urban[[5]] <- c(rep(0, 16), 
                           round2(rgbeta(11, mean = 6.9*gmx, var = se_to_var(5.9, 11), min = 0.006, max = 100)*1000/6)) #30-39
egg_burden_urban[[6]] <- c(rep(0, 6), 
                           round2(rgbeta(7, mean = 13.3*gmx, var = se_to_var(10.1, 7), min = 0.006, max = 100)*1000/6)) #40-49
egg_burden_urban[[7]] <- c(rep(0, 10), 
                           round2(rgbeta(16, mean = 1.8*gmx, var = se_to_var(0.6, 16), min = 0.006, max = 100)*1000/6)) #50-59
egg_burden_urban[[8]] <- c(rep(0, 2), 
                           round2(rgbeta(7, mean = 2.4*gmx, var = se_to_var(1.4, 7), min = 0.006, max = 100)*1000/6)) #60+

length(unlist(egg_burden_urban))
length(unlist(age_cases_urban))


###############################
## Post intervention surveys ##
##############################

# Study X1) Jongsuksuntigul et al 1997 | egg counts
# The impact of a decade long opisthorchiasis control program in Northeastern Thailand
age_cases_J <- list()
egg_burden_J <- list()

age_cases_J[[1]] <- rgbeta(5, mean = 12, var = 3, min = 9.5, max = 14.4) #10-14
age_cases_J[[2]] <- rgbeta(1, mean = 17, var = 3, min = 14.5, max = 19.4) #15-19
age_cases_J[[3]] <- rgbeta(7, mean = 25, var = 10, min = 19.5, max = 29.4) #20-29
age_cases_J[[4]] <- rgbeta(16, mean = 35, var = 10, min = 29.5, max = 39.4) #30-39
age_cases_J[[5]] <- rgbeta(19, mean = 45, var = 10, min = 39.5, max = 49.4) #40-49
age_cases_J[[6]] <- rgbeta(12, mean = 55, var = 10, min = 49.5, max = 59.4) #50-59
age_cases_J[[7]] <- rgbeta(5, mean = 66, var = 20, min = 59.5, max = 80) #60+

egg_burden_J[[1]] <- round2(rgbeta(5, mean = 230, var = 40000, min = 40, max = 999)) #10-14
egg_burden_J[[2]] <- round2(rgbeta(1, mean = 92, var = 92, min = 40, max = 999)) #15-19
egg_burden_J[[3]] <- round2(c(rgbeta(6, mean = 227, var = 40000, min = 40, max = 999), 
                              rgbeta(1, mean = 2000, var = 2000, min = 1000, max = 3000))) #20-29
egg_burden_J[[4]] <- round2(c(rgbeta(13, mean = 281, var = 40000, min = 40, max = 999),
                              rgbeta(3, mean = 2000, var = 6000, min = 1000, max = 3000))) #30-39
egg_burden_J[[5]] <- round2(c(rgbeta(16, mean = 300, var = 40000, min = 40, max = 999),
                              rgbeta(1, mean = 5000, var = 5000, min = 1000, max = 9000),
                              rgbeta(2, mean = 21200, var = 25000, min = 10000, max = 30000))) #40-49
egg_burden_J[[6]] <- round2(c(rgbeta(8, mean = 320, var = 40000, min = 40, max = 999),
                              rgbeta(3, mean = 2500, var = 5000, min = 1000, max = 9000),
                              rgbeta(1, mean = 10500, var = 12000, min = 10000, max = 15000))) #50-59
egg_burden_J[[7]] <- round2(rgbeta(5, mean = 275, var = 40000, min = 40, max = 999)) #60+

######################
## Set up data.frames 
######################
autopsy_data <- data.frame(age=unlist(age_cases), worms=unlist(worms_cases),
                           study_id=1, study_type=1, year=1985)
ramsay_data <- data.frame(age=ramsay$age, worms=ramsay$flukes_recovered,
                          study_id=2, study_type=2, year=1987)
elkins_data <- data.frame(age=unlist(age_cases_E), worms=unlist(worm_burden_E),
                          study_id=2, study_type=2, year=1989)

worm_data <- rbind(autopsy_data, ramsay_data, elkins_data) #merge data.frames

worm_data$age_class <- ifelse(worm_data$age<5, 1,
                  ifelse(worm_data$age<10, 2, 
                  ifelse(worm_data$age<20, 3,
                  ifelse(worm_data$age<30, 4,
                  ifelse(worm_data$age<40, 5,
                  ifelse(worm_data$age<50, 6, 
                  ifelse(worm_data$age<60, 7, 8)))))))

#Same for egg counts
#Pre-intervention egg data
upatham_data <- data.frame(age=unlist(age_cases_U), eggs=unlist(egg_burden_U), 
                           study_id=3, year=1980)
urban_data <- data.frame(age=unlist(age_cases_urban), eggs=unlist(egg_burden_urban),
                         study_id=4, year=1982)
rural_data <- data.frame(age=unlist(age_cases_rural), eggs=unlist(egg_burden_rural),
                         study_id=5, year=1982)

egg_data <- rbind(upatham_data, urban_data, rural_data)

egg_data$age_class <- ifelse(egg_data$age<5, 1,
                          ifelse(egg_data$age<10, 2, 
                          ifelse(egg_data$age<20, 3,
                          ifelse(egg_data$age<30, 4,
                          ifelse(egg_data$age<40, 5,
                          ifelse(egg_data$age<50, 6,
                          ifelse(egg_data$age<60, 7, 8)))))))

#Post-intervention egg data
jong_data <- data.frame(age=unlist(age_cases_J), eggs=unlist(egg_burden_J), 
                        study_id=2, year=1994)

jong_data$age_class <- ifelse(jong_data$age<5, 1,
                          ifelse(jong_data$age<10, 2, 
                          ifelse(jong_data$age<20, 3,
                          ifelse(jong_data$age<30, 4,
                          ifelse(jong_data$age<40, 5,
                          ifelse(jong_data$age<50, 6,
                          ifelse(jong_data$age<60, 7, 8)))))))
#Write tables
write.table(worm_data, "Google Drive/My Drive/Opisthorchis/CCA-evolution/data/worms-age-data.csv",
            sep=",", quote=F, row.names = F)

write.table(jong_data, "Google Drive/My Drive/Opisthorchis/CCA-evolution/data/jong-eggs-age-data.csv",
            sep=",", quote=F, row.names = F)

#write.table(upatham_data, "Google Drive/My Drive/Opisthorchis/CCA-evolution/upatham-eggs-age-data.csv",
#            sep=",", quote=F, row.names = F)

write.table(egg_data, "Google Drive/My Drive/Opisthorchis/CCA-evolution/data/eggs-age-data.csv",
            sep=",", quote=F, row.names = F)

#rm(list=ls(all=T))
#gc()
