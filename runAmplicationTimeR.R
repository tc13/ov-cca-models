# Run AmplicationTimeR - https://github.com/Wedge-lab/AmplificationTimeR 

library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(scales)
library(cmdstanr)

#Meta data file with cellularity
meta <- read.delim("Desktop/MutationTimeR/master.dpclust.txt")
supp <- read.delim("Desktop/MutationTimeR/jusakul-TH.txt", sep="\t")

#Size of chromosomes
chrs.name <- c((1:22), "X")
chrs.length <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23]

#Load input data into lists
muts <- cn <- mult <- list()

for(i in 1:22){
muts[[i]] <- read.table(paste0("Desktop/MutationTimeR/muts/TH_",i,".muts.txt"), sep="\t")
colnames(muts[[i]]) <- c("chr", "start", "end","ref", "alt")
cn[[i]] <- read.table(paste0("Desktop/MutationTimeR/subclones/T_CCA_TH_",i,"_subclones.txt"), header=T)
mult[[i]] <- read.table(paste0("Desktop/MutationTimeR/dpclust3p/T_CCA_TH_",i,".dpclust3p.txt"), header=T)
}

#Gene chromosomes coordinates 
tp53 = data.frame(gene="TP53", chr=17, start=7565097, stop=7590856)
kras = data.frame(gene="KRAS", chr=12, start=25358180, stop=25403863)
cdkn2a = data.frame(gene="CDKN2A", chr=9, start=21967751, stop=21995323)
smad4 = data.frame(gene="SMAD4", chr=18, start=48556583, stop=48611412)
pik3ca = data.frame(gene="PIK3CA", chr=3, start=178866145, stop=178957881)
cndkn2b = data.frame(gene="CDKN2B", chr=9, start=22002902, stop=22009312)
erbb2 = data.frame(gene="ERBB2", chr=17, start=37844347, stop=37884911)
arid1a = data.frame(gene="ARID1A", chr=1, start=27022506, stop=27108595)
apc = data.frame(gene="APC", chr=5, start=112043195, stop=112181936)
rpph1 = data.frame(gene="RPPH1", chr=14, start=20811230, stop=20811570)
kdm5a = data.frame(gene="KDM5A", chr=12, start=389223, stop=498486)
robo1 = data.frame(gene="ROBO1", chr=3, start=78646389, stop=79817148)
fbxw7 = data.frame(gene="FBXW7", chr=4, start=153241696, stop=153457244)
ctnnb1 = data.frame(gene="CTNNB1", chr=3, start=41240996, stop=41281934)
tgfbr2 = data.frame(gene="TGFBR2", chr=3, start=30647994, stop=30735634)
gnas = data.frame(gene="GNAS", chr=20, start=57414773, stop=57486247)
braf = data.frame(gene="BRAF", chr=7, start=140419127, stop=140419127)
fgf23 = data.frame(gene="FGF23", chr=12, start=4477393, stop=4488878)
pten = data.frame(gene="PTEN", chr=10, start=89622395, stop=89731687)
rasa1 = data.frame(gene="RASA1", chr=5, start=86563700, stop=86687748)
erbb4 = data.frame(gene="ERBB4", chr=2, start=212240446, stop=213403565)
egfr = data.frame(gene="EGFR", chr=7, start=55086710, stop=55279321)
brca2 = data.frame(gene="BRCA2", chr=13, start=32889223, stop=32974405)
map2k4 = data.frame(gene="MAP2K4", chr=17, start=11924141, stop=12047147)
msh3 = data.frame(gene="MSH3", chr=5, start=79950471, stop=80172634)
stk11 = data.frame(gene="STK11", chr=19, start=1205777, stop=1228430)
nf1 = data.frame(gene="NF1", chr=17, start=29421945, stop=29704695)
fgfr2 = data.frame(gene="FGFR2", chr=10, start=123237844, stop=123357972)
pik3r1 = data.frame(gene="PIK3R1", chr=5, start=67511584, stop=67597649)
nras = data.frame(gene="NRAS", chr=1, start=115247090, stop=115259392)
rnf43 = data.frame(gene="RNF43", chr=17, start=56429861, stop=56494895)
elf3 = data.frame(gene="ELF3", chr=1, start=201979715, stop=201986311)
avcr2a = data.frame(gene="ACVR2A", chr=2, start=148602086, stop=148688391)
bap1 = data.frame(gene="BAP1", chr=3, start=52435024, stop=52444024)
pbrm1 = data.frame(gene="PBRM1", chr=3, start=52579383, stop=52719929)
arid2 = data.frame(gene="ARID2", chr=12, start=46123489, stop=46301820)
erbb3 = data.frame(gene="ERBB3", chr=12, start=56473892, stop=56497289)
sf3b1 = data.frame(gene="SF3B1", chr=2, start=198254508, stop=198299817)
idh1 = data.frame(gene="IDH1", chr=2, start=209100951, stop=209119795)
rb1 = data.frame(gene="RB1", chr=13, start=48877887, stop=49056026)
fgfr1 = data.frame(gene="FGFR1", chr=8, start=38268661, stop=38326153)
cdkn1b = data.frame(gene="CDKN1B", chr=12, start=12870302, stop=12875303)
asxl1 = data.frame(gene="ASXL1", chr=20, start=30946134, stop=31027122)
lamc2 = data.frame(gene="LAMC2", chr=1, start=183155373, stop=183214035)
tert = data.frame(gene="TERT", chr=5, start=1253282, stop=1295183)

#Single data.frame
oncogenes = rbind(tp53, kras, cdkn2a, cndkn2b, arid1a,
                  erbb2, pik3ca, smad4, apc, rpph1,
                  kdm5a, robo1, fbxw7, ctnnb1, tgfbr2,
                  gnas, braf, fgf23, pten, rasa1, erbb4,
                  egfr, brca2, map2k4, msh3, stk11, nf1,
                  fgfr2, pik3r1, nras, rnf43, avcr2a, bap1,
                  erbb3, arid2, pbrm1, sf3b1, idh1, fgfr1, 
                  cdkn1b, rb1, asxl1, lamc2, tert)

#############################
## Run amplification timer ##
#############################

onco_amp <- list()

for(d in 1:nrow(oncogenes)){
  onco_amp[[oncogenes$gene[d]]] = list()
  for(i in 1:22){
  #Check whole-genome duplication status
  WGD <- ifelse(i %in% c(5, 15, 19), TRUE, FALSE)
  
  #Run for all genes in loop 
  test <- try(time_amplification(
    cn_data = cn[[i]],
    multiplicity_data = mult[[i]],
    mutation_data = muts[[i]],
    muts_type = "All",
    sample_id = paste0("TH",i),
    amplification_chrom = oncogenes$chr[d],
    amplification_start = oncogenes$start[d],
    amplification_stop = oncogenes$stop[d],
    is_WGD = WGD,
    genome = "hg19"), TRUE)
  
  if(inherits(test, "try-error")){NULL
    }else{
      onco_amp[[oncogenes$gene[d]]][[i]] <- test
    }
  } 
} #end of loop

#Merge lists
mergeR <- function(list, gene, clonality="clonal") {
  df <- bind_rows(list)
  time <- df[which(!is.na(df$t_1) & df$clonality_status==clonality), c(1:6,8:11)]
  time$gene <- gene
  return(time)
}

all_driveR <- function(list, gene){
  df <- bind_rows(list)
  df$gene <- gene
  time <- reshape2::melt(df, id.vars=c("sample", "region", "gene", "clonality_status", "num_mutations_used"),
               measure.vars=c("t_1_mean_bootstrap", "t_2_mean_bootstrap", "t_3_mean_bootstrap",
                              "t_4_mean_bootstrap", "t_5_mean_bootstrap", "t_6_mean_bootstrap",
                              "t_7_mean_bootstrap", "t_8_mean_bootstrap", "t_9_mean_bootstrap",
                              "t_10_mean_bootstrap"), na.rm=T, value.name="timing")
  return(time)
}

clonal_list <- list()
subclonal_list <- list()
driver_list <- list()

for(d in 1:nrow(oncogenes)){
  geneX = oncogenes$gene[d]
  clonal_list[[geneX]] = mergeR(onco_amp[[geneX]], geneX)
  test <- try(mergeR(onco_amp[[geneX]], geneX, clonality="subclonal"), silent=T)
  if(inherits(test, "try-error")){NULL
    }else{
  subclonal_list[[geneX]] = mergeR(onco_amp[[geneX]], geneX, clonality="subclonal")
    }
  driver_list[[geneX]] = all_driveR(onco_amp[[geneX]], geneX)
}

## data.frames of t_1 for clonal and subclonal
clonal <- bind_rows(clonal_list)
subclonal <- bind_rows(subclonal_list)

#All driver genes all timepoints
drivers <- bind_rows(driver_list)

#Keep only timings within [0,1]
#subclonal.qc <- subclonal[which(subclonal$t_1>0 & subclonal$t_1<1),]
#clonal.qc <- subclonal[which(clonal$t_1>0 & clonal$t_1<1),]
drivers.qc <- drivers[which(drivers$timing>0 & drivers$timing<1 & drivers$num_mutations_used>=10),] #Remove 112 values
clonal.qc <- drivers.qc[which(drivers.qc$clonality_status=="clonal"),]
subclonal.qc <- drivers.qc[which(drivers.qc$clonality_status=="subclonal"),]

sub_ids <- sort(unique(subclonal.qc$sample))
clonal_ids <- sort(unique(clonal.qc$sample))
  
sub_list <- c_list <- list()
for(i in 1:length(sub_ids)){
  idx <- sub_ids[i] 
  id_num <- as.numeric(gsub("TH", "", idx))
  indx <- which(subclonal.qc$sample == idx)
  mut_idx <- which(subclonal.qc$timing[indx] == min(subclonal.qc$timing[indx]))[1]
  sub_list[[idx]] = data.frame(id=idx, clonality="subclonal", 
                               gene=subclonal.qc$gene[indx][mut_idx], 
                               timing=min(subclonal.qc$timing[indx]),
                               num_mutations = min(subclonal.qc$num_mutations_used[indx]),
                               age_mut=min(subclonal.qc$timing[indx])*supp$Age.at.surgery[id_num],
                               age_surgery=supp$Age.at.surgery[id_num],
                               CCA=supp$Anatomical.subtype[id_num],
                               sex=supp$Sex[id_num]) 
}

for(i in 1:length(clonal_ids)){
  idx <- clonal_ids[i] 
  id_num <- as.numeric(gsub("TH", "", idx))
  indx <- which(clonal.qc$sample == idx)
  mut_idx <- which(clonal.qc$timing[indx] == min(clonal.qc$timing[indx]))[1]
  c_list[[idx]] = data.frame(id=idx, clonality="clonal", 
                               gene=clonal.qc$gene[indx][mut_idx], 
                               timing=min(clonal.qc$timing[indx]),
                               age_mut=min(clonal.qc$timing[indx])*supp$Age.at.surgery[id_num],
                               age_surgery=supp$Age.at.surgery[id_num],
                               num_mutations = min(clonal.qc$num_mutations_used[indx]),
                               CCA=supp$Anatomical.subtype[id_num],
                               sex=supp$Sex[id_num])  
}

sub_df <- bind_rows(sub_list)
clon_df <- bind_rows(c_list)

summary(clon_df$age_mut)

#Get chromosome and location info
drivers.qc$chr <- as.numeric(sapply(drivers.qc$region, function(x) strsplit(x, ":")[[1]][1]))
drivers.qc$start <- as.numeric(sapply(drivers.qc$region, function(x) strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[1]]))
drivers.qc$end <- as.numeric(sapply(drivers.qc$region, function(x) strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[2]]))

#################################
## Beta regression for drivers ##
#################################
## load stan file
beta_reg <- file.path("Google Drive/My Drive/Opisthorchis/CCA-evolution", "beta-regression.stan")
bmod <- cmdstan_model(beta_reg, compile = FALSE)
bmod$compile(force_recompile = TRUE, cpp_options = list(stan_threads = TRUE))

#collate data into a list
comb.indx <- match(gsub("TH", "TH_",drivers.qc$sample), gsub("CCA_", "", supp$Sample.ID))
drivers.qc$sex <- supp$Sex[comb.indx]
drivers.qc$subtype <- supp$Anatomical.subtype[comb.indx]
drivers.qc.2 <- drivers.qc[drivers.qc$variable %in% c("t_1_mean_bootstrap","t_2_mean_bootstrap") & drivers.qc$num_mutations_used>=8,]
  
beta_data_1 <- list(
  M = nrow(drivers.qc.2),
  N = length(unique(drivers.qc.2$sample)),
  G = length(unique(drivers.qc.2$gene)),
  person = match(drivers.qc.2$sample, unique(drivers.qc.2$sample)),
  gene = match(drivers.qc.2$gene, unique(drivers.qc.2$gene)),
  y = drivers.qc.2$timing,
  nX = 3,
  X = as.matrix(data.frame(
    clonality = ifelse(drivers.qc.2$clonality_status=="clonal", 0, 1),
    subtype = ifelse(drivers.qc.2$subtype=="Intrahepatic", 0, 1),
    sex = ifelse(drivers.qc.2$sex=="F",0, 1)))
)


#Fit beta regression 1
bfit1 <- bmod$sample(
  data = beta_data_1,
  seed = 999,
  chains =  4,
  parallel_chains = 4,
  threads_per_chain = 1,
  refresh = 10, # print update every 50 iters
  iter_warmup = 1600,
  iter_sampling = 250
)

print(bfit1$summary(variables = c("beta", "kappa", "alpha")), n=50)

alpha <- bfit1$draws(variables = c("alpha"), 
                    format = "matrix")
beta <- bfit1$draws(variables = c("beta"), 
                    format = "matrix")
beta_mean <- apply(beta, 2, mean)

alpha_CI <- as.data.frame(t(apply(plogis(alpha), 2, 
            FUN= function(x) quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95)))))
alpha_CI$gene <- unique(drivers.qc.2$gene)
alpha_sort <- alpha_CI[order(alpha_CI$`50%`),]
alpha_x <- 1:nrow(alpha_sort)

#Plot coefficients by gene
#pdf("Google Drive/My Drive/Opisthorchis/CCA-evolution/gene-beta-coefficient.pdf",
#    height=7, width=9)
par(mar=c(6,5,1,0), xpd=F)
plot(alpha_sort$`50%`~alpha_x, pch=16, cex=1.6, axes=F, 
     ylab="Chronological time", xlab="", cex.lab=1.4,
     ylim=c(0.32,1))
arrows(x0=alpha_x, x1=alpha_x,
       y0=alpha_sort$`50%`, y1=alpha_sort$`95%`,
       angle=90, length = 0.00)
arrows(x0=alpha_x, x1=alpha_x,
       y0=alpha_sort$`50%`, y1=alpha_sort$`5%`,
       angle=90, length = 0.00)
arrows(x0=alpha_x, x1=alpha_x, lwd=3, col="navy",
       y0=alpha_sort$`50%`, y1=alpha_sort$`75%`,
       angle=90, length = 0.00)
arrows(x0=alpha_x, x1=alpha_x, lwd=3, col="navy",
       y0=alpha_sort$`50%`, y1=alpha_sort$`25%`,
       angle=90, length = 0.00)
axis(2, las=2, lwd=1.2, cex.axis=1.2)
#abline(h=0.5, lty=2)
axis(1, at=alpha_x,labels=alpha_sort$gene, las=2, lwd=1.2)
text(c("Preferentially Early", "Preferentially Late"), 
     cex.main=1.2, x=c(5, 39), y=c(1,1), cex=1.2)
#dev.off()

#Latent period 
latent <- clon_df$age_surgery-clon_df$age_mut
summary(clon_df$age_mut)
summary(latent)

#plot beta coefficients
alpha_mu <- median(apply(alpha,2,mean))
beta_CI<- as.data.frame(t(apply(beta,2, 
      FUN=function(x) quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95)))))
beta_CI$var <- c("clonality:Subclonal", "subtype:Perihiler", "sex:Male")
beta_sort <- beta_CI[order(beta_CI$`50%`),]
beta_x <- c(0.1,0.2,0.3)

par(mar=c(10,5,1,0), xpd=F)
plot(beta_sort$`50%`~beta_x,
     pch=16, cex=1.6, axes=F, 
     ylab="Coefficient (logit)", xlab="", cex.lab=1.4,
     ylim=c(-0.2,0.6), xlim=c(0.1,0.3))
axis(2, las=2, lwd=1.2, cex.axis=1.2)
axis(1, at=beta_x,labels=beta_sort$var, las=2, lwd=0, cex.axis=1.5, line=-1.2)
arrows(x0=beta_x, x1=beta_x,
       y0=beta_sort$`50%`, y1=beta_sort$`95%`,
       angle=90, length = 0.00)
arrows(x0=beta_x, x1=beta_x,
       y0=beta_sort$`50%`, y1=beta_sort$`5%`,
       angle=90, length = 0.00)
arrows(x0=beta_x, x1=beta_x, lwd=3, col="navy",
       y0=beta_sort$`50%`, y1=beta_sort$`75%`,
       angle=90, length = 0.00)
arrows(x0=beta_x, x1=beta_x, lwd=3, col="navy",
       y0=beta_sort$`50%`, y1=beta_sort$`25%`,
       angle=90, length = 0.00)
abline(h=0, lty=2, lwd=1.2)

#Plot
clonal.indx <- which(drivers.qc$clonality_status=="clonal")
subclonal.indx <- which(drivers.qc$clonality_status=="subclonal")
sort(summary.factor(drivers.qc$chr))

#pdf("Google Drive/My Drive/Opisthorchis/CCA-evolution/chr-copy-num-gains-age.pdf",
#        height=6.2, width=7.8)
par(mar=c(5,5,2,3))
plot(1, type="n", ylim=c(0,1), axes=F, ylab="", xlab="", xlim=c(1,22))
abline(v=12, lty=3, lwd=1.2)
abline(v=17, lty=3, lwd=1.2)

points(drivers.qc$timing[clonal.indx]~drivers.qc$chr[clonal.indx],
       pch=16,  col=alpha("seagreen", alpha=0.85), cex=1.4)
points(drivers.qc$timing[clonal.indx]~drivers.qc$chr[clonal.indx], cex=1.4)
points(drivers.qc$timing[subclonal.indx]~drivers.qc$chr[subclonal.indx],
       pch=16, col=alpha("salmon2", alpha=0.85), cex=1.4)
points(drivers.qc$timing[subclonal.indx]~drivers.qc$chr[subclonal.indx], cex=1.4)
axis(1, cex.axis=1.2, lwd=1.2, gap.axis=0.2,
     at=seq(1,22,1), labels=seq(1,22,1))
axis(2, las=2, cex.axis=1.4, lwd=1.2)
title(xlab="Chromosome", cex.lab=1.4)
title(ylab="Chronological time", cex.lab=1.4)
legend("bottomleft", legend=c("Clonal", "Subclonal"), 
       col=c("seagreen", "salmon2"), 
       bty="n", cex=1.3, pch=16, pt.cex=1.8)
#dev.off()

#Amplification of driver genes by fraction of lifetime
pdf("Google Drive/My Drive/Opisthorchis/CCA-evolution/copy-num-gains-age.pdf",
    height=6, width=7)
par(mar=c(5,6,1,1))
hist(drivers.qc$timing, col=alpha("#E41A1C", alpha=0.75), 
     ylim=c(0,45), xlim=c(0,1), axes=F, xaxs="i", yaxs="i",
     main="", 
     xlab="", 
     ylab="", lwd=2.5, breaks=20)
hist(drivers.qc$timing[drivers.qc$clonality_status=="clonal"], 
     col=alpha("white",1), add=T, lwd=3, breaks=20)
hist(drivers.qc$timing[drivers.qc$clonality_status=="clonal"], 
     col=alpha("#4DAF4A",0.8), add=T, lwd=5, breaks=20)

axis(1, cex.axis=1.4, lwd=1.2)
axis(2, las=2, cex.axis=1.4, lwd=1.2, at=c(0, 0.05, 0.1, 0.15)*nrow(drivers.qc), 
     labels=c(0, 0.05, 0.1, 0.15))
title(ylab="Proportion of copy number gains", cex.lab=1.4, line=4)
title(xlab="Chronological time", cex.lab=1.4)
legend(x=c(0.05,0.2), y=c(40,45), legend=c("Clonal", "Subclonal"), 
       col=c(alpha("#4DAF4A",0.8), alpha("#E41A1C", alpha=0.8)), 
       bty="n", cex=1.4, lwd=7)
dev.off()

#Summary of induction period
summary(clon_df$age_mut-2.2)
quantile(clon_df$age_mut-2.2, probs=c(0.05,0.5,0.95))

#Rough latent period
summary(clon_df$age_at_surgery-clon_df$age_mut)

#Random seed for jitter
set.seed(199)
x_jitter <- rnorm(nrow(sub_df), mean = 1, sd = 0.12)

#######################################################
## Combined boxplot - clonal and subclonal mutations ##
#######################################################

clonal_ihep_idx <- which(clon_df$CCA=="Intrahepatic") 
clonal_peri_idx <- which(clon_df$CCA=="Perihilar") 
sub_ihep_idx <- which(sub_df$CCA=="Intrahepatic") 
sub_peri_idx <- which(sub_df$CCA=="Perihilar") 

#Make pdf
pdf("Google Drive/My Drive/Opisthorchis/CCA-evolution/clonal-sub-mutation-age-boxplot.pdf",
   height=7, width=9)

#x jitter values
set.seed(1996)
x_clonal <- runif(n=nrow(clon_df), min=0.05, max=0.25)
x_subclonal <- runif(n=nrow(sub_df), min=0.40, max=0.60)

par(mar=c(1,6,3.5,0), xpd=T)
boxplot(clon_df$age_mut, 
        axes=F, boxwex=0.40, lty=1, lwd=1.6,
        col=alpha("grey60", 0.6),
        ylab="", ylim=c(5,65), at=0.15, xlim=c(0, 0.65))

boxplot(sub_df$age_mut, 
        axes=F, lty=1, lwd=1.6,
        col=alpha("grey60", 0.6),
        boxwex=0.40,
        xlim=c(0, 0.65),
        ylab="", ylim=c(5,65),at=0.5, add=T)

axis(2, las=2, lwd=1.4, cex.axis=1.4)
title(ylab="Patient age, inferred (years)", cex.lab=1.8, line = 3.8)

#Intrahepatic
points(y=clon_df$age_mut[clonal_ihep_idx],
       x=x_clonal[clonal_ihep_idx], 
       pch=16, cex=1.7, col=alpha("skyblue3", 0.8))

points(y=clon_df$age_mut[clonal_ihep_idx],
       x=x_clonal[clonal_ihep_idx], 
       pch=1, cex=1.7)

points(y=sub_df$age_mut[sub_ihep_idx],
       x=x_subclonal[sub_ihep_idx], 
       pch=16, cex=1.7, col=alpha("skyblue3", 0.8))

points(y=sub_df$age_mut[sub_ihep_idx],
       x=x_subclonal[sub_ihep_idx], 
       pch=1, cex=1.7)

#Perihilar
points(y=clon_df$age_mut[clonal_peri_idx],
       x=x_clonal[clonal_peri_idx], 
       pch=16, cex=1.7, col=alpha("seagreen", 0.8))

points(y=clon_df$age_mut[clonal_peri_idx],
       x=x_clonal[clonal_peri_idx], 
       pch=1, cex=1.7)

points(y=sub_df$age_mut[sub_peri_idx],
       x=x_subclonal[sub_peri_idx], 
       pch=16, cex=1.7, col=alpha("seagreen", 0.8))

points(y=sub_df$age_mut[sub_peri_idx],
       x=x_subclonal[sub_peri_idx], 
       pch=1, cex=1.7)

legend(x=c(0.24, 0.45), y=c(60,65),
       title="Tumour subtype", 
       legend=c("Intrahepatic", "Perihiler"), 
       pch=16, col=c(alpha("seagreen",0.9), alpha("skyblue",0.9)),
  bty="n", cex=1.4, pt.cex=1.8)

legend(x=c(0.24, 0.45), y=c(60,65),
       title="", 
       legend=c("",""), 
       pch=1, col=c("black", "black"),
       bty="n", cex=1.4, pt.cex=1.85)

title("Earliest amplified CCA driver mutations", cex.main=1.8, line=1.3)
text(c("Clonal", "Subclonal"), 
     cex.main=1.5, x=c(0.15, 0.5), y=c(66, 66), cex=1.5)

clonal_offset <- c(2.2, 2.2, 2.0, -2.2, 2.2, 2.2, 2.2, -2.2, -2.2, 2.2, 2.2, 2.2, -2.2, -2.2, -2.2, -2.2, -2.2)
clonal_text_y <- clon_df$age_mut + clonal_offset
text(clon_df$gene, x=x_clonal, y=clonal_text_y, cex=0.95)

subclonal_offset <- c(2.2, 2.2, 2.2, -2.0, 2.0, -2.2, 2.2, -1.8, 1.8, 2.2, -2.0, -2.0, -2.2, -2.2)
sub_x_offset <- c(0,0,0,0,0,0,0,0,0.01,0,-0.01,0,0,0)
subclonal_text_y <- sub_df$age_mut + subclonal_offset
subclonal_text_x <- x_subclonal + sub_x_offset
text(sub_df$gene, x=subclonal_text_x, y=subclonal_text_y, cex=0.95)

dev.off()

driver_chr17 <- drivers.qc[drivers.qc$chr==17,]
y_row <- 1:nrow(driver_chr17)

par(mfrow=c(1,1), xpd=F, mar=c(5,5,2,2))
plot(NULL, ylim=c(1,nrow(driver_chr17)+1),
     xlim=c(0, chrs.length[17]), axes=F, ylab="", xlab="")
for(i in 1:nrow(driver_chr17))
  polygon(x=c(driver_chr17$start[i], driver_chr17$end[i]), y=c(y_row[i], y_row[i]+0.5))

#birth age (guess)
summary(2016 - supp$Age.at.surgery) #all born before 1980

