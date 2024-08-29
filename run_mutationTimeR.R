## Run MutationTimeR script on all 22 TH CCA samples

#Relevant packages
list.of.packages <- c(
  "foreach", "doParallel","ranger","palmerpenguins",
  "tidyverse", "kableExtra", "MutationTimeR", "sitools",
  "AmplificationTimeR", "BSgenome.Hsapiens.UCSC.hg19")

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE)
  )
}

#Extra Functions
averagePloidy <- function(bb) {
  c <- if(!is.null(bb$copy_number)) bb$copy_number else bb$total_cn
  sum(width(bb) * c * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Get ploidy per chromosome
ploidyChrom <- function(bb) {
  ploidy.chr <- c()
  for(i in 1:length(chrs.length)){
    bb.chr <- subsetByOverlaps(bb, chrs.GRanges[[i]])
    c <- if(!is.null(bb.chr$copy_number)) bb.chr$copy_number else bb.chr$total_cn
    ploidy.chr[i] <- sum(width(bb.chr) * c * bb.chr$clonal_frequency, na.rm=TRUE) / sum(width(bb.chr) * bb.chr$clonal_frequency, na.rm=TRUE)
  }
  return(ploidy.chr)
}

averageHom <- function(bb){
  sum(width(bb) * (bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

.classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy
classWgd <- function(bb) .classWgd(averagePloidy(bb), averageHom(bb))

fractionGenomeWgdCompatible <- function(bb, min.dist=0.05){
  m <- findMainCluster(bb)
  l <- pmin(bb$time.lo, bb$time - min.dist)
  u <- pmax(bb$time.up, bb$time + min.dist)
  w <- which(l <= m & u >= m)
  avgCi <- weighted.mean(bb$time.up- bb$time.lo, width(bb), na.rm=TRUE)
  sd.wgd <- sqrt(weighted.mean((bb$time[w] - m)^2, width(bb)[w], na.rm=TRUE))
  sd.all <- sqrt(weighted.mean((bb$time - m)^2, width(bb), na.rm=TRUE))
  c(nt.wgd=sum(as.numeric(width(bb))[w]), nt.total=sum(as.numeric(width(bb))[!is.na(bb$time)]), time.wgd=m, n.wgd=length(w), n.all = sum(!is.na(bb$time)), chr.wgd = length(unique(seqnames(bb)[w])), chr.all = length(unique(seqnames(bb)[!is.na(bb$time)])), sd.wgd=sd.wgd, avg.ci=avgCi, sd.all=sd.all) 
}

findMainCluster <- function(bb, min.dist=0.05){
  w <- which(bb$n.snv_mnv > 20 & !is.na(bb$time))
  s <- seq(0,1,0.01)
  l2 <- pmin(bb$time.lo, bb$time - min.dist)[w]
  u2 <- pmax(bb$time.up, bb$time + min.dist)[w]
  l1 <- (l2 +  bb$time[w])/2
  u1 <- (u2+  bb$time[w])/2
  wd <- as.numeric(width(bb)[w])
  o <- sapply(s, function(i) sum(wd * ( (l2 <= i & u2 >=i) + (l1 <= i & u1 >= i))))
  s[which.max(o)]
}

#Set cores
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()

#Meta data file with cellularity
meta <- read.delim("Desktop/MutationTimeR/master.dpclust.txt")
supp <- read.delim("Desktop/MutationTimeR/jusakul-TH.txt", sep="\t")

#Driver mutations
drivers.all <- read.delim("Desktop/MutationTimeR/drivers/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv")
load("Desktop/MutationTimeR/drivers/TableS2_driver_point_mutations_annotation.RData")

#Lists
bb <- muts <- mult <- cl.list <- clusters <- vcfs <- subclones <- list()

#Ploidy by chromosome
chrs.name <- c((1:22), "X")
chrs.length <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23]
chrs.GRanges <- list()

for(i in 1:length(chrs.name)){
  chrs.GRanges[[i]] <- GRanges(seqnames=chrs.name[i],
                               ranges=IRanges(start = 1, end = chrs.length[i]),
                               strand="*")
}

#Read in files
for(i in 1:nrow(meta)){
  cl.list[[i]] <- read.delim(paste0("Desktop/MutationTimeR/clusters/TH_",i,"_2000iters_1000burnin_bestClusterInfo.txt"))
  vcfs[[i]] <- readVcf(paste0("Desktop/MutationTimeR/vcf/TH_",i,".somatic.pass.vcf.gz"))
  subclones[[i]] <- read.delim(paste0("Desktop/MutationTimeR/subclones/T_CCA_TH_",i,"_subclones.txt"))
  muts[[i]] <- read.table(paste0("Desktop/MutationTimeR/muts/TH_",i,".muts.txt"))
  colnames(muts[[i]]) <- c("chr", "start", "end","ref", "alt")
  mult[[i]] <- read.delim(paste0("Desktop/MutationTimeR/dpclust3p/T_CCA_TH_",i,".dpclust3p.txt"))
}

#Update cluster df to get cellular fraction
for(i in 1:nrow(meta)){
  clusters[[i]] <- data.frame(
                  cluster=cl.list[[i]]$cluster.no, 
                  n_ssms=cl.list[[i]]$no.of.mutations, 
                  proportion=cl.list[[i]]$location*meta$cellularity[i])

  #Create bb object 
  bb[[i]] <- GRanges(seqnames=subclones[[i]]$chr, 
                ranges=paste0(subclones[[i]]$startpos,"-", subclones[[i]]$endpos), 
                major_cn=subclones[[i]]$nMaj1_A, minor_cn=subclones[[i]]$nMin1_A, 
                clonal_frequency=meta$cellularity[i])
  
  bb[[i]]$total_cn <- bb[[i]]$major_cn+bb[[i]]$minor_cn
  
}

#Convert columns on vcf
for(i in 1:nrow(meta)){
  #Add t_alt_count
  t_ref_count <- do.call(rbind.data.frame, geno(vcfs[[i]])$AD[,2])[,1]
  t_alt_count <- do.call(rbind.data.frame, geno(vcfs[[i]])$AD[,2])[,2]
  
  info(vcfs[[i]])$t_ref_count <- t_ref_count
  info(vcfs[[i]])$t_alt_count <- t_alt_count
}

#Run mutation timer - in parallel
#Takes about 10 mins with n.boot=100
mt <- foreach(i = 1:22,
              .packages = c("MutationTimeR", "VariantAnnotation", "S4Vectors")
              ) %dopar% {
        mutationTime(vcf=vcfs[[i]], 
                cn=bb[[i]], 
                clusters=clusters[[i]], 
                n.boot=100, 
                gender=meta$sex[i])
}

#stop the cluster
parallel::stopCluster(cl = my.cluster)

#Assess MutationTimeR output
for(i in 1:22){
  table(mt[[i]]$V$CLS)
  vcfs[[i]] <- addMutTime(vcfs[[i]], mt[[i]]$V)
  mcols(bb[[i]]) <- cbind(mcols(bb[[i]]),mt[[i]]$T)
}
#Assess MutationTimeR output
q1 <- sapply(1:22, FUN=function(x) {1-mean(abs(0.5-info(vcfs[[x]])$pMutCNTail)>0.495 , na.rm=TRUE)}) #This should be higher than 0.95
q5 <- sapply(1:22, FUN=function(x) {1-mean(abs(0.5-info(vcfs[[x]])$pMutCNTail)>0.475 , na.rm=TRUE)}) #This should be higher than 0.95
par(mfrow=c(1,1))
boxplot(q5, ylim=c(0.6,1), lty=1, 
        ylab="Fraction of data inside theoretical 95% CI")
abline(h=0.95, lty=3, lwd=1.2)

#Plot VAF and timings
#for(i in 1:22){
i <- 19

#Mutation timer is not very exact within chromosomes
for(i in 1:22){
  try(plotSample(vcfs[[i]], bb[[i]], regions =chrs.GRanges[[17]],
             title=paste0("TH",i)))
}


#pdf("Desktop/TH19-MutationTimerR.pdf", width=10, height=7.8)
#plotSample_hack(vcfs[[i]],bb[[i]])
  plotSample(vcfs[[i]],bb[[i]])
             #title=paste0("TH",i, " | Tumour purity = ", 
              #            round(meta$cellularity[i],digits = 2), 
               #           " | q5 = ", round(q5[i], digits = 2),
                #          " | Patient age = ", supp$Age.at.surgery[i],
                 #         " | Subtype = ", supp$Anatomical.subtype[i])
             
#dev.off()
  #}

#QQ plots
par(mfrow=c(5,5), mar=c(3,3,3,1))
for(i in 1:22){
  n <- nrow(vcfs[[i]])
  qqnorm(qnorm(info(vcfs[[i]])$pMutCNTail[sample(1:n, min(1e4,n))]), 
       main=paste0("TH_",i, " Q5 = ", signif(q5[i],2), ", Q1 = ", signif(q1[i],2)), 
       xlim=c(-5,5), ylim=c(-5,5), pch=16)
  abline(0,1, col='red')
}

boxplot(data=supp, Age.at.surgery~as.factor(supp$Sex))

#Check WGD status
finalPloidy <- sapply(bb, averagePloidy)
finalHom <- sapply(bb, averageHom)
isWgd <- .classWgd(finalPloidy, finalHom)
table(isWgd)

#ad mt columns to bb
for(i in 1:22){
  bb[[i]]$type <- mt[[i]]$T$type
  bb[[i]]$time <- mt[[i]]$T$time
  bb[[i]]$time.lo <- mt[[i]]$T$time.lo
  bb[[i]]$time.up <- mt[[i]]$T$time.up
  bb[[i]]$time.star <- mt[[i]]$T$time.star
  bb[[i]]$n.snv_mnv <- mt[[i]]$T$n.snv_mnv
}

fracGenomeWgdComp <- t(sapply(bb, function(bb) {
  fgw <- try(fractionGenomeWgdCompatible(bb)); 
  if(class(fgw)!='try-error') fgw
  else rep(NA,10)}))

chrOffset <- cumsum(c(0,as.numeric(seqlengths(BSgenome.Hsapiens.UCSC.hg19))))

wgdStar <- factor(rep(1,nrow(fracGenomeWgdComp)), levels=0:3, 
                  labels=c("unlikely","uninformative","likely","very likely"))

wgdStar[fracGenomeWgdComp[,"avg.ci"]<=0.75 & fracGenomeWgdComp[,"nt.total"]/chrOffset[25] >= 0.33 ] <- "likely"
wgdStar[fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"] < 0.66] <- "unlikely"
wgdStar[wgdStar=="likely" & fracGenomeWgdComp[,"nt.wgd"]/fracGenomeWgdComp[,"nt.total"]>0.8 & fracGenomeWgdComp[,"sd.wgd"]<0.1 & fracGenomeWgdComp[,"nt.total"]/chrOffset[25] > 0.5] <- "very likely"
prop.table(table(wgdStar[!isWgd]))
wgdPoss <- !isWgd & 2.5 - 1.5 * finalHom <= finalPloidy
wgdStat <- factor(wgdPoss + 2*isWgd - wgdPoss*isWgd, labels=c("absent","present"))
table(wgdStat, wgdStar)
which(wgdStar %in% c("very likely"))

#Temporal distribution of chromosomal gains
aggregatePerChromosome <- function(bb, isWgd=FALSE){
.aggregateSegments <- function(m){
  #m <- mcols(bb)
  t <- weighted.mean(m$time, m$n.snv_mnv, na.rm=TRUE)
  n <- sum(m$n.snv_mnv[!is.na(m$time)], na.rm=TRUE)
  sd <- sd(m$time, na.rm=TRUE)
  ci <- weighted.mean(m$time.up-m$time.lo, m$n.snv_mnv, na.rm=TRUE)
  w <- sum(m$width[!is.na(m$time)], na.rm=TRUE)
  c(time=t, n=n, sd=sd, ci=ci,w=w)
}
#   if(!isWgd){
s <- split(as.data.frame(bb)[,c("time","time.up","time.lo","n.snv_mnv","width")], seqnames(bb))
r <- t(sapply(s, .aggregateSegments))
r <- r[c(1:22,"X"),]
#   }else{
w <- .aggregateSegments(as.data.frame(bb))
r <- rbind(r,WGD=w)
#   }
return(r)
}

allChrAgg <- simplify2array(mclapply(bb, aggregatePerChromosome, mc.cores=2))
t <- t(allChrAgg[1:23,"time",])
t[allChrAgg[1:23,"w",] < diff(chrOffset)[1:23]*.33] <- NA

par(mfrow=c(6,4))
for(i in 1:23){
  if(all(is.na(t[,i]))){
    plot(NA,NA,xlab="",ylab="", xlim=c(0,1),ylim=c(0,6), bty="n")
  }else{
  try(hist(na.omit(t[,i]), xlim=c(0,1), ylim=c(0,6)))
  }
}
par(mfrow=c(6,4))
for(i in 1:22){
  if(all(is.na(t[i,]))){
    plot(NA,NA,xlab="",ylab="", xlim=c(0,1),ylim=c(0,6), bty="n")
  }else{
    try(hist(na.omit(t[i,]), xlim=c(0,1), ylim=c(0,6)))
  }
}

par(mfrow=c(6,4))
for(i in 1:22){
  if(all(is.na(bb[[i]]$time))){
    plot(NA,NA,xlab="",ylab="", xlim=c(0,1),ylim=c(0,6), bty="n")
  }else{
    try(hist(na.omit(bb[[i]]$time), xlim=c(0,1)))
  }
}


all_times <- c()
earliest_gain <- c()
for(i in 1:22){
  all_times <- c(all_times, na.omit(bb[[i]]$time))
  earliest_gain <- c(earliest_gain, min(bb[[i]]$time, na.rm=T))
}
summary(all_times)
summary(earliest_gain[!is.infinite(earliest_gain)])
summary()

n <- 10
at <- function(x, n){
  if(sum(!is.na(x))<3) return(rep(sum(!is.na(x))/n,n))
  bw=if(sum(!is.na(x))< 6) 0.5 else "nrd0"
  d <- density(x, n=n, from=1/n/2, to=1-1/n/2, bw=bw, na.rm=TRUE)
  d$y/sum(d$y)*d$n
}

allChrCancerHist <- sapply(t(t), at, n=n, simplify="array")
u <- data.frame(WGD=allChrAgg["WGD","time",isWgd])
wgdCancerHist <- at(u$WGD,n=n)
allChrCancerHist <- abind::abind(allChrCancerHist, 
                                 All=sapply(sapply(t(t), as.matrix), 
                                 at, n=n, simplify="array")/23*5, 
                                 WGD=wgdCancerHist, along=2)

#Crazy plot
prgn <- RColorBrewer::brewer.pal(11,"PRGn")
set1 <- RColorBrewer::brewer.pal(9,"Set1")
col <- colorRampPalette(set1[c(4,9,3)])(n)

p <- 0
v <- table(droplevels(donor2type[sample2donor[names(finalSnv)]]))
h <- (allChrCancerHist + p)  / rep(v + p, each=prod(dim(allChrCancerHist)[1:2]))
h <- aperm(h, c(2,3,1))


chr.ploid.mat <- sapply(bb, ploidyChrom)

apply(chr.ploid.mat, 1, median)

par(mfrow=c(4,6))
for(i in 1:length(chrs.length)){
  boxplot(chr.ploid.mat[i,], ylim=c(0,5), lty=1)
  title(main=chrs.name[i])
  title(ylab="Ploidy", cex.lab=1.2, line=2.7)
  abline(h=2, lty=2)
}

#Check ploidy by sample and identify the no-hopers
par(mfrow=c(4,6))
for(i in 1:22){
  boxplot(chr.ploid.mat[,i], ylim=c(0,5), lty=1)
  title(main=paste0("TH",i))
  title(ylab="Ploidy", cex.lab=1.2, line=2.7)
  abline(h=2, lty=2)
}

tumours <- list()
tumour_depth <- c()
for(i in 1:22){
  t.path <- paste("Google Drive/My Drive/Opisthorchis/CCA-evolution/CCA_TH_depth/single_sample/T_CCA_TH_", i, ".bin.1Mb", sep="")
  tumours[[i]] <- read.table(t.path, header=T)
  tumour_depth[i] <- round(median(tumours[[i]]$mean))
}

#Power calculation
power <- function(coverage, purity, ploidy){
  pw = coverage * (purity/((ploidy*purity)+(1-purity)*2))
  return(pw)
}

power.tumour <- c()
for(i in 1:22){
  power.tumour[i] <- power(tumour_depth[i], purity=bb[[i]]$clonal_frequency[1], ploidy=finalPloidy[i])
}

#Attempt to merge vcf TH9 with annotated drivers
# Get timings from amplificationTimer
findOverlaps(vcfs[[1]], finalDrivers)
finalDrivers[c(2106, 2132, 2146, 684)]
cn.idx <- which(!is.na(mt[[1]]$T$time))
bb[[1]][cn.idx]
findOverlaps(bb[[1]], finalDrivers[c(2106,684)])
bb[[1]][c(18,83)]
mt[[1]]$V[c(1427, 5521),]
info(vcfs[[1]])[c(1427, 5521),]
summary.factor(mt[[1]]$T$time)
