
getwd()

library(devtools)
 library(ggplot2)
geno <- read.table("geno.txt", header=FALSE)
dim(geno)
head(geno)
names(geno) <- c("chr", "pos", "ref", "alt", paste0("plant", 1:20))

###Create an F2 population
### Just sample 100 markers from this Mt chr
set.seed(12579)
markers <- sample(1:nrow(geno), size=100)
f <- geno[sort(markers), c("chr", "pos", "ref", "alt", "plant1", "plant18")]
# select just one haplotype
f$plant1 <- gsub("/.*", "", f$plant1)
f$plant18 <- gsub("/.*", "", f$plant18)
# recoding to use -1,0,1
f[f==0] <- -1
# simulate the recombination rate
f$cM <- f$pos/5000
usethis::edit_r_environ()
gh::gh_whoami()
r
#install.packages("devtools")
devtools::install_github("lian0090/simuPoisson")
library(simuPoisson)
#Then, simulate an F2 population with 200 individuals and 100 SNP markers.
pgeno <- t(f[, c("plant1", "plant18")])
pgeno <- apply(pgeno, 2, as.numeric)
map <- f[, c("chr", "pos", "cM")]
f2 <- simuPoisson(pgeno, map$chr, map$cM, 200)
f2 <- as.data.frame(f2)
names(f2) <- paste0(f$chr, "_", f$pos)
#f2=data.frame(ID=paste0("F2_",1:200),f2)
write.table(f2, "f2_geno.csv", sep=",", quote=FALSE)

f2 <- read.csv("f2_geno.csv", header=TRUE)
table(f2[,2])
table(f2[,5])

##Observed allele frequency
# For A1 allele
p <- (53*2+102)/((53+102+45)*2)
# For A2 allele
q <- (45*2+102)/((53+102+45)*2)

##Observed genotype frequency
# For A1A1 genotype
A1A1 <- 53/(53+102+45)
# For A1A2 genotype
A1A2 <- 102/(53+102+45)
# For A2A2 genotype
A2A2 <- 45/(53+102+45)

p^2
2*p*q
q^2
chisq.test(rbind(c(A1A1, A1A2, A2A2), c(p^2, 2*p*q, q^2)))

###tools
#vcftools --vcf TEO_LR_MZ_test2.vcf --freq --out TEO_LR_MZ_test2
#vcftools --vcf TEO_LR_MZ_test2.vcf --counts --out TEO_LR_MZ_test2
#plink --vcf TEO_LR_MZ_test2.vcf --hardy --out TEO_LR_MZ_test2
#TASSEL

####Phenotype data
pheno <- read.csv("f2_pheno.csv")
hist(pheno$height, main="Plant Height", xlab="Value (inch)", breaks=20)

#Combine genotype and phenotype files
gp <- cbind(pheno, f2)

##Let's find out and at a specific Marker Mt_29145:
ggplot(gp, aes(x=as.factor(Mt_29145), y=height, color=as.factor(Mt_29145))) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#fe6f5e"))+
  labs(title="Mt_29145", y="Plant Height", x = "Genotype")+
  theme_classic() +
  guides(color=FALSE) +
  theme(plot.title = element_text(size=20, face = "bold"),
        axis.text=element_text(size=16, face="bold"),
        strip.text.y = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=18, face="bold"),
  )

u <- mean(gp$height) # population mean
# A1A1
h1 <- mean(subset(gp, Mt_29145 == -1)$height)
# A1A2
h12 <- mean(subset(gp, Mt_29145 == 0)$height)
# A2A2
h2 <- mean(subset(gp, Mt_29145 == 1)$height)

a <- (h2 - h1)/2
midpoint <- h1+a
midpoint=(h2 + h1)/2
d <- h12 - midpoint

###The average effect of A1 and A2:
alpha <- a + d*(q - p)
alpha1 <- q*alpha
alpha2 <- -p*alpha

#Breeding value is the value of an individual as a parent!
bv1 = u+alpha1 + alpha1
bv2 = u+alpha2 + alpha2
bv12 = u+alpha1 + alpha2

#Genotypic value and breeding value
plot(c(0, 1, 2), c(h1, h12, h2), xlab="Genotype",ylab="", cex.lab=1.5, xaxt="n",pch=16, col="red",ylim=c(75,85))
     axis(1, at=c(0, 1, 2), labels=c("A1A1", "A1A2", "A2A2"));
     mtext("Breeding Value", side = 4, line = 1, cex=1.5, col="blue");
     mtext("Genotypic Value", side = 2, line = 2, cex=1.5, col="red")
     points(c(0, 1, 2), c(bv2, bv12, bv1), cex=2, col="blue")
     lines(c(0, 1, 2), c(bv2, bv12, bv1), lwd=2, col="blue")

###Additive and dominance variance
     Vp <- var(gp$height) #Phenotypic variance
     Va <- 2*p*q*(a + d*(q - p))^2 #Additive genetic variance
     Vd <- (2*p*q*d)^2 #Dominance genetic variance
     
     Vg <- Va + Vd
     
     h2 <- Va/Vp ##Narrow-sense heritability
     H2 <- Vg/Vp ##Broad-sense heritability
     
     