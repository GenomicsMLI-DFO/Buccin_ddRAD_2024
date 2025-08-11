# Info --------------------------------------------------------------------

# Filtering pipeline (including STACKS version >2 population)
# ALL SAMPLES > 10X (to get a second dataset)
#
# Audrey Bourret
# 2024-06-12
#

# Library -----------------------------------------------------------------
library(parallel) # to detect the number of cores
library(tidyverse)
#library(readxl)

library(vcfR)
library(adegenet)
library(hierfstat)

library(here)

`%nin%` = Negate(`%in%`)

# Internal functions
for(i in 1:length( list.files("./01_Codes/Functions") )){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
}

# Data --------------------------------------------------------------------

stacks.ref.path <- "00_Data/05a_Stacks.denovo.buccinum.M1n1/"
filter.ref.path <-  "00_Data/06a_Filtering.denovo.buccinum.M1n1/"

pop.data <- read_csv(file.path(get.value("info.path"),"Project_Infos_20240220.csv"))

pop.data 

pop.data <- pop.data %>% mutate(Notes_groupe = ifelse(Numero_unique_groupe == "G_23_00984","Groupe_EGSL_Baie_Comeau1_Offshore",
                                                      ifelse(Numero_unique_groupe == "G_23_00989","Groupe_EGSL_Baie_Comeau2_Offshore",       
                                                             Notes_groupe)),
                                Shore = ifelse(str_detect(Notes_groupe, "Offshore"), "Offshore",
                                               ifelse(str_detect(Notes_groupe, "Inshore"), "Inshore", "Problems"                 
                                               )),
                                Site = Notes_groupe %>% str_remove_all("Groupe_|_Offshore|_Inshore|Newfoundland_|EGSL_|Baie_de_Fundy_"),
                                Region = Notes_groupe %>% str_remove("Groupe_") %>% str_remove("_Offshore|_Inshore") %>% str_remove(paste(paste0("_", Site),  collapse = "|" ))
)
pop.data %>% filter(Cat_sample == "Sample") %>% 
group_by(Espece, Region, Shore, Site ) %>% summarise(N = n())

# What is the sample size
pop.data %>% pull(Cat_sample) %>% table()

pop.data %>% filter(Cat_sample == "Sample") %>% 
  group_by(Region_echantillonnage) %>% 
  summarise(N = n()) #%>% write_csv("clipboard")

# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

numCores <- if(get_os() %in% c("os","linux")){
  detectCores() # Utilise le max de coeurs si sur linux
} else 1

numCores

# Coverage ----------------------------------------------------------------

cov.data <- read_csv(file.path("02_Results/00_Stacks/05a_Stacks.denovo/AllIndividuals_DeNovo.buccinum.M1n1_NreadsNloci.csv"))

cov.data 

cov.data %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), fill = Shore)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 10, lty = "dashed", col = "darkgray") +
  facet_wrap(~Region, nrow = 2) +
  labs(x = "Mean coverage") + 
  theme_bw() +
  theme(legend.position = "bottom")

cov.data %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Espece))) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = c(8,10), lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = c(90000,100000), lty = "dashed", col = "darkgray") +
  #  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "none")

cov.data %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = as.factor(Region), col = as.factor(Espece))) +
  geom_boxplot(col = "black") +  
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 10, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")


cov.data %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov)), y = as.numeric(as.character(mean_cov_ns)), col = Region)) +
  geom_point() +  
  #geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  #geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")


cov.data %>% 
  ggplot(aes(x = as.numeric(as.character(n_loci)), y = as.numeric(as.character(mean_cov_ns)), col = Region)) +
  geom_point() +  
  #geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  #geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")


men.cov <- cov.data$mean_cov_ns

quantile(men.cov,probs=c(0.01, 0.05, 0.95, 0.99))

cov.data %>% dplyr::filter(mean_cov_ns < 8) %>% nrow()
cov.data %>% dplyr::filter(mean_cov_ns > 65) %>% nrow() 

cov.data %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), fill = Region)) +
  geom_histogram(binwidth = 5) +  
  #geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = quantile(men.cov, probs=c(0.01, 0.05, 0.95, 0.99)), lty = "dashed", col = "darkgray") +
  geom_vline(xintercept = c(10, 100),  col = "red") +
  
  
  
  #geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")


# List of species to use for next steps 
# I will use all buccinum samples (no duplicate), and also I will perform 
# some stats by Gen_ZONE
ID.coverage <- cov.data %>% filter(n_loci >=90000) %>% 
  mutate(Pop = "NoPop") %>% 
  select(sample, Pop) 

nrow(cov.data)
nrow(ID.coverage)

write.table(ID.coverage, 
            file = file.path(filter.ref.path, "A_90K", "popmap.coverage_90K.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

#read.table(file.path(get.value("info.path"), "popmap.coverage.txt")) %>% View()

# Filtering --------------------------------------------

# Parameters0
r.value       <- 0.75 # Minimum within pop
maf.value     <- 0.01 # Overall MAF
#maf.pop.value <- 0.05

# Filtering #1 : -r and -MAF ---------------------------------------------
if(!file.exists(file.path(filter.ref.path, "A_90K", "01_r75_MAF01"))){
  dir.create(file.path(filter.ref.path, "A_90K", "01_r75_MAF01"))
  print(file.path(filter.ref.path,"A_90K", "01_r75_MAF01"))
}

cmd <- paste("-P", stacks.ref.path, 
             "-M", file.path(filter.ref.path,"A_90K", "popmap.coverage_90K.txt"),
             "--out-path",  file.path(filter.ref.path, "A_90K",  "01_r75_MAF01"),
             "-t", 8,
             "-r", r.value, #              
             "--min-maf", maf.value,
             #"--smooth",
             #"--write-single-snp",
             "--vcf"
)

A <- system2("populations", cmd, stdout=T, stderr=T)
A

# If you want the big file to be ignored, run the following :

cat("*", "!.gitignore", "*.log", sep = "\n",
    file = file.path(filter.ref.path,  "A_90K",  "01_r75_MAF01", ".gitignore"))


# General check : Ho and Fis, He vs Ho outliers ---------------------------------

# Load the post-filtration statistic table

filter.stat <- read.delim(file.path(filter.ref.path, "A_90K",  "01_r75_MAF01", "populations.sumstats.tsv"), 
                          skip=1, sep = "\t", header = T )
names(filter.stat)[1] <- "Locus.ID"

nrow(filter.stat) 

summary(filter.stat)
head(filter.stat)

# For fun, the distribution of the number of SNPs, by scaffold and locus

filter.stat %>% group_by(Locus.ID) %>% summarise(Nsnps = n()) %>% 
  group_by(Nsnps) %>% summarise(Nloc = n()) %>% View()

filter.stat %>% group_by(Chr) %>% summarise(Nloc = length(Locus.ID %>% unique())  ) %>% 
  group_by(Nloc) %>% summarise(Nscaffold = n()) %>% View()

# Check the distribution of Fis and Ho
#scaffold84112,15249,f87087Z15249

filter.stat %>% ggplot(aes(x = Obs.Het, fill=Pop.ID)) +
  geom_histogram() +
  theme_bw()

filter.stat %>% ggplot(aes(x = Fis)) +
  geom_histogram() +
  theme_bw()


filter.stat %>% ggplot(aes(x = Obs.Het, y = Fis, col = Pop.ID)) +
  #geom_point(alpha = 1/100) +
  geom_point()+
  geom_vline(xintercept = 0.5, col = "red", lty = "dashed")+
  geom_hline(yintercept = c(-0.3), col = "red", lty = "dashed") +
  theme_bw()

filter.stat %>% ggplot(aes(x = Obs.Het, y = Exp.Het, col = Fis)) +
  geom_point() +
  geom_vline(xintercept = 0.5, col = "red", lty = "dashed")+
  scale_colour_gradientn(colours=rainbow(4)) +
  geom_abline(slope = 1 ) +
  #geom_hline(yintercept = c(-0.3,0.3), col = "red", lty = "dashed") +
  theme_bw()

filter.stat %>% ggplot(aes(x = Obs.Het, y = Exp.Het)) +
  geom_point(alpha = 0.02) +
  geom_vline(xintercept = c(0.5, 0.6), col = "red", lty = "dashed")+
  scale_colour_gradientn(colours=rainbow(4)) +
  geom_abline(slope = c(1,2), col = "blue")  +
  #geom_hline(yintercept = c(-0.3,0.3), col = "red", lty = "dashed") +
  theme_bw()

# Filtration #2 : Missing data -----------------------------------------

# # Then more filtering with VCF tools

# Verif : individual and SNPs with > 30% missing data

list.files(file.path(filter.ref.path,  "A_90K",  "01_r75_MAF01"))

vcf.path <- file.path(filter.ref.path,  "A_90K",  "01_r75_MAF01", "populations.snps.vcf")

cmd1 <- paste("--vcf", file.path(current.wd, vcf.path), 
              #"--out",  "MAX2.NArm",
              #" --max-alleles", 2,
              #"--max-missing", 0.80,
              #"--missing-site",
              "--missing-indv"
              #"--kept-sites"
              #"--recode"
)

# cmd2 <- paste("--vcf", file.path(current.wd, vcf.path), 
#               #"--out",  "MAX2.NArm",
#               #"--max-alleles", 2,
#               #"--max-missing", 0.70,
#               "--missing-site"
#               #"--missing-indv",
#               #"--kept-sites"
#               #"--recode"
# )


setwd(file.path(filter.ref.path,  "A_90K",  "02_MissingData")) 

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)
A1

cat(file = "populations.filt_buccinum_90K_Individuals_wMissing_A.log",
    "\n", cmd1, "\n",
    A1, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

# Many ind with more 30% missing value

imiss <- read.delim(file.path(filter.ref.path,  "A_90K",  "02_MissingData", "out.imiss"), skip=0, sep = "\t", header = T )
imiss %>% head()
imiss %>% filter(F_MISS >.30)

imiss %>% left_join(pop.data %>% select(INDV = ID_GQ, Region)) %>%  ggplot(aes(x=F_MISS, fill = Region)) + geom_histogram()


graph2.0 <- imiss %>% left_join(pop.data, by = c("INDV" = "ID_GQ")) %>% 
  ggplot(aes(x = Site, y = F_MISS, col =  F_MISS)) + 
  geom_boxplot(col = "blue") +
  geom_jitter(height = 0, alpha = 0.75) +
  geom_hline(yintercept = 0.3, lty = "dashed", col = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
graph2.0

ggsave(filename = file.path(filter.ref.path, "A_90K",  "02_MissingData", "NAbyInd_90K_June2024.png"),
       plot = graph2.0,
       width = 6, height = 4, units = "in")

# From there, what I will do it's to remove duplicated for the SNPs selection process,
# and add them back later when necessary.

imiss %>% left_join(cov.data, by = c("INDV" = "sample")) %>% 
  #left_join(pop.data, by = c("INDV" = "ID_GQ")) %>% 
  ggplot(aes(y = F_MISS, x =  mean_cov_ns, col = Region)) + 
  geom_point()#+
  #facet_grid(~Region)

test.ID <- imiss %>% #left_join(cov.data, by = c("INDV" = "sample")) %>% 
  left_join(pop.data, by = c("INDV" = "ID_GQ")) %>% 
  arrange(F_MISS) %>% 
  mutate(Dup_post_MISSGING = duplicated(Numero_unique_specimen)) 


test.ID %>% dplyr::filter(Numero_unique_specimen %in% test.ID$Numero_unique_specimen[test.ID$Dup_post_MISSGING == T]) %>% 
  
  
  ggplot(aes(x = Numero_unique_specimen, y = F_MISS, col =  Dup_post_MISSGING)) + 
  #geom_boxplot(col = "blue") +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


good.ID <- test.ID  %>% filter(Dup_post_MISSGING == F,
                               F_MISS <=0.3) %>% pull(INDV) %>% as.character()


pop.data %>% filter(ID_GQ %in% good.ID) %>% #select(ID_GQ, Gen_ZONE) %>% 
  group_by(Region, Site) %>% summarise(N = n())


cmd2 <- paste("--vcf", vcf.path, 
              "--recode",
              paste("--indv", good.ID, collapse = " "),
              "--out",  vcf.path %>% str_replace("01_r75_MAF01", "02_MissingData") %>% 
                str_replace("populations.snps.vcf", paste0("populations.snps.", length(good.ID), "ind.vcf"))
)

cmd2

A2 <- system2("vcftools", cmd2, stdout=T, stderr=T)
A2 %>% tail()

cat(file = file.path(filter.ref.path,  "A_90K",  "02_MissingData","VCFtools_RemoveIndMissing.log"),
    "\n", cmd2, "\n",
    A2, # what to put in my file
    append= F, sep = "\n")


list.files(file.path(filter.ref.path,  "A_90K",  "02_MissingData"))

vcf.path <- file.path(filter.ref.path,  "A_90K",  "02_MissingData", "populations.snps.431ind.vcf.recode.vcf")

file.exists(vcf.path)

cmd3 <- paste("--vcf", file.path(current.wd,  vcf.path ), 
              #"--out",  "MAX2.NArm",
              #"--max-alleles", 2,
              #"--max-missing", 0.70,
              "--missing-site"
              #"--missing-indv",
              #"--kept-sites"
              #"--recode"
)

cmd3

setwd(file.path(filter.ref.path,  "A_90K",  "02_MissingData")) 

A3 <- system2("vcftools", cmd3, stdout=T, stderr=T)
A3

cat(file = "ppopulations.filt_buccinum_90K_Loci_wMissing_B.log",
    "\n", cmd3, "\n",
    A3, # what to put in my file
    append= F, sep = "\n")


setwd(current.wd)

# Change the initial 25% missing to 10% missing

lmiss <- read.delim(file.path(filter.ref.path, "A_90K", "02_MissingData", "out.lmiss"), header = T, sep = "\t" )
head(lmiss %>% arrange(desc(F_MISS)))

graph2.1 <- lmiss %>% ggplot(aes(x=F_MISS)) + 
  geom_histogram() + 
  geom_vline(xintercept = 0.10, lty = "dashed", col = "red") +
  theme_bw()


graph2.1  

ggsave(filename = file.path(filter.ref.path,  "A_90K", "02_MissingData", "NAbySNPs_10x_June2024.png"),
       plot = graph2.1,
       width = 4, height = 4, units = "in")


good.LOC <- lmiss %>% filter(F_MISS <= .10) %>% select(CHR,POS)

nrow(good.LOC)

setwd(current.wd)


cmd4 <- paste("--vcf", file.path(filter.ref.path,  "A_90K", "02_MissingData", "populations.snps.431ind.vcf.recode.vcf"), 
              "--recode",
              "--max-missing", "0.9",
              #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
              
              #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
              #"--snps", file.path(current.wd,get.value("filter.ref.path"),"A_GLOBAL_PB_5x", "06c_Indiduals_wMissing","LOC.10.csv"),
              
              "--out", file.path(current.wd, filter.ref.path, "A_90K", "02_MissingData", paste0("populations.", nrow(good.LOC),"snps", ".431ind"))
)

cmd4


A4 <- system2("vcftools", cmd4, stdout=T, stderr=T)
A4

cat(file = file.path(filter.ref.path,  "A_90K", "02_MissingData","VCFtools+RemoveLociMissing_Loci.log"),
    "\n", cmd4, "\n",
    A4, # what to put in my file
    append= F, sep = "\n")

# If you want the big file to be ignored, run the following :

cat("*.vcf", "!.gitignore", sep = "\n",
    file = file.path(filter.ref.path,  "A_90K", "02_MissingData", ".gitignore"))


# Filtration #3 : HW outliers ---------------------------------------

# This part can be long ... (data convertion to hierfstat)

list.files(file.path(filter.ref.path,"A_90K", "02_MissingData"))

# Load the VCF file

vcf.data <- vcfR::read.vcfR(file.path(filter.ref.path,  "A_90K", "02_MissingData", "populations.100393snps.431ind.recode.vcf"))

head(vcf.data)

vcf_field_names(vcf.data , tag = "FORMAT")

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1),
         scaffold.length = sapply(str_split(CHROM, ","), `[`,2))

SCAFFOLD.info$RADloc %>% unique() %>% length()

# Conversion to various R formats
#gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# Computing diversity (Ho, He)

div <- summary(gi.data)

str(div)
 

div.graph <- data.frame(ID = names(div$Hobs),
                        Hobs = div$Hobs,
                        Hexp = div$Hexp)


#div.graph <- div.graph %>% left_join(HW.loc) %>% 
#                           mutate(Npop = ifelse(is.na(Npop), 0, Npop))

div.graph %>% head()

#save(list = c("div", "div.graph"),
#     file = file.path(filter.ref.path, "A_90K", "03_HW", "SNPs_90K_HW.RES.RData"))

load(file.path(filter.ref.path, "A_90K", "03_HW", "SNPs_90K_HW.RES.RData"))

graph3.0 <-  div.graph %>% ggplot(aes(x = Hobs, y = Hexp)) +
  geom_point(alpha = 0.1) +
  scale_colour_distiller(palette = "Spectral") +
  geom_vline(xintercept = 0.6) +
  labs(title = "Ho vs He overall")+ 
  theme_bw()
graph3.0


ggsave(filename = file.path(filter.ref.path,"A_90K", "03_HW", "He_Ho_June2024.png"),
       plot = graph3.0,
       width = 4, height = 4, units = "in")


# Add info om loci

div.graph <- div.graph %>% mutate(Loc = sapply(str_split(ID, ":"), `[`, 1)) 

Nspns.tab <- div.graph %>% group_by(Loc) %>% 
  summarise(N = n())                    

# Keep 1 snp / rad loci
#list.h05.unique <- div.graph %>% filter(Hobs <=0.5)  %>% distinct(Loc, .keep_all = T) 
list.h06.unique <- div.graph %>% filter(Hobs <=0.6)  %>% distinct(Loc, .keep_all = T) 
#list.n13.unique <- div.graph %>% filter(Npop <=13)  %>% distinct(Loc, .keep_all = T) 

# Keep ALL
#list.h05.all  <- div.graph %>% filter(Hobs <=0.5) 
list.h06.all  <- div.graph %>% filter(Hobs <=0.6) 
#list.n13.all <- div.graph %>% filter(Npop <=13) 

#nrow(list.h05.unique)
#nrow(list.h05.all)
nrow(list.h06.unique)
nrow(list.h06.all)

#write.csv(list.h05.unique %>% select(ID), file.path(filter.ref.path,"06c_HW", "Loc.h05.unique.csv"), 
#          row.names = F, quote = F)

#write.csv(list.h05.all %>% select(ID), file.path(filter.ref.path,"06c_HW", "Loc.h05.all.csv"), 
#          row.names = F, quote = F)


write.csv(list.h06.unique %>% select(ID), file.path(filter.ref.path,"A_90K", "03_HW", "Loc.h06.unique.csv"), 
          row.names = F, quote = F)

write.csv(list.h06.all %>% select(ID), file.path(filter.ref.path, "A_90K", "03_HW","Loc.h06.all.csv"), 
          row.names = F, quote = F)

#setwd(current.wd)

# CREATE VCF WITH UNIQUE


cmd3 <- paste("--vcf", file.path(filter.ref.path, "A_90K/02_MissingData/populations.100393snps.431ind.recode.vcf"), 
              "--recode",
              #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
              
              #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
              "--snps", file.path(here::here(),filter.ref.path,"A_90K/03_HW", "Loc.h06.unique.csv"),
              
              "--out", file.path(here::here(),filter.ref.path,"A_90K/03_HW", paste0("populations.",list.h06.unique %>% nrow(),"snps.431ind.H06.single"))
)

cmd3

A3 <- system2("vcftools", cmd3, stdout=T, stderr=T)
A3

cat(file = file.path(here::here(),filter.ref.path,"A_90K/03_HW","VCFtools_SnpWhiteList.H06.log"),
    "\n", cmd3, "\n",
    A3, # what to put in my file
    append= F, sep = "\n")


cmd4 <- paste("--vcf", file.path(filter.ref.path, "A_90K/02_MissingData/populations.100393snps.431ind.recode.vcf"), 
              "--recode",
              #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
              
              #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
              "--snps", file.path(here::here(),filter.ref.path,"A_90K/03_HW", "Loc.h06.all.csv"),
              "--out", file.path(here::here(),filter.ref.path,"A_90K/03_HW", paste0("populations.",list.h06.all %>% nrow(),"snps.431ind.h06.all"))
)

cmd4

A4 <- system2("vcftools", cmd4, stdout=T, stderr=T)
A4

cat(file = file.path(here::here(),filter.ref.path,"A_90K/03_HW","VCFtools_SnpWhiteList.H06.log"),
    "\n", cmd4, "\n",
    A4, # what to put in my file
    append= T, sep = "\n")


# Add a gitignore if necessary

if(!file.exists(file.path(filter.ref.path,"A_90K/03_HW", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.RData", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "A_90K/03_HW", ".gitignore")) 
}

# rapid PCA
file.path(current.wd,filter.ref.path,"A_90K/03_HW") %>% list.files()

vcf.path <- file.path(filter.ref.path,"A_90K/03_HW", "populations.47673snps.431ind.H06.single.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

# Function to create a list of loci, from a genind object

filter.MAF.NA <- function(gi, MAF.trs = 0.5, NA.trs = 0.5){
  # Create vectors for each loci
  MAF.res <- adegenet::minorAllele(gi)
  NA.res  <- na.gi.count(gi)
  
  # Filter by threshold
  MAF.loc <- dimnames(MAF.res[MAF.res >= MAF.trs])[[1]]
  cat("There is", length( MAF.loc), "loci with MAF =", MAF.trs, "\n")
  
  NA.loc <- names(NA.res[NA.res <= NA.trs])
  cat("There is", length(NA.loc), "loci with NA =", NA.trs, "\n")
  
  # LOCI with both conditions
  LOCI.res <- c(MAF.loc, NA.loc)[duplicated(c(MAF.loc, NA.loc)) == T]
  LOCI.res %>% length()
  
  cat("There is", length(LOCI.res), "loci with BOTH MAF =", MAF.trs, "and NA =" , NA.trs, "\n")
  
  return(LOCI.res)
}

LOC.MAF10.NA05 <- filter.MAF.NA(gi.data, MAF.trs = 0.10, NA.trs = 0.05)

pca.test  <- glPca(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05], center = TRUE, scale = FALSE,  
                   parallel = TRUE, n.core =8, nf = 1000)

#save(list = c("pca.test", "gl.data", "gi.data" ),
#     file = here::here(filter.ref.path, "A_90K", "03_HW", "PCA.Rdata"))


load(  here(filter.ref.path, "A_90K", "03_HW", "PCA.Rdata"))

QuickPop::pca_varplot(pca.test)

gPCA <- pca.test %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, shape = Shore)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Site) +
#  stat_ellipse()+
 
  ggforce::geom_mark_ellipse(aes(fill = Region, label = Site, group = Site)) +
   geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw() #+ theme(legend.position = "none")
gPCA


#weird.ID.PCA <- pca.test %>% QuickPop::pca_scoretable(naxe = 5) %>%
#  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
#  dplyr::filter(score.PC1 < -10) %>% pull(ID)

ggsave(filename = file.path(filter.ref.path, "A_90K", "03_HW", "PCA_test.png"),
       plot = gPCA,
       width = 5, height = 4, units = "in")

# Filtration #5 : Too low and too high depth ----------------------------------------------------------

vcf.path <-  file.path(here::here(), filter.ref.path,  "A_90K","03_HW", "populations.98582snps.431ind.h06.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

# Check stats

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
   select(ID, CHROM, POS) %>% 
   mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
          RADloc = sapply(str_split(ID, ":"), `[`,1),
          scaffold.length = sapply(str_split(CHROM, ","), `[`,2))

SCAFFOLD.info %>% head()

length(unique(SCAFFOLD.info$RADloc))
length(unique(SCAFFOLD.info$ID))

# Remove loci with too much depth

gt.tidy <- extract_gt_tidy(vcf.data, format_types = NULL, format_fields = c("DP"))
gt.tidy <- gt.tidy %>% mutate(gt_DP = as.numeric(as.character(gt_DP)))
head(gt.tidy)


vcf.fix <- as.data.frame(vcf.data@fix) %>% mutate(Key = 1:nrow(.))


# Group by individuals


gt.ind <-  gt.tidy %>% group_by(Indiv) %>% summarise(N = n(),
                                                     minDP = min(gt_DP, na.rm = T),
                                                     medianDP = median(gt_DP, na.rm = T),
                                                     Nmin3 = length(Key[gt_DP <3]),
                                                     Nmin5 = length(Key[gt_DP <5])
) 

gt.ind %>% head()

graphDP0 <- gt.ind %>% 
  ggplot(aes(x = Nmin3, y = Nmin5, col = medianDP)) +
  geom_jitter(alpha = 0.5)+
  geom_vline(xintercept = quantile(gt.ind$Nmin3,  c(0.01, 0.5, 0.99)), lty = "dashed") +
  scale_color_viridis_c() +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()
graphDP0


gt.ind %>% left_join(pop.data, by = c("Indiv" = "ID_GQ")) %>% 
  ggplot(aes(x = Site, y = Nmin3, col = medianDP)) +
  geom_boxplot() +
  geom_jitter(height = 0, alpha = 0.5) +
  scale_color_viridis_c() 

gt.ind %>% left_join(pop.data, by = c("Indiv" = "ID_GQ")) %>% 
  ggplot(aes(x = No_plaque_envoi, y = Nmin3, col = medianDP)) +
  geom_boxplot() +
  geom_jitter(height = 0, alpha = 0.5) +
  scale_color_viridis_c() +
  facet_grid(.~No_soumission_GQ, scale = "free") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# group by snps
gt.key <- gt.tidy %>% group_by(Key) %>% summarise(medianDP = median(gt_DP, na.rm = T),
                                                  meanDP = mean(gt_DP, na.rm = T),
                                                  sdDP = sd(gt_DP, na.rm = T),
                                                  maxDP = max(gt_DP, na.rm = T),
                                                  minDP = min(gt_DP, na.rm = T)) %>% 
  left_join(vcf.fix %>% select(Key, ID))

graphDP1 <- gt.key %>% 
  ggplot(aes(x = medianDP, y = meanDP, col = maxDP)) +
  geom_jitter(alpha = 0.5)+
  geom_vline(xintercept = quantile(gt.key$medianDP,  c(0.01, 0.5, 0.99)), lty = "dashed") +
  scale_color_viridis_c() +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()
graphDP1

graphDP2 <- gt.key %>% 
  ggplot(aes(x = medianDP)) +
  geom_histogram(bins = 50)+
  geom_vline(xintercept = quantile(gt.key$medianDP, c(0.01, 0.5, 0.99)), lty = "dashed") +
  scale_y_continuous(trans = "log10") +
  labs(title = "Some SNPs with way too high DP")+
  theme_bw()
graphDP2

quantile(gt.key$medianDP, c(0.001,0.01, 0.5,.99, 0.999))

mean(gt.key$medianDP) - (2* sd(gt.key$medianDP))
mean(gt.key$medianDP) + (2* sd(gt.key$medianDP))

# Final (2 times SD)

LOC.DP.min <- gt.key %>% filter(medianDP < 10) %>% pull(ID)
LOC.DP.max <- gt.key %>% filter(medianDP > 50) %>% pull(ID)

length(LOC.DP.min)
length(LOC.DP.max)

LOC.H06 <- read.csv(file.path(filter.ref.path, "A_90K", "03_HW", "Loc.h06.all.csv"))

LOC.H06.DP <- LOC.H06 %>% filter(ID %nin% c(LOC.DP.min, LOC.DP.max))

nrow(LOC.H06.DP) + length(LOC.DP.min) + length(LOC.DP.max)== nrow(LOC.H06)

write.csv(LOC.H06.DP, file.path(filter.ref.path, "A_90K", "04_DP", "Loc.H06.all.DP.csv"), 
          row.names = F, quote = F)


# Do the filtration

cmd5 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              #paste("--snp", LOC.FILTER$ID[1:2000], collapse = " "),
              
              #"--positions", file.path(get.value("filter.ref.path"),"06c_HeHo_byPlate","Test.tsv"),
              "--snps", file.path(current.wd,filter.ref.path,"A_90K", "04_DP", "Loc.H06.all.DP.csv"),
              "--out", file.path(current.wd,filter.ref.path,"A_90K", "04_DP", paste0("populations.",LOC.H06.DP %>% nrow(),"snps.431ind.H06.all"))
)

cmd5

A5 <- system2("vcftools", cmd5, stdout=T, stderr=T)
A5

cat(file = file.path(current.wd,filter.ref.path,"A_90K", "04_DP","VCFtools_SnpWhiteList.H06.DP.log"),
    "\n", cmd5, "\n",
    A5, # what to put in my file
    append= F, sep = "\n")


if(!file.exists(file.path(filter.ref.path,"A_90K", "04_DP", ".gitignore")) ){
  cat("*.vcf", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "A_90K", "04_DP", ".gitignore")) 
}

# Filtration - Batch effect -----------------------------------------------

file.path(current.wd,filter.ref.path, "A_90K", "04_DP") %>% list.files()

vcf.path <- file.path(filter.ref.path,"A_90K", "04_DP", "populations.97819snps.431ind.H06.all.recode.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

# Stats

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1),
         scaffold.length = sapply(str_split(CHROM, ","), `[`,2))

SCAFFOLD.info %>% head()

length(unique(SCAFFOLD.info$RADloc))
length(unique(SCAFFOLD.info$ID))

gl.data  <- vcfR::vcfR2genlight(vcf.data) 

rm(list = c("vcf.data"))

# RDA on plate

pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% left_join(pop.data) %>% pull(No_plaque_envoi)

pop(gl.data) %>% table()

library(vegan)

## Full model - maybe not the best implementation of NA values (mean)
RDA.plate <- vegan::rda(tab(gl.data) ~ pop(gl.data))

# tab(gl.data.sex)[1:5,1:5]

save(list = c("gl.data", "RDA.plate"),
     file = here(filter.ref.path, "A_90K", "05_BatchEffect", "rda.Rdata"))


load(here::here(filter.ref.path,"A_90K" ,"05_BatchEffect", "rda.Rdata"))

# Group by individuals
anova(RDA.plate)

# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = tab(gl.data) ~ pop(gl.data))
# Df Variance      F Pr(>F)  
# Model      5    162.6 1.0801  0.012 *
#   Residual 425  12798.0                
# ---

RsquareAdj(RDA.plate)$adj.r.squared  # 0.0009302119

res.plate <- as.data.frame(scores(RDA.plate, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDA.plate, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

res.plate.corrected <- as.data.frame(scores(RDA.plate.corrected, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDA.plate.corrected, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

screeplot(RDA.plate)

gg.plate.rda <- res.plate %>% #mutate(Region_test = ifelse(Region_echantillonnage %in% c("JAM", "SLE", "SAN"), Region_echantillonnage, "Other")) %>% 
  ggplot(aes(x = RDA1, y = RDA2, group = No_plaque_envoi)) +
  #facet_wrap(~No_soumission_GQ) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  geom_point(aes( col = No_plaque_envoi)) +
  stat_ellipse()+
  ggtitle("Before correction") +
  #geom_point(data = res2, colour = "grey70", cex = 1) +
  #geom_histogram() +
  #geom_density()+
  labs(x=c("RDA 1")) +  
  theme_bw() #+
#theme(legend.position = "none")
gg.plate.rda


# Source Capblancq

rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

library(robust)
library(qvalue)

# The one chosen was k2 thres 0.05, but I also test k4 thres 10 and the improvement was not high

rdadapt_plate <- rdadapt(RDA.plate, 2)

## P-values threshold after Bonferroni correction
thres_plate <- 0.05/length(rdadapt_plate$p.values)

## Identifying the loci that are below the p-value threshold
outliers_plate <- data.frame(Loci = colnames(tab(gl.data))[which(rdadapt_plate$p.values<thres_plate)], 
                             p.value = rdadapt_plate$p.values[which(rdadapt_plate$p.values<thres_plate)])
outliers_plate  %>% head()  
outliers_plate  %>% nrow()

## List of outlier names

## Formatting table for ggplot
locus_scores <- scores(RDA.plate, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores) %>% 
  mutate(type = ifelse(names %in% outliers_plate$Loci, "Outlier", "Non outlier"))

TAB_var <- as.data.frame(scores(RDA.plate, choices=c(1,2), display="bp")) %>% 
  mutate(ID = row.names(scores(RDA.plate, choices=c(1,2), display="bp")) ) #%>% left_join(final.env.names)
TAB_var  # pull the biplot scores

## Biplot of RDA loci and variables scores
gbiplot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1, alpha = 1/2) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  # geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5) +
  #ggrepel::geom_label_repel(data = TAB_var, aes(label = ID, x=1.1*RDA1, y=1.1*RDA2), size = 2, max.overlaps = 20
  #) +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11) +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.8)),
        strip.background = element_rect(fill = "white"))
gbiplot


ggsave(filename = file.path(filter.ref.path,"A_90K", "05_BatchEffect", "Plate_biplot_20240612.png"),  # radloci different order than Fis.MeanDepthPerSNP.139317snps.png (graph0.1)
       plot = gbiplot,
       width = 9, height = 7, units = "in", bg = "white")

# Filtration
length(outliers_plate$Loci)

LOC.H06.DP <- read.csv(file.path(filter.ref.path,"A_90K","04_DP", "Loc.H06.all.DP.csv"))

LOC.H06.PLATE <- LOC.H06.HW %>% filter(ID %nin% outliers_plate$Loci)

nrow(LOC.H06.PLATE) + length(outliers_plate$Loci) == nrow(LOC.H06.DP)

write.csv(LOC.H06.PLATE, file.path(filter.ref.path,"A_90K", "05_BatchEffect", "Loc.H06.all.PLATE.csv"), 
          row.names = F, quote = F)

vcf.path

cmd6 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              "--snps", file.path(filter.ref.path,"A_90K", "05_BatchEffect", "Loc.H06.all.PLATE.csv"),
              "--out", file.path(filter.ref.path,"A_90K", "05_BatchEffect", paste0("populations.",LOC.H06.PLATE %>% nrow(),"snps.431ind.H06.all"))  # understand why snps in name are wrongly estimated
)
cmd6

A6 <- system2("vcftools", cmd6, stdout=T, stderr=T)
A6

cat(file = file.path(filter.ref.path,"A_90K", "05_BatchEffect","VCFtools_SnpWhiteList.H06.PLATE.log"),
    "\n", cmd6, "\n",
    A6, # what to put in my file
    append= F, sep = "\n")


if(!file.exists(file.path(filter.ref.path, "A_90K", "05_BatchEffect", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.RData", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "A_90K", "05_BatchEffect", ".gitignore")) 
}


# Filtration #7: Sex-linked SNPs ----------------------------------------

list.files(file.path(filter.ref.path, "A_90K", "05_BatchEffect"))

vcf.data <- vcfR::read.vcfR(file.path(filter.ref.path, "A_90K", "05_BatchEffect", "populations.96095snps.431ind.H06.all.recode.vcf"))

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  # how many RAD loci left after removing sex-linked loci
  select(ID, CHROM, POS) %>% 
  mutate(chromosome = sapply(str_split(CHROM, "_"), `[`,4),
         scaffold = sapply(str_split(CHROM, "_"), `[`,6), #%>% str_remove_all("scaffold|contig"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

head(SCAFFOLD.info)
SCAFFOLD.info %>% group_by(chromosome) %>% summarise(N = n())

length(unique(SCAFFOLD.info$RADloc))  # 46828
length(SCAFFOLD.info$RADloc) # 96095

# RDA on sex
gl.data <- vcfR::vcfR2genlight(vcf.data)

pop.data %>% pull(Sexe_visuel) %>% table(useNA = "ifany")

ID.withSex <- pop.data %>% filter(Sexe_visuel %in% c("F","M")) %>% pull(ID_GQ)
ID.withSex %>% length()

gl.data.sex <- gl.data[indNames(gl.data) %in% ID.withSex, ]
pop(gl.data.sex) <- data.frame(ID_GQ = indNames(gl.data.sex)) %>% left_join(pop.data) %>% pull(Sexe_visuel)

pop(gl.data.sex) %>% table()

library(vegan)

## Full model
RDA.sex <- vegan::rda(tab(gl.data.sex) ~ pop(gl.data.sex))
# tab(gl.data.sex)[1:5,1:5]

anova(RDA.sex)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = tab(gl.data.sex) ~ pop(gl.data.sex))
# Df Variance      F Pr(>F)   
# Model      1     44.2 1.1283  0.002 **
#   Residual 298  11683.5                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                 

RsquareAdj(RDA.sex)$adj.r.squared   # 0.0004289733

res.sex <- as.data.frame(scores(RDA.sex, display="sites", scaling=1)) %>% 
  mutate(ID_GQ = dimnames(scores(RDA.sex, display="sites", scaling=1))[[1]]) %>% 
  left_join(pop.data)

# save(list = c("RDA.sex", "res.sex"),
#      file = file.path(filter.ref.path,"A_90K", "06_Sex", "SNPs_SexLinked.RData"))
load(file.path(filter.ref.path, "A_90K", "06_Sex", "SNPs_SexLinked.RData"))

res.sex %>% ggplot(aes(x = RDA1, fill = Sexe_visuel)) +
  #facet_wrap(~RegionAssesment) +
  geom_vline(xintercept = 0) +   geom_hline(yintercept = 0) +
  #geom_point(data = res2, colour = "grey70", cex = 1) +
  geom_histogram() +
  #geom_density()+
  labs(x=c("RDA 1")) +  theme_bw()  

hist(RDA.sex$CCA$v)

sd(RDA.sex$CCA$v)*4

sex.loc <- dimnames(RDA.sex$CCA$v)[[1]][abs(RDA.sex$CCA$v) > 0.013]
sex.loc

length(sex.loc)


LOC.H06.PLATE <- read.csv(file.path(filter.ref.path,"A_90K", "05_BatchEffect", "Loc.H06.all.PLATE.csv")) 

LOC.H06.SEX <- LOC.H06.PLATE %>% filter(ID %nin% c(sex.loc))

nrow(LOC.H06.SEX) + length(sex.loc) == nrow(LOC.H06.PLATE)

write.csv(LOC.H06.SEX, file.path(filter.ref.path, "A_90K", "06_Sex", "Loc.H06.all.SEX.csv"), 
          row.names = F, quote = F)

list.files(file.path(filter.ref.path, "A_90K", "05_BatchEffect"))

vcf.path <-  file.path(file.path(filter.ref.path, "A_90K", "05_BatchEffect", "populations.96095snps.431ind.H06.all.recode.vcf"))

cmd6 <- paste("--vcf", file.path(vcf.path), 
              "--recode",
              "--snps", file.path(filter.ref.path, "A_90K", "06_Sex", "Loc.H06.all.SEX.csv"),
              "--out", file.path(filter.ref.path, "A_90K", "06_Sex", paste0("populations.",LOC.H06.SEX %>% nrow(),"snps.431ind.H06.all"))  # understand why snps in name are wrongly estimated
)
cmd6

A6 <- system2("vcftools", cmd6, stdout=T, stderr=T)
A6

cat(file = file.path(current.wd,filter.ref.path,"A_90K", "06_Sex","VCFtools_SnpWhiteList.H06.SEX.log"),
    "\n", cmd6, "\n",
    A6, # what to put in my file
    append= F, sep = "\n")

if(!file.exists(file.path(filter.ref.path, "A_90K", "06_Sex", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "*.RData", "*.data", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path, "A_90K", "06_Sex", ".gitignore")) 
}


# General Check #8 : Relatedness ----------------------------------------------

current.wd <- here::here()

file.path(filter.ref.path, "A_90K", "06_Sex") %>% list.files()

# In VCF tool, check for duplicated ind.

cmd <- paste("--vcf", file.path(here::here(), filter.ref.path, "A_90K", "06_Sex", "populations.95718snps.431ind.H06.all.recode.vcf"), 
             #sub_indv,
             # --depth,
             "--relatedness2"
             
)


setwd(file.path(filter.ref.path,"A_90K", "07_Relatedness")) 

A <- system2("vcftools", cmd, stdout=T, stderr=T)
A

cat(file = file.path("VCFtools.Relatedness.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

setwd(current.wd)

# Add a gitignore if necessary

if(!file.exists(file.path(filter.ref.path,"A_90K", "07_Relatedness", ".gitignore")) ){
  cat("out.relatedness2", "!.gitignore", sep = "\n",
      file = file.path(filter.ref.path,"A_90K", "07_Relatedness",  ".gitignore")) 
}

related2.all <- read.delim(file.path(filter.ref.path,"A_90K", "07_Relatedness","out.relatedness2"), header = T, sep = "\t" )
related2.all %>% head()

length(related2.all$INDV1 %>% unique)

# Remove t=identical individuals
related2.ok <- related2.all %>% mutate(iden = ifelse(INDV1 == INDV2, 1, 0)) %>% # Remove identique
  filter(iden == 0) %>% 
  left_join(pop.data, by = c("INDV1" = "ID_GQ"))

related2.ok %>% filter(RELATEDNESS_PHI >=1/8) %>% arrange(desc(RELATEDNESS_PHI)) %>% View()

graph.relatedness <- related2.ok %>% 
   ggplot(aes(x =  RELATEDNESS_PHI)) +
  labs(x = "Relatedness coefficient", y = "N observations", title = "Relatedness coefficient distribution")+
  geom_histogram(bins = 1000) +
  geom_vline(xintercept = c(1/2, 1/4, 1/8, 1/16), lty = "dashed")+
  #facet_wrap(~No_plaque) +
  annotate("text", 
           x = c(1/2-0.005, 1/4-0.005, 1/8-0.005, 1/16-0.005), y = 1000, 
           label = c("Individual-self", "Siblings / Parent-offspring", "Half-siblings / Grandparent-grandchild", "First cousins"), 
           angle = 90, hjust = 0, vjust = 0) +
  theme_bw()  + #
  ggforce::facet_zoom(xlim = c(0.125, 0.5), ylim = c(0, 5), zoom.size = 1, horizontal = FALSE)
#  theme(axis.text.x = element_text(angle = 60, hjust = 1))

graph.relatedness 

ggsave(filename = file.path(filter.ref.path,"A_90K", "07_Relatedness", "Relatedness.png"),
       plot = graph.relatedness ,
       width = 7, height = 7, units = "in")


# Final dataset -----------------------------------------------------------

vcf.path <- file.path(filter.ref.path, "A_90K/06_Sex/populations.95718snps.431ind.H06.all.recode.vcf")

file.exists(vcf.path)

vcf.data <- vcfR::read.vcfR(vcf.path)

#gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# Transformer le tout en beau jeu de données :
SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1)
  )

SCAFFOLD.info

SCAFFOLD.info$RADloc %>% unique() %>% length()

na.gi.count <- function(gi){
  res <- apply(tab(gi), MARGIN = 2, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  res <- res[str_ends(names(res), "[.]0")] 
  
  names(res) <- names(res) %>% str_remove("[.]0")
  
  return(res)
  
}

# Check that the order was kept
table(locNames(gi.data) == SCAFFOLD.info$ID )

# ADD NA and MAF info to select the "BEST" SNP / RAD loc

SCAFFOLD.info$NNA <- na.gi.count(gi.data)
SCAFFOLD.info$MAF <- adegenet::minorAllele(gi.data)


SCAFFOLD.info %>% group_by(RADloc) %>% summarise(N = n()) %>% arrange(desc(N)) %>% head()

SCAFFOLD.info$MAF %>% hist()

#SCAFFOLD.info.final <- SCAFFOLD.info %>% 
#  dplyr::filter(between(MAF, 0.01, 0.5)) %>% 
#  mutate(RADloc = as.numeric(RADloc)) %>% 
#  arrange(RADloc, desc(round(MAF, 1)), round(NNA,2)) %>% 
#  distinct(RADloc, .keep_all = T)

SCAFFOLD.info.MAF01 <-  SCAFFOLD.info %>% 
  dplyr::filter(MAF >= 0.01, NNA <=0.10) %>% 
  mutate(RADloc = as.numeric(RADloc)) %>% 
  arrange(RADloc, POS) %>% #head()
  distinct(RADloc, .keep_all = T) # %>% nrow()

SCAFFOLD.info.MAF05 <-  SCAFFOLD.info %>% 
  dplyr::filter(MAF >= 0.05, NNA <=0.10) %>% 
  mutate(RADloc = as.numeric(RADloc)) %>% 
  arrange(RADloc, POS) %>% #head()
  distinct(RADloc, .keep_all = T) # %>% nrow()

SCAFFOLD.info$MAF %>% hist()
SCAFFOLD.info.MAF01$MAF %>% hist()
SCAFFOLD.info.MAF05$MAF %>% hist()


write.csv(SCAFFOLD.info.MAF01%>% select(ID), file.path(filter.ref.path,"A_90K", "08_Final","Loc.H06.DP.unique.final.MAF01.csv"), 
          row.names = F, quote = F)

write.csv(SCAFFOLD.info.MAF01, file.path(filter.ref.path,"A_90K", "08_Final", "Scaffold.info.H06.DP.unique.final.MAF01.csv"), 
          row.names = F, quote = F)

write.csv(SCAFFOLD.info.MAF05%>% select(ID), file.path(filter.ref.path,"A_90K", "08_Final", "Loc.H06.DP.unique.final.MAF05.csv"), 
          row.names = F, quote = F)

write.csv(SCAFFOLD.info.MAF05, file.path(filter.ref.path,"A_90K", "08_Final", "Scaffold.info.H06.DP.unique.final.MAF05.csv"), 
          row.names = F, quote = F)

SCAFFOLD.info.MAF01 <- read.csv(file.path(filter.ref.path,"A_90K", "08_Final", "Scaffold.info.H06.DP.unique.final.MAF01.csv"))
SCAFFOLD.info.MAF01 %>% nrow()

SCAFFOLD.info.MAF05 <- read.csv(file.path(filter.ref.path,"A_90K", "08_Final",  "Scaffold.info.H06.DP.unique.final.MAF05.csv"))
SCAFFOLD.info.MAF05 %>% nrow()


# CREATE VCF WITH UNIQUE

vcf.path <- "./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/06_Sex/populations.95718snps.431ind.H06.all.recode.vcf"

cmd <- paste("--vcf", vcf.path, 
             "--recode",
             "--snps", file.path(here::here(),filter.ref.path, "A_90K", "08_Final", "Loc.H06.DP.unique.final.MAF05.csv"),
             
             "--out", file.path(here::here(),filter.ref.path, "A_90K", "08_Final", paste0("populations.", SCAFFOLD.info.MAF05 %>% nrow(),"snps.","431ind.single.final"))
)

cmd

A <- system2("vcftools", cmd, stdout=T, stderr=T)

tail(A)

cat(file =  file.path(filter.ref.path, "A_90K", "08_Final", "VCF.filter.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= T, sep = "\n")

# # Reload, and save as Rdata
# 
# vcf.path <- file.path( "./00_Data/06b_Filtering.ref",  "C_10X_samples_2024", "07_Final", "populations.53384snps.522indwArctogadus.H06.DP.single.final.recode.vcf")
# vcf.data <- vcfR::read.vcfR(vcf.path)
# 
# gl.final  <- vcfR::vcfR2genlight(vcf.data) 
# gi.final  <- vcfR::vcfR2genind(vcf.data) 
# 
# save(list = c("gl.final", "gi.final"),
#      file = file.path("./00_Data/06b_Filtering.ref", "C_10X_samples_2024", "07_Final", "populations.53384snps.522indwArctogadus.adegenet.Rdata"))
# 


if(!file.exists(file.path("./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final", ".gitignore")) ){
  cat("*.vcf", "*.Rdata", "!.gitignore", sep = "\n",
      file = file.path("./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final", ".gitignore")) 
}

