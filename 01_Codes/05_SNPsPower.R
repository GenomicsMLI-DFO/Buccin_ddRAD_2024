# Info --------------------------------------------------------------------

# Test how many markers needed to detect structure
#  431 individuals
#
# Audrey Bourret
# 2024-07-17
#

# Library -----------------------------------------------------------------

library(tidyverse)
library(ggpubr)

library(adegenet)

`%nin%` = Negate(`%in%`)

library(QuickPop)

# Data --------------------------------------------------------------------

pop.data <- read_csv(file.path("00_Data/00_FileInfos/" ,"Project_Infos_20240220.csv"))

pop.data 

pop.data <- pop.data %>% mutate(Notes_groupe = ifelse(Numero_unique_groupe == "G_23_00984","Groupe_EGSL_Baie_Comeau1_Offshore",
                                                      ifelse(Numero_unique_groupe == "G_23_00989","Groupe_EGSL_Baie_Comeau2_Offshore",       
                                                             Notes_groupe)),
                                Shore = ifelse(str_detect(Notes_groupe, "Offshore"), "Offshore",
                                               ifelse(str_detect(Notes_groupe, "Inshore"), "Inshore", "Problems"                 
                                               )),
                                Site = Notes_groupe %>% str_remove_all("Groupe_|_Offshore|_Inshore|Newfoundland_|EGSL_|Baie_de_Fundy_"),
                                Region = Notes_groupe %>% str_remove("Groupe_") %>% str_remove("_Offshore|_Inshore") %>% str_remove(paste(paste0("_", Site),  collapse = "|" )),
                                RegionGen = ifelse(Site %in% c("Bic1", "Bic2"), "EGSL-2a", 
                                                   ifelse(Site %in% c("Iles_de_la_Madeleine", "Mingan"), "EGSL-2b",
                                                          ifelse(Region %in% c("EGSL"), "EGSL-1", Region)
                                                   ))
)


pop.data %>% filter(Cat_sample == "Sample") %>% 
  mutate(Geno = ifelse(ID_GQ %in% indNames(gl.data) %>% str_remove("_rep"),
                       T,F)) %>% 
  group_by(Espece, Region, Shore, Site ) %>% summarise(N = n(),
                                                       Ngenotyped = length(ID_GQ[Geno == T]))


other.SP <- c("S_23_01903", "S_23_01910", "S_23_01913", "S_23_01922")

pop.data %>% dplyr::filter(ID_GQ %in% other.SP) %>% View()

# Genetic Data ------------------------------------------------------------

vcf.path <- file.path("./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.recode.vcf")

vcf.data <- vcfR::read.vcfR(vcf.path)

gl.data  <- vcfR::vcfR2genlight(vcf.data) 

pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(pop.data) %>% pull(Site)

table(pop(gl.data))


# Subset for nGSL / EGSL-1

ind.ESGL1 <- pop.data %>% dplyr::filter(RegionGen == "EGSL-1") %>% pull(ID_GQ)
length(ind.ESGL1)

gl.EGSL1 <- gl.data[indNames(gl.data) %in% ind.ESGL1, ]
gl.EGSL1


nLoc(gl.data)



# Reduce the N of SNPs ----------------------------------------------------

gl.data.sub <- list()
gl.EGSL1.sub <- list()

n.loc <-  c(20000, 15000, 10000, 5000, 2000, 1000, 500, 100)

for(i in n.loc){
 
 randsel1 <- sample(1:nLoc(gl.data), i, replace = FALSE)
 randsel2 <- sample(1:nLoc(gl.EGSL1), i, replace = FALSE) 
 
 gl.data.sub[paste0("Loc_", i)] <- gl.data[, randsel1]
 gl.EGSL1.sub[paste0("Loc_", i)] <- gl.EGSL1[, randsel2]
 
}

# PCA ---------------------------------------------------------------------

pca.all <- list()
pca.EGSL1 <- list()

for(i in n.loc){

  print(i)
  pca.all[[paste0("Loc_", i)]]   <- glPca( gl.data.sub[[paste0("Loc_", i)]] , center = TRUE, scale = FALSE,  
                    parallel = TRUE, n.core =16, nf = 100)
  
  #pca.EGSL1[[paste0("Loc_", i)]]   <- glPca( gl.EGSL1.sub[[paste0("Loc_", i)]] , center = TRUE, scale = FALSE,  
  #                                       parallel = TRUE, n.core =16, nf = 100)

  }


pca.EGSL1[[paste0("Loc_", i)]]

pca.fig <- list()
pca.fig.all <- list()

for(i in n.loc){

pca.fig[[paste0("Loc_", i)]] <-pca.EGSL1[[paste0("Loc_", i)]] %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, shape = Shore)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
  stat_ellipse() +
  #ggforce::geom_mark_ellipse(aes(fill = RegionGen, label = Site, group = Site), 
  #                            alpha = 0.25,  label.fontsize = 8, label.fontface = "plain", size = 0.2, con.size = 0.2 ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen)) +  
  scale_shape_manual(values = c(24,21)) +  
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  
    #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data), vjust=1, hjust=0) +
  
  labs(title = paste("N SNPs:",  i),
       x = paste0("PC1 (", QuickPop::pca_var(pca.EGSL1[[paste0("Loc_", i)]])$p.eig[1] %>% round(3) *100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(pca.EGSL1[[paste0("Loc_", i)]])$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()


}



fig1 <- ggpubr::ggarrange(pca.fig$Loc_20000, pca.fig$Loc_10000, pca.fig$Loc_5000,
                          pca.fig$Loc_2000, pca.fig$Loc_1000, pca.fig$Loc_500, #pca.fig$Loc_100,
                  common.legend = T)

fig1


ggsave(plot = fig1,
       filename = "02_Results/02_PowerTest/Test_SNPsPower.png",
       height = 8, width = 10, bg = "white")



for(i in n.loc){
  
  pca.fig.all[[paste0("Loc_", i)]] <-pca.all[[paste0("Loc_", i)]] %>% QuickPop::pca_scoretable(naxe = 5) %>%
    left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
    ggplot(aes(x = score.PC1, y = score.PC2)) +
    geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
    geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
    
    ggforce::geom_mark_ellipse(aes(fill = RegionGen, colour = RegionGen,label = NULL, group = Site), 
                               alpha = 0.10,  label.fontsize = 6, label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.3, expand = unit(2, "mm"), label.colour = "inherit" ) +
    geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
    scale_shape_manual(values = c(24,21)) +  
    scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
    scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
    
    labs(title = paste("N SNPs:",  i),
         x = paste0("PC1 (", QuickPop::pca_var(pca.all[[paste0("Loc_", i)]])$p.eig[1] %>% round(3) *100, "%)"),
         y = paste0("PC2 (", QuickPop::pca_var(pca.all[[paste0("Loc_", i)]])$p.eig[2] %>% round(3) *100, "%)")) +
    theme_bw()
  
  
}


fig2 <- ggpubr::ggarrange(pca.fig.all$Loc_20000, pca.fig.all$Loc_5000,
                          pca.fig.all$Loc_2000, pca.fig.all$Loc_1000, pca.fig.all$Loc_500, pca.fig.all$Loc_100,
                          common.legend = T)

fig2


ggsave(plot = fig2,
       filename = "02_Results/02_PowerTest/Test_SNPsPower_ALL_20250725.png",
       height = 8, width = 10, bg = "white")


# Save the random datasets ------------------------------------------------

#save(list = c("pca.fig", "pca.EGSL1", "pca.fig.all", "pca.all", "gl.EGSL1.sub", "gl.data.sub" ),
#     file = here("02_Results/02_PowerTest", "NlocSubsets.Rdata"))

load(here("02_Results/02_PowerTest", "NlocSubsets.Rdata"))
