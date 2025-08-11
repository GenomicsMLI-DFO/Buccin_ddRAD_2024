# Info --------------------------------------------------------------------

# Population genetics analysis of Buccinum undatum
#  431 individuals
#
# Audrey Bourret
# 2024-06-13
#

# Library -----------------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(ggforce)
library(Hmisc)

library(vcfR)
library(adegenet)
library(hierfstat)

`%nin%` = Negate(`%in%`)

library(QuickPop)

# Mapping

library(terra)
library(tidyterra)
library(rnaturalearth)
library(sf)

# install.packages("devtools")
#devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)
library(ggnewscale)
library(ggrepel)

# install.packages("shadowtext") 
library(shadowtext)


# Data --------------------------------------------------------------------

pop.data <- read_csv(file.path("00_Data/00_FileInfos/","Project_Infos_20240220.csv"))

pop.data 

pop.data <- pop.data %>% mutate(Notes_groupe = ifelse(Numero_unique_groupe == "G_23_00984","Groupe_EGSL_Baie_Comeau1_Offshore",
                                                      ifelse(Numero_unique_groupe == "G_23_00989","Groupe_EGSL_Baie_Comeau2_Offshore",       
                                                             Notes_groupe)),
                                Shore = ifelse(str_detect(Notes_groupe, "Offshore"), "Subtidal",
                                               ifelse(str_detect(Notes_groupe, "Inshore"), "Intertidal", "Problems"                 
                                               )),
                                Site = Notes_groupe %>% str_remove_all("Groupe_|_Offshore|_Inshore|Newfoundland_|EGSL_|Baie_de_Fundy_"),
                                Region = Notes_groupe %>% str_remove("Groupe_") %>% str_remove("_Offshore|_Inshore") %>% str_remove(paste(paste0("_", Site),  collapse = "|" )),
                               
                               
                                RegionGen = ifelse(Site %in% c("Bic1", "Bic2"), "SSLE", 
                                            ifelse(Site %in% c("Iles_de_la_Madeleine", "Mingan"), "GSL",
                                            ifelse(Region %in% c("EGSL"), "NSLE", 
                                            ifelse(Region %in% c("Baie_de_Fundy"), "Bay of Fundy",        
                                                   
                                                   str_replace_all(Region, "_", " ")))
                                                   )),
                                            Site = ifelse(Site == "GP", "Green's Point",
                                                          ifelse(Site == "SA", "Saint Andrews",
                                                                 ifelse(Site == "GL", "Gannet Lighthouse",
                                                                        ifelse(Site == "NL", "3L", str_replace_all(Site, "_", " " )))))
                                            )


other.SP <- c("S_23_01903", "S_23_01910", "S_23_01913", "S_23_01922")

# Genetic Data ------------------------------------------------------------

vcf.path <- file.path("./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.recode.vcf")

vcf.data <- vcfR::read.vcfR(vcf.path)

# Extract genotypes
gt <- extract.gt(vcf.data)

# Double-check missing rate
missing.rate <- sum(is.na(gt)) / length(gt) * 100
cat("% of missing SNPs :", round(missing.rate, 2), "%\n")

# gl object
gl.data  <- vcfR::vcfR2genlight(vcf.data) 

as.data.frame(gl.data[1:10, 1:10])

pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
                left_join(pop.data) %>% pull(Site)

nInd(gl.data)

table(pop(gl.data)) 

#gi object

gi.data  <- vcfR::vcfR2genind(vcf.data) 

pop(gi.data) <- data.frame(ID_GQ = indNames(gi.data)) %>% 
  left_join(pop.data) %>% pull(Site)



# Map ---------------------------------------------------------------------

admin <- terra::vect(rnaturalearth::ne_countries(scale = "large", continent = "north america", returnclass	  = "sf"))

pop <- terra::vect(pop.data %>% 
                     dplyr::filter(#Cat_sample == "Sample",
                                   #Espece == "undatum",
                                   ID_GQ %in% indNames(gl.data) # to keep only those with gen data
                                   ) %>% 
                     group_by(Site, RegionGen, Annee_echantillonnage, Longitude_echantillonnage_DD, Latitude_echantillonnage_DD) %>% summarise(N = length(ID_GQ)) 
                   ,
                   geom = c("Longitude_echantillonnage_DD", "Latitude_echantillonnage_DD"), keepgeom = T, crs = "+proj=longlat")

plot(crop(admin, ext(pop)))
plot(pop, add = T)

### Helper function to download bathymetry data ###
# function based on that from Michael Sumner Github gist (related to Pull Request for {terra})
# https://gist.github.com/mdsumner/aaa6f1d2c1ed107fbdd7e83f509a7cf3
get_elev <- function(x, method = "bilinear", maxcell = 25e6, silent = TRUE) {
  if (terra::ncell(x) > maxcell) {
    stop("number of cells in x exceeds 'maxcell'")
  }
  
  src <- "/vsicurl/https://gebco2023.s3.valeria.science/gebco_2023_land_cog.tif"
  # src <- "/vsicurl/https://public.services.aad.gov.au/datasets/science/GEBCO_2021_GEOTIFF/GEBCO_2021.tif"
  
  terra::project(terra::rast(src), x, by_util = TRUE, method = method)
}

# Load bathymetry
bathym <- get_elev(rast(ext(st_bbox(pop)[c('xmin','xmax','ymin','ymax')])+5,
                        crs = 'EPSG:4326',
                        res = 0.004166667)
)

# If interested, create contours as desired level(s)
isobath <- terra::as.contour(bathym, levels = c(-50, -250)) %>%
  st_as_sf() %>%
  st_cast("LINESTRING")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)

annotations <- data.frame(
  name = c("Bay of Fundy", "Scotian Shelf", 
            "Gulf of\nSt. Lawrence", "St. Lawrence\nEstuary", "Newfoundland"),
  x = c(250000, 550000, 500000, 200000,  900000),
  y = c(120000, 85000,  520000, 650000, 600000)
)

gg.fig1a <- ggplot(pop) + 
  #geom_spatraster(data = log(-bathym), alpha = 0.75,  show.legend = F) +
  #cmocean::scale_fill_cmocean(name = "gray", direction = -1) +
  geom_sf(data= isobath, aes(col = factor(level)), alpha = 0.75) + 
#  scale_color_manual(values = c( "cyan" , "deepskyblue3" , "blue4"), breaks =c(-50, -150, -300)  ,name = "Depth", ) +
  scale_color_manual(values = c( "deepskyblue" , "blue4"), breaks =c(-50, -250)  ,name = "Depth" ) +
  
    new_scale_color() +
    geom_sf(data = admin, fill="gray95", size=0.5, alpha = 1) +
  
  
  geom_rect(aes(xmin = -90000, xmax = 80000, ymin = 450000, ymax = 600000),
            color = "black", fill = NA, size = 0.5) +

  geom_shadowtext(data = annotations, aes(x = x, y = y, label = name),
            size = 3, bg.color = "white", fontface = "italic", color = c(rep("blue2", 4), "black")) +  
  
  geom_sf(data = pop, alpha = .5, size = 5,aes(col =  RegionGen
                                      #size = N, 
                                      #shape = Shore 
                                      )) +  
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  #scale_colour_manual(values = c("firebrick1",  "dodgerblue","limegreen", "black"), guide = "none") +
  #scale_color_viridis_c(name = "Depth (m)", direction = -1, option = "plasma")+
  # Map limits
  geom_sf_text_repel(data = pop %>% filter(RegionGen %in% c("Bay of Fundy", "GSL", "Newfoundland")),
                     aes(label = Site), pch = 16) +

  #scale_shape_manual(values = c(17, 16)) +
  #  coord_sf(xlim = c(-50000, 1100000), ylim = c(-20000, 1000000), crs = sf::st_crs("EPSG:6622")) +
  coord_sf(xlim = c(-50000, 1200000), ylim = c(-50000, 800000), crs = sf::st_crs("EPSG:6622")) +
  # Others
  #facet_grid(. ~ Category) +
  xlab("Longitude") + ylab("Latitude") +
  #ggtitle("ALL samples") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=10),
        strip.background = element_rect(fill = "white")#,
        #legend.position = "top"
        )  

gg.fig1a

gg.fig1b <- ggplot(pop) + 
  #geom_spatraster(data = log(-bathym), alpha = 0.75,  show.legend = F) +
  #cmocean::scale_fill_cmocean(name = "gray", direction = -1) +
  geom_sf(data= isobath, aes(col = factor(level)), alpha =0.75) + 
#  scale_color_manual(values = c( "cyan" , "deepskyblue3" , "blue4"), breaks =c(-50, -150, -300)  ,name = "Depth", ) +
  scale_color_manual(values = c( "deepskyblue" , "blue4"), breaks =c(-50, -250)  ,name = "Depth", ) +
  
    new_scale_color() +
  
   
  geom_sf(data = admin, fill="gray95", size=0.5, alpha = 1) +
  
  geom_shadowtext(data = annotations, aes(x = 0, y = 540000, label = "St. Lawrence Estuary"),
                  size = 3, bg.color = "white", fontface = "italic", color = "blue2") +  
  
  
  geom_sf(data = pop, alpha = .5, size = 5, aes(col =  RegionGen#, 
                                      #size = N, 
                                      #shape = Shore 
                                      )) +  
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
 
  
   
  #  scale_colour_manual(values = c("firebrick1",  "dodgerblue","limegreen", "black"), guide = "none") +
  #scale_color_viridis_c(name = "Depth (m)", direction = -1, option = "plasma")+
  # Map limits
  geom_sf_text(
    aes(label = Site), nudge_y = 4000) +
  scale_shape_manual(values = c(17, 16)) +
  #  coord_sf(xlim = c(-50000, 1100000), ylim = c(-20000, 1000000), crs = sf::st_crs("EPSG:6622")) +
  coord_sf(xlim = c(-80000, 90000), ylim = c(470000, 590000), crs = sf::st_crs("EPSG:6622")) +
  # Others
  #facet_grid(. ~ Category) +
  xlab("Longitude") + ylab("Latitude") +
  #ggtitle("ALL samples") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=10),
        strip.background = element_rect(fill = "white")#,
        #legend.position = "top"
  )  

gg.fig1b

gg.fig1 <- ggpubr::ggarrange(gg.fig1a, gg.fig1b, common.legend = T, align = "hv",
                             labels = LETTERS, legend = "none")
gg.fig1

ggsave(plot = gg.fig1,
       filename = "02_Results/01_PopStruct/Map_20250725.png",
       height = 4.5, width = 12, bg = "white")

# PCA ---------------------------------------------------------------------

pca.all  <- glPca(gl.data, center = TRUE, scale = FALSE,  
                  parallel = TRUE, n.core =16)
#save(list = c("pca.all" ),
#     file = here("02_Results/01_PopStruct", "PCA.Rdata"))

load( here::here("02_Results/01_PopStruct", "PCA.Rdata"))

gPCA.var <- QuickPop::pca_varplot(pca.all, 20)

gPCA.var

gPCA.1 <- pca.all %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +

  ggforce::geom_mark_ellipse(aes(fill = RegionGen, colour = RegionGen,label = ifelse(Site %in% c("Gannet Lighthouse"), as.character(Site), NA), group = Site), 
                             alpha = 0.10,  label.fontsize = 6, label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.3, expand = unit(2, "mm"), label.colour = "inherit" ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  scale_shape_manual(values = c(24,21)) +
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.all)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.all)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) + theme(legend.title = element_blank())
gPCA.1

gPCA.2 <- pca.all %>% QuickPop::pca_scoretable(naxe = 5) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
  
  ggforce::geom_mark_ellipse(aes(fill = RegionGen, colour = RegionGen, label = ifelse(Site %in% c("Gannet Lighthouse", "Iles de la Madeleine", "Saint Andrews", "Green's Point"), as.character(Site), NA),, group = Site), 
                             alpha = 0.10,  label.fontsize = 6,  label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.2, expand = unit(2, "mm"), label.colour = "inherit" ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
  scale_shape_manual(values = c(24,21)) +
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC3 (", QuickPop::pca_var(pca.all)$p.eig[3] %>% round(3) *100, "%)"),
    y = paste0("PC4 (", QuickPop::pca_var(pca.all)$p.eig[4] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.2


gPCA.3 <- pca.all %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC5, y = score.PC6)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
  
  #ggforce::geom_mark_ellipse(data = . %>% dplyr::filter(Site %in% "Green's Point"), 
  #                        aes(fill = RegionGen, colour = RegionGen, label = ifelse(Site %in% c("Gannet Lighthouse", "Green's Point"), as.character(Site), NA), group = Site), 
   #                       alpha = 0.10,  label.fontsize = 6, label.fontface = "bold", size = 0.2, con.size = 0.2 ,label.colour = "inherit") +  
  
  ggforce::geom_mark_ellipse(#data = . %>% dplyr::filter(Site %nin% "Green's Point"),
                             aes(fill = RegionGen, colour = RegionGen, label = ifelse(Site %in% c("Iles de la Madeleine",  "Mingan", "Green's Point", "3Ps", "3L"), as.character(Site), NA),, group = Site), 
                             alpha = 0.10,  label.fontsize = 6,label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.2, expand = unit(2, "mm"),label.colour = "inherit" ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
  scale_shape_manual(values = c(24,21)) +
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  
    #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC5 (", QuickPop::pca_var(pca.all)$p.eig[5] %>% round(3) *100, "%)"),
    y = paste0("PC6 (", QuickPop::pca_var(pca.all)$p.eig[6] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.3

gPCA.4 <- pca.all %>% QuickPop::pca_scoretable(naxe = 10) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC7, y = score.PC8)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +

  ggforce::geom_mark_hull(data = . %>% dplyr::filter(Site %in% "Green's Point"), 
                          aes(colour= RegionGen, fill= RegionGen, label = ifelse(Site %in% c("Gannet Lighthouse", "Green's Point"), as.character(Site), NA), group = Site), 
                          alpha = 0.15,  label.fontsize = 6, label.fontface = "bold", label.buffer = unit(2, "mm"), size = 0.2, con.size = 0.2, label.colour = "inherit", concavity = -1000 ) +  
  
  ggforce::geom_mark_ellipse(data = . %>% dplyr::filter(Site %nin% "Green's Point"), 
                             aes(colour= RegionGen, fill= RegionGen, label = ifelse(Site %in% c("Gannet Lighthouse", "Green's Point"), as.character(Site), NA), group = Site), 
                             alpha = 0.10,  label.fontsize = 6, label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.2, expand = unit(2, "mm"), label.colour = "inherit") +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +

#  ggforce::geom_mark_ellipse(aes(fill = RegionGen, label = ifelse(Site %in% c("Gannet Lighthouse", "Green's Point"), as.character(Site), NA),, group = Site), 
#                             alpha = 0.25,  label.fontsize = 8, label.fontface = "bold", size = 0.2, con.size = 0.2 ) +
#  geom_point(alpha = 0.5, size = 1) + 
  

  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC7 (", QuickPop::pca_var(pca.all)$p.eig[7] %>% round(3) *100, "%)"),
    y = paste0("PC8 (", QuickPop::pca_var(pca.all)$p.eig[8] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) #+ theme(legend.position = "none")
gPCA.4

gPCA <- ggpubr::ggarrange(gPCA.1, gPCA.2,gPCA.3,gPCA.4,
                          nrow = 2, ncol = 2, common.legend = T, legend = "top",
                          labels = LETTERS)

gPCA

ggsave(filename = file.path("02_Results/01_PopStruct/PCA_20250725.png"),
       plot = gPCA,
       width = 9, height = 8, units = "in", bg = "white")



# Admixture ---------------------------------------------------------------

# A2

vcf.path

cmd1 <- paste("--vcf", vcf.path,
              #"--recode",
              "--plink-tped",
              "--out", "./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final")


cmd1

A1 <- system2("vcftools", cmd1, stdout=T, stderr=T)



# Correct the map file to remove chromosomes (change to 0)
tped.plink <- read.delim("00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.tped",
                                header = F)
tped.plink[1:6, 1:6]
tped.plink$V1 <- 0

write.table(tped.plink, 
            file = "00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.tped",
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


cmd2a <- paste("--tfam", "./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.tfam",
               "--tped", "./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.tped",
               "--make-bed",
               "--out", "./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final")



A2a <- system2("/home/genyoda/Documents/Programs/plink_linux_x86_64_20210606/plink", cmd2a, stdout=T, stderr=T)
A2a

# Original SNP

bed.file <- file.path(here::here(),"./00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.bed" )
file.exists(bed.file)
fam.file <- bed.file %>% str_replace(".bed", ".fam")
fam <- read.table(fam.file)
file.exists(bed.file)


for(k in 1:20){
  
  print(k)  
  
  setwd(file.path(here::here(), "02_Results/01_PopStruct/Admixture/") ) 
  
  cmd <- paste("--cv", # to perform cross-validation in the log file 
               bed.file,
               k, # the number of K
               "-j8"#
  )
  
  A <- system2("admixture", cmd, stdout = T, stderr = T) 
  
  cat(file = paste0("Bundatum.k",k, ".log"),
      "\n", cmd, "\n",
      A, # what to put in my file
      append= F, sep = "\n")
  
  setwd(here::here())
  
}

# Cross-validation results:

CV.all.res <- data.frame(k = 1:20,
                         CV = NA,
                         stringsAsFactors = F)


for(i in 1:nrow(CV.all.res)){
  # Which k
  k <- CV.all.res[i, "k"]
  
  # Extract from the log file
  temp <- readLines(file.path("./02_Results/01_PopStruct/Admixture/", paste0("Bundatum.k",k, ".log")))
  CV.temp <- temp %>% str_subset("CV error")
  CV <- sapply(str_split(CV.temp, ":"), `[`, 2) %>% str_remove_all(" ")
  
  # Add to my data.frame
  CV.all.res[i, "CV"] <- CV
  
}

CV.all.res$CV <- as.numeric(as.character(CV.all.res$CV))

CV.all.res %>% arrange(CV)

plot(CV.all.res$CV)

gg.CV.all <- CV.all.res %>% mutate(color = ifelse(k %in% c(2,3,4,10), "red", "black")) %>% 
  ggplot(aes(x = factor(k), y = CV)) + 
  geom_point(size = 2, aes(col = color)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "K", y = "Cross-validation error") +
  theme_bw() +
  theme(legend.position = "none")
gg.CV.all

ggsave(filename = file.path(here::here(), "02_Results/01_PopStruct/", "Admixture.ALL.CV.png"), 
       plot = gg.CV.all,
       height = 3.5, width = 4, units = "in")   


gg.struct <- ggpubr::ggarrange(gPCA.var + labs(title = NULL), 
                               gg.CV.all,
                               labels = LETTERS) 
gg.struct 


ggsave(filename = file.path(here::here(), "02_Results/01_PopStruct/", "FigS2_20250304.png"), 
       plot = gg.struct ,
       height = 3.5, width = 7, units = "in")   




k.all <- 15

Q.all.k2.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.2.Q"))
Q.all.k3.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.3.Q"))
Q.all.k4.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.4.Q"))
Q.all.k5.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.5.Q"))
Q.all.k6.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.6.Q"))
Q.all.k7.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.7.Q"))
Q.all.k8.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.8.Q"))
Q.all.k9.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.9.Q"))
Q.all.k10.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.10.Q"))
Q.all.k11.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.11.Q"))
Q.all.k12.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.12.Q"))
Q.all.k13.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.13.Q"))
Q.all.k14.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.14.Q"))
Q.all.k15.res <-  read.table(file.path(here::here(), "02_Results/01_PopStruct/Admixture/populations.23405snps.431ind.single.final.15.Q"))



Q.all.res <- bind_rows(  
  cbind(fam$V1, Q.all.k15.res, K = 15),
  cbind(fam$V1, Q.all.k14.res, K = 14),
  cbind(fam$V1, Q.all.k13.res, K = 13),
  cbind(fam$V1, Q.all.k12.res, K = 12),
  cbind(fam$V1, Q.all.k11.res, K = 11),
  cbind(fam$V1, Q.all.k10.res, K = 10),
  cbind(fam$V1, Q.all.k9.res, K = 9),
  cbind(fam$V1, Q.all.k8.res, K = 8),
  cbind(fam$V1, Q.all.k7.res, K = 7),
  cbind(fam$V1, Q.all.k6.res, K = 6),
  cbind(fam$V1, Q.all.k5.res, K = 5),
  cbind(fam$V1, Q.all.k4.res, K = 4),
  cbind(fam$V1, Q.all.k3.res, K = 3),
  cbind(fam$V1, Q.all.k2.res, K = 2)
)
names(Q.all.res) <- c("ID_GQ", paste0("Q", 1:k.all), "K")

names(Q.all.k3.res) <- c("V2", "V1", "V3")
names(Q.all.k4.res) <- c("V1", "V2", "V4", "V3")
names(Q.all.k10.res) <- c("V4", "V1",   "V5","V6", "V7", "V8", "V3", "V9", "V10", "V2" )

 
 Q.all.res <- bind_rows(  
   cbind(fam$V1, Q.all.k10.res, K = 10),
   cbind(fam$V1, Q.all.k4.res, K = 4),
   cbind(fam$V1, Q.all.k3.res, K = 3),
   cbind(fam$V1, Q.all.k2.res, K = 2)
 )
 
 names(Q.all.res) <- c("ID_GQ", paste0("Q", 1:10), "K")

head(Q.all.res)

#Q.all.res <- cbind(fam.all$V1, Q.all.res)

#reorder(ID, Qvalue, FUN = function(x) min(x))

max.df <- Q.all.res %>% dplyr:::filter(K == 4) %>% mutate(MaxQ = Q1 ) %>% dplyr::select(ID_GQ, MaxQ)

gg.str.all <- Q.all.res %>% 
  

  pivot_longer(cols =  paste0("Q", 1:10), names_to = "Group", values_to = "Q") %>% 
  # mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
  #  dplyr::filter(ID_GQ %nin% bad.samples) %>% 
  dplyr::filter(K %in% c(2,3,4,10),
                Q != 0) %>%
  left_join(max.df) %>% 
  left_join(pop.data) %>% 
  mutate(Site = factor(Site, levels = c( "Green's Point",  "Saint Andrews", "Gannet Lighthouse",
                                          "Iles Penchees", "Baie des Bacon", "Pointe a Emile", "Forestville", "Portneuf",  "Baie Comeau1", "Baie Comeau2",
                                         "Bic1", "Bic2", "Iles de la Madeleine","Mingan", "3L", "3Ps"
                                         
                                         ))) %>% 
  
  
    ggplot(aes(x = reorder(ID_GQ, MaxQ, median), y = Q, fill = Group)) + 
  
  geom_col() +
  facet_grid(K ~ Site, space = "free", scale = "free") +
  #facet_grid(K ~ paste(RegionGen, Site, sep = ": "), space = "free", scale = "free") +
  #scale_fill_viridis_d(option = "A") +
  scale_fill_brewer(palette = "Paired") +
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
gg.str.all


ggsave(plot = gg.str.all,
       filename = "02_Results/01_PopStruct/Admixture.ALL_20250304.png",
       width = 8, height = 7)

ggsave(plot = gg.str.all,
       filename = "02_Results/01_PopStruct/Admixture.ALL_20250304.pdf",
       device = cairo_pdf,
       width = 8, height = 7)

ggsave(plot = gg.str.all,
       filename = "02_Results/01_PopStruct/Admixture.ALL_20250304.svg",
       width = 8, height = 7)

ggsave(plot = gg.str.all,
       filename = "02_Results/01_PopStruct/Admixture.ALL_20250304.eps",
       width = 8, height = 7)



gg.str.k15 <- Q.all.res %>% 
  
  
  pivot_longer(cols =  paste0("Q", 1:15), names_to = "Group", values_to = "Q") %>% 
  # mutate(Group = factor(Group, levels = c("Q1", "Q2", "Q4", "Q3"))) %>% 
  #  dplyr::filter(ID_GQ %nin% bad.samples) %>% 
  dplyr::filter(K %in% c(15),
                Q != 0) %>%
  left_join(max.df) %>% 
  left_join(pop.data) %>% 
  mutate(Site = factor(Site, levels = c( "Green's Point",  "Saint Andrews", "Gannet Lighthouse",
                                         "Iles Penchees", "Baie des Bacon", "Pointe a Emile", "Forestville", "Portneuf",  "Baie Comeau1", "Baie Comeau2",
                                         "Bic1", "Bic2", "Iles de la Madeleine","Mingan", "3L", "3Ps"
                                         
  ))) %>% 
  
  
  ggplot(aes(x = reorder(ID_GQ, MaxQ, median), y = Q, fill = Group)) + 
  
  geom_col() +
  facet_grid(K ~ Site, space = "free", scale = "free") +
  #facet_grid(K ~ paste(RegionGen, Site, sep = ": "), space = "free", scale = "free") +
  #scale_fill_viridis_d(option = "A") +
  #scale_fill_brewer(palette = "Paired") +
  labs(y="Membership probability") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.text.y = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "white", colour  = "white"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt") )
gg.str.k15


gg.struct2 <- ggpubr::ggarrange(ggpubr::ggarrange(gPCA.var + labs(title = NULL), 
                                                  gg.CV.all,
                                                  labels = LETTERS),
                                gg.str.k15, labels = c("", "C"),
                                nrow = 2, heights = c(2,3))

gg.struct2 


ggsave(filename = file.path(here::here(), "02_Results/01_PopStruct/", "FigS2_20250619.png"), 
       plot = gg.struct2,
       height = 7, width = 7, units = "in")   


# Fst ---------------------------------------------------------------------

library(dartR)

table(pop(gl.data))

res.FST <- dartR::gl.fst.pop(gl.data, nboots = 999, percent = 95, nclusters = 20)

write_csv(res.FST$Bootstraps,
          "02_Results/01_PopStruct/FSTwBootstrap.csv")


res.FST.Bootstraps <- read_csv("02_Results/01_PopStruct/FSTwBootstrap.csv")

res.FST.Bootstraps$`p-value` %>% max()
res.FST.Bootstraps$Fst %>% mean()
res.FST.Bootstraps$Fst %>% min()
res.FST.Bootstraps$Fst %>% max()

# Compared

# Functions

table.fst <- function(fst){
  res <-  fst %>% dplyr::select(Population1, Population2, "Lower bound CI limit", "Upper bound CI limit", "p-value", "Fst")
  
  return(res)
  
}

heat.fst <- function(fst){
  res <- bind_rows(table.fst(fst),
                   table.fst(fst) %>% mutate(Pop1.int = Population2,
                                             Pop2.int = Population1,
                                             Population1 = Pop1.int,
                                             Population2 = Pop2.int) %>% 
                     dplyr::select(-c(Pop1.int, Pop2.int))
  )
  
  return(res)         
  
}


heat.graph.panel <- function(heat){
  graph <-  heat %>%  
    ggplot(aes(x=Population1, y=Population2, fill=Fst, shapp = Sign)) + 
    geom_tile(colour=  "gray") +
    geom_point(aes(x = Population1, y=Population2, pch = Sign), size = 3)+
    #scale_y_discrete(limits=rev) +
    #scale_x_discrete(limits=rev) +
    scale_fill_gradient(low = "ivory1", high = "red3")+
    # scale_fill_gradient(low = "ivory1", high = "dodgerblue2")+
    scale_shape_manual(values = c("","*"), guide = "none") +
    #scale_fill_distiller(palette = "Spectral") +
    labs(x = NULL, y = NULL) +
    #theme_bw() +
    facet_grid(SFA2~SFA1, space = "free", scale = "free") +
    theme_minimal() + 
    theme(#axis.text.x = element_blank(),
      strip.text = element_text(angle = 0),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "cm"),
      panel.border = element_rect(fill = NA, colour = "black"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      #legend.background = element_rect(colour = "black", fill = "grey95", size = 0.2),
      legend.position = "right",
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))
  print(graph)
  return(graph)
}


g0.fst <- res.FST.Bootstraps %>% heat.fst() %>%  
  left_join(pop.data %>% dplyr::select(Population1 = Site, Region1 = RegionGen)  %>% distinct()) %>% 
  left_join(pop.data %>% dplyr::select(Population2 = Site, Region2 = RegionGen)  %>% distinct()) %>% 
  mutate(Sign = ifelse(`p-value` <= 0.05, "*", ""),
         Region1 = factor(Region1, levels = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland" )), 
  Region2 = factor(Region2, levels = c("Newfoundland",  "GSL", "SSLE", "NSLE", "Bay of Fundy")),
  Population1 = factor(Population1, levels = c( "Green's Point",  "Saint Andrews", "Gannet Lighthouse",
                                         "Iles Penchees", "Baie des Bacon", "Pointe a Emile", "Forestville", "Portneuf",  "Baie Comeau1", "Baie Comeau2",
                                         "Bic1", "Bic2", "Iles de la Madeleine","Mingan", "3L", "3Ps")) ,
  Population2 = factor(Population2, levels = c( "Green's Point",  "Saint Andrews", "Gannet Lighthouse",
                                                "Iles Penchees", "Baie des Bacon", "Pointe a Emile", "Forestville", "Portneuf",  "Baie Comeau1", "Baie Comeau2",
                                                "Bic1", "Bic2", "Iles de la Madeleine","Mingan", "3L", "3Ps"))                                        
                                         )%>% 
  
    ggplot(aes(x=Population1, y=Population2, fill=Fst)) + 
  geom_tile(colour=  "gray") +
  # geom_point(aes(x = Population1, y=Population2), size = 3)+
  geom_text(aes(label = paste0(format(round(Fst, dig = 3 ), scientific=F)))) +
  #scale_y_discrete(limits=rev) +
  #scale_x_discrete(limits=rev) +
  scale_fill_gradient(low = "ivory1", high = "red3")+
  # scale_fill_gradient(low = "ivory1", high = "dodgerblue2")+
  scale_shape_manual(values = c("","*"), guide = "none") +
  #scale_fill_distiller(palette = "Spectral") +
  labs(x = NULL, y = NULL) +
  #theme_bw() +
  facet_grid(Region2~Region1, space = "free", scale = "free") +
  theme_minimal() + 
  theme(#axis.text.x = element_blank(),
    strip.text = element_text(angle = 0),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "cm"),
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    #legend.background = element_rect(colour = "black", fill = "grey95", size = 0.2),
    #legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle =  90, vjust = 0.5, hjust = 1))

#g0.fst <-  g1.fst + scale_y_discrete(limits=rev) + scale_x_discrete(limits=rev)  +    scale_fill_viridis_c(option = "C")
g0.fst + scale_y_discrete(limits=rev) + scale_x_discrete(limits=rev)  #    scale_fill_viridis_c(option = "C")
g0.fst

ggsave(plot = g0.fst,
       filename = "02_Results/01_PopStruct/Fst_20250722.png",
       width = 12, height = 7)


res.FST.Bootstraps %>% dplyr::filter(`p-value` > 0.01) %>% select(Fst, `p-value`)

# PCadapt -----------------------------------------------------------------

library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")


# Convertion to plink .bed format

SCAFFOLD.info <- vcf.data@fix %>% as.data.frame() %>%  
  select(ID, CHROM, POS) %>% 
  mutate(scaffold = sapply(str_split(CHROM, ","), `[`,1) %>% str_remove("scaffold"),
         RADloc = sapply(str_split(ID, ":"), `[`,1),
         scaffold.length = sapply(str_split(CHROM, ","), `[`,2))

SCAFFOLD.info$RADloc %>% unique() %>% length()



# Read .bed in PCAadapt
pcadapt.genotype  <- read.pcadapt(file.path("00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.bed"),
                                      type = "bed")

pcadapt.snp <- read.delim(file.path("00_Data/06a_Filtering.denovo.buccinum.M1n1/A_90K/08_Final/populations.23405snps.431ind.single.final.bim" ),
                              header = F) %>% pull(V2)

length(pcadapt.snp) 
nLoc(gl.data)

# Run pcadapt

K.init <- 20

pcadapt.test   <- pcadapt(pcadapt.genotype, K =K.init)
# Check screeplot

plot(pcadapt.test, option = "screeplot") 

# Check structure

plot(pcadapt.test, option = "scores") 

pcadapt.k3 <- pcadapt(pcadapt.genotype , K = 3)

plot(pcadapt.k3, option = "manhattan")
hist(pcadapt.k3$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(pcadapt.k3, option = "qqplot")

# Statistics
#x$pvalues 
alpha <- 0.05

qval.k3 <- qvalue::qvalue(pcadapt.k3$pvalues)$qvalues
outliers.k3 <- which(qval.k3 < alpha)
length(outliers.k3)
#outliers.high.sfa.k2 <- c(which(is.na(qval.sfa.k2)), which(qval.sfa.k2 < 3.5e-6)) 

pcadapt.outliers <- pcadapt.snp[outliers.k3]

bonf.k3 <- p.adjust(pcadapt.k3$pvalues,method="bonferroni")
outliers.bonf.k3 <- which(bonf.k3 < alpha)
length(outliers.bonf.k3)

pcadapt.k3$pvalues

res.pcadapt <- data.frame(ID = pcadapt.snp,
                              pvalues = pcadapt.k3$pvalues,
                              qvalues = qval.k3) %>% 
#  left_join(SCAFFOLD.info %>% dplyr::select(ID, CHROM, POS)) %>% 
  mutate(Outliers.qvalue = ifelse(ID %in% pcadapt.snp[outliers.k3], T, F),
         Outliers.bonferroni = ifelse(ID %in% pcadapt.snp[outliers.bonf.k3], T, F),
         Cat.SNP = ifelse(Outliers.bonferroni == T, "Strong outlier",
                          ifelse(Outliers.qvalue == T, "Outlier", "Neutral")))

man.pcadapt <-  res.pcadapt %>% 
  ggplot(aes(x = ID, y = -log10(pvalues), col =Outliers.bonferroni  )) +
  geom_point(alpha = 0.5)+
  scale_colour_manual(values = c("gray", "darkorange"))  +
  labs(x = "Position") +
#  facet_grid(. ~ as.numeric(CHROM2 %>% str_remove("Chr")), space = "free", scale = "free") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        axis.ticks.x =  element_blank(), 
        panel.spacing = unit(0, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
man.pcadapt


pca.neutral <- glPca(gl.data[ , locNames(gl.data) %nin%  pcadapt.snp[outliers.k3] ] , center = TRUE, scale = FALSE,  
                         parallel = TRUE, n.core = 20, nf = 10)

pca.outlier <- glPca(gl.data[ , locNames(gl.data) %in%  pcadapt.snp[outliers.bonf.k3] ] , center = TRUE, scale = FALSE,  
                         parallel = TRUE, n.core = 20, nf = 10)

QuickPop::pca_varplot(pca.neutral, 10)
QuickPop::pca_varplot(pca.outlier, 10)

#save(list = c("pca.neutral", "pca.outlier"),
#     file = "./02_Results/01_PopStruct/PCA_neutral_outlier_20240621.Rdata")

load("./02_Results/01_PopStruct/PCA_neutral_outlier_20240621.Rdata")

gg.pca.neutral1 <- pca.neutral%>% QuickPop::pca_scoretable(naxe = 10) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
  
  ggforce::geom_mark_ellipse(aes(fill = RegionGen, colour = RegionGen,label = ifelse(Site %in% c("Gannet Lighthouse"), as.character(Site), NA), group = Site), 
                             alpha = 0.10,  label.fontsize = 6, label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.3, expand = unit(2, "mm"), label.colour = "inherit" ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch= 21) + 
  scale_shape_manual(values = c(24,21)) +  
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  labs(title = "Putative neutral SNPs (n = 22,523 SNPs)",
    #subtitle = "Divided by sample category",
    x = paste0("PC1 (", QuickPop::pca_var(pca.neutral)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.neutral)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) + theme(legend.title = element_blank())
gg.pca.neutral1 


gg.pca.neutral2 <- pca.neutral%>% QuickPop::pca_scoretable(naxe = 10) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC3, y = score.PC4)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
  
  ggforce::geom_mark_ellipse(aes(fill = RegionGen, colour = RegionGen,label = ifelse(Site %in% c("Iles de la Madeleine", "Mingan"), as.character(Site), NA), group = Site), 
                             alpha = 0.10,  label.fontsize = 6, label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.3, expand = unit(2, "mm"), label.colour = "inherit" ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
  scale_shape_manual(values = c(24,21)) +  
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  labs(title = "Putative neutral SNPs (n = 22,523 SNPs)",
       #subtitle = "Divided by sample category",
       x = paste0("PC3 (", QuickPop::pca_var(pca.neutral)$p.eig[3] %>% round(3) *100, "%)"),
       y = paste0("PC4 (", QuickPop::pca_var(pca.neutral)$p.eig[4] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) + theme(legend.title = element_blank())
gg.pca.neutral2 

gg.pca.outlier <- pca.outlier %>% QuickPop::pca_scoretable(naxe = 10) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2)) +
  geom_hline(yintercept = 0, col = "darkgray", size = 0.5) +
  geom_vline(xintercept = 0, col = "darkgray", size = 0.5) +
  
  ggforce::geom_mark_ellipse(aes(fill = RegionGen, colour = RegionGen,label = ifelse(Site %in% c("Gannet Lighthouse"), as.character(Site), NA), group = Site), 
                             alpha = 0.10,  label.fontsize = 6, label.buffer = unit(2, "mm"), label.fontface = "bold", size = 0.2, con.size = 0.3, expand = unit(2, "mm"), label.colour = "inherit" ) +
  geom_point(alpha = 0.5, size = 1.5, aes(fill = RegionGen), pch = 21) + 
  scale_color_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  scale_fill_manual(values = cols, breaks = c("Bay of Fundy", "NSLE", "SSLE", "GSL", "Newfoundland")) +
  
  scale_shape_manual(values = c(24,21)) +  
  labs(title = "Outlier SNPs (n = 284 SNPs)",
       #subtitle = "Divided by sample category",
       x = paste0("PC1 (", QuickPop::pca_var(pca.outlier)$p.eig[1] %>% round(3) *100, "%)"),
       y = paste0("PC2 (", QuickPop::pca_var(pca.outlier)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw(base_size = 8) + theme(legend.title = element_blank())
gg.pca.outlier

gg.pca.outlier + facet_grid(~RegionGen)


# Figure PCADAPT
gg.comp.pca <- ggpubr::ggarrange(gg.pca.neutral1, gg.pca.neutral2, gg.pca.outlier,
                  man.pcadapt,
  labels = LETTERS,
  nrow = 2, ncol = 2, common.legend = T, align = "hv"
)

gg.comp.pca


ggsave(filename = file.path("02_Results/01_PopStruct/PCA_neutralVSoutlier_20250725.png"),
       plot = gg.comp.pca,
       width = 8, height = 7, units = "in", bg = "white")

