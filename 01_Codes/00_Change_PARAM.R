# Info --------------------------------------------------------------------

# This file allowed to set/modify folder structure
#
# Audrey Bourret
# 2021-01-12
#

# Library -----------------------------------------------------------------

# Internal functions
if(length(list.files("./01_Codes/Functions"))>=1){ 
  for(i in 1:length( list.files("./01_Codes/Functions") )){
    source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
  }
}

PARAM <- data.frame(x = character(), value = character())

# Specific functions ------------------------------------------------------

"add.param<-" <- function(PARAM, ..., value){
  PARAM <- rbind(PARAM, value)
  names(PARAM) <- c("x", "value")
  PARAM$x     <- as.character(PARAM$x)
  PARAM$value <- as.character(PARAM$value)
  PARAM
  
}


# Path --------------------------------------------------------------------

# 00_Data
add.param(PARAM) <- c("info.path", "./00_Data/00_FileInfos")

add.param(PARAM) <- c("raw.path", "./00_Data/01_Raw") 
#add.param(PARAM) <- c("cutadapt.path", "./00_Data/02_Cutadapt") 
add.param(PARAM) <- c("trimmo.path", "./00_Data/02_Trimmomatic") 

add.param(PARAM) <- c("demulti.path", "./00_Data/03a_Demultiplex") 
add.param(PARAM) <- c("align.path", "./00_Data/03b_Align") # External disk
#add.param(PARAM) <- c("align.mtDNA.path", "./00_Data/03c_Align_mtDNA") # External disk

#add.param(PARAM) <- c("samples.path", "./00_Data/04_Samples") # External disk

add.param(PARAM) <- c("test.denovo.path", "./00_Data/04a_Test.denovo") # External disk
add.param(PARAM) <- c("test.ref.path", "./00_Data/04b_Test.ref") # External disk

#add.param(PARAM) <- c("ustacks.denovo.path", "./00_Data/05a_Stacks.denovo") # External disk
add.param(PARAM) <- c("stacks.denovo.path", "./00_Data/05a_Stacks.denovo") # External disk
add.param(PARAM) <- c("stacks.ref.path", "./00_Data/05b_Stacks.ref") # External disk


add.param(PARAM) <- c("filter.denovo.path", "./00_Data/06a_Filtering.denovo")
add.param(PARAM) <- c("filter.ref.path", "./00_Data/06b_Filtering.ref")

add.param(PARAM) <- c("ref.genome.path", "./00_Data/99_REF_Genome")

# Results - Stacks

add.param(PARAM) <- c("fastqc.path", "./02_Results/00_Stacks/01_FastQC")

add.param(PARAM) <- c("fastqc.raw.path", "./02_Results/00_Stacks/01_FastQC/01_Raw")
add.param(PARAM) <- c("fastqc.trim.path", "./02_Results/00_Stacks/01_FastQC/02_Trimmomatic")
add.param(PARAM) <- c("fastqc.demulti.path", "./02_Results/00_Stacks/01_FastQC/03_Demultiplex")


add.param(PARAM) <- c("cutadapt.log", "./02_Results/00_Stacks/02_Cutadapt")
add.param(PARAM) <- c("trimmo.log", "./02_Results/00_Stacks/02_Trimmomatic")
add.param(PARAM) <- c("demulti.log", "./02_Results/00_Stacks/03_Demultiplex")
add.param(PARAM) <- c("test.denovo.log", "./02_Results/00_Stacks/04a_Test.denovo")
add.param(PARAM) <- c("test.ref.log", "./02_Results/00_Stacks/04b_Test.ref")
add.param(PARAM) <- c("stacks.log", "./02_Results/00_Stacks/05a_Stacks.denovo")
add.param(PARAM) <- c("stacks.ref.log", "./02_Results/00_Stacks/05b_Stacks.ref")

#add.param(PARAM) <- c("filter.data", "./02_Results/00_Stacks/06_Filtering")

#add.param(PARAM) <- c("results", "./02_Results/07_Analysis")
#add.param(PARAM) <- c("admix.res", "./02_Results/07_Analysis/07a_Admixture")
#add.param(PARAM) <- c("fastStruct.res", "./02_Results/07_Analysis/07b_FastStructure")

# Update file before continue to ensure next section will do OK!
if(dir.exists("./00_Data/00_FileInfos") == F){
  dir.create("./00_Data/00_FileInfos")
}

write.csv2(PARAM, file = file.path("./00_Data/00_FileInfos", "Options.csv"), row.names=F)

# Files -------------------------------------------------------------------

# External programs -------------------------------------------------------

# Data --------------------------------------------------------------------


# Save Parameters ---------------------------------------------------------

PARAM

write.csv2(PARAM, file = file.path("./00_Data/00_FileInfos", "Options.csv"), row.names=F)

get.value("raw.path")



