# Info --------------------------------------------------------------------

# RAD-seq pipeline with STACKS version >2, up-to Gstacks
# NS data
# Trimmomatics + with rad check on both sides
# With and without ref genome
#
# Audrey Bourret
# 2024-02-21
#


# Library -----------------------------------------------------------------

library(parallel)
library(tidyverse)

# Internal functions
for(i in 1:length( list.files("./01_Codes/Functions") )){
  source(file.path("./01_Codes/Functions",  list.files("./01_Codes/Functions")[i]))  
}

`%nin%` = Negate(`%in%`)

# Add pythom env to this specific project
Sys.setenv(PATH = paste(c("/home/genyoda/Documents/PythonVenv/GenoBaseEnv/bin",
                          Sys.getenv("PATH")),
                        collapse = .Platform$path.sep))

system2("multiqc", "--help")

# Check folders -----------------------------------------------------------
# To run if you need the folder skeleton 

#auto.folder()

# Data --------------------------------------------------------------------

pop.info <- read.table(file.path(get.value("info.path"), "popmap.txt"))
names(pop.info) <- c("Sample", "POP")
pop.info


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
pop.data %>% group_by(Region, Shore, Site ) %>% summarise(N = n())

# Define initial working directory (just in case something goes wrong)
current.wd <- getwd()

numCores <- if(get_os() %in% c("os","linux")){
  detectCores() # Utilise le max de coeurs si sur linux
} else 1

cat("There's", numCores, "cores available in this computer", sep = " ")

# FastQC ------------------------------------------------------------------

# Took around 30 min by plate (NS and HI)
# Now with an overwrite function

fastqc <- function(folder.in, folder.out, overwrite = F, nthread, fastq = "fastq") {
#@folder.in : where to find the data
#@folder.out : where to put the data
#@overwrite : if T, will remove everything within the folder  
    
  if(get_os() %in% c("os","linux")){ # to run only on os and not windows
    
    cat("Performing a FastQC analysis\n")
    
    # Remove old files or not (addition June 2020)
    if(isTRUE(overwrite)){
      file.remove(list.files(folder.out, full.name = T, pattern =fastq))
      files.to.use <- list.files(folder.in, full.names = T, pattern =fastq) %>% 
        str_remove(".md5") %>% 
        str_subset(paste0("1.",fastq,"|2.", fastq)) %>% unique()
    }
    
    if(isTRUE(overwrite == F)){
      previous.analysis <- list.files(folder.out, pattern = ".html") %>% str_replace(paste0("_",fastq,".html"), paste0(".",fastq,".gz"))
      
      files.to.use <- list.files(folder.in, full.names = T, pattern = fastq) %>% 
        str_remove(".md5") %>% unique() %>% 
        str_subset(paste(previous.analysis, collapse = "|"), negate = T) %>% 
        str_subset(paste0("R1.",fastq,"|R2.", fastq))
      
    }
    
    cat("Results could be find here:", folder.out ,"\n")
    
    mclapply(files.to.use,
             FUN = function(x){
               
               cmd <- paste("--outdir", folder.out, x, 
                            "-t", 1)
               system2("fastqc", cmd) 
              
             } ,
             mc.cores = nthread
    )
    

  } else {cat("Cannot perform FastQC on windows yet -- sorry!!")}
} # End of my function

# Test a multiqc function

multiqc <- function(folder.in){
  # Multi QC aggregation - run pretty fast ...
  for(s in c("1_fastqc", "2_fastqc")){
    print(s)  
    cmd <- paste(paste(list.files(folder.in, full.names = T) %>% 
                         #str_subset(l) %>% 
                         str_subset(paste0("",s)) %>% 
                         str_subset(".zip"), collapse = " "),
                 "--outdir", file.path(folder.in, "MultiQC_report"),
                 "--filename", paste0("multiqc_report_", s, ".html"),
                 "-f" # to rewrite on previous data
    )
    
    system2("multiqc", cmd)
    
  } 
  
}



# Run FastQC - change for overwrite T for the first time to or it will not work
fastqc(folder.in = get.value("raw.path"), folder.out = get.value("fastqc.raw.path"), overwrite = T, nthread = 20, fastq = "fastq")


# Multi QC aggregation - run pretty fast ...

multiqc(folder.in = get.value("fastqc.raw.path"))

# CHECK on the R2 side that the adapter "CGG" is present. 
# IF not, the PART2 of trimmomatic allows to remove 5 pb on this side

# Trimmomatics ------------------------------------------------------------

# In paired-end
# To remove the Illumina adapter and to cut 3pb in R2

# Check that you can reach the program
trimmomatic.path <- "/home/genyoda/Documents/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
system2("java", paste("-jar", trimmomatic.path, "PE", "-version"), stdout=T, stderr=T) 

# Check which files you want to use

files.to.use <- list.files(get.value("raw.path"), full.names = T) %>% 
                           str_subset("R1") %>% 
                           str_subset(".md5", negate = T)

files.to.use


#  REMOVE Illumina adaptors
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Names of the important files
             file1 <- x
             
             # Original file
             file2 <- file1 %>% str_replace("_R1", "_R2")
             # Version with 3 pb removed
             #file2 <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
             #                    str_replace("_R1.fastq.gz", "_R2_HC3.fastq.gz")
             
             fileout <- file1 %>%  str_replace(get.value("raw.path"), get.value("trimmo.path")) %>% 
                                   str_replace("_R1.fastq.gz", ".fastq.gz")   
             
             logout <- file1 %>% str_replace(get.value("raw.path"), get.value("trimmo.log")) %>% 
                                 str_replace("_R1.fastq.gz", ".log")
             
             # The command

             cmd1 <- paste("-jar",
                           trimmomatic.path, "PE",
                           "-threads", 1,
                           #"-trimlog", logout,
                           file1, file2,
                           "-baseout", fileout,
                           "ILLUMINACLIP:00_Data/00_FileInfos/adapters/TruSeq3-PE-2.fa:2:30:10:8:TRUE",
                           #"HEADCROP:1-99",
                           sep = " ")
             
             A1 <- system2("java", cmd1, stdout=T, stderr=T) # to R console
             A1
             
             # save a file log
             cat(file = logout,
                 "STEP1 - REMOVE ILLUMINA ADAPTORS", cmd1, A1,
                 #"STEP1 - CUT 3 pb ON R2 PAIRED because of quality drop on Novaseq",cmd1, A1, # what to put in my file
                 append= F, sep = "\n\n")
             
           } ,
           mc.cores = 6
  )
  

# SOME POST ANALYSIS FILE MANIPULATION

# Rename the files

old.name <-list.files(get.value("trimmo.path"), full.names = T)
old.name

new.name <- old.name %>% 
  str_replace("_1P.fastq", "_R1.fastq") %>% 
  str_replace("_2P.fastq", "_R2.fastq") 
  
new.name

file.rename(from = old.name,
            to = new.name)

# Then remove ALL unecessary files

file.to.remove <- list.files(get.value("trimmo.path"), full.names = T, pattern = "1U.fastq.gz|2U.fastq.gz|HC3.fastq.gz")
sum(file.size(file.to.remove))/  1024 / 1024 / 1024

file.remove(file.to.remove)


# Compute stats by plate/RunSeq
# Will work for paired reads

trim.data <- data.frame(ID = character(),
                        RunSeq = character(),
                        No_plaque = character(),
                        Ntotal = numeric(),
                        Nsurvival = numeric(),
                        stringsAsFactors = F)

for(x in  list.files(get.value("trimmo.log"), pattern = ".log")){

  ID <- x %>% str_remove(".log")
  ID.int = ID %>% str_replace("[.][A-Z][:digit:][:digit:][:digit:]---[A-Z][:digit:][:digit:][:digit:][.]Sebaste_", "___")
  RunSeq = sapply(str_split(ID.int, "___"), `[`, 1)
  No_plaque = sapply(str_split(ID.int, "___"), `[`, 2) %>% str_remove("Plaque-")
  temp <- readLines(file.path(get.value("trimmo.log"), x ))
  temp <- temp %>% str_subset(pattern = "Input Read Pairs")
  Ntotal    <- sapply(str_split(temp, " "), `[`, 4)
  Nsurvival <- sapply(str_split(temp, " "), `[`, 7)

  trim.data <- bind_rows(trim.data,
                        data.frame(ID = ID,
                                   RunSeq = RunSeq,
                                   No_plaque = No_plaque,
                                   Ntotal = as.numeric(Ntotal),
                                   Nsurvival = as.numeric(Nsurvival),
                                   stringsAsFactors = F))
  
}

trim.data

write_csv(x = trim.data, file = file.path(get.value("trimmo.log"), "Summary_Nreads.csv") )

gg.trim <- trim.data %>% ggplot(aes(x = No_plaque, y = Nsurvival, fill = RunSeq)) +
           geom_bar(stat = "identity") +
           labs(title = "N reads by plates post trimmomatic") +
           theme_bw() + 
           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = perc_surv)) +
  geom_histogram() +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = Nsurvival, fill = No_plaque)) +
  geom_histogram() +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

gg.trim <- trim.data %>% mutate(perc_surv = Nsurvival/Ntotal) %>% 
  ggplot(aes(x = Ntotal, y = Nsurvival, col = No_plaque)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(title = "N reads by plates post trimmomatic") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.trim

#ggsave(filename = file.path(get.value("trimmo.log"), "Summary_Nreads.png"),
#       plot = gg.trim,
#       width = 5, height = 4, units = "in")


# If you want the log file to be ignore, run the following :

cat("*.log", "!.gitignore", sep = "\n",
    file = file.path(get.value("trimmo.log"), ".gitignore"))


# Run FastQC - change for overwrite T first the first time to or it will not work
# THIS PART SHOULD BE BENCHMARKED TOO

fastqc(folder.in = get.value("trimmo.path"), folder.out = get.value("fastqc.trim.path"), 
       overwrite = T, nthread = 12)

multiqc(get.value("fastqc.trim.path"))


# From the fastqc html file, I extracted the sequence length distribution

fastqc_length <-read_tsv("02_Results/00_Stacks/01_FastQC/02_Trimmomatic/MultiQC_report/fastqc_sequence_length_distribution_plot_R1.tsv")

names(fastqc_length)[1] <- "Sequence_length"


fastqc_length.tidy <- fastqc_length %>% pivot_longer(-c(Sequence_length), names_to = "Plate", values_to = "Nread") %>% 
  group_by(Plate) %>%  arrange(Plate) %>% 
  mutate(Total_read = sum(Nread),
  Cumulative_read = cumsum(Nread),
  Inv_cumulative_read = Total_read - Cumulative_read + Nread,
  P_read = Cumulative_read / Total_read * 100,
  Inv_p_read = 100 * Inv_cumulative_read / Total_read,
  N_nuc = Nread  * Sequence_length,
  Max_nuc = max(Sequence_length) * Total_read,
  Inv_Sequence_length = max(Sequence_length) - Sequence_length,
  Inv_cumulative_nucl = Inv_cumulative_read * Sequence_length,
  P_nucl = Inv_cumulative_nucl / max(Inv_cumulative_nucl) * 100
  #  Inv_cumulative_nucl = Max_nuc - Cumulative_nucl
  #P_nucl =  (N_nuc / Max_nuc * 100)
  ) 


fastqc_length.tidy %>% View()

fastqc_length.tidy %>%   ggplot(aes(x = Sequence_length, y = Inv_p_read, group = Plate)) +
  geom_line() 
  

gg.length <- fastqc_length.tidy %>% pivot_longer(-c(Sequence_length, Plate)) %>% 
  dplyr::filter(name %in% c("Inv_p_read", "P_nucl")) %>% 
  #ungroup() %>% 
  group_by(name, Sequence_length) %>% 
  summarise(mean = mean(value),
         sd = sd(value)) %>% 
  ggplot(aes(x = Sequence_length, y = mean, col = name)) +
  geom_hline(yintercept = c(80,85,90,95), col = "gray", lty = "dashed") +
  geom_vline(xintercept = c(100), col = "darkred", cex = 1.5) +
  geom_line()+ geom_point()+
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd)) +
  labs(y = "Mean %") +
  ggtitle("Relation between % reads and % pb to maximize ddRAD analysis")+
  #facet_grid(.~name) +
  theme_bw()

gg.length 
  
ggsave(filename = "02_Results/00_Stacks/02_Trimmomatic/Sequence_length.png", plot = gg.length)

fastqc_length.tidy %>%   ggplot(aes(x = Sequence_length, y = Inv_p_read, group = Plate)) +
  geom_line(alpha = 0.2) 

fastqc_length.tidy %>%   ggplot(aes(x = Sequence_length, y = P_nucl, group = Plate)) +
  geom_line(alpha = 0.2) 


# Demultiplex -------------------------------------------------------------

# Barcode files are important for thise step
# They couldbe created with the 01_Create_PopMap.R

# TAKE CARE, all individuals with similar names will be collapse
length(pop.info$Sample)
length(pop.info$Sample %>% unique())

# If you specify the same output directory for two differents
# process_radtags runs, the second run will overwrite identical filenames
# from the first run. Output the data into separate directories, and then
# concatenate the shared samples together after the runs complete

sub.dir <- c("NS.LH00487_0017.001")
sub.dir

# Before running, check that the right barcode file is found - here an example
# not in a loop

files.to.use <-  list.files(get.value("trimmo.path"), full.names = T, pattern = "_R1")# %>% 
 # str_subset(sub.dir[1])
files.to.use

barcode <- list.files(get.value("info.path"), full.names = T, pattern = "barcodes.txt") %>% 
  str_subset(paste0(files.to.use[1] %>% str_remove(get.value("trimmo.path")) %>% 
               str_remove("_R1.fastq.gz") %>% 
               str_remove(sub.dir[1]) %>%
               str_remove("[.][A-Z][:digit:][:digit:][:digit:]---[A-Z][:digit:][:digit:][:digit:][.]") %>% 
               str_remove("/"),
               "_barcodes.txt")) 
barcode



for(i in sub.dir){
  
  files.to.use <-  list.files(get.value("trimmo.path"), full.names = T, pattern = "_R1") %>% 
    str_subset(i) #%>% str_subset("P02|P03")
  
  # files.to.use
  
  cat(paste("\nWorking with the run:", i),
      paste("There are" , length(files.to.use) , "plates within this run"),
      sep= "\n")
  
  # Créer un sous-dossier s'il n'existe pas
  if(file.exists(file.path(get.value("demulti.path"), i)) == F) {
    
    cat(paste("\nCreating a new directory:", file.path(get.value("demulti.path"), i)),
        sep= "\n")
    
    dir.create(file.path(get.value("demulti.path"), i), recursive = T)}
  
  
# Parallel version of process_radtag 

  
mclapply(files.to.use,
         FUN = function(x){
           
           # Files
           file1 <- x
           file2 <- file1 %>% str_replace("_R1", "_R2")
           
           barcode <- list.files(get.value("info.path"), full.names = T, pattern = "barcodes.txt") %>% 
             str_subset(paste0(x %>% str_remove(get.value("trimmo.path")) %>% 
                          str_remove("_R1.fastq.gz") %>% 
                          str_remove(i) %>%
                            str_remove("[.][A-Z][:digit:][:digit:][:digit:]---[A-Z][:digit:][:digit:][:digit:][.]") %>% 
                            str_remove("/"),
                          "_barcodes.txt")) 
           
           plate <- barcode %>% str_remove(get.value("info.path")) %>% 
                                str_remove("_barcodes.txt") %>% 
                                str_remove("/")
           
           # Créer un sous-dossier s'il n'existe pas
           if(file.exists(file.path(get.value("demulti.path"), i, plate)) == F) {
             
             #cat(paste("\nCreating a new directory:", file.path(get.value("demulti.path"), i)),
             #     sep= "\n")
             
             dir.create(file.path(get.value("demulti.path"), i, plate), recursive = T)}
     
           # Command
           cmd <- paste("--paired",
                        "-1", file1,
                        "-2", file2,
                        "-o", file.path(get.value("demulti.path"), i, plate),
                        "--inline_null",   
                        "-b", barcode,
                        "--renz_1", "pstI", # Check RAD site on R1 
                        "--renz_2", "mspI", # CGG on R2 
                         "-E", "phred33",
                        "--filter-illumina",
                        "-c", # clean
                        "-r", #rescue
                        "-q", #check quality
                        "-t", 90,  #truncate at 150 - 22pb adaptors - 8 (max barcode) = 111 (all reads within sample must be the same length)
                        "-i", "gzfastq"
           )
           
           A <- system2("process_radtags", cmd, stdout=T, stderr=T)
           A
           # save a log file 
           log.file <- file1 %>% str_replace(get.value("trimmo.path"), get.value("demulti.log")) %>% 
             str_replace(".fastq.gz","_summary.log") %>% 
             str_remove("_R1")
           
           cat(file = log.file,
               cmd, "\n\n",
               A, # what to put in my file
               append= F, sep = "\n")
           
           # Detailed log file
           file.rename(from = list.files(file.path(get.value("demulti.path"), i, plate), full.names = T, pattern = ".log"),
                       to = log.file %>% str_replace("summary", "detailed")
                       )
           
         } ,
         mc.cores = 15
)

gc()

}


# If you want the log/csv file to be ignore, run the following :

cat("*", "!.gitignore", "!*.png", "!AllIndividuals_Nreads.csv", sep = "\n",
    file = file.path(get.value("demulti.log"), ".gitignore"))


# Extract data for each summary_detailed
# You must change the str_subset code that is the index for each project

for(x in list.files( get.value("demulti.log"), pattern = "detailed", full.names = T)){
  
  data <- readLines(x) %>% 
    str_subset("S_") %>% 
    str_split(pattern = "\t")
  
  cat("\n",length(data), "samples retrieved in", x %>% str_remove(get.value("demulti.log")) %>% str_remove("/"))
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  # If there is a RUN column, it's because some individuals come from more than one sample
  names(data) <- c("Barcode", "Filename", "Total", "NoRadTag", "LowQuality", "Retained", "Pct_Retained", "Pct_Total_Reads")
  
  data <- data %>% mutate(Total = as.numeric(as.character(Total)),
                          NoRadTag = as.numeric(as.character(NoRadTag)),
                          LowQuality = as.numeric(as.character(LowQuality)),
                          Retained = as.numeric(as.character(Retained)),
                          Run = x %>% str_remove(get.value("demulti.log")) %>% 
                            str_remove("_log_detailed.txt")  %>% 
                            str_remove("/")
  )
  write_csv(data, x %>% str_replace("detailed.log", "Nreads.csv"))
  
}

# Create one big log file 

Nreads.data <- data.frame()

for(x in list.files( get.value("demulti.log"), pattern = "Nreads", full.names = T) %>% str_subset(".csv") %>% str_subset("NS.")){
  data.int <- read_csv(x)
  
  Nreads.data <- bind_rows(Nreads.data, data.int)
  
}

# Compute N read removed
Nreads.data$Removed <- Nreads.data$Total - Nreads.data$Retained  

# Add a column for the run name
Nreads.data$RunSeq <- paste(sapply(str_split(Nreads.data$Run, "[.]"),`[`,1),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,2),
                            sapply(str_split(Nreads.data$Run, "[.]"),`[`,3),
                            sep = ".")


# Save the result
 write_csv(Nreads.data, file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))

Nreads.data <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))

trim.data <- read_csv(file = file.path(get.value("trimmo.log"), "Summary_Nreads.csv") )

Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ")) 

# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)
Nreads.data$Filename %>% unique() %>% length()

# Stats by plate
Nreads.data %>% group_by(No_plaque_envoi) %>% summarise(Retained = sum(Retained)) %>%
  group_by() %>% 
  summarise(mean = mean(Retained)/2,
            sd = sd(Retained)/2,
            min = min(Retained)/2,
            max = max(Retained)/2)

# Impact of overall process

gg.radtag <- Nreads.data %>% group_by(No_plaque_envoi, RunSeq) %>% 
                #mutate(RunSeq = paste()) 
                summarise(Total = sum(Total)/2,
                          Retained = sum(Retained)/2) %>% 
                left_join(trim.data %>% select(RunSeq, No_plaque_envoi = No_plaque, Nsurvival)) %>% 
                mutate(perc_with_barcode = Total/Nsurvival,
                       perc_keep = Retained / Nsurvival,
                       perc_keep_over_barcode = Retained / Total) %>% 
                select(RunSeq, No_plaque_envoi, perc_with_barcode, perc_keep,  perc_keep_over_barcode) %>% 
                pivot_longer(c(perc_with_barcode, perc_keep,  perc_keep_over_barcode), names_to = "Cat", values_to = "Perc") %>% 
  ggplot(aes(x = No_plaque_envoi, y = Perc, col = RunSeq, group = RunSeq)) +
  geom_point(position= position_dodge(width = 1/4)) +
  labs(title = "% reads retained post demultiplexing") +
  facet_grid(Cat ~ . ) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

gg.radtag

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Preads_byPlates.png"),
       plot = gg.radtag,
       width = 8, height = 6, units = "in")


gg.retained.plate <- Nreads.data %>% #gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(RunSeq, No_plaque_envoi) %>% summarise(Retained = sum(Retained)/2,
                                            Removed = sum(Removed)/2) %>% 
  ggplot(aes(y = Retained, x = No_plaque_envoi, fill = RunSeq)) + 
  geom_bar(stat = "identity")+
  theme_bw() + 
  labs(title = "N reads retained post demultiplexing") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
gg.retained.plate

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byPlates.png"),
       plot = gg.retained.plate,
       width = 5, height = 4, units = "in")


# Check that we have all the data - all with unique data

head(Nreads.data)
Nreads.data %>% nrow() / length(sub.dir)
Nreads.data$Filename %>% unique() %>% length()

# Stats by IND
Nreads.data %>% group_by(Filename) %>% summarise(Retained = sum(Retained)) %>%
                group_by() %>% 
                summarise(mean = mean(Retained)/2,
                          sd = sd(Retained)/2,
                          min = min(Retained)/2,
                          max = max(Retained)/2)

gg.retained.ind <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(Filename, No_plaque_envoi, Cat) %>% summarise(N = sum(N)/2) %>% 
  ggplot(aes(x = Filename, y = N, fill= Cat)) + 
  geom_bar(stat= "identity") +
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(y = "N reads", x = "Samples", title = "N reads by ind") +
  facet_wrap(~ No_plaque_envoi, scale = "free_x") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")
gg.retained.ind

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byInd.png"),
       plot = gg.retained.ind,
       width = 8, height = 8, units = "in")

gg.retained.pop <- Nreads.data %>% gather(Removed, Retained, key = "Cat", value = "N") %>% #head()
  group_by(Filename, Site, Shore, Region, Espece) %>% summarise(N = sum(N)/2) %>% 
  ggplot(aes(x = Site, y = N)) + 
  geom_boxplot(fill = NA) +
  geom_jitter(height = 0, aes(col = Espece))+
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000),
  #                   labels = c(1:7))+
  labs(y = "N reads", x = "Samples", title = "N reads by ind") +
  facet_grid(~ Region, scale = "free_x", space = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle  = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")
gg.retained.pop 

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byPop.png"),
       plot = gg.retained.pop,
       width = 8, height = 6, units = "in")


gg.retained.ind2.a <- Nreads.data %>% 
  mutate(barcode_length = nchar(BC),
         first.nuc = str_sub(BC, 1,1)) %>% 
  group_by(Filename, barcode_length, first.nuc) %>% summarise(Retained = sum(Retained)/2) %>%
  #  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained)) +
  scale_y_continuous(trans = "log10") +
  geom_bar(stat = "identity") +
  #facet_grid(barcode_length ~ first.nuc, scale = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())
gg.retained.ind2.a


gg.retained.ind2.b <- Nreads.data %>% 
  mutate(barcode_length = nchar(BC),
         first.nuc = str_sub(BC, 1,1)) %>% 
  group_by(Filename, barcode_length, first.nuc) %>% summarise(Retained = sum(Retained)/2) %>%
  #  arrange(Retained) %>% 
  ggplot(aes(x = reorder(Filename, Retained), y = Retained)) +
  scale_y_continuous(trans = "log10") +
  geom_bar(stat = "identity") +
  facet_grid(barcode_length ~ first.nuc, scale = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())
gg.retained.ind2.b


gg.retained.ind2 <- ggpubr::ggarrange(gg.retained.ind2.a, gg.retained.ind2.b,
                                      nrow = 2, heights = c(1:3))
gg.retained.ind2

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byInd_withBarcodes.png"),
       plot = gg.retained.ind2,
       width = 8, height = 8, units = "in")


gg.barcodes <- Nreads.data %>% 
   group_by(No_plaque_envoi, BC, Filename,RunSeq) %>% summarise(Total = sum(Total)/2) %>% 
  ggplot(aes(x = BC, y = Total)) + 
  geom_boxplot() +
  geom_jitter(aes(col =RunSeq), cex = 1, alpha = 0.5) +
 # facet_wrap(~Barcode.y)
  #scale_y_continuous(breaks = c(1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000),
  #                   labels = c(1:10))+
  scale_y_continuous(trans = "log10") +
  labs(y = "N total reads (log)", x = "Barcodes") +
 # facet_grid(. ~ Espece, scale = "free_x", space = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "bottom")
gg.barcodes

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode.png"),
       plot = gg.barcodes,
       width =8, height = 5, units = "in")

# This is THE graph
gg.barcodes.summary <- Nreads.data %>% 
  group_by(BC, Filename) %>% summarise(Total =mean(Total)/2, 
                                              Retained =mean(Retained)/2,
                                              NoRadTag = mean(NoRadTag)/2,
                                              LowQuality = mean(LowQuality)/2) %>%
  pivot_longer(c(Total, Retained, NoRadTag, LowQuality), names_to = "Cat", values_to = "N") %>% 
  mutate(barcode_length = nchar(BC),
         first.nuc = str_sub(BC, 1,1)) %>% 
  ggplot(aes(x = barcode_length, y =N, group = barcode_length) )+ 
  geom_violin() +
  geom_jitter(alpha = 0.5, aes(col = first.nuc), cex = 1) +
  #geom_violin()+
  # facet_wrap(~Barcode.y)
  scale_y_continuous(trans = "log10") +
  labs(y = "N reads (log)", x = "Barcodes length (4-8)") +
  facet_grid(Cat ~ ., scale = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ),
        legend.position = "bottom")

gg.barcodes.summary

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode_Summary.png"),
       plot = gg.barcodes.summary,
       width =6, height = 6, units = "in")

gg.barcodes.summary2 <- Nreads.data %>% 
  group_by(BC, Filename) %>% summarise(Total =mean(Total)/2, 
                                              Retained =mean(Retained)/2,
                                              NoRadTag = mean(NoRadTag)/2,
                                              LowQuality = mean(LowQuality)/2) %>%
  pivot_longer(c(Total, Retained, NoRadTag, LowQuality), names_to = "Cat", values_to = "N") %>% 
  mutate(barcode_length = nchar(BC),
         first.nuc = str_sub(BC, 1,1)) %>% 
  group_by(first.nuc, barcode_length) %>% summarise(Nretained = mean(N[Cat == "Retained"])) %>% 
  arrange(Nretained) %>% 
  ggplot(aes(x = factor(barcode_length), y = first.nuc, fill = Nretained)) +
  geom_bin2d() + theme_bw()
gg.barcodes.summary2

ggsave(filename = file.path(get.value("demulti.log"), "Summary_Nreads_byBarcode_Summary2.png"),
       plot = gg.barcodes.summary2,
       width =5, height = 4, units = "in")


## Create one file with up to 3 files by individual
## UPDATE: sould work with up to 3 subdir 
sub.dir <- c("NS.LH00487_0017.001")
sub.dir
plates <- file.path(get.value("demulti.path"), sub.dir[1]) %>% list.files()
plates

for(p in plates){
  
  files.to.use <- list.files(file.path(get.value("demulti.path"), sub.dir[1], p), full.names = T)
  
  cat(paste("\nWorking with the plate:", p),
      paste("There are" , length(files.to.use) , "files to process"),
      sep= "\n")
  
  mclapply(files.to.use,
           FUN = function(x){
             
             # Files
             file1 <-  x
             
             if(length(sub.dir) >= 2){
               file2 <- file1 %>% str_replace(sub.dir[1], sub.dir[2])  
             }
             
             if(length(sub.dir) >= 3){
               file3 <- file1 %>% str_replace(sub.dir[1], sub.dir[3])  
             }
             
             file.join <- file1 %>% str_remove(file.path(sub.dir[1], p))# %>% 
             #str_remove("[:digit:][:digit:][:digit:][:digit:][.][:digit:][:digit:][:digit:][.][A-Z][:digit:][:digit:][:digit:][.]Panomics-" )
             
             if(length(sub.dir) == 1){
               
               file.rename(from = file1, to = file.join) 
               
             }
             
             if(length(sub.dir) == 2){
               # Cat command - this is so simple !
               cmd <- paste(file1, file2, ">", file.join)
               system2("cat", cmd, stdout=T, stderr=T)
             }
             
             if(length(sub.dir) == 3){
               # Cat command - this is so simple !
               cmd <- paste(file1, file2, file3, ">", file.join)
               system2("cat", cmd, stdout=T, stderr=T)
             }
             
             
           } ,
           mc.cores = 20
  )
  
  gc()
  
}

# Number of files observed
list.files(file.path(get.value("demulti.path")), pattern = ".fq.gz") %>% length()

# Number of files expected
384 * length(plates)

1152/384

# fastqc(folder.in = file.path("00_Data/03a_Demultiplex/NS.1754.001/Plaque-1/"), folder.out = get.value("fastqc.demulti.path"), 
#        overwrite = T, nthread = 20, fastqc = "fq")
# 
# multiqc(get.value("fastqc.demulti.path"))

# Select individuals for tests --------------------------------------------

# Select between 10-20 representative individuals
# Will be used both for testing alignement and stacks

Nreads.data <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))
Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ")) 

Nreads.data %>% group_by(Espece, Cat_sample) %>% summarise(N = n())

Nreads.data %>% group_by(Espece, Cat_sample, Region_echantillonnage) %>% summarise(N = n())

Nreads.data %>% nrow()

Nreads.data  %>% ggplot(aes(x = Region_echantillonnage, y = Retained, col =  Espece)) + 
  geom_boxplot() +
  geom_jitter(height = 0) +
  facet_grid(.~Espece, space ="free", scale = "free")+
  theme_bw()

Nreads.data  %>% dplyr::filter(Espece == "undatum") %>% 
  ggplot(aes(x = Numero_unique_groupe, y = Retained, col =  Region_echantillonnage)) + 
  geom_boxplot() +
  geom_jitter(height = 0) +
  geom_hline(yintercept =  1000000, lty = "dashed")+
  facet_grid(.~Espece, space ="free", scale = "free")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5 ))


# Create the list of individuals for testing parameters :

hist(Nreads.data$Retained)
#summary(Nreads.data$Retained.by.sample)
summary(Nreads.data)
Nreads.data %>% pull(Cat_sample) %>% unique()

test.ID <- Nreads.data %>% filter(Cat_sample == "Sample",
                                  Espece== "undatum") %>% 
            mutate(Pop = paste(Espece, Region_echantillonnage, Mois_echantillonnage, sep = "_")) %>% 
           group_by(Filename, Pop) %>% 
           summarise(Retained = sum(Retained)) %>% 
  filter(Retained >= quantile(Retained, probs = 0.25),
         Retained <= quantile(Retained, probs = 0.75)) %>% 
  group_by(Pop) %>% #summarise((N = n()))
  sample_n(2) %>% pull(Filename)  

test.ID 

#write.table(pop.info %>% select(Sample, POP) %>% 
#                          filter(Sample %in% test.ID) %>% 
#                          mutate(POP = "noPOP"), 
#             file = file.path(get.value("info.path"), "popmap.test_samples.txt"),
#             quote = FALSE, sep = "\t",
#             row.names = F, col.names = F)

# Test ID 2

previous.test.ID <- read.table(file.path(get.value("info.path"), "popmap.test_samples.txt"))
names(previous.test.ID) <- c("Sample", "POP")


test.SP <- Nreads.data %>% filter(Cat_sample == "Sample",
                                  Espece!= "undatum") %>% 
  mutate(Pop = paste(Espece, Region_echantillonnage, Mois_echantillonnage, sep = "_")) %>% 
  group_by(Filename, Pop) %>% 
  summarise(Retained = sum(Retained)) %>% 
  filter(Retained >= quantile(Retained, probs = 0.25),
         Retained <= quantile(Retained, probs = 0.75)) %>% 
  #group_by(Pop) %>% #summarise((N = n()))
  #sample_n(2) %>% 
  pull(Filename)  

test.SP


#write.table(rbind(previous.test.ID,
#              pop.info %>% select(Sample, POP) %>% 
#              filter(Sample %in% test.SP) %>% 
#              mutate(POP = "noPOP")), 
#            file = file.path(get.value("info.path"), "popmap.test_species.txt"),
#            quote = FALSE, sep = "\t",
#            row.names = F, col.names = F)

# Redone, but with individual from only one specie

test.ID <- Nreads.data %>% filter(Cat_sample == "Sample",
                                  Espece== "undatum",
                                  Filename %in% GOOD.SP) %>% 
  mutate(Pop = paste(Espece, Region_echantillonnage, Mois_echantillonnage, sep = "_")) %>% 
  group_by(Filename, Pop) %>% 
  summarise(Retained = sum(Retained)) %>% 
  #filter(Retained >= quantile(Retained, probs = 0.25),
  #       Retained <= quantile(Retained, probs = 0.75)) %>% 
  group_by(Pop) %>% #summarise((N = n()))
  sample_n(2) %>% pull(Filename)  

test.ID 

#write.table(pop.info %>% select(Sample, POP) %>% 
#            filter(Sample %in% test.ID) %>% 
#            mutate(POP = "noPOP"), 
#            file = file.path(get.value("info.path"), "popmap.test_undatum.txt"),
#            quote = FALSE, sep = "\t",
#            row.names = F, col.names = F)
 
 # BWA - index the reference genome ----------------------------------------
 
 # Check that 
 A <- system2("bwa", "", stdout=T, stderr=T)
 A
 list.files(get.value("ref.genome.path"))
 
 # Babylonia_areolata
 
 cmd <- paste("index",
              "-p",  file.path(get.value("ref.genome.path"),"Babylonia_areolata_GCA_011634625.1_BA9_1.0_genomic"), 
              "-a", "bwtsw",
              file.path(get.value("ref.genome.path"),"Babylonia_areolata_GCA_011634625.1_BA9_1.0_genomic.fna")
 )
 
 
 A <- system2("bwa", cmd, stdout=T, stderr=T)
 A
 # save a log file 
 
 cat(file = file.path(get.value("ref.genome.path"), "Babylonia_areolata.Feb2024.log" ),
     cmd, "\n\n",
     A, # what to put in my file
     append= T, sep = "\n")

 
 # BWA - Align reads to the reference genome -------------------------------
 
 # This part take times, can start on a subset (and test Stack at the same times)
 # 2-3 min by file
 
 TEST.ID <- read.table( file.path(get.value("info.path"), "popmap.test_samples.txt")) %>% pull(V1) %>% paste(collapse = "|")
 TEST.ID <- read.table( file.path(get.value("info.path"), "popmap.test_species.txt")) %>% pull(V1) %>% paste(collapse = "|")
 
 # ATTENTION - ICI POUR PAIRED 
 # Test version
 demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>%
                  str_subset(".rem.", negate = T) %>%  
                  str_subset(TEST.ID)
 # Complete version
 demulti.files <- list.files(get.value("demulti.path"), pattern = ".1.fq.gz", full.names = T) %>% 
   str_subset(".rem.", negate = T)# %>%  str_subset(TEST.ID, negate = T)

demulti.files %>% length()
  

 #demulti.files[1:10]
 
  mclapply(demulti.files,
          FUN = function(x){
          # How I will rename all this : 
          file.R1 <- x
          file.R2 <- x %>% str_replace(".1.fq.gz", ".2.fq.gz")
          file.bam <- x %>% str_replace(get.value("demulti.path"), get.value("align.path")) %>% 
                            str_replace(".1.fq.gz", ".bam")
          file.sort.bam <- file.bam %>% str_replace(".bam", ".sorted.bam")
          stat.tsv <- file.bam %>% str_replace(".bam", ".stat.tsv")
          
          if(!file.exists(stat.tsv)){
          
          # DO THE ALIGMENT  
          cmd1 <- paste("mem",
                        "-t", 16,
                        "-M",

                        file.path(get.value("ref.genome.path"),"Babylonia_areolata_GCA_011634625.1_BA9_1.0_genomic"), # the index ref genome
                        file.R1,
                        file.R2,
                        "2> /dev/null",
                        "| samtools", "view", "-Sb", 
                        "--threads", 15, # Number of additional threads
                        "-o", file.bam#,
          )# the file
          
          A <- system2("bwa", cmd1, stdout=T, stderr=T)
          
          # Sort the file 
          cmd2 <- paste("sort",
                        "--threads", 15,
                        #"-n",
                        "-O", "BAM",
                        # the new file (sorted)
                        "-o", file.sort.bam,
                        # the bam file to sort
                        file.bam
                        
          )
          
          A <- system2("samtools", cmd2, stdout=T, stderr=T)
          
          # Compute stats
          
          cmd3 <- paste("flagstat",
                        "--threads", 15,
                        "-O", "tsv",
                        file.sort.bam,
                        ">",
                        stat.tsv
          )
          
          A <- system2("samtools", cmd3, stdout=T, stderr=T)
          }
          
          },
          mc.cores = 2
 ) 
 
  
  
  # Remove unsorted files that are unecessary 
  
  files.to.remove <- list.files(get.value("align.path"), pattern = ".bam", full.names = T)  %>% 
    str_subset(pattern = "sorted", negate = T)
  
  length(files.to.remove)
  files.to.remove[1:20]
  
  for(x in files.to.remove){
    
    file.remove(x)
  }
  
  
  # If you want the big files to be ignore, run the following :
  
  cat("*.bam", "*.tsv", "!.gitignore", sep = "\n",
      file = file.path(get.value("align.path"), ".gitignore"))
  
  

 # Compute the aligned reads -----------------------------------------------

 map.res <- data.frame(ID = character(),
                       total = numeric(),
                       secondary = numeric(),
                       supplementary = numeric(),
                       duplicates = numeric(),
                       mapped = numeric(),
                       mapped_perc = numeric(),
                       paired = numeric(),
                       read1 = numeric(),
                       read2 = numeric(),
                       properly_paired = numeric(),
                       properly_paired_perc = numeric(),
                       twith_itself_mate_mapped = numeric(),
                       singletons = numeric(),
                       singletons_perc = numeric(),
                       twith_mate_mapped_diff_chr = numeric(),
                       twith_mate_mapped_diff_chr_HMQ = numeric(),
                       #Nmappedprim = numeric(),
                        stringsAsFactors = F)
 
 files.to.use <- list.files(get.value("align.path"), pattern = ".stat.tsv", full.names = T)  
  
for(x in seq_along(files.to.use)){
  
   ID <-   files.to.use[x] %>% str_remove(get.value("align.path")) %>% 
                               str_remove(".stat.tsv") %>% 
                               str_remove("/")
  
   temp <-  sapply(str_split(readLines(files.to.use[x]), "\t"), `[`,1)  
   map.res[x,] <- c(ID, temp %>% str_remove("%"))
   
   
}  

for(x in 2:ncol(map.res)){
  map.res[,x] <- as.numeric(as.character(map.res[,x]))
  
}  
    
nrow(map.res)
head(map.res)  
summary(map.res)  
map.res <-  map.res %>% left_join(pop.data, by = c("ID" = "ID_GQ"))

graph1.0 <- map.res %>% #filter(ID != "Dp_2066") %>% 
  ggplot(aes(x = mapped_perc, col = Espece)) + 
  #scale_x_continuous(limits= c(97,100)) +
  #geom_boxplot(col = "black") +
  #geom_jitter(width = 0, height = 0.25, cex = 2) +
  geom_density(alpha = 0.5) +
  #  geom_vline(xintercept = 96, lty = "dashed") +
  #facet_grid(Espece ~ . , scales = "free")
  labs(x = "Percentage of reads mapped", y = "Density") +
  theme_bw()  

graph1.0

ggsave(filename = file.path(get.value("demulti.log"), "Alignment_density.png"),
       plot = graph1.0,
       width =5, height = 4, units = "in")

graph1.1 <-  map.res %>%  ggplot(aes(y = mapped_perc, x = Espece)) + 
   geom_boxplot() +
   geom_jitter(aes(col = factor(Annee_echantillonnage)),height = 0, alpha = 0.5) +
   #facet_grid(Espece ~ . , scales = "free")
   labs(y = "Percentage of reads mapped") +
  geom_hline(yintercept = 96, lty = "dashed") +
     theme_bw() + 
   theme(legend.position = "bottom",
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

graph1.1 

ggsave(filename = file.path(get.value("demulti.log"), "Alignment_boxplot.png"),
       plot = graph1.1,
       width =5, height = 4, units = "in")
 
graph1.2 <- map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
   ggplot(aes(x=total, y=mapped_perc)) +
   geom_point(alpha = 0.5)+
  geom_hline(yintercept = 96, lty = "dashed") +
  #scale_x_continuous(trans = "log10") + 
  labs(y = "Percentage of reads mapped", x = "N read total") +
   theme_bw() 
 
graph1.2

map.res %>%  ggplot(aes(x = mapped, fill = Espece)) + 
  geom_histogram()+
  #facet_grid(Espece ~ . , scales = "free")
  geom_vline(xintercept = 100000, lty = "dashed") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

map.res %>%  ggplot(aes(x = mapped, y = total, col = mapped_perc)) + 
  geom_point()+
  #facet_grid(Espece ~ . , scales = "free")
  geom_vline(xintercept = c(100000), lty = "dashed") +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

graph1.3 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
   ggplot(aes(y = mapped_perc, x = Espece)) + 
   geom_boxplot(alpha = 0.5) +
   geom_hline(yintercept = 96, lty = "dashed") +
   geom_jitter(aes(col = total),height = 0, alpha = 0.75) +
   scale_color_distiller(palette = "Spectral", trans = "log10") +
   #facet_wrap(~ Gen_ZONE)
   labs(y = "Percentage of reads mapped") +
   theme_bw()  +
   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

graph1.3 
 

graph1.4 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped, col = No_plaque_envoi)) +
  geom_density() +
  #geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
  #           lty = "dashed" )+
  #facet_wrap(~ Gen_ZONE)+
  labs(x = "N reads mapped") +
  theme_bw()
graph1.4

graph1.5 <- map.res %>% #filter(Espece_revision %in% c("Pb"))  %>% 
  ggplot(aes(x = mapped)) +
  geom_density() +
  #geom_vline(xintercept = map.res %>% filter(Espece_revision %in% c("Pb")) %>% pull(mapped) %>% median, 
  #           lty = "dashed" )+
  labs(x = "N reads mapped") +
  facet_wrap(~ Espece)+
  theme_bw()
graph1.5
 

## Overall popmap

# Select those that will be used based on alignment
write.table(map.res %>% filter(mapped >= 100000) %>% select(ID) %>% 
                               mutate(Pop = "NoPop"), 
            file = file.path(get.value("info.path"), "popmap_good_align.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

 map.res %>% #group_by(Espece_revision) %>% 
             summarise(MedianNread = median(total),
                       MedianPerc = median(mapped_perc),
                       MinPerc = min(mapped_perc),
                       MaxPerc = max(mapped_perc))
 
plot(map.res[,c("total", "secondary", "mapped", "paired", "properly_paired", "singletons", "twith_itself_mate_mapped", "twith_mate_mapped_diff_chr")])

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=secondary, col=factor(No_plaque))) +
  geom_smooth()+
  #geom_point() +
  #facet_wrap(~Espece_revision)+
  theme_bw()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=secondary, y=singletons, col=factor(No_plaque))) +
  #geom_smooth()+
  geom_point() +
 # facet_wrap(~Espece_revision)+
  theme_bw()


map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=total, y=singletons, col=factor(No_plaque))) +
  geom_point()

map.res %>% #filter(Espece_revision %in% c("Pb", "Pm")) %>% 
  ggplot(aes(x=twith_itself_mate_mapped, y=singletons, col=factor(No_plaque))) +
  geom_point()


hist(map.res$total) 
 
quantile(map.res$mapped, c(0.01,0.05,0.1, 0.5))# %>% summary()
  
# map.res %>%  %>% summarise(Mean = mean(Permapped))
 
head(map.res)

#write_csv(map.res, file = file.path(get.value("demulti.log"), "Mapping_Results.csv"))
 

# Stacks - Testing parameters - Ref-genome --------------------------------



# This part will run the ref_map.pl pipeline with a subset of samples,
# just to be sure that everything is alright

cmd <- paste("--samples", get.value("align.path"),
             "--popmap", file.path(get.value("info.path"), "popmap.test_species.txt"),
             "-o", get.value("test.ref.path"),

             "-T", 20,
             "-X", "\"populations:-r 0.80 --vcf --write-single-snp\"",
             "-X", "\"gstacks:-S .sorted.bam\"" # en espérant que ça fonctionne ...
)

A <- system2("ref_map.pl", cmd, stdout=T, stderr=T)
A

# Stats

# Extract statistics

  cmd1 <- paste(file.path(get.value( "test.ref.path"), "populations.log.distribs"),
               "samples_per_loc_postfilters")
  
  res1 <- system2("stacks-dist-extract", cmd1, stdout = T)
  
  data.int <- res1[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
  names(data.int) <- c("n_samples", "n_loci")
  
  n_loci <- as.numeric(as.character(data.int$n_loci)) %>%  sum() 
  
  # N SNP / locus
  
  cmd2 <- paste(file.path(get.value( "test.ref.path"), "populations.log.distribs"),
               "snps_per_loc_postfilters")
  
  res2 <- system2("stacks-dist-extract", cmd2, stdout = T)
  
  data.int <- res2[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data.int <-  data.frame(matrix(unlist(data.int), nrow= length(data.int), byrow = T))
  names(data.int) <- c("n_snps", "n_loci")
  
  data.int$n_snps <- as.numeric(as.character(data.int$n_snps))
  data.int$n_loci <- as.numeric(as.character(data.int$n_loci))  
 
  no_div <- data.int %>% filter(n_snps == 0) %>% pull(n_loci)  %>% as.numeric()
  n_loci_poly <- n_loci - no_div

 cat("\nThere is", n_loci_poly, "polymorphic loci (r80) out of", n_loci, "loci")

 
vcf.path <- ("00_Data/04b_Test.ref/populations.snps.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)
 
library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
                left_join(pop.data) %>% pull(Espece)

table(pop(gl.data))

plot(gl.data)



  res <- apply(tab(gl.data,NA.method = c("asis") ), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })

pop.data %>% dplyr::filter(ID_GQ %in% names( res[res>0.5])) %>% dplyr::select(ID_GQ, Espece)
  
    
#library(remotes)
#remotes::install_github("biodray/QuickPop")

pca.test  <- adegenet::glPca(gl.data[indNames(gl.data) %nin%  names( res[res>0.30]),], center = TRUE, scale = FALSE,  
                   parallel = TRUE, nf = 10)

gPCA <- pca.test %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col =Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
#  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
 # facet_wrap(~Espece) +
  #scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data)), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA

ggsave(filename = file.path(get.value("test.ref.log"), "PCA_test_StackRef.png"),
       plot = gPCA,
       width = 5, height = 4, units = "in")

# If you want the big files to be ignore, run the following :
cat("catalog.*", "*.tsv", "*.vcf", "!.gitignore", sep = "\n",
     file = file.path(get.value("test.ref.path"), ".gitignore"))
 

# Gstakcs - Ref - discovered SNPs -----------------------------------------

cmd <- paste("-I", get.value("align.path"),
             "-M", file.path(get.value("info.path"), "popmap_good_align.txt"),
             "-O", get.value("stacks.ref.path"),
             "-S", ".sorted.bam",
             "-t", 8)

A <- system2("gstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.ref.log"), "gstacks.ref.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")



# Check the distribution of ...


# N SNP / locus

cmd <- paste(file.path(get.value("stacks.ref.path"), "gstacks.log.distribs"),
             "effective_coverages_per_sample")

res <- system2("stacks-dist-extract", cmd, stdout = T)

data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
names(data) <- c("sample", "n_loci", "n_used_fw_reads", "mean_cov", "mean_cov_ns")


data <- data %>% left_join(pop.data, by = c("sample" = "ID_GQ"))

data %>% View()


data %>% group_by(Annee_echantillonnage) %>% summarise(N = n())
data %>% group_by(Region_echantillonnage) %>% summarise(N = n())

graph2.1 <- data %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), fill = Espece)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  facet_wrap(~Region_echantillonnage, nrow = 2) +
  labs(x = "Mean coverage") + 
  theme_bw() +
  theme(legend.position = "bottom")

graph2.1

ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageByRegion_May2022.png"),
       plot = graph2.1,
       width = 7.5, height = 4, units = "in")

graph2.2 <- data %>%# left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  ggplot(aes(x = as.numeric(as.character(mean_cov_ns)), y = as.numeric(as.character(n_loci)), col = as.factor(Espece))) +
  geom_point(alpha = 0.5) +
 # scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  geom_vline(xintercept = 20, lty = "dashed", col = "darkgray") +
  #geom_hline(yintercept = 60000, lty = "dashed", col = "darkgray") +
#  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Mean coverage", y = "N loci") + 
  theme_bw() +
  theme(legend.position = "none")
graph2.2

ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageVSloci_March2024.png"),
       plot = graph2.2,
       width = 7, height = 4, units = "in")

graph2.3 <- data %>% #left_join(pop.data, by = c("sample" = "ID_GQ")) %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = as.factor(Annee_echantillonnage), col = as.factor(Espece))) +
  geom_boxplot(col = "black") +  
  geom_jitter(height = 0, alpha = 0.5) +
  #scale_y_continuous(breaks = c(0, 60000, 120000, 180000, 240000, 300000, 360000, 420000, 480000)) +
  #geom_vline(xintercept = 5, lty = "dashed", col = "darkgray") +
  geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(y = "Mean coverage", y = "") + 
  theme_bw() +
  theme(legend.position = "none")

graph2.3

ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageVsYear_May2022.png"),
       plot = graph2.3,
       width = 7, height = 4, units = "in")


graph2.4 <- data %>% left_join(map.res %>% select(ID, total), by = c("sample" = "ID")) %>% 
  ggplot(aes(y = as.numeric(as.character(mean_cov_ns)), x = as.numeric(as.character(total)), col = as.factor(Espece))) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  geom_hline(yintercept = 5, lty = "dashed", col = "darkgray") +
  #geom_hline(yintercept = 60000, lty = "dashed", col = "darkgray") +
  #  facet_wrap(~Gen_ZONE, nrow = 2) +
  labs(x = "Nreads obtained", y = "Coverage") + 
  theme_bw() +
  theme(legend.position = "none")
graph2.4

ggsave(filename = file.path(get.value("stacks.ref.log"), "CoverageVsNreads_May2022.png"),
       plot = graph2.4,
       width = 7, height = 4, units = "in")


data %>% mutate(cat_cov = ifelse( as.numeric(as.character(mean_cov_ns)) <3, "LOW", "OK")) %>% group_by(cat_cov) %>% summarise(N = n())

# Save the results

readr::write_csv(data, file.path(get.value("stacks.ref.log"), "AllIndividuals_RefGen_NreadsNloci.csv"))

# Do a quick filtration equivalent r80


ID.coverage <- data %>% dplyr::filter(mean_cov_ns >= 20) %>% 
  mutate(Pop = "NoPop") %>% 
  select(sample, Pop) 
nrow(ID.coverage)

write.table(ID.coverage, 
            file = file.path(get.value("stacks.ref.path"), "popmap.coverage_20x.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

cmd <- paste("-P", get.value("stacks.ref.path"), 
             "-M", file.path(get.value("stacks.ref.path"), "popmap.coverage_20x.txt"),
             "--out-path",  file.path(get.value("stacks.ref.path")),
             "-t", 8,
             "-r", 0.8, #              
            # "--min-maf", maf.value,
             #"--smooth",
             "--write-single-snp",
             "--vcf"
)

A <- system2("populations", cmd, stdout=T, stderr=T)
A


vcf.path <- ("00_Data/05b_Stacks.ref/populations.snps.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(pop.data) %>% pull(Espece)

pop(gi.data) <- data.frame(ID_GQ = indNames(gi.data)) %>% 
  left_join(pop.data) %>% pull(Espece)


table(pop(gl.data))

plot(gl.data)


res <- apply(tab(gl.data,NA.method = c("asis") ), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
freq.na <- n.na / length(l)
return(freq.na)
})

pop.data %>% dplyr::filter(ID_GQ %in% names( res[res>0.90])) %>% dplyr::select(ID_GQ, Espece)

hist(res)

#library(remotes)
#remotes::install_github("biodray/QuickPop")

pca.test  <- adegenet::glPca(gl.data[indNames(gl.data) %nin%  names( res[res>0.50]),], center = TRUE, scale = FALSE,  
                             parallel = TRUE, nf = 10)

gPCA <- pca.test %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
#  left_join(df.g5, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
   facet_wrap(~Espece) +
  #scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data)), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA

ggsave(filename = file.path(get.value("test.ref.log"), "PCA_test_StackRef_all.png"),
       plot = gPCA,
       width = 5, height = 4, units = "in")


pop.data$Espece %>% table()

pca.test %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter(score.PC1 > -1) %>% View()


# library(adegenet)

Buccinum.clust.k4    <- snapclust(gi.data, k = 4, hybrids = F)
Buccinum.clust.k5    <- snapclust(gi.data, k = 5, hybrids = F)

compoplot(MenFas.clust.k2, n.col = 2, col.pal = hybridpal())
compoplot(MenFas.clust.k2.bc, n.col = 2, col.pal = hybridpal())

table(Buccinum.clust.k4$group) # 11
table(Buccinum.clust.k5$group) # 28


df.g4 <- data.frame(ID_GQ = names( Buccinum.clust.k4$group),
                    K = Buccinum.clust.k4$group)

df.g5 <- data.frame(ID_GQ = names( Buccinum.clust.k5$group),
                    K = Buccinum.clust.k5$group)

GOOD.SP <- df.g5 %>% dplyr::filter(K == 1,
                                   ID_GQ %in% names(res[res<=0.50])) %>% pull(ID_GQ)
length(GOOD.SP)



pca.test2  <- adegenet::glPca(gl.data[indNames(gl.data) %in% GOOD.SP,], center = TRUE, scale = FALSE,  
                             parallel = TRUE, nf = 10)

gPCA <- pca.test2 %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  #left_join(df.g5, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC2, y = score.PC3, col =Region_echantillonnage)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #facet_wrap(~Numero_unique_groupe) +
  #scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data)), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test2)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test2)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA


# Stacks - Testing parameters - De-novo - Buccinum ----------------------------------------

# parameters
m.para <- c(3,4,5) # 3 is the recommand value .... but can test 4 and 5
M.para <- c(1:8) # Maybe test the 0
n.delta <- c(-1, 0, 1)# n = M

get.value("test.denovo.path")
get.value("test.denovo.log")

file.exists(get.value("test.denovo.path"))
file.exists(get.value("test.denovo.log"))
file.exists(file.path(get.value("info.path"), "popmap.test_undatum.txt"))
#dir.create(get.value("test.denovo.path"))

# WILL RUN ONLY IF THE FILE DOESN'T EXIST!!!
# Check if it's paired or not ...

for(m in m.para){
  for(M in M.para){
    for(nn in n.delta){
      
      # last missing parameters
      n.para <- M + nn
      
      # Directory name
      test.dir <-  file.path(get.value("test.denovo.path"), paste0("stack.m",m,".M",M,".n",n.para))
      
      # Run the code only if the directory doesn't exist
      if(file.exists(test.dir) == F){
        
        cat("\nProcessing: ", test.dir, "\n", sep="")
        
        dir.create(test.dir)
        
        # Run denovo_map.pl
        
        cmd <- paste("--samples", get.value("demulti.path"),
                     "--popmap", file.path(get.value("info.path"), "popmap.test_undatum.txt"),
                     "-o", test.dir,
                     "-M", M,
                     "-n", n.para,
                     "-m", m, # default = 3
                     "--paired",
                     "-T", 8,
                     "-X", "\"populations:-r 0.80\""
        )
        
        A <- system2("denovo_map.pl", cmd, stdout=T, stderr=T)
        
        # save a log file
        
        cat(file = file.path(get.value("test.denovo.log"), paste0("stack.m",m,".M",M,".n",n.para,".log" )),
            cmd, "\n\n",
            A, # what to put in my file
            append= F, sep = "\n")
      } # Close the file exist condition
      
    } # close nn (aka n)
  } # Close M
} # Close m


# Extract statistics

nloci.res <- data.frame(m = character(), M = character(), n = character(),
                        r80 = character(), n_loci = numeric())

nsnps.res <- data.frame(m = character(), M = character(), n = character(),
                        r80 = character(), n_snps = numeric(), n_loci = numeric())



for(x in list.files(get.value("test.denovo.path"), full.names = T)){
  
  # Extract parameters
  info.param <- x %>% str_remove(get.value("test.denovo.path")) %>%
    str_remove("/stack.")
  
  m.param <- sapply(str_split(info.param, "\\."), `[`)[[1]] %>% str_remove("m")
  M.param <- sapply(str_split(info.param, "\\."), `[`)[[2]] %>% str_remove("M")
  n.param <- sapply(str_split(info.param, "\\."), `[`)[[3]] %>% str_remove("n")
  
  # N loci
  
  cmd <- paste(file.path(x, "populations.log.distribs"),
               "samples_per_loc_postfilters")
  
  res <- system2("stacks-dist-extract", cmd, stdout = T)
  
  data <- res[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  names(data) <- c("n_samples", "n_loci")
  
  nloci.res <- bind_rows(nloci.res,
                         data.frame(m = m.param, M = M.param, n = n.param,
                                    r80 = "yes", n_loci = as.numeric(as.character(data$n_loci)) %>%  sum() )
  )
  
  # N SNP / locus
  
  cmd <- paste(file.path(x, "populations.log.distribs"),
               "snps_per_loc_postfilters")
  
  res <- system2("stacks-dist-extract", cmd, stdout = T)
  
  data <- res[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  names(data) <- c("n_snps", "n_loci")
  
  data$n_snps <- as.numeric(as.character(data$n_snps))
  data$n_loci <- as.numeric(as.character(data$n_loci))
  
  data$m   <- m.param
  data$M   <- M.param
  data$n   <- n.param
  data$r80 <- "yes"
  
  nsnps.res <- bind_rows(nsnps.res, data)
  
}

nloci.res <- nloci.res %>% left_join(nsnps.res %>% filter(n_snps == 0) %>% select(m,n,M, no_div = n_loci),
                                     by = c("m", "M", "n")) %>%
  mutate(n_loci_poly = n_loci - no_div)

nloci.res
nsnps.res %>% head()

nsnps.res %>% group_by( m, M, n, r80) %>% summarise(n_loc_total = sum(n_loci),
                                                    n_loc_poly = sum(n_loci[n_snps > 0])) 

nloci.res %>% head()

# Le meilleur paramètre
nloci.res %>% select(m, M, n, n_loci_poly) %>% arrange(desc(n_loci_poly)) %>% head(10)

graph1.1 <- nloci.res %>% gather(n_loci, n_loci_poly, key = "loci", value = "n_loci") %>%
  mutate(n.rel = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1")),
         grouped = paste(loci, n.rel)) %>%
  ggplot(aes(x = M, y = n_loci, group = grouped, col = n.rel)) +
  geom_line(aes(lty = loci), show.legend = F) +
  geom_point() +
  labs(y = "No. of loci\nshared by 80% of samples",
       colour = "n") +
  facet_grid(.~m, labeller = label_both) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme_bw()

graph1.1

graph1.2 <- nsnps.res %>% mutate(n_snps_tot = n_snps * n_loci) %>%
  group_by(m, M, n) %>%
  summarise(n_snps = sum(n_snps_tot)) %>%
  mutate(n.rel = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1"))) %>%
  ggplot(aes(x = M, y = n_snps, group = n.rel, col = n.rel)) +
  geom_line() +
  geom_point() +
  labs(y = "No. of snps\nshared by 80% of samples", colour = "n") +
  facet_grid(.~m, labeller = label_both) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme_bw()

graph1.2

graph1.3 <- nsnps.res %>% mutate(n = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1"))) %>%
  filter(m == 3, n %in% c("M + 1", "M", "M - 1")) %>%
  group_by(m, M, n) %>%
  mutate(n_snps = ifelse(n_snps > 10, 11, n_snps),
         p_loci = n_loci / sum(n_loci)) %>%
  ggplot(aes(x = n_snps, y = p_loci, fill = M)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Percentage of loci") +
  facet_grid(n ~ m, labeller = label_both)+
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme_bw()

graph1.3

ggsave(filename = file.path(get.value("test.denovo.log"), "NLoci_test.png"),
       plot = graph1.1,
       width = 7.5, height = 4, units = "in")

ggsave(filename = file.path(get.value("test.denovo.log"), "NSnps_test.png"),
       plot = graph1.2,
       width = 7.5, height = 4, units = "in")

ggsave(filename = file.path(get.value("test.denovo.log"), "NSnpsByLocus_test.png"),
       plot = graph1.3,
       width = 7.5, height = 4, units = "in")


# Stacks - Testing parameters - De-novo - ALL ----------------------------------------

Nreads.data <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))
Nreads.data <- Nreads.data %>% left_join(pop.data, by = c("Filename" = "ID_GQ")) 


test.ALL <- Nreads.data %>% filter(Cat_sample == "Sample") %>% 
  group_by(Filename, Espece) %>% 
  summarise(Retained = sum(Retained)) %>% 
  filter(Retained >= quantile(Retained, probs = 0.25),
         Retained <= quantile(Retained, probs = 0.75)) %>% 
  group_by(Espece) %>% #summarise((N = n()))
sample_n(3) %>% 
pull(Filename)  

test.ALL

Nreads.data %>% dplyr::filter(Espece == "terraenovae")

#write.table(pop.info %>% select(Sample, POP) %>% 
#            filter(Sample %in% test.ALL) %>% 
#            mutate(POP = "noPOP"), 
#            file = file.path(get.value("info.path"), "popmap.test_ALL.txt"),
#            quote = FALSE, sep = "\t",
#            row.names = F, col.names = F)


# parameters
m.para <- c(3) # 3 is the recommand value .... but can test 4 and 5
M.para <- c(1:8) # Maybe test the 0
n.delta <- c(0)# n = M



file.exists("./00_Data/04a_Test.denovo.ALL")
file.exists(get.value("test.denovo.log"))
file.exists(file.path(get.value("info.path"), "popmap.test_ALL.txt"))
#dir.create(get.value("test.denovo.path"))

# WILL RUN ONLY IF THE FILE DOESN'T EXIST!!!
# Check if it's paired or not ...

for(m in m.para){
  for(M in M.para){
    for(nn in n.delta){
      
      # last missing parameters
      n.para <- M + nn
      
      # Directory name
      test.dir <-  file.path("./00_Data/04a_Test.denovo.ALL", paste0("stack.m",m,".M",M,".n",n.para))
      
      # Run the code only if the directory doesn't exist
      if(file.exists(test.dir) == F){
        
        cat("\nProcessing: ", test.dir, "\n", sep="")
        
        dir.create(test.dir)
        
        # Run denovo_map.pl
        
        cmd <- paste("--samples", get.value("demulti.path"),
                     "--popmap", file.path(get.value("info.path"), "popmap.test_ALL.txt"),
                     "-o", test.dir,
                     "-M", M,
                     "-n", n.para,
                     "-m", m, # default = 3
                     "--paired",
                     "-T", 8,
                     "-X", "\"populations:-r 0.80\""
        )
        
        A <- system2("denovo_map.pl", cmd, stdout=T, stderr=T)
        
        # save a log file
        
        cat(file = file.path(get.value("test.denovo.log"), paste0("stack.m",m,".M",M,".n",n.para,"_ALL.log" )),
            cmd, "\n\n",
            A, # what to put in my file
            append= F, sep = "\n")
      } # Close the file exist condition
      
    } # close nn (aka n)
  } # Close M
} # Close m


# Extract statistics

nloci.res <- data.frame(m = character(), M = character(), n = character(),
                        r80 = character(), n_loci = numeric())

nsnps.res <- data.frame(m = character(), M = character(), n = character(),
                        r80 = character(), n_snps = numeric(), n_loci = numeric())



for(x in list.files("./00_Data/04a_Test.denovo.ALL", full.names = T)){
  
  # Extract parameters
  info.param <- x %>% str_remove("./00_Data/04a_Test.denovo.ALL") %>%
    str_remove("/stack.")
  
  m.param <- sapply(str_split(info.param, "\\."), `[`)[[1]] %>% str_remove("m")
  M.param <- sapply(str_split(info.param, "\\."), `[`)[[2]] %>% str_remove("M")
  n.param <- sapply(str_split(info.param, "\\."), `[`)[[3]] %>% str_remove("n")
  
  # N loci
  
  cmd <- paste(file.path(x, "populations.log.distribs"),
               "samples_per_loc_postfilters")
  
  res <- system2("stacks-dist-extract", cmd, stdout = T)
  
  data <- res[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  names(data) <- c("n_samples", "n_loci")
  
  nloci.res <- bind_rows(nloci.res,
                         data.frame(m = m.param, M = M.param, n = n.param,
                                    r80 = "yes", n_loci = as.numeric(as.character(data$n_loci)) %>%  sum() )
  )
  
  # N SNP / locus
  
  cmd <- paste(file.path(x, "populations.log.distribs"),
               "snps_per_loc_postfilters")
  
  res <- system2("stacks-dist-extract", cmd, stdout = T)
  
  data <- res[c(-1, -2)] %>% str_split(pattern = "\t")
  
  data <-  data.frame(matrix(unlist(data), nrow= length(data), byrow = T))
  names(data) <- c("n_snps", "n_loci")
  
  data$n_snps <- as.numeric(as.character(data$n_snps))
  data$n_loci <- as.numeric(as.character(data$n_loci))
  
  data$m   <- m.param
  data$M   <- M.param
  data$n   <- n.param
  data$r80 <- "yes"
  
  nsnps.res <- bind_rows(nsnps.res, data)
  
}

nloci.res <- nloci.res %>% left_join(nsnps.res %>% filter(n_snps == 0) %>% select(m,n,M, no_div = n_loci),
                                     by = c("m", "M", "n")) %>%
  mutate(n_loci_poly = n_loci - no_div)

nloci.res
nsnps.res %>% head()

nsnps.res %>% group_by( m, M, n, r80) %>% summarise(n_loc_total = sum(n_loci),
                                                    n_loc_poly = sum(n_loci[n_snps > 0])) 

nloci.res %>% head()

# Le meilleur paramètre
nloci.res %>% select(m, M, n, n_loci_poly) %>% arrange(desc(n_loci_poly)) %>% head(10)

graph1.1 <- nloci.res %>% gather(n_loci, n_loci_poly, key = "loci", value = "n_loci") %>%
  mutate(n.rel = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1")),
         grouped = paste(loci, n.rel)) %>%
  ggplot(aes(x = M, y = n_loci, group = grouped, col = n.rel)) +
  geom_line(aes(lty = loci), show.legend = F) +
  geom_point() +
  labs(y = "No. of loci\nshared by 80% of samples",
       colour = "n") +
  facet_grid(.~m, labeller = label_both) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme_bw()

graph1.1

graph1.2 <- nsnps.res %>% mutate(n_snps_tot = n_snps * n_loci) %>%
  group_by(m, M, n) %>%
  summarise(n_snps = sum(n_snps_tot)) %>%
  mutate(n.rel = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1"))) %>%
  ggplot(aes(x = M, y = n_snps, group = n.rel, col = n.rel)) +
  geom_line() +
  geom_point() +
  labs(y = "No. of snps\nshared by 80% of samples", colour = "n") +
  facet_grid(.~m, labeller = label_both) +
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme_bw()

graph1.2

graph1.3 <- nsnps.res %>% mutate(n = ifelse(M==n, "M", ifelse(n>M, "M + 1", "M - 1"))) %>%
  filter(m == 3, n %in% c("M + 1", "M", "M - 1")) %>%
  group_by(m, M, n) %>%
  mutate(n_snps = ifelse(n_snps > 10, 11, n_snps),
         p_loci = n_loci / sum(n_loci)) %>%
  ggplot(aes(x = n_snps, y = p_loci, fill = M)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Percentage of loci") +
  facet_grid(n ~ m, labeller = label_both)+
  guides(colour = guide_legend(title.hjust = 0.5))+
  theme_bw()

graph1.3

ggsave(filename = file.path(get.value("test.denovo.log"), "NLoci_test.png"),
       plot = graph1.1,
       width = 7.5, height = 4, units = "in")

ggsave(filename = file.path(get.value("test.denovo.log"), "NSnps_test.png"),
       plot = graph1.2,
       width = 7.5, height = 4, units = "in")

ggsave(filename = file.path(get.value("test.denovo.log"), "NSnpsByLocus_test.png"),
       plot = graph1.3,
       width = 7.5, height = 4, units = "in")



# ALL - m3 M4 n4 ----------------------------------------------------------

# # Ustacks - Build loci de novo --------------------------------------------

files.to.use <- list.files(get.value("demulti.path"), full.names = T, pattern = ".1.fq.gz") %>% str_subset(".rem", negate = T)
length(files.to.use)


# Attention - doit rouler des ID  "-i" unique, ce qui peut être problématique
duplicated(files.to.use) %>% table()

# Version 8 parallel * 8 threads

mclapply(seq_along(files.to.use),
         FUN = function(x) {
           
           cat("Processing: ", files.to.use[x], "\n", sep="")
           
           cmd <- paste("-t", "gzfastq",
                        "-f", files.to.use[x],
                        "-o", "00_Data/05a_Stacks.denovo.ALL.M4n4/",
                        "-i", x, # Sequential number - THIS IS IMPORTANT
                        "-m", 3,
                        "--name", files.to.use[x] %>% str_remove(get.value("demulti.path")) %>%
                          str_remove(".1.fq.gz") %>%
                          str_remove("/"),
                        "-M", 4,
                        "-p", 8
           )
           
           A <- system2("ustacks", cmd, stdout=T, stderr=T)
           
           # save a log file
           
           if(file.exists(get.value("stacks.log")) == F){
             dir.create(get.value("stacks.log"))
           }
           
           cat(file = file.path(get.value("stacks.log"), paste0(files.to.use[x] %>% str_remove(get.value("demulti.path")) %>%
                                                                  str_remove(".fq.gz") %>%
                                                                  str_remove("/"),
                                                                ".ustacks.m4M4.denovo.log")),
               "\n", cmd, "\n",
               A, # what to put in my file
               append= F, sep = "\n")
           
         },
         mc.cores = 8
) #)



# Read coverage

Coverage.res <- data.frame(Sample = character(),
                           Coverage.mean = numeric(),
                           Coverage.SD = numeric(),
                           Nreads = numeric(),
                           Preads = numeric(),
                           stringsAsFactors = F)

for(x in list.files( get.value("stacks.log"), pattern = "ustacks.m4M4.denovo.log", full.names = T)){
  
  res <- readLines(x) %>% str_subset("Final coverage") %>%
    str_remove("Final coverage: ") %>%
    str_split("; ")
  
  
  if(length(res) == 0){
    
    res <- c(x %>% str_remove(get.value("stacks.log")) %>% str_remove(".ustacks.m4M4.denovo.log") %>% str_remove("[/]") %>% str_remove("[.]1") ,
             0,
             0,
             0,
             0
    )
  }
  if(length(res) == 1){
    res.reads <- res[[1]][4] %>% str_remove("n_reads=") %>%
      str_split("[(]")
    
    res <- c(x %>% str_remove(get.value("stacks.log")) %>% str_remove(".ustacks.m4M4.denovo.log") %>% str_remove("[/]") %>% str_remove("[.]1") ,
             res[[1]][1] %>% str_remove("mean="),
             res[[1]][2] %>% str_remove("stdev="),
             res.reads[[1]][1],
             res.reads[[1]][2] %>% str_remove("%[)]")
    )
  }
  Coverage.res[nrow(Coverage.res)+1,] <- res
  
}

Coverage.res %>% head(30)
pop.data %>% head()

Coverage.res <- Coverage.res %>% mutate(#Sample = str_remove(Sample, "\\.1"),
  Coverage.mean = as.numeric(as.character(Coverage.mean)),
  Coverage.SD = as.numeric(as.character(Coverage.SD)),
  Nreads = as.numeric(as.character(Nreads)),
  Preads = as.numeric(as.character(Preads))) %>%
  left_join(pop.data, by = c("Sample" = "ID_GQ"))

summary(Coverage.res)

Coverage.res %>% View()

# Basic stats overall

Coverage.res %>%
  summarise(Mean.all = mean(Coverage.mean),
            Mean.5x = mean(Coverage.mean[Coverage.mean >= 5]),
            Mean.10x = mean(Coverage.mean[Coverage.mean >= 10]),
            n.all =n(),
            n.5x = length(Coverage.mean[Coverage.mean >= 5]),
            n.10x = length(Coverage.mean[Coverage.mean >= 10]) #,
            #n.60x = length(Coverage.mean[Coverage.mean >= 60]),
            #n.83.4x = length(Coverage.mean[Coverage.mean >= 83.4]),
            #Mean.60x = mean(Coverage.mean[Coverage.mean >= 60]),
            #Mean.83.4x = mean(Coverage.mean[Coverage.mean >= 83.4])
  ) #%>% View()

Coverage.res %>% ggplot(aes(x = Coverage.mean, fill = Espece)) +
  geom_histogram() +
  geom_vline(xintercept = 5) +
  facet_wrap(~Espece, scale = "free_y")
scale_x_continuous(limits = c(0,50))


Coverage.res %>% ggplot(aes(x = Nreads, y = Coverage.mean, col = Espece)) +
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = 10)

# Check N by pop

Coverage.res %>% group_by(Notes_groupe) %>% summarise( n.all =n(),
                                                       n.5x = length(Coverage.mean[Coverage.mean >= 5]),
                                                       n.10x = length(Coverage.mean[Coverage.mean >= 10])
                                                       ) %>% 
                  mutate(P.10X = n.10x / n.all)                 



# Cstacks - Assemble catalog - ALL ---------------------------------------------

#Create a list of random samples for the catalog (we want around 40-200 representative individuals)

# Coverage between 20 and 40
# No duplicates
# Not from plaque 2 HI (cause this is a duplicate plate)
Coverage.res$Espece_revision %>% unique()

ID.cov <- c(# Most sample site
  Coverage.res %>% filter(Cat_sample == "Sample",
                          Espece != "undatum",
                          Coverage.mean >= 5) %>%  pull(Sample),
  # Very low sample size site
  Coverage.res %>% filter(Cat_sample == "Sample",
                          Espece == "undatum",
                          Coverage.mean >= 5,
                          Coverage.mean <= 35) %>%
    group_by(Notes_groupe ) %>% #summarize(N = n()) %>% View()
    sample_n(5) %>% pull(Sample) %>% unique()
)

write.table(pop.info %>% filter(Sample %in% ID.cov) %>%
              select(Sample, POP),
            file = file.path(get.value("info.path"), "popmap_catalog_M4n4_ALL.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Run cstacks

file.prefix <- list.files("00_Data/05a_Stacks.denovo.ALL.M4n4/", full.names = T, pattern = ".snps.tsv.gz") %>%
  str_subset(paste(paste0(ID.cov, ".snps.tsv.gz"), collapse = "|")) %>%
  str_remove(".snps.tsv.gz")

ID.cov

#if(file.exists(file.path("00_Data/05a_Stacks.denovo.ALL.Mn1/", "Catalog_ALL")) == F){
#  dir.create(file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/", "Catalog_ALL"))
#}

cmd <- paste(#"-P", get.value("ustacks.path"),
  "-o", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/"),
  #"-M", file.path(get.value("info.path"), "popmap_catalog_PbPm.txt"),
  "-s", file.prefix,
  "-n", 4,
  "-p", 16
)

A <- system2("cstacks", cmd, stdout=T, stderr=T)
A %>% tail()
# save a log file

cat(file = file.path(get.value("stacks.log"), "cstacks.ALL.M4n4.denovo.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

# Sstacks - Match to catalog - ALL ----------------------------------------------

#Create the list of individuals to use

Coverage.res %>% nrow()

ID.cov <- Coverage.res %>% filter(between(Coverage.mean,5 ,35)) %>%
  pull(Sample)

ID.cov %>% length()

file.prefix <- list.files("00_Data/05a_Stacks.denovo.ALL.M4n4", full.names = T, pattern = ".snps.tsv.gz") %>%
  str_subset(paste(paste0(ID.cov, ".snps.tsv.gz"), collapse = "|")) %>%
  str_remove(".snps.tsv.gz")

file.prefix

# J'ai essayé de le faire rouler avec mcapply 8 * 8 mais ça prend trop de mémoire et ça plante

write.table(pop.data %>% filter(ID_GQ %in% ID.cov) %>%
              select(ID_GQ) %>% mutate(POP = "nopop"),
            
            file = file.path(get.value("info.path"), "popmap_sstacks_ALL_M4n4.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Run sstacks

cmd <- paste(#"-P", get.value("ustacks.path"),
  "-o", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/"),
  "-c",  file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/"),
  paste("-s", file.prefix, collapse = " "),
  # "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
  "-p", 20
)


A <- system2("sstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "sstacks.denovo.ALL.M4n4.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Tsv2BAM - Transpose - ALL -----------------------------------------------------

# Run tsv2bam

# *Notice - when running with 96 cores, it failed for not apparent reason. So I'm running it
# with only 1 core. Maybe something in between may be OK...


# The Catalog should be in the P with the samples

cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/" ),
             #"-o", get.value("cstacks.path"),
             #paste("-s", file.prefix, collapse = " "),
             "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M4n4.txt"),
             "-t", 20,
             "-R", get.value("demulti.path")  # for paired-end read
)

A <- system2("tsv2bam", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "tsv2bam.denovo.ALL.M4n4.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Gstacks - incorporate PE - ALL ------------------------------------------------

# Was run on genyoda because it need more than 400 gb of RAM
# Then send back 4 files, gstacks* and catalog*


cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/"),
             
             "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M4n4.txt"),
             "-t", 8
)

cmd


A <- system2("gstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "gstacks.denovo.ALL.M4n4.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Info we can get from the use of gstacks - denovo

cmd <- paste(file.path("00_Data/05a_Stacks.denovo.ALL.M4n4", "gstacks.log.distribs"),
             "effective_coverages_per_sample")

res <- system2("stacks-dist-extract", cmd, stdout = T)

cov.data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

cov.data <-  data.frame(matrix(unlist(cov.data), nrow= length(cov.data), byrow = T), stringsAsFactors = F)
names(cov.data) <- res[3] %>% str_split(pattern = "\t") %>% unlist()

cov.data %>% head()

cov.data$mean_cov <- as.numeric(as.character(cov.data$mean_cov))
cov.data$mean_cov_ns <- as.numeric(as.character(cov.data$mean_cov_ns))
cov.data$n_used_fw_reads <- as.numeric(as.character(cov.data$n_used_fw_reads))
cov.data$n_loci <- as.numeric(as.character(cov.data$n_loci))

# Add info about raw reads
RawReads <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))
head(RawReads)

cov.data <- cov.data %>% left_join(RawReads %>% select(sample = Filename, Raw = Retained)) %>% 
  left_join(pop.data %>% select(sample = ID_GQ, Espece, Notes_groupe))
cov.data %>% head()

#write_csv(cov.data, file.path(get.value("stacks.log"),"AllIndividuals_DeNovo.ALL.M4n4_NreadsNloci.csv"))

cov.data <- read_csv(file.path(get.value("stacks.log"),"AllIndividuals_DeNovo.ALL.M4n4_NreadsNloci.csv"))

graph1.0 <- cov.data %>%
  ggplot(aes(x = Raw/2, y = n_used_fw_reads, col = Espece)) +
  geom_point() +
  #facet_wrap(~Espece) %>% 
  theme_bw()
graph1.0

graph1.1 <- cov.data %>%
  ggplot(aes(x = mean_cov_ns, y = n_loci, col = Espece)) +
  geom_point() +
  theme_bw()
graph1.1


graph1.2 <- cov.data %>%
  ggplot(aes(x = mean_cov_ns)) +
  geom_histogram()+
  theme_bw()
graph1.2

pdf(file.path(get.value("stacks.log"),"BasicGraph_PostGstacks_DeNovo.ALL.M4n4_2024-06-08.pdf"))
print(graph1.0)
print(graph1.1)
print(graph1.2)
dev.off()


# Test population - ALL --------------------------------------------

# Parameters
r.value       <- 0.75 # Minimum within pop
maf.value     <- 0.05 # Overall MAF
#maf.pop.value <- 0.05

write.table(pop.data %>% filter(ID_GQ %in% ID.cov) %>%
              mutate(Espece = ifelse(Espece == "undatum", "undatum", "other")) %>% 
              select(ID_GQ,  Espece), 
            
            file = file.path(get.value("info.path"), "popmap_ALL_M1n1_popSP.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/"),
             "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M4n4.txt/population_R75_NoPop"),
             "--out-path", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/"),
             "-t", 8,
             "-r", r.value, #              
             "--min-maf", maf.value,
             #"--smooth",
             "--write-single-snp",
             "--vcf"
)

A <- system2("populations", cmd, stdout=T, stderr=T)

cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),
             "-M", file.path(get.value("info.path"), "popmap_ALL_M1n1_popSP.txt"),
             "--out-path", file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/population_R75_SP"),
             "-p 2",
             "-t", 8,
             "-r", r.value, #              
             "--min-maf", maf.value,
             #"--smooth",
             "--write-single-snp",
             "--vcf"
)

A <- system2("populations", cmd, stdout=T, stderr=T)
A


# rapid PCA
file.path(current.wd,filter.ref.path,"06c_HW") %>% list.files()

#vcf.path <- file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/population_R75_NoPop/populations.snps.vcf")
#vcf.data <- vcfR::read.vcfR(vcf.path)

vcf.path <- file.path("00_Data/05a_Stacks.denovo.ALL.M4n4/population_R75_SP/populations.snps.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)


library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# subset

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

LOC.MAF10.NA10 <- filter.MAF.NA(gi.data, MAF.trs = 0.05, NA.trs = 0.10)


pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(pop.data) %>% pull(Espece)

table(pop(gl.data))

count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}


na.info <- data.frame(ID_GQ = indNames(gl.data),
                      NNA = count.ind.na.gl(gl.data))


hist(na.info$NNA)


na.info %>% left_join(pop.data) %>% 
  ggplot(aes(x = NNA, fill = Espece)) +
  geom_histogram() +
  facet_grid(~Espece)


no.loci.id <-  na.info %>% dplyr::filter(NNA ==1) %>% pull(ID_GQ)
no.loci.id


pca.sp  <- glPca(gl.data[indNames(gl.data) %nin% no.loci.id,], center = TRUE, scale = FALSE,  
                 parallel = TRUE, n.core =8, nf = 1000)

#save(list = c("pca.sp", "gl.data"),
#     file = file.path("00_Data/05a_Stacks.denovo.ALL.M4n4", "PCA.Rdata"))

load(file.path("00_Data/05a_Stacks.denovo.ALL.M4n4", "PCA.Rdata"))

gPCA.sp <- pca.sp %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter(NNA < 0.75) %>% 
  #dplyr::filter(Espece == "undatum") %>% 
  mutate(Espece =  paste("B.", Espece)) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, fill = Espece, shape = Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +

  #facet_wrap(~Espece) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  scale_shape_manual(values = c(22:24,21), name = "Species")+
  scale_fill_discrete(name = "Species") +
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  labs(#title = paste("All snps:",  nLoc(gl.final)),
     x = paste0("PC1 (", QuickPop::pca_var(pca.sp)$p.eig[1] %>% round(3) *100, "%)"),
     y = paste0("PC2 (", QuickPop::pca_var(pca.sp)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw() + theme()
gPCA.sp

ggsave(plot = gPCA.sp, file = "02_Results/00_Stacks/05a_Stacks.denovo/PCA_species_M4N4.png")



gPCA.sp.full.b <- pca.sp %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  #dplyr::filter(NNA < 0.75) %>% 
  #dplyr::filter(Espece == "undatum") %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, fill =NNA * 100, shape = paste("B.", Espece))) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_fill_distiller(palette = "Spectral", name = "% NA") +
  #facet_wrap(~Espece) +
  #  stat_ellipse(aes(col = Espece))+
  scale_shape_manual(values = c(22:24,21), name = "Species")+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.sp)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.sp)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw() 
gPCA.sp.full.b

ggsave(plot = gPCA.sp.full.b, file = "02_Results/00_Stacks/05a_Stacks.denovo/PCA_species_M4N4_NAall.png")

gPCA.sp.big <- ggpubr::ggarrange(gPCA.sp + guides(fill=guide_legend(nrow=2, byrow=TRUE),
                                                  shape=guide_legend(nrow=2, byrow=TRUE)) + theme(legend.position = "bottom"), 
                                 gPCA.sp.full.b + guides(shape = FALSE) + theme(legend.position = "bottom"),
                                 nrow = 1, labels = LETTERS, align = "hv")

gPCA.sp.big

ggsave(plot = gPCA.sp.big, file = "02_Results/00_Stacks/05a_Stacks.denovo/PCA_species_M4N4_BIG.png",
       width = 7.5, height = 5, bg = "white")

ID.to.sequence <- pca.sp %>% QuickPop::pca_scoretable(naxe = 6) %>% 
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter( score.PC1 >2, score.PC2 >1) %>% 
  pull(ID)

pca.sp %>% QuickPop::pca_scoretable(naxe = 6) %>% 
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter(ID %nin% ID.to.sequence,
                Espece!= "undatum") %>% 
  View()


# ALL - m3 M1 n1 ----------------------------------------------------------

# # Ustacks - Build loci de novo --------------------------------------------

# Analyze only the 1.fq.gz
files.to.use <- list.files(get.value("demulti.path"), full.names = T, pattern = ".1.fq.gz") %>% str_subset(".rem", negate = T)
length(files.to.use)


# Attention - doit rouler des ID  "-i" unique, ce qui peut être problématique
duplicated(files.to.use) %>% table()

# Version 8 parallel * 8 threads

mclapply(seq_along(files.to.use),
         FUN = function(x) {

           cat("Processing: ", files.to.use[x], "\n", sep="")

           cmd <- paste("-t", "gzfastq",
                        "-f", files.to.use[x],
                        "-o", "00_Data/05a_Stacks.denovo.ALL.M1n1/",
                        "-i", x, # Sequential number - THIS IS IMPORTANT
                        "-m", 3,
                        "--name", files.to.use[x] %>% str_remove(get.value("demulti.path")) %>%
                          str_remove(".1.fq.gz") %>%
                          str_remove("/"),
                        "-M", 1,
                        "-p", 8
           )

           A <- system2("ustacks", cmd, stdout=T, stderr=T)

           # save a log file

           if(file.exists(get.value("stacks.log")) == F){
             dir.create(get.value("stacks.log"))
           }

           cat(file = file.path(get.value("stacks.log"), paste0(files.to.use[x] %>% str_remove(get.value("demulti.path")) %>%
                                                                  str_remove(".fq.gz") %>%
                                                                  str_remove("/"),
                                                                ".ustacks.m3M1.denovo.log")),
               "\n", cmd, "\n",
               A, # what to put in my file
               append= F, sep = "\n")

         },
         mc.cores = 8
) #)


# Read coverage

Coverage.res <- data.frame(Sample = character(),
                           Coverage.mean = numeric(),
                           Coverage.SD = numeric(),
                           Nreads = numeric(),
                           Preads = numeric(),
                           stringsAsFactors = F)

for(x in list.files( get.value("stacks.log"), pattern = "ustacks.m3M1.denovo.log", full.names = T)){

  res <- readLines(x) %>% str_subset("Final coverage") %>%
    str_remove("Final coverage: ") %>%
    str_split("; ")

  
  if(length(res) == 0){
    
      res <- c(x %>% str_remove(get.value("stacks.log")) %>% str_remove(".ustacks.m3M1.denovo.log") %>% str_remove("[/]") %>% str_remove("[.]1") ,
           0,
           0,
           0,
           0
  )
}
  if(length(res) == 1){
  res.reads <- res[[1]][4] %>% str_remove("n_reads=") %>%
    str_split("[(]")

  res <- c(x %>% str_remove(get.value("stacks.log")) %>% str_remove(".ustacks.m3M1.denovo.log") %>% str_remove("[/]") %>% str_remove("[.]1") ,
           res[[1]][1] %>% str_remove("mean="),
           res[[1]][2] %>% str_remove("stdev="),
           res.reads[[1]][1],
           res.reads[[1]][2] %>% str_remove("%[)]")
  )
  }
  Coverage.res[nrow(Coverage.res)+1,] <- res

}

Coverage.res %>% head(30)
pop.data %>% head()

Coverage.res <- Coverage.res %>% mutate(#Sample = str_remove(Sample, "\\.1"),
                                        Coverage.mean = as.numeric(as.character(Coverage.mean)),
                                        Coverage.SD = as.numeric(as.character(Coverage.SD)),
                                        Nreads = as.numeric(as.character(Nreads)),
                                        Preads = as.numeric(as.character(Preads))) %>%
  left_join(pop.data, by = c("Sample" = "ID_GQ"))

summary(Coverage.res)

Coverage.res %>% View()

# Basic stats overall

Coverage.res %>%
  summarise(Mean.all = mean(Coverage.mean),
            Mean.5x = mean(Coverage.mean[Coverage.mean >= 5]),
            Mean.10x = mean(Coverage.mean[Coverage.mean >= 10]),
            n.all =n(),
            n.5x = length(Coverage.mean[Coverage.mean >= 5]),
            n.10x = length(Coverage.mean[Coverage.mean >= 10]) #,
            #n.60x = length(Coverage.mean[Coverage.mean >= 60]),
            #n.83.4x = length(Coverage.mean[Coverage.mean >= 83.4]),
            #Mean.60x = mean(Coverage.mean[Coverage.mean >= 60]),
            #Mean.83.4x = mean(Coverage.mean[Coverage.mean >= 83.4])
  ) #%>% View()


Coverage.res %>% ggplot(aes(x = Coverage.mean, fill = Espece)) +
  geom_histogram() +
  geom_vline(xintercept = 5) +
  facet_wrap(~Espece, scale = "free_y")
scale_x_continuous(limits = c(0,50))


Coverage.res %>% ggplot(aes(x = Nreads, y = Coverage.mean, col = Espece)) +
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = 10)

# Check N by pop

Coverage.res %>% group_by(Notes_groupe) %>% summarise( n.all =n(),
                                              n.5x = length(Coverage.mean[Coverage.mean >= 5]),
                                              n.10x = length(Coverage.mean[Coverage.mean >= 10])) #%>% write_csv("test_25oct.csv")



# Cstacks - Assemble catalog - ALL ---------------------------------------------

#Create a list of random samples for the catalog (we want around 40-200 representative individuals)

# Coverage between 20 and 40
# No duplicates
# Not from plaque 2 HI (cause this is a duplicate plate)
Coverage.res$Espece_revision %>% unique()

ID.cov <- c(# Most sample site
  Coverage.res %>% filter(Cat_sample == "Sample",
                          Espece != "undatum",
                          Coverage.mean >= 5) %>%  pull(Sample),
  # Very low sample size site
  Coverage.res %>% filter(Cat_sample == "Sample",
                          Espece == "undatum",
                          Coverage.mean >= 5,
                          Coverage.mean <= 35) %>%
    group_by(Notes_groupe ) %>% #summarize(N = n()) %>% View()
    sample_n(5) %>% pull(Sample) %>% unique()
)

write.table(pop.info %>% filter(Sample %in% ID.cov) %>%
              select(Sample, POP),
            file = file.path(get.value("info.path"), "popmap_catalog_M1n1_ALL.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Run cstacks

file.prefix <- list.files("00_Data/05a_Stacks.denovo.ALL.M1n1/", full.names = T, pattern = ".snps.tsv.gz") %>%
   str_subset(paste(paste0(ID.cov, ".snps.tsv.gz"), collapse = "|")) %>%
  str_remove(".snps.tsv.gz")

ID.cov

if(file.exists(file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/", "Catalog_ALL")) == F){
  dir.create(file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/", "Catalog_ALL"))
}

cmd <- paste(#"-P", get.value("ustacks.path"),
             "-o", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/", "Catalog_ALL"),
             #"-M", file.path(get.value("info.path"), "popmap_catalog_PbPm.txt"),
             "-s", file.prefix,
             "-n", 1,
             "-p", 16
)

A <- system2("cstacks", cmd, stdout=T, stderr=T)

# save a log file

cat(file = file.path(get.value("stacks.log"), "cstacks.ALL.M1n1.denovo.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

# Sstacks - Match to catalog - ALL ----------------------------------------------

#Create the list of individuals to use

Coverage.res %>% nrow()

ID.cov <- Coverage.res %>% filter(between(Coverage.mean,5 ,35)) %>%
  pull(Sample)

ID.cov %>% length()

file.prefix <- list.files("00_Data/05a_Stacks.denovo.ALL.M1n1", full.names = T, pattern = ".snps.tsv.gz") %>%
     str_subset(paste(paste0(ID.cov, ".snps.tsv.gz"), collapse = "|")) %>%
     str_remove(".snps.tsv.gz")

file.prefix

# J'ai essayé de le faire rouler avec mcapply 8 * 8 mais ça prend trop de mémoire et ça plante

write.table(pop.data %>% filter(ID_GQ %in% ID.cov) %>%
            select(ID_GQ) %>% mutate(POP = "nopop"),
            
            file = file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Run sstacks

cmd <- paste(#"-P", get.value("ustacks.path"),
             "-o", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),
             "-c",  file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/", "Catalog_ALL"),
             paste("-s", file.prefix, collapse = " "),
            # "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
             "-p", 20
)


A <- system2("sstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "sstacks.denovo.ALL.M1n1.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

#list.files(pattern = "matches.tsv.gz", full.names = T)
 
# old.name <- list.files(file.path("00_Data/05a_Stacks.denovo.ALL.M1n1", "Catalog_ALL"), pattern = "matches.tsv.gz", full.names = T)
# old.name
# 
# new.name <-  old.name %>% str_remove("Catalog_ALL/")
# new.name
# 
# file.rename(from = old.name,
#             to = new.name)

# Tsv2BAM - Transpose - ALL -----------------------------------------------------

# Run tsv2bam

# *Notice - when running with 96 cores, it failed for not apparent reason. So I'm running it
# with only 1 core. Maybe something in between may be OK...


# The Catalog should be in the P with the samples

cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/" ),
             #"-o", get.value("cstacks.path"),
             #paste("-s", file.prefix, collapse = " "),
              "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
             "-t", 20,
             "-R", get.value("demulti.path")  # for paired-end read
)

A <- system2("tsv2bam", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "tsv2bam.denovo.ALL.M1n1.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Gstacks - incorporate PE - ALL ------------------------------------------------

# Was run on genyoda because it need more than 400 gb of RAM
# Then send back 4 files, gstacks* and catalog*


cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),

             "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
             "-t", 8
)

cmd


A <- system2("gstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "gstacks.denovo.ALL.M1n1.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Info we can get from the use of gstacks - denovo

cmd <- paste(file.path("00_Data/05a_Stacks.denovo.ALL.M1n1", "gstacks.log.distribs"),
             "effective_coverages_per_sample")

res <- system2("stacks-dist-extract", cmd, stdout = T)

cov.data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

cov.data <-  data.frame(matrix(unlist(cov.data), nrow= length(cov.data), byrow = T), stringsAsFactors = F)
names(cov.data) <- res[3] %>% str_split(pattern = "\t") %>% unlist()

cov.data %>% head()

cov.data$mean_cov <- as.numeric(as.character(cov.data$mean_cov))
cov.data$mean_cov_ns <- as.numeric(as.character(cov.data$mean_cov_ns))
cov.data$n_used_fw_reads <- as.numeric(as.character(cov.data$n_used_fw_reads))
cov.data$n_loci <- as.numeric(as.character(cov.data$n_loci))

# Add info about raw reads
RawReads <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))
head(RawReads)

cov.data <- cov.data %>% left_join(RawReads %>% select(sample = Filename, Raw = Retained)) %>% 
                         left_join(pop.data %>% select(sample = ID_GQ, Espece, Notes_groupe))
cov.data %>% head()

#write_csv(cov.data, file.path(get.value("stacks.log"),"AllIndividuals_DeNovo.ALL.M1n1_NreadsNloci.csv"))

cov.data <- read_csv(file.path(get.value("stacks.log"),"AllIndividuals_DeNovo.ALL.M1n1_NreadsNloci.csv"))

graph1.0 <- cov.data %>%
  ggplot(aes(x = Raw/2, y = n_used_fw_reads, col = Espece)) +
  geom_point() +
  #facet_wrap(~Espece) %>% 
  theme_bw()
graph1.0

graph1.1 <- cov.data %>%
  ggplot(aes(x = mean_cov_ns, y = n_loci, col = Espece)) +
  geom_point() +
  theme_bw()
graph1.1


graph1.2 <- cov.data %>%
  ggplot(aes(x = mean_cov_ns)) +
  geom_histogram()+
  theme_bw()
graph1.2

pdf(file.path(get.value("stacks.log"),"BasicGraph_PostGstacks_DeNovo.ALL.M1n1_2024-06-04.pdf"))
print(graph1.0)
print(graph1.1)
print(graph1.2)
dev.off()

# Test population - ALL --------------------------------------------

# Parameters
r.value       <- 0.75 # Minimum within pop
maf.value     <- 0.05 # Overall MAF
#maf.pop.value <- 0.05

# Filtering #1 : -r and -MAF ---------------------------------------------

write.table(pop.data %>% filter(ID_GQ %in% ID.cov) %>%
            select(ID_GQ,  Espece),
            
            file = file.path(get.value("info.path"), "popmap_ALL_M1n1_popSP.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)


cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),
             "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
             "--out-path", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),
             "-t", 8,
             "-r", r.value, #              
             "--min-maf", maf.value,
             #"--smooth",
             "--write-single-snp",
             "--vcf"
)

A <- system2("populations", cmd, stdout=T, stderr=T)

cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),
             "-M", file.path(get.value("info.path"), "popmap_ALL_M1n1_popSP.txt"),
             "--out-path", file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/"),
             "-p 4",
             "-t", 8,
             "-r", r.value, #              
             "--min-maf", maf.value,
             #"--smooth",
             "--write-single-snp",
             "--vcf"
)



A <- system2("populations", cmd, stdout=T, stderr=T)
A


# rapid PCA
file.path(current.wd,filter.ref.path,"06c_HW") %>% list.files()

vcf.path <- file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/populations.snps.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)

vcf.path <- file.path("00_Data/05a_Stacks.denovo.ALL.M1n1/Population_r75_overall/populations.snps.vcf")
vcf.data <- vcfR::read.vcfR(vcf.path)


library(adegenet)
gl.data  <- vcfR::vcfR2genlight(vcf.data) 
gi.data  <- vcfR::vcfR2genind(vcf.data) 

# subset

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

LOC.MAF10.NA10 <- filter.MAF.NA(gi.data, MAF.trs = 0.10, NA.trs = 0.10)


pop(gl.data) <- data.frame(ID_GQ = indNames(gl.data)) %>% 
  left_join(pop.data) %>% pull(Espece)

table(pop(gl.data))

count.ind.na.gl <- function(gl){
  res <- apply(tab(gl,  NA.method = c("asis")), MARGIN = 1, FUN = function(l){   n.na <- length(l[is.na(l) == T])
  freq.na <- n.na / length(l)
  return(freq.na)
  })
  return(res)
  
}


na.info <- data.frame(ID_GQ = indNames(gl.data),
                      NNA = count.ind.na.gl(gl.data))


hist(na.info$NNA)


na.info %>% left_join(pop.data) %>% 
  ggplot(aes(x = NNA, fill = Espece)) +
  geom_histogram() +
  facet_grid(~Espece)


no.loci.id <-  na.info %>% dplyr::filter(NNA ==1) %>% pull(ID_GQ)
no.loci.id


pca.sp  <- glPca(gl.data[indNames(gl.data) %nin% no.loci.id,], center = TRUE, scale = FALSE,  
                   parallel = TRUE, n.core =8, nf = 1000)

save(list = c("pca.sp", "gl.data"),
file = file.path("00_Data/05a_Stacks.denovo.ALL.M1n1", "PCA.Rdata"))


gPCA.sp <- pca.sp %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  #dplyr::filter(NNA < 0.30) %>% 
  #dplyr::filter(Espece == "undatum") %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col =NNA, pch = Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_color_distiller(palette = "Spectral") +
  #facet_wrap(~Espece) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.sp)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.sp)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA.sp



gPCA <- pca.test %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  #dplyr::filter(NNA < 0.30) %>% 
  #dplyr::filter(Espece == "undatum") %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col =NNA, pch = Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_color_distiller(palette = "Spectral") +
  #facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA

gPCA.red <- pca.test.red %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  #dplyr::filter(NNA < 0.30) %>% 
  #dplyr::filter(Espece == "undatum") %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = NNA, pch = Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  scale_color_distiller(palette = "Spectral") +
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test.red)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test.red)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA.red



Bad.SP1 <-  pca.test  %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter(score.PC2 <= -5,
                Espece ==  "undatum") %>% pull(ID) 


Bad.SP2 <-  pca.test.red  %>% QuickPop::pca_scoretable(naxe = 6) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  dplyr::filter(score.PC1 >= 1,
                score.PC2 <= -1,
                Espece ==  "undatum") %>% pull(ID) 

table(Bad.SP1 == Bad.SP2) 
Bad.NA <- na.info %>% left_join(pop.data) %>%
               dplyr::filter(NNA >0.3,
                Espece ==  "undatum") %>% pull(ID_GQ) 

na.info %>% left_join(pop.data) %>%
  ggplot(aes(x = Notes_groupe, y = NNA))+
   geom_boxplot() + geom_jitter(height = 0, aes(col = Espece)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


na.info %>% left_join(Coverage.res, by = c("ID_GQ" = "Sample")) %>% 
  ggplot(aes(x= Coverage.mean, y = NNA, col = Espece)) +
  geom_point()

  
ggsave(filename = file.path("00_Data/05a_Stacks.denovo.ALL.M1n1","PCA_test.png"),
       plot = gPCA,
       width = 5, height = 4, units = "in")

gPCA.NA <- pca.test %>% QuickPop::pca_scoretable(naxe = 3) %>%
  left_join(pop.data, by = c("ID" = "ID_GQ")) %>% 
  left_join(na.info, by = c("ID" = "ID_GQ")) %>% 
  ggplot(aes(x = score.PC1, y = score.PC2, col = NNA, pch = Espece)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #facet_wrap(~Region_echantillonnage) +
  #  stat_ellipse(aes(col = Espece))+
  geom_point(alpha = 0.5, size = 3) +  
  scale_color_distiller(palette = "Spectral") +
  #  scale_colour_manual(name = "Region", values = c("black","blue", "darkorange","red", "magenta"))+    
  # annotate("text",  x=-Inf, y = Inf, label = paste("Test snps:",  nLoc(gl.data[, locNames(gl.data) %in% LOC.MAF10.NA05])), vjust=1, hjust=0) +
  
  labs(#title = paste("All snps:",  nLoc(gl.final)),
    x = paste0("PC1 (", QuickPop::pca_var(pca.test)$p.eig[1] %>% round(3) *100, "%)"),
    y = paste0("PC2 (", QuickPop::pca_var(pca.test)$p.eig[2] %>% round(3) *100, "%)")) +
  theme_bw()
gPCA.NA


# Buccinum - m3 M1 n1 ----------------------------------------------------------

# # Ustacks - Build loci de novo --------------------------------------------

# Analyze only the 1.fq.gz
files.to.use <- list.files(get.value("demulti.path"), full.names = T, pattern = ".1.fq.gz") %>% str_subset(".rem", negate = T)
length(files.to.use)


# Attention - doit rouler des ID  "-i" unique, ce qui peut être problématique
duplicated(files.to.use) %>% table()

# Version 8 parallel * 8 threads

mclapply(seq_along(files.to.use),
         FUN = function(x) {
           
           cat("Processing: ", files.to.use[x], "\n", sep="")
           
           cmd <- paste("-t", "gzfastq",
                        "-f", files.to.use[x],
                        "-o", "00_Data/05a_Stacks.denovo.buccinum.M1n1/",
                        "-i", x, # Sequential number - THIS IS IMPORTANT
                        "-m", 3,
                        "--name", files.to.use[x] %>% str_remove(get.value("demulti.path")) %>%
                          str_remove(".1.fq.gz") %>%
                          str_remove("/"),
                        "-M", 1,
                        "-p", 8
           )
           
           A <- system2("ustacks", cmd, stdout=T, stderr=T)
           
           # save a log file
           
           if(file.exists(get.value("stacks.log")) == F){
             dir.create(get.value("stacks.log"))
           }
           
           cat(file = file.path(get.value("stacks.log"), paste0(files.to.use[x] %>% str_remove(get.value("demulti.path")) %>%
                                                                  str_remove(".fq.gz") %>%
                                                                  str_remove("/"),
                                                                ".ustacks.m3M1.denovo.log")),
               "\n", cmd, "\n",
               A, # what to put in my file
               append= F, sep = "\n")
           
         },
         mc.cores = 8
) #)


# Read coverage

Coverage.res <- data.frame(Sample = character(),
                           Coverage.mean = numeric(),
                           Coverage.SD = numeric(),
                           Nreads = numeric(),
                           Preads = numeric(),
                           stringsAsFactors = F)

for(x in list.files( get.value("stacks.log"), pattern = "ustacks.m3M1.denovo.log", full.names = T)){
  
  res <- readLines(x) %>% str_subset("Final coverage") %>%
    str_remove("Final coverage: ") %>%
    str_split("; ")
  
  
  if(length(res) == 0){
    
    res <- c(x %>% str_remove(get.value("stacks.log")) %>% str_remove(".ustacks.m3M1.denovo.log") %>% str_remove("[/]") %>% str_remove("[.]1") ,
             0,
             0,
             0,
             0
    )
  }
  if(length(res) == 1){
    res.reads <- res[[1]][4] %>% str_remove("n_reads=") %>%
      str_split("[(]")
    
    res <- c(x %>% str_remove(get.value("stacks.log")) %>% str_remove(".ustacks.m3M1.denovo.log") %>% str_remove("[/]") %>% str_remove("[.]1") ,
             res[[1]][1] %>% str_remove("mean="),
             res[[1]][2] %>% str_remove("stdev="),
             res.reads[[1]][1],
             res.reads[[1]][2] %>% str_remove("%[)]")
    )
  }
  Coverage.res[nrow(Coverage.res)+1,] <- res
  
}

Coverage.res %>% head(30)
pop.data %>% head()

Coverage.res <- Coverage.res %>% mutate(#Sample = str_remove(Sample, "\\.1"),
  Coverage.mean = as.numeric(as.character(Coverage.mean)),
  Coverage.SD = as.numeric(as.character(Coverage.SD)),
  Nreads = as.numeric(as.character(Nreads)),
  Preads = as.numeric(as.character(Preads))) %>%
  left_join(pop.data, by = c("Sample" = "ID_GQ"))

summary(Coverage.res)

Coverage.res %>% View()

# Basic stats overall

Coverage.res %>%
  summarise(Mean.all = mean(Coverage.mean),
            Mean.5x = mean(Coverage.mean[Coverage.mean >= 5]),
            Mean.10x = mean(Coverage.mean[Coverage.mean >= 10]),
            n.all =n(),
            n.5x = length(Coverage.mean[Coverage.mean >= 5]),
            n.10x = length(Coverage.mean[Coverage.mean >= 10]) #,
            #n.60x = length(Coverage.mean[Coverage.mean >= 60]),
            #n.83.4x = length(Coverage.mean[Coverage.mean >= 83.4]),
            #Mean.60x = mean(Coverage.mean[Coverage.mean >= 60]),
            #Mean.83.4x = mean(Coverage.mean[Coverage.mean >= 83.4])
  ) #%>% View()


Coverage.res %>% ggplot(aes(x = Coverage.mean, fill = Espece)) +
  geom_histogram() +
  geom_vline(xintercept = 5) +
  facet_wrap(~Espece, scale = "free_y")
scale_x_continuous(limits = c(0,50))


Coverage.res %>% ggplot(aes(x = Nreads, y = Coverage.mean, col = Espece)) +
  geom_point(alpha = 0.5)+
  geom_hline(yintercept = 10)

# Check N by pop

Coverage.res %>% group_by(Notes_groupe) %>% summarise( n.all =n(),
                                                       n.5x = length(Coverage.mean[Coverage.mean >= 5]),
                                                       n.10x = length(Coverage.mean[Coverage.mean >= 10])) #%>% write_csv("test_25oct.csv")



# Cstacks - Assemble catalog - ALL ---------------------------------------------

#Create a list of random samples for the catalog (we want around 40-200 representative individuals)

# Coverage between 20 and 40
# No duplicates
# Not from plaque 2 HI (cause this is a duplicate plate)
Coverage.res$Espece_revision %>% unique()

ID.cov <-   Coverage.res %>% filter(Cat_sample == "Sample",
                          Espece == "undatum",
                          # The none undatum!!
                          Sample %nin% c("S_23_01903", "S_23_01910", "S_23_01913", "S_23_01922"),
                          Coverage.mean >= 10,
                          Coverage.mean <= 30) %>%
    group_by(Notes_groupe ) %>% #summarize(N = n()) %>% View()
    sample_n(10) %>% pull(Sample) %>% unique()


write.table(pop.info %>% filter(Sample %in% ID.cov) %>%
              select(Sample, POP),
            file = file.path(get.value("info.path"), "popmap_catalog_M1n1_Buccinum.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Run cstacks

file.prefix <- list.files("00_Data/05a_Stacks.denovo.buccinum.M1n1/", full.names = T, pattern = ".snps.tsv.gz") %>%
  str_subset(paste(paste0(ID.cov, ".snps.tsv.gz"), collapse = "|")) %>%
  str_remove(".snps.tsv.gz")

length(ID.cov) == length(file.prefix)

cmd <- paste(#"-P", get.value("ustacks.path"),
  "-o", file.path("00_Data/05a_Stacks.denovo.buccinum.M1n1"),
  #"-M", file.path(get.value("info.path"), "popmap_catalog_PbPm.txt"),
  "-s", file.prefix,
  "-n", 1,
  "-p", 16
)

A <- system2("cstacks", cmd, stdout=T, stderr=T)

# save a log file

cat(file = file.path(get.value("stacks.log"), "cstacks.ALL.M1n1.denovo.buccinum.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

# Sstacks - Match to catalog - ALL ----------------------------------------------

#Create the list of individuals to use

Coverage.res %>% nrow()

ID.cov <- Coverage.res %>% filter(between(Coverage.mean, 5 ,35),
                                          Espece == "undatum",
                                          # The none undatum!!
                                          Sample %nin% c("S_23_01903", "S_23_01910", "S_23_01913", "S_23_01922")) %>% 
  pull(Sample)

ID.cov %>% length()

file.prefix <- list.files("00_Data/05a_Stacks.denovo.buccinum.M1n1", full.names = T, pattern = ".snps.tsv.gz") %>%
  str_subset(paste(paste0(ID.cov, ".snps.tsv.gz"), collapse = "|")) %>%
  str_remove(".snps.tsv.gz")

file.prefix

# J'ai essayé de le faire rouler avec mcapply 8 * 8 mais ça prend trop de mémoire et ça plante

write.table(pop.data %>% filter(ID_GQ %in% ID.cov) %>%
              select(ID_GQ) %>% mutate(POP = "nopop"),
            
            file = file.path(get.value("info.path"), "popmap_sstacks_buccinum_M1n1.txt"),
            quote = FALSE, sep = "\t",
            row.names = F, col.names = F)

# Run sstacks

cmd <- paste(#"-P", get.value("ustacks.path"),
  "-o", file.path("00_Data/05a_Stacks.denovo.buccinum.M1n1/"),
  "-c",  file.path("00_Data/05a_Stacks.denovo.buccinum.M1n1/"),
  paste("-s", file.prefix, collapse = " "),
  # "-M", file.path(get.value("info.path"), "popmap_sstacks_ALL_M1n1.txt"),
  "-p", 20
)


A <- system2("sstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "sstacks.denovo.buccinum.M1n1.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")

#list.files(pattern = "matches.tsv.gz", full.names = T)

# old.name <- list.files(file.path("00_Data/05a_Stacks.denovo.ALL.M1n1", "Catalog_ALL"), pattern = "matches.tsv.gz", full.names = T)
# old.name
# 
# new.name <-  old.name %>% str_remove("Catalog_ALL/")
# new.name
# 
# file.rename(from = old.name,
#             to = new.name)

# Tsv2BAM - Transpose - ALL -----------------------------------------------------

# Run tsv2bam

# *Notice - when running with 96 cores, it failed for not apparent reason. So I'm running it
# with only 1 core. Maybe something in between may be OK...


# The Catalog should be in the P with the samples

cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.buccinum.M1n1/" ),
             #"-o", get.value("cstacks.path"),
             #paste("-s", file.prefix, collapse = " "),
             "-M", file.path(get.value("info.path"), "popmap_sstacks_buccinum_M1n1.txt"),
             "-t", 20,
             "-R", get.value("demulti.path")  # for paired-end read
)

A <- system2("tsv2bam", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "tsv2bam.denovo.buccinum.M1n1.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Gstacks - incorporate PE - ALL ------------------------------------------------

# Was run on genyoda because it need more than 400 gb of RAM
# Then send back 4 files, gstacks* and catalog*


cmd <- paste("-P", file.path("00_Data/05a_Stacks.denovo.buccinum.M1n1/"),
             
             "-M", file.path(get.value("info.path"), "popmap_sstacks_buccinum_M1n1.txt"),
             "-t", 8
)

cmd


A <- system2("gstacks", cmd, stdout=T, stderr=T)

cat(file = file.path(get.value("stacks.log"), "gstacks.denovo.buccinum.M1n1.log"),
    "\n", cmd, "\n",
    A, # what to put in my file
    append= F, sep = "\n")


# Info we can get from the use of gstacks - denovo

cmd <- paste(file.path("00_Data/05a_Stacks.denovo.buccinum.M1n1", "gstacks.log.distribs"),
             "effective_coverages_per_sample")

res <- system2("stacks-dist-extract", cmd, stdout = T)

cov.data <- res[c(-1, -2, -3)] %>% str_split(pattern = "\t")

cov.data <-  data.frame(matrix(unlist(cov.data), nrow= length(cov.data), byrow = T), stringsAsFactors = F)
names(cov.data) <- res[3] %>% str_split(pattern = "\t") %>% unlist()

cov.data %>% head()

cov.data$mean_cov <- as.numeric(as.character(cov.data$mean_cov))
cov.data$mean_cov_ns <- as.numeric(as.character(cov.data$mean_cov_ns))
cov.data$n_used_fw_reads <- as.numeric(as.character(cov.data$n_used_fw_reads))
cov.data$n_loci <- as.numeric(as.character(cov.data$n_loci))

# Add info about raw reads
RawReads <- read_csv(file.path(get.value("demulti.log"),"AllIndividuals_Nreads.csv"))
head(RawReads)

cov.data <- cov.data %>% left_join(RawReads %>% select(sample = Filename, Raw = Retained)) %>% 
  left_join(pop.data %>% select(sample = ID_GQ, Espece, Notes_groupe, Region, Site, Shore))
cov.data %>% head()

#write_csv(cov.data, file.path(get.value("stacks.log"),"AllIndividuals_DeNovo.buccinum.M1n1_NreadsNloci.csv"))

cov.data <- read_csv(file.path(get.value("stacks.log"),"AllIndividuals_DeNovo.buccinum.M1n1_NreadsNloci.csv"))

graph1.0 <- cov.data %>%
  ggplot(aes(x = Raw/2, y = n_used_fw_reads, col = Notes_groupe)) +
  geom_point() +
  #facet_wrap(~Espece) %>% 
  theme_bw()
graph1.0

# Check which one are "outliers"

cov.data %>%
 # dplyr::filter(Region =="Newfoundland") %>% 
  ggplot(aes(x = Raw/2, y = n_used_fw_reads, col = Notes_groupe)) +
  geom_point() +
  geom_text( 
    data=cov.data %>% filter(#Region =="Newfoundland", 
                             n_used_fw_reads < 1250000,
                             Raw > 5000000), # Filter data first
    aes(label=sample),
    nudge_y = 150000,
    nudge_x = 500000
  ) +
  theme_bw()

  cov.data$Region


cov.data %>%
  ggplot(aes(x = Raw/2, y = Site, col = Region)) +
  geom_violin()+
  geom_jitter(width = 0) +
  #facet_wrap(~Espece) %>% 
  theme_bw()

cov.data %>%
  ggplot(aes(x =  n_loci, y = Site, col = Region)) +
  geom_violin()+
  geom_jitter(width = 0) +
  #facet_wrap(~Espece) %>% 
  theme_bw()

cov.data %>% dplyr::filter(n_used_fw_reads < 10000)

pop.data%>% group_by(Notes_groupe, Numero_unique_groupe) %>% summarise(N = n())

graph1.1 <- cov.data %>%
  ggplot(aes(x = mean_cov_ns, y = n_loci, col = Notes_groupe)) +
  geom_point() +
  theme_bw()
graph1.1


graph1.2 <- cov.data %>%
  ggplot(aes(x = mean_cov_ns)) +
  geom_histogram()+
  theme_bw()
graph1.2

pdf(file.path(get.value("stacks.log"),"BasicGraph_PostGstacks_DeNovo.buccin.M1n1_2024-06-12.pdf"))
print(graph1.0)
print(graph1.1)
print(graph1.2)
dev.off()

