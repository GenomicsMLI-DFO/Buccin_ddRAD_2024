#############
## Code to plot treemix output for Buccinum undatum
## June 2025
#############
#library(R.utils)
library(here)
#here::i_am("Buccinum_treemix_plots.R")


#################
#Step 0: treemix runs on Unix systems. Our input file for treemix
# is called "treemix_buccinum_input.gz"
#Here we will call various output files and plot results

#Now read in treemix files run on linux computer
#code to plot

source("01_Codes/Functions/plotting_funcs.R")
#Plot unrooted tree with no migration events
plot_tree(here("02_Results/04_treemix_output/", "treemix_output", "out_stem0"), mbar=FALSE) #this one looks good! But remove legend from plot since no migration events


#plot residuals
plot_resid(
  here("02_Results/04_treemix_output/","treemix_output", "out_stem0"),
  here("02_Results/04_treemix_output/","treemix_output", "treemix_populations.txt")
)
#Poor fits include: Iles de la madeleine poor fit to tree, bai comeau 1 and 2 and ille penche/bic

#now plot rooted version with one migration event
plot_tree(here("02_Results/04_treemix_output/","treemix_output","out_stem1")) #adds event from bic across egsl
plot_resid(
  here("02_Results/04_treemix_output/","treemix_output","out_stem1"),
  here("02_Results/04_treemix_output/","treemix_output", "treemix_populations.txt") )#still have mageleine as poor fit as well as NL and parts of ESL

#now plot rooted version with two migration events
plot_tree(here("02_Results/04_treemix_output/","treemix_output","out_stem2")) #Nonsensicle event from 3l --> BoF
plot_resid(
  here("02_Results/04_treemix_output/","treemix_output","out_stem2"),
  here("02_Results/04_treemix_output/","treemix_output","treemix_populations.txt")) #seems to lead to poorer fit + mingan very poor fit

#now plot rooted version with three migration events
plot_tree(here("02_Results/04_treemix_output/","treemix_output","out_stem3")) #2 across EGSL and then one weird from Mingan --> BoF
plot_resid(
  here("02_Results/04_treemix_output/","treemix_output","out_stem3"),
  here("02_Results/04_treemix_output/","treemix_output","treemix_populations.txt")) #Mingan poor fit

#TAKE-HOME: migration events do not seem plausible
#but need to look at f3 statistics for events.
#Doing so shows lack of support, so we conclude
#tree without historical migration is most appropriate fit
