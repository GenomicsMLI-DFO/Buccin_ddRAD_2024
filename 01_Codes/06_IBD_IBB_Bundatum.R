################################################
#
# Code to test for IBD and IBB in B. undatum
#
################################################


#load libraries
library(vegan)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(here)


#Part 1: IBD within EGSL northern shore
#7 sites; 21 unique comparisons

#Read in matrices of pairwise FST and geog distance (euclidean)
#Here, ordered: Iles Penchees -	Baier des Bacon	Point a Emilie	Forestville	Portneuf	Bai Comeau 1	Bai Comeau 2

fst <- read.table(here( "02_Results/03_IBD_IBB_Analysis/" , "matrices", "fst_egsl.txt"), header = FALSE)
dim(fst) #7 x 7

geog <- read.table(here("02_Results/03_IBD_IBB_Analysis/" ,"matrices", "geographic_northern_egsl.txt"), header=T)
geog[upper.tri(geog, diag = T)] <- NA #Matrix is symmetrical so convert the upper redundant triangle and the diagonal to NA
geog
dim(geog) #7 x 7

#Linearize FST and natural log distance

fstLin <- fst/(1 - fst)
ln_geog <- log(geog)

#test corr b/w lin fst and log distance
mant_linlog <- mantel (ln_geog, fstLin, permutation=999)
mant_linlog #Mantel statistic r:  0.8512  ; p = 0.002

#Can plot and do linear model for comparison with next part.

y = as.dist(fstLin)[1:length(as.dist(fstLin))]#y= linearized Fst
x_log = as.dist(ln_geog)[1:length(as.dist(ln_geog))] #x= geographic distance
#The above two lines of code are unravelling the matrices into numeric vectors that can then be compared and plotted against each other in a regression

ibd_data <- as.data.frame(cbind(x_log, y))

# Fit the model
ibd_reg <- lm(y ~ x_log, data = ibd_data)
summary(ibd_reg)

#Plot

# Create prediction data
pred_data1 <- data.frame(
  x_log = seq(min(ibd_data$x_log), max(ibd_data$x_log), length.out = 100)
)

# Get predictions and confidence intervals
predictions <- predict(ibd_reg, newdata = pred_data1, interval = "confidence")
pred_data1$y_pred  <- predictions[, "fit"]
pred_data1$y_lower <- predictions[, "lwr"]
pred_data1$y_upper <- predictions[, "upr"]


p1 <- ggplot() +
  # Data points (main points with legend)
  geom_point(data = ibd_data,
             aes(x = x_log, y = y),
             size = 3, alpha = 0.7, color= "#E69F00") +
  # Confidence bands (no legend)
  geom_ribbon(data = pred_data1,
              aes(x = x_log, y = y_pred,
                  ymin = y_lower, ymax = y_upper),
              alpha = 0.2, show.legend = FALSE,fill= "#E69F00" ) +
  # Predicted trendline (no additional legend)
  geom_line(data = pred_data1,
            aes(x = x_log, y = y_pred),
            color = "#E69F00", linewidth = 1.2) +
  theme_minimal() +
  labs(x = "Ln (Geographic distance km) ",
       y = expression(Linearized~F[ST]~"("~F[ST]~"/(1-"~F[ST]~")"~")"))+
  scale_y_continuous(limits = c( -0.00142, 0.0286)) +  # Set y-axis limit; comes from limits on next plot
  scale_color_manual(values = c("#E69F00")) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),  # Remove gridlines
        axis.line = element_line(color = "black", size = .5),  # Add solid axis lines
        axis.ticks = element_line(color = "black", size = .5))  # Keep tick marks
p1

### Part 2: Look at effect of adding (z) depth barrier

#Read in dataframe that has all pairwise values within EGSL (including "southern" shore) as
#well as BoF, though we will drop bof for now...
#Here a barrier is a deep-water barrier b/w sites

#read in geog distance and then take natural log
bar_df <- read.csv(here("02_Results/03_IBD_IBB_Analysis/" ,"matrices", "barrier_data_egsl_bof.csv"), header=T)
dim(bar_df); head(bar_df)
bar_df$ln_geog <- log(bar_df$geog)

#Let's remove Bof for now; focus within EGSL
bar_egsl = bar_df[1:36, 1:7]
str(bar_egsl) #here, barrier is integer

# Ensure 'barrier' is treated as a factor in the original data and the new data
bar_egsl$barrier <- as.factor(bar_egsl$barrier)
str(bar_egsl)

#Finally, linearize FST
bar_egsl$fstLin <- bar_egsl$FST/(1 - bar_egsl$FST)

# Fit the model
model2 <- lm(fstLin ~ ln_geog * barrier, data = bar_egsl)
summary(model2) #all terms and interaction significant

# Create prediction data for barrier = 0
pred_data_0 <- data.frame(
  ln_geog = seq(min(bar_egsl$ln_geog[bar_egsl$barrier == 0]), max(bar_egsl$ln_geog[bar_egsl$barrier == 0]), length.out = 100),
  barrier = factor(0, levels = levels(bar_egsl$barrier))  # barrier = 0
)

# Create prediction data for barrier = 1
pred_data_1 <- data.frame(
  ln_geog = seq(min(bar_egsl$ln_geog[bar_egsl$barrier == 1]), max(bar_egsl$ln_geog[bar_egsl$barrier == 1]), length.out = 100),
  barrier = factor(1, levels = levels(bar_egsl$barrier))  # barrier = 1
)

# Combine both prediction data sets
pred_data <- rbind(pred_data_0, pred_data_1)

# Get predictions with confidence intervals
predictions <- predict(model2, newdata = pred_data, interval = "confidence")

# Add predictions and confidence intervals to pred_data
pred_data$FST_pred  <- predictions[, "fit"]
pred_data$FST_lower <- predictions[, "lwr"]
pred_data$FST_upper <- predictions[, "upr"]

# Plot with separate trendlines and confidence bands
p2 <- ggplot() +
  # Data points (main points with legend)
  geom_point(data = bar_egsl,
             aes(x = ln_geog, y = fstLin, color = barrier),
             size = 3, alpha = 0.7) +
  # Confidence bands (no legend)
  geom_ribbon(data = pred_data,
              aes(x = ln_geog, y = FST_pred,
                  ymin = FST_lower, ymax = FST_upper, fill = barrier),
              alpha = 0.2, show.legend = FALSE) +
  # Predicted trendlines (no additional legend; same color as points)
  geom_line(data = pred_data,
            aes(x = ln_geog, y = FST_pred, color = barrier),
            linewidth = 1.2) +
  theme_minimal() +
  labs(x = "Ln (Geographic distance km)",
       y = expression(Linearized~F[ST]~"("~F[ST]~"/(1-"~F[ST]~")"~")"),
       color = "Barrier") +
  scale_color_manual(values = c("#E69F00", "#8CCF8C"), labels = c("No Barrier", "Barrier")) +
  scale_fill_manual(values = c("#E69F00", "#8CCF8C")) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),  # Remove gridlines
        axis.line = element_line(color = "black", size = .5),
        axis.ticks = element_line(color = "black", size = 0.5),
        legend.title = element_blank())  # Add solid axis lines
p2
plot_build <- ggplot_build(p2)

# Check the y-axis limits
plot_build$layout$panel_scales_y[[1]] #-0.00142 -- 0.0278

#combine plots


# Assuming p1 and p2 are already created
# Combine the plots with labels A) and B) in the top-left corners
combined_plot <- p1 + p2 +
  plot_layout(ncol = 2) +  # Arrange the plots side by side
  plot_annotation(tag_levels = "A")

# Show the combined plot
combined_plot

#Last thing is to do a partial mantel for just the EGSL sites

#As reminder, here are objects for regular 2 matrix mantel
#mant_linlog <- mantel (ln_geog, fstLin, permutation=999)
#need to create barrier matrix in same layout and expand other ones to include bic sites
#order of sites:
##Iles Penchees -	Baier des Bacon	Point a Emilie	Forestville	Portneuf	Bai Comeau 1	Bai Comeau 2 BIC1 BIC2

geog_pm <- read.table(here("02_Results/03_IBD_IBB_Analysis/" ,"matrices", "geog_egsl1-2a.txt"), header=F)
dim(geog_pm) #9 x 9
geog_pm_ln <- log(geog_pm)

fst_pm <- read.table(here("02_Results/03_IBD_IBB_Analysis/" ,"matrices", "fst_egsl1-2a.txt"), header=F)
dim(fst_pm)  #9 x 9
fst_pm_lin <- fst_pm/(1 - fst_pm)

barrier_pm <-  read.table(here("02_Results/03_IBD_IBB_Analysis/" ,"matrices", "barrier_egsl1-2a.txt"), header=F)
dim(barrier_pm) #9x9


mantel(fst_pm_lin, geog_pm_ln, method="pearson", permutations=999 )
#here relationship isn't as strong; r = 0.3692; p = 0.032 because cross-region
#comparisons included
mantel.partial(fst_pm_lin, geog_pm_ln, barrier_pm, method="pearson", permutations=999)
#partial mantel higher once account for barrier #r = 0.7653; p = 0.024
