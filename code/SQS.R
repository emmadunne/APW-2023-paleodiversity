## This script will walk you through two coverge-based methods of sampling standardisation,
## (a) Shareholder Quorum Subsampling (SQS) via iNEXT and (b) squares, following the  
## methods outlined in Dunne et al. (2018) and Allen et al. (2020)


## Packages used in this script:
library(tidyverse)

require(devtools)
install_version("iNEXT", version = "2.0.20")
library(iNEXT)



# (a) Shareholder Quorum Subsampling ------------------------------------------

## Shareholder Quorum Subsampling (SQS) uses rank-order abundance to estimate diversity 
## through subsampling at different degrees of sampling coverage, or 'quorum levels', 
## and estimates diversity using a metric called Good's u. 

## We will use the package iNEXT, which estimates diversity using Hill numbers via 
## subsampling with the equations of Chao and Jost (2012), and also implements extrapolation, 
## using the Chao1 estimator - for more info, see the main citation:
citation("iNEXT")
## Or the paper describing the package: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12613


## Start off by loading your occurrence data (if not already in R form earlier):
occ_data <- read_csv("./03_sampling/PBDB_pseudos.csv")

## Pare this down to just the columns we need, while creating a new data object so we don't overwrite the original:
genus_data <- subset(occ_data, select=c(genus, accepted_name, occurrence_no, collection_no,
                                        early_interval, late_interval, min_ma, max_ma))


## To get your data in the right shape for iNEXT, you'll need:
##   1. a dataframe of interval names
##   2. total number of occurrences in each interval (i.e. sampling units)
##   3. genus incidence frequencies for each interval

## In the , the first entry of each list object is the total number of 
## sampling units, followed by the taxa incidence frequencies
## For example: Interval_A : 150 99 96 80 74 68 60 54 46 45
## = there are 150 taxa in Interval_A, 99 in the first collection, 96 in the second, etc.

## 1. We've already got our intervals information from earlier, load that now if its not already in R:
intervals <- read_csv("./data/intervals_Car_Tor.csv")

## 2 + 3. Get the genus incidence frequencies for each interval
freq_data <- lapply(1:nrow(intervals), function(i) {
  tmp <- genus_data %>% filter(max_ma >= intervals[i,"max_ma"] & min_ma <= intervals[i,"min_ma"]) %>% 
    count(., genus) %>% arrange(desc(n)) %>% 
    add_row(n = sum(.$n), .before = 1) %>%
    select(n)
  freq_raw <- as.numeric(tmp$n)
  freq_raw
})
names(freq_data) <- intervals$interval_name # give each list element its correct interval name
str(freq_data) # call up the list to see a summary of its structure
freq_data[[7]] # check the data for the Norian



## Now that we've got our data in the right shape, let's do some estimates!

## Create a vector of quorum levels that we want to compute
## 0.4 is considered the 'standard', but the fashion now is to plot multiple quorum levels
quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)

## And create a new list to store output
estD_output <- list() 

## estimateD() is the main function for running SQS in iNEXT

for(i in 1:length(quorum_levels)) {
  # estimateD function in iNEXT, running over each quorum level
  estD_tmp <- estimateD(freq_data, datatype = "incidence_freq", base = "coverage", level = quorum_levels[i])
  # filter to the diversity estimates (order = 1):
  estD_tmp <- filter(estD_tmp, order == 1)
  # organise the output:
  estD_tmp$quorum_level <- quorum_levels[i]
  estD_tmp$mid_ma <- intervals$mid_ma
  # add the output to the newly created list
  estD_output[[i]] <- estD_tmp
}

## The output fis a 'list' object so we'll need to convert it into a dataframe and clean it up before plotting
estD_plotting <- bind_rows(estD_output) # binds rows of a list

## Ensure that the quorum level column is being treated as a 'factor' to avoid errors while plotting:
estD_plotting$quorum_level <- as.factor(estD_plotting$quorum_level)

## Create a colour gradient for as many colours as you have quorum levels:
teal_gradient <- scales::seq_gradient_pal("turquoise", "darkslategrey", "Lab")(seq(0, 1, length.out = 4))

## Set your interval boundaries:
int_boundaries <- c(237.0, 228.0, 208.5, 201.3, 199.3, 190.8, 182.7, 174.1)

iNEXT_plot <- ggplot(estD_plotting, aes(x = mid_ma, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
  ## Each quorum level is called individually to be plotted:
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.3), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[1], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.4), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[2], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.5), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[3], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.6), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[4], alpha = 0.2) +
  ## Set our line and point sizes (and shapes):
  geom_line(linewidth = 1) +
  geom_point(aes(pch = method), size = 4.5) +
  scale_shape_manual(values=c(15, 16, 17)) +
  ## Add our colours, theme, and axes labels:
  scale_colour_manual(values = teal_gradient) +
  #scale_x_reverse(breaks = int_boundaries) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness") +
  theme_minimal()
iNEXT_plot # Call the plot to the plots tab


## Save a copy of the plot to the plots folder
ggsave("./plots/iNEXT_gen.pdf", plot = iNEXT_plot, 
       width = 30, height = 18, units = "cm")



inc.data <- iNEXT(freq_data, q = 0, datatype = "incidence_freq")
cov_rare <- inc.data$iNextEst
for(i in 1:length(cov_rare)) {
  cov_rare[[i]]$stage_int <- names(cov_rare)[i]
}

cov_rare <- do.call(rbind, cov_rare) %>% as_tibble() #convert to tibble for ease of plotting
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name[1:4]), "Period"] <- "Jurassic"
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name[5:7]), "Period"] <- "Triassic"

cov_rare_plot <- ggplot(data = cov_rare, aes(x = SC, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = stage_int, colour = Period, lty = method)) +
  geom_line(size = 1) +
  scale_linetype_manual(values=c("dotted", "solid", "dotdash")) +
  scale_colour_manual(values = c("#67A599","#F04028")) +
  #geom_point(data = cov_rare, aes(x = SC, y = qD, pch = Period, colour = Period), size = 3, inherit.aes = F) +
  theme(panel.background = element_blank(),
        legend.position="none",
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=14)) +
  labs(x = "Coverage", y = "Species richness") +
  scale_x_continuous(limits = c(0, 1.05), expand=c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 150), expand=c(0,0), breaks = seq(0, 210, 30))
cov_rare_plot
cov_rare_plot + geom_dl(data=cov_rare1, aes(label=stage_int),method=list("last.points",rot=30))







