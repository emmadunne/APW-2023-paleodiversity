# ******************************************************
#
#   Analytical Paleobiology Workshop 2023
#
#   Module 2: Paleodiversity analyses
#   Day 4 | Thursday, August 24th
#
#   Emma Dunne (emma.dunne@fau.de)
#   Lisa Schnetz (lisa.schnetz@gmail.com)
# ______________________________________________________
#
#   3. Sampling standardisation
# 
# ******************************************************


## In this script we will explore 2 methods of sampling standardisation:
##    1. Coverage-based rarefaction (similar to Shareholder Quroum Subsampling, SQS)
##    2. Squares, another coverage-based method of sampling standardisation
## Before getting into these methods, we will first explore our data using
##    Good's u and Rarefaction/Extrapolation curves curves via the iNEXT package



# 0. Packages used in this script -----------------------------------------

library(tidyverse)

require(devtools)
install_version("iNEXT", version = "2.0.20")# the most up-to-date version of iNEXT is still a little glitchy, so we'll use this one for now
library(iNEXT)




# Good's u ----------------------------------------------------------------

## We will start off by examining Good's u, which is ued by the traditional SQS 
##    method as a measure of sample 'coverage'. In iNEXT the equations of Chao & Jost
##    do this job - see notes below for more information

## Start off by loading your *cleaned* occurrence data (if not already in R form earlier):
#occ_data <- read_csv("./data/PBDB_pseudos.csv")

## Pare this down to just the columns we need, while creating a new data object so 
##  we don't overwrite the original. We'll use genus level data for these analyses:
genus_data <- subset(occ_data, select=c(genus, accepted_name, occurrence_no, collection_no,
                                        early_interval, late_interval, min_ma, max_ma))

## And load your data information if its not already in R:
#intervals <- read_csv("./data/intervals_Car_Tor.csv")


N <- sum(Individuals) # N
n1 <- length(which(Individuals==1))
u <- 1- n1/N # Good's u


## To compute Good's u for each interval, we need to know the frequencies of each taxon (genus):
tax_freq <- lapply(1:nrow(intervals), function(i) {
  tmp <- genus_data %>% filter(max_ma >= intervals[i,"max_ma"] & min_ma <= intervals[i,"min_ma"]) %>% 
    count(., genus) %>% arrange(desc(n)) %>% # count no. genera in each interval
    select(n)
  freq_raw <- as.numeric(tmp$n)
  freq_raw
})
names(tax_freq) <- intervals$interval_name # give each list element its correct interval name
str(tax_freq) # call up the list to see a summary of its structure


## This is the function to calculate Good's u. 
## You only need to run it once, but you can tweak some parts of it, depending on what you are interested in
goodsU <- function(occvec) {
  n <- sum(occvec) # sum of all occurrences
  f1 <- sum(occvec == 1); f2 <- sum(occvec == 2) # sum of singletons vs. doubletons
  # out <- 1 - (f1 / sum(occvec)) # singleton-only coverage estimator
  out <- 1 - length(which(occvec==1)) / n
  #out <- 1 - (f1/n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2)) # Chao and Jost (2012) coverage estimator using singletons and doubletons
  out[is.nan(out) | is.infinite(out)] <- NA
  return(out)
}

## Now let's compute Good's u for each interval:
intervals_Gu <- lapply(tax_freq, goodsU)
intervals_Gu <- as.data.frame(do.call(rbind, intervals_Gu)) # convert to dataframe 
intervals_Gu # pull results up the the console

## Good's u runs from 0 (totally incomplete) to 1 (very 'complete')
## How 'complete' are your intervals?



# 2. Data set-up for iNEXT ------------------------------------------------

## Now let's explore with iNEXT. This is a very comprehensive and flexible
##    package - we will only be scratching the surface here. Take a look at the
##    publicaito and associated documentation/tuturials for more info
citation("iNEXT") # https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12613

## First, we need to get our data set up properly:
##   1. a dataframe of interval names
##   2. total number of occurrences in each interval (i.e. sampling units)
##   3. taxon incidence frequencies for each interval

## In the incidence frequencies, the first entry of each list object is the total number of 
##    sampling units, followed by the taxa incidence frequencies
##    For example: Interval_A : 150 99 96 80 74 68 60 54 46 45
##      = there are 150 taxa in Interval_A, 99 in the first collection, 96 in the second, etc.

## 1. We've already got our intervals information from earlier, load that now if its not already in R:
#intervals <- read_csv("./data/intervals_Car_Tor.csv")

## 2 + 3. Genus incidence frequencies for each interval
## This loop computes incidence frequences for each interval in the intervals data:
freq_data <- lapply(1:nrow(intervals), function(i) {
  tmp <- genus_data %>% dplyr::filter(max_ma >= intervals[i,"max_ma"] & min_ma <= intervals[i,"min_ma"]) %>% 
    count(., genus) %>% arrange(desc(n)) %>% # count no. genera in each interval
    add_row(n = sum(.$n), .before = 1) %>% # find total number of genera
    select(n)
  freq_raw <- as.numeric(tmp$n) 
  freq_raw
})
names(freq_data) <- intervals$interval_name # give each list element its correct interval name
str(freq_data) # call up the list to see a summary of its structure
freq_data[[6]] # check the data for the Norian

## Now we're ready to go!




# 3. Rarefaction/Extrapolation curves -------------------------------------


## The main function of the package is iNEXT() (which gives the package its name) and means
##    "Interpolation and extrapolation of Hill number with order q"
##    This function has a much longer and more complex description that you can find
##    in the original literature, but for our purposes here, it is computing diversity 
##    estimates, sample coverage estimates and related statistics for sample sizes

## Compute iNEXT for our incidence frequency data:
inc_data <- iNEXT(freq_data, q = 0, datatype = "incidence_freq")

## This returns the "iNEXT" object including different output lists. We need $iNextEst,
##    which shows size- and coverage-based diversity estimates along with related
##    statistics for a series of rarefied and extrapolated samples
cov_rare <- inc_data$iNextEst # extract these data from the output lists

## The function ggiNEXT() extends ggplot2 to help plot the output.
?ggiNEXT # put up the help tab

## It allows 3 types of curves - see the documentation for more information
##  1. Sample-size-based R/E curve (type=1) - plots diversity estimates with confidence 
##      intervals (if se=TRUE) as a function of sample size up to double
ggiNEXT(inc_data, type=1, facet.var="site", se=TRUE) + theme_minimal()

##  2. Sample completeness curve (type=2) with confidence intervals - plots 
##      the sample coverage with respect to sample size 
ggiNEXT(inc_data, type=2) + theme_minimal()

##  3. Coverage-based R/E curve (type=3) - plots the diversity estimates with 
##      confidence intervals as a function of sample coverage up to the maximum coverage
ggiNEXT(inc_data, type=3) + theme_minimal()


## To make a publication-worthy Coverage-based R/E curve, we can plot all of the intervals together
## First let's go some tidying up so we can colour the Triassic and Jurassic separently:
for(i in 1:length(cov_rare)) {
  cov_rare[[i]]$stage_int <- names(cov_rare)[i] # add stage_int
}
cov_rare <- do.call(rbind, cov_rare) %>% as_tibble() # convert to tibble for ease of plotting
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name[1:4]), "Period"] <- "Jurassic" # mark Jurassic stages
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name[5:7]), "Period"] <- "Triassic" # mark Jurassic stages

## Finally, set up the plot in ggplot
cov_rare_plot <- ggplot(data = cov_rare, aes(x = SC, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = stage_int, colour = Period, lty = method)) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(values=c("dotted", "solid", "dotdash")) +
  # Add hex codes for the periods (see here: https://github.com/crimeacs/Geochronological_Colors/blob/master/Geochronological_scale_HEX.pdf)
  scale_colour_manual(values = c("#34b2c9","#812B92")) +
  theme_minimal() +
  labs(x = "Coverage", y = "Species richness") +
  scale_x_continuous(limits = c(0, 1.05), expand=c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 80), expand=c(0,0), breaks = seq(0, 80, 20))
cov_rare_plot



# 4. Shareholder Quorum Subsampling (SQS) ---------------------------------


## Shareholder Quorum Subsampling (SQS) uses rank-order abundance to estimate diversity 
##    through subsampling at different degrees of sampling coverage, or 'quorum levels', 
##    and estimates diversity using a metric called Good's u. 

## We will use the package iNEXT, which estimates diversity using Hill numbers via 
##    subsampling with the equations of Chao and Jost (2012), and also implements extrapolation, 
##    using the Chao1 estimator - for more info, see Hsieh et al. (2016):


## Start by creating a vector of quorum levels that we want to compute
## 0.4 is considered the 'standard', but the fashion now is to plot multiple quorum levels
quorum_levels <- round(seq(from = 0.5, to = 0.8, by = 0.1), 1)

## And create a new list to store output
estD_output <- list() 


## estimateD() is the main function for running SQS in iNEXT
?estimateD # pull up more information to the help tab

## Compute estimateD() for each interval
## This make take a few moments to run, depending on the speed of your computer
for(i in 1:length(quorum_levels)) {
  # estimateD function running over each quorum level:
  estD_tmp <- estimateD(freq_data, datatype = "incidence_freq", base = "coverage", level = quorum_levels[i])
  # filter to the correct diversity estimates (order = 1):
  estD_tmp <- filter(estD_tmp, order == 1)
  # organise the output:
  estD_tmp$quorum_level <- quorum_levels[i]
  estD_tmp$mid_ma <- intervals$mid_ma
  # add the output to the newly created list:
  estD_output[[i]] <- estD_tmp
}

## The output is a 'list' object so we'll need to convert it into a dataframe and clean it up before plotting
estD_plotting <- bind_rows(estD_output) # binds rows of a list

## Ensure that the quorum level column is being treated as a 'factor' to avoid errors while plotting:
estD_plotting$quorum_level <- as.factor(estD_plotting$quorum_level)

## Create a colour gradient for as many colours as you have quorum levels:
col_gradient <- scales::seq_gradient_pal("#B7E3B6", "#0A6B09", "Lab")(seq(0, 1, length.out = 4))

## Set your interval boundaries:
int_boundaries <- c(237.0, 228.0, 208.5, 201.3, 199.3, 190.8, 182.7, 174.1)

## Construct the plot in ggplot
iNEXT_plot <- ggplot(estD_plotting, aes(x = mid_ma, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
  ## Each quorum level is called individually to be plotted:
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.3), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = col_gradient[1], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.4), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = col_gradient[2], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.5), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = col_gradient[3], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.6), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = col_gradient[4], alpha = 0.2) +
  ## Set our line and point sizes (and shapes):
  geom_line(linewidth = 1) +
  geom_point(aes(pch = method), size = 4.5) +
  scale_shape_manual(values=c(15, 16, 17)) +
  ## Add our colours, theme, and axes labels:
  scale_colour_manual(values = col_gradient) +
  scale_x_reverse(breaks = int_boundaries) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness") +
  theme_minimal()
iNEXT_plot # Call the plot to the plots tab


## Save a copy of the plot to the plots folder
ggsave("./plots/iNEXT_gen.pdf", plot = iNEXT_plot, 
       width = 30, height = 18, units = "cm")





# 5.Squares ------------------------------------------------------------------

## Squares is an extrapolator and coverage-based approach which performs well when the 
##    rank abundance distributions of samples are particularly strongly skewed 
##    (i.e. there are many rare taxa, as is often the case with fossil data)
## In other words, the Chao1 estimator from estimateD() in iNEXT above, assumes that  
##    abundance distributions are uniform. Squares attempts to correct this by not  
##    making this assumption
## You can find out more information in the original publication: Alroy (2020):
##    https://doi.org/10.1017/pab.2019.40


## Just like we did for the analysis above, we need to get the genus incidence frequencies 
##    for each interval, but this time we don't need the totals at the beginning of the string
freq_data_sq <- lapply(1:nrow(intervals), function(i) {
  tmp <- genus_data %>% filter(max_ma >= intervals[i,"max_ma"] & min_ma <= intervals[i,"min_ma"]) %>% 
    count(., genus) %>% arrange(desc(n)) %>% 
    select(n)
  freq_raw <- as.numeric(tmp$n)
  freq_raw
})
names(freq_data_sq) <- intervals$interval_name # give each list element its correct interval name
str(freq_data_sq) # call up the list to see a summary of its structure
freq_data_sq[[6]] # check the data for the Norian


## Using these data, you can plot rank abundances for each interval
## Let's look at the Norian:
barplot(freq_data_sq[["Norian"]], log = "y", # log the y-axis
        col = "#e94196", border = 0, space = 0,
        xlab = "Rank order", ylab = "Abundance")


## Now let's get estimating diversity using squares
## The code here is adapted from Dr Bethany Allen's tutorial which you can find here:
##    https://github.com/bethany-j-allen/sampling_bias_workshop/blob/master/Code/3_SQSAndSquares.R

## Start by making an empty vector to store the squares estimates in:
squares_list <- vector("numeric", length = 0)

## Then loop through each stage:
for(i in 1:length(freq_data_sq)) {
  # filter the rank-abundance string of the stage of interest:
  count_list <- freq_data_sq[[i]]
  # Compute counts:
  genus_count <- length(count_list) # total number of genera in the string:
  sing_count <- sum(count_list == 1) # umber of singletons (the number of genera in the string with an abundance of 1):
  ind_count <- sum(count_list) # number of individuals (the sum of the string):
  # Find the sum of the string squared (this is where the name comes from):
  sum_nsq <- sum(count_list^2)
  # Put input parameters into the squares equation:
  squares <- genus_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*genus_count)))
  # Because the equation is a fraction, you can get an estimate of infinity if all taxa are singletons. 
  # As a safeguard against that, if the estimate comes out as infinite, we replace it with the observed genus count:
  if(squares == Inf){squares <- genus_count}
  # Finally add the squares estimate to the end of the vector:
  squares_list <- append(squares_list, squares)
}

## Label the squares estimates with their midpoints:
to_plot <- data.frame(squares_list, intervals$mid_ma)
## Rename the mid_ma column to keep things tidy:
to_plot <- rename(to_plot, "mid_ma" = "intervals.mid_ma")

## And finally, plot these data!
squares_plot <- ggplot(to_plot, aes(x = mid_ma, y = squares_list)) +
  geom_line(colour = "#e94196", size = 1) + 
  geom_point(colour = "#e94196", size = 4, shape = 16) +
  scale_x_reverse(breaks = int_boundaries) + 
  labs(x = "Time (Ma)", y = "Squares genus diversity") +
  theme_minimal()
squares_plot # call to plot tab

## Save a copy of the plot to the plots folder
ggsave("./plots/squares_gen.pdf", plot = sauares_plot, 
       width = 20, height = 14, units = "cm")


## How do your results for coverage-rarified richness and squares compare/differ?


