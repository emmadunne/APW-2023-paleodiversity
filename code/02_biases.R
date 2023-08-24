# ******************************************************
#
#   Analytical Paleobiology Workshop 2023
#
#   Module 2: Paleodiversity analyses
#   Day 4 | Thursday, August 24th
#
#   Emma Dunne (emma.dunne@fau.de)
# ______________________________________________________
#
#   2. Exploring fossil record biases
# 
# ******************************************************


# 0. Packages used in this script -----------------------------------------

library(tidyverse)
library(geoscale) # for plotting with the geological time scale on the x-axis (uses base R syntax)
library(viridis) # for colour scales
library(vegan) # for diversity metrics

select <- dplyr::select # ensure the select function is coming from dplyr



# 1(a). Sampling proxy counts ---------------------------------------------------


## Let's explore sampling patterns!
## First we'll calculate counts of sampling proxies and plot these alongside raw diversity

# Taxa per interval 
count_taxa <- vector("numeric") # create empty vector for the loop below to populate
for (i in 1:nrow(intervals)) { # for-loop to count each taxon that appears in each interval
  out <- subset(occ_data, max_ma > intervals[i,]$min_ma & min_ma < intervals[i,]$max_ma) # uses our intervals dataframe
  count_taxa[i] <- (length(unique(out$accepted_name)))
  print(count_taxa[i])
}

# Collections per interval
count_colls <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- subset(occ_data, max_ma > intervals[i,]$min_ma & min_ma < intervals[i,]$max_ma)
  count_colls[i] <- (length(unique(out$collection_no)))
  print(count_colls[i])
}

# Formations per interval
count_formations <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- subset(occ_data, max_ma > intervals[i,]$min_ma & min_ma < intervals[i,]$max_ma)
  count_formations[i] <- (length(unique(out$formation)))
  print(count_formations[i])
}


## For equal-area gird cells, I would recommend the package 'icosa' (Kocsis, 2017)
## For more info see: http://cran.nexr.com/web/packages/icosa/vignettes/icosaIntroShort.pdf


## Gather the proxy information together in a new dataframe for plotting:
proxy_counts <- data.frame(intervals$interval_name, intervals$mid_ma, count_taxa, count_colls, count_formations)
## Rename the columns for ease:
proxy_counts <- rename(proxy_counts, 
                       "interval_name" = "intervals.interval_name", 
                       "mid_ma" = "intervals.mid_ma")

## Finally, convert all zero's to NAs for plotting 
## This means that the plots won't register zero and instead will leave gaps 
## where there is no sampling instead - this gives a more realistic picture
proxy_counts[proxy_counts == 0] <- NA 



# 1(b). Sampling plots ----------------------------------------------------------


## Let's get plotting these patterns!
## Below are two options depending on your preferred aesthetic
##    1. via ggplot
##    2. using geoscale to plot a geological time scale on the x-axis (similar to 'deeptime')


## Option 1: Plotting using ggplot
##_________________________________

## Set interval boundaries for the dotted lines on the plot
## We'll also use this vector again, so its handy to have :)
int_boundaries <- c(237.0, 228.0, 208.5, 201.3, 199.3, 190.8, 182.7, 174.1)

## Set up your ggplot layers (first layer goes on the bottom, etc):
proxy_plot <- ggplot() + 
  # Formations (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_formations), colour = "orangered3", linewidth = 1.2, linetype = "dashed")  +
  geom_point(data = proxy_counts, aes(mid_ma, count_formations), colour = "orangered3", size = 4, shape = 16) +
  # Collections (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_colls), colour = "peru", linewidth = 1.2, linetype = "dashed")  +
  geom_point(data = proxy_counts, aes(mid_ma, count_colls), colour = "peru", size = 5, shape = 16) +
  # Taxa (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_taxa), colour = 'black', linewidth = 1.2)  +
  geom_point(data = proxy_counts, aes(mid_ma, count_taxa), colour = "black", size = 4, shape = 16) +
  # Add a minimal theme - but you can make your own custom themes too!
  theme_minimal() + 
  labs(x = "Time (Ma)", y = "Sampling proxy counts") +
  # Make sure to reverse the x-axis to match geological time!
  scale_x_reverse(breaks = int_boundaries) +
  # And tidy up our y-axis with even breaks that match the totals in our dataframe:
  scale_y_continuous(breaks = seq(0, 320, 20))
## Call the finished plot to the RStudio plots tab:
proxy_plot

## Set dimensions and save plot (as pdf) to the plots folder
#dir.create("./plots") # create new folder if one doesn't already exist
ggsave(plot = proxy_plot,
       width = 20, height = 15, dpi = 500, units = "cm", 
       filename = "./plots/sampling_proxies.pdf", useDingbats=FALSE)



## Option 2: Plotting using geoscale
##__________________________________

## NB: this is base R syntax, not ggplot syntax, so its a little different!
## Take a look at the helpfile for the geoscale package for more info on each function
## RStudio > Help panel (bottom right) > Search "geoscalePlot"

## Before you make the plot, set up parameters for exporting a PDF of the plot to:
pdf("sampling_proxies_geoscale.pdf", width = 9, height = 7) 

## Set up the base of the plot with the timescale:
geoscalePlot(proxy_counts$mid_ma, proxy_counts$count_taxa, # ages and main data points
             units = c("Period", "Age"), # which intervals to show in your timeline  (can also add Epoch)
             tick.scale = "Epoch", # resolution of the tick marks on timescale
             boxes = "Age", # option to include grey boxes for individual time bins
             abbrev = c("Age"), # option to abbreviate names of geological units
             lty = 1, pch = NA, 
             cex.age = 1, cex.ts = 1, cex.pt = 1, # size of numbers, text, and points on the scale bar
             age.lim = c(237.0, 174.1), # oldest and youngest ages of the entire time interval you're plotting (i.e. temporal range)
             data.lim = c(0, 320), # range of the data you're plotting - check it matches the largest number in your dataframe!
             ts.col = TRUE, # include colours in the timescale? 
             ts.width = 0.17, # space taken up by plotting the time scale (value = 0-1)
             label = "Counts", # label for y-axis
             direction ="horizontal", # orientation of the plot
             erotate = 0) # numerical value for the rotation for the temporal units (default is 90)

# Now add the different lines and points separately:
# Formations (as dots and a line):
lines(proxy_counts$mid_ma, proxy_counts$count_formations, col="peru", lty = 2, lwd = 3)
points(proxy_counts$mid_ma, proxy_counts$count_formations, pch =16, cex = 1.5, col = "peru")
# Collections (as dots and a line):
lines(proxy_counts$mid_ma, proxy_counts$count_colls, col="orange", lty = 3, lwd = 3)
points(proxy_counts$mid_ma, proxy_counts$count_colls, pch =16, cex = 1.5, col = "orange")
# Taxa (as dots and a line):
lines(proxy_counts$mid_ma, proxy_counts$count_taxa, col="black", lty = 1, lwd = 3)
points(proxy_counts$mid_ma, proxy_counts$count_taxa, pch =16, cex = 1.5, col = "black")

# Finally, add a legend
legend('topright', legend = c("Formations", "Collections", "Taxa"), 
       col = c("peru", "orange", "black"), 
       pch = c(16), box.lty = 1, pt.cex=2, bg= "white")

dev.off() ## Turn off graphic device, to trigger the export of the pdf
## Sometimes this doesn't look good in the RStudio window so you might need 
## to navigate to the plots folder to see it properly



# 2. Collections per latitude ------------------------------------------------

## Yesterday, you had a look with alpha diversity ('local richness') with Wolfgang
##    When visualised, alpha diversity can provide more insights into sampling patterns
##    especially as it adds a spatial element as opposed to just temporal patterns.
## Let's create a plot to see where the collections across time and paleolatitude.
##    We'll also colour our plot according to the number of taxa in each collection

## There is evidence to suggest that alpha diversity is not as strongly affected by sampling biases
##    as gamma (or 'global') diversity. For a more sophisticated way to calculate alpha diversity by
##    treating taxonomically indeterminate occurrences as valid, see the method described in 
##    Close et al. (2019) - code available here: https://github.com/emmadunne/local_richness

## Let's get our data set up:
lat_data <- occ_data # rename object to keep the original separate

## Create new column for mid_ma
lat_data$mid_ma <- (lat_data$max_ma + lat_data$min_ma)/2 

## Next, we'll need to count the number of taxa per collection (i.e. their frequency):
taxa_freqs <- count(lat_data, collection_no)

## Subset lat_data to only the columns we need:
lat_data <- lat_data %>% 
  select(collection_no, paleolat, paleolng, mid_ma) %>% 
  distinct() %>% na.omit()

## Add add the frequency information:
lat_data <- left_join(taxa_freqs, lat_data, by = "collection_no")

## Before we plot, let's order the frequencies and remove any NAs that have crept in:
lat_data <- lat_data %>% arrange(n) %>% na.omit()

## Take a look:
View(lat_data)


## Set up our ggplot layers
lat_plot <- ggplot(data = lat_data, aes(x = mid_ma, y = paleolat, colour = n)) +
  geom_vline(xintercept = int_boundaries, lty = 2, col = "grey90") +
  geom_hline(yintercept = 0, colour = "grey10") +
  scale_color_viridis(trans = "log", breaks = c(1, 2, 12), direction = -1, option = "D") + # set the break= to match your richness data
  #scale_y_continuous(labels = function(x) format(x, width = 5), limits = c(-70, 70), breaks = seq(from = -60, to = 60, by = 20)) +
  scale_x_reverse(breaks = int_boundaries) + 
  theme_minimal() + 
  theme(legend.direction = "vertical", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Palaeolatitude (º)") +
  geom_point(size = 4, alpha = 0.5) # (alpha sets point transparency)
lat_plot # call to plot window


## Set dimensions and save plot (as pdf)
ggsave(plot = lat_plot,
       width = 20, height = 10, dpi = 500, units = "cm", 
       filename = "./plots/lat_alpha_div.pdf", useDingbats=FALSE)



# 3. Regression -----------------------------------------------------------

## We can also apply some simple statistics to better quantify the
##    relationship between raw richness and sampling proxies/effort
## Let's do some simple regression plots:

## Raw richness vs. collections
reg_colls <- ggplot(proxy_counts_NC, aes(x=count_taxa, y=count_colls)) + 
  geom_point(shape=17, size = 6, colour = "orange")+
  geom_smooth(method=lm, colour = "orange4", fill = "orange1")  +
  theme_minimal()
reg_colls

## Raw richness vs. formations
reg_forms <- ggplot(proxy_counts_NC, aes(x=count_taxa, y=count_formations)) + 
  geom_point(shape=16, size = 5, colour = "#30C430")+
  geom_smooth(method=lm, colour = "#0A6B09", fill = "#B0ECB0") +
  theme_minimal()
reg_forms


## Let's quantify these relationships through a liner model:

## Raw richness vs. collections
lm_colls = lm(count_colls ~ count_taxa, proxy_counts_NC)
summary(lm_colls) # summary of results

## Raw richness vs. formations
lm_forms = lm(count_formations ~ count_taxa, proxy_counts_NC)
summary(lm_forms)

## We can test the senstivity of the data to well-sampled intervals
##    In this case, we could remove the Norian and/or the Carnian and re-run the analyses

# Take out row 6, which is the Norian
proxy_counts_N <- proxy_counts[-6,]
# Take out row 6 + 7 (Norian + Carnian)
proxy_counts_NC <- proxy_counts[-c(6,7),]

## Are there any changes to your results?



# 4. Collector's curve -------------------------------------------------------

## Collector's curves, also known as a species accumulation curves, can also tell us about sampling: 
##    they display the cumulative number of taxa as a function of the cumulative effort (i.e. sampling proxy)

## First, get our data into shape:
## Table the number of each species per collection
abun_data <- table(occ_data$collection_no, occ_data$accepted_name) # abundance
#abun_data[abun_data > 0] <- 1 # if we want a presence/absence table, we can add this step

## Turn this into a matrix:
abun_matrix <- matrix(abun_data, ncol = length(unique(occ_data$accepted_name)))
colnames(abun_matrix) <- colnames(abun_data) # add the column names back in for when we need to check anything
rownames(abun_matrix) <- rownames(abun_data) # same for the row names

## Using the vegan package, we can make the collectors or species accumulation curve
## Check out the help file for specaccum() to find out more about the methods
sp_accum <- specaccum(abun_matrix, method = "collector")
summary(sp_accum) # bring up the summary

## Plot the curve in base R - you can make this pretty in ggplot if you prefer!
plot(sp_accum, ci.type = "poly", col = "#0E6A8A", lwd = 2, ci.lty = 0, ci.col = "#5CBCDD")





# 5. Modern world map --------------------------------------------------------

## Finally, let's explore our data on a modern world map and see if we can spot
##    any geographic (and even socio-economic) patterns...

## First, let's pear down or occurrence data to only keep the info we need for making the map
locality_info <- occ_data %>% 
  dplyr::select(collection_name, lat, lng, early_interval, late_interval, max_ma, min_ma) %>% 
  distinct(collection_name, .keep_all = TRUE) %>% 
  na.omit()

## Grab a world map for ggplot to work with:
world_map <- map_data("world")
ggplot() + geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region)) 

## Let's make it pretty and add our data
modern_map <- ggplot() + 
  geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region), 
           color = "grey80", fill = "grey90", size = 0.1) +
  geom_point(data = locality_info, aes(lng, lat), alpha = 0.3, size = 4, colour = "#9B1999") +
  theme_void() + theme(legend.position = "none")
modern_map

## And save as a .pdf
ggsave(plot = modern_map,
       width = 8, height = 5, dpi = 600, 
       filename = "./plots/Modern_map.pdf", useDingbats=FALSE)


