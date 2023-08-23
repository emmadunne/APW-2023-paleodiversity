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
#   1. Accessing databases in R 
#       & manipulating imported data
# 
# ******************************************************


# 0. Packages used in this script -----------------------------------------

library(tidyverse) 
library(divDyn)
library(sepkoski)


## Clear R's environment before starting so you're working with a clean slate:
rm(list = ls())

## If you've been using a lot of different packages, some function names might be masked;
## this step ensures that the function 'select' is coming from the dplyr package (part of tidyverse)
select <- dplyr::select



# 1. Fetching data from packages -------------------------------------------


### (a) sepkoski

## This package allows easy access to Sepkoski's fossil marine animal genera 
## compendium (Sepkoski, 2002), ported from Shanan Peters' online database.
## More information here: https://github.com/LewisAJones/sepkoski

## Accessing the datasets
data("sepkoski_raw") # Sepkoski's raw fossil marine animal genera compendium (Sepkoski, 2002)
data("sepkoski") # Sepkoski's compendium with first and last appearance intervals updated to be consistent with stages from the International Geological Time Scale 2022
data("interval_table") # a table linking intervals in Sepkoski's compendium with the International Geological Time Scale 2022.

## Let's look at the data...
View(sepkoski_raw) # opens a new tab in RStudio

## What variables have we got?
glimpse(sepkoski_raw) # dplyr (tidyyverse function)
str(sepkoski_raw) # base R function

## Let's plot Sepkoski's famous curve
sepkoski_curve()

## Take a look at the help file to customise the plot
?sepkoski_curve
sepkoski_curve(fill = TRUE)



### (b) divDyn

## Let's check out the corals dataset within divDyn
## (You will learn a lot more about this package tomorrow with Adam :) )

data(corals) # 'attach' the data so that R can 'see' it

## Let's make it easier to view:
View(corals) # opens a new tab in RStudio

## We can also view the variables using these different functions:
glimpse(corals) # dplyr (tidyverse) function

## Let's look at the types of variables we have:
?corals # pull up the help file to look at the variables

class(corals$genus) # what type of data is in the GENUS column?
## What type of data is in...
class(corals$max_ma)
class(corals$ecologyBoth)




# 2. Importing PBDB download ----------------------------------------------

## Any dataset can be imported into R from a file
## Let's download some data from the PBDB using the Download Generator 
## Where are you going to store this file?

## Now let's import it:
pbdb_data_raw <- read.csv("./data/pbdb_data.csv", skip = 18) 

## Take a look inside:
View(pbdb_data_raw)
glimpse(pbdb_data_raw)




# 3. PBDB via URL ---------------------------------------------------------

## The Paleobiology Database data can accessed through an API request
## Note that this requires an internet connection

## First, choose a taxonomic group and time interval and create new objects:
taxon_group <- "Pseudosuchia" # Taxon group
start_interval <- "Carnian" # Interval to start at
stop_interval <- "Toarcian" # Interval to stop at

## Create an API request form the PBDB and store this URL as an object
## A list of API options can be found here: https://paleobiodb.org/data1.2/
URL <- paste0("https://paleobiodb.org/data1.2/occs/list.csv?base_name=", # occurrence data, as a .csv
              taxon_group, "&interval=", start_interval, ",", stop_interval, # use our inputs from above
              "&show=full&pres=regular") # any additional columns we want 

## Then use this to load the data into R:
occ_data_raw <- as_tibble(read.csv(URL, header = TRUE, stringsAsFactors = FALSE))

## Take a peep:
glimpse(occ_data_raw) # view columns
View(occ_data_raw) # open as new tab

## It's good practice to save copies of your data as you go:
write_csv(occ_data_raw, "./data/PBDB_pseudos_24_08_23.csv")



# 4. Cleaning occurrence data ---------------------------------------------

## Raw occurrence data is imperfect, especially if you have not curated it yourself 
## 'Cleaning' is a very  important step to ensure you don't include unnecessary info

## Let's go through step by step and out some of the noise:

## Remove 'super-generic' identifications, so that we only retain occurrences to species- and genus-level
occ_data_raw2 <- filter(occ_data_raw, (identified_rank %in% c("species","genus")))

## Remove occurrences with “aff.”, “ex. gr.”, “sensu lato”, “informal”, or quotation marks in their identified names
## For PBDB data, why don't we use the "accepted_name" column for this?
occ_data_raw3 <- occ_data_raw2 %>% 
  filter(!grepl("cf\\.|aff\\.|\\?|ex\\. gr\\.|sensu lato|informal|\\\"", identified_name)) 

## Remove ichnotaxa (trace fossils) so that only regular taxa remain
## We can do this via the pres_mode column and entries marked as 'trace' or 'soft', and those with no genus name
occ_data_raw4 <- occ_data_raw3[occ_data_raw3$pres_mode != "trace", ] # trace taxa
occ_data_raw5 <- occ_data_raw4[!grepl("soft",occ_data_raw4$pres_mode), ] # 'soft' preservation

## Remove entries without a genus name - in the PBDB these can be errors
occ_data_raw6 <- occ_data_raw5[occ_data_raw5$genus != "", ]

## Finally, filter the data so any duplicate taxon names or collection numbers are eliminated:
occ_data <- distinct(occ_data_raw6, accepted_name, collection_no, .keep_all = TRUE)

## Take a look at the end result - How much has our data been reduced by?
length(unique(occ_data_raw$occurrence_no)) # start
length(unique(occ_data$occurrence_no)) # finish

## It's important to check your dataset for errors. In the case of fossil occurrence data, 
## errors can come from taxonomy, statigraphy, geography, etc.




# 5. Intervals data -------------------------------------------------------


## We can also grab time intervals data from the PBDB API. We'll do this here to help us 
## plotting later. However, these intervals data are quite out of date

## Download names and ages of time intervals from the PBDB:
intervals_all <- read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")
View (intervals_all) # take a look

## Make a vector of stage names that we are interested in:
interval_names <- c("Carnian", "Norian", "Rhaetian", # Late Triassic
                    "Hettangian", "Sinemurian", "Pliensbachian", "Toarcian") # Early Jurassic

## Select these intervals from the full PBDB intervals dataset:
intervals <- filter(intervals_all, interval_name %in% interval_names)

## Pare this down to just the 3 columns we'll need:
intervals <- select(intervals, interval_name, early_age, late_age)

## For ease of use later, let's rename the age columns to match the occurrence data:
intervals <- rename(intervals, "max_ma" = "early_age", "min_ma" = "late_age")

## And finally, calculate the midpoint for each interval and add it to a new (4th) column
intervals$mid_ma <- (intervals$min_ma + intervals$max_ma)/2

## Take a peep:
View(intervals) # open as new tab

## Save a copy as a .csv file - Note: your file path will differ!
write_csv(intervals, "./data/intervals_Car_Tor.csv")



