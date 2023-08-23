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

## Raw occurrence data is imperfect, especially if you have not curated it
## yourself 

