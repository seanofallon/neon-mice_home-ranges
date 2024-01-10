####################################### LOAD AND CLEAN NEON SMALL MAMMAL DATA
# This script processes data products downloaded directly from the NEON website, specifically
# the files needed for the analysis performed in 'Uncovering multiple influences on space use
# by deer mice using NEON data' by S. O'Fallon, K. Mabry, and N. Pinter-Wollman.

# The INPUT for this script is the zip folder that downloads when you select the "Download
# Data" button at NSF-NEON site: https://data.neonscience.org/data-products/DP1.10072.001
# Because this download/processing step takes some time, the outputs of this script that are 
# needed as inputs for later steps in this study are included as csv files in this R Project/
# Github repo - so if you already have 'pejn_caps.csv' and 'mnka-bp.csv' in your working
# directory, you can safely skip to 'get_pe_home_ranges.R' to replicate the O'Fallon et al. study.

# However, if you'd like to start from scratch, or work with a different set of data than was
# used in that study, this script generates useful OUTPUTS for working with NEON small mammal
# capture data:
# mam_plotNight_nodups.csv - 
# mam_trapNight_nodups.csv
# captures_w_id.csv
# mnka-bp.csv
# pejn_caps.csv

# The script should work regardless of the filters you choose to put on your download, but for
# this project, all sites were selected

### Clear env, load packages
# clear
rm(list=ls())

# load packages
library(dplyr)
library(neonUtilities)
library(neonOS)

### Load data
# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

# # Just need to do this once to get csvs saved locally - then comment out
# # "Stacks" data downloaded in zip folder into usable format (note that you need path to where you downloaded)
# stackByTable("C:/Users/seano/Downloads/NEON_count-small-mammals.zip")

# # read csvs downloaded from NEON website
# categoricalCodes_10072 <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/categoricalCodes_10072.csv")
# issueLog_10072 <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/issueLog_10072.csv")
# mam_identificationHistory <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/mam_identificationHistory.csv")
mam_perplotnight <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/mam_perplotnight.csv")
mam_pertrapnight <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/mam_pertrapnight.csv")
# mam_voucher <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/mam_voucher.csv")
# validation_10072 <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/validation_10072.csv")
variables_10072 <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/NEON_count-small-mammals/stackedFiles/variables_10072.csv")

### This commented chunk removes duplicate data entries, which takes a long time given the size of the objects.
### To ensure this step is only needed once, write.csv is called to save the objects as files locally - so 
### this only needs to be done once. 
# ### Remove duplicate data entries
# # 1.Check the perplotnight table by nightuid using standard removeDups function
# mam_plotNight_nodups <- neonOS::removeDups(data=mam_perplotnight,
#                                            variables=variables_10072,
#                                            table='mam_perplotnight')
# # remove flagged entries - they now have duplicateRecordQF = 1 or 2
# mam_plotNight_nodups <- filter(mam_plotNight_nodups, duplicateRecordQF < 1)
# 
# # # write mam_plotNight_nodups to csv, so get_home_ranges script can run without these slow steps being repeated
# # write.csv(mam_plotNight_nodups, file = "C:/Users/seano/Desktop/Projects/NEONSmallMammals/mam_plotNight_nodups.csv")
# 
# 
# # 2. Filter out multiple captures of untagged individuals in a single trap
# # (trapStatus = 4) before running the removeDups function on the
# # mam_pertrapnight data.
# mam_trapNight_multipleCaps <- mam_pertrapnight %>%
#   filter(trapStatus == "4 - more than 1 capture in one trap" &
#            is.na(tagID) & is.na(individualCode))
# 
# # This data subset contains no multiple captures so no further filtering is necessary
# 
# #Check the pertrapnight table using standard removeDups function
# mam_trapNight_nodups <- neonOS::removeDups(data=mam_pertrapnight,
#                                            variables=variables_10072,
#                                            table='mam_pertrapnight')
# 
# # how many rows w/o dups?
# nrow(mam_trapNight_nodups)
# # filter out rows flagged by removeDups
# mam_trapNight_nodups <- filter(mam_trapNight_nodups, duplicateRecordQF < 1)
# 
# # write mam_trapNight_nodups to csv in data folder bc ridding duplicates takes forever
# write.csv(mam_trapNight_nodups, file = "C:/Users/seano/Desktop/Projects/NEONSmallMammals/mam_trapNight_nodups.csv")

# load in perPlotNight and perTrapNight tables w/ dups removed
mam_plotNight_nodups <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/mam_plotNight_nodups.csv")
mam_trapNight_nodups <- read.csv("C:/Users/seano/Desktop/Projects/NEONSmallMammals/mam_trapNight_nodups.csv")

# filter to just traps that captured a tagged animal (or animal that was then tagged)
captures_w_id <- filter(mam_trapNight_nodups, !is.na(tagID))
captures_w_id <- filter(captures_w_id, tagID != "")

### Ensure that trapStatus is correct by making sure all rows w/ positive trapStatus has tagID associated
# makes object w/ rows that have tagID but no "capture" status
trapStatusErrorCheck <- captures_w_id %>%
  filter(!is.na(tagID)) %>%
  filter(!grepl("capture",trapStatus))

nrow(trapStatusErrorCheck)
# There are 13 records that have a tagID without a captured trapStatus

# makes an object w/ rows that have NA in tagID but a "capture" status
tagIDErrorCheck <- captures_w_id %>%
  filter(is.na(tagID)) %>%
  filter(grepl("capture",trapStatus))

nrow(tagIDErrorCheck)
# There are 0 records that lack a tagID but are marked with a captured trapStatus

# filter out records where a tagID was recorded w/o a "capture" status
captures_w_id <- filter(captures_w_id, grepl("capture",trapStatus))

# write captures_w_id as csv to data folder bc its a convenient df
write.csv(captures_w_id, file = "C:/Users/seano/Desktop/Projects/NEONSmallMammals/captures_w_id.csv")


######################### CALCULATE PEROMYSCUS MNKA FOR ALL PLOTS AT ALL SAMPLING EVENTS
# filter just PEMA captures
pe_caps <- filter(captures_w_id, taxonID %in% c('PEMA','PELEPEMA','PESP',
                                                'PEBO','PEGO','PEGOPELE',
                                                'PEKE','PEKEPEMA','PELE',
                                                'PEMAPEBO','PETR','PEAT'))

# make ref df w/ each tagID and # captures (rows) to use to filter main object
unq_ids <- unique(pe_caps$tagID)

# combine pe_caps and plotNight table to associate eventID's w/ tag
pejn_caps<-neonOS::joinTableNEON(mam_plotNight_nodups,
                                   pe_caps, name1 = "mam_perplotnight",
                                   name2 = "mam_pertrapnight")
pejn_caps<-filter(pejn_caps,tagID!='NA')

#Filter trap dataset to just the capture records of target taxa and a few core
# fields needed for the analyses.
coreFields <- c("nightuid","plotID","collectDate.x","tagID","eventID")


caps_for_mnka <- pejn_caps %>%
  dplyr::select(all_of(coreFields)) %>%
  rename('collectDate' = 'collectDate.x')

caps_for_mnka<-na.omit(caps_for_mnka)

# Generate a column of all of the unique tagIDs included in the dataset
uTags <- caps_for_mnka %>% dplyr::select(tagID) %>% distinct()

#create empty data frame to populate
known_alive <- slice(caps_for_mnka,0)

# For each tagged individual, add a record for each night of trapping done on
# the plots on which it was captured between the first and last dates of capture
for (i in uTags$tagID){
  indiv <- caps_for_mnka %>% filter(tagID == i)
  firstCap <- as.Date(min(indiv$collectDate), "%Y-%m-%d", tz = "UTC")
  lastCap <- as.Date(max(indiv$collectDate), "%Y-%m-%d", tz = "UTC")
  possibleDates <- seq(as.Date(firstCap), as.Date(lastCap), by="days")
  potentialNights <- mam_plotNight_nodups %>%
    filter(as.character(collectDate) %in%
             as.character(possibleDates) & plotID %in% indiv$plotID) %>%
    dplyr::select(nightuid,plotID, collectDate, eventID) %>%
    mutate(tagID=i)
  allnights <- left_join(potentialNights, indiv)
  known_alive <- bind_rows(known_alive, allnights)
}

# # Make fxn to calculate MNKA at each plot at each sampling bout
mnka_per_plot <- function(capture_data) {
  mnka_by_plot_bout <- capture_data %>% group_by(eventID,plotID) %>%
    summarize(n=n_distinct(tagID))
  return(mnka_by_plot_bout)
}

# call fxn on capsNew (which has records of all nights that each individuals was known to be alive)
MNKA.bp<-mnka_per_plot(capture_data = known_alive)
head(MNKA.bp)

write.csv(as.data.frame(MNKA.bp), file = "C:/Users/seano/Desktop/Projects/NEONSmallMammals/all-mnka-bp.csv")

# 32 captures did not get eventID, so must be removed
pejn_caps<-pejn_caps[!is.na(pejn_caps[,'eventID']),]

# clean up vars


write.csv(pejn_caps, file = "C:/Users/seano/Desktop/Projects/NEONSmallMammals/pejn_caps.csv")

