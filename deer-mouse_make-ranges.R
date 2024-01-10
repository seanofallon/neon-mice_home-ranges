# This script uses NSF-NEON small mammal box trapping capture data to calculate home ranges
# for Peromyscus maniculatus (PEMA) individuals, and assign traits and conditions to those individuals
# for the analysis performed in 'Uncovering multiple influences on space use by deer mice using
# NEON data' by S. O'Fallon, N. Pinter-Wollman, and K. Mabry. For details on NEON data collection
# and data products, refer to their website: https://data.neonscience.org/

# INPUTS:
# pejn_caps.csv - NEON-recorded data for each PEMA capture w/ tagID in their 'small mammal box
#                   trapping' protocol. Used to calculate home ranges and assign traits/conditions.
# all-mnka-bp.csv - Minimum number (of PEMA) known alive, representing population density, at each 
#               small mammal sampling plot at each sampling event.

# OUTPUTS:
# ao_ranges - dataframe containing home range area calculations and trait/condition assignments
#             for all eligible PEMA individuals analyzed in this study. Used as input for statistical
#             analysis and data visualization in 'deer-mouse_stats-figs.R'


##### Start up: clear env, load packages, read in data
# clear env
rm(list=ls())

# load packages
library(dplyr)
library(neonUtilities)
library(neonOS)
library(raster)
library(adehabitatHR)
library(stringr)
library(ggpubr)

# read in pejn_caps and mnka-bp (if these are in your working directory, run as is. Otherwise,
# change file name to the full path name to where you've saved pejn_caps.csv and mnka-bp.csv)
pejn_caps<-read.csv("pejn_caps.csv")
mnka.bp<-read.csv("all-mnka-bp.csv")

# define functions needed later:
# Function that converts letters to numbers, to deal w/ NEON coordinates:
LETTER2num <- function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}

### Reduce pejn_caps df to peromyscus captured at least 5 times
# find tagIDs w/ at least 5 caps
recap_ids <- pejn_caps %>% count(tagID) %>% filter(n >= 5)
# filter capture df to only tags captured >4x
pe_recaps <- filter(pejn_caps, tagID %in% recap_ids$tagID)


### Set outlier mice to be removed - found in later steps to have data entry errors
# remove unrealistic weights caused by data entry/measurement errors - replace w/ NA
pe_recaps["weight"][pe_recaps["weight"] > 50] <- NA
# remove unrealistic hindfoot lengths caused by data entry/measurement errors - replace w/ NA
pe_recaps["hindfootLength"][pe_recaps["hindfootLength"] > 28] <- NA

# remove 'NEON.MAM.D01.R1761', who had errors in pregnancyStatus
pe_recaps["weight"][pe_recaps["tagID"] == 'NEON.MAM.D01.R1761'] <- NA

# remove "NEON.MAM.D08.R1921", who had an unrealistically low final capture weight
pe_recaps["weight"][pe_recaps["tagID"] == "NEON.MAM.D08.R1921"] <- NA

# pull out pregnant females to plot body condition by weight
nopreg_recaps <- pe_recaps[pe_recaps$pregnancyStatus != 'pregnant',]


# add column w/ size category based on TaxonID
pe_recaps$size <- with(pe_recaps, ifelse(taxonID %in% c('PEMA','PELE','PELEPEMA'),
                                         'small', 'big'))
nopreg_recaps$size <- with(nopreg_recaps, ifelse(taxonID %in% c('PEMA','PELE','PELEPEMA'),
                                                 'small', 'big'))

## boxplot of big vs. small for all caps
bigVsmall<-boxplot(weight ~ size, data = pe_recaps)

t.test(weight ~ size, data = pe_recaps)
t.test(weight ~ size, data = nopreg_recaps)


### Add col w/ lifestage based on weight to restrict dataset to adult mice
for (i in 1:nrow(pe_recaps)) {
  if (pe_recaps[i,'size'] == 'small') {
    if (is.na(pe_recaps[i,'weight']) == T) {
      pe_recaps[i,'lifeStageByWeight']<-'unknown'
    }else if (pe_recaps[i,'weight']<16) {
      pe_recaps[i,'lifeStageByWeight']<-'subadult'
    }else if (pe_recaps[i,'weight']>16) {
      pe_recaps[i,'lifeStageByWeight']<-'adult'
    } else
      pe_recaps[i,'lifeStageByWeight']<-'unknown'
  }
  else if (pe_recaps[i,'size'] == 'big') {
    if (is.na(pe_recaps[i,'weight']) == T) {
      pe_recaps[i,'lifeStageByWeight']<-'unknown'
    }else if (pe_recaps[i,'weight']<19) {
      pe_recaps[i,'lifeStageByWeight']<-'subadult'
    }else if (pe_recaps[i,'weight']>19) {
      pe_recaps[i,'lifeStageByWeight']<-'adult'
    } else
      pe_recaps[i,'lifeStageByWeight']<-'unknown'
  }
}

# loop goes through all indivs in recap_ids, replacing all lifeStageByWeight
# after first 'adult' w/ 'adult' (mice stay adults after growing to adulthood)
for (i in 1:nrow(recap_ids)) {
  imouse<-filter(pe_recaps,tagID==recap_ids$tagID[i]) # take all caps for ith mouse
  if ('adult' %in% imouse$lifeStageByWeight) { # if we ever assign 'adult' to imouse
    # assign date of first 'adult' assignment to first.adult
    first.adult<-min(imouse[imouse$lifeStageByWeight=='adult','collectDate.x'])
    # put 'adult' for lifeStageByWeight in imouse for all caps after 1st 'adult' assignment
    imouse[imouse$collectDate.x >= first.adult,'lifeStageByWeight'] <- 'adult'
    # apply to pe_recaps
    pe_recaps[pe_recaps$tagID==recap_ids$tagID[i],'lifeStageByWeight']<-
      imouse$lifeStageByWeight
  }
}

### Filter pe_recaps to mice captured at least 5x as an adult - our dataset
# remove non-adult captures
ao_recaps<-filter(pe_recaps,lifeStageByWeight=='adult')
# find tagIDs w/ at least 5 adult caps
ao_ids <- ao_recaps %>% count(tagID) %>% filter(n >= 5)
# filter capture df to only tags captured >4x
ao_recaps <- filter(ao_recaps, tagID %in% ao_ids$tagID)


##### Calculate home ranges
# set up df to hold home ranges
ao_ranges <- data.frame(tagID=character(),
                        UD=numeric())

# make objects for HR loop
unq_id=unique(ao_recaps$tagID) # all IDs - to dictate # of iterations

# get HR of all individuals:
for(i in 1:length(unq_id)){
  ix=unq_id[i]
  ao_ranges[i,'tagID']<-ix
  traps=ao_recaps$trapCoordinate[ao_recaps$tagID==ix] # take the column with the trap coordinates of a single individual from the NEON data
  splt= data.table::tstrsplit(traps, split="") # take the NEON coordinates and split into a letter element (1) and a number element (2)
  x_coor=10*unlist(lapply(splt[[1]],FUN=LETTER2num)) # turn A,B,C etc to numbers and multiply by 10 to get units in meters
  if(max(x_coor)>200){ # some grid locations were marked as X so need to get rid of those - missing data
    x_coor[x_coor==240]=NA  
  }
  # the tstrsplit function above creates 2 elements for single-digit numbered grid cells (like A1) and 3 elements for 2-digit numbered grid cells (like A10) - the only one two digit one is 10 so we deal with it separately 
  if(length(splt)==2){ # if there are no traps at row 10
    y_coor=10*as.numeric(unlist(splt[[2]])) # multiplied by 10 to get units in meters
  }
  if(length(splt)==3){# if there are traps at row 10
    y_coor=10*as.numeric(unlist(splt[[2]])) # multiplied by 10 to get units in meters
    y_coor[!is.na(splt[[3]])]=100
  }
  
  coor=cbind(x_coor,y_coor) # bind xy into a single variable
  coor=coor[!is.na(coor[,1]),] #remove NAs - if X was entered instead of a letter, this will remove those entries
  coor=coor[!is.na(coor[,2]),] #remove NAs - if X was entered instead of a number, this will remove those entries
  coor = as.data.frame(coor)
  for_samp = dim(coor)[1]
  # get Home range using KDE:
  if(dim(unique(coor))[1]>2 & nrow(coor)>4){ # need more than 4 points to calculate a HR and at least 3 need to be unique locations
    coor = sp::SpatialPoints(coor) # convert coordinates to the format that adehabitatHR likes
    # get area based on kernel
    ud_est= kernelUD(coor) # make the UD of the animal
    ao_ranges[i,'UD']=as.numeric(kernel.area(ud_est, percent = 50, unin = "m", unout = "m2")) # get the area of the 50% UD of the animal
  }else{
    ao_ranges[i,'UD'] = NA
  }
}


# remove NAs from ao_ranges
ao_ranges <- na.omit(ao_ranges)


### Attach vars to ao_ranges for modeling
for (i in 1:nrow(ao_ranges)) { # for each range calculated:
  # take adult caps of ith mouse
  imouse <- filter(ao_recaps,tagID==ao_ranges$tagID[i])
  ### SPECIES
  ao_ranges[i,'taxonID']<-imouse$taxonID[1]
  ### SIZE CATEGORY
  ao_ranges[i,'size']<-imouse$size[1]
  ### SEX
  # if mouse was ever pregnant, assign it F in ao_ranges
  if ('pregnant' %in% imouse$pregnancyStatus) {
    ao_ranges[i,'sex']<-'F'}
  # if mouse assigned M > F, assign it M in ao_ranges
  else if ('M' > 'F' %in% imouse$sex) {
    ao_ranges[i,'sex']<-'M'}
  # if mouse assigned F > M, assign it F in ao_ranges
  else if ('M' < 'F' %in% imouse$sex) {
    ao_ranges[i,'sex']<-'F'}
  # that should cover all cases - but if it doesn't, assign U for unknown
  else {
    ao_ranges[i,'sex']<-'U'}
  ### SITE
  ao_ranges[i,'site'] <- imouse[1,"siteID"]
  ### PLOT ID (needed for MNKA assignment...)
  ao_ranges[i,'plotID'] <- imouse[1,"plotID"]
  ### VEG TYPE
  # if imouse caps occured in forests/woodlands, assign to vegType 'forest'
  if (imouse$nlcdClass[1]=='deciduousForest' | imouse$nlcdClass[1]=='mixedForest'
      | imouse$nlcdClass[1]=='evergreenForest' | imouse$nlcdClass[1]=='woodyWetlands') {
    ao_ranges[i,'vegType']<-'forest'}
  # if imouse caps occured in grassy habitat, assign to vegType 'grassland'
  else if (imouse$nlcdClass[1]=='cultivatedCrops' | imouse$nlcdClass[1]=='grasslandHerbaceous'
           | imouse$nlcdClass[1]=='pastureHay') {
    ao_ranges[i,'vegType']<-'grassland'}
  # if imouse caps occured in shrubby habitat, assign to vegType 'shrubland'
  else if (imouse$nlcdClass[1]=='shrubScrub') {
    ao_ranges[i,'vegType']<-'shrubland'}
  # that should cover everything, but in case of errors, assign U for unknown
  else  {
    ao_ranges[i,'vegType']<-'U'
  }
  ### HINDFOOT LENGTH
  ao_ranges[i,'meanHindfootLength']<-mean(imouse[imouse$pregnancyStatus!='pregnant' & imouse$pregnancyStatus!='unknown',]$hindfootLength,na.rm=T)  
  ### WEIGHT
  ao_ranges[i,'meanWeight']<-mean(imouse[imouse$pregnancyStatus!='pregnant' & imouse$pregnancyStatus!='unknown',]$weight,na.rm=T)
  ### RANGE OF DATES COLLECTED
  ao_ranges[i,'firstCollected']<-min(imouse$collectDate.x)
  ao_ranges[i,'lastCollected']<-max(imouse$collectDate.x)
  ### YEAR
  ao_ranges[i,'year']<-str_trunc(ao_ranges[i,'lastCollected'],4,side='right',ellipsis = '')
  ### LATITUDE
  ao_ranges[i,'latitude'] <- imouse[1,'decimalLatitude']
  ### MIN NUMBER KNOWN ALIVE (PEROMYSCUS DENSITY)
  mouse_alive<-pejn_caps[pejn_caps$collectDate.x >= ao_ranges$firstCollected[i] &
                           pejn_caps$collectDate.x <= ao_ranges$lastCollected[i] &
                           pejn_caps$site == ao_ranges$site[i] &
                           pejn_caps$plotID == ao_ranges$plotID[i],]
  
  all_mnka4mouse<-mnka.bp[mnka.bp$eventID %in% unique(mouse_alive$eventID) &
                            mnka.bp$plotID == unique(mouse_alive$plotID),'n']
  ao_ranges[i,'meanMNKA']<-mean(all_mnka4mouse,na.rm=T)
  ao_ranges[i,'maxMNKA']<-max(all_mnka4mouse,na.rm=T)
}


### Calculate bodyCondition based on residuals of relationship btwn meanHindfootLength and meanWeight
# first, eliminate NAs from ao_ranges as they cannot be included in lm
ao_ranges<-na.omit(ao_ranges)

# check - is there a different relationship btwn weight and foot length in big vs. small species?
sVb_test <- lm(meanWeight ~ meanHindfootLength * size, ao_ranges)
summary(sVb_test)
# yes there is - need to calculate separately

# calculate relationship for big species
mw.by.mhl.big <- lm(meanWeight ~ meanHindfootLength, ao_ranges[ao_ranges$size=='big',])
# calculate relationship for small species
mw.by.mhl.small <- lm(meanWeight ~ meanHindfootLength, ao_ranges[ao_ranges$size=='small',])

### Make Figure S2 w/ saved plots
par(mfrow = c(1, 3))

boxplot(weight ~ size, data = pe_recaps,
        ylab = 'Weight', xlab = 'Size')

plot(meanWeight ~ meanHindfootLength, ao_ranges[ao_ranges$size=='big',],
     ylab = 'Mean Weight', xlab = 'Mean Hindfoot Length')
abline(mw.by.mhl.big)

plot(meanWeight ~ meanHindfootLength, ao_ranges[ao_ranges$size=='small',],
     ylab = 'Mean Weight', xlab = 'Mean Hindfoot Length')
abline(mw.by.mhl.small)


# residuals
mw.by.mhl.small$residuals
mw.by.mhl.big$residuals

aob_ranges <- ao_ranges[ao_ranges$size=='big',]
nrow(aob_ranges)

for (i in 1:nrow(aob_ranges)) {
  aob_ranges[i,'bodyCondition'] <- mw.by.mhl.big$residuals[i]
}

aos_ranges <- ao_ranges[ao_ranges$size=='small',]
nrow(aos_ranges)

for (i in 1:nrow(aos_ranges)) {
  aos_ranges[i,'bodyCondition'] <- mw.by.mhl.small$residuals[i]
}

ao_ranges<-rbind(aos_ranges, aob_ranges)

# it can be helpful to save ao_ranges so that it can be loaded and used in the analysis script
save(ao_ranges, file = "all_ao_ranges.RData") # might need to add specification for directory in which to save in the filename
