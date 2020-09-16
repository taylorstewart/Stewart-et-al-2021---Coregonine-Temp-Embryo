#################################################################
##
##
##
#################################################################

## CLEAR THE ENVIRONMENT FIRST ==================================

rm(list = ls(all.names = TRUE))


## LOAD PACKAGES ================================================

library(dplyr)           # data manipulation
library(magrittr)        # for %<>%
library(stringr)         # for logical pattern matching
library(data.table)      # for fread() - super fast load-in of data


## MAKE A LIST OF FILES TO LOAD =================================

length.files <- list.files('/Users/Taylor/CloudStation/Cisco-Climate-Change/Development-Pictures/LakeOntario-Cisco-2019',
                        full.names = TRUE, recursive = TRUE, pattern = '\\.txt$')
family.data <- read_excel("data/LakeOntario-Cisco-2019.xlsx", sheet = "FamilyPlateData") %>% 
  dplyr::select(cross = Cross, trt = Treatment, plate.no = PlateNo)

## READ IN BATCH TEXT FILES =====================================

larval.length.trt <- do.call(rbind, lapply(length.files, function(i) {
  print(i)
  
  raw <- fread(file = i) %>% 
    rename(well.id = 1, length.pix = 2) %>% 
    mutate(calib = ifelse(str_detect(well.id, "PL") == TRUE, "Calibration", "Length"),
           date = as.POSIXct(paste0(substr(i, 138, 141), "-", substr(i, 142, 143), "-", substr(i, 144, 145)), format = "%Y-%m-%d"),
           trt = as.numeric(substr(i, 150, 150)),
           plate = as.numeric(substr(i, 157, 159)),
           rep = as.numeric(substr(i, 164, 164))) %>% 
    dplyr::select(date, trt, plate, rep, well.id, length.pix, calib)
  
  calibrartion <- raw %>% filter(calib == "Calibration") %>% 
    group_by(date, trt, plate, rep) %>% 
    summarise(calib.pix = mean(length.pix))
  
  data <- raw %>% filter(calib == "Length") %>% 
    dplyr::select(-calib) %>% 
    left_join(calibrartion) %>% 
    mutate(row = substr(well.id, 1, 1),
           column = substr(well.id, 2, 2),
           plate.no = paste0(trt, "-", plate, row),
           well.no = paste0(trt, "-", plate, row, "-", column),
           length.mm = length.pix/calib.pix) %>% 
    left_join(family.data) %>% 
    dplyr::select(date, cross, trt, plate, plate.no, well.no, length.pix, calib.pix, length.mm) %>% 
    mutate(trt = ifelse(trt == 1, 2.0, 4.5))
}))

write.csv(larval.length.trt, "data/LakeOntario-Cisco-2019-Lengths.csv", row.names = FALSE)
