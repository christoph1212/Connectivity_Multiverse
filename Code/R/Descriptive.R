## Descriptive Analysis
# This script is part of the study "Robustness of EEG functional brain networks 
# associated with fluid intelligence: A multiverse analysis of connectivity and 
# thresholding methods". It performs descriptive analyses
#
# Written by: Christoph Fruehlinger
# Last edit: July 2026

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

## Load relevant packages for further analysis
if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, apaTables, tidyverse, effsize, ggsignif, ggpubr, psychometric, 
               stringr, ggstance, ggh4x)

options(scipen=999)
rm(list = ls())
set.seed(42)

data <- read.csv("Data/Connectivity/Graph_data.csv", header = TRUE)
main_path_conn   <- "imcoh"
main_path_thresh <- "dens"
data <- data %>%
  mutate(SetCond = paste(Run, Condition, sep = "_"))
# Select main path columns
pattern        <- paste('[a-z0-9]+', main_path_conn, '[a-z0-9]+', main_path_thresh, sep = '_')
main_path_cols <- grep(pattern, colnames(data), value = TRUE)
main_path_data <- data[, c('ID', 'SetCond', main_path_cols)]

subs <- main_path_data %>% 
  filter(SetCond == "first_eyes_closed") %>% 
  pull("ID") %>% 
  unique() %>%
  sort()

## Descriptive Analysis

# Add Sociodemographic Data

socio_table = read.csv("Data/SocioDemographics.txt", header = T, na.strings = c("", "NA"))

## Convert columns to factors and rename levels
socio_table$Gender <- as.factor(socio_table$Gender)
socio_table$HormonalContraceptives <- as.factor(socio_table$HormonalContraceptives)
socio_table$Ethnicity <- as.factor(socio_table$Ethnicity)
socio_table$MaritalStatus <- as.factor(socio_table$MaritalStatus)
socio_table$Occupancy <- as.factor(socio_table$Occupancy)
socio_table$HighestDegree <- as.factor(socio_table$HighestDegree)
socio_table$monthlyBruttoIncome <- as.factor(socio_table$monthlyBruttoIncome)

levels(socio_table$Gender)[levels(socio_table$Gender) == 1] <- "Female"
levels(socio_table$Gender)[levels(socio_table$Gender) == 2] <- "Male"
levels(socio_table$HormonalContraceptives)[levels(socio_table$HormonalContraceptives) == 1] <- "Yes"
levels(socio_table$HormonalContraceptives)[levels(socio_table$HormonalContraceptives) == 2] <- "No"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 1] <- "European"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 2] <- "Arabic"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 3] <- "African"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 4] <- "Asian"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 5] <- "Hisp"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 6] <- "other"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 1] <- "Single"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 2] <- "Relationship"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 3] <- "legal Partnership"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 4] <- "Divorced"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 5] <- "Widowed"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 1] <- "Working"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 2] <- "Student"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 3] <- "unemployed"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 1] <- "Hauptschule"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 2] <- "Realschule"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 3] <- "Abitur"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 4] <- "Hochschulabschluss"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 5] <- "Promotion"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 1] <- "<500"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 2] <- "500-1000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 3] <- "1000-2000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 4] <- "2000-3000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 5] <- ">3000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 6] <- "no info"

# Filter included participants
socio_fullSample <- socio_table %>% 
  filter(ID %in% subs)

# Calculate descriptive statistics of the sample demographics
# absolute Values
table(socio_fullSample$HighestDegree, useNA = "always")
table(socio_fullSample$Occupancy, useNA = "always")
table(socio_fullSample$Ethnicity, useNA = "always")
table(socio_fullSample$Gender, useNA = "always")

# relative Values
prop.table(table(socio_fullSample$HighestDegree, useNA = "always"))
prop.table(table(socio_fullSample$Occupancy, useNA = "always"))
prop.table(table(socio_fullSample$Ethnicity, useNA = "always"))
prop.table(table(socio_fullSample$Gender, useNA = "always"))

# Calculate McDonald's Omega
CC_omega <- omega(main_path_data %>%
                     dplyr::select(contains("cc")))
Pathl_omega <- omega(main_path_data %>%
                    dplyr::select(contains("pathl")))
Eglob_omega <- omega(main_path_data %>%
                       dplyr::select(contains("eglob")))
Eloc_omega <- omega(main_path_data %>%
                       dplyr::select(contains("eloc")))
SW_omega <- omega(main_path_data %>%
                       dplyr::select(contains("smallworld")))


message(paste0("\n", "McDonald's Omega (across Frequency Bands):", 
           "\n", "CC:    ", round(CC_omega$omega.tot, 3),
           "\n", "Pathl: ", round(Pathl_omega$omega.tot, 3),
           "\n", "Eglob: ", round(Eglob_omega$omega.tot, 3), 
           "\n", "Eloc:  ", round(Eloc_omega$omega.tot, 3),
           "\n", "SW:    ", round(SW_omega$omega.tot, 3)))


IST_full <- read.csv("Data/IST_table.csv", sep = ";", header = T)

IST_full$fluid <- rowSums(IST_full[,2:21])

IST_filtered <- IST_full %>%
  filter(ID %in% subs)

gf_omega <- omega(IST_filtered[,2:21])

message(paste0("\n", "Fluid Intelligence McDonald's Omega:", "\n",
           "gf: ", round(gf_omega$omega.tot, 3)))
