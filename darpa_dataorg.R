# Data organization

# Script that organizes count tables from feature counts into proper order, pairs the with a col data file and saves 
# them to easily identifiable names. 

library("glue")
library("readr")
library("writexl")
library("tidyr")
library("dplyr")
library("tibble")

# ------------- A549 (flu) high MOI, miRNA -------------------

AH_cd <- read.csv("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/miRNA/AH_cd.csv", stringsAsFactors = FALSE)
AH_cd$Time <- as.factor(AH_cd$Time)

AH_cd$X <- as.character(AH_cd$X)
AH_cd[5,1]="6_Mock1_1"
AH_cd[6,1]="6_Mock1_2"
AH_cd[15,1]="12_Mock1_1"
AH_cd[16,1]="12_Mock1_2"
AH_cd[25,1]="24_Mock1_1"
AH_cd[26,1]="24_Mock1_2"

AH_cd <- data.frame(AH_cd, row.names = 1)

AH_miRNA_df <- read.csv("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/miRNA/AH_miRNAs_full_new.csv", 
                        stringsAsFactors = FALSE)
AH_miRNA_df <- AH_miRNA_df[-c(659),]
rownames(AH_miRNA_df) <- AH_miRNA_df$X.miRNA
AH_miRNA_df <- AH_miRNA_df[,-c(1,65:124)]
AH_miRNA_df <- AH_miRNA_df[,-c(1:3)]



AH_miRNA_df <- as.data.frame(sapply(AH_miRNA_df, as.numeric), row.names = row.names(AH_miRNA_df))

# create a temporary variable to store the row names
col_names <- row.names(AH_cd)
# set the column names of the count data to match the row names of the colData
names(AH_miRNA_df) <- col_names
# clean up
rm(col_names)

cd <- AH_cd
df <- AH_miRNA_df
rm(AH_cd)
rm(AH_miRNA_df)

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_countdata/AH_mirna_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/AH_mirna_coldata.csv")

file_prefix <- "AH_miRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$Time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(Time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)
rm(cd_time)

# ------------- A549 high MOI, mRNA -------------------


# set working directory to the one containing the column data
# read in the colData
# change the Time column to a factor
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/old_col_data")
AHm_cd <- read.csv("A549_high_cd_full.csv", row.names=1)
AH_cd$Time <- as.factor(AHm_cd$Time)
AHm_cd <- rbind(AHm_cd, "6_H1N1_2"=c("H1N1",6,2))
AHm_cd <- AHm_cd[c(1:5,36,6:35),]
AHm_cd <- AHm_cd[c(1,2,13,14,25,26,5,6,17,18,29,30,7,8,19,20,31,32,3,4,15,16,27,28,9,10,21,22,33,34,11,12,23,24,35,36),]

# read in the full count table, set the gene IDs to be the row names and make strings
# be actual strings and not factors
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/full_tables")
AHm_df <- read.table("AH_Flu.count.txt", row.names = 1, stringsAsFactors = FALSE)

AHm_df <- AHm_df[,-c(1:5)] # get rid of the first 5 columns
AHm_df <- AHm_df[,c(5,6,1:4,11,12,7:10,17,18,13:16,35,36,31:34,23,24,19:22,29,30,25:28)] # reorder to the same order as the df file

AHm_df <- AHm_df[-c(1),]

# convert the counts to numbers
AHm_df <- as.data.frame(sapply(AHm_df, as.numeric), row.names = row.names(AHm_df))

# create a temporary variable to store the row names
col_names <- row.names(AHm_cd)
# set the column names of the count data to match the row names of the colData
names(AHm_df) <- col_names
# clean up
rm(col_names)

cd <- AHm_cd
df <- AHm_df

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_countdata/AH_Flu_mRNA_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/AH_Flu_mRNA_coldata.csv")

file_prefix <- "AH_Flu_mRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$Time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(Time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)
rm(AHm_df)
rm(AHm_cd)

# ------------- HBEpC (Flu) low MOI, mRNA -------------------

# set working directory to the one containing the column data
# read in the colData
# change the Time column to a factor
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/old_col_data")
cd <- read.csv("HBEpC_low_cd.csv", row.names=1)
cd$Time <- as.factor(cd$Time)

# read in the full count table, set the gene IDs to be the row names and make strings
# be actual strings and not factors
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/full_tables")
df <- read.table("HBEpC_low_Flu_counts.txt", row.names = 1, stringsAsFactors = FALSE)

df <- df[-c(1), -c(1:5)] # get rid of the first row and the first 5 columns

# convert the counts to numbers
df <- as.data.frame(sapply(df, as.numeric), row.names = row.names(df))

# create a temporary variable to store the row names
col_names <- row.names(cd)
# set the column names of the count data to match the row names of the colData
names(df) <- col_names
# clean up
rm(col_names)

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_countdata/HBEL_Flu_mRNA_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/HBEL_Flu_mRNA_coldata.csv")

file_prefix <- "HBEL_Flu_mRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$Time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(Time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)

# ------------- A549 (flu) low MOI, mRNA -------------------

# set working directory to the one containing the column data
# read in the colData
# change the Time column to a factor
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/old_coldata")
cd <- read.csv("A549_high_cd_full.csv", row.names=1)
cd$Time <- as.factor(cd$Time)
cd <- rbind(cd, "6_H1N1_2"=c("H1N1",6,2))
cd <- cd[c(1:5,36,6:35),]
cd <- cd[c(1,2,13,14,25,26,5,6,17,18,29,30,7,8,19,20,31,32,3,4,15,16,27,28,9,10,21,22,33,34,11,12,23,24,35,36),]

# read in the full count table, set the gene IDs to be the row names and make strings
# be actual strings and not factors
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/full_tables")
df <- read.table("AH_Flu.count.txt", row.names = 1, stringsAsFactors = FALSE)

df <- df[,-c(1:5)] # get rid of the first 5 columns
df <- df[,c(5,6,1:4,11,12,7:10,17,18,13:16,35,36,31:34,23,24,19:22,29,30,25:28)] # reorder to the same order as the df file

df <- df[-c(1),]

# convert the counts to numbers
df <- as.data.frame(sapply(df, as.numeric), row.names = row.names(df))

# create a temporary variable to store the row names
col_names <- row.names(cd)
# set the column names of the count data to match the row names of the colData
names(df) <- col_names
# clean up
rm(col_names)

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_countdata/AL_Flu_mRNA_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/AL_Flu_mRNA_coldata.csv")

file_prefix <- "AL_Flu_mRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$Time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(Time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)

# ------------- A549 (zika) high MOI, mRNA -------------------

# set working directory to the one containing the column data
# read in the colData
# change the Time column to a factor
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/old_col_data")
cd <- read.csv("A549_high_ZV_cd.csv", row.names=1)
cd$Time <- as.factor(cd$Time)
cd <- cd[c(3,4,9,10,15,16,5,6,11,12,17,18,21,22,27,28,33,34,23,24,29,30,35,36,1,2,19,20,7,8,25,26,13,14,31,32),]

# read in the full count table, set the gene IDs to be the row names and make strings
# be actual strings and not factors
setwd("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/full_tables")
df <- read.table("AH_Zika.count.txt", row.names = 1, stringsAsFactors = FALSE)

df <- df[,-c(1:5)] # get rid of the first 5 columns
df <- df[,c(35,36,31:34,23,24,19:22,29,30,25:28,5,6,1:4,11,12,17,18,7:10,13:16)] # reorder to the same order as the df file

df <- df[-c(1),]

# convert the counts to numbers
df <- as.data.frame(sapply(df, as.numeric), row.names = row.names(df))

# create a temporary variable to store the row names
col_names <- row.names(cd)
# set the column names of the count data to match the row names of the colData
names(df) <- col_names
# clean up
rm(col_names)

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_countdata/AH_zika_mRNA_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/AH_zika_mRNA_coldata.csv")

file_prefix <- "AH_zika_mRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$Time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(Time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Flu_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)

rm(df)
rm(cd)
rm(AH_cd)
rm(AH_miRNA_df)
rm(cd_time)

#------------------------- HBEpC (Flu) high MOI, miRNA ------------------------------

setwd("~/Dropbox/Ghedin_Lab/DARPA_seq/count_Tables")
df <- read.csv("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_Tables/HBEpC_mature_miRNAexpressed.csv", 
                        stringsAsFactors = FALSE)
df <- df[,-c(35:64)]
df <- df[-c(659),]
rownames(df) <- df$X.miRNA
df <- df[,-c(1:4)]

# convert the counts to numbers
df <- as.data.frame(sapply(df, as.numeric), row.names = row.names(df))

cd <- read.csv("./HBEpC_miRNA_cd.csv")
rownames(cd) <- glue("{cd$virus}_{cd$time}-{cd$replicate}")

# create a temporary variable to store the row names
col_names <- row.names(cd)
# set the column names of the count data to match the row names of the colData
names(df) <- col_names
# clean up
rm(col_names)

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_seq/count_Tables/processed_countdata/HBEH_miRNA_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_seq/count_Tables/processed_coldata/HBEH_miRNA_coldata.csv")

file_prefix <- "HBEH_miRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)

rm(df)
rm(cd)
rm(cd_time)

#------------------------- HFF1 (Zika) high MOI, miRNA ------------------------------

setwd("~/Dropbox/Ghedin_Lab/DARPA_seq/count_Tables")
df <- read.csv("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_Tables/HFF1_mature_miRNAsexpressed.csv", 
               stringsAsFactors = FALSE)
df <- df[,-c(35:64)]
df <- df[-c(659),]
rownames(df) <- df$X.miRNA
df <- df[,-c(1:4)]

# convert the counts to numbers
df <- as.data.frame(sapply(df, as.numeric), row.names = row.names(df))

cd <- read.csv("./HFF1_miRNA_cd.csv")
rownames(cd) <- glue("{cd$virus}_{cd$time}-{cd$replicate}")

# create a temporary variable to store the row names
col_names <- row.names(cd)
# set the column names of the count data to match the row names of the colData
names(df) <- col_names
# clean up
rm(col_names)

df <- df %>% rownames_to_column("rownames")
cd <- cd %>% rownames_to_column("rownames")

write_excel_csv(df, path = "~/Dropbox/Ghedin_Lab/DARPA_seq/count_Tables/processed_countdata/HFF1H_miRNA_counts.csv")
write_excel_csv(cd, path = "~/Dropbox/Ghedin_Lab/DARPA_seq/count_Tables/processed_coldata/HFF1H_miRNA_coldata.csv")

file_prefix <- "HFF1H_miRNA"

# Create a vector of the Time column of the column data
timepoints <- as.vector(cd$time)
# Find the unique time points within this column
timepoints <- (unique(timepoints))

# make seperate cd files for each time point
for (i in seq_along(timepoints))
{
  # filter the column data for rows that contain that time point
  # place this new data frame into a list
  cd_time <- cd %>%
    filter(time==timepoints[i])
  fn <- glue("{file_prefix}_{timepoints[i]}_cd.csv")
  write_excel_csv(cd_time, path = glue("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_Tables/processed_coldata/{fn}"))
}
rm(i)

rm(df)
rm(cd)
rm(cd_time)
rm(file_prefix)
rm(fn)
rm(timepoints)

#------------------------- A549 low MOI, miRNA ------------------------------

setwd("~/Dropbox/Ghedin_Lab/DARPA_Seq")
df = read.csv(file = "./count_Tables/full_tables/AL_miRNA_expressed_all.csv", 
                            stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
df <- counts_table[,1:88]
rownames(df) <- df[,1]


cd <- read.csv("./count_Tables/full_tables/AL_miRNA_coldata/csv")
rownames(cd) <- glue("{cd$virus}_{cd$time}-{cd$replicate}")

# create a temporary variable to store the row names
col_names <- row.names(cd)
# set the column names of the count data to match the row names of the colData
names(df) <- col_names
# clean up
rm(col_names)