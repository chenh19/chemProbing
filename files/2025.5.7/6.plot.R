library(tidyverse)
library(ggplot2)
library(expss)

# create folder
if (dir.exists("./3.analysis/6.plot/")==FALSE){
  dir.create("./3.analysis/6.plot/")
}
if (dir.exists("./3.analysis/6.plot/absolute")==FALSE){
  dir.create("./3.analysis/6.plot/absolute")
}
if (dir.exists("./3.analysis/6.plot/relative")==FALSE){
  dir.create("./3.analysis/6.plot/relative")
}


# plot absolute values
csvs <- list.files(path='./3.analysis/5.parse', pattern='*.mpileup.csv', full.names = TRUE)
for (csv in csvs) {

  df <- read.csv(csv, header = TRUE)
  df <- filter(df, !is.na(Mut_percentage))
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  filename <- gsub("./3.analysis/5.parse/WV2JYX_","",csv)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  pdf(paste0("./3.analysis/6.plot/absolute/", filename, ".pdf"), width = 15, height = 3)
  print(
    ggplot(df, aes(x=Position, y=Mut_percentage)) +
      geom_col(fill = "steelblue", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(0, 15)) +
      labs(title = filename, x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}


# plot relative values (SL5 37 degrees)
SL5_37Cs <- csvs[grepl("37C_15min", csvs)]
internal_ctrl <- read.csv("./3.analysis/5.parse/WV2JYX_2_02_SL5_DMSO_37C_15min_Mn.mpileup.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SARS-CoV2-5UTR-Amplicon")
internal_ctrl$Mut_percentage[internal_ctrl$Mut_percentage < 0] <- 0

for (SL5_37C in SL5_37Cs) {
  
  df <- read.csv(SL5_37C, header = TRUE)
  df <- filter(df, Region == "SARS-CoV2-5UTR-Amplicon")
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  difference <- add_columns(df, internal_ctrl, by=c("Region", "Position", "Ref_base"))
  colnames(difference) <- c("Region", "Position", "Ref_base", "Mut_percentage", "Mut_percentage_ctrl")
  difference <- mutate(difference, Difference = Mut_percentage - Mut_percentage_ctrl)
  
  filename <- gsub("./3.analysis/5.parse/WV2JYX_","",SL5_37C)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  pdf(paste0("./3.analysis/6.plot/relative/", filename, ".pdf"), width = 15, height = 3)
  print(
    ggplot(difference, aes(x=Position, y=Difference)) +
      geom_col(fill = "firebrick", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(-5, 10)) +
      labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}


# plot relative values (SL5 80 degrees)
SL5_80Cs <- csvs[grepl("80C_15min", csvs)]
internal_ctrl <- read.csv("./3.analysis/5.parse/WV2JYX_8_08_SL5_DMSO_80C_15min_Mn.mpileup.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SARS-CoV2-5UTR-Amplicon")
internal_ctrl$Mut_percentage[internal_ctrl$Mut_percentage < 0] <- 0

for (SL5_80C in SL5_80Cs) {
  
  df <- read.csv(SL5_80C, header = TRUE)
  df <- filter(df, Region == "SARS-CoV2-5UTR-Amplicon")
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  difference <- add_columns(df, internal_ctrl, by=c("Region", "Position", "Ref_base"))
  colnames(difference) <- c("Region", "Position", "Ref_base", "Mut_percentage", "Mut_percentage_ctrl")
  difference <- mutate(difference, Difference = Mut_percentage - Mut_percentage_ctrl)
  
  filename <- gsub("./3.analysis/5.parse/WV2JYX_","",SL5_80C)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  pdf(paste0("./3.analysis/6.plot/relative/", filename, ".pdf"), width = 15, height = 3)
  print(
    ggplot(difference, aes(x=Position, y=Difference)) +
      geom_col(fill = "firebrick", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(-5, 10)) +
      labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}


# plot relative values (SMN2 80 degrees)
SMN2_80Cs <- csvs[grepl("80C_5min", csvs)]
internal_ctrl <- read.csv("./3.analysis/5.parse/WV2JYX_13_13_Ctrl_DMSO_80C_5min_Mn.mpileup.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SMN2-Exon7-Amplicon")
internal_ctrl$Mut_percentage[internal_ctrl$Mut_percentage < 0] <- 0

for (SMN2_80C in SMN2_80Cs) {
  
  df <- read.csv(SMN2_80C, header = TRUE)
  df <- filter(df, Region == "SMN2-Exon7-Amplicon")
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  difference <- add_columns(df, internal_ctrl, by=c("Region", "Position", "Ref_base"))
  colnames(difference) <- c("Region", "Position", "Ref_base", "Mut_percentage", "Mut_percentage_ctrl")
  difference <- mutate(difference, Difference = Mut_percentage - Mut_percentage_ctrl)
  
  filename <- gsub("./3.analysis/5.parse/WV2JYX_","",SMN2_80C)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  pdf(paste0("./3.analysis/6.plot/relative/", filename, ".pdf"), width = 15, height = 3)
  print(
    ggplot(difference, aes(x=Position, y=Difference)) +
      geom_col(fill = "firebrick", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(-5, 10)) +
      labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}


# cleanup
rm(list = ls())
