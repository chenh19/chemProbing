library(tidyverse)
library(ggplot2)
library(expss)

# create folder
if (dir.exists("./3.analysis/8.plot/")==FALSE){
  dir.create("./3.analysis/8.plot/")
}
if (dir.exists("./3.analysis/8.plot/absolute")==FALSE){
  dir.create("./3.analysis/8.plot/absolute")
}
if (dir.exists("./3.analysis/8.plot/relative")==FALSE){
  dir.create("./3.analysis/8.plot/relative")
}


# plot absolute values
csvs <- list.files(path='./3.analysis/7.parse', pattern='*.mpileup.csv', full.names = TRUE)
for (csv in csvs) {

  df <- read.csv(csv, header = TRUE)
  df <- filter(df, !is.na(Mut_percentage))
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  filename <- gsub("./3.analysis/7.parse/WV2JYX_","",csv)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  filename <- gsub("-cut-trimmed", "", filename)
  filename <- gsub("-short-trimmed", "", filename)
  
  p <- ggplot(df, aes(x=Position, y=Mut_percentage)) +
      geom_col(fill = "steelblue", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(0, 15)) +
      labs(title = filename, x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  # ggsave(paste0("./3.analysis/8.plot/absolute/", filename, ".pdf"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
  ggsave(paste0("./3.analysis/8.plot/absolute/", filename, ".jpg"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
}


## plot relative values (SL5 37 degrees)
SL5_37Cs <- csvs[grepl("37C_15min", csvs)]
internal_ctrl <- read.csv("./3.analysis/7.parse/WV2JYX_2_02_SL5_DMSO_37C_15min_Mn_cut_trimmed.mpileup.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed")
internal_ctrl$Mut_percentage[internal_ctrl$Mut_percentage < 0] <- 0

summary <- c("avg_A","avg_C","avg_G","avg_U")

for (SL5_37C in SL5_37Cs) {
  
  df <- read.csv(SL5_37C, header = TRUE)
  df <- filter(df, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed")
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  difference <- add_columns(df, internal_ctrl, by=c("Region", "Position", "Ref_base"))
  colnames(difference) <- c("Region", "Position", "Ref_base", "Mut_percentage", "Mut_percentage_ctrl")
  difference <- mutate(difference, Difference = Mut_percentage - Mut_percentage_ctrl)
  
  a <- filter(difference, Ref_base == "A" | Ref_base == "a")
  a <- round(mean(a$Difference),4)
  c <- filter(difference, Ref_base == "C" | Ref_base == "c")
  c <- round(mean(c$Difference),4)
  g <- filter(difference, Ref_base == "G" | Ref_base == "g")
  g <- round(mean(g$Difference),4)
  t <- filter(difference, Ref_base == "T" | Ref_base == "t")
  t <- round(mean(t$Difference),4)
  summary <- data.frame(summary,c(a,c,g,t))
  
  filename <- gsub("./3.analysis/7.parse/WV2JYX_","",SL5_37C)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  filename <- gsub("-cut-trimmed", "", filename)
  filename <- gsub("-short-trimmed", "", filename)
  
  p <- ggplot(difference, aes(x=Position, y=Difference)) +
      geom_col(fill = "firebrick", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(-5, 10)) +
      labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  # ggsave(paste0("./3.analysis/8.plot/relative/", filename, ".pdf"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
  ggsave(paste0("./3.analysis/8.plot/relative/", filename, ".jpg"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
}
names <- gsub("WV2JYX_", "", basename(SL5_37Cs))
names <- gsub(".mpileup.csv", "", names)
colnames(summary) <- c("Base", names)
summary_long <- pivot_longer(summary, cols = 2:7, names_to = "Group", values_to = "Value")
summary_long$Group <- factor(summary_long$Group, levels = c("Base",
                                                            "1_01_SL5_DMSO_37C_15min_Mg_cut_trimmed",
                                                            "2_02_SL5_DMSO_37C_15min_Mn_cut_trimmed",
                                                            "3_03_SL5_C12_5mM_37C_15min_Mn_cut_trimmed",
                                                            "4_04_SL5_C15_5mM_37C_15min_Mn_cut_trimmed",
                                                            "5_05_SL5_C12_0.5mM_37C_15min_Mn_cut_trimmed",
                                                            "6_06_SL5_C15_0.5mM_37C_15min_Mn_cut_trimmed"
                                                            
))
summary_long$Group <- gsub('_cut_trimmed', '', summary_long$Group)
sum1 <- summary_long


## plot relative values (SL5 80 degrees)
SL5_80Cs <- csvs[grepl("80C_15min", csvs)]
internal_ctrl <- read.csv("./3.analysis/7.parse/WV2JYX_8_08_SL5_DMSO_80C_15min_Mn_cut_trimmed.mpileup.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed")
internal_ctrl$Mut_percentage[internal_ctrl$Mut_percentage < 0] <- 0

summary <- c("avg_A","avg_C","avg_G","avg_U")

for (SL5_80C in SL5_80Cs) {
  
  df <- read.csv(SL5_80C, header = TRUE)
  df <- filter(df, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed")
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  difference <- add_columns(df, internal_ctrl, by=c("Region", "Position", "Ref_base"))
  colnames(difference) <- c("Region", "Position", "Ref_base", "Mut_percentage", "Mut_percentage_ctrl")
  difference <- mutate(difference, Difference = Mut_percentage - Mut_percentage_ctrl)
  
  a <- filter(difference, Ref_base == "A" | Ref_base == "a")
  a <- round(mean(a$Difference),4)
  c <- filter(difference, Ref_base == "C" | Ref_base == "c")
  c <- round(mean(c$Difference),4)
  g <- filter(difference, Ref_base == "G" | Ref_base == "g")
  g <- round(mean(g$Difference),4)
  t <- filter(difference, Ref_base == "T" | Ref_base == "t")
  t <- round(mean(t$Difference),4)
  summary <- data.frame(summary,c(a,c,g,t))
  
  filename <- gsub("./3.analysis/7.parse/WV2JYX_","",SL5_80C)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  filename <- gsub("-cut-trimmed", "", filename)
  filename <- gsub("-short-trimmed", "", filename)
  
  p <- ggplot(difference, aes(x=Position, y=Difference)) +
      geom_col(fill = "firebrick", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(-5, 10)) +
      labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  # ggsave(paste0("./3.analysis/8.plot/relative/", filename, ".pdf"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
  ggsave(paste0("./3.analysis/8.plot/relative/", filename, ".jpg"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
}
names <- gsub("WV2JYX_", "", basename(SL5_80Cs))
names <- gsub(".mpileup.csv", "", names)
colnames(summary) <- c("Base", names)
summary_long <- pivot_longer(summary, cols = 2:7, names_to = "Group", values_to = "Value")
summary_long$Group <- factor(summary_long$Group, levels = c("Base",
                                                            "7_07_SL5_DMSO_80C_15min_Mg_cut_trimmed",
                                                            "8_08_SL5_DMSO_80C_15min_Mn_cut_trimmed",
                                                            "9_09_SL5_C12_5mM_80C_15min_Mn_cut_trimmed",
                                                            "10_10_SL5_C15_5mM_80C_15min_Mn_cut_trimmed",
                                                            "11_11_SL5_C12_0.5mM_80C_15min_Mn_cut_trimmed",
                                                            "12_12_SL5_C15_0.5mM_80C_15min_Mn_cut_trimmed"
))
summary_long$Group <- gsub('_cut_trimmed', '', summary_long$Group)
sum2 <- summary_long


## plot base summary
all <- rbind(sum1,sum2)
all$Group <- factor(all$Group, levels = c("Base",
                                          "1_01_SL5_DMSO_37C_15min_Mg",
                                          "2_02_SL5_DMSO_37C_15min_Mn",
                                          "3_03_SL5_C12_5mM_37C_15min_Mn",
                                          "4_04_SL5_C15_5mM_37C_15min_Mn",
                                          "5_05_SL5_C12_0.5mM_37C_15min_Mn",
                                          "6_06_SL5_C15_0.5mM_37C_15min_Mn",
                                          "7_07_SL5_DMSO_80C_15min_Mg",
                                          "8_08_SL5_DMSO_80C_15min_Mn",
                                          "9_09_SL5_C12_5mM_80C_15min_Mn",
                                          "10_10_SL5_C15_5mM_80C_15min_Mn",
                                          "11_11_SL5_C12_0.5mM_80C_15min_Mn",
                                          "12_12_SL5_C15_0.5mM_80C_15min_Mn"
                                                            
))
colors <- c("avg_A" = "#1F77B4", 
            "avg_C" = "#FF7F0E", 
            "avg_G" = "#2CA02C",
            "avg_U" = "#D62728")
p <- ggplot(all, aes(x = Group, y = Value, fill = Base)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  scale_fill_manual(values = colors) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  theme_minimal() +
  labs(title = "Summary (with reads trimming)", x = NULL, y = "Mutation Rate %") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10))
  )
# ggsave(paste0("./3.analysis/8.plot/relative/base_summary_37C.pdf"), plot = p, width = 8, height = 3, units = "in", dpi = 1200)
ggsave(paste0("./3.analysis/8.plot/relative/base_summary.jpg"), plot = p, width = 10, height = 5, units = "in", dpi = 1200)



# cleanup
rm(list = ls())
