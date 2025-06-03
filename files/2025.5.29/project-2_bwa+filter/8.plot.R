library(tidyverse)
library(ggplot2)
library(expss)

# create folder
if (dir.exists("./3.analysis/8.plot/")==FALSE){
  dir.create("./3.analysis/8.plot/")
}
if (dir.exists("./3.analysis/8.plot/3.absolute_mut")==FALSE){
  dir.create("./3.analysis/8.plot/3.absolute_mut")
}
if (dir.exists("./3.analysis/8.plot/4.relative_mut")==FALSE){
  dir.create("./3.analysis/8.plot/4.relative_mut")
}
if (dir.exists("./3.analysis/8.plot/4.relative_mut/summary")==FALSE){
  dir.create("./3.analysis/8.plot/4.relative_mut/summary")
}
if (dir.exists("./3.analysis/8.plot/5.base_summary")==FALSE){
  dir.create("./3.analysis/8.plot/5.base_summary")
}

################################################################################

# plot absolute values
csvs <- list.files(path='./3.analysis/7.parse', pattern='*.mpileup.csv', full.names = TRUE)
for (csv in csvs) {
  
  df <- read.csv(csv, header = TRUE)
  df <- filter(df, !is.na(Mut_percentage))
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  filename <- gsub("./3.analysis/7.parse/","",csv)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  p <- ggplot(df, aes(x=Position, y=Mut_percentage)) +
    geom_col(fill = "steelblue", width = 0.8) +
    scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
    coord_cartesian(ylim = c(0, 10)) +
    labs(title = filename, x = "Position", y = "Mutation (%)") +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  #ggsave(paste0("./3.analysis/8.plot/3.absolute_mut/", filename, ".pdf"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
  ggsave(paste0("./3.analysis/8.plot/3.absolute_mut/", filename, ".jpg"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
}

################################################################################

# calculate mean mut percentage
csvs <- list.files(path='./3.analysis/7.parse', pattern='*.mpileup.csv', full.names = TRUE)
list <- read.csv("./list.csv", header = FALSE)
colnames(list) <- c("Sample", "Group")
groups <- unique(list$Group)

for (group in groups) {
  
  #group <- groups[[1]]
  
  replicate <- filter(list, Group == group)
  replicate <- replicate$Sample
  
  merge <- add_columns(read.csv(paste0("./3.analysis/7.parse/", replicate[[1]], ".mpileup.csv"), header = TRUE), 
                       read.csv(paste0("./3.analysis/7.parse/", replicate[[2]], ".mpileup.csv"), header = TRUE), 
                       by=c("Region", "Position", "Ref_base"))
  merge <- add_columns(merge, 
                       read.csv(paste0("./3.analysis/7.parse/", replicate[[3]], ".mpileup.csv"), header = TRUE), 
                       by=c("Region", "Position", "Ref_base"))
  merge <- mutate(merge,
                  Mean_mut_percentage = round((Mut_percentage + Mut_percentage_1 + Mut_percentage_2) / 3, 4))
  
  merge <- merge[c("Region", "Position", "Ref_base", "Mean_mut_percentage")]
  merge$Mean_mut_percentage[merge$Mean_mut_percentage < 0] <- 0
  
  write.table(merge, file=paste0("./3.analysis/8.plot/4.relative_mut/summary/", group, ".mean.csv"), sep=",", row.names = FALSE)
  
}

################################################################################

# Project-2 - part-1

csvs <- list.files(path='./3.analysis/8.plot/4.relative_mut/summary', pattern='*.mean.csv', full.names = TRUE)
SL5_30mins <- csvs[grepl("37C-30min", csvs)]
internal_ctrl <- read.csv("./3.analysis/8.plot/4.relative_mut/summary/18_rna_DMSO_37C-30min_Mn.mean.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed+3")
internal_ctrl$Mean_mut_percentage[internal_ctrl$Mean_mut_percentage < 0] <- 0

summary <- c("avg_A","avg_C","avg_G","avg_U")

for (SL5_30min in SL5_30mins) {
  
  df <- read.csv(SL5_30min, header = TRUE)
  df <- filter(df, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed+3")
  df$Mean_mut_percentage[df$df$Mean_mut_percentage < 0] <- 0
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
  
  filename <- gsub("./3.analysis/8.plot/4.relative_mut/summary/","",SL5_30min)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mean.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  p <- ggplot(difference, aes(x=Position, y=Difference, fill=Ref_base)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("a" = "#1F77B4", "c" = "#FF7F0E", "g" = "#2CA02C", "t" = "#D62728")) +
    scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
    scale_y_continuous(limits = c(-2.5, 5)) +
    labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  #ggsave(paste0("./3.analysis/8.plot/4.relative_mut/", filename, ".pdf"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
  ggsave(paste0("./3.analysis/8.plot/4.relative_mut/", filename, ".jpg"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
}
names <- gsub(".mean.csv", "", basename(SL5_30mins))
colnames(summary) <- c("Base", names)
summary_long <- pivot_longer(summary, cols = 2:5, names_to = "Group", values_to = "Value")
summary_long$Group <- factor(summary_long$Group, levels = c("Base",
                                                            "17_rna_DMSO_37C-30min_Mg",
                                                            "18_rna_DMSO_37C-30min_Mn",
                                                            "23_rna_C34-Squarate-1mM_37C-30min_Mn",
                                                            "24_rna_Squarate-only-1mM_37C-30min_Mn"
                                                            
))
sum1 <- summary_long

################################################################################

# Project-2 - part-2

csvs <- list.files(path='./3.analysis/8.plot/4.relative_mut/summary', pattern='*.mean.csv', full.names = TRUE)
SL5_2hs <- csvs[grepl("37C-2h", csvs)]
internal_ctrl <- read.csv("./3.analysis/8.plot/4.relative_mut/summary/20_rna_DMSO_37C-2h_Mn.mean.csv", 
                          header = TRUE)
internal_ctrl <- filter(internal_ctrl, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed+3")
internal_ctrl$Mean_mut_percentage[internal_ctrl$Mean_mut_percentage < 0] <- 0

summary <- c("avg_A","avg_C","avg_G","avg_U")

for (SL5_2h in SL5_2hs) {
  
  df <- read.csv(SL5_2h, header = TRUE)
  df <- filter(df, Region == "SARS-CoV2-5UTR-Amplicon-short-trimmed+3")
  df$Mean_mut_percentage[df$df$Mean_mut_percentage < 0] <- 0
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
  
  filename <- gsub("./3.analysis/8.plot/4.relative_mut/summary/","",SL5_2h)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mean.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  p <- ggplot(difference, aes(x=Position, y=Difference, fill=Ref_base)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("a" = "#1F77B4", "c" = "#FF7F0E", "g" = "#2CA02C", "t" = "#D62728")) +
    scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
    scale_y_continuous(limits = c(-2.5, 5)) +
    labs(title = paste0(filename, " (relative to DMSO-Mn-Ctrl)"), x = "Position", y = "Mutation (%)") +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  #ggsave(paste0("./3.analysis/8.plot/4.relative_mut/", filename, ".pdf"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
  ggsave(paste0("./3.analysis/8.plot/4.relative_mut/", filename, ".jpg"), plot = p, width = 15, height = 3, units = "in", dpi = 1200)
}
names <- gsub(".mean.csv", "", basename(SL5_2hs))
colnames(summary) <- c("Base", names)
summary_long <- pivot_longer(summary, cols = 2:5, names_to = "Group", values_to = "Value")
summary_long$Group <- factor(summary_long$Group, levels = c("Base",
                                                            "19_rna_DMSO_37C-2h_Mg",
                                                            "20_rna_DMSO_37C-2h_Mn",
                                                            "21_rna_C34-NM-1mM_37C-2h_Mn",
                                                            "22_rna_NM-only-1mM_37C-2h_Mn"
                                                            
))
sum2 <- summary_long

################################################################################

## plot base summary
all <- rbind(sum1,sum2)
all$Group <- factor(all$Group, levels = c("Base",
                                          "19_rna_DMSO_37C-2h_Mg",
                                          "20_rna_DMSO_37C-2h_Mn",
                                          "21_rna_C34-NM-1mM_37C-2h_Mn",
                                          "22_rna_NM-only-1mM_37C-2h_Mn",
                                          "17_rna_DMSO_37C-30min_Mg",
                                          "18_rna_DMSO_37C-30min_Mn",
                                          "23_rna_C34-Squarate-1mM_37C-30min_Mn",
                                          "24_rna_Squarate-only-1mM_37C-30min_Mn"
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
  labs(title = "Summary", x = NULL, y = "Mutation Rate %") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10))
  )
# ggsave(paste0("./3.analysis/8.plot/5.base_summary/base_summary.pdf"), plot = p, width = 8, height = 5, units = "in", dpi = 1200)
ggsave(paste0("./3.analysis/8.plot/5.base_summary/base_summary.jpg"), plot = p, width = 8, height = 5, units = "in", dpi = 1200)

################################################################################

# cleanup
rm(list = ls())
