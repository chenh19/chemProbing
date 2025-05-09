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
  
  filename <- gsub("./3.analysis/5.parse/","",csv)
  filename <- gsub("_", "-", filename)
  filename <- gsub(".mpileup.csv","",filename)
  filename <- paste0(filename, "_", region)
  
  pdf(paste0("./3.analysis/6.plot/absolute/", filename, ".pdf"), width = 15, height = 10)
  print(
    ggplot(df, aes(x=Position, y=Mut_percentage)) +
      geom_col(fill = "steelblue", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 5)) +
      labs(title = filename, x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}


# plot relative values
RM <- csvs[grepl("RM1", csvs)]
ctrl <- csvs[grepl("DMSO", csvs)]
ctrl <- lapply(ctrl, read.csv)
internal_ctrl <- add_columns(ctrl[[1]], ctrl[[2]], by=c("Region", "Position", "Ref_base"))
internal_ctrl <- add_columns(internal_ctrl, ctrl[[3]], by=c("Region", "Position", "Ref_base"))
internal_ctrl <- mutate(internal_ctrl, Mut_percentage_avg = round((Mut_percentage + Mut_percentage_1 + Mut_percentage_2)/3, 4))
internal_ctrl <- internal_ctrl[c("Region", "Position", "Ref_base", "Mut_percentage_avg")]
internal_ctrl$Mut_percentage_avg[internal_ctrl$Mut_percentage_avg < 0] <- 0
summary <- c("avg_A","avg_C","avg_G","avg_U")

for (rm in RM) {
  
  df <- read.csv(rm, header = TRUE)
  df$Mut_percentage[df$Mut_percentage < 0] <- 0
  region <- unique(df$Region)
  
  difference <- add_columns(df, internal_ctrl, by=c("Region", "Position", "Ref_base"))
  difference <- mutate(difference, Difference = Mut_percentage - Mut_percentage_avg)
  
  a <- filter(difference, Ref_base == "A" | Ref_base == "a")
  a <- round(mean(a$Difference),4)
  c <- filter(difference, Ref_base == "C" | Ref_base == "c")
  c <- round(mean(c$Difference),4)
  g <- filter(difference, Ref_base == "G" | Ref_base == "g")
  g <- round(mean(g$Difference),4)
  t <- filter(difference, Ref_base == "T" | Ref_base == "t")
  t <- round(mean(t$Difference),4)
  summary <- data.frame(summary,c(a,c,g,t))
  
  filename <- gsub("./3.analysis/5.parse/","",rm)
  filename <- gsub("_001.mpileup.csv","",filename)
  filename <- gsub("_", "-", filename)
  filename <- paste0(filename, "_", region)
  
  pdf(paste0("./3.analysis/6.plot/relative/", filename, ".pdf"), width = 15, height = 3)
  print(
    ggplot(difference, aes(x=Position, y=Difference)) +
      geom_col(fill = "firebrick", width = 0.8) +
      scale_x_continuous(breaks = seq(0, max(df$Position), by = 50)) +
      scale_y_continuous(limits = c(-5.5, 9.5)) +
      labs(title = paste0(filename, " (relative to DMSO-average)"), x = "Position", y = "Mutation (%)") +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))
  )
  dev.off()
}
names <- gsub(".mpileup.csv", "", basename(RM))
colnames(summary) <- c("Base", names)
summary_long <- pivot_longer(summary, cols = 2:5, names_to = "Group", values_to = "Value")
summary_long$Group <- factor(summary_long$Group, levels = c("Base",
                                                            "RM12-1_S4_L001_001",
                                                            "RM12-2_S1_L001_001",
                                                            "RM15-1_S8_L001_001",
                                                            "RM15-2_S45_L001_001"
                                                            
))
colors <- c("avg_A" = "#1F77B4", 
            "avg_C" = "#FF7F0E", 
            "avg_G" = "#2CA02C",
            "avg_U" = "#D62728")
pdf(paste0("./3.analysis/6.plot/relative/base_summary.pdf"), width = 5.5, height = 6)
print(
  ggplot(summary_long, aes(x = Group, y = Value, fill = Base)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
    scale_fill_manual(values = colors) +
    coord_cartesian(ylim = c(-1, 3)) +
    theme_minimal() +
    labs(title = "Summary", x = NULL, y = "Mutation Rate %") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 65, hjust = 1),
      axis.title.x = element_text(margin = margin(t = 10))
    )
)
dev.off()


# cleanup
rm(list = ls())
