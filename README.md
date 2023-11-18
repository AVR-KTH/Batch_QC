---
title: "Quality control"
author: "Andrea Villanueva Raisman"
output: html_document
---

```{r}
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(factoextra)
library(Rtsne)
library(umap)
library(pheatmap)
library(broom)
library(DEqMS)
library(ggvenn)
library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(modelr)
```
# Introduction

This is an exploratory analysis with the objective of checking for missing values, outliers or bias. Regarding missing values, R’s default way of handling them is to drop them. However, the command na.action = na.warn informs about these instances.  This exploratory analysis is performed on 17 plates. Each of these plates has 88 randomized plasma samples and 3 replicates of in-house plasma. The exploratory analysis begins with the quality control using these in-house plasma replicates. 

By looking at the normalized area of the peptide intensity peaks, we can assess if they follow a normal distribution, which they do. By looking at the boxplot, we can find outliers. 
```{r}
#Missing values warning
options(na.action = na.warn)

#Read in data
data <- read_csv2('Peptide Ratio Results.csv')
#Get rid of QTags
data <- subset(data, Protein != "QTag")
#Transform from character to numeric class
data$`Ratio To Standard` <- as.numeric(sub(",", ".", data$`Ratio To Standard`))
data$`Normalized Area` <- as.numeric(sub(",", ".", data$`Normalized Area`))
data$`Protein Abundance` <- as.numeric(sub(",", ".", data$`Protein Abundance`))
data$`DotProductLightToHeavy` <- as.numeric(sub(",", ".", data$`DotProductLightToHeavy`))

#Histogram of normalized area
ggplot(data, aes(x = log10(`Normalized Area`))) +
  geom_histogram(fill = 'grey50', color = 'grey30') +
  theme_bw()

# Boxplot
ggplot(data, aes(x = `Replicate Name`, y = log10(`Ratio To Standard`))) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

```
If we filter the peptides whose log10(`Ratio To Standard`) is above 5, we see that well A12 in sample 25 and B12 in sample 30

```{r}
data %>%
  filter(log10(`Ratio To Standard`) > 5)

filtered_data <- data %>%
  filter(`DotProductLightToHeavy` > 0.85)

#Histogram of normalized area
ggplot(filtered_data, aes(x = log10(`Normalized Area`))) +
  geom_histogram(fill = 'grey50', color = 'grey30') +
  theme_bw()

# Boxplot
ggplot(filtered_data, aes(x = `Replicate Name`, y = log10(`Ratio To Standard`))) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


```
# Peptide quantification inter-plate comparisons: APOA1 and VTCN

APOA1 and VTNC are used to compare between plates because they are very stable, proteins that are commonly measured in blood. The stability of APOA1 and VTNC makes them suitable candidates for comparisons between plates or as reference proteins in certain experimental settings. Being both stable and well-characterized makes them suitable references to ensure consistency and reliability in this experiment, especially when comparing results across plates.

```{r}
#Separation of plate (batch) and wells
data2 <- separate(data, `Replicate Name`, into = c("Plate", "Well"), sep = "_")

#APOA1
apoa1_data <- subset(data2, Protein == "sp|P02647|APOA1_HUMAN")

#Visualize ratio to standard
ggplot(apoa1_data, aes(x = `Plate`, y = `Ratio To Standard` )) +
  geom_point(aes(color = `Peptide`)) +
  labs(title = "APOA1 Ratio to Standard", x = "Plate", y = "Ratio to standard")

#VTNC
vtnc_data <- subset(data2, Protein == "sp|P04004|VTNC_HUMAN")

#Visualize ratio to standard
ggplot(vtnc_data, aes(x = `Plate`, y = `Ratio To Standard` )) +
  geom_point(aes(color = `Peptide`)) +
  labs(title = "VTNC Ratio to Standard", x = "Plate", y = "Ratio to standard")
```
```{r}
# Extract last three characters from 'Well'
data2$Well <- substr(data2$Well, nchar(data2$Well) - 2, nchar(data2$Well))

# APOA1
apoa1_pep <- subset(data2, Protein == "sp|P02647|APOA1_HUMAN") %>%
  subset(Peptide == "THLAPYSDELR")
  
# Visualize ratio to standard for APOA1
ggplot(apoa1_pep, aes(x = Plate, y = `Ratio To Standard`, color = Well, group = Well)) +
  geom_point() +
  geom_line() +
  labs(title = "APOA1-THLAPYSDELR Ratio to Standard", x = "Plate", y = "Ratio to standard") +
  theme_minimal()

# VTNC
vtnc_pep <- subset(data2, Protein == "sp|P04004|VTNC_HUMAN") %>%
  subset(Peptide == "DVWGIEGPIDAAFTR")

# Visualize ratio to standard for VTNC
ggplot(vtnc_pep, aes(x = Plate, y = `Ratio To Standard`, color = Well, group = Well)) +
  geom_point() +
  geom_line() +
  labs(title = "VTNC-DVWGIEGPIDAAFTR Ratio to Standard", x = "Plate", y = "Ratio to standard") +
  theme_minimal()

```
# Correlation between APOA1 and VTNC
To see if both peptides follow the same pattern, we can make a plot and see the correlation of their peptides's median or the correlation of the peptide with the smallest inter-plate CV.
```{r}
# Calculate median 'Ratio To Standard' for each 'Peptide' within 'Protein'
apoa1_median <- apoa1_data %>%
  group_by(interaction(Plate, Well)) %>%
  summarize(Median_Ratio = median(`Ratio To Standard`))

vtnc_median <- vtnc_data %>%
  group_by(interaction(Plate, Well)) %>%
  summarize(Median_Ratio = median(`Ratio To Standard`))

correlation_data_1 <- merge(apoa1_median, vtnc_median, by = "interaction(Plate, Well)")

# Perform linear regression to get residuals
lm_model <- lm(Median_Ratio.x ~ Median_Ratio.y, data = correlation_data_1)
residuals <- abs(residuals(lm_model))

# Find indices of samples with largest residuals
outliers <- order(residuals, decreasing = TRUE)[1:5] 

# Create a logical vector to indicate outliers
correlation_data_1$outlier <- row.names(correlation_data_1) %in% row.names(correlation_data_1[outliers, ])

# Plot correlation between median ratios with labels for top outliers
ggplot(correlation_data_1, aes(x = Median_Ratio.x, y = Median_Ratio.y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = subset(correlation_data_1, outlier), aes(label = `interaction(Plate, Well)`), hjust = -0.2, vjust = -0.5) +
  labs(title = "Correlation of Median Ratios", x = "APOA1 Median Ratio", y = "VTNC Median Ratio") +
  theme_minimal()

```


```{r}
# Calculate CV for each 'Peptide' within 'Protein'
apoa1_cv <- apoa1_data %>%
  group_by(Peptide) %>%
  summarize(Interplate_CV = sd(`Ratio To Standard`) / mean(`Ratio To Standard`) * 100)

vtnc_cv <- vtnc_data %>%
  group_by(Peptide) %>%
  summarize(Interplate_CV = sd(`Ratio To Standard`) / mean(`Ratio To Standard`) * 100)

# Find the smallest inter-plate CV for each protein
apoa1_min_cv <- apoa1_cv %>%
  filter(Interplate_CV == min(Interplate_CV)) %>%
  pull(Interplate_CV)

vtnc_min_cv <- vtnc_cv %>%
  filter(Interplate_CV == min(Interplate_CV)) %>%
  pull(Interplate_CV)

# Extract the peptide with the smallest CV
min_cv_peptides <- bind_rows(
  apoa1_cv %>% filter(Interplate_CV == apoa1_min_cv),
  vtnc_cv %>% filter(Interplate_CV == vtnc_min_cv)
)

# Calculate median 'Ratio To Standard' for the smallest CV peptides
apoa1_smallest_cv <- apoa1_data %>%
  filter(Peptide %in% min_cv_peptides$Peptide) %>%
  group_by(interaction(Plate, Well)) %>%
  summarize(Ratio = mean(`Ratio To Standard`))

vtnc_smallest_cv <- vtnc_data %>%
  filter(Peptide %in% min_cv_peptides$Peptide) %>%
  group_by(interaction(Plate, Well)) %>%
  summarize(Ratio = mean(`Ratio To Standard`))

correlation_data <- merge(apoa1_smallest_cv, vtnc_smallest_cv, by = "interaction(Plate, Well)")

# Perform linear regression to get residuals
lm_model <- lm(Ratio.x ~ Ratio.y, data = correlation_data)
residuals <- abs(residuals(lm_model))

# Find indices of samples with largest residuals
outliers <- order(residuals, decreasing = TRUE)[1:5] 

# Create a logical vector to indicate outliers
correlation_data$outlier <- row.names(correlation_data) %in% row.names(correlation_data[outliers, ])

# Plot correlation between median ratios with labels for top outliers
ggplot(correlation_data, aes(x = Ratio.x, y = Ratio.y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = subset(correlation_data, outlier), aes(label = `interaction(Plate, Well)`), hjust = -0.2, vjust = -0.5) +
  labs(title = "Correlation of ratios for peptide with smallest CV", x = "APOA1-THLAPYSDELR  ratio", y = "VTNC-DVWGIEGPIDAAFTR  ratio") +
  theme_minimal()
```
