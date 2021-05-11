# Uniformity Plate Analysis
# All wells contain the total binding (TOTB) control which generates the maximum raw data signal, but represents vehicle (0% biological activity).

#Set up -------------------------------------------------
library(tidyverse)
library(lubridate)
library(viridis)
library(ggpubr)
library(rstatix)

#rescale_mean -------------------------------------------
rescale_mean <- function(x) {
  x / mean(x)
}

# rescale_med ------------------------------------------------
rescale_med <- function(x) {
  x / median(x)
}

# outlier_robust Calculate outliers using box plot criteria----
outlier_robust <- function(x) {
  x > quantile(x, 0.75) + (1.5 * IQR(x)) | x < quantile(x, 0.25) - (1.5 * IQR(x))
}

# Prepare data -------------------

  UniformData <- read_csv('Data/AllData.csv')%>%
  filter(PlateMap == 'Uniform') %>%
  select(-PlateMap) %>%
  mutate(Assay = as_factor(Assay),
         PlateType = as.character(PlateType),
         Well = as_factor(Well),
         PlateId = as.character(PlateId)) %>%
  group_by(PlateId,ExpTime, Assay, PlateType) %>%
  nest() %>%
  mutate(ExpTime = mdy_hms(ExpTime),
         ExpDate = as_date(ExpTime)) %>%
  unnest(cols= c(data)) %>%
  ungroup() %>%
  group_by(PlateId) %>%
  mutate(MedData = rescale_med(Data)) %>%
  ungroup() %>%
  pivot_longer(cols = contains('Data'),
               names_to = 'Scale',
               values_to = 'Data') %>%
  mutate(Scale = if_else(Scale == 'Data', 'Raw', 'Median'),
         Scale = as_factor(Scale),
         Scale = fct_rev(Scale)
  ) %>%
  ungroup()

ggplot(UniformData, aes(x = PlateId, y = Data)) +
  geom_boxplot(aes(color = PlateType)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, alpha = 0.5) +
  labs(x = "Plate",
       y = 'Data') +
  scale_fill_viridis(discrete=TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  facet_grid(Scale ~ Assay, scales = 'free') +
  labs(title = 'TOTB Plate Uniformity Data')
