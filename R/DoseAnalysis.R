# Dose Response Data Analysis
# All wells contain the total binding (TOTB) control which generates the maximum raw data signal, but represents vehicle (0% biological activity).

#Set up -------------------------------------------------
library(tidyverse)
library(lubridate)
library(viridis)
library(ggpubr)
library(rstatix)

# outlier_robust Calculate outliers using box plot criteria----
outlier_robust <- function(x) {
  x > quantile(x, 0.75) + (1.5 * IQR(x)) | x < quantile(x, 0.25) - (1.5 * IQR(x))
}

# Prepare data -------------------
DRPlateMap <- read_csv('Data/Platemap.csv')

PlateId2  <- data.frame(PlateId2 = 1:28)

DoseData <- read_csv('Data/AllData.csv')%>%
  filter(PlateMap == 'Dose') %>%
  select(-PlateMap) %>%
  mutate(Assay = as_factor(Assay),
         PlateType = as.character(PlateType),
         Well = as_factor(Well),
         PlateId = as.character(PlateId)) %>%
  left_join(DRPlateMap) %>%
  mutate(Sample = fct_cross(Cmpd, as_factor(Conc))) %>%
  group_by(PlateId,ExpTime, Assay, PlateType) %>%
  nest() %>%
  mutate(ExpTime = mdy_hms(ExpTime),
         ExpDate = as_date(ExpTime),
         Run = as_factor(if_else(ExpDate < mdy('12/1/2019'), 'Run1', 'Run2'))) %>%
  arrange(Assay, Run) %>%
  bind_cols(PlateId2)%>%
  mutate(AssayPlate = if_else(Assay == 'Tgt1', PlateId2 - 14, PlateId2 - 0),
         AssayPlate = as_factor(AssayPlate))%>%
  select(-PlateId2) %>%
  unnest(cols= c(data)) %>%
  ungroup()

#  Normalize to Plate controls -------------------
CtrlData <- DoseData %>%
  filter(Cmpd == 'TOTB' | Cmpd == 'NSB') %>%
  group_by(Assay, PlateId, Cmpd) %>%
  summarise(Mean = mean(Data),
            Median = median(Data)) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with('M'), names_to = 'NormStatType', values_to = 'StatVal') %>%
  mutate(NormStat = paste(Cmpd, NormStatType, sep = '_')) %>%
  select(-Cmpd, -NormStatType) %>%
  pivot_wider(names_from = starts_with('Norm'), values_from = StatVal)

DoseData <- DoseData %>%
  left_join(CtrlData) %>%
  mutate(PctActMean = 100 * (1 - ((NSB_Mean - Data) / (NSB_Mean - TOTB_Mean))),
         PctActMedian = 100 * (1 - ((NSB_Median - Data) / (NSB_Median - TOTB_Median)))
  ) %>%
  select(-starts_with('NSB'), -starts_with('TOTB'))

# Plate Controls and QC ---------------------------
PlateQCData <- DoseData %>%
  filter(Cmpd == 'TOTB' | Cmpd == 'NSB') %>%
  select(-Conc, -Sample) %>%
  pivot_longer(cols = starts_with('Pct'), names_to = 'Scale', values_to = 'Activity') %>%
  mutate(Scale = if_else(str_detect(Scale, 'Mean'), 'Mean', 'Median')
         )

ggplot(PlateQCData, aes(x = AssayPlate, y = Activity, fill = Cmpd)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, alpha = 0.5) +
  labs(x = "Plate",
       y = 'Data') +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  facet_grid(Scale ~ Assay, scales = 'free') +
  labs(title = 'Dose Response Plate Controls')

PlateQCData <- PlateQCData %>%
  group_by(Assay, Run, AssayPlate, Cmpd, Scale) %>%
  summarise(Avg = mean(Activity),
            Med = median(Activity),
            StDev = sd(Activity),
            MAD = mad(Activity)) %>%
  ungroup() %>%
  mutate(ZLim = if_else(Cmpd == 'NSB', Avg - 3 * StDev, Avg + 3 * StDev),
         ZLimRob = if_else(Cmpd == 'NSB', Med - 3 * MAD, Med + 3 * MAD)) %>%
  pivot_longer(cols = 6:11, names_to = 'StatType', values_to = 'StatVal') %>%
  pivot_wider(names_from = Cmpd | StatType, names_sep = '_', values_from = StatVal) %>%
  mutate(Z = (NSB_ZLim - TOTB_ZLim) / (NSB_Avg - TOTB_Avg),
         Zrob = (NSB_ZLimRob - TOTB_ZLimRob) / (NSB_Med - TOTB_Med)
  )

ZComp <- PlateQCData %>%
  select(Assay, Scale, AssayPlate, Z, Zrob)%>%
  filter(Scale =='Mean') %>%
  pivot_longer(cols = 4:5, names_to = 'ZFactor', values_to = 'Zval')

ZCompTgt1 <- ZComp %>%
  filter(Assay == 'Tgt1') %>%
  t_test(Zval ~ ZFactor, paired = TRUE) %>%
  add_significance() %>%
  add_xy_position(x = 'ZFactor')

ggpaired(filter(ZComp, Assay == 'Tgt1'), x = 'ZFactor', y = 'Zval',
         id = 'AssayPlate',
         line.color = 'AssayPlate',
         order = c('Z', 'Zrob'),
         ylab = 'Z Value', xlab = 'Z Type') +
  stat_pvalue_manual(ZCompTgt1, tip.length = 0) +
  labs(title = 'Tgt1 Z Factor Comparison',
       subtitle = get_test_label(ZCompTgt1, detailed= TRUE)) +
  theme(legend.position = 'right')

ZCompTgt2 <- ZComp %>%
  filter(Assay == 'Tgt2') %>%
  t_test(Zval ~ ZFactor, paired = TRUE) %>%
  add_significance() %>%
  add_xy_position(x = 'ZFactor')

ggpaired(filter(ZComp, Assay == 'Tgt2'), x = 'ZFactor', y = 'Zval',
         id = 'AssayPlate',
         order = c('Z', 'Zrob'),
         line.color = 'AssayPlate',
         ylab = 'Z Value', xlab = 'Z Type') +
  stat_pvalue_manual(ZCompTgt2, tip.length = 0) +
  labs(title = 'Tgt2 Z Factor Comparison',
       subtitle = get_test_label(ZCompTgt2, detailed= TRUE)) +
  theme(legend.position = 'right')

# Cmpd Data - remove all wells without test compounds and put into tidy format ----------

CmpdData <- DoseData %>%
  filter(Conc != 0) %>%
  select(-Data, -Row, -Column) %>%
  group_by(Assay, AssayPlate) %>%
  pivot_longer(cols = starts_with('Pct'), names_to = 'Scale', values_to = 'Activity') %>%
  mutate(Scale = if_else(str_detect(Scale, 'Mean'), 'Mean', 'Median')) %>%
  ungroup()

PlateSmplData <- CmpdData %>%
  group_by(Assay, Run, AssayPlate, Sample, Scale) %>%
  mutate(Activity = if_else(Scale == 'Mean', mean(Activity), median(Activity))) %>%
  ungroup()

SampleSummData <- PlateSmplData %>%
  group_by(Assay, Sample, Scale) %>%
  mutate(Activity = if_else(Scale == 'Mean', mean(Activity), median(Activity))) %>%
  ungroup()

# Activite Well determination ------------------------

# Plate Activity Limits are based on the within plate TOTB controls correspond to the Mean/Median PctAct + 3*(SD/Mad) consistent with the normalization mode (Scale). This represents a real time assessment of Active wells on the plate independent of sample replication information and could be used to determine activity for either individual wells or aggregated sample replicates on a sample plate.

PlateActLims <- PlateQCData %>%
  select(Assay, AssayPlate, Scale, starts_with('TOTB_Z')) %>%
  mutate(ActLim = if_else(Scale == 'Mean', TOTB_ZLim, TOTB_ZLimRob)) %>%
  select(-starts_with('TOTB'))

# Phil stopping here. This is older code that shows different approaches, but needs to be updated to the current R objects above. The ggplot templates might be useful for some figures. Once I focused on one-way ANOVA Sample~Activity the direct calculation of the residuals was quicker for all the entire data set and I didn't have to subset the data by the Scale variable, but you might be better with purrr map functions than I am.
# # Model Residuals================
# #  All data 64 wells each sample -------
#
# CmpdOutlierRate <- CmpdData %>%
#   group_by(Assay) %>%
#   summarise(Wells = n(), Outliers = sum(AssayOutlier), OutlierRate = Outliers/Wells) %>%
#   ungroup()
#
# #All Data MSDs (64 Observatuibs/Sample) --------------------
# SampleResiduals <- CmpdData %>%
#   group_by(Assay,Sample) %>%
#   mutate(TrueAct = if_else(Scale == 'Mean', mean(Activity), median(Activity)),
#          Residual = Activity - TrueAct) %>%
#   ungroup()
#
# ggplot(SampleResiduals, aes(x= Residual)) +
#   geom_histogram() +
#   geom_vline(xintercept = quantile(SampleResiduals$Residual, 0.95), linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = quantile(SampleResiduals$Residual, 0.05), linetype = 'dashed', color = 'blue') +
#   labs(title = 'Sample Residuals by Assay and Scale',
#        x = 'Residuals') +
#   theme_minimal() +
#   facet_grid(Scale ~ Assay)
#
# SampleMSD <- SampleResiduals %>%
#   group_by(Assay, Scale) %>%
#   summarise(SDResiduals = sd(Residual),
#             MSD = 2 * sqrt(2) * SDResiduals,
#             N = n()) %>%
#   ungroup()%>%
#   mutate(Summarization = 'All')
#
# # Samples summarized by plate and 16 plates/sample  -------------
# SamplePlateResiduals <- CmpdData %>%
#   group_by(Assay, Scale, ReadId, Sample) %>%
#   summarise(MeanAct = mean(Activity), MedAct = median(Activity)) %>%
#   ungroup() %>%
#   mutate(Activity = if_else(Scale == 'Mean', MeanAct, MedAct)) %>%
#   select(-starts_with('M')) %>%
#   group_by(Assay,Sample, Scale) %>%
#   mutate(TrueAct = if_else(Scale == 'Mean', mean(Activity), median(Activity)),
#          Residual = Activity - TrueAct) %>%
#   ungroup()
#
# ggplot(SamplePlateResiduals, aes(x= Residual)) +
#   geom_histogram() +
#   geom_vline(xintercept = quantile(SamplePlateResiduals$Residual, 0.95), linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = quantile(SamplePlateResiduals$Residual, 0.05), linetype = 'dashed', color = 'blue') +
#   labs(title = 'Sample:Plate Residuals by Assay and Scale',
#        x = 'Residuals') +
#   theme_minimal() +
#   facet_grid(Scale ~ Assay)
#
# SamplePlateMSD <- SamplePlateResiduals %>%
#   group_by(Assay, Scale) %>%
#   summarise(SDResiduals = sd(Residual),
#             MSD = 2 * sqrt(2) * SDResiduals,
#             N = n()) %>%
#   ungroup() %>%
#   mutate(Summarization = 'Plate')
#
# # Samples summarized by Run and Plate Order and 2 plates/sample  -------------
# SamplePlateRunResiduals <- CmpdData %>%
#   group_by(Assay, Scale, Run, RunSample) %>%
#   summarise(MeanAct = mean(Activity), MedAct = median(Activity)) %>%
#   ungroup() %>%
#   mutate(Activity = if_else(Scale == 'Mean', MeanAct, MedAct)) %>%
#   select(-starts_with('M')) %>%
#   group_by(Assay, RunSample, Scale) %>%
#   mutate(TrueAct = if_else(Scale == 'Mean', mean(Activity), median(Activity)),
#          Residual = Activity - TrueAct) %>%
#   ungroup()
#
# ggplot(SamplePlateRunResiduals, aes(x= Residual)) +
#   geom_histogram() +
#   geom_vline(xintercept = quantile(SamplePlateRunResiduals$Residual, 0.95), linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = quantile(SamplePlateRunResiduals$Residual, 0.05), linetype = 'dashed', color = 'blue') +
#   labs(title = 'Sample:Plate:Run Residuals by Assay and Scale',
#        x = 'Residuals') +
#   theme_minimal() +
#   facet_grid(Scale ~ Assay)
#
# SamplePlateRunMSD <- SamplePlateResiduals %>%
#   group_by(Assay, Scale) %>%
#   summarise(SDResiduals = sd(Residual),
#             MSD = 2 * sqrt(2) * SDResiduals,
#             N = n()) %>%
#   ungroup() %>%
#   mutate(Summarization = 'Plate and Order')
#
# # MSD Summary Table -------------------
#
# SummMSDTable <- SampleMSD %>%
#   bind_rows(SamplePlateMSD) %>%
#   bind_rows(SamplePlateRunMSD)
#
# #ANOVAs ---------------------
# #All Data MSDs (64 Observatuibs/Sample) --------------------
# # Tgt1 ============================================
# # 2-way ANOVA -------------------
# Tgt1All.aov<- aov(Activity ~ Sample * Scale, data = filter(CmpdData, Assay == 'Tgt1'))
#
# Tgt1AllMSD <- 2 * sqrt(2) * sd(residuals(object = Tgt1All.aov))
#
# plot(Tgt1All.aov, 2)
#
# A <- quantile(Tgt1All.aov$residuals, 0.025)
# B <- quantile(TTgt1All.aov$residuals, 0.975)
#
# ggplot(Tgt1All.aov, aes(x= Tgt1All.aov$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt1 All Data',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Tgt1All.aov))
#
# cat('MSD = ', signif(Tgt1AllMSD, 3))
#
# # 1-way ANOVA (Mean) =========================
# Tgt1AllMean.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt1', Scale == 'Mean'))
#
#
# Tgt1AllMeanMSD <- 2 * sqrt(2) * sd(residuals(object = Tgt1AllMean.aov))
#
# plot(Tgt1AllMean.aov, 2)
#
# A <- quantile(Tgt1AllMean.aov$residuals, 0.025)
# B <- quantile(Tgt1AllMean.aov$residuals, 0.975)
#
# ggplot(Tgt1AllMean.aov, aes(x= Tgt1AllMean.aov$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt1 All Data (Means)',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Tgt1AllMean.aov))
#
# cat('MSD = ', signif(Tgt1AllMeanMSD, 3))
#
# # 1-way ANOVA Median ============
# Tgt1AllMed.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt1', Scale == 'Median'))
#
# Tgt1AllMedMSD <- 2 * sqrt(2) * sd(residuals(object = Tgt1AllMed))
#
# plot(Tgt1AllMed, 2)
#
# A <- quantile(Tgt1AllMed$residuals, 0.025)
# B <- quantile(Tgt1AllMed$residuals, 0.975)
#
# ggplot(Tgt1AllMed, aes(x= Tgt1AllMed$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt1 All Data (Medians)',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Tgt1AllMed))
#
# cat('MSD = ', signif(Tgt1AllMedMSD, 3))
#
# # Tgt2 ANOVAS ------------------
# # 2-way ANOVA ---------------
# Tgt2All.aov<- aov(Activity ~ Sample * Scale, data = filter(CmpdData, Assay == 'Tgt2', Scale != 'Raw'))
#
# Tgt2AllMSD <- 2 * sqrt(2) * sd(residuals(object = Tgt2All.aov))
#
# plot(Tgt2All.aov, 2)
#
# A <- quantile(Tgt2All.aov$residuals, 0.025)
# B <- quantile(Tgt2All.aov$residuals, 0.975)
#
# ggplot(Tgt2All.aov, aes(x= Tgt2All.aov$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt2 All Data',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Tgt2All.aov))
#
# cat('MSD = ', signif(Tgt2AllMSD, 3))
#
# # 1-way ANOVA (Mean) =========================
# Tgt2AllMean.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt2', Scale == 'Mean'))
#
# Tgt2AllMeanMSD <- 2 * sqrt(2) * sd(residuals(object = Tgt2AllMean.aov))
#
# plot(Tgt2AllMean.aov, 2)
#
# A <- quantile(Tgt2AllMean.aov$residuals, 0.025)
# B <- quantile(Tgt2AllMean.aov$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt2 All Data (Means)',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Tgt2AllMean.aov))
#
# cat('MSD = ', signif(Tgt2AllMeanMSD, 3))
#
# # 1-way ANOVA Median ============
# Tgt2AllMed.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt2', Scale == 'Median'))
#
# Tgt2AllMedMSD <- 2 * sqrt(2) * sd(residuals(object = Tgt2AllMed.aov))
#
# plot(Tgt2AllMed.aov, 2)
#
# A <- quantile(Tgt2AllMed.aov$residuals, 0.025)
# B <- quantile(Tgt2AllMed.aov$residuals, 0.975)
#
# ggplot(Tgt2AllMed.aov, aes(x= Tgt2AllMed.aov$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt2 All Data (Medians)',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Tgt2AllMed.aov))
#
# cat('MSD = ', signif(Tgt2AllMedMSD, 3))
#
# # Plate median MSDs ------------------------------
# # Samples summarized by plate and 16 plates/sample  -------------
# CmpdSampleData <- CmpdData %>%
#   filter(Scale != 'Raw') %>%
#   group_by(Assay, ReadId, Sample) %>%
#   mutate(Activity = if_else(Scale =='Mean', mean(Activity), median(Activity))) %>%
#   ungroup() %>%
#   select(-(13:19)) %>%
#   unique()
#
# # ANOVAs -----------
# # Tgt1 Sample--------------------
# Tgt1Sample.aov<- aov(Activity ~ Sample * Scale, data = filter(CmpdData, Assay == 'Tgt1'))
#
# Temp <- Tgt1Sample.aov
#
# Tgt1SampleMSD <- 2 * sqrt(2) * sd(residuals(object = Temp))
#
# plot(Temp, 2)
#
# A <- quantile(Temp$residuals, 0.025)
# B <- quantile(Temp$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt1 Plate Summarized Data',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Temp))
#
# cat('MSD = ', signif(Tgt1SampleMSD, 3))
#
# # 1-way ANOVA (Means) --------------
# Tgt1SampleMean.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt1', Scale == 'Mean'))
#
# Temp <- Tgt1SampleMean.aov
#
# Tgt1SampleMeanMSD <- 2 * sqrt(2) * sd(residuals(object = Temp))
#
# plot(Temp, 2)
#
# A <- quantile(Temp$residuals, 0.025)
# B <- quantile(Temp$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt1 Plate Summarized Data (Means',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Temp))
#
# cat('MSD = ', signif(Tgt1SampleMeanMSD, 3))
#
# # 1-way ANOVA Medians -------------------
# Tgt1SampleMedian.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt1', Scale == 'Median'))
#
# Temp <- Tgt1SampleMedian.aov
#
# Tgt1SampleMedianMSD <- 2 * sqrt(2) * sd(residuals(object = Temp))
#
# plot(Temp, 2)
#
# A <- quantile(Temp$residuals, 0.025)
# B <- quantile(Temp$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt1 Plate Summarized Data (Medians',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Temp))
#
# cat('MSD = ', signif(Tgt1SampleMedianMSD, 3))
#
# # Tgt2 ------------------------
# # 2-way ANOVA --------------------
# Tgt2Sample.aov<- aov(Activity ~ Sample * Scale, data = filter(CmpdData, Assay == 'Tgt2'))
#
# Temp <- Tgt2Sample.aov
#
# Tgt2SampleMSD <- 2 * sqrt(2) * sd(residuals(object = Temp))
#
# plot(Temp, 2)
#
# A <- quantile(Temp$residuals, 0.025)
# B <- quantile(Temp$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt2 Plate Summarized Data',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Temp))
#
# cat('MSD = ', signif(Tgt2SampleMSD, 3))
#
# # 1-way ANOVA Means --------------
# Tgt2SampleMean.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt2', Scale == 'Mean'))
#
# Temp <- Tgt2SampleMean.aov
#
# Tgt2SampleMeanMSD <- 2 * sqrt(2) * sd(residuals(object = Temp))
#
# plot(Temp, 2)
#
# A <- quantile(Temp$residuals, 0.025)
# B <- quantile(Temp$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt2 Plate Summarized Data (Means)',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Temp))
#
# cat('MSD = ', signif(Tgt2SampleMeanMSD, 3))
#
# # 1-way ANOVA Median -------------------------------
# Tgt2SampleMedian.aov<- aov(Activity ~ Sample, data = filter(CmpdData, Assay == 'Tgt2', Scale == 'Median'))
#
# Temp <- Tgt2SampleMedian.aov
#
# Tgt2SampleMedianMSD <- 2 * sqrt(2) * sd(residuals(object = Temp))
#
# plot(Temp, 2)
#
# A <- quantile(Temp$residuals, 0.025)
# B <- quantile(Temp$residuals, 0.975)
#
# ggplot(Temp, aes(x= Temp$residuals)) +
#   geom_histogram() +
#   geom_vline(xintercept = A, linetype = 'dashed', color = 'blue') +
#   geom_vline(xintercept = B, linetype = 'dashed', color = 'blue') +
#   labs(title = 'Model Tgt2 Plate Summarized Data (Medians)',
#        x = 'Residuals') +
#   theme_minimal()
#
# print(summary(Temp))
#
# cat('MSD = ', signif(Tgt2SampleMedianMSD, 3))
#
#
# # Dose Response Curves -----------------------------------------------
#
# ggplot(filter(CmpdData, ReadId == '05', Scale == 'Median'), aes(x = Conc, y = Activity, shape = DRSample, color = DRSample)) +
#   geom_point() +
#   geom_smooth(method = 'loess', se = FALSE) +
#   labs(title = 'Tgt2 Dose Curves (representative plate',
#        x = 'Concentration',
#        y = 'Activity') +
#   scale_x_log10() +
#   theme_minimal()
#
# ggplot(filter(CmpdData, ReadId == '21', Scale == 'Median'), aes(x = Conc, y = Activity, shape = DRSample, color = DRSample)) +
#   geom_point() +
#   geom_smooth(method = 'loess', se = FALSE) +
#   labs(title = 'Tgt1 Dose Curves (representative plate',
#        x = 'Concentration',
#        y = 'Activity') +
#   scale_x_log10() +
#   theme_minimal()
#
# # Residuals vs Activity --------------------
# CmpdData <- CmpdData %>%
#   group_by(ReadId, Sample, Scale) %>%
#   mutate(SumAct = if_else(Scale == 'Median', median(Activity), mean(Activity)),
#          Residual = SumAct - Activity,
#          Lim = if_else(Scale == 'Median', 2 * mad(Residual), 2 * sd(Residual)),
#          Extreme = abs(Residual > 50)
#   )%>%
#   ungroup()
#
# ggplot(CmpdData, aes(x = SumAct, y = Residual)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth() +
#   labs(x = "Activity",
#        y = 'Residual') +
#   theme_minimal() +
#   facet_grid(Assay ~ Scale, scales = 'free') +
#   labs(title = 'Residials vs. Sample Summarized Activity')
#
# ggplot(SamplePlateResiduals, aes(x = TrueAct, y = Residual)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth() +
#   labs(x = "Activity",
#        y = 'Residual') +
#   theme_minimal() +
#   facet_grid(Assay ~ Scale, scales = 'free') +
#   labs(title = 'Residials vs. Sample Summarized Activity')
