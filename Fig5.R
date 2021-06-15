# Recreate Phil's plots in the same style as other figures.

T1EstTrueMeanN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt1'), aes(x=true.mean.est, y=PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 1)' )

T2EstTrueMeanN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt2'), aes(x=true.mean.est, y=PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 1)' )

T1EstTrueMedN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt1'), aes(x=true.mean.est, y=PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 1)' )

T2EstTrueMedN1 <- ggplot (filter(CmpdData3, Assay == 'Tgt2'), aes(x=true.mean.est, y=PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 1)' )

T1EstTrueMeanN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt1'), aes(x=true.mean.est, y=mean.PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 4)' )

T2EstTrueMeanN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt2'), aes(x=true.mean.est, y=mean.PctAct.Mean)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Mean)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 4)' )

T1EstTrueMedN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt1'), aes(x=true.mean.est, y=median.PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt1 (N = 4)' )

T2EstTrueMedN4 <- ggplot (filter(SummPerPlate2, Assay == 'Tgt2'), aes(x=true.mean.est, y=median.PctAct.Median)) +
  geom_point() +
  labs(x = 'Estimated True Pct. Activity',
       y = 'Pct. Activity (Median)') +
  theme_classic() +
  labs(title = 'Tgt2 (N = 4)' )

T1Scatter <- (T1EstTrueMeanN1/T1EstTrueMedN1/T1EstTrueMeanN4/T1EstTrueMedN4)

T2Scatter <- (T2EstTrueMeanN1/T2EstTrueMedN1/T2EstTrueMeanN4/T2EstTrueMedN4)

Fig5 <- (T1Scatter | T2Scatter) +
  plot_annotation(title = 'Figure 5. Measured Data vs. Estimated Truth', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 5.jpg', plot = Fig5, height = 8, width = 6, units = 'in', dpi = 300)

