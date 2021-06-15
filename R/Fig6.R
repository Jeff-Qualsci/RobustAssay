# Figure 6 ROC curves

T1N1ppv = ggplot(plot.Indiv.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +
  theme(legend.position = 'none')

T1N4ppv <- ggplot(plot.results.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +
  theme(legend.position = 'none')

T1N1sel <- ggplot(plot.Indiv.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +
  theme(legend.position = 'none')

T1N4sel <- ggplot(plot.results.Tgt1, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt1, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +
  theme(legend.position = 'none')

T2N1ppv = ggplot(plot.Indiv.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic()

T2N4ppv <- ggplot(plot.results.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("PPV or NPV") +
  theme_classic() +


T2N1sel <- ggplot(plot.Indiv.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=1/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic()

T2N4sel <- ggplot(plot.results.Tgt2, aes (x=cutoff, y=Value, color=Result, linetype=Scale)) +
  geom_line(lwd=1.1) + theme_bw() + ylim (0.8, 1) +
  ggtitle ("Tgt2, N=4/Sample") + xlab("Activity Cutoff, %") + ylab("Sensitivity or Specificity") +
  theme_classic() +


T1ROC <- (T1N1sel/T1N1ppv/T1N4sel/T1N4ppv)

T2ROC <- (T2N1sel/T2N1ppv/T2N4sel/T2N4ppv)

Fig6 <- (T1ROC | T2ROC) +
  plot_annotation(title = 'Figure 6. ROC Analysis', tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave('Figures/Weidner Fig 6.jpg', plot = Fig6, height = 8, width = 6, units = 'in', dpi = 300)


