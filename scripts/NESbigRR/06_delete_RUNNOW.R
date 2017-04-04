#low plots
#coi1.cam
#[7]
jpeg(paste("plots/BcAt_MAF20_noogs_lowTR_", names(HEM.plotdata[7]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[7]))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste(names(HEM.plotdata[7]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[7]), sep="")), colour = "blue", linetype="longdash") +
  geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[7]), sep="")), colour = "blue", linetype="longdash") +
  geom_text(aes(0,get(paste("TH95pos_", names(HEM.plotdata[7]), sep="")), label = "95% Threshold", vjust = 1.5, hjust = .05), col = "blue")+
  geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[7]), sep="")), linetype="longdash") +
  geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[7]), sep="")), linetype="longdash") +
  geom_text(aes(0,get(paste("TH99pos_", names(HEM.plotdata[7]), sep="")), label = "99% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  expand_limits(y=0)+
  theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
dev.off()

#col0.cam
#[8]
jpeg(paste("plots/BcAt_MAF20_noogs_lowTR_", names(HEM.plotdata[8]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[8]))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste(names(HEM.plotdata[8]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[8]), sep="")), colour = "blue", linetype="longdash") +
  geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[8]), sep="")), colour = "blue", linetype="longdash") +
  geom_text(aes(0,get(paste("TH95pos_", names(HEM.plotdata[8]), sep="")), label = "95% Threshold", vjust = 1.5, hjust = .05), col = "blue")+
  geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[8]), sep="")), linetype="longdash") +
  geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[8]), sep="")), linetype="longdash") +
  geom_text(aes(0,get(paste("TH99pos_", names(HEM.plotdata[8]), sep="")), label = "99% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  expand_limits(y=0)+
  theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
dev.off()

#npr1.cam
#[9]
jpeg(paste("plots/BcAt_MAF20_noogs_lowTR_", names(HEM.plotdata[9]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[9]))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste(names(HEM.plotdata[9]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[9]), sep="")), colour = "blue", linetype="longdash") +
  geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[9]), sep="")), colour = "blue", linetype="longdash") +
  geom_text(aes(0,get(paste("TH95pos_", names(HEM.plotdata[9]), sep="")), label = "95% Threshold", vjust = 1.5, hjust = .05), col = "blue")+
  geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[9]), sep="")), linetype="longdash") +
  geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[9]), sep="")), linetype="longdash") +
  geom_text(aes(0,get(paste("TH99pos_", names(HEM.plotdata[9]), sep="")), label = "99% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  expand_limits(y=0)+
  theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
dev.off()

#anac055.cam
#[6]
jpeg(paste("plots/BcAt_MAF20_noogs_lowTR_", names(HEM.plotdata[6]), ".ManhattanPlot.jpg", sep=""), width=8, height=4, units='in', res=600)
ggplot(HEM.plotdata, aes(x=Index, y=HEM.plotdata[6]))+
  theme_bw()+
  geom_point(aes(color = factor(Chrom)))+
  labs(list(y="SNP Effect Estimate", x="Chromosome position", title=paste(names(HEM.plotdata[6]))))+
  guides(col = guide_legend(nrow = 8, title="Chromosome"))+
  geom_hline(yintercept=get(paste("TH95pos_", names(HEM.plotdata[6]), sep="")), colour = "blue", linetype="longdash") +
  geom_hline(yintercept=get(paste("TH95neg_", names(HEM.plotdata[6]), sep="")), colour = "blue", linetype="longdash") +
  geom_text(aes(0,get(paste("TH95pos_", names(HEM.plotdata[6]), sep="")), label = "95% Threshold", vjust = 1.5, hjust = .05), col = "blue")+
  geom_hline(yintercept=get(paste("TH99pos_", names(HEM.plotdata[6]), sep="")), linetype="longdash") +
  geom_hline(yintercept=get(paste("TH99neg_", names(HEM.plotdata[6]), sep="")), linetype="longdash") +
  geom_text(aes(0,get(paste("TH99pos_", names(HEM.plotdata[6]), sep="")), label = "99% Threshold", vjust = 1.5, hjust = .05), col = "black")+
  expand_limits(y=0)+
  theme(legend.position="none")+
  scale_x_continuous(name="Chromosome", breaks = c(1677889, 5253114, 9013367, 11074212, 13595791, 17206983, 20036067, 22404724, 24429409, 26804549, 28608225, 30154184, 31914256, 34033137, 35838514, 38953687), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))
dev.off()