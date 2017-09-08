#quick plot

#HEM.plotdata <- HEM.plotdata[which(abs(HEM.plotdata$LA2093) > 3e-05),]




for (i in c(7)){
jpeg(paste("paper/plots/FigR6/bw_Sl_LesionSize_trueMAF20_NA10_lowTR_", names(HEM.plotdata[i]), ".ManhattanPlot.jpg", sep=""), width=8, height=5, units='in', res=600)
plot(ggplot(HEM.plotdata, aes(x=Index, y=100*HEM.plotdata[,i]))+
         theme_bw()+
         colScale+
         geom_point(aes(color = factor(Chrom)))+
         labs(list(y=expression(paste("Estimated Effect Size (",mm^{2},")")), title=paste("Lesion Size on ", names(HEM.plotdata[i]))))+
                guides(col = guide_legend(nrow = 8, title="Chromosome"))+
                geom_hline(yintercept=100*get(paste("TH95pos_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
                geom_hline(yintercept=100*get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), colour = "black", lty=3) +
                geom_text(aes(0,100*get(paste("TH95neg_", names(HEM.plotdata[i]), sep="")), label = "95% Threshold", vjust = 2.2, hjust=.05), col = "black")+
                geom_hline(yintercept=100*get(paste("TH99pos_", names(HEM.plotdata[i]), sep="")), lty=2) +
                geom_hline(yintercept=100*get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), lty=2) +
                geom_text(aes(0,100*get(paste("TH99neg_", names(HEM.plotdata[i]), sep="")), label = "99% Threshold", vjust = 5.2, hjust=.05), col = "black")+
                theme(legend.position="none")+
                theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
                #NA10 chromosomes
                scale_x_continuous(name="Chromosome", breaks = c(1674920, 5242762, 8999082, 11057880, 13580364, 17181767, 20009659, 22371413, 24388631, 26765466, 28559407, 30104939, 31864287, 33980281, 35775028, 38876579), labels = c("1", "2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14", "15", "16"))+
       geom_vline(xintercept=1099438, lty=2)+
       geom_vline(xintercept=1836245, lty=2)+
       geom_vline(xintercept=31154483, lty=2)+
       geom_vline(xintercept=33853054, lty=2)+
       geom_vline(xintercept=36555407, lty=2)+
       geom_vline(xintercept=41350409, lty=2)+
       geom_vline(xintercept=5703231, lty=2)+
       geom_vline(xintercept=6589294, lty=2)+
       geom_vline(xintercept=7955289, lty=2)+
       geom_vline(xintercept=11188054, lty=2)+
       geom_vline(xintercept=22332692, lty=2)+
       geom_vline(xintercept=24530790, lty=2)
)
#dev.off()
}
              

