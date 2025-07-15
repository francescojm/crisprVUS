

clc<-read.xlsx('../../data/raw/CL_tissue_ctype_colors.xlsx',sheet = 2,rowNames = TRUE)

load('_allHits.RData')

hist(allHits$percBasalEXP_of_ps_cl,border=FALSE,col='darkcyan',main=paste('Basal expression percentile of the hosting gene\nin cell line(s) arbouring the DAMs'),ylim=c(0,100),xlab='th')

availableExpValues<-length(which(!is.na(allHits$percBasalEXP_of_ps_cl) & allHits$percBasalEXP_of_ps_cl> -Inf))


pp<-round(100*length(which(allHits$avgBasalEXP_fpkm_in_ps_cl>=1))/availableExpValues,2)

print(paste(pp, '% of DAMs are in genes expressed in the cell line in which the DAM is observed'))

pdf('exploration/figures/percExpressedDAMs.pdf',5,5)
pie(c(pp,100-pp),col=c('darkcyan','gray'),border=FALSE,labels = c('Expressed','not Expressed'))
dev.off()

print(paste(round(100*length(which(allHits$percBasalEXP_of_ps_cl>50))/availableExpValues,2), '% of DAMs are in genes whith a basal expression over the 50th percentile of the cell line in which the DAM is observed'))
print(paste(round(100*length(which(allHits$percBasalEXP_of_ps_cl>80))/availableExpValues,2), '% of DAMs are in genes whith a basal expression over the 80th percentile of the cell line in which the DAM is observed'))
print(paste(round(100*length(which(allHits$percBasalEXP_of_ps_cl==100))/availableExpValues,2), '% of DAMs are in genes whith a basal expression over the highest percentile of the cell line in which the DAM is observed'))

pdf('exploration/figures/DAMexpPercentiles.pdf',5,5)
hist(allHits$percBasalEXP_of_ps_cl,border=FALSE,col='darkcyan',main=paste('Basal expression percentile of the hosting gene\nin cell line(s) arbouring the DAMs'),ylim=c(0,100),xlab='th')
dev.off()


pdf('exploration/figures/essentialityDAMmatching.pdf',5,5)
hist(100*allHits$matching,border=FALSE,col='gray',main='')
dev.off()

length(which(allHits$matching==1))/length(allHits$matching)





pdf('exploration/figures/essentialityDAMmatchingHGp.pdf',5,6)
plot(-log10(allHits$hypTest_p),bg=adjustcolor("blue", alpha.f = 0.3),col=NA,pch=21,cex=2,frame.plot=FALSE,ylab='-log10(HG p)')
abline(h= -log10(0.05),lty=2)
dev.off()



cdg<-read.table('../../data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv',sep='\t',header=TRUE)
cdg<-unique(cdg$SYMBOL)

#allHits<-allHits[which(!is.element(allHits$GENE,cdg)),]


load('summary_byvar.RData')


library(ggplot2)
library(ggrepel)

set.seed(123)

# Build dataframe
df <- data.frame(
  x = 1:nrow(allHits),
  y = allHits$medFitEff,
  label = paste(allHits$GENE,allHits$var)
)

# Add transparent fill colours
df$fill <- adjustcolor(clc[allHits$ctype, 1], alpha.f = 0.6)

# Top 20 most negative fitness effects
top_points <- df[order(df$y), ][1:50, ]


utype<-unique(allHits$ctype)
xpos<-match(utype,allHits$ctype)

manual_labels <- data.frame(
  x = xpos,
  y = rep(-0.3,length(xpos)),
  label = utype,
  colours=clc[unique(allHits$ctype), 1]
)

pdf('All_Hits_summary.pdf',14,6)
ggplot(df, aes(x = x, y = y)) +
  geom_point(
    aes(fill = fill),
    shape = 21,
    colour = "black",
    size = 3,
    stroke = 0.5
  ) +
  geom_text_repel(
    data = top_points,
    aes(label = label),
    nudge_y = -0.1, 
    size=3.5,
    box.padding = 0.8,
    point.padding = 0.5,
    min.segment.length = 0.01,
    segment.color = "grey50",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +
  geom_text(
    data = manual_labels,
    aes(x = x, y = y, label = label,color=colours),
    size = 3.5,
    fontface = "italic",    # optional
    hjust = 0,
    angle = 90# or 0.5 or 1 depending on alignment
  ) +
  scale_fill_identity() +
  theme_minimal()
dev.off()


plot(manual_labels$x,manual_labels$y)
text(manual_labels$x,manual_labels$y,srt=90,labels = unique(allHits$ctype),col = clc[unique(allHits$ctype), 1])
