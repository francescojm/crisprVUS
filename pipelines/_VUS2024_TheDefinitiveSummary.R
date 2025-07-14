
clc<-read.xlsx('../../data/raw/CL_tissue_ctype_colors.xlsx',sheet = 2,rowNames = TRUE)

load('_allHits.RData')

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
