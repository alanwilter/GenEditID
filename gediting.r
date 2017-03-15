library(ggplot2)

# read data
protein_abundance_data <- read.csv('data/20170118-GED-EXP-00000001-ICW.csv')
protein_abundance_data

# calculate relative protein abundance, ratio 800 to 700
protein_abundance_data['ratio800to700'] = protein_abundance_data['Channel800']/protein_abundance_data['Channel700']
protein_abundance_data

# plot ratio against cell line
protein_abundance_plot <- ggplot(protein_abundance_data, aes(x = Target, y = ratio800to700)) +
  geom_violin() +
  geom_point(aes(colour=Target)) +
  theme_bw() + 
  ggtitle('InCellWestern') +
  xlab('Cell line') +
  ylab('Relative protein abundance (ratio 800 to 700 nm)') +
  theme(axis.text.x  = element_text(angle=15, vjust=0.5, size=10))

print(protein_abundance_plot)

