

library(ggplot2)
library(ggthemes)


bc_sizes = read.csv('bc_sizes.csv')
colnames(bc_sizes) = c('size')
bc_sizes = bc_sizes / 2

p = ggplot(bc_sizes, aes(bc_sizes$size)) +
  geom_histogram(binwidth=1) +
  theme_tufte() +
  xlab('Read Cloud Size') +
  ylab('Count') +
  #scale_y_log10() +
  xlim(0, 50)

print(p)