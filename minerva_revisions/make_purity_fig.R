library(ggplot2)
library(gridExtra)

d1 = read.csv('d1.test_3.eval.csv')
d1.base = read.csv('d1.baseline.eval.csv')
d2 = read.csv('d2.test_2.eval.csv')
d2.base = read.csv('d2.baseline.eval.csv')

d1 = d1[d1$size > 1,]
d2 = d2[d2$size > 1,]
d1.base = d1.base[d1.base$size >= 30,]
d2.base = d2.base[d2.base$size >= 30,]

d1$enhanced = 'enhanced'
d2$enhanced = 'enhanced'
d1.base$enhanced = 'standard'
d2.base$enhanced = 'standard'

d1 = rbind(d1, d1.base)
d2 = rbind(d2, d2.base)

p1 = ggplot(d1, aes(x=purity,group=enhanced,fill=as.factor(enhanced))) + 
  geom_density(position="identity", alpha=0.5) +
  scale_fill_discrete(name="Kind") +
  theme_bw() +
  ylab('Density') + 
  xlab('Purity') + 
  theme(text = element_text(size=30))
p2 = ggplot(d2, aes(x=purity,group=enhanced,fill=as.factor(enhanced))) + 
  geom_density(position="identity", alpha=0.5) +
  scale_fill_discrete(name="Kind") +
  theme_bw() +
  ylab('Density') + 
  xlab('Purity') + 
  theme(text = element_text(size=30))

h1 = ggplot(d1, aes(x=entropy,group=enhanced,fill=as.factor(enhanced))) + 
  geom_density(position="identity", alpha=0.5) +
  scale_fill_discrete(name="Kind") +
  theme_bw() +
  ylab('Density') + 
  xlab('Entropy') +
  xlim(0, 6) +
  theme(text = element_text(size=30))
h2 = ggplot(d2, aes(x=entropy,group=enhanced,fill=as.factor(enhanced))) + 
  geom_density(position="identity", alpha=0.5) +
  scale_fill_discrete(name="Kind") +
  theme_bw() +
  ylab('Density') + 
  xlab('Entropy') +
  xlim(0, 6) +
  theme(text = element_text(size=30))

p = grid.arrange(
  p1, h1,
  p2, h2,
  ncol=2
)

png('purity_figure_raw.png', width=1920, height=1080)
plot(p)
dev.off()