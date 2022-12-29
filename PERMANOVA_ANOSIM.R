# vegan package learning
# install.packages("vegan")
library(vegan)
data(dune, dune.env) # example
?dune # look data description

# Homogeneity test
groupDune <- as.factor(dune.env$Management)
dune.dist <- vegdist(dune, method = 'euclidean') 
modDune <- betadisper(dune.dist, groupDune)
anova(modDune)

# PERMANOVA
adonis2(dune ~ Management, data = dune.env, method = 'bray')
adonis2(dune ~ Management, data = dune.env, method = 'euclidean')

# ANOSIM
dune.dist <- vegdist(dune, method = 'euclidean') 
dune.ano <- with(dune.env, anosim(dune.dist, Management)) 
dune.ano
summary(dune.ano) 
plot(dune.ano)



