library(GGally)
library(network)
library(sna)

df.connectivity <- fread(io$paga.connectivity) %>%
  matrix.please %>% .[opts$celltypes,opts$celltypes]

df.coordinates <- fread(io$paga.coordinates) %>% 
  matrix.please %>% .[opts$celltypes,]

# Parse data
df.connectivity[df.connectivity<0.25] <- 0

# Create network
net.paga = network(df.connectivity)

# Define coordinates
# x = gplot.layout.fruchtermanreingold(net.paga, NULL)
net.paga %v% "x" = df.coordinates[, 1]
net.paga %v% "y" = df.coordinates[, 2]