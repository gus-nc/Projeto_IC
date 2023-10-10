#-----------------------------------------------------------------------------
# Plotting the preliminary matrices 
#-----------------------------------------------------------------------------

setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC

###-----------------------------------------------------------------------------
# Using igraph package
library("igraph")
# Import data as matrix for use in igraph
links_eu_vol = as.matrix(read.csv("Diet_E_Vol_Xincomp.csv", header = T, row.names = 1))
links_mt_vol = as.matrix(read.csv("Diet_M_Vol_Xincomp.csv", header = T, row.names = 1))
links_eu_freq = as.matrix(read.csv("Diet_E_Freq_Xincomp.csv", header = T, row.names = 1))
links_mt_freq = as.matrix(read.csv("Diet_M_Freq_Xincomp.csv", header = T, row.names = 1))

# Convert network in a igraph object
netw_eu_vol = graph_from_incidence_matrix(links_eu_vol, weighted = TRUE)
E(netw_eu_vol)$width = log10(E(netw_eu_vol)$weight)
netw_mt_vol = graph_from_incidence_matrix(links_mt_vol, weighted = TRUE)
E(netw_mt_vol)$width = log10(E(netw_mt_vol)$weight)
netw_eu_freq = graph_from_incidence_matrix(links_eu_freq, weighted = TRUE)
E(netw_eu_freq)$width = E(netw_eu_freq)$weight
netw_mt_freq = graph_from_incidence_matrix(links_mt_freq, weighted = TRUE)
E(netw_eu_freq)$width = E(netw_eu_freq)$weight


# Plot the networks for volumes matrices
par(mfrow=c(1,2))
plot(netw_eu_vol, vertex.label=NA, vertex.size=8, layout=layout.bipartite)
plot(netw_mt_vol, vertex.label=NA, vertex.size=8, layout=layout.bipartite)

# Plot the networks for frequency matrices
par(mfrow=c(1,2))
plot(netw_eu_freq, vertex.label=NA, vertex.size=8, layout=layout.bipartite)
plot(netw_mt_freq, vertex.label=NA, vertex.size=8, layout=layout.bipartite)
###-----------------------------------------------------------------------------
# Using bipartite package
# Load required packages for bipartite
library(permute)
library(lattice)
library(vegan)
library(statnet.common)
library(network)
library(sna)
library(bipartite) # Main package

#A matrix > higher trophic level species=anurans (columns) and lower trophic level species=prey (rows).
# Read all the interaction network matrices into a list
linksb_eu_vol = t(as.matrix(read.csv("Diet_E_Vol_Xincomp.csv", header = T, row.names = 1)))
linksb_mt_vol = t(as.matrix(read.csv("Diet_M_Vol_Xincomp.csv", header = T, row.names = 1)))
linksb_eu_freq = t(as.matrix(read.csv("Diet_E_Freq_Xincomp.csv", header = T, row.names = 1)))
linksb_mt_freq = t(as.matrix(read.csv("Diet_M_Freq_Xincomp.csv", header = T, row.names = 1)))


webs = list(linksb_eu_vol, linksb_mt_vol, linksb_eu_freq, linksb_eu_freq)
names(webs) = c("vol_eu", "vol_mt", "freq_eu", "freq_mt") 

# Visualize the observed networks from the datasets
#For Eucalipto
visweb(webs$freq_eu, type = "none")
visweb(webs$W_eu, type = "nested") #It kinda works but doesn't seem nested
visweb(webs$W_eu, type = "diagoal") #It kinda works
# For Mata Atlantica
visweb(webs$W_mt, type = "none")
visweb(webs$W_mt, type = "nested") #It kinda works but doesn't seem nested
visweb(webs$W_mt, type = "diagoal") #It kinda works

# Plot the Network
plotweb(webs$freq_eu, y.width.high = 0.05, y.width.low = 0.05, text.rot=90, low.y = 0.8, 
        high.y = 1.3, col.low = "green", col.high = "blue")

plotweb(webs$W_mt, y.width.high = 0.05, y.width.low = 0.05, text.rot=90, low.y = 0.8, 
        high.y = 1.3, col.low = "green", col.high = "blue")

