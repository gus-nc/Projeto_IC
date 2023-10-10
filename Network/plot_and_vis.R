###-----------------------------------------------------------------------------
### Plotting and Visualizing the matrices with the bipartite package

# Load required packages for bipartite
library(permute)
library(lattice)
library(vegan)
library(statnet.common)
library(network)
library(sna)
library(bipartite) # Main package

###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----------------------------------------------------------------------------
# A matrix as > higher trophic level species=anurans (columns) 
# and lower trophic level species=prey (rows).

# Read all the interaction network matrices into a list
links_eu_vol = t(as.matrix(read.csv("Diet_E_Vol_Xincomp.csv", header = T, row.names = 1)))
links_mt_vol = t(as.matrix(read.csv("Diet_M_Vol_Xincomp.csv", header = T, row.names = 1)))
links_eu_freq = t(as.matrix(read.csv("Diet_E_Freq_Xincomp.csv", header = T, row.names = 1)))
links_mt_freq = t(as.matrix(read.csv("Diet_M_Freq_Xincomp.csv", header = T, row.names = 1)))

webs = list(links_eu_vol, links_mt_vol, links_eu_freq, links_mt_freq) # create a list w/ the matrices
envs = c("vol_eu", "vol_mt", "freq_eu", "freq_mt") # names of each env and type combination
names(webs) = envs # name the list

# Filter the Eucalyptus matrices
col_sum_eu = colSums(webs$freq_eu) # Only for freq bc it reveals the discrete number of preys
cols_to_remove_eu = which(col_sum_eu < 5) # Min of 5 prey items
if (length(cols_to_remove_eu) > 0) { # Check if is a int(0)
  webs$vol_eu = webs$vol_eu[, -cols_to_remove_eu] # Remove columns
  webs$freq_eu = webs$freq_eu[, -cols_to_remove_eu] # Remove columns
}

rows_to_remove_eu = which(rowSums(webs$freq_eu) == 0) # Check if there are items without a predator
if (length(rows_to_remove_eu) > 0) { # Check if is a int(0)
  webs$vol_eu = webs$vol_eu[-rows_to_remove_eu,] # Remove lines
  webs$freq_eu = webs$freq_eu[-rows_to_remove_eu,] # Remove lines
}
# Filter the Atlantic Forest matrices
col_sum_mt = colSums(webs$freq_mt) # Only for freq bc it reveals the discrete number of preys
cols_to_remove_mt = which(col_sum_mt < 5)  # Min of 5 prey items
if (length(cols_to_remove_mt) > 0) { # Check if is a int(0)
  webs$vol_mt = webs$vol_mt[, -cols_to_remove_mt] # Remove columns
  webs$freq_mt = webs$freq_mt[, -cols_to_remove_mt] # Remove columns
}

rows_to_remove_mt = which(rowSums(webs$freq_mt) == 0) # Check if there are items without a predator
if (length(rows_to_remove_mt) > 0) { # Check if is a int(0)
  webs$vol_mt = webs$vol_mt[-rows_to_remove_mt,] # Remove lines
  webs$freq_mt = webs$freq_mt[-rows_to_remove_mt,] # Remove lines
}
###-----------------------------------------------------------------------------
### Visualize the observed networks from the datasets as incidence matrices
# Re-order the interaction matrices
for (env in envs) {
  row_order = order(rowSums(webs[[env]])) # Calculate row sums
  webs[[env]] = webs[[env]][row_order, ] # Use the index vector to reorder the rows of the matrix
}

#For Eucalipto Compar
par(mfrow=c(1,2))
visweb(webs$freq_eu, plotsize = 6, labsize = 1) #It kinda works but doesn't seem nested
visweb(webs$vol_eu, plotsize = 6, labsize = 1)
# For Mata Atlantica Compar
par(mfrow=c(1,2))
visweb(webs$freq_mt,  plotsize = 9, labsize = 1)
visweb(webs$vol_mt, plotsize = 9, labsize = 1)

# For Frequency
par(mar=c(0.1,0.1,0.1,0.1))
layout(matrix(c(1,2,1,3),ncol=2),heights=c(0.5,3,0.5,3))
plot.new()
text(0.5,0.5,"Eucalyptus Vs Atlântic Forest Frequency Matrices",cex=2,font=2)  # Add title 
visweb(webs$freq_eu, type = "nested", plotsize = 6, labsize = 1) #It kinda works but doesn't seem nested
visweb(webs$freq_mt, type = "nested", plotsize = 6, labsize = 1)

# For Volume
par(mar=c(0.1,0.1,0.1,0.1))
layout(matrix(c(1,2,1,3),ncol=2),heights=c(0.5,3,0.5,3))
plot.new()
text(0.5,0.5,"Eucalyptus Vs Atlântic Forest Volume Matrices",cex=2,font=2)  # Add title 
visweb(webs$vol_eu, type = "nested", plotsize = 6, labsize = 1) #It kinda works but doesn't seem nested
visweb(webs$vol_mt, type = "nested", plotsize = 6, labsize = 1)

###-----------------------------------------------------------------------------
# Visualize the observed networks from the datasets as incidence matrices

# Plot the Network Frequency
plotweb(webs$freq_eu, y.width.high = 0.05, y.width.low = 0.05, text.rot=90, low.y = 0.8, 
        high.y = 1.3, col.low = "green", col.high = "blue")

plotweb(webs$freq_mt, y.width.high = 0.05, y.width.low = 0.05, text.rot=90, low.y = 0.8, 
        high.y = 1.3, col.low = "green", col.high = "blue")

# Plot the Network Volume
plotweb(webs$vol_eu, y.width.high = 0.05, y.width.low = 0.05, text.rot=90, low.y = 0.8, 
        high.y = 1.3, col.low = "green", col.high = "blue")

plotweb(webs$vol_mt, y.width.high = 0.05, y.width.low = 0.05, text.rot=90, low.y = 0.8, 
        high.y = 1.3, col.low = "green", col.high = "blue")

