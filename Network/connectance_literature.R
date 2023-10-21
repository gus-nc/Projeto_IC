###
### Calculating Connectance relationships with literature data Ceron et al (2019)
###-----Set a Directory---------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----Import and Filter the Matrices------------------------------------------
# Import the matrices and create a webs object
# A matrix as > higher trophic level species=anurans (columns) 
# and lower trophic level species=prey (rows).

# Read all the interaction network matrices into a list
links_eu_freq = t(as.matrix(read.csv("Diet_E_Freq_Xincomp.csv", header = T, row.names = 1)))
links_mt_freq = t(as.matrix(read.csv("Diet_M_Freq_Xincomp.csv", header = T, row.names = 1)))

webs = list(links_eu_freq, links_mt_freq) # create a list w/ the matrices
envs = c("freq_eu", "freq_mt") # names of each env and type combination
names(webs) = envs # name the list

# Filter the Eucalyptus matrices
cols_to_remove_eu = which(colSums(webs$freq_eu) < 5) # Only for freq bc it reveals the discrete number of preys
# Min of 5 prey items p/ species
if (length(cols_to_remove_eu) > 0) { # Check if is a int(0)
  webs$vol_eu = webs$vol_eu[, -cols_to_remove_eu] # Remove columns
  webs$freq_eu = webs$freq_eu[, -cols_to_remove_eu] # Remove columns
}

rows_to_remove_eu = which(rowSums(webs$freq_eu) == 0) # Check if there are items without a predator
if (length(rows_to_remove_eu) > 0) { # Check if is a int(0)
  webs$vol_eu = webs$vol_eu[-rows_to_remove_eu,] # Remove columns
  webs$freq_eu = webs$freq_eu[-rows_to_remove_eu,] # Remove columns
}
# Filter the Atlantic Forest matrices
cols_to_remove_mt = which(colSums(webs$freq_mt) < 5) # Only for freq bc it reveals the discrete number of preys
# Min of 5 prey items
if (length(cols_to_remove_mt) > 0) { # Check if is a int(0)
  webs$vol_mt = webs$vol_mt[, -cols_to_remove_mt] # Remove columns
  webs$freq_mt = webs$freq_mt[, -cols_to_remove_mt] # Remove columns
}

rows_to_remove_mt = which(rowSums(webs$freq_mt) == 0) # Check if there are items without a predator
if (length(rows_to_remove_mt) > 0) { # Check if is a int(0)
  webs$vol_mt = webs$vol_mt[-rows_to_remove_mt,] # Remove columns
  webs$freq_mt = webs$freq_mt[-rows_to_remove_mt,] # Remove columns
}

###------Network Size and Connectance-------------------------------------------
# 1 Network Size
ind_size = list ()
for (env in envs) {
  ind_size[[env]] = dim(webs[[env]])[1]+dim(webs[[env]])[2] # Here defined as the number of nodes in the matrix
}
print(ind_size)

# 2 Connectance
ind_connect = list()
for (env in envs) {
  ind_connect[[env]] = sum(webs[[env]] > 0)/prod(dim(webs[[env]])) # Total links / number of possible links
}
print(ind_connect)
###------Import Ceron´s data base-----------------------------------------------
library(dplyr)
library(tidyverse)

# Produce one big fat matrix
data = read.csv("Network/Ceron_Sup_Augusto.csv", header=T, dec= ".", row.names="id")

plot(data$size, data$com, 
     xlab = "Network Size",
     ylab = "Connectance (C)",
     ylim = c(0.2,1),
     pch = 19,   # Set the point character (you can change this)
     col = "blue" # Set the point color (you can change this)
)
points(ind_size$freq_eu, ind_connect$freq_eu, pch = 19, col = "red")
points(ind_size$freq_mt, ind_connect$freq_mt, pch = 19, col = "red")
