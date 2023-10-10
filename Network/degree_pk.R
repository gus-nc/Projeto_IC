###
### Calculating indices related to the P(k) degree distribution
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----------------------------------------------------------------------------
# Import the matrices and create a webs object
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
###-----------------------------------------------------------------------------
#1 Network Size
ind_size = list ()
for (env in envs) {
  ind_size[[env]] = dim(webs[[env]])[1]+dim(webs[[env]])[2] # Here defined as the number of nodes in the matrix
}
print(ind_size)

#2 Connectance
ind_connect = list()
for (env in envs) {
  ind_connect[[env]] = sum(webs[[env]] > 0)/prod(dim(webs[[env]])) # Total links / number of possible links
}
print(ind_connect)

#3 Simple Degree distribution
# Predator degree 
pred_k = list()
pred_hi = list()
for (env in envs) {
  pred_k[[env]][["distr"]] = colSums(webs[[env]][,colnames(webs[[env]])] > 0) # Create a list inside envlist() with the values
  pred_k[[env]][["mean"]] = mean(pred_k[[env]][["distr"]])
  pred_k[[env]][["max"]] = max(pred_k[[env]][["distr"]])
  pred_k[[env]][["min"]] = min(pred_k[[env]][["distr"]])
  pred_k[[env]][["sd"]] = sd(pred_k[[env]][["distr"]])
  pred_hi[[env]] = hist(pred_k[[env]][["distr"]], 0:max(pred_k[[env]][["distr"]]), plot = FALSE) # Histogram
}

# Prey item degree
prey_k = list()
prey_hi = list()
for (env in envs) {
  prey_k[[env]][["distr"]] = rowSums(webs[[env]][,colnames(webs[[env]])] > 0) # Create a list inside envlist() with the values
  prey_hi[[env]] = hist(prey_k[[env]][["distr"]], 0:max(prey_k[[env]][["distr"]]), plot = FALSE) # Histogram
}

# Plot all the histograms
# Eucalyptus Vs Atlântic Forest Predator P(k)
barplot(pred_hi[["freq_eu"]]$density, names.arg = c(1:max(pred_k[["freq_eu"]][["distr"]])),
        xlab = "Predator degree", ylab = "Frequency")
barplot(pred_hi[["freq_mt"]]$density, names.arg = c(1:max(pred_k[["freq_mt"]][["distr"]])),
        xlab = "Predator degree", ylab = "Frequency")

# Eucalyptus Vs Atlântic Forest Prey P(k)
barplot(prey_hi[["freq_eu"]]$density, names.arg = c(1:max(prey_k[["freq_eu"]][["distr"]])),
        xlab = "Prey degree", ylab = "Frequency")
barplot(prey_hi[["freq_mt"]]$density, names.arg = c(1:max(prey_k[["freq_mt"]][["distr"]])),
        xlab = "Prey degree", ylab = "Frequency")
