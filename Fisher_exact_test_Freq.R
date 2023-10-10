#-------------------------------------------------------------------------------
# Script to Produce Fisher's exact test for 
# each species that occur in both environments
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Augusto's Windows PC
setwd("/Volumes/GoogleDrive/Other computers/Meu laptop/Code/Projeto_IC/") # Augusto's MacOS Laptop
###-----------------------------------------------------------------------------
# 1. With Incomplete Prey Items

# Import tables .csv files
links_eu = read.csv("Diet_E_Freq.csv", header = T, row.names = 1)
links_mt = read.csv("Diet_M_Freq.csv", header = T, row.names = 1)

links_eu = t(as.matrix(links_eu)) #Eucalyptus matrix
links_mt = t(as.matrix(links_mt)) #Mata matrix

# Remove species with less than 5 items
col_sum_eu = colSums(links_eu) # Check the sum for each column for eucalyptus
col_sum_mt = colSums(links_mt) # Check the sum for each column for mata

cols_to_remove_eu = which(col_sum_eu < 5) # Identify columns with sum < 5 for eu
cols_to_remove_mt = which(col_sum_mt < 5) # Identify columns with sum < 5 for mt

links_eu = links_eu[, -cols_to_remove_eu] # Remove columns for eu
links_mt = links_mt[, -cols_to_remove_mt] # Remove columns for mt



# Retrieve species that occur in both environments
species = Reduce(intersect, list(colnames(links_eu), colnames(links_mt))) # Identify common elements in the columns

mlst = list() # Create a list for logging each species matrix

for (sp in species) {
  eucalyptus = links_eu[,sp]  # Extract the column respective to the species for eucalyptus
  mata = links_mt[,sp] # Extract the column respective to the species for mata atl?ntica
  mlst[[sp]] = matrix(c(eucalyptus, mata), # Generate matrix eu/mt for X-squared tests
                      nrow = 29,
                      byrow = FALSE)
  
  colnames(mlst[[sp]]) = c("eu","mt") # Naming the columns 
}
# Remove the lines that sum 0, because they do not work for the statistic test
mlstnonzero = list() # list for logging modified species matrix

for (sp in species) {
  r_sum = rowSums(mlst[[sp]]) # Check the sum for each row
  l = c() # Create a vector to store the rows with sum 0
  for (i in 1:nrow(mlst[[sp]])) { # Iterate through each row
    if (r_sum[i] == 0) { 
      l = c(l, i) # If row = 0 add the row number to a vector
    }
  }
  mlstnonzero[[sp]] = mlst[[sp]][-l,] # Create a new matrix without rows with sum 0
}



#Producing the statistics table
lstfish = list() # Create a list to store the chisq result for each species
df = data.frame(matrix(nrow = 1, ncol = 0)) # Data frame to plot results from chisq.test()
rownames(df) = "p.value" 

for (sp in species) {
  lstfish[[sp]] = fisher.test(mlstnonzero[[sp]], # Monte Carlo correction applied
                      alternative = "two.sided", simulate.p.value = TRUE, B = 2000) 
  df[,sp] = unlist(lstfish[[sp]]["p.value"]) # Add a new column for each species
}

df = as.data.frame(t(df)) # Transpose the data frame>matrix>dataframe again
print(df)



# View formated data.frame
library(gt)
library(gtExtras)
df |> 
  gt(rownames_to_stub = TRUE) |> 
  tab_header(
    title = "Fisher's Exact Test",
    subtitle = "Co-occuring species in Mata Atlântica and Eucalyptus"
  ) |> 
  fmt_number(columns = "p.value", decimals = 4)

###-------------------------------------------------------------------------------
# 2. Without Incomplete Preys
###-------------------------------------------------------------------------------
# Import tables .csv files
links_eu_x = read.csv("Diet_E_Freq_Xincomp.csv", header = T, row.names = 1)
links_mt_x = read.csv("Diet_M_Freq_Xincomp.csv", header = T, row.names = 1)

links_eu_x = t(as.matrix(links_eu_x)) #Eucalyptus matrix
links_mt_x = t(as.matrix(links_mt_x)) #Mata matrix

# Remove species with less than 5 items
col_sum_eu_x = colSums(links_eu_x) # Check the sum for each column for eucalyptus
col_sum_mt_x = colSums(links_mt_x) # Check the sum for each column for mata

cols_to_remove_eu_x = which(col_sum_eu_x < 5) # Identify columns with sum < 5 for eu
cols_to_remove_mt_x = which(col_sum_mt_x < 5) # Identify columns with sum < 5 for mt

links_eu_x <- links_eu_x[, -cols_to_remove_eu_x] # Remove columns for eu
links_mt_x <- links_mt_x[, -cols_to_remove_mt_x] # Remove columns for mt

#-------------------------------------------------------------------------------
# Retrieve species that occur in both environments
species_x = Reduce(intersect, list(colnames(links_eu_x), colnames(links_mt_x))) # Identify common elements in the columns

mlst_x = list() # Create a list for logging each species matrix

for (sp in species_x) {
  eucalyptus_x = links_eu_x[,sp]  # Extract the column respective to the species for eucalyptus
  mata_x = links_mt_x[,sp] # Extract the column respective to the species for mata atl?ntica
  mlst_x[[sp]] = matrix(c(eucalyptus_x, mata_x), # Generate matrix eu/mt for X-squared tests
                      nrow = 24,
                      byrow = FALSE)
  
  colnames(mlst_x[[sp]]) = c("eu","mt") # Naming the columns 
}
# Remove the lines that sum 0, because they do not work for the statistic test
mlstnonzero_x = list() # list for logging modified species matrix

for (sp in species_x) {
  r_sum_x = rowSums(mlst_x[[sp]]) # Check the sum for each row
  l_x = c() # Create a vector to store the rows with sum 0
  for (i in 1:nrow(mlst_x[[sp]])) { # Iterate through each row
    if (r_sum_x[i] == 0) { 
      l_x = c(l_x, i) # If row = 0 add the row number to a vector
    }
  }
  mlstnonzero_x[[sp]] = mlst_x[[sp]][-l_x,] # Create a new matrix without rows with sum 0
}

#-------------------------------------------------------------------------------
#Producing the statistics table
lstfish_x = list() # Create a list to store the chisq result for each species
df_x = data.frame(matrix(nrow = 1, ncol = 0)) # Data frame to plot results from chisq.test()
rownames(df_x) = "p.value" 

for (sp in species_x) {
  lstfish_x[[sp]] = fisher.test(mlstnonzero_x[[sp]], # Monte Carlo correction applied
                              alternative = "two.sided", simulate.p.value = TRUE, B = 2000) 
  df_x[,sp] = unlist(lstfish_x[[sp]]["p.value"]) # Add a new column for each species
}

df_x = as.data.frame(t(df_x)) # Transpose the data frame>matrix>dataframe again
print(df_x)

library(gt)
library(gtExtras)

# View formated data.frame
df_x |> 
  gt(rownames_to_stub = TRUE) |> 
  tab_header(
    title = "Fisher's Exact Test",
    subtitle = "Co-occuring species in Mata Atlântica and Eucalyptus"
  ) |> 
  fmt_number(columns = "p.value", decimals = 4)

