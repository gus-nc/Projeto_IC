#-------------------------------------------------------------------------------
# Script to Produce Chi-square test of independence for 
# each species that occur in both environments
#-------------------------------------------------------------------------------

setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Augusto's Laptop

# Import tables .csv files
links_eu = read.csv("Diet_E_Vol_Xincomp.csv", header = T, row.names = 1)
links_mt = read.csv("Diet_M_Vol_Xincomp.csv", header = T, row.names = 1)

links_eu = t(as.matrix(links_eu)) #Eucalyptus matrix
links_mt = t(as.matrix(links_mt)) #Mata matrix

# Remove species with less than 5 items
freq_eu_x =t(as.matrix((read.csv("Diet_E_Freq_Xincomp.csv", header = T, row.names = 1))))
freq_mt_x =t(as.matrix((read.csv("Diet_M_Freq_Xincomp.csv", header = T, row.names = 1))))

col_sum_eu = colSums(freq_eu_x) # Check the sum for each column for eucalyptus
col_sum_mt = colSums(freq_mt_x) # Check the sum for each column for mata

cols_to_remove_eu = which(col_sum_eu < 5) # Identify columns with sum < 5 for eu
cols_to_remove_mt = which(col_sum_mt < 5) # Identify columns with sum < 5 for mt

links_eu <- links_eu[, -cols_to_remove_eu] # Remove columns for eu
links_mt <- links_mt[, -cols_to_remove_mt] # Remove columns for mt

#-------------------------------------------------------------------------------

# Retrieve species that occur in both environments 
species = Reduce(intersect, list(colnames(links_eu), colnames(links_mt))) # Identify common elements in the columns
mlst = list() # Create a list for logging each species matrix

for (sp in species) {
  eucalyptus = links_eu[,sp]  # Extract the column respective to the species for eucalyptus
  mata = links_mt[,sp] # Extract the column respective to the species for mata atl?ntica
  mlst[[sp]] = matrix(c(eucalyptus, mata), # Generate matrix eu/mt for X-squared tests
                      nrow = nrow(links_eu),
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
lstchisq = list() # Create a list to store the chisq result for each species
df = data.frame(matrix(nrow = 3, ncol = 0)) # Create a data frame to plot the meaningful results from chisq.test()
rownames(df) = c("statistic", "parameter", "p.value") 

for (sp in species) {
  lstchisq[[sp]] = chisq.test(mlstnonzero[[sp]], correct = FALSE)  # No continuity correction for 2 x 2
  x = c(unlist(lstchisq[[sp]]["statistic"]),unlist(lstchisq[[sp]]["parameter"]),
        unlist(lstchisq[[sp]]["p.value"]))
  df[,sp] = x # Add a new column for each species
}

df = as.data.frame(t(df)) # Transpose the matrix
print(df)

library(gt)
library(gtExtras)

# View formated data.frame
df |> 
  gt(rownames_to_stub = TRUE) |> 
  tab_header(
    title = "Chisquare Test of Independence",
    subtitle = "Co-occuring species in Mata Atlântica and Eucalyptus"
  ) |> 
  fmt_number(columns = c("p.value","statistic"), decimals = 3)

