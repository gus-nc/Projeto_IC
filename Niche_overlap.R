###
### Calculating Niche Overlap between predator species
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----------------------------------------------------------------------------
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

# Filter the Atlantic Forest matrices
cols_to_remove_mt = which(colSums(webs$freq_mt) < 5) # Only for freq bc it reveals the discrete number of preys
# Min of 5 prey items
if (length(cols_to_remove_mt) > 0) { # Check if is a int(0)
  webs$vol_mt = webs$vol_mt[, -cols_to_remove_mt] # Remove columns
  webs$freq_mt = webs$freq_mt[, -cols_to_remove_mt] # Remove columns
}

###-----------------------------------------------------------------------------
### Retrieve species that occur in both environments
species_x = Reduce(intersect, list(colnames(webs[["vol_eu"]]), colnames(webs[["vol_mt"]]))) # Identify common elements in the columns

# Volume
mlst_vol = list() # Create a list for logging each species matrix
for (sp in species_x) {
  eucalyptus_vol = webs[["vol_eu"]][,sp]  # Extract the column respective to the species for eucalyptus
  mata_vol = webs[["vol_mt"]][,sp] # Extract the column respective to the species for mata atlantica
  mlst_vol[[sp]] = matrix(c(eucalyptus_vol, mata_vol), # Generate matrix eu/mt for X-squared tests
                        nrow = nrow(webs[["vol_mt"]]),
                        byrow = FALSE)
  
  colnames(mlst_vol[[sp]]) = c("eu","mt") # Naming the columns 
}
# Frequency
mlst_freq = list() # Create a list for logging each species matrix
for (sp in species_x) {
  eucalyptus_freq = webs[["freq_eu"]][,sp]  # Extract the column respective to the species for eucalyptus
  mata_freq = webs[["freq_mt"]][,sp] # Extract the column respective to the species for mata atlantica
  mlst_freq[[sp]] = matrix(c(eucalyptus_freq, mata_freq), # Generate matrix eu/mt for X-squared tests
                          nrow = nrow(webs[["freq_mt"]]),
                          byrow = FALSE)
  
  colnames(mlst_freq[[sp]]) = c("eu","mt") # Naming the columns 
}

# Remove the lines that sum 0, the lines removed for vol and freq should be the same
for (sp in species_x) {
  rows_to_remove_vol = which(rowSums(mlst_vol[[sp]]) == 0) # Check the sum for each row volume
  rows_to_remove_freq = which(rowSums(mlst_freq[[sp]]) == 0) # Check the sum for each row volume
  mlst_vol[[sp]] = mlst_vol[[sp]][-rows_to_remove_vol,] # Create a new matrix without rows with sum 0
  mlst_freq[[sp]] = mlst_freq[[sp]][-rows_to_remove_freq,] # Create a new matrix without rows with sum 0
}

###-----------------------------------------------------------------------------
# Calculate Niche Overlap Pianka (1986)
over_vol = list() # Create a list for the niche overlap w/ volume
over_freq = list() # Create a list for the niche overlap w/ frequency
v_i = c() # Vector for storing the pij*pik results
v_eu = c() # Vector for storing the pij*pij results [eucalyptus]
v_mt = c() # Vector for storing the pik*pik results [mata]
for (sp in species_x) { # For each species
  for (l in 1:nrow(mlst_vol[[sp]])) { # go through each row in vol and retrieve
    pij_pik = (as.vector(mlst_vol[[sp]][l,1])/
               sum(mlst_vol[[sp]][,1]))*
            (as.vector(mlst_vol[[sp]][l,2])/
               sum(mlst_vol[[sp]][,2]))
    pij_2 = (as.vector(mlst_vol[[sp]][l,1])/
              sum(mlst_vol[[sp]][,1]))^2
    pik_2 = (as.vector(mlst_vol[[sp]][l,2])/
              sum(mlst_vol[[sp]][,2]))^2
    v_i = append(v_i, pij_pik) # append this value for each row to a vector
    v_eu = append(v_eu, pij_2) # append this value for each row to a vector
    v_mt = append(v_mt, pik_2) # append this value for each row to a vector
  }
  over_vol[[sp]]= sum(v_i)/sqrt(sum(v_eu)*sum(v_mt)) # Calculate the niche overlap w/ volume for sp
  
  for (l in 1:nrow(mlst_freq[[sp]])) { # go through each row in freq and retrieve
    pij_pik = (as.vector(mlst_freq[[sp]][l,1])/
                 sum(mlst_freq[[sp]][,1]))*
      (as.vector(mlst_freq[[sp]][l,2])/
         sum(mlst_freq[[sp]][,2]))
    pij_2 = (as.vector(mlst_freq[[sp]][l,1])/
               sum(mlst_freq[[sp]][,1]))^2
    pik_2 = (as.vector(mlst_freq[[sp]][l,2])/
               sum(mlst_freq[[sp]][,2]))^2
    v_i = append(v_i, pij_pik) # append this value for each row to a vector
    v_eu = append(v_eu, pij_2) # append this value for each row to a vector
    v_mt = append(v_mt, pik_2) # append this value for each row to a vector
  }
  over_freq[[sp]]= sum(v_i)/sqrt(sum(v_eu)*sum(v_mt)) # Calculate the niche overlap w/ freq for sp
}
###-----------------------------------------------------------------------------
#Producing the statistics table
df = data.frame(matrix(nrow = 3, ncol = 2)) 
rownames(df) = species_x 
colnames(df) = c("Volume", "Frequency")

for (sp in species_x) {
  df[sp,"Volume"] = over_vol[[sp]]
  df[sp,"Frequency"] = over_freq[[sp]]
}

print(df)

library(gt)
library(gtExtras)

# View formated data.frame
df |> 
  gt(rownames_to_stub = TRUE) |> 
  tab_header(
    title = "Niche Overlap ",
    subtitle = "Co-occuring species in Mata Atlântica and Eucalyptus"
  ) |> 
  fmt_number(columns = colnames(df), decimals = 4)

