###
### Calculating Nestedness for the networks
###-----Set Directory-----------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----Organize Observed-------------------------------------------------------
# Import the matrices and create a webs list
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

# Re-order the interaction matrices
for (env in envs) {
  row_order = order(rowSums(webs[[env]]), decreasing = TRUE) # Calculate row sums
  webs[[env]] = webs[[env]][row_order, ] # Use the index vector to reorder the rows of the matrix
  col_order = order(colSums(webs[[env]]), decreasing = TRUE) # Calculate col sums
  webs[[env]] = webs[[env]][,col_order] # Use the index vector to reorder the cols of the matrix
}

###-----Obs wNODF---------------------------------------------------------------
### Calculate the brute value of WNODF for each env empirical network
WNODFc = list() # List for storing the sum values for cols
WNODFr = list() # List for storing the sum values for rowws
WNODF_obs = list() # List for storing the final value for each matrix

for (env in envs) { # For each environment
  # WNODFc values for cols
  WNODFc_temp = 0 # Create a temporary object for WNODFc
  for (i in 1:(ncol(webs[[env]])-1)) { # for every i:n-1
    for (j in (i+1):ncol(webs[[env]])) { # and for every i+1:n
      if(sum(webs[[env]][,i]) <= sum(webs[[env]][,j])) { # if F(ci) <= F(cj)
        N_paired = 0
      }
      else {
        k_ij = 0 # number of 1's in ci = cj, if cj != 0?
        N_j = sum(webs[[env]][,j] > 0) # Cells > 0 in cj
        for (r in 1:nrow(webs[[env]])) { # for every row in t cols sum k_ij
          if (webs[[env]][r,i] > webs[[env]][r,j] && webs[[env]][r,j] != 0){
            k_ij = k_ij+1
          }
        }
        N_paired = k_ij/N_j
      }
      WNODFc_temp = WNODFc_temp + N_paired 
    }
  }
  WNODFc[[env]] = WNODFc_temp*100
  # WNODFr 
  WNODFr_temp = 0
  for (i in 1:(nrow(webs[[env]])-1)) {
    for (j in (i+1):nrow(webs[[env]])) {
      if(sum(webs[[env]][i,]) <= sum(webs[[env]][j,])) {
        N_paired = 0
      }
      else {
        k_ij = 0
        N_j = sum(webs[[env]][j,] > 0)
        for (c in 1:ncol(webs[[env]])) {
          if (webs[[env]][i,c] > webs[[env]][j,c] && webs[[env]][j,c] != 0){
            k_ij = k_ij+1
          }
        }
        N_paired = k_ij/N_j
      }
      WNODFr_temp= WNODFr_temp + N_paired
    }
  }
  WNODFr[[env]] = WNODFr_temp*100
  # WNODF
  m = dim(webs[[env]])[1]
  n = dim(webs[[env]])[2]
  WNODF_obs[[env]] = 2*(WNODFc[[env]]+WNODFr[[env]])/ # The final average of WNODF for the env matrix
    (m*(m-1)+(n*(n-1)))
}
  
library(dplyr)
library(tidyverse)
how_much = factorial(20)/(factorial(10)*factorial(10))

# Produce one big fat matrix
data = read.csv("Data_filtered.csv", header=T, dec= ",", row.names=NULL)

permutations = 1000 # number of permutations
trilhas = sort(unique(data$"Trilha")) # unique sample points

###-----Frequency wNODF---------------------------------------------------------
rep_freq = 0
for (n in permutations) {
  rep_freq = rep_freq + 1
  trail_samp = sample(trilhas, 10)
  
  # Produce an horizonral version of the matrix
  data_filt = dplyr::select(data, c(1,3,6,7,11)) %>%  #selects relevant columns
    group_by(Espécie, Classificação.2) %>% #groups by codigo and item
    dplyr::summarise(Quantidade = sum(Quantidade)) %>% #sum quantity of each item
    pivot_wider(names_from = Classificação.2,
              values_from = Quantidade) %>%  
    as.data.frame() %>% #changes dataset format to Codigo-row oriented
    replace(is.na(.), 0) #removes NA and replaces with 0
  
  row.names(data_filt) = data_filt$Espécie
  data_filt = dplyr:: select(data_filt, !Espécie)
  data_filt = data_filt[,!names(data_filt) 
                            %in% c("Não Identificado", "Semente", "Semente de capim", "Apenas substrato", "Vazio")]

  # Filter the Matrice
  cols_to_remove = which(colSums(data_filt) < 5) # Min of 5 prey items
  if (length(cols_to_remove) > 0) { # Check if is a int(0)
    data_filt = data_filt[, -cols_to_remove] # Remove columns
  }
  rows_to_remove = which(rowSums(data_filt) == 0) # Check if there are items without a predator
  if (length(rows_to_remove) > 0) { # Check if is a int(0)
    data_filt = data_filt[-rows_to_remove,] # Remove columns
  }
  
  # Re-order the interaction matrices
  row_order = order(rowSums(data_filt), decreasing = TRUE) # Calculate row sums
  data_filt = data_filt[row_order, ] # Use the index vector to reorder the rows of the matrix
  col_order = order(colSums(data_filt), decreasing = TRUE) # Calculate col sums
  data_filt = data_filt[,col_order] # Use the index vector to reorder the cols of the matrix
  
}
###-----Volume wNODF------------------------------------------------------------
##------p_value AND ggplot------------------------------------------------------