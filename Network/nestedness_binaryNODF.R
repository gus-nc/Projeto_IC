###
### Calculating Nestedness for the networks according to NODF Almeida 2008
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

#Create a Binary Network
for (env in envs) {
  webs[[env]] = ifelse(webs[[env]] > 0, 1, webs[[env]])
}

# Re-order the interaction matrices
for (env in envs) {
  row_order = order(rowSums(webs[[env]]), decreasing = TRUE) # Calculate row sums
  webs[[env]] = webs[[env]][row_order, ] # Use the index vector to reorder the rows of the matrix
  col_order = order(colSums(webs[[env]]), decreasing = TRUE) # Calculate col sums
  webs[[env]] = webs[[env]][,col_order] # Use the index vector to reorder the cols of the matrix
}
###-----Export the matrices-----------------------------------------------------
for (env in envs) {
  write.table(webs[[env]], paste( "Network/Aninhado_3/", env, ".csv", sep = ""), 
              sep = ",", eol = "\n", col.names = FALSE, row.names = FALSE)
}
###-----Observed NODF-----------------------------------------------------------
### Calculate the brute value of NODF for each env empirical network
NODFc = list() # List for storing the sum values for cols
NODFr = list() # List for storing the sum values for rowws
NODF_obs = list() # List for storing the final value for each matrix

for (env in envs) { # For each environment
  # NODFc values for cols
  NODFc_temp = 0 # Create a temporary object for NODFc
  for (i in 1:(ncol(webs[[env]])-1)) { # for every col i:n-1
    for (j in (i+1):ncol(webs[[env]])) { # and for every col i+1:n
      if(sum(webs[[env]][,i]) <= sum(webs[[env]][,j])) { # if F(ci) <= F(cj)
        N_paired = 0 # The Nested Value for this pair of cols
      }
      else {
        PO = 0 # number of 1's in ci = cj, if cj != 0?
        N_j = sum(webs[[env]][,j] == 1) # Cells > 0 in cj
        for (r in 1:nrow(webs[[env]])) { # for every row in t cols sum PO
          if (webs[[env]][r,i] == 1 && webs[[env]][r,j] == 1){
            PO = PO+1
          }
        }
        N_paired = PO/N_j # The Nested Value for this pair of cols
      }
      NODFc_temp = NODFc_temp + N_paired 
    }
  }
  NODFc[[env]] = NODFc_temp*100
  # NODFr values for rows
  NODFr_temp = 0 # Create a temporary object for NODFr
  for (i in 1:(nrow(webs[[env]])-1)) { # for every col i:n-1
    for (j in (i+1):nrow(webs[[env]])) { # and for every col i+1:n
      if(sum(webs[[env]][i,]) <= sum(webs[[env]][j,])) { # if F(ci) <= F(cj)
        N_paired = 0 # The Nested Value for this pair of cols
      }
      else {
        PO = 0 # number of 1's in ci = cj, if cj != 0?
        N_j = sum(webs[[env]][j,] == 1)  #Cells > 0 in cj
        for (c in 1:ncol(webs[[env]])) { # for every row in t cols sum PO
          if (webs[[env]][i,c] == 1 && webs[[env]][j,c] == 1){
            PO = PO+1
          }
        }
        N_paired = PO/N_j # The Nested Value for this pair of cols
      }
      NODFr_temp= NODFr_temp + N_paired
    }
  }
  NODFr[[env]] = NODFr_temp*100
  # NODF
  m = dim(webs[[env]])[1] # Number of rows
  n = dim(webs[[env]])[2] # Number of cols
  NODF_obs[[env]] = 2*(NODFc[[env]]+NODFr[[env]])/ # The final average of NODF for the env matrix
    (m*(m-1)+(n*(n-1)))
}
