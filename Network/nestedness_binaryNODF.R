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

###-----  
library(dplyr)
library(tidyverse)
factorial(20)/(factorial(10)*factorial(10)) # How much possible permutations

# Produce one big fat matrix
data = read.csv("Data_filtered.csv", header=T, dec= ",", row.names=NULL)

permutations = 10000 # number of permutations
trilhas = sort(unique(data$"Trilha")) # unique sample points
###-----Permutation for NODF---------------------------------------------------------
NODF_freq = c()
for (n in 1:permutations) {
  trail_samp = sample(trilhas, 10)
  
  # Produce an horizonral version of the matrix
  data_filt = filter(data, Trilha %in% trail_samp) %>%
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
  
  # Make a Binary Matrix
  data_filt = ifelse(as.matrix(data_filt) > 0, 1, as.matrix(data_filt))
  
  # Re-order the interaction matrices
  row_order = order(rowSums(data_filt), decreasing = TRUE) # Calculate row sums
  data_filt = data_filt[row_order, ] # Use the index vector to reorder the rows of the matrix
  col_order = order(colSums(data_filt), decreasing = TRUE) # Calculate col sums
  data_filt = data_filt[,col_order] # Use the index vector to reorder the cols of the matrix
  
  ### Calculate NODF Frequency permutation networks
  
  # NODFc values for cols
  NODFc_temp = 0 # Create a temporary object for NODFc
  for (i in 1:(ncol(data_filt)-1)) { # for every i:n-1
    for (j in (i+1):ncol(data_filt)) { # and for every i+1:n
      if(sum(data_filt[,i]) <= sum(data_filt[,j])) { # if F(ci) <= F(cj)
        N_paired = 0        
      }
      else {
        PO = 0 # number of 1's in ci = cj, if cj != 0?
        N_j = sum(data_filt[,j] == 1) # Cells > 0 in cj
        for (r in 1:nrow(data_filt)) { # for every row in t cols sum PO
          if (data_filt[r,i] == 1 && data_filt[r,j] == 1){
            PO = PO+1            }
        }
        N_paired = PO/N_j
      }
      NODFc_temp = NODFc_temp + N_paired 
    } 
  }
  NODFc = NODFc_temp*100
  # NODFr 
  NODFr_temp = 0
  for (i in 1:(nrow(data_filt)-1)) {
    for (j in (i+1):nrow(data_filt)) {
      if(sum(data_filt[i,]) <= sum(data_filt[j,])) {
        N_paired = 0
      }
      else {
        PO = 0
        N_j = sum(data_filt[j,] == 1)
        for (c in 1:ncol(data_filt)) {
          if (data_filt[i,c] == 1 && data_filt[j,c] == 1){
            PO = PO+1            }
        }
        N_paired = PO/N_j
      }
      NODFr_temp= NODFr_temp + N_paired
    }
  }
  NODFr = NODFr_temp*100
  # NODF
  m = dim(data_filt)[1]
  n = dim(data_filt)[2]
  NODF_freq = append(NODF_freq, 2*(NODFc+NODFr)/(m*(m-1)+(n*(n-1)))) 
}
###------P values and Histograms------------------------------------------------
p_value = list()
p_value[["freq_mt"]] = sum(NODF_freq >= NODF_obs[["freq_mt"]]) / permutations
p_value[["freq_eu"]] = sum(NODF_freq >= NODF_obs[["freq_eu"]]) / permutations

# Plot the p_values
f_path = "D:/Drive/Other computers/Meu laptop/LAB_VERT/Dieta_MV/Results/raw/NODF_permutation_p_values.txt"

sink(f_path)
cat("Permutations:", permutations, "\n")
for (env in envs) {
  cat(env, "\n")
  cat("Observed Test Statistic:", NODF_obs[[env]], "\n")
  cat("P-Value:", p_value[[env]], "\n\n")
}
sink() # Close the sink to stop redirecting output


library(ggplot2)
# png(filename = "C:/Users/nunes/Documentos/LAB_VERT/Dieta_MV/Results/raw/perm_freq_NODF.png"
#, width = 1000, height = 800) # LAPTOP
png(filename = "D:/Drive/Other computers/Meu laptop/LAB_VERT/Dieta_MV/Results/raw/perm_binary_NODF.png"
    , width = 1000, height = 800) # PC
# Data for the histogram
plot_data <- data.frame(Statistic = NODF_freq)

# Create histogram for Frequency Permutation
ggplot(plot_data, aes(x = Statistic)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black") +
  geom_vline(xintercept = NODF_obs[["freq_eu"]], color = "red", linetype = "dashed") +
  geom_vline(xintercept = NODF_obs[["freq_mt"]], color = "red", linetype = "dashed") +
  labs(title = "NODF Distribution",
       x = "NODF",
       y = "Frequency") +
  theme_minimal()
dev.off()
