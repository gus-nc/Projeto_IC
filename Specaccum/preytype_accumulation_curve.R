##-----------------------------------------------------------------------------
### Crate Prey Types acumulation curves for 
### each environment
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Augusto's Windows Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----------------------------------------------------------------------------
### Produce individual matrices for each species with >= 5 items
library(dplyr)
library(tidyverse)

# Produce one big fat matrix
data_filt = read.csv("Data_filtered.csv", header=T, dec= ",", row.names=NULL)

# Produce an horizonral version of the matrix
data_filt = dplyr::select(data_filt, 1:7) %>%  #selects relevant columns
  group_by(Código, Ambiente, Classificação.2) %>% #groups by codigo and item
  dplyr::summarise(Classificação.2, Quantidade = sum(Quantidade)) %>% #sum quantity of each item
  distinct() %>% #only unique values
  pivot_wider(names_from = Classificação.2,
              values_from = Quantidade) %>%  as.data.frame() %>% #changes dataset format to Codigo-row oriented
  replace(is.na(.), 0) #removes NA and replaces with 0

# Remove unwanted columns
data_filt = data_filt[,!names(data_filt) 
                                        %in% c("Não Identificado", "Semente", "Semente de capim", "Apenas substrato", "Vazio")]

# Create 1 matrice for each environment
env_matrices = list() # Create a list for each environment matrix

envs = c("mata", "eucalipto", "total")
for (env in envs) {
  if (env != "total") {
    env_matrices[[env]] = data_filt[data_filt[,2] == env,]
  } else {
    env_matrices[[env]] = data_filt
  }
  row.names(env_matrices[[env]]) = env_matrices[[env]]$Código
  env_matrices[[env]] = env_matrices[[env]][,-c(1,2)]
  cols_to_remove = which(colSums(env_matrices[[env]]) == 0) # Identify columns with sum == 0
  if (length(cols_to_remove) > 0) {
    env_matrices[[env]]= env_matrices[[env]][, -cols_to_remove] # Remove columns
  }
}
###-----------------------------------------------------------------------------
### Produce a specaccum curve for each environment
library("vegan")

### Create individual curves for each environment
curves = list()
for (env in envs) {
  curves[[env]] = specaccum(env_matrices[[env]], method = "exact")
}

# Extract the relevant data from the specaccum object
x = list()
y = list()
sd = list()
for (env in envs) {
  x[[env]] = c(1:nrow(env_matrices[[env]])) # Retrieves the number of individuals
  y[[env]] = curves[[env]]$richness # Extract the expected richness at a given number of sites
  sd[[env]] = curves[[env]]$sd # Extract the expected sd at a given number of sites
  
}

cl = c("black", "darkgray", "gray")

plot(curves[["total"]], ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
      ci.lty=0, xlab = "Individuals", ylab = "Prey Types",
      col = cl[1], 
      xlim = c(1,250), ylim=c(0,30))
par(new = TRUE)

plot(curves[["eucalipto"]], ci.type = "line", lty=2, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Individuals", ylab = "Prey Types",
     col = cl[2], 
     xlim = c(1,250), ylim=c(0,30))
par(new = TRUE)

plot(curves[["mata"]], ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Individuals", ylab = "Prey Types",
     col = cl[3], 
     xlim = c(1,250), ylim=c(0,30))
par(new = TRUE)

title(main = "Prey Types Accumulation Curve") # Add title 
legend("bottomright", legend = c("Total", "Eucalyptus", "Atlântic Forest"),
       col= cl, lty = c(1,2,1),
       lwd = 2, bty = "n") # Add lines legend

#-------------------------------------------------------------------------------