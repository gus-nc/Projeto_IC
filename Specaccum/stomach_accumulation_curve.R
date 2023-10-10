###-----------------------------------------------------------------------------
### Crate stomach content acumulation curves for 
### each species in each environment
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Augusto's Windows PC
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----------------------------------------------------------------------------
### Produce individual matrices for each species with >= 5 items
library(dplyr)
library(tidyverse)

# Produce one big fat matrix
data_filt = read.csv("Data_filtered.csv", header=T, dec= ",", row.names=NULL)

# Detect unique species
species = sort(unique(data_filt$"Espécie"))
# Remove species with less than 5 items
item_sum = sapply(species, function(sp) { # Create the sum for each species
  filt = list(); # Create a list to append the data.frames for each specie
  filt = subset(data_filt, Espécie == sp); # Subset the data.frame based on column "Espécie"
  sum(filt$Quantidade) # Sum the "Quanridaede"column
}, simplify=FALSE, USE.NAMES=TRUE) # Prevent simplify and force USE.NAMES

species_to_keep = names(which(item_sum >= 5)) # Define species with >= 5 items


# Produce an horizonral version of the matrix
data_filt = dplyr::select(data_filt, 1:7) %>%  #selects relevant columns
  group_by(Código, Classificação.2) %>% #groups by codigo and item
  dplyr::summarise(Espécie, Classificação.2, Quantidade = sum(Quantidade)) %>% #sum quantity of each item
  distinct() %>% #only unique values
  pivot_wider(names_from = Classificação.2,
              values_from = Quantidade) %>%  as.data.frame() %>% #changes dataset format to Codigo-row oriented
  replace(is.na(.), 0) #removes NA and replaces with 0

data_filt_sp = list()
for (sp in species_to_keep) { # to generate a species-specific "matrix"
  data_filt_sp[[sp]] = data_filt %>% 
  filter(Espécie == sp) %>% 
  dplyr::select(-2)
  row.names(data_filt_sp[[sp]]) = data_filt_sp[[sp]]$Código
  data_filt_sp[[sp]] = dplyr:: select(data_filt_sp[[sp]], !Código)
  data_filt_sp[[sp]] = data_filt_sp[[sp]][,!names(data_filt_sp[[sp]]) 
              %in% c("Não Identificado", "Semente", "Semente de capim", "Apenas substrato", "Vazio")]
  cols_to_remove = which(colSums(data_filt_sp[[sp]]) == 0) # Identify columns with sum == 0
  print.default(cols_to_remove)
  if (length(cols_to_remove) > 0) {
  data_filt_sp[[sp]] = data_filt_sp[[sp]][, -cols_to_remove] # Remove columns
  }
  rows_to_remove = which(rowSums(data_filt_sp[[sp]]) == 0) # Identify rows with sum == 0
  if (length(rows_to_remove) > 0) {
  data_filt_sp[[sp]] = data_filt_sp[[sp]][-rows_to_remove, ] # Remove rows
  }
  #print(head(data_filt_sp))
  #write.csv(data_filt_sp, paste0("data_filt_", sp, ".csv"))
}

###-----------------------------------------------------------------------------
### Produce each species specaccum curve
library("vegan")

### Create individual curves for every species
curves = list()
for (sp in species_to_keep) {
  curves[[sp]] = specaccum(data_filt_sp[[sp]], method = "exact")
}

# Extract the relevant data from the specaccum object
x = list()
y = list()
sd = list()
for (sp in species_to_keep) {
  x[[sp]] = c(1:nrow(data_filt_sp[[sp]])) # Retrieves the number of individuals
  y[[sp]] = curves[[sp]]$richness # Extract the expected richness at a given number of sites
  sd[[sp]] = curves[[sp]]$sd # Extract the expected sd at a given number of sites
  
}

for (sp in species_to_keep) {
  plot(curves[[sp]], ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
       ci.lty=0, xlab = "Individuals", ylab = "Prey Types", 
       xlim = c(1,20), ylim=c(0,15))
  #points(x[[sp]], y[[sp]], pch = 16, col = "black", # Plot the richness values 
        # cex = 1, xlim = c(1:nrow(data_filt_sp[[sp]])), ylim=c(0,30))
  par(new = TRUE)
}

title(main = "Prey Items Accumulation Curve") # Add title 
#legend("bottomright", legend = c("Total", "Eucalyptus", "Atlântic Forest"),
     #  col= c("black", "gray", "gray"),
     #  lwd = 2, bty = "n") # Add lines legend

#legend("bottomright", legend = c("Total", "Eucalyptus", "Atlântic Forest"), 
      # lwd = 0, bty = "n", pch = c(16, 15, 17), 
      # title = "Habitats") # Add rest of the legend
png("total_items_curves.png",  width = 600, height = 600)
#-------------------------------------------------------------------------------