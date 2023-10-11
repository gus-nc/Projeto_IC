###
### Calculating the niche breadth predators
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Laptop
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
###-----------------------------------------------------------------------------
#Read all the interaction network matrices into a list
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
  webs$vol_eu = webs$vol_eu[-rows_to_remove_eu,] # Remove lines
  webs$freq_eu = webs$freq_eu[-rows_to_remove_eu,] # Remove lines
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
  webs$vol_mt = webs$vol_mt[-rows_to_remove_mt,] # Remove lines
  webs$freq_mt = webs$freq_mt[-rows_to_remove_mt,] # Remove lines
}




###-----------------------------------------------------------------------------
### Retrieve species that occur in both environments
species_x = Reduce(intersect, list(colnames(webs[["vol_eu"]]), colnames(webs[["vol_mt"]]))) # Identify common elements in the columns

# Prey availability in the environment as the sum for all species
sum_spp = list()
for (env in envs) {
  sum_spp[[env]] = colSums(webs[[env]])
}
# Prey proportion in the matrix per species
relative_webs = list()
for (env in envs) { # Convert each matrix into a relative abundance matrix
  relative_webs[[env]] = sweep(webs[[env]], MARGIN = 2, STATS = sum_spp[[env]], FUN = "/")
}

# Calculate Niche Breadth Pianka (1986) for every env and spp
nb_spp = list() # create a list for the NB values for each env
for (env in envs) { # for each env
  for (j in 1:ncol(webs[[env]])) # look for each species
   nb_spp[[env]] = append(nb_spp[[env]], 1/(sum(sapply(relative_webs[[env]][,j],
                                        function(x) x^2)))) # append the value and spp name
  names(nb_spp[[env]]) = colnames(webs[[env]])
}

###-----------------------------------------------------------------------------
###Producing the table for all species and all environments
spp_total = sort(unique(c(colnames(webs[["vol_mt"]]),colnames(webs[["vol_eu"]])))) # all species in the system
df = data.frame(matrix(nrow = length(spp_total), ncol = 4)) # Create an empty df
rownames(df) =  spp_total # all species ass rows
colnames(df) = envs # environments as columns

for (env in envs) { # for every environment
  for (i in 1:nrow(df)) { # look for each species
    if (rownames(df[i,]) %in% names(nb_spp[[env]])) { # see if the species is in the env
      df[i,env] = as.numeric(nb_spp[[env]][rownames(df[i,])])
      }
    else {
      df[i,env] = NA
      }
    }
}

print(df)

# Calculate the mean NB value for each environment
mean_nb_spp = c()
for (env in envs) {
  mean_nb_spp[env] = mean(df[[env]], na.rm = TRUE)
}
mean_nb_spp
###-----------------------------------------------------------------------------
### Plotting and saving as an image
library(gt)
library(gtExtras)

### View formated data.frame
tab_1 = df |> 
  gt(rownames_to_stub = TRUE) |> 
  tab_header(
    title = "Niche Breadth",
    subtitle = "Co-occuring species in Mata Atlântica and Eucalyptus"
  ) |> 
  fmt_number(columns = colnames(df), decimals = 4) |> 
  tab_style(style = cell_text(align = "center", v_align = "middle"), locations = cells_body()) |> 
  sub_missing( columns = everything(), rows = everything(),
    missing_text = "-") 
tab_1

#Save as PNG
gtsave(tab_1, filename = "Niche_Breadth.png", path = "D:/Drive/Other computers/Meu laptop/LAB_VERT/Dieta_MV/Results/")