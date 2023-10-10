###
### Calculating and Plotting the eligibility of predators
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

# Filter 
rows_to_remove_eu = which(rowSums(webs$freq_eu) == 0) # Check if there are items without a predator
if (length(rows_to_remove_eu) > 0) { # Check if is a int(0)
  webs$vol_eu = webs$vol_eu[-rows_to_remove_eu,] # Remove lines
  webs$freq_eu = webs$freq_eu[-rows_to_remove_eu,] # Remove lines
}
rows_to_remove_mt = which(rowSums(webs$freq_mt) == 0) # Check if there are items without a predator
if (length(rows_to_remove_mt) > 0) { # Check if is a int(0)
  webs$vol_mt = webs$vol_mt[-rows_to_remove_mt,] # Remove lines
  webs$freq_mt = webs$freq_mt[-rows_to_remove_mt,] # Remove lines
}

# Re-order the interaction matrices
for (env in envs) {
  row_order = order(rowSums(webs[[env]]), decreasing = TRUE) # Calculate row sums
  webs[[env]] = webs[[env]][row_order, ] # Use the index vector to reorder the rows of the matrix
}

# Prey availability in the environment as the sum for all species
avail = list()
for (env in envs) {
  avail[[env]] = rowSums(webs[[env]])
}

# Prey proportion in the matrix per species
elec_matrix = list()
for (env in envs) { # Convert each matrix into a relative abundance matrix
  relative_web = sweep(webs[[env]], MARGIN = 1, STATS = avail[[env]], FUN = "/")
  col_sum = colSums(relative_web)
  elec_matrix[[env]] = sweep(relative_web, MARGIN = 2, STATS = col_sum, FUN = "/")
}
###-----------------------------------------------------------------------------
### Plot the Barplots in LAPTOP

# Plot the availability graphs
for (env in envs) {
  png(filename = paste("C:/Users/nunes/Documentos/LAB_VERT/Dieta_MV/Results/raw/eligibility/", 
                       env, "_avail.png", sep = ""), width = 500, height = 250)
  par(mar = c(2,6,3,1))
  barplot(avail[[env]], ylab = "Abundance of Prey \n (Frequency)", 
        ylim = c(0, max(avail[[env]])+50), names.arg = "", cex.axis = 1.5, lwd =2, cex.lab=1.5)
  bar_midpoints = barplot(avail[[env]], plot = FALSE)
  axis(side = 1, at = bar_midpoints, labels = FALSE, lwd = 2)
  dev.off()
}

# Plot the respective species graphs
for (env in envs) {
  png(filename = paste("C:/Users/nunes/Documentos/LAB_VERT/Dieta_MV/Results/raw/eligibility/", 
                       env, "_elig.png", sep = ""), width = 500, height = 3000)
  layout(c(1:dim(webs[[env]])[2]),12)
  for (j in 1:dim(webs[[env]])[2]) {
    par(mar = c(2,5,2,2))
    barplot(elec_matrix[[env]][,j], width = , ylim = c(0, max(elec_matrix[[env]][,j])+0.1), 
             names.arg = "", cex.axis = 1.5, lwd =2, axis.lty= 0)
    mtext(colnames(webs[[env]])[j], side = 3, line = 0, at = max(par("usr")[1], par("usr")[2]), cex = 1.5, font = 3, adj =1)
  
    bar_midpoints = barplot(elec_matrix[[env]][,j], plot = FALSE)
    axis(side = 1, at = bar_midpoints-1.2, labels = FALSE, lwd = 2)
    
    axis(side = 1, at = bar_midpoints+1.2, labels = FALSE, lwd = 2)
    }
  # mtext("Y-Axis Label", side = 2, line = 3, cex = 3)
  dev.off()
}


###-----------------------------------------------------------------------------
### Plot the Barplots in PC ###

# Plot the availability graphs
for (env in envs) {
  png(filename = paste("D:/Drive/Other computers/Meu laptop/LAB_VERT/Dieta_MV/Results/raw/eligibility/", 
                       env, "_avail.png", sep = ""), width = 600, height = 400)
  par(mar = c(4,6,2,2))
  barplot(avail[[env]], ylab = "Abundance of Prey \n (Frequency)", 
          ylim = c(0, max(avail[[env]])+50), names.arg = "", cex.axis = 1.5, lwd =2, cex.lab=1.5)
  bar_midpoints = barplot(avail[[env]], plot = FALSE)
  axis(side = 1, at = bar_midpoints-1.2, labels = FALSE, lwd = 2) # Bold axis
  
  axis(side = 1, at = bar_midpoints+1.2, labels = FALSE, lwd = 2) # Bold axis
  axis(side = 1, at = bar_midpoints, labels = names(avail[[env]]), cex.axis = 1.5, las = 2)
  dev.off()
}

# Plot the respective species graphs
for (env in envs) { # For each environment 
  max_spp = apply(elec_matrix[[env]][,1:dim(webs[[env]])[2]], # Find the maximum proportion for each species
                  MARGIN = 2, FUN = "max", simplify = TRUE)
  n = 0 # Counter object 
  for (sp in max_spp) { # Transform max proportions < 0.4 into 0.4
    n = n+1
    if (sp < 0.4) {
      max_spp[n] = 0.4
      }  
  }
  
  png(filename = paste("D:/Drive/Other computers/Meu laptop/LAB_VERT/Dieta_MV/Results/raw/eligibility/", # Open a png device for plotting
                       env, "_elig.png", sep = ""), width = 600, height = 300*sum(max_spp)) # The total height as a sum of the max_spp proportions
  layout(matrix(c(1:dim(webs[[env]])[2]),ncol=1),heights=max_spp) # Create the layout 
  
  for (j in 1:dim(webs[[env]])[2]) { # create the barplot for each species
    par(mar = c(4,5,0.5,2))
    barplot(elec_matrix[[env]][,j], width = , ylim = c(0, max(elec_matrix[[env]][,j])+0.1), 
            names.arg = "", cex.axis = 1.5, lwd =2, axis.lty= 0)
    mtext(colnames(webs[[env]])[j], side = 3, line = 0, at = max(par("usr")[1], par("usr")[2]), cex = 1.5, font = 3, adj =1)
    
    bar_midpoints = barplot(elec_matrix[[env]][,j], plot = FALSE)
    axis(side = 1, at = bar_midpoints-1.2, labels = FALSE, lwd = 2) # Bold axis
    
    axis(side = 1, at = bar_midpoints+1.2, labels = FALSE, lwd = 2) # Bold axis
  }
  axis(side = 1, at = bar_midpoints, labels = rownames(webs[[env]]), cex.axis = 1.5, las = 2)
  # mtext("Y-Axis Label", side = 2, line = 3, cex = 3)
  dev.off()
}
      