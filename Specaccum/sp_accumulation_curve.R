###-----------------------------------------------------------------------------
### Create species acumulation curves
###-----------------------------------------------------------------------------

setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Augusto's Windows PC
setwd("D:/Drive/Other computers/Meu laptop/Code/Projeto_IC") # PC
library("vegan")

# Import tables .csv files
spp_t = read.csv("Specacum_Spp_T.csv", header=T, row.names=1)
spp_eu = read.csv("Specacum_Spp_Eu.csv", header=T, row.names=1)
spp_mt = read.csv("Specacum_Spp_Mt.csv", header=T, row.names=1)

# Defining the colors
# blue50 = rgb(173, 216, 230, max=255, alpha=125)
# red50 = rgb(250, 182, 193, max=255, alpha=125)
# green50 = rgb(144, 238, 144, max=255, alpha=125)


### Create individually the curve data for all the sampling
curve_t <- specaccum(spp_t, method = "exact")

# Extract the relevant data from the specaccum object
x_t = c(1:20) # For the total sampling there are 20 sites
y_t = curve_t$richness # Extract the expected richness at a given number of sites
sd_t = curve_t$sd # Extract the expected sd at a given number of sites

plot(curve_t, ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Sites", ylab = "Species Richness", 
     xlim=c(1,20), ylim=c(0,30))

points(x_t, y_t, pch = 16, col = "black", # Plot the richness values 
       cex = 1, xlim=c(1,20), ylim=c(0,30))

# arrows(x,y-sd,x,y+sd, code=3, length=0.02, angle = 90)

### Create individually  the curve data for eucalyptus sampling
curve_eu <- specaccum(spp_eu, method = "exact")

# Extract the relevant data from the specaccum object
x_eu = c(1:10) # For the total sampling there are 10 sites
y_eu = curve_eu$richness # Extract the expected richness at a given number of sites
sd_eu = curve_eu$sd # Extract the expected sd at a given number of sites

plot(curve_eu, ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Sites", ylab = "Species Richness", 
     xlim=c(1,20), ylim=c(0,30))

points(x_eu, y_eu, pch = 15, col = "black", # Plot the richness values 
       cex = 1, xlim=c(1,20), ylim=c(0,30))

### Create individually  the curve data for Mata Atl?ntica sampling
curve_mt <- specaccum(spp_mt, method = "exact")

# Extract the relevant data from the specaccum object
x_mt = c(1:10) # For the total sampling there are 20 sites
y_mt = curve_mt$richness # Extract the expected richness at a given number of sites
sd_mt = curve_mt$sd # Extract the expected sd at a given number of sites

plot(curve_mt, ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Sites", ylab = "Species Richness", 
     xlim=c(1,20), ylim=c(0,30))

points(x_mt, y_mt, pch = 17, col = "black", # Plot the richness values 
       cex = 1, xlim=c(1,20), ylim=c(0,30))

### Plot all curves in the same graph
plot(curve_t, ci.type = "line", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Sites", ylab = "Species Richness", 
     xlim=c(1,20), ylim=c(0,30))
points(x_t, y_t, pch = 16, col = "black", # Plot the richness values 
       cex = 1, xlim=c(1,20), ylim=c(0,30))
par(new = TRUE)
plot(curve_eu, ci.type = "line", col="gray", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Sites", ylab = "Species Richness", 
     xlim=c(1,20), ylim=c(0,30))
points(x_eu, y_eu, pch = 15, col = "black", # Plot the richness values 
       cex = 1, xlim=c(1,20), ylim=c(0,30))
par(new = TRUE)
plot(curve_mt, ci.type = "line", col= "gray", lty=1, lwd=2, # plot the sample-based regression 
     ci.lty=0, xlab = "Sites", ylab = "Species Richness", 
     xlim=c(1,20), ylim=c(0,30))
points(x_mt, y_mt, pch = 17, col = "black", # Plot the richness values 
       cex = 1, xlim=c(1,20), ylim=c(0,30))

title(main = "Species Accumulation Curve") # Add title 
legend("bottomright", legend = c("Total", "Eucalyptus", "Atlântic Forest"),
       col= c("black", "gray", "gray"),
       lwd = 2, bty = "n") # Add lines legend

legend("bottomright", legend = c("Total", "Eucalyptus", "Atlântic Forest"), 
      lwd = 0, bty = "n", pch = c(16, 15, 17), 
      title = "Habitats") # Add rest of the legend

#-------------------------------------------------------------------------------

