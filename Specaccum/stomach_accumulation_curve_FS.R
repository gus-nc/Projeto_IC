###-----------------------------------------------------------------------------
### Crate stomach content acumulation curves for 
### each species in each environment
###-----------------------------------------------------------------------------
setwd("C:/Users/nunes/Documentos/Code/Projeto_IC") # Augusto's Windows PC
setwd("/Volumes/GoogleDrive/Other computers/Meu laptop/Code/Projeto_IC/") # Augusto's MacOS Laptop
###-----------------------------------------------------------------------------
library(dplyr)
Sys.getlocale(category="LC_ALL")

# Produce one big fat matrix
data_filt = read.csv("Data_filtered.csv", header=T, dec= ",", row.names=NULL) %>% # reads dataset 
  dplyr::select(1:7) %>%  #selects relevant columns
  group_by(Código, Classificação.2) %>% #groups by codigo and item
  dplyr::summarise(Espécie, Classificação.2, Quantidade = sum(Quantidade)) %>% #sum quantity of each item
  distinct() %>% #only unique values
  pivot_wider(names_from = Classificação.2,
              values_from = Quantidade) %>%  as.data.frame() %>% #changes dataset format to Codigo-row oriented
  replace(is.na(.), 0) #removes NA and replaces with 0
  
## to generate a species-specific "matrix"
data_filt_sp = data_filt %>%
  filter(Espécie == "Aplastodiscus albosignatus")

write.csv(data_filt, "data_filt_FS.csv")
