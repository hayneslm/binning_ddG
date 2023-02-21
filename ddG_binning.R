install.packages("tidyr")

library(dplyr)
library(readr)
library(tidyr)

nbins = 100

ddg_eq02 = read_tsv("3q02_EvoEF_BM_ddG.xls", 
                    col_names = c("Pos","AA","Mut","ddG"))
ddg_eq02_bins = ddg_eq02 %>% mutate(bin = cut(ddG, breaks = nbins-1))
ddg_eq02_parsed = NULL
for (i in c(1:nbins)) {
  bin = i
  temp = ddg_eq02_bins %>% filter(bin == i) %>% nrow()
  temp2 = c(bin, temp)
  ddg_eq02_parsed = rbind(ddg_eq02_parsed, temp2)
}

  
  
