install.packages("tidyr")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# calculate bins based on ddG data
nbins = 50

ddg_eq02 = read_tsv("3q02_EvoEF_BM_ddG.xls", 
                    col_names = c("Pos","AA","Mut","ddG"))
ddg_eq02_cut = ddg_eq02 %>% mutate(bin = cut(ddG, breaks = nbins-1))
ddg_eq02_parsed = table(ddg_eq02_cut$bin) %>% as.data.frame()
write_tsv(ddg_eq02_parsed,"eq02_all.xls")

plot(ddg_eq02_parsed)

ddg_dvn1 = read_tsv("1dvn_EvoEF_BM_ddG.xls",
                    col_names = c("Pos","AA","Mut","ddG"))
ddg_dvn1_cut = ddg_dvn1 %>% mutate(bin = cut(ddG, breaks = nbins-1))
ddg_dvn1_parsed = table(ddg_dvn1_cut$bin) %>% as.data.frame()
write_tsv(ddg_dvn1_parsed,"dvn1_all.xls")

plot(ddg_dvn1_parsed)

########################
# Acive - Latent ddG
AminL = full_join(ddg_eq02,ddg_dvn1, by = c("Pos", "AA", "Mut")) %>% 
  na.omit() %>% mutate(AmL = ddG.x - ddG.y)
AminL_cut = AminL %>% mutate(bin = cut(AmL, breaks = nbins-1))
AminL_table = table(AminL_cut$bin) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(AminL_table,"AminL_all.xls")
plot(AminL_table)

#Latent-Active ddG
LminA = AminL %>% rename(LmA = AmL) %>% mutate(LmA = ddG.y-ddG.x)
LminA_cut = LminA %>% mutate(bin = cut(LmA, breaks = nbins-1))
LminA_table = table(LminA_cut$bin) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(LminA_table,"LminA_all.xls")
plot(LminA_table)

# Active plus Latent ddG
AplusL = full_join(ddg_eq02,ddg_dvn1, by = c("Pos", "AA", "Mut")) %>% 
  na.omit() %>% mutate(ApL = ddG.x + ddG.y)
AplusL_cut = AplusL %>% mutate(bin = cut(ApL, breaks = nbins-1))
AplusL_table = table(AplusL_cut$bin) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(AplusL_table,"AplusL_all.xls")
plot(AplusL_table)
