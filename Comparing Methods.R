library(dplyr)

# Import and format all files
data_0h = read_tsv("filtered_missense_0h.xls")

ddg_eq02 = read_tsv("3q02_EvoEF_BM_ddG.xls", 
                    col_names = c("Pos","AA","Mut","ddG")) %>% 
  mutate(mutation = paste(AA,Pos,Mut, sep = ""))

EVE_score = read_csv("PAI1_human_eve.csv") %>% filter(mutations != "wt") %>% 
  separate(mutations, into = c("AA","Rest"), sep = "(?<=[A-Z])(?=[0-9])") %>%   
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(Pos = as.double(Pos)) %>% 
  mutate(Pos = Pos-23) %>% 
  mutate(mutation = paste(AA,Pos,Mut, sep = "")) %>% filter(Pos>0)

AM_score = read_tsv("PAI1_alphamissense.txt", col_names = F) %>% 
  select(!X1) %>% rename(mutations = X2, AMscore = X3, type = X4) %>% 
  separate(mutations, into = c("AA","Rest"), sep = "(?<=[A-Z])(?=[0-9])") %>%   
  separate(Rest, into= c("Pos","Mut"), sep = "(?<=[0-9])(?=[A-Z])") %>% 
  mutate(Pos = as.double(Pos)) %>% 
  mutate(Pos = Pos-23) %>% 
  mutate(mutation = paste(AA,Pos,Mut, sep = "")) %>% filter(Pos>0)
 

Pooled = full_join(data_0h, ddg_eq02, by = c("AA","Pos","Mut","mutation")) %>% 
  full_join(., EVE_score, by = c("AA","Pos","Mut","mutation")) %>% 
  full_join(., AM_score, by = c("AA","Pos","Mut","mutation")) %>% na.omit()

