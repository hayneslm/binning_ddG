nbins = 50

data_0h = read_tsv("filtered_missense_0h.xls")
data_0h_func = data_0h %>% filter(log2FoldChange > 0)
data_0h_nf = data_0h %>%  filter(log2FoldChange <= 0)

AplusL_func = full_join(data_0h_func, AplusL_cut) %>% na.omit() %>% select(bin)
AplusL_func_table = table(AplusL_func) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(AplusL_func_table, "AplusL_func.xls")
plot(AplusL_func_table)


AplusL_nf = full_join(data_0h_nf, AplusL_cut) %>% na.omit() %>% select(bin)
AplusL_nf_table = table(AplusL_nf) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(AplusL_nf_table, "AplusL_nf.xls")
plot(AplusL_nf_table)

############################
all_48h = read_tsv("filtered_missense_48h.xls")
func_0h_list = data_0h_func %>% select(mutation)
func48 = left_join(all_48h, func_0h_list) %>% na.omit()

func_48h_active = full_join(AplusL_cut, func48) %>% 
  filter(log2FoldChange > 0) %>% na.omit() %>% 
  select(bin)
table48h_active = table(func_48h_active) %>% as.data.frame()%>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(table48h_active, "ApL_func_active.xls")
plot(table48h_active)

func_48h_latent = full_join(AplusL_cut, func48) %>% 
  filter(log2FoldChange <= 0) %>% na.omit() %>% 
  select(bin)
table48h_latent = table(func_48h_latent) %>% as.data.frame()%>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(table48h_latent, "ApL_func_latent.xls")
plot(table48h_latent)
