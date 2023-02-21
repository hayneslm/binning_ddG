nbins = 50

data_0h = read_tsv("filtered_missense_0h.xls")
data_0h_func = data_0h %>% filter(log2FoldChange > 0)
data_0h_nf = data_0h %>%  filter(log2FoldChange <= 0)

AminL_func = full_join(data_0h_func, AminL_cut) %>% na.omit() %>% select(bin)
AminL_func_table = table(AminL_func) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(AminL_func_table, "AminL_func.xls")
plot(AminL_func_table)


AminL_nf = full_join(data_0h_nf, AminL_cut) %>% na.omit() %>% select(bin)
AminL_nf_table = table(AminL_nf) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(AminL_nf_table, "AminL_nf.xls")
plot(AminL_nf_table)

############################
all_48h = read_tsv("filtered_missense_48h.xls")
func_0h_list = data_0h_func %>% select(mutation)
func48 = left_join(all_48h, func_0h_list) %>% na.omit()

func_48h_active = full_join(AminL_cut, func48) %>% 
  filter(log2FoldChange > 0) %>% na.omit() %>% 
  select(bin)
table48h_active = table(func_48h_active) %>% as.data.frame()%>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(table48h_active, "AmL_func_active.xls")
plot(table48h_active)

func_48h_latent = full_join(AminL_cut, func48) %>% 
  filter(log2FoldChange <= 0) %>% na.omit() %>% 
  select(bin)
table48h_latent = table(func_48h_latent) %>% as.data.frame()%>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(table48h_latent, "AmL_func_latent.xls")
plot(table48h_latent)
