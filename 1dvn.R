nbins = 50

data_0h = read_tsv("filtered_missense_0h.xls")
data_0h_func = data_0h %>% filter(log2FoldChange > 0)
data_0h_nf = data_0h %>%  filter(log2FoldChange <= 0)

dvn1_func = full_join(data_0h_func, ddg_dvn1_cut) %>% na.omit() %>% select(bin)
dvn1_func_table = table(dvn1_func) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(dvn1_func_table, "dvn1_func.xls")
# plot(dvn1_func_table)


dvn1_nf = full_join(data_0h_nf, ddg_dvn1_cut) %>% na.omit() %>% select(bin)
dvn1_nf_table = table(dvn1_nf) %>% as.data.frame() %>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(dvn1_nf_table, "dvn1_nf.xls")
# plot(dvn1_nf_table)

############################
all_48h = read_tsv("filtered_missense_48h.xls")
func_0h_list = data_0h_func %>% select(mutation)
func48 = left_join(all_48h, func_0h_list) %>% na.omit()

func_48h_active = full_join(ddg_dvn1_cut, func48) %>% 
  filter(log2FoldChange > 0) %>% na.omit() %>% 
  select(bin)
table48h_active = table(func_48h_active) %>% as.data.frame()%>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(table48h_active, "dvn1_func_active.xls")
plot(table_48h)

func_48h_latent = full_join(ddg_dvn1_cut, func48) %>% 
  filter(log2FoldChange <= 0) %>% na.omit() %>% 
  select(bin)
table48h_latent = table(func_48h_latent) %>% as.data.frame()%>% 
  mutate(fraction = Freq/sum(Freq))
write_tsv(table48h_latent, "dvn1_func_latent.xls")
plot(table_48h)

