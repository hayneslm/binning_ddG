nbins = 50

data_0h = read_tsv("filtered_missense_0h.xls")
data_0h_func = data_0h %>% filter(log2FoldChange > 0)
data_0h_nf = data_0h %>%  filter(log2FoldChange <= 0)

eq02_func = full_join(data_0h_func, ddg_eq02_cut) %>% na.omit() %>% select(bin)
eq02_func_table = table(eq02_func)
plot(eq02_func_table)


eq02_nf = full_join(data_0h_nf, ddg_eq02_cut) %>% na.omit() %>% select(bin)
eq02_nf_table = table(eq02_nf)
plot(eq02_nf_table)

############################
all_48h = read_tsv("filtered_missense_48h.xls")
func_0h_list = data_0h_func %>% select(mutation)
func48 = left_join(all_48h, func_0h_list) %>% na.omit()

func_48h_parsed2 = full_join(ddg_eq02_cut, func48) %>% na.omit() %>% 
  select(bin)
table48h = table(func_48h_parsed2) %>% as.data.frame()
plot(table_48h)




nf_0h_list = data_0h_nf %>% select(mutation)
nf48 = left_join(all_48h, nf_0h_list) %>% na.omit()

nf_48h_parsed2 = full_join(ddg_eq02_cut, nf48) %>% na.omit() %>% 
  select(bin)
table48h_nf = table(nf_48h_parsed2) %>% as.data.frame()
plot(table48h_nf)



