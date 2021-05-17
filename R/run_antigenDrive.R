
## Analyze 30 k sampled samples including 
## 
## i) CD8+ sorted T-LGLL from pb and bm; 
## ii) MNC sorted T-LGLL from multiple time points; 
## iii) CD8+ sorted RA
## iv) CD8+ sorted healthy
## v) CD8+ sorted melanoma
## 
## and vi) MNC sorted healthy
## and vii) sc-T-LGLL

## Get the gliph files; groups i-v here
filename <- c(list.files("src/bash/gliph2/30k/", pattern = "cluster.csv", full.names = T), list.files("src/bash/gliph2/emerson/", pattern = "cluster.csv", full.names = T), list.files("src/bash/gliph2/lgll_sc/", pattern = "cluster.csv", full.names = T))
gliph_df <- lapply(filename, FUN = function(x){message(x); readGliphFile(x) %>% filter(number_unique_cdr3 > 1)}) %>% rbindlist() %>% mutate(key = paste0(TcRb, Sample))

## Get the largest clones; define these by
## a ) taking the biggest clone
## b ) taking the > 5 % clone
filenames   <- c(list.files("results/manuscript/gliph2/30k/", full.names = T), list.files("results/manuscript/gliph2/emerson/", full.names = T), list.files("results/manuscript/gliph2/lgll_sc/", full.names = T))
input_files <- lapply(filenames, FUN = function(x){fread(x) %>% mutate(name = extractFileName(x))}) %>% rbindlist()
colnames(input_files) <- make.names(colnames(input_files))

biggest_clone <- input_files %>% group_by(name) %>% top_n(n=1, wt=count) %>% mutate(key = paste0(CDR3b, name))
over5_clone   <- input_files %>% filter(count>0.05) %>% mutate(key = paste0(CDR3b, name))

samples_with_clone <- over5_clone %>% pull(name) %>% unique()
sample_list         <- biggest_clone %>% pull(name) %>% unique()


## Get gliph with the largest clones
sample_df <- fread("results/manuscript/gliph2/sample_list.txt")

gliph_biggest_df <- gliph_df %>% filter(key %in% biggest_clone$key) %>% left_join(sample_df, by = c("Sample" = "file")) 
gliph_over_df    <- gliph_df %>% filter(key %in% over5_clone$key) %>% left_join(sample_df, by = c("Sample" = "file"))  %>% filter(Freq > 0.05)

gliph_biggest_df2 <- gliph_biggest_df %>% dplyr::select(pattern:cluster_size_score,TcRb:Freq,type.y,sort,key.y)
colnames(gliph_biggest_df2)[c(17,19)] <- c("sample type", "sample key")
fwrite(gliph_biggest_df2, "results/manuscript/gliph2/gliph_biggest.txt")

antigen_driven_biggest_sample <- gliph_biggest_df %>% pull(Sample) %>% unique()
antigen_driven_over_sample    <- gliph_over_df %>% pull(Sample) %>% unique()

sample_df <- sample_df %>% mutate(antigen_driven_biggest_sample = file %in% antigen_driven_biggest_sample, 
                                  over_sample = file %in% samples_with_clone)
sample_df$antigen_driven_over_sample    <- sample_df$file %in% antigen_driven_over_sample
sample_df$antigen_driven_biggest_sample <- sample_df$file %in% antigen_driven_biggest_sample

sample_df %>% 
  group_by(key, sort, antigen_driven_biggest_sample) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) %>% filter(antigen_driven_biggest_sample) %>% 
  
  ggplot(aes(reorder(key,prop), prop)) + geom_bar(stat = "identity",fill = "darkslategray") + geom_hline(yintercept = 0.5, linetype = "dotted") + ylim(c(0,1)) + facet_grid(~sort, scales = "free_x", space = "free_x") +
  ggpubr::rotate_x_text(45) + labs(x = "", y = "prop of \nantigen driven samples") + facets_nice
ggsave("results/manuscript/gliph2/bar_antigen_drive.pdf", width = 5, height = 5)

