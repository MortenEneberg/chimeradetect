

Rlibpath<- "/user_data/men/sepseq/R_lib/x86_64-pc-linux-gnu-library/3.6"
.libPaths(c(Rlibpath, .libPaths()))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(gggenes)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))

sample_names<-c("dummy114")

read_sam <- function(batch) {
  file<-paste(batch, ".sam", sep = "")
  sam_read<-read.csv2(file = file, header = F, sep = "\t")
colnames(sam_read)<- c("barcode", "flag", "read", "start_position", "mapq", "cigar", "drop1", "drop2", "drop3", "drop4", "drop5", "align_metrics") 
sam_read<-sam_read %>%    
  select(barcode, flag, read, start_position, mapq, cigar, align_metrics) %>%
    filter(flag != 4) %>%
    separate(align_metrics, into = c(NA, NA, "align_score", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), sep = ":") %>%
    separate(align_score, into = c("align_score", NA)) %>%
    mutate(batch = batch)%>%
    mutate(end_position = start_position+24) # differentiate for adapters and barcodes!
  
  read_lengths<-read.csv2(file = paste(batch, "_read_lengths.tsv", sep=""), sep = "\t", header = F)
  colnames(read_lengths)<-c("read", "length")
  read_lengths<-read_lengths%>%
    separate(read, into = c("read", NA), sep=" ")
  
  sam_read<-sam_read %>%
    left_join(read_lengths, by = "read")
  
  return(sam_read)
}


combine_sam_files<-function(names){
  combined_sam<-data.frame()
  for (i in 1:length(names)) {
    sam_batch<-read_sam(names[i])
    combined_sam<-combined_sam %>%
      bind_rows(
        sam_batch
      )
  }
  return(combined_sam)
}

sam_files<-combine_sam_files(sample_names)

sam_bar_sum<-sam_files %>%
  filter(barcode != "seq_adp") %>%
  group_by(read, batch) %>%
  summarise(n_barcode=length(barcode),
            unique_barcode=length(unique(barcode)))

sam_adp_sum<-sam_files %>%
  filter(barcode == "seq_adp") %>%
  group_by(read, batch) %>%
  summarise(n_adp=length(barcode))

sam_sum <- sam_bar_sum %>%
  left_join(sam_adp_sum, by = c("read", "batch"))



ggplot(sam_sum, aes(x = n_barcode)) +
  geom_histogram(binwidth=1,boundary=-0.5) +
  scale_x_continuous(breaks=1:max(sam_sum$n_barcode)) +
  theme_minimal()+
  facet_grid(rows = vars(batch))


sam_sum_counts <- sam_sum %>%
  group_by(batch) %>%
  count(n_barcode, unique_barcode) %>%
  mutate_all(~replace_na(.,0))

ggplot(sam_sum_counts, aes(x=n_barcode, y=unique_barcode, fill = n)) + 
  geom_tile() +
  theme_minimal()+
  facet_grid(rows=vars(batch))+
  scale_x_continuous(breaks=1:max(sam_sum_counts$n_barcode)) +
  scale_y_continuous(breaks=1:9) +
  geom_text(aes(label = n), color = "white", size = 3)


##Visualising the structure of barcodes for 
sam_many_barcodes <- sam_files %>%
  filter(read %in% filter(sam_sum, n_barcode > 15)$read)

colors <- colorRampPalette(brewer.pal(8, "Set2"))(17)

ggplot(sam_many_barcodes, aes(xmin=start_position, xmax=end_position, y=read, fill = barcode, label = barcode)) +
  facet_wrap(~read, scales = "free", ncol = 1) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "left", size = 0.01) +
  scale_fill_manual(values=colors) +
  theme_genes()+
  theme(legend.position = "top",
        axis.text.y = element_blank())
