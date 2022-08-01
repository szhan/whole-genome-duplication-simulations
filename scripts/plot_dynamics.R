require(cowplot)
require(dplyr)
require(ggplot2)
require(readr)
require(viridis)


in_dir <- "data/Sol/"
title_text <- "Solanaceae-like BiSSE parameters"
out_pdf_file <- "bisse_sims_solanaceae_like.pdf"

s_pref <- 'a'
s_list <- 1:1
e_pref <- 'b'
e_list <- 1:1
q_pref <- 'c'
q_list <- 1:5

sims <- tibble()
for(i in s_list){
  for(j in e_list){
    for(k in q_list){
      csv_file <- paste0(in_dir,
                         s_pref, i, '_',
                         e_pref, j, '_',
                         q_pref, k, '.csv')
      print(csv_file)
      tmp <- as_tibble(read.csv(csv_file))
      sims <- bind_rows(sims, tmp)
    }
  }
}

min_wgds <- floor(min(sims$nbr_transitions))
max_wgds <- ceiling(max(sims$nbr_transitions))

set.seed(123456)

d <- sims %>%
  ggplot(aes(y = nbr_transitions,
             x = q10)) +
  geom_hex(color = "white",
           bins = 10) +
  scale_fill_viridis_c() +
  ylim(min_wgds, max_wgds) +
  ggtitle(title_text) +
  ylab("Number of diploid-to-polyploid transitions\n(mean across all lineages)") +
  xlab("Rate of polyploid-to-diploid transition") +
  labs(fill = "Number of\nsimulations") +
  theme_minimal() +
  theme(title = element_text(size = 14),
        panel.background = element_rect(fill = "grey92",
                                        colour = NA),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        axis.text.x  = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
  )

pdf(out_pdf_file, useDingbats = FALSE, width = 8.5, height = 7)
d
dev.off()
