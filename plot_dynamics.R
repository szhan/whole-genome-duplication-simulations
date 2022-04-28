require(cowplot)
require(dplyr)
require(ggplot2)
require(viridis)


in_dir <- "data/sim8/"

s_pref <- 'x'
s_list <- 1:5
e_pref <- 'y'
e_list <- 1:5
q_pref <- 'z'
q_list <- 1:5

out_pdf_file <- "rate_dynamics.pdf"

p_x <- list()

pth <- 0
for(i in s_list){
  for(j in e_list){
    sims <- tibble()
    for(k in q_list){
      csv_file <- paste0(in_dir, s_pref, i, '_', e_pref, j, '_', q_pref, k, '.csv')
      tmp <- as_tibble(read.csv(csv_file))
      tmp <- tmp[tmp$replicate == 1, ] # Take only replicate 1
      print(tmp)
      sims <- bind_rows(sims, tmp)
    }
    sims <- sims %>%
      mutate(mean_nbr_transitions = trunc(mean_nbr_transitions))
    
    max_wgds <- ceiling(max(sims$mean_nbr_transitions))
    rates_q10 <- sort(unique(sims$q10))
    wgds_by_q10 <- expand.grid(y = max_wgds:1,
                               x = rates_q10)
    wgds_by_q10$z <- rep(0, nrow(wgds_by_q10))
    for(h in 1:nrow(sims)){
      wgds_by_q10[wgds_by_q10$y == sims[h,]$mean_nbr_transitions &
                    wgds_by_q10$x == sims[h,]$q10, ]$z <- 
        wgds_by_q10[wgds_by_q10$y == sims[h,]$mean_nbr_transitions &
                      wgds_by_q10$x == sims[h,]$q10, ]$z + 1
    }
    
    p_tmp <- wgds_by_q10 %>%
      ggplot(aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_y_continuous(breaks = max_wgds:1,
                         labels = max_wgds:1) +
      scale_x_continuous(breaks = rates_q10,
                         labels = rates_q10) +
      scale_fill_viridis(discrete = FALSE) +
      theme_minimal() +
      theme(axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.text.y  = element_text(size = 8),
            axis.text.x  = element_text(size = 6),
            legend.position = "none")
    
    pth <- pth + 1
    p_x[[pth]] <- p_tmp
  }
}


pdf(out_pdf_file, useDingbats = FALSE)
grid_size <- floor(sqrt(length(p_x)))^2
plot_grid(
  p_x[[1]],
  nrow = grid_size,
  ncol = grid_size
)
dev.off()
