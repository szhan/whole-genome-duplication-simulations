require(cowplot)
require(dplyr)
require(ggplot2)


p_x <- list()

pth <- 0
for(i in 1:5){
  for(j in 1:5){
    sims <- tibble()
    for(k in 1:5){
      csv_file <- paste0("s", i, "_e", j, "_q", k, ".csv")
      tmp <- as_tibble(read.csv(csv_file))
      sims <- bind_rows(sims, tmp)
    }
    sims <- sims %>%
      mutate(median_nbr_transitions = trunc(median_nbr_transitions))
    
    max_wgds <- max(sims$median_nbr_transitions)
    rates_q10 <- sort(unique(sims$q10))
    wgds_by_q10 <- expand.grid(y = max_wgds:1,
                               x = rates_q10)
    wgds_by_q10$z <- rep(0, nrow(wgds_by_q10))
    for(h in 1:nrow(sims)){
      wgds_by_q10[wgds_by_q10$y == sims[h,]$median_nbr_transitions &
                    wgds_by_q10$x == sims[h,]$q10, ]$z <- 
        wgds_by_q10[wgds_by_q10$y == sims[h,]$median_nbr_transitions &
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
    print(pth)
    p_x[[pth]] <- p_tmp
    p_x <- c(p_x, p_tmp)
  }
}


pdf("rate_dynamics.pdf", useDingbats = FALSE)
plot_grid(p_x[[1]] , p_x[[2]] , p_x[[3]] , p_x[[4]] , p_x[[5]],
          p_x[[6]] , p_x[[7]] , p_x[[8]] , p_x[[9]] , p_x[[10]],
          p_x[[11]], p_x[[12]], p_x[[13]], p_x[[14]], p_x[[15]],
          p_x[[16]], p_x[[17]], p_x[[18]], p_x[[19]], p_x[[20]],
          p_x[[21]], p_x[[22]], p_x[[23]], p_x[[24]], p_x[[25]],
          nrow = 5,
          ncol = 5)
dev.off()
