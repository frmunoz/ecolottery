get_number_of_gens <- function(given_size, pool, nbrep = 100, prob = 0, d = 1,
                               gens = NULL, limit.sim = F, coeff.lim.sim = 1,
                               sigm = 0.1, filt = NULL, prob.death = NULL,
                               method.dist = "euclidean", plot_gens = FALSE) {
  
  if (requireNamespace("changepoint", quietly = TRUE)) {
    
    # Simulate starting community from given pool
    start_com <- pool[sample(1:nrow(pool), given_size),]
    
    # Gives the number of generation to test
    if (is.null(gens)) {
      gens <- 100*given_size
    }
    
    if (is.null(filt)) {
      filt <- function(x) return(1)
    }
    
    gens_conv_cpt <- c()
    nb_sp_gen <- c()
    
    # Loop to simulate the communities and determine changepoint
    for (i in 1:nbrep) {
      final <- forward(initial = start_com, prob = prob, d = d, gens = gens,
                       keep = F, pool = pool, limit.sim = limit.sim, 
                       coeff.lim.sim = coeff.lim.sim, sigm = sigm, filt = filt,
                       prob.death = prob.death, method.dist = method.dist,
                       plot_gens = plot_gens)
      
      
      nb_sp <- final$sp_t
      nb_sp_gen <- rbind(nb_sp_gen, data.frame(gens = 1:gens, rich = nb_sp,
                                               stringsAsFactors = F))
      
      cpt_bic <- changepoint::cpt.mean(nb_sp, penalty = "BIC")
      
      gens_conv_cpt <- c(gens_conv_cpt, changepoint::cpts(cpt_bic))
    }
    
    # Average number of generations
    nb_sp_gen <- data.frame(apply(nb_sp_gen, 2,
                                  function(x) tapply(x, nb_sp_gen$gens, mean)))
    
    # Make a plot to check if it is really convergent
    diagnostic_plot <- ggplot(nb_sp_gen, aes_string("gens", "rich")) +
      geom_point(alpha = 0.5) +
      geom_vline(xintercept = max(gens_conv_cpt), colour = "red",
                 linetype = "dashed") +
      labs(x = "Generation number", y = "Specific richness") +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")) +
      theme_bw()
    
    return(list(n = max(gens_conv_cpt), plot = diagnostic_plot))
  
  } else {
    stop("Install package 'changepoint' to compute number of generations")
  }
  
}
