get_number_of_gens <- function(given_size, pool, nbrep = 5, prob = 1, d = 1,
                               gens = NULL, limit.sim = FALSE,
                               coeff.lim.sim = 1, sigm = 0.1, filt = NULL,
                               prob.death = NULL, method.dist = "euclidean",
                               plot_gens = FALSE) {
  
    if (is.character(pool)) {
      pool <- data.frame(id = 1:length(pool), sp = pool, trait = rep(NA, length(pool)), stringsAsFactors = FALSE)
      if (limit.sim | !is.null(filt)) {
        cat("No trait information provided in the regional pool\n") 
        pool[, 3] <- runif(nrow(pool))
        cat(paste0("Random trait values attributed to individuals of the ", "regional pool\n"))
        colnames(pool) <- c("id", "sp", "trait")
      }
    }
    
    # At least 100 total replacements in expectation
    if(is.null(gens)) gens <- 100*given_size/d
    gens = min(100*given_size/d, gens, na.rm=T) 
    
    # Simulate starting community from given pool
    start_com <- pool[sample(1:nrow(pool), given_size),]
    
    if (is.null(filt)) {
      filt <- function(x) return(1)
    }
    
    gens_conv_cpt <- c()
    nb_sp_gen <- c()
    
    # Loop to simulate the communities and determine changepoint
    for (i in 1:nbrep) {
      final <- forward(initial = start_com, prob = prob, d = d, gens = gens,
                       keep = FALSE, pool = pool, limit.sim = limit.sim, 
                       coeff.lim.sim = coeff.lim.sim, sigm = sigm, filt = filt,
                       prob.death = prob.death, method.dist = method.dist,
                       plot_gens = plot_gens)
      
      nb_sp_gen <- rbind(nb_sp_gen, data.frame(gens = 1:gens, rich = final$sp_t,
                                               stringsAsFactors = FALSE))
      
      # Customized non-linear regression
      final_sp <- final$sp_t[length(final$sp_t)]
      init_sp <- length(unique(start_com$sp))
      changePoint <- function(t, spf, Tval, init_sp) {
        init_sp*exp(log(spf/init_sp)*t/Tval)*as.numeric(t<=Tval)+spf*as.numeric(t>Tval)
      }
      sqerror <- function (par, x, t) {
        sum((x - changePoint(t, par[1], par[2], init_sp))^2)
      }
      sp.fit <- optim(par = c(final_sp, median(data$t)), fn = sqerror, x = data$sp, t = data$t, lower=c(1,1), 
                      upper=c(given_size,gens), method = "L-BFGS-B")
      gens_conv_cpt <- c(gens_conv_cpt, sp.fit$par[2])
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
    
    return(list(n = gens_conv_cpt, plot = diagnostic_plot))
  
}
