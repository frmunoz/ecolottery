# Function to show the distribution of trait in the community vs. distribution
# of trait in the regional pool, or the abundances in the local community vs.
# in the regional pool.
plot_comm <- function(x, type = "trait", seltrait = 1, main = NULL) 
# x should be the output of coalesc or forward
# seltrait is the index of the trait to be plotted (in case of multiple traits)
{
  
  # Check parameters
  if (!is.numeric(seltrait) | seltrait <= 0){
    stop(paste0("seltrait must be the number of the column of desired trait",
                " in the community data.frame."))
  }
  
  if (!is.null(main)){
    if (!is.character(main)){
      stop("main must be a character string describing the title of the plot.")
    }
  }
  
  if(type == "abund")
    stop("Type abund is no longer supported, see options in ?abund")
  
  if(type != "trait")
    # Compute community and regional pool abundances by species
    ab <- abund(x) 
  
  switch(type,
    "trait" =
     {
       # Color vectors
       metaCol <- rgb(1,0,0,0.2)
       localCol <- rgb(0,0,1,0.2)
       # data.frame of trait pool and local values
       data <- data.frame(
         level = c(rep("pool", nrow(x$pool)), rep("comm", nrow(x$com))),
         trait = c(x$pool[, seltrait + 2], x$com[, seltrait + 2]))
       # Plot
       if (requireNamespace("ggplot2", quietly = TRUE)) {
         
         ggplot(data, aes_string(x = "trait"), main = main) +
         geom_density(aes_string(group = "level", fill = "level"),
                      alpha = 0.5) +
         scale_fill_manual(values = c("pool" = metaCol, "comm" = localCol)) +
         ggtitle(main)
         
      } else {
         tmin <- min(c(x$pool[, seltrait + 2],x$com[, seltrait + 2]))
         tmax <- max(c(x$pool[, seltrait + 2],x$com[, seltrait + 2]))
         h1 <- hist(x$pool[, seltrait + 2], plot = FALSE, breaks = 10)
         h1$density <- h1$counts/sum(h1$counts)*100
         h2 <- hist(x$com[, seltrait + 2], plot = FALSE, breaks = 10)
         h2$density <- h2$counts/sum(h2$counts)*100
         plot(h1,
              xlim = c(0.9*tmin,1.1*tmax), main = main, xlab = "trait",
              ylim = c(0,max(c(h1$density, h2$density))),
              col = metaCol, freq = FALSE)
         plot(h2, col = localCol, add = TRUE, freq = FALSE)
      }
     },
    "locreg" =
     {
       # Plot graphics
       plot(ab$pool[rownames(ab$com),"relab"],
            ab$com$relab,
            main = main,
            xlab = "Regional abundance",
            ylab = "Local abundance",
            log = "xy")
       abline(0,1)
     },
    "sad" =
    {
      #Fits log-series distribution to abundance data
      ab.com.ls <- sads::fitsad(ab$com$ab, "ls")
      
      #Show species abundance distributions
      sads::ppsad(ab.com.ls)
    },
    "rad" =
    {
      #Fits geometric series to abundance data
      ab.com.gs <- sads::fitrad(ab$com$ab, "gs")
      
      # Show rank abundance distributions
      sads::pprad(ab.com.gs)
    },
    warning("Need to choose a type of plot to display (trait, locreg, sad or rad)"))
}
