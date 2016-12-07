# Function to show the distribution of trait in the community vs. distribution
# of trait in the regional pool, or the abundances in the local community vs.
# in the regional pool.
plot_comm <- function(x, type = "trait", seltrait = 1, main = NULL) 
# x should be the output of coalesc or forward
# seltrait is the index of the trait to be plotted (in case of multiple traits)
{
  
  # Check parameters
  if (!is.numeric(seltrait) | seltrait <= 0){
    stop("seltrait must be the number of the column of desired trait in the community data.frame.")
  }
  
  if (!is.null(main)){
    if (!is.character(main)){
      stop("main must be a character string describing the title of the plot.")
    }
  }
  
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
         ggplot2::ggplot(data, ggplot2::aes_string(x = "trait"), main = main) +
         ggplot2::geom_density(ggplot2::aes_string(group = "level", fill = "level"), alpha = 0.5) +
         ggplot2::scale_fill_manual(values = c("pool" = metaCol, "comm" = localCol)) +
         ggplot2::ggtitle(main)
       } else {
         tmin <- min(c(x$pool[, seltrait + 2],x$com[, seltrait + 2]))
         tmax <- max(c(x$pool[, seltrait + 2],x$com[, seltrait + 2]))
         h1 <- hist(x$pool[, seltrait + 2], plot=F, breaks = 10)
         h1$density=h1$counts/sum(h1$counts)*100
         h2 <- hist(x$com[, seltrait + 2], plot=F, breaks = 10)
         h2$density=h2$counts/sum(h2$counts)*100
         plot(h1, xlim=c(0.9*tmin,1.1*tmax), ylim=c(0,max(c(h1$density, h2$density))), main=main, xlab="trait", col=metaCol, freq=F)
         plot(h2, col=localCol, add=T, freq=F)
      }
     },
    "abund" =
     {
       # Compute community and regional pool abundances by species
       ab <- abund(x) 
       
       # Plot graphics
       plot(ab$pool[rownames(ab$com),"relab"],
            ab$com$relab,
            main = main,
            xlab = "Regional abundance",
            ylab = "Local abundance",
            log = "xy")
       abline(0,1)
     },
    warning("Need to choose a type of plot to display (trait or abund)"))
}
