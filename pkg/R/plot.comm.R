# Plotting of regional and local trait distributions

plot.comm <- function(x, seltrait=1,main=NULL) 
# x should be the output of coalesc or forward
# seltrait is the index of the trait to be plotted (in case of multiple traits)
{
  require(ggplot2)
  # Color vectors
  metaCol <- rgb(1,0,0,0.2); localCol <- rgb(0,0,1,0.2)
  
  # Dataframe of trait pool and local values
  data <- data.frame(level=c(rep("pool", nrow(x$pool)),rep("comm", nrow(x$com))),
                     trait=c(x$pool[, seltrait+2],x$com[, seltrait+2]))
  
  # Plot
  ggplot(data, aes(trait),main=main) +
    geom_density(aes(group = level, fill = level), alpha = 0.5) +
    scale_fill_manual(values = c(metaCol, localCol)) +
    ggtitle(main) + 
    theme_classic()
}
 
