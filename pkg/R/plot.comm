# Plotting of regional and local trait distributions

plot.comm <- function(x,seltrait=1) 
# x should be the output of coalesc or forward
# seltrait is the index of the trait to be plotted (in case of multiple traits)
{
  require(ggplot2)
  metaCol <- rgb(1,0,0,0.2); localCol <- rgb(0,0,1,0.2)
  data <- data.frame(level=c(rep("pool",nrow(x$pool)),rep("comm",nrow(x$com))),trait=c(x$pool[,seltrait+2],x$com[,seltrait+2]))
  qplot(trait, data = data, geom = "density" , group = level, fill = level, alpha=.3)
}
 
