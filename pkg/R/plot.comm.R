# Plotting of regional and local trait or abundance distributions

plot.comm <- function(x, type="trait", seltrait=1,main=NULL) 
# x should be the output of coalesc or forward
# seltrait is the index of the trait to be plotted (in case of multiple traits)
{
  switch(type,
    "trait"=
     {
       require(ggplot2)
       # Color vectors
       metaCol <- rgb(1,0,0,0.2); localCol <- rgb(0,0,1,0.2)
       # Dataframe of trait pool and local values
       data <- data.frame(level=c(rep("pool", nrow(x$pool)),rep("comm", nrow(x$com))),
                     trait=c(x$pool[, seltrait+2],x$com[, seltrait+2]))
       # Plot
       ggplot(data, aes(trait),main=main) + geom_density(aes(group = level, fill = level), alpha = 0.5) +
        scale_fill_manual(values = c(metaCol, localCol)) +
        ggtitle(main) + 
        theme_classic()
     },
    "abund"=
     {
      ab <- abund(x) 
      plot(ab$reg[names(ab$loc)]/sum(ab$reg),ab$loc/sum(ab$loc),main=main,xlab="Regional abundance", ylab="Local abundance",log="xy"); abline(0,1)
     },
    warning("Need to choose a type of plot to display (trait or abund)"))
}
 
