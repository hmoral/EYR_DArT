###########################################
## Sliding widow plots                    #
###########################################

# writen by Anders Goncalves da Silva (https://github.com/andersgs)

# USAGE
# This script has to be used "manually"
# It does NOT take any input from the command line

##################
# Load libraries #
##################
library(ggplot2)
##################
# Read function  #
##################
source("window_summary.R")
##################
# set directory  #
##################
setwd("./")

##################################
# Read input parameters and data #
##################################
ata <- read.table("Example_Anders_Chr1A.txt", 
                   stringsAsFactors = F,
                   header = T)

tmp_mean <- window_summary(data = data, size = 5000000, offset = 500000)
tmp_sd <- window_summary(data = data, size = 5000000, offset = 500000, FUN = sd)
setkeyv(tmp_mean, c("start", "end"))
tmp <- tmp_mean[tmp_sd[,list(start,end,stat)]]
setnames(tmp, old = c("stat", "i.stat"), new = c("mean", "sd"))

p1 <- ggplot(data, aes(x = Position, y = HetObsv)) + 
  geom_point(color = "gray50") +
  geom_line(data = tmp, aes(x = start, y = mean), size = 1) +
  geom_point(data = tmp, aes(x = start, y = mean, size = N), color = 'red') +
  geom_ribbon(data = tmp, aes(x = start, 
                              y = mean, 
                              ymin = mean - sd, 
                              ymax = mean + sd),
              alpha = 0.3) +
  ylab("Observed He")
#ggsave(filename = "window_plot.png", plot = p1)
