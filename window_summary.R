###########################################
## Sliding widow plots function           #
###########################################

# writen by Anders Goncalves da Silva (https://github.com/andersgs)

#' window_summary: Summarize SNP statistics along a chromosome given a certain
#' window size, and offset
#' 
#' @parm data A data.table or data.frame containing the position of SNPs
#'            and the statistic(s) of interest
#' @param size An integer indicating the size of window in base-pairs
#' @param offset An integer indicating the offset for overlap between
#'               windows. A value of 0 indicates no overlap (default: 0)
#' @param pos_name Name of the column that has positions
#' @param stat_name Name of the column that has the position of interest
#' @param FUN Name of the function to apply (e.g., mean)
#'               
#'               

window_summary <- function(data, 
                           size, 
                           offset = 0, 
                           pos_name = "Position",
                           stat_name = "HetObsv", 
                           FUN = mean) {
  # test of dependencies
  miss_libs = c()
  if(!require(data.table)) {
    miss_libs = c(miss_libs, "data.table")
  }
  
  if(length(miss_libs > 0)) {
    print("Could not find some libraries that are needed to run properly.")
    print("To install the necessary libraries, type:")
    libs_install <- paste("\"", miss_libs, "\"", sep = "", collapse = ",")
    cat(paste("install.packages(c(", libs_install, "))\n", sep = ""))
    return("Please install necessary packages to continue")
  }
  
  # make it into a data.table
  if(!data.table::is.data.table(data)) {
    data <- data.table::data.table(data)
  }
  
  # get distribution of inter-SNP distances
  #dist_dist <- data[,get(pos_name)[2:nrow(data)] - get(pos_name)[1:(nrow(data)-1)]]
  
  #create windows
  pos_range <- range(data[[pos_name]])
  print(pos_range)
  if(offset == 0) {
    offset = size + 1
  }
  start_pos <- seq(pos_range[1], pos_range[2] - size, offset)
  end_pos <- seq(pos_range[1] + size, pos_range[2], offset)
  # fixing the ends as seq will not go all the way if the offset is not a 
  # multiple. It means the last bin will be shorter
  if(end_pos[length(end_pos)] < pos_range[2]) {
    start_pos <- c(start_pos, end_pos[length(end_pos)] + 1)
    end_pos <- c(end_pos, pos_range[2])
  }
  wins <- data.table(start = start_pos, end = end_pos)
  # group snps by windows
  wins_dt <- apply(wins, 1, function(win) {
    ix <- which(data[,get(pos_name)] >= win[1] &
            data[,get(pos_name)] <= win[2])
    if(nrow(data[ix,]) == 0)
      NULL
    else{
      data.table(start = win[1], 
                 end = win[2],
                 N = length(ix),
                 stat = data[ix,FUN(get(stat_name))])
    }
  })
  ix <- sapply(wins_dt, is.null)
  wins_dt <- data.table::rbindlist(wins_dt[which(!(ix))])
  return(wins_dt)
}
