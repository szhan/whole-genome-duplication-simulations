require(dplyr)
require(diversitree)
require(optparse)
require(parallel)
require(stringr)


# Example command:
# Rscript --no-save simulate_histories.R \
#   --spec0=1 --spec1=0.4 --ext0=0.5 --ext1=0.50 --q01=0.5 --q10=0.65 \
#   --ntaxa=10000 --nreps=100 --out_file=s1_e3_q4.csv

option_list = list(
  make_option(c("--spec0"),
              help = "Diploid speciation rate",
              type = "double",
              metavar = "double",
              default = NULL),
  make_option(c("--spec1"),
              help = "Polyploid speciation rate",
              type = "double",
              metavar = "double",
              default = NULL),
  make_option(c("--ext0"),
              help = "Diploid extinction rate",
              type = "double",
              metavar = "double",
              default = NULL),
  make_option(c("--ext1"),
              help = "Polyploid extinction rate",
              type = "double",
              metavar = "double",
              default = NULL),
  make_option(c("--q01"),
              help = "Diploid-to-polyploid transition rate",
              type = "double",
              metavar = "double",
              default = NULL),
  make_option(c("--q10"),
              help = "Polyploid-to-diploid transition rate",
              type = "double",
              metavar = "double",
              default = NULL),
  make_option(c("--ntaxa"),
              help = "Number of extant tip taxa",
              type = "integer",
              metavar = "integer",
              default = NULL),
  make_option(c("--nreps"),
              help = "Number of replicates",
              type = "integer",
              metavar = "integer",
              default = NULL),
  make_option(c("--out_file"),
              help = "Output file name",
              type = "character",
              metavar = "character",
              default = NULL)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print("Input parameters")
print(paste0("spec0     : ", opt$spec0))
print(paste0("spec1     : ", opt$spec1))
print(paste0("ext0      : ", opt$ext0))
print(paste0("ext1      : ", opt$ext1))
print(paste0("q01       : ", opt$q01))
print(paste0("q10       : ", opt$q10))
print(paste0("ntaxa     : ", opt$ntaxa))
print(paste0("nreps     : ", opt$nreps))
print(paste0("out_file  : ", opt$out_file))


# count number of 0->1 transitions in a binary sequence
count_transitions_dp <- function(x){
  return(str_count(paste0(x, collapse = ""), "01"))
}


get_number_transitions <- function(phy, his){
  nodes <- names(his$node.state)[-1]
  tips  <- names(his$tip.state)
  
  trans_per_lineage <- rep(0, length(tips))
  trans_overall     <- 0
  
  for(i in 1:length(nodes)){
    his_df      <- data.frame(his$history[nodes[i]])
    
    if(nrow(his_df) == 1)
      next
    
    nbr_trans   <- count_transitions_dp(his_df[, 2])
    list_desc   <- get.descendants(nodes[i],
                                   phy,
                                   tips.only = TRUE)
    
    trans_per_lineage[list_desc] <- trans_per_lineage[list_desc] + nbr_trans
    trans_overall <- trans_overall + nbr_trans * length(list_desc)
  }
  
  for(i in 1:length(tips)){
    his_df      <- data.frame(his$history[tips[i]])
    
    if(nrow(his_df) == 1)
      next
    
    nbr_trans   <- count_transitions_dp(his_df[, 2])
    
    trans_per_lineage[i] <- trans_per_lineage[i] + nbr_trans
    trans_overall <- trans_overall + nbr_trans
  }
  
  #mean_trans   <- trans_overall / length(his$tip.state)
  median_trans <- median(trans_per_lineage)
  
  return(median_trans)
}


is_wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}


# wrapper for mclapply
simulate_bisse_tree <- function(pars){
  stopifnot(length(pars) == 7)
  stopifnot(all(is.double(pars[1:6])))
  stopifnot(is_wholenumber(pars[7]))
  
  phy   <- trees(pars = pars[1:6], # rate parameters
                 max.taxa = pars[7], # number of extant taxa
                 n    = 1, # number of simulations
                 x0   = 0, # initial root state
                 type = "bisse",
                 include.extinct = FALSE)
  
  his   <- history.from.sim.discrete(phy[[1]], 0:1)
  median_nbr_transitions <- get_number_transitions(phy[[1]], his)
  nbr_tips_polyploid <- sum(his$tip.state)
  
  # keep model parameters and simulated tree plus evolutionary history
  random_word <- paste0(sample(letters, 20, TRUE), collapse = "")
  out_save_file <- paste0(random_word, ".RData")
  save(phy, his, pars, median_nbr_transitions,
       opt, .Random.seed, # global variables
       file = out_save_file)
  
  # TODO: there should be a better way to return a vector
  # that can be fed into mclapply
  return(paste0(median_nbr_transitions, ",", nbr_tips_polyploid))
}


results <- mclapply(rep(list(c(opt$spec0, # diploid speciation rate
                               opt$spec1, # polyploid speciation rate
                               opt$ext0,  # diploid extinction rate
                               opt$ext1,  # polyploid extinction rate
                               opt$q01,   # diploid-to-polyploid transition
                               opt$q10,   # polyploid-to-diploid transition
                               opt$ntaxa  # number of extant taxa
                               )),
                        opt$nreps),       # number of simulations
                    simulate_bisse_tree,
                    mc.cores = detectCores())

output <- tibble(replicate = 1:opt$nreps,
                 ntaxa     = rep(opt$ntaxa, opt$nreps),
                 spec0     = rep(opt$spec0, opt$nreps),
                 spec1     = rep(opt$spec1, opt$nreps),
                 ext0      = rep(opt$ext0,  opt$nreps),
                 ext1      = rep(opt$ext1,  opt$nreps),
                 q01       = rep(opt$q01,   opt$nreps),
                 q10       = rep(opt$q10,   opt$nreps),
                 median_nbr_transitions = unlist(results))


write.csv(x         = output,
          file      = opt$out_file,
          row.names = FALSE,
          quote     = FALSE)
