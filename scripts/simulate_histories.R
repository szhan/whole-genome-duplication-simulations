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


# Count number of 0->1 or 1->0 transitions in a string of 0s and 1s.
# TODO: Processing strings is highly inefficient.
count_transitions <- function(x, reverse = FALSE){
  if(reverse){
    pattern <- "10"
  }else{
    pattern <- "01" # TRUE
  }
  nbr_transitions <- str_count(paste0(x, collapse = ""), pattern)
  return(nbr_transitions)
}


# Count number of 0->1 transitions from the root to each tip,
# and return the number of transitions across the tips.
get_number_transitions_per_lineage <- function(phy, # 'phylo' object
                                               his){
  nodes <- names(his$node.state)[-1] # Exclude root node (listed as first)
  tips <- names(his$tip.state)
  
  # Tally up transitions from root to each tip
  trans_per_path_01 <- rep(0, length(tips))
  trans_per_path_10 <- rep(0, length(tips))
  
  # Number of transitions in the history of an internal node is shared by
  # its descendants at the tips.
  for(i in 1:length(nodes)){
    tmp_df <- data.frame(his$history[nodes[i]])
    
    if(nrow(tmp_df) == 1)
      next
    
    nbr_trans_01 <- count_transitions(tmp_df[, 2], reverse = FALSE)
    nbr_trans_10 <- count_transitions(tmp_df[, 2], reverse = TRUE)
    list_desc <- get.descendants(nodes[i], # Internal node ID
                                 phy,
                                 tips.only = TRUE)
    
    trans_per_path_01[list_desc] <- trans_per_path_01[list_desc] + nbr_trans_01
    trans_per_path_10[list_desc] <- trans_per_path_10[list_desc] + nbr_trans_10
  }
  
  for(i in 1:length(tips)){
    tmp_df <- data.frame(his$history[tips[i]])
    
    if(nrow(tmp_df) == 1)
      next
    
    nbr_trans_01 <- count_transitions(tmp_df[, 2], reverse = FALSE)
    nbr_trans_10 <- count_transitions(tmp_df[, 2], reverse = TRUE)
    
    trans_per_path_01[list_desc] <- trans_per_path_01[list_desc] + nbr_trans_01
    trans_per_path_10[list_desc] <- trans_per_path_10[list_desc] + nbr_trans_10
  }
  
  names(trans_per_path_01) <- tips
  names(trans_per_path_10) <- tips # Not returned
  
  # It is not possible to have more reverse than forward transitions
  # if root state is set to 0.
  # TODO: Investigate why sometimes these assertions are false.
  #stopifnot(all(trans_per_path_01 >= trans_per_path_10))
  #stopifnot(max(trans_per_path_01 - trans_per_path_10) == 0 ||
  #          max(trans_per_path_01 - trans_per_path_10) == 1)
  
  result <- list()
  result[[1]] <- trans_per_path_01
  result[[2]] <- trans_per_path_10
  names(result) <- c("trans_per_path_01",
                     "trans_per_path_10")
  
  return(result)
}


# Helper function
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}


# Wrapper for mclapply()
simulate_bisse_tree <- function(pars){
  stopifnot(length(pars) == 7)
  stopifnot(all(is.double(pars[1:6])))
  stopifnot(is_wholenumber(pars[7]))
  
  phy <- trees(pars = pars[1:6], # Rate parameters
               max.taxa = pars[7], # Nbr of extant taxa
               n = 1, # Nbr of simulations
               x0 = 0, # Initial root state
               type = "bisse",
               include.extinct = FALSE)[[1]]
  
  his <- history.from.sim.discrete(phy, 0:1)
  nbr_trans_per_lineage <- get_number_transitions_per_lineage(phy, his)
  mean_nbr_trans_01 <- mean(nbr_trans_per_lineage$trans_per_path_01)
  nbr_tips_0 <- sum(his$tip.state == 0) # Not kept or returned
  nbr_tips_1 <- sum(his$tip.state == 1) # Not kept or returned
  
  # Keep model parameters and simulated tree plus character history.
  random_word <- paste0(sample(letters, 20, TRUE), collapse = "")
  out_save_file <- paste0(random_word, ".RData")
  save(phy,
       his,
       pars,
       nbr_trans_per_lineage,
       #opt,
       .Random.seed,
       file = out_save_file)
  
  return(mean_nbr_trans_01)
}


results <- mclapply(rep(list(c(opt$spec0, # Diploid speciation rate
                               opt$spec1, # Polyploid speciation rate
                               opt$ext0,  # Diploid extinction rate
                               opt$ext1,  # Polyploid extinction rate
                               opt$q01,   # Diploid-to-polyploid transition rate
                               opt$q10,   # Polyploid-to-diploid transition rate
                               opt$ntaxa  # Nbr of extant taxa
                               )),
                        opt$nreps),       # Nbr of simulations
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
                 nbr_transitions = unlist(results))

write.csv(x         = output,
          file      = opt$out_file,
          row.names = FALSE,
          quote     = FALSE)
