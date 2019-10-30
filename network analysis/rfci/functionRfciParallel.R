rfci_parallel<- function (suffStat, indepTest, alpha, labels, p, skel.method = c("parallel"), 
          mem.efficient = FALSE, fixedGaps = NULL, fixedEdges = NULL, 
          NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10), conservative = FALSE, 
          maj.rule = FALSE, verbose = FALSE, num.cores = detectCores()) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  if (conservative && maj.rule) 
    stop("Can only choose one of conservative or majority rule RFCI")
  if (verbose) 
    cat("Compute Skeleton\n================\n")
  num_workers <- num.cores
  if (num_workers < 2) {
    stop("The number of cores is insufficient to run parallel-rfci")
  }
  workers <- NULL
  skel <- skeleton_parallel(suffStat, indepTest, alpha, labels = labels, 
                            method = skel.method, workers = workers, num_workers = num_workers, 
                            fixedGaps = fixedGaps, fixedEdges = fixedEdges, mem.efficient = mem.efficient, 
                            NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  sk.A <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  u.t <- find.unsh.triple(sk.A, check = FALSE)
  r.v. <- rfci.vStruc(suffStat, indepTest, alpha, sepset, sk.A, 
                      unshTripl = u.t$unshTripl, unshVect = u.t$unshVect, conservative = (conservative || 
                                                                                            maj.rule), version.unf = c(1, 1), maj.rule = maj.rule, 
                      verbose = verbose)
  A <- r.v.$amat
  sepset <- r.v.$sepset
  if (Sys.info()[["sysname"]] == "Windows") {
    stopCluster(workers)
  }
  if (verbose) 
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), 
        "\n")
  res <- udag2apag(A, suffStat, indepTest, alpha, sepset, rules = rules, 
                   unfVect = r.v.$unfTripl, verbose = verbose)
  Amat <- res$graph
  colnames(Amat) <- rownames(Amat) <- labels
  new("fciAlgo", amat = Amat, call = cl, n = integer(0), max.ord = as.integer(skel@max.ord), 
      max.ordPDSEP = 0L, n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0, 
      sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", 
                                                               p))
}