#modified version of the skeleton-parallel function
skeleton_parallel <- function (suffStat, indepTest, alpha, labels, p, method = c("parallel"), 
          mem.efficient = FALSE, workers, num_workers, m.max = Inf, 
          fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, verbose = FALSE) 
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
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
  }
  else if (!identical(dim(fixedGaps), c(p, p))) 
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps))) 
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p))) 
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (fixedEdges != t(fixedEdges)) 
    stop("fixedEdges must be symmetric")
  pval <- NULL
  sepset <- lapply(seq_p, function(.) vector("list", p))
  pMax <- matrix(-Inf, nrow = p, ncol = p)
  diag(pMax) <- 1
  done <- FALSE
  ord <- 0
  n.edgetests <- numeric(1)
  edge_test_xy <- function(x, y) {
    G_xy <- TRUE
    num_tests_xy <- 0
    pMax_xy <- pMax[x, y]
    sepset_xy <- NULL
    done_xy <- TRUE
    if (G_xy && !fixedEdges[y, x]) {
      nbrsBool <- G.l[[x]]
      nbrsBool[y] <- FALSE
      nbrs <- seq_p[nbrsBool]
      length_nbrs <- length(nbrs)
      if (length_nbrs >= ord) {
        if (length_nbrs > ord) 
          done_xy <- FALSE
        S <- seq_len(ord)
        repeat {
          num_tests_xy <- num_tests_xy + 1
          pval <- indepTest(x, y, nbrs[S], suffStat)
          if (is.na(pval)) 
            pval <- as.numeric(NAdelete)
          if (pMax_xy < pval) 
            pMax_xy <- pval
          if (pval >= alpha) {
            G_xy <- FALSE
            sepset_xy <- nbrs[S]
            break
          }
          else {
            nextSet <- getNextSet(length_nbrs, ord, S)
            if (nextSet$wasLast) 
              break
            S <- nextSet$nextSet
          }
        }
      }
    }
    list(G_xy, sepset_xy, num_tests_xy, pMax_xy, done_xy)
  }
  edge_test <- function(i) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    num_tests_i <- 0
    G_i <- TRUE
    pMax_xy <- pMax[x, y]
    pMax_yx <- pMax[y, x]
    sepset_xy <- NULL
    sepset_yx <- NULL
    done_i <- TRUE
    res_x <- edge_test_xy(x, y)
    G_i <- res_x[[1]]
    sepset_xy <- res_x[[2]]
    num_tests_i <- num_tests_i + res_x[[3]]
    pMax_xy <- res_x[[4]]
    done_i <- done_i & res_x[[5]]
    if (G_i) {
      if (ord == 0) {
        num_tests_i <- num_tests_i + 1
      }
      else {
        res_y <- edge_test_xy(y, x)
        G_i <- res_y[[1]]
        sepset_yx <- res_y[[2]]
        num_tests_i <- num_tests_i + res_y[[3]]
        pMax_yx <- res_y[[4]]
        done_i <- done_i & res_y[[5]]
      }
    }
    list(i, G_i, sepset_xy, sepset_yx, num_tests_i, pMax_xy, 
         pMax_yx, done_i)
  }
  edge_tests <- function(l) {
    res <- vector("list", length(l))
    for (k in 1:length(l)) {
      res[[k]] <- edge_test(l[[k]])
    }
    res
  }
  total_mem <- function() {
    tryCatch({
      if (Sys.info()[["sysname"]] == "Linux") {
        total <- (as.numeric(gsub("[^0-9]*([0-9]*)", 
                                  "\\1", system("egrep '^MemFree:' /proc/meminfo", 
                                                intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                   "\\1", system("egrep '^Cached:' /proc/meminfo", 
                                                                                                 intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                                                                    "\\1", system("egrep '^Inactive:' /proc/meminfo", 
                                                                                                                                                  intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                                                                                                                     "\\1", system("egrep '^Buffers:' /proc/meminfo", 
                                                                                                                                                                                                   intern = TRUE))))/1000
        return(total)
      }
      else if (Sys.info()[["sysname"]] == "Windows") {
        total <- (as.numeric(gsub("[^0-9]*([0-9]*)", 
                                  "\\1", system("wmic OS get FreePhysicalMemory /Value", 
                                                intern = TRUE))[3]))/1000
        return(total)
      }
      else if (Sys.info()[["sysname"]] == "Darwin") {
        total <- 4096 * (as.numeric(gsub("[^0-9]*([0-9]*)", 
                                         "\\1", system("vm_stat | grep 'Pages free'", 
                                                       intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                          "\\1", system("vm_stat | grep 'Pages inactive'", 
                                                                                                        intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                                                                           "\\1", system("vm_stat | grep 'Pages speculative'", 
                                                                                                                                                         intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                                                                                                                            "\\1", system("vm_stat | grep 'Pages purgeable'", 
                                                                                                                                                                                                          intern = TRUE))))/1e+06
        return(total)
      }
      else {
        total <- (as.numeric(gsub("[^0-9]*([0-9]*)", 
                                  "\\1", system("vmstat -s | grep 'free memory'", 
                                                intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                   "\\1", system("vmstat -s | grep 'inactive memory'", 
                                                                                                 intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                                                                    "\\1", system("vmstat -s | grep 'buffer memory'", 
                                                                                                                                                  intern = TRUE))) + as.numeric(gsub("[^0-9]*([0-9]*)", 
                                                                                                                                                                                     "\\1", system("vmstat -s | grep 'swap cache'", 
                                                                                                                                                                                                   intern = TRUE))))/1000
        return(total)
      }
    }, error = function(e) {
      return(1024)
    }, warning = function(e) {
      return(1024)
    })
  }
  parallel_threshold <- 100
  if (mem.efficient) {
    mem_per_test <- 2
    tests_per_batch <- as.integer(total_mem()/mem_per_test)
  }
  start_total <- proc.time()
  while (!done && any(G) && ord <= m.max) {
    n.edgetests[ord1 <- ord + 1L] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[, 1]), ]
    ind <- subset(ind, ind[, 1] < ind[, 2])
    remainingEdgeTests <- nrow(ind)
    G.l <- split(G, gl(p, p))
    if (!mem.efficient) {
      tests_per_batch <- remainingEdgeTests
    }
    
    if (Sys.info()[["sysname"]] == "Linux") {
      workers <- makeCluster(num_workers, type = "FORK")
    }
    else{
      if(!require('doSNOW')){
        stop('Package doSNOW is required in Window! (function: skeleton_parallel)')
      }
      workers <- makePSOCKcluster(num_workers)
      registerDoSNOW(workers)
      clusterExport(workers, list = c('indepTest', 'alpha', 'bicTest', 'singleBicTest'), envir = environment());
      eval(suffStat)
      clusterEvalQ(workers, library(pcalg))
      #clusterEvalQ(workers, source('independenceTest.R'))
    }
    
    for (j in seq(1, remainingEdgeTests, by = tests_per_batch)) {
      l <- min(remainingEdgeTests, j + tests_per_batch - 
                 1)
      if (l - j + 1 < num_workers) {
        num_workers <- l - j + 1
      }
      res <- NULL
      if (l - j + 1 < parallel_threshold) {
        res <- lapply(j:l, edge_test)
      }
      else{
        res <- do.call("c", clusterApply(workers, clusterSplit(workers, j:l), edge_tests))
      }
      
      for (p_obj in res) {
        i <- p_obj[[1]]
        x <- ind[i, 1]
        y <- ind[i, 2]
        n.edgetests[ord1] <- n.edgetests[ord1] + p_obj[[5]]
        pMax[x, y] <- p_obj[[6]]
        pMax[y, x] <- p_obj[[7]]
        G[x, y] <- G[y, x] <- p_obj[[2]]
        if (!p_obj[[2]]) {
          if (!is.null(p_obj[[3]])) 
            sepset[[x]][[y]] <- p_obj[[3]]
          if (!is.null(p_obj[[4]])) 
            sepset[[y]][[x]] <- p_obj[[4]]
        }
        done <- done & p_obj[[8]]
      }
    }
    ord <- ord + 1
    stopCluster(workers)
  }
  total_t = proc.time() - start_total
  #cat("n=", suffStat[[2]], ",p=", p, "\n", sep = "")
  #cat("Num CI Tests=", n.edgetests, ",Total CI Tests=", sum(unlist(n.edgetests)), 
  #    ",Total Time=", total_t[3], "\n", sep = " ")
  for (i in 1:(p - 1)) {
    for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i, 
                                                        j], pMax[j, i])
  }
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}