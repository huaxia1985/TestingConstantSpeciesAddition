#modified tree.geosse function to record all time periods that have Ntaxa number of endemic species in the phylogeny.
make.tree.classe <- function(pars, k, max.taxa=Inf, max.t=Inf, x0,
                             single.lineage=TRUE,Ntaxa) {
  # The other models don't require k, but this function is hidden away,
  # so no worry about passing k in rather than recomputing it.
  
  if ( x0 < 1 || x0 > k )
    stop("x0 must be an integer in [1,k]")
  
  # arrange the parameters in a list with elements:
  #   lambda = lambda_ijk array, mu = mu vector, q = q_ij array, 
  #   nstates = number of states
  pars.list <- inflate.pars.classe(pars, k)
  
  # for drawing samples below, it's nicer to have 0 than NA for the
  # non-applicable speciation rates
  pars.list$lambda[which(is.na(pars.list$lambda))] <- 0
  
  # row i is all states != i, i.e., states that can be transitioned to
  to <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])),
               k, k-1, TRUE)
  
  # pars is a "k x k+1" matrix giving, for a lineage in state row i,
  # the rate at which speciation, extinction, or anagenetic transition
  # to each other state happens.  This approach loses speciation info
  # (retained in pars.list$lambda) and requires an extra sample() call
  # within the speciation "if" below, but it makes the indices less
  # heinous.
  pars <- cbind(rowSums(pars.list$lambda), pars.list$mu, 
                matrix(pars[-seq_len(k*k*(k+1)/2+k)], k, k-1, TRUE))
  # r.i = total rate at which something happens to a lineage in state i
  r.i <- rowSums(pars)
  
  extinct <- FALSE
  split   <- FALSE
  parent <- 0
  n.i <- rep(0, k)    # number of lineages in state i at this time
  len <- 0            # branch lengths
  t <- 0              # time elapsed
  hist <- list()      # history of transitions
  nl <- matrix(c(0,0,1,0),1,4) #number of lineages through time
  infolist <- list()
  histlist <- list()
  tlist <- numeric()
  
  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0] <- 1
    start <- 0
  } else {
    ##states <- rep(x0, 2)
    ##n.taxa <- lineages <- n.i[x0] <- 2
    stop("Nope.")
  }
  
  while ( sum(states[lineages]==3) <= max.taxa && n.taxa > 0 ) {
        
    # When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt
    
    # Stop if it happens too late.
    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }
    
    len[lineages] <- len[lineages] + dt
    
    # What state does the event happen to?
    state <- sample(k, 1, FALSE, r.n/r.tot)
    
    # What lineage with that state gets the event?
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]
    
    # What event happens?  1 = speciation, 2 = extinction, 
    #    type>2 = transition (type & to provide new state)
    type <- sample(k+1, 1, FALSE, pars[state,])
    
    if ( type == 1 ) {                      # Speciation
      if ( sum(states[lineages]==3) == max.taxa )
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      
      # get daughter states from indices of lamda_ijk 
      lam <- pars.list$lambda[state,,]
      j <- sample(k*k, 1, FALSE, lam)
      s.daught <- c((j-1) %% k + 1, (j-1) %/% k + 1)
      states[new.i] <- s.daught
      n.i[state] <- n.i[state] - 1
      n.i[s.daught[1]] <- n.i[s.daught[1]] + 1
      n.i[s.daught[2]] <- n.i[s.daught[2]] + 1
      
      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      n.taxa <- n.taxa + 1
      lineages <- which(!split & !extinct)
      
    } else if ( type == 2 ) {               # Extinction
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1
      
    } else {                                # Transition (anagenetic)
      states[lineage] <- state.new <- to[state, type - 2]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1,-1)
      hist[[length(hist)+1]] <- c(lineage, t, state, state.new)
    }
    if (n.i[3]!=nl[length(nl[,4]),4]) {nl <- rbind(nl,c(t,n.i))}
    if (n.i[3]==Ntaxa) {
    	info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     start=start, state=states, extinct=extinct,
                     split=split)
 		hist1 <- as.data.frame(do.call(rbind, hist))
  		if ( nrow(hist1) == 0 ) {hist1 <- as.data.frame(matrix(NA, 0, 4))}
  		names(hist1) <- c("idx", "t", "from", "to")
  		hist1$x0 <- info$start[match(hist1$idx, info$idx)]
 		hist1$tc <- hist1$t - hist1$x0
     	infolist <- c(infolist,list(info))
    	histlist <- c(histlist,list(hist1))
    	tlist <- c(tlist,t)
    }
  }
  if (n.taxa==0) {
    info <- "extinct"
  } else {
 	info <- list(infolist=infolist,histlist=histlist,tlist=tlist,nl=nl)
  }
  info
}

#Other functions called in the tree.geosse function.
stationary.freq.classe.ev <- function(pars, k) {
  A <- projection.matrix.classe(pars, k)
  ## continuous time, so the dominant eigenvalue is the largest one
  ## return its eigenvector, normalized
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  evA$vectors[,i] / sum(evA$vectors[,i])
}

projection.matrix.classe <- function(pars, k) {
  A <- matrix(0, nrow=k, ncol=k)

  nsum <- k*(k+1)/2
  kseq <- seq_len(k)
  pars.lam <- pars[seq(1, nsum*k)]
  pars.mu <- pars[seq(nsum*k+1, (nsum+1)*k)]
  pars.q <- pars[seq((nsum+1)*k+1, length(pars))]

  ## array indices of lambda's in parameter vector
  idx.lam <- cbind(rep(kseq, each=nsum), rep(rep(kseq, times=seq(k,1,-1)), k),
                   unlist(lapply(kseq, function(i) i:k)))
  ## transpose of matrix indices of q's in parameter vector
  idx.q <- cbind(unlist(lapply(kseq, function(i) (kseq)[-i])), 
                 rep(kseq, each=k-1))

  ## take care of off-diagonal elements
  for (n in seq_len(nsum*k)) {
    ## add this lambda to A[daughter states, parent state]
    ## (separate steps in case the daughter states are the same)
    r <- idx.lam[n,]
    A[r[2], r[1]] <- A[r[2], r[1]] + pars.lam[n]
    A[r[3], r[1]] <- A[r[3], r[1]] + pars.lam[n]
  }
  A[idx.q] <- A[idx.q] + pars.q

  ## fix the diagonal elements
  diag(A) <- 0
  diag(A) <- -colSums(A) + unlist(lapply(kseq, function(i) 
                 sum(pars.lam[seq((i-1)*nsum+1, i*nsum)]) - pars.mu[i]))
  A
}

pars.ge.to.cl <- function(pars.ge)
{
  if (is.null(names(pars.ge)))
    names(pars.ge) <- default.argnames.geosse()
  pars.cl <- rep(0, 27)
  names(pars.cl) <- default.argnames.classe(3)
  pars.cl['lambda222'] <- pars.cl['lambda112'] <- pars.ge['sA']
  pars.cl['lambda333'] <- pars.cl['lambda113'] <- pars.ge['sB']
  pars.cl['lambda123'] <-  pars.ge['sAB']
  pars.cl['mu2'] <- pars.cl['q13'] <- pars.ge['xA']
  pars.cl['mu3'] <- pars.cl['q12'] <- pars.ge['xB']
  pars.cl['q21'] <- pars.ge['dA']
  pars.cl['q31'] <- pars.ge['dB']
  pars.cl
}
