library(ape)
library(deSolve)
# load a list called tree that include 300 randomly chosen trees from Crottini et al. 2012
# copy all the functions included in the folder
result <- vector("list",300)
# for the ith tree
for (i in 1:300) {
	dt <- 10 # set time step
	rr <- 0.5 # set prior parameter
	w <- rep(1,7)
	# enter tip numbers of each subtree
	native <- c(3:17,19:21,44:52,55,58,59,62)
	Gerrhosauridae <- cbind(c(18,45,44),c(18,2,17))
	Scinicinae <- cbind(c(7,31,37,43,52),c(36,31,0,53,31))
	Trachylepis <- cbind(c(3,40,51),c(13,60,0))
	Iguania <- cbind(c(9,10,28),c(8,1,32))
	Chamaeleonidae <- cbind(c(14:17,30,41,42),c(44,18,2,29,20,41,17))
	Boidae <- cbind(c(11,13,22,2),c(2,2,53,1))
	Pseudoxyrhophiinae <- cbind(c(47:49,23),c(80,80,80,6))
	Psammophiinae <- cbind(c(29,50),c(50,1))
	Typhopidae <- cbind(c(21,46,62:65),c(1,15,0,124,18,102))
	Phelsuma_Lygo <- cbind(c(5,6,55:59),c(34,0,24,38,38,24,24))
	Homo_Blae <- cbind(c(12,25),c(4,3))
	Paroedura <- cbind(c(4,39),c(19,NA))
	Urop_afro<- cbind(c(20,24),c(19,NA))
	Hemidactylus <- cbind(c(8,19,27),c(1,1,1))
	phy <- tree[[i]]
	tot <- phy$edge.length[1]
	parent <- min(which(phy$edge[,1]==67))
	while (phy$edge[parent,2]>66) {
		tot <- tot+phy$edge.length[parent]
		parent <- min(which(phy$edge[,1]==phy$edge[parent,2]))
	}
	tot <- tot+phy$edge.length[parent]
	N <- length(phy$tip.label)
	# extract subtrees
	Gerrhosauridae <- extractcom(phy,Gerrhosauridae,category=2,N,native)
	Scinicinae <- extractcom(phy,Scinicinae,category=2,N,native)
	Trachylepis <- extractcom(phy,Trachylepis,category=2,N,native)
	Iguania <- extractcom(phy,Iguania,category=2,N,native)
	Chamaeleonidae <- extractcom(phy,Chamaeleonidae,category=2,N,native)
	Boidae <- extractcom(phy,Boidae,category=2,N,native)
	Pseudoxyrhophiinae <- extractcom(phy,Pseudoxyrhophiinae,category=2,N,native)
	Psammophiinae <- extractcom(phy,Psammophiinae,category=2,N,native)
	Typhopidae <- extractcom(phy,Typhopidae,category=2,N,native)
	Phelsuma_Lygo <- extractcom(phy,Phelsuma_Lygo,category=2,N,native)
	Homo_Blae <- extractcom(phy,Homo_Blae,category=2,N,native)
	Paroedura <- extractcom(phy,Paroedura,category=2,N,native)
	Urop_afro <- extractcom(phy,Urop_afro,category=2,N,native)
	Hemidactylus <- extractcom(phy,Hemidactylus,category=2,N,native)
	# function to calculate the posterior probability of the tree
	posterior <- function(pars) {
		out <- try(likcal(pars,Gerrhosauridae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Scinicinae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Trachylepis,category=1,is.AB=F,samplet=NA,samplef=1/60,root=2,samplem=1)+
		likcal(pars,Iguania,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Chamaeleonidae,category=2,is.AB=F,samplet=NA,samplef=1,root=NA,samplem=5/8)+
		likcal(pars,Boidae,category=1,is.AB=F,samplet=NA,samplef=1/53,root=1,samplem=1)+
		likcal(pars,Pseudoxyrhophiinae,category=1,is.AB=F,samplet=c(3/80,1/6),samplef=1,root=1,samplem=1)+
		likcal(pars,Psammophiinae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Typhopidae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=3/4)+
		likcal(pars,Phelsuma_Lygo,category=2,is.AB=F,samplet=c(3/24,2/38),samplef=1,root=NA,samplem=1)+
		likcal(pars,Homo_Blae,category=2,is.AB=F,samplet=NA,samplef=0,root=NA,samplem=1)+
		likcal(pars,Paroedura,category=3,is.AB=F,samplet=NA,samplef=NA,root=NA,samplem=NA)+
		likcal(pars,Urop_afro,category=4,is.AB=F,samplet=NA,samplef=NA,root=NA,samplem=NA)+
		likcal(pars,Hemidactylus,category=2,is.AB=T,samplet=NA,samplef=0,root=NA,samplem=1)+
		sum(dexp(pars, rr, log=TRUE)),silent=T)
		if ((inherits(out, "try-error") || (!is.finite(out))||(is.na(out))))
     	 100000
   	else -out
	}
	# calculate the maximum likelihood estimate as the initial value of MCMC
	x.init <- nlm(posterior,p=c(0.1,0.1,10,0.1,0.1,0.1,0.1),hessian=F)$estimate
	posterior <- function(pars) {
		out <- try(likcal(pars,Gerrhosauridae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Scinicinae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Trachylepis,category=1,is.AB=F,samplet=NA,samplef=1/60,root=2,samplem=1)+
		likcal(pars,Iguania,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Chamaeleonidae,category=2,is.AB=F,samplet=NA,samplef=1,root=NA,samplem=5/8)+
		likcal(pars,Boidae,category=1,is.AB=F,samplet=NA,samplef=1/53,root=1,samplem=1)+
		likcal(pars,Pseudoxyrhophiinae,category=1,is.AB=F,samplet=c(3/80,1/6),samplef=1,root=1,samplem=1)+
		likcal(pars,Psammophiinae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=1)+
		likcal(pars,Typhopidae,category=1,is.AB=F,samplet=NA,samplef=1,root=1,samplem=3/4)+
		likcal(pars,Phelsuma_Lygo,category=2,is.AB=F,samplet=c(3/24,2/38),samplef=1,root=NA,samplem=1)+
		likcal(pars,Homo_Blae,category=2,is.AB=F,samplet=NA,samplef=0,root=NA,samplem=1)+
		likcal(pars,Paroedura,category=3,is.AB=F,samplet=NA,samplef=NA,root=NA,samplem=NA)+
		likcal(pars,Urop_afro,category=4,is.AB=F,samplet=NA,samplef=NA,root=NA,samplem=NA)+
		likcal(pars,Hemidactylus,category=2,is.AB=T,samplet=NA,samplef=0,root=NA,samplem=1)+
		sum(dexp(pars, rr, log=TRUE)),silent=T)
		if ((inherits(out, "try-error") || (!is.finite(out))||(is.na(out))))
     	 -100000
   	else out
	}
	# start MCMC
	y.init <- posterior(x.init)
	if ( y.init == -100000 )
    stop("Starting point must have finite probability")
    output <- matrix(NA,100,length(x.init)+1)
	for ( ii in c(1:100)) {
      for (j in seq_along(x.init)) {
      	z <- y.init - rexp(1)
      	L <- x.init
      	R <- x.init
      	while (L[j]<0)
      	   {u <- runif(1) * w[j]
  			L[j] <- L[j] - u
  			R[j] <- R[j] + (w[j]-u)}
  		while ( posterior(L) > z )
    		{L[j] <- L[j] - w[j]}
  		while ( posterior(R) > z )
    		{R[j] <- R[j] + w[j]}
    	L<-max(0,L[j])
    	R<-R[j]
    	xs <- x.init
    	repeat {
    		xs[j] <- runif(1, L, R)
    		ys <- posterior(xs)
    		if ( ys > z )
      		break
    		if ( xs[j] < x.init[j] )
      		L <- xs[j]
    		else
      		R <- xs[j]
      	}
      	x.init <- xs
      	y.init <- ys
      }
      output[ii,] <- c(x.init,y.init)
  	}
  	for (j in seq_along(x.init)) {
	w[j] <- diff(quantile(output[,j],c(0.025, 0.975)))
	}
	output <- matrix(NA,2000,length(x.init)+1)
	for ( ii in c(1:2000)) {
      for (j in seq_along(x.init)) {
      	z <- y.init - rexp(1)
      	L <- x.init
      	R <- x.init
      	while (L[j]<0)
      	   {u <- runif(1) * w[j]
  			L[j] <- L[j] - u
  			R[j] <- R[j] + (w[j]-u)}
  		while ( posterior(L) > z )
    		{L[j] <- L[j] - w[j]}
  		while ( posterior(R) > z )
    		{R[j] <- R[j] + w[j]}
    	L<-max(0,L[j])
    	R<-R[j]
    	xs <- x.init
    	repeat {
    		xs[j] <- runif(1, L, R)
    		ys <- posterior(xs)
    		if ( ys > z )
      		break
    		if ( xs[j] < x.init[j] )
      		L <- xs[j]
    		else
      		R <- xs[j]
      	}
      	x.init <- xs
      	y.init <- ys
    }
    output[ii,] <- c(x.init,y.init)
	coltest <- matrix(NA,1000,3)
	colprob <- matrix(NA,1000,floor(tot/dt))
	colnllprob <- matrix(NA,1000,floor(tot/dt))
    # for the jth set of parameter values
    for (j in 1001:2000) {
      x.init<-as.numeric(output[j,1:7])
      # calculate null expectations for colonization times
      colnll <- colnllcal(x.init,tot,dt)
      # infer probability distribution of colonization time for the tree
      coltime <- vector("list",14)
      coltime[[1]] <- colcal(x.init,Gerrhosauridae,category=2,is.AB=F,samplet=NA,samplef=1,samplem=1,dt,is.root=F)
      coltime[[2]] <- colcal(x.init,Scinicinae,category=2,is.AB=F,samplet=NA,samplef=1,samplem=1,dt,is.root=F)
      coltime[[3]] <- colcal(x.init,Trachylepis,category=2,is.AB=F,samplet=NA,samplef=1/60,samplem=1,dt,is.root=F)
      coltime[[4]] <- colcal(x.init,Iguania,category=2,is.AB=F,samplet=NA,samplef=1,samplem=1,dt,is.root=F)
      coltime[[5]] <- colcal(x.init,Chamaeleonidae,category=2,is.AB=F,samplet=NA,samplef=1,samplem=5/8,dt,is.root=T)
      coltime[[6]] <- colcal(x.init,Boidae,category=2,is.AB=F,samplet=NA,samplef=1/53,samplem=1,dt,is.root=F)
      coltime[[7]] <- colcal(x.init,Pseudoxyrhophiinae,category=2,is.AB=F,samplet=c(3/80,1/6),samplef=1,samplem=1,dt,is.root=F)
      coltime[[8]] <- colcal(x.init,Psammophiinae,category=3,is.AB=F,samplet=NA,samplef=NA,samplem=1,dt,is.root=F)
      coltime[[9]] <- colcal(x.init,Typhopidae,category=2,is.AB=F,samplet=NA,samplef=1,samplem=3/4,dt,is.root=F)
      coltime[[10]] <- colcal(x.init,Phelsuma_Lygo,category=2,is.AB=F,samplet=c(3/24,2/38),samplef=1,samplem=1,dt,is.root=F)
      coltime[[11]] <- colcal(x.init,Homo_Blae,category=3,is.AB=F,samplet=NA,samplef=0,samplem=1,dt,is.root=F)
      coltime[[12]] <- colcal(x.init,Paroedura,category=3,is.AB=F,samplet=NA,samplef=NA,samplem=NA,dt,is.root=F)
      coltime[[13]] <- colcal(x.init,Urop_afro,category=4,is.AB=F,samplet=NA,samplef=NA,samplem=NA,dt,is.root=F)
      coltime[[14]] <- colcal(x.init,Hemidactylus,category=2,is.AB=T,samplet=NA,samplef=0,samplem=1,dt,is.root=F)
      for (ii in 1:14) {
      coltime[[ii]][[1]][coltime[[ii]][[1]]<0] <- 0
      coltime[[ii]][[1]] <- coltime[[ii]][[1]]/sum(coltime[[ii]][[1]])
      coltime[[ii]][[4]] <- sapply(1:ifelse(is.matrix(coltime[[ii]][[4]]),dim(coltime[[ii]][[4]])[1],length(coltime[[ii]][[4]])),setprob,coltime[[ii]],10)
	}
	prob <- numeric(floor(tot/dt))
	for (ii in c(8,11,12,13)) {
		prob[1:length(coltime[[ii]][[4]])] <- prob[1:length(coltime[[ii]][[4]])]+coltime[[ii]][[1]]
	}
	for (ii in c(1:7,9,10,14)) {
		prob1 <- numeric(length(colnll))
		prob1[coltime[[ii]][[4]]] <- prob1[coltime[[ii]][[4]]]+as.numeric(t(matrix(rep(coltime[[ii]][[1]],dim(coltime[[ii]][[4]])[1]),length(coltime[[ii]][[1]]),dim(coltime[[ii]][[4]])[1])))
		prob <- prob+prob1/sum(prob1)
	}
	prob <- prob/14
	colprob[j,]<- prob
	colnllprob[j,] <- colnll
	# test if the inferred colonization times are random draws from the null expectations
	coltest[j,] <- colsigtest(colnll,prob)
    # results include: parameter values, colonisation time test, inferred colonisation times, null expectations for colonisation times
   	result[[i]] <- list(ParameterList=output,ColonisationResult=coltest,InferredColonisationTime=colprob,NullColonisationTime=colnllprob)
 }
