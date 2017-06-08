#this code simulates phylogenies using gsa to approximate the null distribution of colonization times.
#different states can be defined as follows to use our analytical solution implemented in calnullcal function to derive the null distribution.
#states when only consider endemic taxa
states <- c(1,0,0,1-1/100,1-1/mean(nB),1-1/mean(nAB))
#states when only consider widespread taxa
states <- c(0,0,1,1-1/100,1-1/mean(nB),1-1/mean(nAB))
#states when consider both endemic and widespread taxa
states <- c(100/(100+mean(nAB)),0,mean(nAB)/(100+mean(nAB)),1-1/100,1-1/mean(nB),1-1/mean(nAB))

max.taxa <- 200
max.t <- Inf
par(mfrow=c(3,3))
for (z in c(1:9)) {
pars <- parslist[[z]]
pars.cl <- pars.ge.to.cl(pars)
colontime <- numeric()
treelist <- vector("list",100)
numbsvector <- list()
timesvector <- list()
timeperiodvector <- list()
timeoverallvector <- vector()
sumtimes <- 0
for (i in 1:100) {
info <- make.tree.classe(pars.cl, 3, max.taxa, max.t, x0=sample(3, 1, FALSE, stationary.freq.classe.ev(pars.cl, 3)),Ntaxa=100)
while (is.character(info)) {
  info <- make.tree.classe(pars.cl, 3, max.taxa, max.t, x0=sample(3, 1, FALSE, stationary.freq.classe.ev(pars.cl, 3)),Ntaxa=100)
}
treelist[[i]] <- info
numbs <- info$nl
timeperiods <- which(numbs[,4]==100)
temp <- vector()
	if (length(timeperiods)>0){
		for (i in 1:length(timeperiods)){
			temp <- c(temp, numbs[timeperiods[i]+1,1]-numbs[timeperiods[i],1])
		}	
	}
numbsvector <- c(numbsvector, list(numbs))
timesvector <- c(timesvector, list(temp))
timeperiodvector <- c(timeperiodvector, list(timeperiods))
timeoverallvector <- c(timeoverallvector, sum(temp))
sumtimes <- sumtimes + sum(temp)
}
timemax <- sumtimes/100
nAB <- numeric()
nB <- numeric()
for (j in 1:length(treelist)) {
    timeperiods <- timeperiodvector[[j]]
    if (length(timeperiods) > 0) {
         tree <- treelist[[j]]
         timelength <- timesvector[[j]]
         timeoverall <- timeoverallvector[[j]]
         numbs <- numbsvector[[j]]
         rsamp <- runif(1, min = 0, max = 1)
         if (rsamp <= (timeoverall/timemax)) {
             r <- runif(1, min = 0, max = timeoverall)
             timepassed <- 0
             jtest <- 0
             while (r > timepassed) {
                  jtest <- jtest + 1
                  timepassed <- timepassed + timelength[jtest]
             }
             cutinterval <- timeperiods[jtest]
             cuttime <- numbs[cutinterval, 1] + runif(1, min = 0, max = timelength[jtest])
			 tmp <- abs(cuttime-tree$tlist)
			 tmp <- which(tmp==min(tmp))
			 info <- tree$infolist[[tmp]]
			 info$state <- info$state - 1
			 attr(info,"t") <- tree$tlist[[tmp]]
			 attr(info,"hist") <- tree$histlist[[tmp]]
			 phy <- me.to.ape.bisse(info[-1,], info$state[1])
			 phy <- prune(phy)
			 tmp <- attr(info,"t")-extractcol(phy$hist,phy)
			 colontime <- c(colontime,tmp)
			 nAB <- c(nAB,sum(phy$tip.state==0))
			 nB <- c(nB,sum(phy$tip.state==1))
         }
     }
}    
obs <- hist(colontime,plot=F)
x <- obs$breaks
obs <- obs$density/sum(obs$density)
exp <- colnllcal(states,pars,max(x),x[2]-x[1])
plot(x[-1],exp,ylim=c(0,1))
points(x[-1],obs,pch=4)
}

#this function extract colonization times from simulated phylogenies. 
extractcol <- function (hist,phy) {
	nodename <- unique(hist$name2)
	#when only consider endemic taxa
	tipname <- phy$tip.label[phy$tip.state==2]
	#when only consider widespread taxa
	tipname <- phy$tip.label[phy$tip.state==0]
	#when consider both endemic and widespread taxa
	tipname <- phy$tip.label[phy$tip.state!=1]
	is.tip <- grepl('sp',nodename)
	removetip <- is.element(nodename[is.tip],tipname)
	nodename <- nodename[-which(removetip==F)]
	out <- numeric()
	for (i in nodename) {
		idx <- which(hist$name2==i)
		col <- c(0,which((hist$from[idx]==2)*(hist$to[idx]==1)==1))
		emg <- c(0,which((hist$from[idx]==1)*(hist$to[idx]==2)==1))
		if (max(col)>max(emg)) {
			out <- c(out,hist$t[idx[max(col)]])
		} else {
			out <- c(out,NA)
		}
	}
	names(out) <- nodename
	is.tip <- grepl('sp',nodename)
	node <- nodename[!is.tip]
	library(gtools)
	node <- mixedsort(node,decreasing=T)
	for (j in node) {
		tree <- extract.clade(phy,j)
		tip2 <- intersect(tree$tip.label,tipname)
		if (length(tip2) > 0) {
		existtip2 <- intersect(nodename,tip2)
		misstip2 <- setdiff(tip2,existtip2)
		existtip2 <- existtip2[is.na(out[existtip2])]
		if (length(misstip2)==0 && length(existtip2)==0) {
			out[j] <- NA	
		}
		if (length(misstip2)>0 && !is.na(out[j])) {
			nodename <- c(nodename,misstip2)
			out2 <- as.numeric(replicate(length(misstip2),out[j]))
			names(out2) <- misstip2
			out <- c(out,out2)
		}
		if (length(existtip2)>0 && !is.na(out[j])) {
			out[existtip2] <- out[j]
			out[j] <- NA
		}
		} else {
			out[j] <- NA
		}
	}
	out <- out[!is.na(out)]
	out <- out[out>0]
}


