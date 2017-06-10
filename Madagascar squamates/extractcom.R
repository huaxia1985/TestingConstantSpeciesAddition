extractcom <- function(phy,tips,category,N,native) {
	aa<-order(phy$edge[,1],decreasing=F)
	n<-length(phy$edge[,1])/2
	phy$edge<-phy$edge[aa,]
	phy$edge.length<-phy$edge.length[aa]
	phy <- cbind(phy$edge[seq(from=1,by=2,length.out=n),1],matrix(phy$edge[,2],ncol=2,byrow=T),matrix(phy$edge.length,ncol=2,byrow=T))
	idx <- which(is.element(phy[,c(2,3)],tips[,1]))
	if (category==1) {
		idx <- unique(c(idx[idx<=n],idx[idx>n]-n))
		idx <- c((min(idx)-1):max(idx))
		edge <- phy[idx,]
		if (sum(edge[1,2:3] %in% edge[,1])<2) {
			edge <- edge[-1,]
		}
		n <- dim(edge)[1]
		edge <- cbind(edge,matrix(NA,n,4))
		for (i in n:1) {
			if (edge[i,2]<=N) {
				if (is.element(edge[i,2],native)) {
					edge[i,6] <- 1
				}
				else {
					edge[i,6] <- 0
				}
			}
			else {
				edge[i,6] <- which(edge[,1]==edge[i,2])
			}
			tmp <- which(tips[,1]==edge[i,2])
			if (length(tmp)>0) {edge[i,8] <- tips[tmp,2]}
			if (edge[i,3]<=N) {
				if (is.element(edge[i,3],native)) {
					edge[i,7] <- 1
				}
				else {
					edge[i,7] <- 0
				}
			}
			else {
				edge[i,7] <- which(edge[,1]==edge[i,3])
			}
			tmp <- which(tips[,1]==edge[i,3])
			if (length(tmp)>0) {edge[i,9] <- tips[tmp,2]}
		}
	}
	if (category==2) {
		idx <- unique(c(idx[idx<=n],idx[idx>n]-n))
		idx <- c((min(idx)-1):max(idx))
		edge <- phy[idx,]
		if (sum(edge[1,2:3] %in% edge[,1])==0) {
			edge <- edge[-1,]
			edge<-rbind(phy[c(which(phy[,2]==edge[1]),which(phy[,3]==edge[1])),],edge)
		}
		while (sum(edge[1,2:3] %in% edge[,1])==2) {
			edge<-rbind(phy[c(which(phy[,2]==edge[1]),which(phy[,3]==edge[1])),],edge)
		}
		n <- dim(edge)[1]
		edge <- cbind(edge,matrix(NA,n,4))
		for (i in n:1) {
			if (edge[i,2]<=N) {
				if (is.element(edge[i,2],native)) {
					edge[i,6] <- 1
				}
				else {
					edge[i,6] <- 0
				}
			}
			else {
				tmp <- which(edge[,1]==edge[i,2])
				edge[i,6] <- ifelse(length(tmp)==0,n+1,tmp)
			}
			tmp <- which(tips[,1]==edge[i,2])
			if (length(tmp)>0) {edge[i,8] <- tips[tmp,2]}
			if (edge[i,3]<=N) {
				if (is.element(edge[i,3],native)) {
					edge[i,7] <- 1
				}
				else {
					edge[i,7] <- 0
				}
			}
			else {
				tmp <- which(edge[,1]==edge[i,3])
				edge[i,7] <- ifelse(length(tmp)==0,n+1,tmp)
			}
			tmp <- which(tips[,1]==edge[i,3])
			if (length(tmp)>0) {edge[i,9] <- tips[tmp,2]}
		}
	}
	matrix(edge,dim(edge)[1],9)
}