probtwo <- function (RES) {
a <- dim(RES[[4]])
ll <- sapply(c(1:a[2]),function (i) max(RES[[4]][,i]))
RES[[2]] <- array(NA,ll+1)
RES[[2]][RES[[4]]+1] <- RES[[1]]
prob <- array(NA,ll+1)
if (ll[2]>2&&ll[3]>2) {
	prob[1,2:(ll[2]-1),2:(ll[3]-1)] <- RES[[2]][1,2:(ll[2]-1),2:(ll[3]-1)]-RES[[2]][1,3:ll[2],3:ll[3]]
}
if (ll[2]>2) {
	prob[1,2:(ll[2]-1),ll[3]] <- RES[[2]][1,2:(ll[2]-1),ll[3]]-RES[[2]][2,3:ll[2],(ll[3]+1)]
}
if (ll[3]>2) {
	prob[1,ll[2],2:(ll[3]-1)] <- RES[[2]][1,ll[2],2:(ll[3]-1)]-RES[[2]][2,(ll[2]+1),3:ll[3]]
}
prob[1,ll[2],ll[3]] <- RES[[2]][1,ll[2],ll[3]]-RES[[2]][2,(ll[2]+1),(ll[3]+1)]
if (ll[1]>1) {
	if (ll[2]>2) {
		prob[2:ll[1],2:(ll[2]-1),ll[3]+1] <- RES[[2]][2:ll[1],2:(ll[2]-1),ll[3]+1]-RES[[2]][3:(ll[1]+1),3:ll[2],(ll[3]+1)]
	}
	if (ll[3]>2) {
		prob[2:ll[1],ll[2]+1,2:(ll[3]-1)] <- RES[[2]][2:ll[1],ll[2]+1,2:(ll[3]-1)]-RES[[2]][3:(ll[1]+1),(ll[2]+1),3:ll[3]]
	}
	prob[2:ll[1],ll[2]+1,ll[3]+1] <- RES[[2]][2:ll[1],ll[2]+1,ll[3]+1]-RES[[2]][3:(ll[1]+1),(ll[2]+1),(ll[3]+1)]
	prob[2:ll[1],ll[2],ll[3]+1] <- RES[[2]][2:ll[1],ll[2],ll[3]+1]-RES[[2]][3:(ll[1]+1),ll[2]+1,ll[3]+1]
	prob[2:ll[1],ll[2]+1,ll[3]] <- RES[[2]][2:ll[1],ll[2]+1,ll[3]]-RES[[2]][3:(ll[1]+1),ll[2]+1,ll[3]+1]
}
if (ll[2]>2) {
	prob[ll[1]+1,2:(ll[2]-1),ll[3]+1] <- RES[[2]][ll[1]+1,2:(ll[2]-1),ll[3]+1]
}
if (ll[3]>2) {
	prob[ll[1]+1,ll[2]+1,2:(ll[3]-1)] <- RES[[2]][ll[1]+1,ll[2]+1,2:(ll[3]-1)]
}
prob[ll[1]+1,ll[2]+1,ll[3]+1] <- RES[[2]][ll[1]+1,ll[2]+1,ll[3]+1]
prob[ll[1]+1,ll[2],ll[3]+1] <- RES[[2]][ll[1]+1,ll[2],ll[3]+1]
prob[ll[1]+1,ll[2]+1,ll[3]] <- RES[[2]][ll[1]+1,ll[2]+1,ll[3]]
prob <- as.numeric(prob[RES[[4]]+1])
prob[prob<0]=0
prob <- prob/sum(prob)
prob
}