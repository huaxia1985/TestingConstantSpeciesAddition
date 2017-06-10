colcal <- function(pars,edge,category,is.AB,samplet,samplef,samplem,dt,is.root) {
	func0 <- function (time,state1,par) {
		with(as.list(c(state1,par)), {
			sA <- par[1]
			sB <- par[2]
			sAB <- par[3]
			xA <- par[4]
			xB <- par[5]
			dA <- par[6]
			dB <- par[7]
			deA <- -(sA+dA+xA)*eA+xA+dA*eAB+sA*eA*eA
			deB <- -(sB+dB+xB)*eB+xB+dB*eAB+sB*eB*eB
			deAB <- -(sA+sB+sAB+xA+xB)*eAB+xA*eB+xB*eA+sA*eAB*eA+sB*eAB*eB+sAB*eA*eB
			dnA <- -(sA+dA+xA)*nA+dA*nAB+2*sA*nA*eA
			dnB <- -(sB+dB+xB)*nB+0+2*sB*nB*eB
			dnAB <- -(sA+sB+sAB+xA+xB)*nAB+xA*nB+xB*nA+sA*(eA*nAB+eAB*nA)+sB*(eB*nAB+eAB*nB)+sAB*(eA*nB+eB*nA)
			return(list(c(dnA,dnB,dnAB,deA,deB,deAB)))
		})
	}
	func <- function (time,state1,par) {
		with(as.list(c(state1,par)), {
			sA <- par[1]
			sB <- par[2]
			sAB <- par[3]
			xA <- par[4]
			xB <- par[5]
			dA <- par[6]
			dB <- par[7]
			deA <- -(sA+dA+xA)*eA+xA+dA*eAB+sA*eA*eA
			deB <- -(sB+dB+xB)*eB+xB+dB*eAB+sB*eB*eB
			deAB <- -(sA+sB+sAB+xA+xB)*eAB+xA*eB+xB*eA+sA*eAB*eA+sB*eAB*eB+sAB*eA*eB
			dnA <- -(sA+dA+xA)*nA+dA*nAB+2*sA*nA*eA
			dnB <- -(sB+dB+xB)*nB+dB*nAB+2*sB*nB*eB
			dnAB <- -(sA+sB+sAB+xA+xB)*nAB+xA*nB+xB*nA+sA*(eA*nAB+eAB*nA)+sB*(eB*nAB+eAB*nB)+sAB*(eA*nB+eB*nA)
			return(list(c(dnA,dnB,dnAB,deA,deB,deAB)))
		})
	}
	dtcal <- function (state,pars,dt,t,t1,is.root) {
		if (t1>=t) {
			pp <- numeric(2)
			state1 <- matrix(NA,2,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func0,pars)
			pp[2] <- log(sum(tmp[2,2:4]))
			state1[2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			t1 <- t1-t
		} else {
			Tt <- trunc((t-t1)/dt)
			pp <- numeric(Tt+3)
			state1 <- matrix(NA,Tt+3,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- pp[1]+log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func0,pars)
			pp[Tt+3] <- pp[Tt+3]+log(sum(tmp[2,2:4]))
			state1[Tt+3,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			if (t1>0) {
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t1),func0,pars)
				pp[2:(Tt+2)] <- pp[2:(Tt+2)]+log(sum(tmp[2,2:4]))
				state <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(t1,t),func,pars)
				pp[2] <- pp[2]+log(sum(tmp[2,2:4]))
				state1[2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			}
			if (Tt!=0) {
				for (i in 1:Tt) {
					tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,i*dt),func0,pars)
					pp[i+2] <- pp[i+2]+ifelse(sum(tmp[2,2:4])<=0,-1000,log(sum(tmp[2,2:4])))
					state2 <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
					tmp <- ode(c(nA=state2[1],nB=state2[2],nAB=state2[3],eA=state2[4],eB=state2[5],eAB=state2[6]),c(i*dt,t),func,pars)
					pp[i+2] <- pp[i+2]+ifelse(sum(tmp[2,2:4])<=0,-1000,log(sum(tmp[2,2:4])))
					state1[i+2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				}
			}
				t1 <- dt-(t-Tt*dt-t1)
			if (is.na(state1[2,1])) {
				pp <- pp[-2]
				state1 <- state1[-2,]
			}
		}
		list(pp,state1,t1)
	}
		dtcal2 <- function (state,pars,dt,t,t1) {
		if (t1>=t) {
			pp <- numeric(1)
			state1 <- matrix(NA,1,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			t1 <- t1-t
		} else {
			Tt <- trunc((t-t1)/dt)
			pp <- numeric(Tt+2)
			state1 <- matrix(NA,Tt+2,6)
			tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
			pp[1] <- pp[1]+log(sum(tmp[2,2:4]))
			state1[1,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			if (t1>0) {
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t1),func0,pars)
				pp[2:(Tt+2)] <- pp[2:(Tt+2)]+log(sum(tmp[2,2:4]))
				state <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(t1,t),func,pars)
				pp[2] <- pp[2]+log(sum(tmp[2,2:4]))
				state1[2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
			}
			if (Tt!=0) {
				for (i in 1:Tt) {
					tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,i*dt),func0,pars)
					pp[i+2] <- pp[i+2]+ifelse(sum(tmp[2,2:4])<=0,-1000,log(sum(tmp[2,2:4])))
					state2 <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
					tmp <- ode(c(nA=state2[1],nB=state2[2],nAB=state2[3],eA=state2[4],eB=state2[5],eAB=state2[6]),c(i*dt,t),func,pars)
					pp[i+2] <- pp[i+2]+ifelse(sum(tmp[2,2:4])<=0,-1000,log(sum(tmp[2,2:4])))
					state1[i+2,] <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
				}
			}
				t1 <- dt-(t-Tt*dt-t1)
			if (is.na(state1[2,1])) {
				pp <- pp[-2]
				state1 <- state1[-2,]
			}
		}
		list(pp,state1,t1)
	}
		odecal <- function (state,t,pars) {
		tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
		pp <- sum(tmp[2,2:4])
		state <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
		list(pp,state)
	}
	if (category==2) {
		n <- dim(edge)[1]
		state <- vector("list",n)
		res <- vector("list",n)
		for (i in n:2) {
			if (edge[i,6]<=1 && edge[i,7]<=1) {
				ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1-samplef,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
				ext <- ext[2,5:7]
				ext <- ext+(1-ext)*(1-samplem)
				if (edge[i,6]==0 && edge[i,7]==0) {
				if (edge[i,8]==0||edge[i,9]==0) {
					tmp <- ode(c(nA=0,nB=2/sum(edge[i,8:9]),nAB=0,eA=0,eB=1-2/sum(edge[i,8:9]),eAB=0),c(0,edge[i,4]),func,pars)
					tmp[2,5:7] <- ext
					tmp2 <- tmp
				} else if (length(samplet)>1) {
					tmp <- ode(c(nA=0,nB=samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars)
					tmp2 <- tmp
				} else {
					tmp <- ode(c(nA=0,nB=1/edge[i,8],nAB=0,eA=0,eB=1-1/edge[i,8],eAB=0),c(0,edge[i,4]),func,pars)
					tmp[2,5:7] <- ext
					tmp2 <- ode(c(nA=0,nB=1/edge[i,9],nAB=0,eA=0,eB=1-1/edge[i,9],eAB=0),c(0,edge[i,4]),func,pars)
				}
					state[[i]] <- as.numeric(c(tmp[2,2]*tmp2[2,2]*pars[1],tmp[2,3]*tmp2[2,3]*pars[2],(tmp[2,2]*tmp2[2,4]+tmp2[2,2]*tmp[2,4])*pars[1]/2+(tmp[2,4]*tmp2[2,3]+tmp2[2,4]*tmp[2,3])*pars[2]/2+(tmp[2,2]*tmp2[2,3]+tmp2[2,2]*tmp[2,3])*pars[3]/2,tmp[2,5:7]))
					state[[i]][1:3] <- state[[i]][1:3]/sum(state[[i]][1:3])
					state[[i]] <- matrix(state[[i]],1,6)
				} else if (edge[i,6]==1 && edge[i,7]==1) {
				if (edge[i,8]==0||edge[i,9]==0) {
					tmp <- dtcal(c(2/sum(edge[i,8:9]),0,0,1-2/sum(edge[i,8:9]),0,0),pars,dt,edge[i,4],0)
					tmp[[2]][,4] <- ext[1]
					tmp[[2]][,5] <- ext[2]
					tmp[[2]][,6] <- ext[3]
					tmp2 <- tmp
				} else if (length(samplet)>1) {
					tmp <- dtcal(c(samplet[1],0,0,1-samplet[1],1-samplet[2],0),pars,dt,edge[i,4],0)
					tmp2 <- tmp
				} else {
					tmp <- dtcal(c(1/edge[i,8],0,0,1-1/edge[i,8],0,0),pars,dt,edge[i,4],0)
					tmp[[2]][,4] <- ext[1]
					tmp[[2]][,5] <- ext[2]
					tmp[[2]][,6] <- ext[3]
					tmp2 <- dtcal(c(1/edge[i,9],0,0,1-1/edge[i,9],0,0),pars,dt,edge[i,4],0)
				}
					ll <- length(tmp[[1]])
					idx <- sapply(c(1:ll),function(j) rep(j,ll))
					idx <- cbind(as.numeric(idx),as.numeric(t(idx)))
					sta <- matrix(NA,dim(idx)[1],6)
					sta[,1] <- tmp[[2]][idx[,1],1]*tmp2[[2]][idx[,2],1]*pars[1]
					sta[,2] <- tmp[[2]][idx[,1],2]*tmp2[[2]][idx[,2],2]*pars[2]
					sta[,3] <- (tmp[[2]][idx[,1],3]*tmp2[[2]][idx[,2],1]+tmp[[2]][idx[,1],1]*tmp2[[2]][idx[,2],3])*pars[1]/2+(tmp[[2]][idx[,1],2]*tmp2[[2]][idx[,2],3]+tmp[[2]][idx[,1],3]*tmp2[[2]][idx[,2],2])*pars[2]/2+(tmp[[2]][idx[,1],1]*tmp2[[2]][idx[,2],2]+tmp[[2]][idx[,1],2]*tmp2[[2]][idx[,2],1])*pars[3]/2
					sta[,4] <- tmp[[2]][1,4]
					sta[,5] <- tmp[[2]][1,5]
					sta[,6] <- tmp[[2]][1,6]
					res[[i]] <- vector("list",5)
					res[[i]][[1]] <- tmp[[1]][idx[,1]]+tmp2[[1]][idx[,2]]+log(rowSums(sta[,1:3]))
					sta[,1:3] <- sta[,1:3]/rowSums(sta[,1:3])
					res[[i]][[2]] <- sta
					res[[i]][[3]] <- tmp[[3]]
					res[[i]][[4]] <- idx
					res[[i]][[5]] <- c(seq(ll,ll*(ll-1),ll),(ll*(ll-1)+1):(ll*ll))
				} else {
				if (length(samplet)==1) {
					res[[i]] <- dtcal(c(1/edge[i,which(edge[i,6:7]==1)+7],0,0,1-1/edge[i,which(edge[i,6:7]==1)+7],0,0),pars,dt,edge[i,4],0)
					tmp <- ode(c(nA=0,nB=1/as.numeric(edge[i,which(edge[i,6:7]==0)+7]),nAB=0,eA=0,eB=1-1/as.numeric(edge[i,which(edge[i,6:7]==0)+7]),eAB=0),c(0,edge[i,4]),func,pars)
					res[[i]][[2]][,4] <- ext[1]
					res[[i]][[2]][,5] <- ext[2]
					res[[i]][[2]][,6] <- ext[3]
				} else {
					res[[i]] <- dtcal(c(samplet[1],0,0,1-samplet[1],1-samplet[2],0),pars,dt,edge[i,4],0)
					tmp <- ode(c(nA=0,nB=samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars)
				}
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- res[[i]][[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
				}
			} else if (edge[i,6]==1 || edge[i,7]==1) {
				a <- ifelse (edge[i,6]>1,edge[i,6],edge[i,7])
				if (length(samplet)==1) {
					if (!is.AB) {
						res[[i]] <- dtcal(c(1/ifelse(edge[i,6]>1,edge[i,9],edge[i,8]),0,0,1-1/ifelse(edge[i,6]>1,edge[i,9],edge[i,8]),0,0),pars,dt,ifelse(edge[i,6]>1,edge[i,5],edge[i,4]),0)
					} else {
						res[[i]] <- dtcal(c(0,0,1/2,0,1,1/2),pars,dt,ifelse(edge[i,6]>1,edge[i,5],edge[i,4]),0)
					}
					res[[i]][[2]][,4] <- ext[1]
					res[[i]][[2]][,5] <- ext[2]
					res[[i]][[2]][,6] <- ext[3]
				} else {
					res[[i]] <- dtcal(c(samplet[1],0,0,1-samplet[1],1-samplet[2],0),pars,dt,ifelse(edge[i,6]>1,edge[i,5],edge[i,4]),0)
				}
				if (is.null(res[[a]])) {
					tmp <- state[[a]]
					tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,ifelse(edge[i,6]>1,edge[i,4],edge[i,5])),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- res[[i]][[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					res[[i]][[1]] <- res[[i]][[1]]*rowSums(state[[i]][,1:3])
					state[[i]][,1:3] <- state[[i]][,1:3]/rowSums(state[[i]][,1:3])
				} else if (is.null(state[[a]])) {
					res1 <- apply(matrix(res[[a]][[2]][res[[a]][[5]],],length(res[[a]][[5]]),6),1,dtcal,pars,dt,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),res[[a]][[3]])
					res2 <- apply(matrix(res[[a]][[2]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),6),1,odecal,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),pars)
					RES <- vector("list",5)
					RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[5]]][j]),as.numeric(res[[a]][[1]][res[[a]][[5]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[5]])),function (j) res1[[j]][[1]][z]))))
					RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
					RES[[3]] <- res1[[3]][[3]]
					RES[[4]] <- rbind(cbind(0,matrix(res[[a]][[4]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),dim(res[[a]][[4]])[2])),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],],2,rep,length(res1[[1]][[1]]))))
					RES[[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[5]])+1:length(res[[a]][[5]])
					if (edge[a,6]==0||edge[a,7]==0) {
						RES[[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						RES[[5]] <- which(RES[[4]][,1]==max(RES[[4]][,1]))
					} else {
						if (edge[a,6]>1) {
							if (is.null(res[[edge[a,6]]])) {
								RES[[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
								RES[[5]] <- which(RES[[4]][,1]==max(RES[[4]][,1]))
							}
						}
						if (edge[a,7]>1) {
							if (is.null(res[[edge[a,7]]])) {
								RES[[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						 		RES[[5]] <- which(RES[[4]][,1]==max(RES[[4]][,1]))
							}
						}
					}
					idx <- rbind(cbind(rep(c(1:(length(RES[[1]])-length(RES[[5]]))),length(res[[i]][[1]])),as.numeric(sapply(c(1:length(res[[i]][[1]])),function(j) rep(j,length(RES[[1]])-length(RES[[5]]))))),cbind(rep((length(RES[[1]])-length(RES[[5]])+1):length(RES[[1]]),length(res[[i]][[1]])),as.numeric(sapply(c(1:length(res[[i]][[1]])),function(j) rep(j,length(RES[[5]]))))))
					sta <- matrix(NA,dim(idx)[1],6)
					sta[,1] <- RES[[2]][idx[,1],1]*res[[i]][[2]][idx[,2],1]*pars[1]
					sta[,2] <- RES[[2]][idx[,1],2]*res[[i]][[2]][idx[,2],2]*pars[2]
					sta[,3] <- (RES[[2]][idx[,1],3]*res[[i]][[2]][idx[,2],1]+RES[[2]][idx[,1],1]*res[[i]][[2]][idx[,2],3])*pars[1]/2+(RES[[2]][idx[,1],2]*res[[i]][[2]][idx[,2],3]+RES[[2]][idx[,1],3]*res[[i]][[2]][idx[,2],2])*pars[2]/2+(RES[[2]][idx[,1],1]*res[[i]][[2]][idx[,2],2]+RES[[2]][idx[,1],2]*res[[i]][[2]][idx[,2],1])*pars[3]/2
					sta[,4] <- res[[i]][[2]][1,4]
					sta[,5] <- res[[i]][[2]][1,5]
					sta[,6] <- res[[i]][[2]][1,6]
					tmp <- vector("list",5)
					tmp[[1]] <- RES[[1]][idx[,1]]+res[[i]][[1]][idx[,2]]+log(rowSums(sta[,1:3]))
					sta[,1:3] <- sta[,1:3]/rowSums(sta[,1:3])
					tmp[[2]] <- sta
					tmp[[3]] <- res[[i]][[3]]
					tmp[[4]] <- rbind(cbind(as.numeric(sapply(c(1:length(res[[i]][[1]])),function(j) rep(j,length(RES[[1]])-length(RES[[5]])))),apply(RES[[4]][-RES[[5]],],2,rep,length(res[[i]][[1]]))),cbind(as.numeric(sapply(c(1:length(res[[i]][[1]])),function(j) rep(j,length(RES[[5]])))),apply(RES[[4]][RES[[5]],],2,rep,length(res[[i]][[1]]))))
					tmp[[5]] <- (length(RES[[1]])-length(RES[[5]]))*(length(res[[i]][[1]])-1)+1:(length(RES[[1]])-length(RES[[5]])+length(res[[i]][[1]])*length(RES[[5]]))
					res[[i]] <- tmp
				} else {
					ll <- length(res[[a]][[1]])
					res1 <- dtcal(state[[a]][ll,],pars,dt,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),res[[a]][[3]])
					res2 <- apply(matrix(state[[a]][-ll,],(ll-1),6),1,odecal,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),pars)
					RES <- vector("list",3)
					RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][j]),(res[[a]][[1]][ll]+res1[[1]]))
					RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
					l1 <- length(res[[i]][[1]])
					l2 <- length(RES[[1]])
					idx <- cbind(as.numeric(sapply(c(1:l1),function(j) rep(j,l2))),as.numeric(t(sapply(c(1:l2),function(j) rep(j,l1)))))
					sta <- matrix(NA,dim(idx)[1],6)
					sta[,1] <- res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],1]*pars[1]
					sta[,2] <- res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],2]*pars[2]
					sta[,3] <- (res[[i]][[2]][idx[,1],3]*RES[[2]][idx[,2],1]+res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],3])*pars[1]/2+(res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],3]+res[[i]][[2]][idx[,1],3]*RES[[2]][idx[,2],2])*pars[2]/2+(res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],2]+res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],1])*pars[3]/2
					sta[,4] <- res[[i]][[2]][1,4]
					sta[,5] <- res[[i]][[2]][1,5]
					sta[,6] <- res[[i]][[2]][1,6]
					res[[i]][[1]] <- res[[i]][[1]][idx[,1]]+RES[[1]][idx[,2]]+log(rowSums(sta[,1:3]))
					sta[,1:3] <- sta[,1:3]/rowSums(sta[,1:3])
					res[[i]][[2]] <- sta
					res[[i]][[4]] <- idx
					res[[i]][[5]] <- c(seq(l2,l2*(l1-1),l2),(l2*(l1-1)+1):(l1*l2))
				}
				if (i==3&&(edge[i,2]==17||edge[i,2]==21||edge[i,3]==17||edge[i,3]==21)) {
					ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,ifelse(edge[i,6]==1,edge[i,4],edge[i,5])),func,pars),warning=function(w) w)
					ext <- ext[2,5:7]
					res[[i]][[2]][,4] <- ext[1]
					res[[i]][[2]][,5] <- ext[2]
					res[[i]][[2]][,6] <- ext[3]
				}
			} else if (edge[i,6]==0 || edge[i,7]==0) {
				a <- ifelse (edge[i,6]>1,edge[i,6],edge[i,7])
				if (length(samplet)==1) {
					tmp <- ode(c(nA=0,nB=1/ifelse(edge[i,6]>1,edge[i,9],edge[i,8]),nAB=0,eA=0,eB=1-1/ifelse(edge[i,6]>1,edge[i,9],edge[i,8]),eAB=0),c(0,ifelse(edge[i,6]>1,edge[i,5],edge[i,4])),func,pars)
				} else {
					tmp <- ode(c(nA=samplet[1],nB=samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars)
				}
				tmp <- tmp[2,2:4]/sum(tmp[2,2:4])
				if (is.null(res[[a]])) {
					tmp1 <- state[[a]]
					tmp1 <- ode(c(nA=tmp1[1],nB=tmp1[2],nAB=tmp1[3],eA=tmp1[4],eB=tmp1[5],eAB=tmp1[6]),c(0,ifelse(edge[i,6]>1,edge[i,4],edge[i,5])),func,pars)
					state[[i]] <- matrix(tmp1[2,2:7],1,6)
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
				} else if (is.null(state[[a]])) {
					res1 <- apply(matrix(res[[a]][[2]][res[[a]][[5]],],length(res[[a]][[5]]),6),1,dtcal,pars,dt,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),res[[a]][[3]])
					res2 <- apply(matrix(res[[a]][[2]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),6),1,odecal,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),pars)
					res[[i]] <- vector("list",5)
					res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[5]]][j]),as.numeric(res[[a]][[1]][res[[a]][[5]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[5]])),function (j) res1[[j]][[1]][z]))))
					res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
					res[[i]][[3]] <- res1[[3]][[3]]
					res[[i]][[4]] <- rbind(cbind(0,matrix(res[[a]][[4]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),dim(res[[a]][[4]])[2])),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],],2,rep,length(res1[[1]][[1]]))))
					res[[i]][[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[5]])+1:length(res[[a]][[5]])
					if (edge[a,6]==0||edge[a,7]==0) {
						res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
					} else {
						if (edge[a,6]>1) {
							if (is.null(res[[edge[a,6]]])) {
								res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
							}
						}
						if (edge[a,7]>1) {
							if (is.null(res[[edge[a,7]]])) {
								res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
								res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
							}
						}
					}
					res[[i]][[2]][,3] <- (res[[i]][[2]][,3]*tmp[1]+res[[i]][[2]][,1]*tmp[3])*pars[1]/2+(res[[i]][[2]][,2]*tmp[3]+res[[i]][[2]][,3]*tmp[2])*pars[2]/2+(res[[i]][[2]][,1]*tmp[2]+res[[i]][[2]][,2]*tmp[1])*pars[3]/2
					res[[i]][[2]][,1] <- res[[i]][[2]][,1]*tmp[1]*pars[1]
					res[[i]][[2]][,2] <- res[[i]][[2]][,2]*tmp[2]*pars[2]
				} else {
					ll <- length(res[[a]][[1]])
					res1 <- dtcal(state[[a]][ll,],pars,dt,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),res[[a]][[3]])
					res2 <- apply(matrix(state[[a]][-ll,],(ll-1),6),1,odecal,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),pars)
					res[[i]] <- vector("list",3)
					res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][j]),(res[[a]][[1]][ll]+res1[[1]]))
					res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
					res[[i]][[3]] <- res1[[3]]
					state[[i]] <- res[[i]][[2]]
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
					res[[i]][[2]] <- state[[i]]
				}
			} else {
				if (is.null(res[[edge[i,6]]])&&is.null(res[[edge[i,7]]])) {
					res1 <- state[[edge[i,6]]]
					res1 <- ode(c(nA=res1[1],nB=res1[2],nAB=res1[3],eA=res1[4],eB=res1[5],eAB=res1[6]),c(0,edge[i,4]),func,pars)
					tmp <- state[[edge[i,7]]]
					tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,edge[i,5]),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					state[[i]] <- matrix(as.numeric(c(res1[2,2:4]/sum(res1[2,2:4]),res1[2,5:7])),1,6)
					state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
					state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
					state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
					state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
				} else if (is.null(res[[edge[i,6]]])||is.null(res[[edge[i,7]]])) {
					a <- ifelse (!is.null(res[[edge[i,6]]]),edge[i,6],edge[i,7])
					b <- ifelse (!is.null(res[[edge[i,6]]]),edge[i,7],edge[i,6])
					tmp <- state[[b]]
					tmp <- ode(c(nA=tmp[1],nB=tmp[2],nAB=tmp[3],eA=tmp[4],eB=tmp[5],eAB=tmp[6]),c(0,ifelse(!is.null(res[[edge[i,6]]]),edge[i,5],edge[i,4])),func,pars)
					tmp <- as.numeric(tmp[2,2:4]/sum(tmp[2,2:4]))
					if (is.null(state[[a]])) {
						res1 <- apply(matrix(res[[a]][[2]][res[[a]][[5]],],length(res[[a]][[5]]),6),1,dtcal,pars,dt,ifelse(is.null(state[[edge[i,6]]]),edge[i,4],edge[i,5]),res[[a]][[3]])
						res2 <- apply(matrix(res[[a]][[2]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),6),1,odecal,ifelse(is.null(state[[edge[i,6]]]),edge[i,4],edge[i,5]),pars)
						res[[i]] <- vector("list",5)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[5]]][j]),as.numeric(res[[a]][[1]][res[[a]][[5]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[5]])),function (j) res1[[j]][[1]][z]))))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
						res[[i]][[3]] <- res1[[3]][[3]]
						res[[i]][[4]] <- rbind(cbind(0,matrix(res[[a]][[4]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),dim(res[[a]][[4]])[2])),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[5]])+1:length(res[[a]][[5]])
					if (edge[a,6]==0||edge[a,7]==0) {
						res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
					} else {
						if (edge[a,6]>1) {
							if (is.null(res[[edge[a,6]]])) {
								res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
							}
						}
						if (edge[a,7]>1) {
							if (is.null(res[[edge[a,7]]])) {
								res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
							}
						}
					}
						res[[i]][[2]][,3] <- (res[[i]][[2]][,3]*tmp[1]+res[[i]][[2]][,1]*tmp[3])*pars[1]/2+(res[[i]][[2]][,2]*tmp[3]+res[[i]][[2]][,3]*tmp[2])*pars[2]/2+(res[[i]][[2]][,1]*tmp[2]+res[[i]][[2]][,2]*tmp[1])*pars[3]/2
						res[[i]][[2]][,1] <- res[[i]][[2]][,1]*tmp[1]*pars[1]
						res[[i]][[2]][,2] <- res[[i]][[2]][,2]*tmp[2]*pars[2]
					} else {
						ll <- length(res[[a]][[1]])
						res1 <- dtcal(state[[a]][ll,],pars,dt,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),res[[a]][[3]])
						res2 <- apply(matrix(state[[a]][-ll,],(ll-1),6),1,odecal,ifelse(edge[i,6]>1,edge[i,4],edge[i,5]),pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][j]),(res[[a]][[1]][ll]+res1[[1]]))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						res[[i]][[3]] <- res1[[3]]
						state[[i]] <- res[[i]][[2]]
						state[[i]][,3] <- (state[[i]][,3]*tmp[1]+state[[i]][,1]*tmp[3])*pars[1]/2+(state[[i]][,2]*tmp[3]+state[[i]][,3]*tmp[2])*pars[2]/2+(state[[i]][,1]*tmp[2]+state[[i]][,2]*tmp[1])*pars[3]/2
						state[[i]][,1] <- state[[i]][,1]*tmp[1]*pars[1]
						state[[i]][,2] <- state[[i]][,2]*tmp[2]*pars[2]
						state[[i]][,1:3] <- state[[i]][,1:3]/sum(state[[i]][,1:3])
						res[[i]][[2]] <- state[[i]]
					}
				} else {
					if (!is.null(state[[edge[i,6]]])&&!is.null(state[[edge[i,7]]])) {
						a <- edge[i,6]
						ll <- length(res[[a]][[1]])
						res1 <- dtcal(state[[a]][ll,],pars,dt,edge[i,4],res[[a]][[3]])
						res2 <- apply(matrix(state[[a]][-ll,],(ll-1),6),1,odecal,edge[i,4],pars)
						res[[i]] <- vector("list",3)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][j]),(res[[a]][[1]][ll]+res1[[1]]))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						res[[i]][[3]] <- res1[[3]]
						a <- edge[i,7]
						ll <- length(res[[a]][[1]])
						res1 <- dtcal(state[[a]][ll,],pars,dt,edge[i,5],res[[a]][[3]])
						res2 <- apply(matrix(state[[a]][-ll,],(ll-1),6),1,odecal,edge[i,5],pars)
						RES <- vector("list",2)
						RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][j]),(res[[a]][[1]][ll]+res1[[1]]))
						RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						l1 <- length(res[[i]][[1]])
						l2 <- length(RES[[1]])
						idx <- cbind(as.numeric(sapply(c(1:l1),function(j) rep(j,l2))),as.numeric(t(sapply(c(1:l2),function(j) rep(j,l1)))))
						sta <- matrix(NA,dim(idx)[1],6)
						sta[,1] <- res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],1]*pars[1]
						sta[,2] <- res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],2]*pars[2]
						sta[,3] <- (res[[i]][[2]][idx[,1],3]*RES[[2]][idx[,2],1]+res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],3])*pars[1]/2+(res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],3]+res[[i]][[2]][idx[,1],3]*RES[[2]][idx[,2],2])*pars[2]/2+(res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],2]+res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],1])*pars[3]/2
						sta[,4] <- res[[i]][[2]][1,4]
						sta[,5] <- res[[i]][[2]][1,5]
						sta[,6] <- res[[i]][[2]][1,6]
						res[[i]][[1]] <- res[[i]][[1]][idx[,1]]+RES[[1]][idx[,2]]+log(rowSums(sta[,1:3]))
						sta[,1:3] <- sta[,1:3]/rowSums(sta[,1:3])
						res[[i]][[2]] <- sta
						#res[[i]][[3]] <- res[[i]][[3]]
						res[[i]][[4]] <- idx
						res[[i]][[5]] <- c(seq(l2,l2*(l1-1),l2),(l2*(l1-1)+1):(l1*l2))
					} else if (is.null(state[[edge[i,6]]])&&is.null(state[[edge[i,7]]])) {
						a <- edge[i,6]
						res1 <- apply(matrix(res[[a]][[2]][res[[a]][[5]],],length(res[[a]][[5]]),6),1,dtcal,pars,dt,edge[i,4],res[[a]][[3]])
						res2 <- apply(matrix(res[[a]][[2]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),6),1,odecal,edge[i,4],pars)
						res[[i]] <- vector("list",5)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[5]]][j]),as.numeric(res[[a]][[1]][res[[a]][[5]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[5]])),function (j) res1[[j]][[1]][z]))))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
						res[[i]][[3]] <- res1[[3]][[3]]
						res[[i]][[4]] <- rbind(cbind(0,matrix(res[[a]][[4]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),dim(res[[a]][[4]])[2])),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[5]])+1:length(res[[a]][[5]])
						if (edge[a,6]==0||edge[a,7]==0) {
						res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
					} else {
						if (edge[a,6]>1) {
							if (is.null(res[[edge[a,6]]])) {
								res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
							}
						}
						if (edge[a,7]>1) {
							if (is.null(res[[edge[a,7]]])) {
								res[[i]][[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- which(res[[i]][[4]][,1]==max(res[[i]][[4]][,1]))
							}
						}
					}
						a <- edge[i,7]
						res1 <- apply(matrix(res[[a]][[2]][res[[a]][[5]],],length(res[[a]][[5]]),6),1,dtcal,pars,dt,edge[i,5],res[[a]][[3]])
						res2 <- apply(matrix(res[[a]][[2]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),6),1,odecal,edge[i,5],pars)
						RES <- vector("list",5)
						RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[5]]][j]),as.numeric(res[[a]][[1]][res[[a]][[5]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[5]])),function (j) res1[[j]][[1]][z]))))
						RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
						RES[[3]] <- res1[[3]][[3]]
						RES[[4]] <- rbind(cbind(0,res[[a]][[4]][-res[[a]][[5]],]),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],],2,rep,length(res1[[1]][[1]]))))
						RES[[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[5]])+1:length(res[[a]][[5]])
						if (edge[a,6]==0||edge[a,7]==0) {
						RES[[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						RES[[5]] <- which(RES[[4]][,1]==max(RES[[4]][,1]))
					} else {
						if (edge[a,6]>1) {
							if (is.null(res[[edge[a,6]]])) {
								RES[[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
								RES[[5]] <- which(RES[[4]][,1]==max(RES[[4]][,1]))
							}
						}
						if (edge[a,7]>1) {
							if (is.null(res[[edge[a,7]]])) {
								RES[[4]] <- rbind(res[[a]][[4]][-res[[a]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[a]][[4]][,1])-1),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
						 		RES[[5]] <- which(RES[[4]][,1]==max(RES[[4]][,1]))
							}
						}
					}
						idx <- rbind(cbind(rep(c(1:(length(res[[i]][[1]])-length(res[[i]][[5]]))),(length(RES[[1]])-length(RES[[5]]))),as.numeric(sapply(c(1:(length(RES[[1]])-length(RES[[5]]))),function(j) rep(j,length(res[[i]][[1]])-length(res[[i]][[5]]))))),cbind(rep(c(1:(length(res[[i]][[1]])-length(res[[i]][[5]]))),length(RES[[5]])),as.numeric(sapply(c((length(RES[[1]])-length(RES[[5]])+1):length(RES[[1]])),function(j) rep(j,length(res[[i]][[1]])-length(res[[i]][[5]]))))),cbind(rep((length(res[[i]][[1]])-length(res[[i]][[5]])+1):length(res[[i]][[1]]),length(RES[[1]])-length(RES[[5]])),as.numeric(sapply(c(1:(length(RES[[1]])-length(RES[[5]]))),function(j) rep(j,length(res[[i]][[5]]))))),cbind(rep((length(res[[i]][[1]])-length(res[[i]][[5]])+1):length(res[[i]][[1]]),length(RES[[5]])),as.numeric(sapply(c((length(RES[[1]])-length(RES[[5]])+1):length(RES[[1]])),function(j) rep(j,length(res[[i]][[5]]))))))
						sta <- matrix(NA,dim(idx)[1],6)
						sta[,1] <- res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],1]*pars[1]
						sta[,2] <- res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],2]*pars[2]
						sta[,3] <- (res[[i]][[2]][idx[,1],3]*RES[[2]][idx[,2],1]+res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],3])*pars[1]/2+(res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],3]+res[[i]][[2]][idx[,1],3]*RES[[2]][idx[,2],2])*pars[2]/2+(res[[i]][[2]][idx[,1],1]*RES[[2]][idx[,2],2]+res[[i]][[2]][idx[,1],2]*RES[[2]][idx[,2],1])*pars[3]/2
						sta[,4] <- RES[[2]][1,4]
						sta[,5] <- RES[[2]][1,5]
						sta[,6] <- RES[[2]][1,6]
						RES1 <- vector("list",5)
						RES1[[1]] <- res[[i]][[1]][idx[,1]]+RES[[1]][idx[,2]]+log(rowSums(sta[,1:3]))
						sta[,1:3] <- sta[,1:3]/rowSums(sta[,1:3])
						RES1[[2]] <- sta
						RES1[[3]] <- RES[[3]]
						RES1[[4]] <- cbind(res[[i]][[4]][idx[,1],],RES[[4]][idx[,2],])
						RES1[[5]] <- c(((length(res[[i]][[1]])-length(res[[i]][[5]]))*(length(RES[[1]])-length(RES[[5]]))+1):dim(idx)[1])
						res[[i]] <- RES1
					} else {
						a <- ifelse (is.null(state[[edge[i,6]]]),edge[i,6],edge[i,7])
						b <- ifelse (is.null(state[[edge[i,6]]]),edge[i,7],edge[i,6])
						ll <- length(res[[b]][[1]])
						res1 <- dtcal(state[[b]][ll,],pars,dt,ifelse(is.null(state[[edge[i,6]]]),edge[i,5],edge[i,4]),res[[a]][[3]])
						res2 <- apply(matrix(state[[a]][-ll,],(ll-1),6),1,odecal,ifelse(is.null(state[[edge[i,6]]]),edge[i,5],edge[i,4]),pars)
						tmp <- vector("list",3)
						tmp[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][j]),(res[[a]][[1]][ll]+res1[[1]]))
						tmp[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
						tmp[[3]] <- res1[[3]]
						res1 <- apply(matrix(res[[a]][[2]][res[[a]][[5]],],length(res[[a]][[5]]),6),1,dtcal,pars,dt,ifelse(is.null(state[[edge[i,6]]]),edge[i,4],edge[i,5]),res[[a]][[3]])
						res2 <- apply(matrix(res[[a]][[2]][-res[[a]][[5]],],length(res[[a]][[1]])-length(res[[a]][[5]]),6),1,odecal,ifelse(is.null(state[[edge[i,6]]]),edge[i,4],edge[i,5]),pars)
						res[[i]] <- vector("list",5)
						res[[i]][[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[a]][[1]][-res[[a]][[5]]][j]),as.numeric(res[[a]][[1]][res[[a]][[5]]]+sapply(c(1:length(res1[[1]][[1]])),function (z) sapply(c(1:length(res[[a]][[5]])),function (j) res1[[j]][[1]][z]))))
						res[[i]][[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1[[1]][[1]])), function (j) t(rbind(sapply(c(1:length(res1)),function (i) res1[[i]][[2]][j,]))))))
						res[[i]][[3]] <- res1[[3]][[3]]
						res[[i]][[4]] <- rbind(cbind(0,res[[a]][[4]][-res[[a]][[5]],]),cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])),function (j) rep(j,length(res[[a]][[5]])))),apply(res[[a]][[4]][res[[a]][[5]],],2,rep,length(res1[[1]][[1]]))))
						res[[i]][[5]] <- length(res2)+(length(res1[[1]][[1]])-1)*length(res[[a]][[5]])+1:length(res[[a]][[5]])
						idx <- rbind(cbind(rep(c(1:(length(res[[i]][[1]])-length(res[[i]][[5]]))),length(tmp[[1]])),as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[1]])-length(res[[i]][[5]]))))),cbind(rep((length(res[[i]][[1]])-length(res[[i]][[5]])+1):length(res[[i]][[1]]),length(tmp[[1]])),as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[5]]))))))
						sta <- matrix(NA,dim(idx)[1],6)
						sta[,1] <- res[[i]][[2]][idx[,1],1]*tmp[[2]][idx[,2],1]*pars[1]
						sta[,2] <- res[[i]][[2]][idx[,1],2]*tmp[[2]][idx[,2],2]*pars[2]
						sta[,3] <- (res[[i]][[2]][idx[,1],3]*tmp[[2]][idx[,2],1]+res[[i]][[2]][idx[,1],1]*tmp[[2]][idx[,2],3])*pars[1]/2+(res[[i]][[2]][idx[,1],2]*tmp[[2]][idx[,2],3]+res[[i]][[2]][idx[,1],3]*tmp[[2]][idx[,2],2])*pars[2]/2+(res[[i]][[2]][idx[,1],1]*tmp[[2]][idx[,2],2]+res[[i]][[2]][idx[,1],2]*tmp[[2]][idx[,2],1])*pars[3]/2
						sta[,4] <- tmp[[2]][1,4]
						sta[,5] <- tmp[[2]][1,5]
						sta[,6] <- tmp[[2]][1,6]
						RES <- vector("list",5)
						RES[[1]] <- res[[i]][[1]][idx[,1]]+tmp[[1]][idx[,2]]+log(rowSums(sta[,1:3]))
						sta[,1:3] <- sta[,1:3]/rowSums(sta[,1:3])
						RES[[2]] <- sta
						RES[[3]] <- tmp[[3]]
						RES[[4]] <- rbind(cbind(as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[1]])-length(res[[i]][[5]])))),apply(res[[i]][[4]][-res[[i]][[5]],],2,rep,length(tmp[[1]]))),cbind(as.numeric(sapply(c(1:length(tmp[[1]])),function(j) rep(j,length(res[[i]][[5]])))),apply(res[[i]][[4]][res[[i]][[5]],],2,rep,length(tmp[[1]]))))
						RES[[5]] <- (length(res[[i]][[1]])-length(res[[i]][[5]]))*(length(tmp[[1]])-1)+1:(length(res[[i]][[1]])-length(res[[i]][[5]])+length(tmp[[1]])*length(res[[i]][[5]]))
						res[[i]] <- RES	
					}
				}
			}
		}
		if (!is.root) {res[[2]][[2]][,5] <- 1}
		res1 <- apply(matrix(res[[2]][[2]][res[[2]][[5]],],length(res[[2]][[5]]),6),1,dtcal2,pars,dt,ifelse(edge[1,6]==2,edge[1,4],edge[1,5]),res[[a]][[3]])
		res2 <- apply(matrix(res[[2]][[2]][-res[[2]][[5]],],length(res[[2]][[1]])-length(res[[2]][[5]]),6),1,odecal,ifelse(edge[1,6]==2,edge[1,4],edge[1,5]),pars)
		RES <- vector("list",5)
		RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+res[[2]][[1]][-res[[2]][[5]]][j]),as.numeric(sapply(c(1:length(res1)),function (j) res[[2]][[1]][res[[2]][[5]][j]]+res1[[j]][[1]])))
		RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),Reduce(rbind,lapply(c(1:length(res1)),function (i) as.matrix(res1[[1]][[2]]))))
		RES[[1]] <- exp(RES[[1]])*RES[[2]][,2]
		RES[[4]] <- rbind(cbind(0,res[[2]][[4]][-res[[2]][[5]],]),cbind(rep(c(1:length(res1[[1]][[1]])),length(res[[2]][[5]])),apply(res[[2]][[4]][res[[2]][[5]],],2,function (i) rep(i,rep(length(res1[[1]][[1]]),length(res[[2]][[5]]))))))
		if (edge[2,6]==0||edge[2,7]==0) {
			RES[[4]] <- rbind(res[[2]][[4]][-res[[2]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[2]][[4]][,1])-1),function (j) rep(j,length(res[[2]][[5]])))),apply(res[[2]][[4]][res[[2]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
		} else {
			if (edge[2,6]>1) {
				if (is.null(res[[edge[2,6]]])) {
					RES[[4]] <- rbind(res[[2]][[4]][-res[[2]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[2]][[4]][,1])-1),function (j) rep(j,length(res[[2]][[5]])))),apply(res[[2]][[4]][res[[2]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
				}
			}
			if (edge[2,7]>1) {
				if (is.null(res[[edge[2,7]]])) {
				RES[[4]] <- rbind(res[[2]][[4]][-res[[2]][[5]],],cbind(as.numeric(sapply(c(1:length(res1[[1]][[1]])+max(res[[2]][[4]][,1])-1),function (j) rep(j,length(res[[2]][[5]])))),apply(res[[2]][[4]][res[[2]][[5]],-1],2,rep,length(res1[[1]][[1]]))))
				}
			}
		}
		ll <- dim(RES[[4]])[2]
		Tt <- vector("list",ll)
		A <- c(which(edge[,6]==1),which(edge[,7]==1))
		A <- A[order(A,decreasing=T)]
		B <- which(RES[[4]][1,]==1)
		B<- B[order(B,decreasing=T)]
		b <- B[!duplicated(A)]
		a <- unique(A)
		for (i in 1:length(b)) {
			Tt[[b[i]]] <- ifelse(edge[a[i],6]==1,edge[a[i],4],edge[a[i],5])
			if (sum(A==a[i])>1) {
				Tt[[b[i]-1]] <- Tt[[b[i]]]
				parent <- c(which(edge[,6]==a[i]),which(edge[,7]==a[i]))
				son <- a[i]
				if (parent>1) {
				while (edge[parent,7]==0||edge[parent,6]==0||(ifelse(edge[parent,6]>1,is.null(res[[edge[parent,6]]]),F))||(ifelse(edge[parent,7]>1,is.null(res[[edge[parent,7]]]),F))) {
					Tt[[b[i]-2]] <- cbind(Tt[[b[i]-2]],ifelse(edge[parent,6]==son,edge[parent,4],edge[parent,5]))
					son <- parent
					parent <- c(which(edge[,6]==parent),which(edge[,7]==parent))
					if (parent==1) {break}
				}
				}
				Tt[[b[i]-2]] <- cbind(Tt[[b[i]-2]],ifelse(edge[parent,6]==son,edge[parent,4],edge[parent,5]))
			} else {
				if (ifelse(edge[a[i],6]>1,!is.null(res[[edge[a[i],6]]]),F)||ifelse(edge[a[i],7]>1,!is.null(res[[edge[a[i],7]]]),F)) {
					parent <- c(which(edge[,6]==a[i]),which(edge[,7]==a[i]))
					son <- a[i]
					if (parent>1) {
					while (edge[parent,7]==0||edge[parent,6]==0||(ifelse(edge[parent,6]>1,is.null(res[[edge[parent,6]]]),F))||(ifelse(edge[parent,7]>1,is.null(res[[edge[parent,7]]]),F))) {
						Tt[[b[i]-1]] <- cbind(Tt[[b[i]-1]],ifelse(edge[parent,6]==son,edge[parent,4],edge[parent,5]))
						son <- parent
						parent <- c(which(edge[,6]==parent),which(edge[,7]==parent))
						if (parent==1) {break}
					}
					}
					Tt[[b[i]-1]] <- cbind(Tt[[b[i]-1]],ifelse(edge[parent,6]==son,edge[parent,4],edge[parent,5]))
				} else {
					parent <- c(which(edge[,6]==a[i]),which(edge[,7]==a[i]))
					son <- a[i]
					if (parent>1) {
					while (edge[parent,7]==0||edge[parent,6]==0||(ifelse(edge[parent,6]>1,is.null(res[[edge[parent,6]]]),F))||(ifelse(edge[parent,7]>1,is.null(res[[edge[parent,7]]]),F))) {
						Tt[[b[i]]] <- cbind(Tt[[b[i]]],ifelse(edge[parent,6]==son,edge[parent,4],edge[parent,5]))
						son <- parent
						parent <- c(which(edge[,6]==parent),which(edge[,7]==parent))
						if (parent==1) {break}
					}
					}
					Tt[[b[i]]] <- cbind(Tt[[b[i]]],ifelse(edge[parent,6]==son,edge[parent,4],edge[parent,5]))
				}
			}
		}
		if (length(Tt)==9) {Tt[[1]]<-ifelse(edge[1,6]==2,edge[1,4],edge[1,5])}
		if (sum(edge[,6:7]==1)==2) {prob <- probtwo(RES)}
		if (sum(edge[,6:7]==1)==3) {prob <- probthree(RES)}
		if (sum(edge[,6:7]==1)==4) {prob <- probfour(RES)}
		if (sum(edge[,6:7]==1)==5) {prob <- probfive(RES)}
		a <- dim(RES[[4]])
		ll <- sapply(c(1:a[2]),function (i) max(RES[[4]][,i]))
	}
	if (category==3) {
		sample.1 <- as.numeric(ifelse(is.na(edge[2,8]),edge[2,9],edge[2,8]))
		RES <- dtcal(c(1/sample.1,0,0,1-1/sample.1,0,0),pars,dt,edge[2,4],0)
		if (!is.na(samplef)) {
			tmp <- ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1,eAB=0),c(0,edge[2,4]),func,pars)
			RES[[2]][,4] <- tmp[2,5]
			RES[[2]][,5] <- tmp[2,6]
			RES[[2]][,6] <- tmp[2,7]
		}
		ll <- length(RES[[1]])
		res1 <- dtcal2(RES[[2]][ll,],pars,dt,ifelse(edge[1,6]==2,edge[1,4],edge[1,5]),RES[[3]])
		res2 <- apply(matrix(RES[[2]][-ll,],(ll-1),6),1,odecal,ifelse(edge[1,6]==2,edge[1,4],edge[1,5]),pars)
		RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+RES[[1]][j]),(RES[[1]][ll]+res1[[1]]))
		RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
		RES[[3]] <- res1[[3]]
		RES[[1]] <- RES[[1]]*rowSums(RES[[2]][,1:3])
		Tt <- c(edge[2,4],ifelse(edge[1,6]==2,edge[1,4],edge[1,5]))
		prob <- exp(RES[[1]][1:(length(RES[[1]])-1)])-exp(RES[[1]][2:length(RES[[1]])])
		prob <- c(prob,exp(RES[[1]][length(RES[[1]])]))
		ll <- length(RES[[1]])
		RES[[4]] <- c(1:ll)
	}
	if (category==4) {
		sample.1 <- as.numeric(ifelse(is.na(edge[2,8]),edge[2,9],edge[2,8]))
		RES <- dtcal(c(1/sample.1,0,0,1-1/sample.1,1,0),pars,dt,edge[2,4],0)
		tmp <- ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1,eAB=0),c(0,edge[2,4]),func,pars)
		RES[[2]][,4] <- tmp[2,5]
		RES[[2]][,5] <- tmp[2,6]
		RES[[2]][,6] <- tmp[2,7]
		ll <- length(RES[[1]])
		res1 <- dtcal2(RES[[2]][ll,],pars,dt,ifelse(edge[1,6]==2,edge[1,4],edge[1,5]),RES[[3]])
		res2 <- apply(matrix(RES[[2]][-ll,],(ll-1),6),1,odecal,ifelse(edge[1,6]==2,edge[1,4],edge[1,5]),pars)
		RES[[1]] <- c(sapply(c(1:length(res2)),function (j) log(res2[[j]][[1]])+RES[[1]][j]),(RES[[1]][ll]+res1[[1]]))
		RES[[2]] <- rbind(t(sapply(c(1:length(res2)),function (i) res2[[i]][[2]])),res1[[2]])
		RES[[3]] <- res1[[3]]
		RES[[1]] <- RES[[1]]*rowSums(RES[[2]][,1:3])
		Tt <- c(edge[2,4],ifelse(edge[1,6]==2,edge[1,4],edge[1,5]))
		prob <- exp(RES[[1]][1:(length(RES[[1]])-1)])-exp(RES[[1]][2:length(RES[[1]])])
		prob <- c(prob,exp(RES[[1]][length(RES[[1]])]))
		ll <- length(RES[[1]])
		RES[[4]] <- c(1:ll)
	}
	list(prob,ll,Tt,RES[[4]])
}