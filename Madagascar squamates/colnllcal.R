colnllcal <- function(pars,tot,dt) {
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
	dtcal2 <- function (state,pars,dt,t,t1) {
		if (t1>=t) {
			pp <- numeric(2)
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
					pp[i+2] <- pp[i+2]+log(sum(tmp[2,2:4]))
					state2 <- as.numeric(c(tmp[2,2:4]/sum(tmp[2,2:4]),tmp[2,5:7]))
					tmp <- ode(c(nA=state2[1],nB=state2[2],nAB=state2[3],eA=state2[4],eB=state2[5],eAB=state2[6]),c(i*dt,t),func,pars)
					pp[i+2] <- pp[i+2]+log(sum(tmp[2,2:4]))
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
	RES <- dtcal2(c(402/403,0,1/403,1,1,1),pars,dt,tot,0)
	prob <- exp(RES[[1]])*RES[[2]][,2]
	prob <- prob[-length(prob)]-prob[-1]
	prob[prob<0] <- 0
	prob <- prob/sum(prob)
}