#this function calculate the null distribution of tip branch lengths
tipnullcal <- function(states,pars,tot,dt) {
  func0 <- function (time,state1,par) {
    with(as.list(c(state1,par)), {
      sB <- par[1]
      sA <- par[2]
      sAB <- par[3]
      xB <- par[4]
      xA <- par[5]
      dB <- par[6]
      dA <- par[7]
      deA <- -(sA+dA+xA)*eA+xA+dA*eAB+sA*eA*eA
      deB <- -(sB+dB+xB)*eB+xB+dB*eAB+sB*eB*eB
      deAB <- -(sA+sB+sAB+xA+xB)*eAB+xA*eB+xB*eA+sA*eAB*eA+sB*eAB*eB+sAB*eA*eB
      dnA <- -(sA+dA+xA)*nA+dA*nAB+2*0*nA*eA
      dnB <- -(sB+dB+xB)*nB+dB*nAB+2*sB*nB*eB
      dnAB <- -(sA+sB+sAB+xA+xB)*nAB+xA*nB+xB*nA+0*(eA*nAB+eAB*nA)+0*eB*nAB+0*eAB*nB+0*eA*nB+0*eB*nA
      return(list(c(dnA,dnB,dnAB,deA,deB,deAB)))
    })
  }
  func <- function (time,state1,par) {
    with(as.list(c(state1,par)), {
      sB <- par[1]
      sA <- par[2]
      sAB <- par[3]
      xB <- par[4]
      xA <- par[5]
      dB <- par[6]
      dA <- par[7]
      deA <- -(sA+dA+xA)*eA+xA+dA*eAB+sA*eA*eA
      deB <- -(sB+dB+xB)*eB+xB+dB*eAB+sB*eB*eB
      deAB <- -(sA+sB+sAB+xA+xB)*eAB+xA*eB+xB*eA+sA*eAB*eA+sB*eAB*eB+sAB*eA*eB
      dnA <- -(sA+dA+xA)*nA+dA*nAB+2*sA*nA*eA
      dnB <- -(sB+dB+xB)*nB+dB*nAB+2*sB*nB*eB
      dnAB <- -(sA+sB+sAB+xA+xB)*nAB+xA*nB+xB*nA+sA*(eA*nAB+eAB*nA)+sB*(eB*nAB+eAB*nB)+sAB*(eA*nB+eB*nA)
      return(list(c(dnA,dnB,dnAB,deA,deB,deAB)))
    })
  }
  dtcal2 <- function (state,pars,dt,t) {
      Tt <- trunc(t/dt)
      state1 <- matrix(NA,Tt+1,6)
      tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,t),func,pars)
      state1[1,] <- as.numeric(c(tmp[2,2:4],tmp[2,5:7]))
        for (i in 1:Tt) {
          tmp <- ode(c(nA=state[1],nB=state[2],nAB=state[3],eA=state[4],eB=state[5],eAB=state[6]),c(0,i*dt),func0,pars)
          state2 <- as.numeric(c(tmp[2,2:4],tmp[2,5:7]))
          tmp <- ode(c(nA=state2[1],nB=state2[2],nAB=state2[3],eA=state2[4],eB=state2[5],eAB=state2[6]),c(i*dt,t),func,pars)
          state1[i+1,] <- as.numeric(c(tmp[2,2:4],tmp[2,5:7]))
        }
     list(state1)
  }
  RES <- dtcal2(states,pars,dt,tot)
  prob <- RES[[1]][,2]*RES[[1]][,2]/(RES[[1]][,1]+RES[[1]][,2]+RES[[1]][,3])+RES[[1]][,1]*RES[[1]][,1]/(RES[[1]][,1]+RES[[1]][,2]+RES[[1]][,3])+RES[[1]][,3]*RES[[1]][,3]/(RES[[1]][,1]+RES[[1]][,2]+RES[[1]][,3])
  prob <- prob[-length(prob)]-prob[-1]
  prob[prob<0] <- 0
  prob/sum(prob)
}