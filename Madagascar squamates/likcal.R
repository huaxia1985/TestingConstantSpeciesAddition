likcal <- function (pars,edge,category,root,is.AB,samplet,samplef,samplem) {
	lg <- 0
	func <- function (time,state,par) {
		with(as.list(c(state,par)), {
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
if (category==1) {
	n <- dim(edge)[1]
	prob <- matrix(NA,n,6)
 	for (i in n:2) {
		if (edge[i,6]<=1 && edge[i,7]<=1) {
			ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1-samplef,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
			ext <- ext[2,5:7]
			ext <- ext+(1-ext)*(1-samplem)
			if (edge[i,6]==0 && edge[i,7]==0) {
				if (edge[i,8]==0||edge[i,9]==0) {
					res1 <- tryCatch(ode(c(nA=0,nB=2/sum(edge[i,8:9]),nAB=0,eA=0,eB=1-2/sum(edge[i,8:9]),eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2+log(samplem)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],ext))
				} else if (length(samplet)>1) {
					res1 <- tryCatch(ode(c(nA=0,nB=samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2+log(samplem)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],res1[2,5:7]))
				} else {
					res1 <- tryCatch(ode(c(nA=0,nB=1/edge[i,8],nAB=0,eA=0,eB=1-1/edge[i,8],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=0,nB=1/edge[i,9],nAB=0,eA=0,eB=1-1/edge[i,9],eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)*2
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,ext))
					}
			}
			else if (edge[i,6]==1 && edge[i,7]==1) {
				if (edge[i,8]==0||edge[i,9]==0) {
					res1 <- tryCatch(ode(c(nA=2/sum(edge[i,8:9]),nB=0,nAB=0,eA=1-2/sum(edge[i,8:9]),eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2+ifelse(is.element(46,edge[i,2:3]),1,2)*log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],ext))
				} else if (length(samplet)>1) {
					res1 <- tryCatch(ode(c(nA=samplet[1],nB=0,nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],res1[2,5:7]))
				} else {
					res1 <- tryCatch(ode(c(nA=1/edge[i,8],nB=0,nAB=0,eA=1-1/edge[i,8],eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=1/edge[i,9],nB=0,nAB=0,eA=1-1/edge[i,9],eB=0,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,ext))
				}
			}
			else {
				if (length(samplet)==1) {
					res1 <- tryCatch(ode(c(nA=as.numeric((edge[i,6]==1)/edge[i,8]),nB=as.numeric((edge[i,6]==0)/edge[i,8]),nAB=0,eA=as.numeric(ifelse(edge[i,6]==1,1-1/edge[i,8],0)),eB=as.numeric(ifelse(edge[i,6]==0,1-1/edge[i,8],0)),eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=as.numeric((edge[i,7]==1)/edge[i,9]),nB=as.numeric((edge[i,7]==0)/edge[i,9]),nAB=0,eA=as.numeric(ifelse(edge[i,7]==1,1-1/edge[i,9],0)),eB=as.numeric(ifelse(edge[i,7]==0,1-1/edge[i,9],0)),eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,ext))
				} else {
					res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)*samplet[1],nB=(edge[i,6]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=(edge[i,7]==1)*samplet[1],nB=(edge[i,7]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
				}
				}
		}
		else {
			if (edge[i,6]<=1) {
				if (length(samplet)==1) {
					if (i<=2) {samplem<-1}
					res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)/edge[i,8],nB=(edge[i,6]==0)/edge[i,8],nAB=0,eA=(edge[i,6]==1)*(1-1/edge[i,8]),eB=(edge[i,6]==0)*(1-1/edge[i,8]),eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
				} else {
					res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)*samplet[1],nB=(edge[i,6]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
				}
				res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res2[2,5:7]))
				if (i==3&&edge[i,2]==21) {
					ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					prob[i,4:6] <- ext[2,5:7]
				}
			}
			else if (edge[i,7]<=1) {
				if (length(samplet)==1) {
					if (i<=2) {samplem<-1}
					res1 <- tryCatch(ode(c(nA=(edge[i,7]==1)/edge[i,9],nB=(edge[i,7]==0)/edge[i,9],nAB=0,eA=(edge[i,7]==1)*(1-1/edge[i,9]),eB=(edge[i,7]==0)*(1-1/edge[i,9]),eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
				} else {
					res1 <- tryCatch(ode(c(nA=(edge[i,7]==1)*samplet[1],nB=(edge[i,7]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
				}
				res2 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res2[2,5:7]))
				if (i==3&&edge[i,3]==21) {
					ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					prob[i,4:6] <- ext[2,5:7]
				}
			}
			else {
				res1 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res2[2,5:7]))
				}
			}
			}
		}
	lik <- ifelse (root==1,sum(prob[i,2:3]),prob[i,2])
if (category==2) {
	n <- dim(edge)[1]
	prob <- matrix(NA,n,6)
 	for (i in n:2) {
		if (edge[i,6]<=1 && edge[i,7]<=1) {
			ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1-samplef,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
			ext <- as.numeric(ext[2,5:7])
			ext <- ext+(1-ext)*(1-samplem)
			if (edge[i,6]==0 && edge[i,7]==0) {
				if (edge[i,8]==0||edge[i,9]==0) {
					res1 <- tryCatch(ode(c(nA=0,nB=2/sum(edge[i,8:9]),nAB=0,eA=0,eB=1-2/sum(edge[i,8:9]),eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2+log(samplem)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],ext))
				} else if (length(samplet)>1) {
					res1 <- tryCatch(ode(c(nA=0,nB=samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],res1[2,5:7]))
				} else {
					res1 <- tryCatch(ode(c(nA=0,nB=1/edge[i,8],nAB=0,eA=0,eB=1-1/edge[i,8],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=0,nB=1/edge[i,9],nAB=0,eA=0,eB=1-1/edge[i,9],eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,ext))
					}
			}
			else if (edge[i,6]==1 && edge[i,7]==1) {
				if (edge[i,8]==0||edge[i,9]==0) {
					res1 <- tryCatch(ode(c(nA=2/sum(edge[i,8:9]),nB=0,nAB=0,eA=1-2/sum(edge[i,8:9]),eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2+log(samplem)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],ext))
					if (is.element(5,edge[i,2:3])) {
						prob[i,4:6] <- res1[2,5:7]
					}
				} else if (length(samplet)>1) {
					res1 <- tryCatch(ode(c(nA=samplet[1],nB=0,nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)*2
					res1[2,2:4] <- res1[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]^2*pars[1],res1[2,3]^2*pars[2],res1[2,2]*res1[2,4]*pars[1]+res1[2,4]*res1[2,3]*pars[2]+res1[2,2]*res1[2,3]*pars[3],res1[2,5:7]))
				} else {
					res1 <- tryCatch(ode(c(nA=1/edge[i,8],nB=0,nAB=0,eA=1-1/edge[i,8],eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=1/edge[i,9],nB=0,nAB=0,eA=1-1/edge[i,9],eB=0,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,ext))
				}
			}
			else {
				if (length(samplet)==1) {
					res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)/edge[i,8],nB=(edge[i,6]==0)/edge[i,8],nAB=0,eA=(edge[i,6]==1)*(1-1/edge[i,8]),eB=(edge[i,6]==0)*(1-1/edge[i,8]),eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=(edge[i,7]==1)/edge[i,9],nB=(edge[i,7]==0)/edge[i,9],nAB=0,eA=(edge[i,7]==1)*(1-1/edge[i,9]),eB=(edge[i,7]==0)*(1-1/edge[i,9]),eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,ext))
				} else {
					res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)*samplet[1],nB=(edge[i,6]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)
					res1[2,2:4] <- res1[2,2:4]/tmp
					res2 <- tryCatch(ode(c(nA=(edge[i,7]==1)*samplet[1],nB=(edge[i,7]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res2[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)
					res2[2,2:4] <- res2[2,2:4]/tmp
					prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res1[2,5:7]))
				}
				}
		}
		else {
			if (edge[i,6]<=1) {
				if (length(samplet)==1) {
					if (i<=3) {samplem<-1}
					if (!is.AB) {res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)/edge[i,8],nB=(edge[i,6]==0)/edge[i,8],nAB=0,eA=(edge[i,6]==1)*(1-1/edge[i,8]),eB=(edge[i,6]==0)*(1-1/edge[i,8]),eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)} 
					else {res1 <- tryCatch(ode(c(nA=0,nB=0,nAB=1/2,eA=0,eB=1,eAB=1/2),c(0,edge[i,4]),func,pars),warning=function(w) w)}
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
				} else {
					res1 <- tryCatch(ode(c(nA=(edge[i,6]==1)*samplet[1],nB=(edge[i,6]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)
					res1[2,2:4] <- res1[2,2:4]/tmp
				}
				res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=ifelse(is.AB,1,prob[edge[i,7],5]),eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res2[2,5:7]))
				if (i==3&&edge[i,2]==17) {
					ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,4]),func,pars),warning=function(w) w)
					prob[i,4:6] <- ext[2,5:7]
				}
			}
			else if (edge[i,7]<=1) {
				if (length(samplet)==1) {
					if (i<=3) {samplem<-1}
					if (!is.AB) {res1 <- tryCatch(ode(c(nA=(edge[i,7]==1)/edge[i,9],nB=(edge[i,7]==0)/edge[i,9],nAB=0,eA=(edge[i,7]==1)*(1-1/edge[i,9]),eB=(edge[i,7]==0)*(1-1/edge[i,9]),eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)}
					else {res1 <- tryCatch(ode(c(nA=0,nB=0,nAB=1/2,eA=0,eB=1/2,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)}
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)+log(samplem)
					res1[2,2:4] <- res1[2,2:4]/tmp
				} else {
					res1 <- tryCatch(ode(c(nA=(edge[i,7]==1)*samplet[1],nB=(edge[i,7]==0)*samplet[2],nAB=0,eA=1-samplet[1],eB=1-samplet[2],eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					tmp <- try(sum(res1[2,2:4]),silent=T)
					if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
					lg <- lg+log(tmp)
					res1[2,2:4] <- res1[2,2:4]/tmp
				}
				res2 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=ifelse(is.AB,1,prob[edge[i,6],5]),eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res2[2,5:7]))
				if (i==3&&!is.AB&&length(samplet)==1) {
					ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,edge[i,5]),func,pars),warning=function(w) w)
					prob[i,4:6] <- ext[2,5:7]
				}
			} 
			else {
				res1 <- tryCatch(ode(c(nA=prob[edge[i,6],1],nB=prob[edge[i,6],2],nAB=prob[edge[i,6],3],eA=prob[edge[i,6],4],eB=prob[edge[i,6],5],eAB=prob[edge[i,6],6]),c(0,edge[i,4]),func,pars),warning=function(w) w)
				tmp <- try(sum(res1[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res1[2,2:4] <- res1[2,2:4]/tmp
				res2 <- tryCatch(ode(c(nA=prob[edge[i,7],1],nB=prob[edge[i,7],2],nAB=prob[edge[i,7],3],eA=prob[edge[i,7],4],eB=prob[edge[i,7],5],eAB=prob[edge[i,7],6]),c(0,edge[i,5]),func,pars),warning=function(w) w)
				tmp <- try(sum(res2[2,2:4]),silent=T)
				if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
				lg <- lg+log(tmp)
				res2[2,2:4] <- res2[2,2:4]/tmp
				prob[i,] <- as.numeric(c(res1[2,2]*res2[2,2]*pars[1],res1[2,3]*res2[2,3]*pars[2],(res1[2,2]*res2[2,4]+res2[2,2]*res1[2,4])*pars[1]/2+(res1[2,4]*res2[2,3]+res2[2,4]*res1[2,3])*pars[2]/2+(res1[2,2]*res2[2,3]+res2[2,2]*res1[2,3])*pars[3]/2,res2[2,5:7]))
				#if (length(samplet)>1&&i==2) {
				#	tmp <- is.element(edge[,2],c(5:6))
				#	tmp <- edge[tmp,4]+edge[2,which(edge[2,2:3]==edge[tmp,1])+3]
				#	ext <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=0,eAB=0),c(0,tmp),func,pars),warning=function(w) w)
				#	prob[i,4:6] <- ext[2,5:7]
				#}
			}
			}
		}
	res <- tryCatch(ode(c(nA=prob[2,1],nB=prob[2,2],nAB=prob[2,3],eA=prob[2,4],eB=ifelse(is.element(edge[2,2],c(5:6)),1,prob[2,5]),eAB=prob[2,6]),c(0,edge[i,c(4,5)[which(edge[1,6:7]==2)]]),func,pars),silent=T)
	if (inherits(res,"simpleWarning")) {stop("negative probability")}
	lik <- res[2,2]
	}
if (category==3) {
	sample.1 <- as.numeric(ifelse(is.na(edge[2,8]),edge[2,9],edge[2,8]))
	res1 <- tryCatch(ode(c(nA=1/sample.1,nB=0,nAB=0,eA=1-1/sample.1,eB=0,eAB=0),c(0,edge[2,4]),func,pars),silent=T)
	tmp <- try(sum(res1[2,2:4]),silent=T)
	if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
	lg <- lg+log(tmp)
	res1[2,2:4] <- res1[2,2:4]/tmp
	res <- tryCatch(ode(c(nA=as.numeric(res1[2,2]),nB=as.numeric(res1[2,3]),nAB=as.numeric(res1[2,4]),eA=as.numeric(res1[2,5]),eB=1,eAB=as.numeric(res1[2,7])),c(0,edge[1,c(4,5)[which(edge[1,6:7]==2)]]),func,pars),silent=T)
	if (inherits(res,"simpleWarning")) {stop("negative probability")}
	lik <- res[2,3]
}
if (category==4) {
	sample.1 <- as.numeric(ifelse(is.na(edge[2,8]),edge[2,9],edge[2,8]))
	res1 <- tryCatch(ode(c(nA=1/sample.1,nB=0,nAB=0,eA=1-sample.1,eB=1,eAB=0),c(0,edge[2,4]),func,pars),silent=T)
	tmp <- try(sum(res1[2,2:4]),silent=T)
	if (tmp<=0||inherits(tmp,"try-error")) {stop("negative probability")}
	lg <- lg+log(tmp)
	res1[2,2:4] <- res1[2,2:4]/tmp
	tmp <- tryCatch(ode(c(nA=1,nB=0,nAB=0,eA=0,eB=1,eAB=0),c(0,edge[2,4]),func,pars),silent=T)
	res <- tryCatch(ode(c(nA=as.numeric(res1[2,2]),nB=as.numeric(res1[2,3]),nAB=as.numeric(res1[2,4]),eA=as.numeric(tmp[2,5]),eB=1,eAB=as.numeric(tmp[2,7])),c(0,edge[1,c(4,5)[which(edge[1,6:7]==2)]]),func,pars),silent=T)
	if (inherits(res,"simpleWarning")) {stop("negative probability")}
	lik <- res[2,3]
}
log(lik)+lg
}
