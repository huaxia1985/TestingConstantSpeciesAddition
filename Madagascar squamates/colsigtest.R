colsigtest <- function (nll,obs) {
	medX <- min(which(cumsum(obs)>=0.5))
	pmed <- sum(sapply(c(1:medX),function (i) sapply(c(0:15),function (j) sapply(c(0:15), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(31-j-z)*factorial(31)/factorial(j)/factorial(z)/factorial(31-j-z)))))
	quadX <- min(which(cumsum(obs)>=0.25))
	pquad <- sum(sapply(c(1:quadX),function (i) sapply(c(0:7),function (j) sapply(c(0:23), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(31-j-z)*factorial(31)/factorial(j)/factorial(z)/factorial(31-j-z)))))
	quad3X <- min(which(cumsum(obs)>=0.75))
	pquad3 <- sum(sapply(c(1:quad3X),function (i) sapply(c(0:23),function (j) sapply(c(0:7), function (z) ifelse(i>1,sum(nll[1:(i-1)]),0)^j*(1-sum(nll[1:i]))^z*nll[i]^(31-j-z)*factorial(31)/factorial(j)/factorial(z)/factorial(31-j-z)))))
	c(pmed,pquad,pquad3)
}