setprob <- function (sample,colinput,dt) {
	if (is.matrix(colinput[[4]])) {
		sample <- colinput[[4]][sample,]
	} else {
		sample <- colinput[[4]][sample]
	}
	sampleold <- sample
	Tt <- colinput[[3]]
	ll <- colinput[[2]]
	if (length(ll)==1) {
		if (length(Tt)>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt)) {
				if (Tt[i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[i]+aa*dt
				}
				if (is.element(sample,idx)) {
					if (i!=1) {
					sample <- floor(sum(Tt[1:(i-1)])/dt)+which(idx==sample)
					}
					break
				}
			}
		}
		out <- sample
	}
	if (length(ll)==3) {
	if (sample[2]<ll[2]) {
		if (length(Tt[[2]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[2]])) {
				if (Tt[[2]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[2]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[2]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[2]][i]+aa*dt
				}
				if (is.element(sample[2],idx)) {
					if (i!=1) {
					sample[2] <- floor(sum(Tt[[2]][1:(i-1)])/dt)+which(idx==sample[2])
					}
					break
				}
			}
		}
	} 
	if (sample[3]<ll[3]) {
		if (length(Tt[[3]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[3]])) {
				if (Tt[[3]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[3]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[3]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[3]][i]+aa*dt
				}
				if (is.element(sample[3],idx)) {
					if (i!=1) {
					sample[3] <- floor(sum(Tt[[3]][1:(i-1)])/dt)+which(idx==sample[3])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[2]])<dt,dt-sum(Tt[[2]]),sum(Tt[[2]])-floor(sum(Tt[[2]])/dt)*dt)
	if (sample[2]==ll[2]||sample[3]==ll[3]) {
		if (length(Tt[[1]])>1) {
			b <- 0
			for (i in 1:length(Tt[[1]])) {
				if (Tt[[1]][i]<t1) {
					idx <- b+1
					b<-b+1
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[1]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[1]][i]+aa*dt
				}
				if (is.element(sample[1],idx)) {
					if (i!=1) {
					sample[1] <- floor(sum(Tt[[1]][1:(i-1)])/dt)+which(idx==sample[1])
					}
					break
				}
			}
		}
		sample[1] <- sample[1]+floor(sum(Tt[[2]])/dt)
	} 
	if (sampleold[2]==ll[2]&&sampleold[3]==ll[3]) {
		out<-c(sample[1],sample[1])
	} else if (sampleold[2]<ll[2]&&sampleold[3]==ll[3]) {
		out<-c(sample[1],sample[2])
	} else if (sampleold[2]==ll[2]&&sampleold[3]<ll[3]) {
		out<-c(sample[1],sample[3])
	} else {
		out<-c(sample[2],sample[3])
	}
	}
	if (length(ll)==5) {
	if (sample[4]<ll[4]) {
		if (length(Tt[[4]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[4]])) {
				if (Tt[[4]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[4]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[4]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[4]][i]+aa*dt
				}
				if (is.element(sample[4],idx)) {
					if (i!=1) {
					sample[4] <- floor(sum(Tt[[4]][1:(i-1)])/dt)+which(idx==sample[4])
					}
					break
				}
			}
		}
	}
	if (sample[5]<ll[5]) {
		if (length(Tt[[5]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[5]])) {
				if (Tt[[5]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[5]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[5]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[5]][i]+aa*dt
				}
				if (is.element(sample[5],idx)) {
					if (i!=1) {
					sample[5] <- floor(sum(Tt[[5]][1:(i-1)])/dt)+which(idx==sample[5])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[4]])<dt,dt-sum(Tt[[4]]),sum(Tt[[4]])-floor(sum(Tt[[4]])/dt)*dt)
	if (sample[3]<ll[3]&&sample[3]>0) {
		if (length(Tt[[3]])>1) {
			b <- 0
			for (i in 1:length(Tt[[3]])) {
				if (Tt[[3]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[3]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[3]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[3]][i]+aa*dt
				}
				if (is.element(sample[3],idx)) {
					if (i!=1) {
					sample[3] <- floor(sum(Tt[[3]][1:(i-1)])/dt)+which(idx==sample[3])
					}
					break
				}
			}
		}
		sample[3] <- sample[3]+floor(sum(Tt[[4]])/dt)
	} 
	if (sample[2]<ll[2]) {
		if (length(Tt[[2]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[2]])) {
				if (Tt[[2]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[2]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[2]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[2]][i]+aa*dt
				}
				if (is.element(sample[2],idx)) {
					if (i!=1) {
					sample[2] <- floor(sum(Tt[[2]][1:(i-1)])/dt)+which(idx==sample[2])
					}
					break
				}
			}
		}
	} 
	t1 <- ifelse(sum(Tt[[2]])<dt,dt-sum(Tt[[2]]),sum(Tt[[2]])-floor(sum(Tt[[2]])/dt)*dt)
	if (sample[2]==ll[2]||sample[3]==ll[3]) {
		if (length(Tt[[1]])>1) {
			b <- 0
			for (i in 1:length(Tt[[1]])) {
				if (Tt[[1]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[1]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[1]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[1]][i]+aa*dt
				}
				if (is.element(sample[1],idx)) {
					if (i!=1) {
					sample[1] <- floor(sum(Tt[[1]][1:(i-1)])/dt)+which(idx==sample[1])
					}
					break
				}
			}
		}
		sample[1] <- sample[1]+floor(sum(Tt[[2]])/dt)
	} 
	if (sampleold[2]==ll[2]) {
		if (sampleold[3]==ll[3]) {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[1],sample[1],sample[1])
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[1],sample[4],sample[1])
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				out<-c(sample[1],sample[1],sample[5])
			}
		} else if (sampleold[3]==0) {
			out<-c(sample[1],sample[4],sample[5])
		} else {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[1],sample[3],sample[3])
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[1],sample[4],sample[3])
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				out<-c(sample[1],sample[3],sample[5])
			}
		}
	} else {
		if (sampleold[3]==ll[3]) {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[2],sample[1],sample[1])
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[2],sample[4],sample[1])
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				out<-c(sample[2],sample[1],sample[5])
			}
		} else if (sampleold[3]==0) {
			out<-c(sample[2],sample[4],sample[5])
		} else {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[2],sample[3],sample[3])
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				out<-c(sample[2],sample[4],sample[3])
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				out<-c(sample[2],sample[3],sample[5])
			}
		}
	}
	}
	if (length(ll)==7) {
	if (sample[6]<ll[6]) {
		if (length(Tt[[6]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[6]])) {
				if (Tt[[6]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[6]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[6]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[6]][i]+aa*dt
				}
				if (is.element(sample[6],idx)) {
					if (i!=1) {
					sample[6] <- floor(sum(Tt[[6]][1:(i-1)])/dt)+which(idx==sample[6])
					}
					break
				}
			}
		}
	}
	if (sample[7]<ll[7]) {
		if (length(Tt[[7]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[7]])) {
				if (Tt[[7]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1<-t1-Tt[[7]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[7]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[7]][i]+aa*dt
				}
				if (is.element(sample[7],idx)) {
					if (i!=1) {
					sample[7] <- floor(sum(Tt[[7]][1:(i-1)])/dt)+which(idx==sample[7])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[6]])<dt,dt-sum(Tt[[6]]),sum(Tt[[6]])-floor(sum(Tt[[6]])/dt)*dt)
	if (sample[5]<ll[5]&&sample[5]>0) {
		if (length(Tt[[5]])>1) {
			b <- 0
			for (i in 1:length(Tt[[5]])) {
				if (Tt[[5]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[5]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[5]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[5]][i]+aa*dt
				}
				if (is.element(sample[5],idx)) {
					if (i!=1) {
					sample[5] <- floor(sum(Tt[[5]][1:(i-1)])/dt)+which(idx==sample[5])
					}
					break
				}
			}
		}
		sample[5] <- sample[5]+floor(sum(Tt[[6]])/dt)
	} 
	if (sample[4]<ll[4]) {
		if (length(Tt[[4]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[4]])) {
				if (Tt[[4]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[4]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[4]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[4]][i]+aa*dt
				}
				if (is.element(sample[4],idx)) {
					if (i!=1) {
					sample[4] <- floor(sum(Tt[[4]][1:(i-1)])/dt)+which(idx==sample[4])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[4]])<dt,dt-sum(Tt[[4]]),sum(Tt[[4]])-floor(sum(Tt[[4]])/dt)*dt)
	if (sample[3]<ll[3]&&sample[3]>0) {
		if (length(Tt[[3]])>1) {
			b <- 0
			for (i in 1:length(Tt[[3]])) {
				if (Tt[[3]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[3]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[3]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[3]][i]+aa*dt
				}
				if (is.element(sample[3],idx)) {
					if (i!=1) {
					sample[3] <- floor(sum(Tt[[3]][1:(i-1)])/dt)+which(idx==sample[3])
					}
					break
				}
			}
		}
		sample[3] <- sample[3]+floor(sum(Tt[[4]])/dt)
	}
	if (sample[2]<ll[2]) {
		if (length(Tt[[2]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[2]])) {
				if (Tt[[2]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[2]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[2]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[2]][i]+aa*dt
				}
				if (is.element(sample[2],idx)) {
					if (i!=1) {
					sample[2] <- floor(sum(Tt[[2]][1:(i-1)])/dt)+which(idx==sample[2])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[2]])<dt,dt-sum(Tt[[2]]),sum(Tt[[2]])-floor(sum(Tt[[2]])/dt)*dt)
	if (sample[2]==ll[2]||sample[3]==ll[3]) {
		if (length(Tt[[1]])>1) {
			b <- 0
			for (i in 1:length(Tt[[1]])) {
				if (Tt[[1]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[1]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[1]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[1]][i]+aa*dt
				}
				if (is.element(sample[1],idx)) {
					if (i!=1) {
					sample[1] <- floor(sum(Tt[[1]][1:(i-1)])/dt)+which(idx==sample[1])
					}
					break
				}
			}
		}
		sample[1] <- sample[1]+floor(sum(Tt[[2]])/dt)
	} 
	if (sampleold[2]==ll[2]) {
		if (sampleold[3]==ll[3]) {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[1],sample[1],sample[1])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[1],sample[6],sample[1])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[1],sample[1],sample[7])
				} else {
					out<-c(sample[1],sample[1],sample[6],sample[7])
				}
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[4],sample[1],sample[1])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[4],sample[6],sample[1])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[4],sample[1],sample[7])
				} else {
					out<-c(sample[1],sample[4],sample[6],sample[7])
				}
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[1],sample[5],sample[5])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[1],sample[6],sample[5])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[1],sample[5],sample[7])
				} else {
					out<-c(sample[1],sample[1],sample[6],sample[7])
				}
			}
		} else if (sampleold[3]==0) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[4],sample[5],sample[5])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[4],sample[6],sample[5])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[4],sample[5],sample[7])
				} else {
					out<-c(sample[1],sample[4],sample[6],sample[7])
				}
		} else {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[3],sample[3],sample[3])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[3],sample[6],sample[3])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[3],sample[3],sample[7])
				} else {
					out<-c(sample[1],sample[3],sample[6],sample[7])
				}
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[4],sample[3],sample[3])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[4],sample[6],sample[3])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[4],sample[3],sample[7])
				} else {
					out<-c(sample[1],sample[4],sample[6],sample[7])
				}
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[3],sample[5],sample[5])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[1],sample[3],sample[6],sample[5])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[1],sample[3],sample[5],sample[7])
				} else {
					out<-c(sample[1],sample[3],sample[6],sample[7])
				}
			}
		}
	} else {
		if (sampleold[3]==ll[3]) {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[1],sample[1],sample[1])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[1],sample[6],sample[1])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[1],sample[1],sample[7])
				} else {
					out<-c(sample[2],sample[1],sample[6],sample[7])
				}
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[4],sample[1],sample[1])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[4],sample[6],sample[1])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[4],sample[1],sample[7])
				} else {
					out<-c(sample[2],sample[4],sample[6],sample[7])
				}
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[1],sample[5],sample[5])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[1],sample[6],sample[5])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[1],sample[5],sample[7])
				} else {
					out<-c(sample[2],sample[1],sample[6],sample[7])
				}
			}
		} else if (sampleold[3]==0) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[4],sample[5],sample[5])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[4],sample[6],sample[5])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[4],sample[5],sample[7])
				} else {
					out<-c(sample[2],sample[4],sample[6],sample[7])
				}
		} else {
			if (sampleold[4]==ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[3],sample[3],sample[3])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[3],sample[6],sample[3])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[3],sample[3],sample[7])
				} else {
					out<-c(sample[2],sample[3],sample[6],sample[7])
				}
			} else if (sampleold[4]<ll[4]&&sampleold[5]==ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[4],sample[3],sample[3])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[4],sample[6],sample[3])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[4],sample[3],sample[7])
				} else {
					out<-c(sample[2],sample[4],sample[6],sample[7])
				}
			} else if (sampleold[4]==ll[4]&&sampleold[5]<ll[5]) {
				if (sampleold[6]==ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[3],sample[5],sample[5])
				} else if (sampleold[6]<ll[6]&&sampleold[7]==ll[7]) {
					out<-c(sample[2],sample[3],sample[6],sample[5])
				} else if (sampleold[6]==ll[6]&&sampleold[7]<ll[7]) {
					out<-c(sample[2],sample[3],sample[5],sample[7])
				} else {
					out<-c(sample[2],sample[3],sample[6],sample[7])
				}
			}
		}
	}
	}
	if (length(ll)==9) {
	if (sample[8]<ll[8]) {
		if (length(Tt[[8]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[8]])) {
				if (Tt[[8]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[8]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[8]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[8]][i]+aa*dt
				}
				if (is.element(sample[8],idx)) {
					if (i!=1) {
					sample[8] <- floor(sum(Tt[[8]][1:(i-1)])/dt)+which(idx==sample[8])
					}
					break
				}
			}
		}
	}
	if (sample[9]<ll[9]) {
		if (length(Tt[[9]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[9]])) {
				if (Tt[[9]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[9]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[9]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[9]][i]+aa*dt
				}
				if (is.element(sample[9],idx)) {
					if (i!=1) {
					sample[9] <- floor(sum(Tt[[9]][1:(i-1)])/dt)+which(idx==sample[9])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[8]])<dt,dt-sum(Tt[[8]]),sum(Tt[[8]])-floor(sum(Tt[[8]])/dt)*dt)
	if (sample[7]<ll[7]&&sample[7]>0) {
		if (length(Tt[[7]])>1) {
			b <- 0
			for (i in 1:length(Tt[[7]])) {
				if (Tt[[7]][i]<t1) {
					idx <- b+1
					b <- b+1
					t1 <- t1-Tt[[7]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[7]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[7]][i]+aa*dt
				}
				if (is.element(sample[7],idx)) {
					if (i!=1) {
					sample[7] <- floor(sum(Tt[[7]][1:(i-1)])/dt)+which(idx==sample[7])
					}
					break
				}
			}
		}
		sample[7] <- sample[7]+floor(sum(Tt[[8]])/dt)
	} 
	if (sample[6]<ll[6]) {
		if (length(Tt[[6]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[6]])) {
				if (Tt[[6]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[6]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[6]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[6]][i]+aa*dt
				}
				if (is.element(sample[6],idx)) {
					if (i!=1) {
					sample[6] <- floor(sum(Tt[[6]][1:(i-1)])/dt)+which(idx==sample[6])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[6]])<dt,dt-sum(Tt[[6]]),sum(Tt[[6]])-floor(sum(Tt[[6]])/dt)*dt)
	if (sample[5]<ll[5]&&sample[5]>0) {
		if (length(Tt[[5]])>1) {
			b <- 0
			for (i in 1:length(Tt[[5]])) {
				if (Tt[[5]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[5]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[5]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[5]][i]+aa*dt
				}
				if (is.element(sample[5],idx)) {
					if (i!=1) {
					sample[5] <- floor(sum(Tt[[5]][1:(i-1)])/dt)+which(idx==sample[5])
					}
					break
				}
			}
		}
		sample[5] <- sample[5]+floor(sum(Tt[[6]])/dt)
	}
	if (sample[3]<ll[3]) {
		if (length(Tt[[3]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[3]])) {
				if (Tt[[3]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[3]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[3]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[3]][i]+aa*dt
				}
				if (is.element(sample[3],idx)) {
					if (i!=1) {
					sample[3] <- floor(sum(Tt[[3]][1:(i-1)])/dt)+which(idx==sample[3])
					}
					break
				}
			}
		}
	}
	if (sample[4]<ll[4]) {
		if (length(Tt[[4]])>1) {
			b <- 0
			t1 <- 0
			for (i in 1:length(Tt[[4]])) {
				if (Tt[[4]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[4]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[4]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[4]][i]+aa*dt
				}
				if (is.element(sample[4],idx)) {
					if (i!=1) {
					sample[4] <- floor(sum(Tt[[4]][1:(i-1)])/dt)+which(idx==sample[4])
					}
					break
				}
			}
		}
	}
	t1 <- ifelse(sum(Tt[[3]])<dt,dt-sum(Tt[[3]]),sum(Tt[[3]])-floor(sum(Tt[[3]])/dt)*dt)
	if (sample[2]<ll[2]&&sample[2]>0) {
		if (length(Tt[[2]])>1) {
			b <- 0
			for (i in 1:length(Tt[[2]])) {
				if (Tt[[2]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[2]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[2]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[2]][i]+aa*dt
				}
				if (is.element(sample[7],idx)) {
					if (i!=1) {
					sample[2] <- floor(sum(Tt[[2]][1:(i-1)])/dt)+which(idx==sample[2])
					}
					break
				}
			}
		}
		sample[2] <- sample[2]+floor(sum(Tt[[3]])/dt)
	}
	t1 <- ifelse(sum(Tt[[3]])<dt,dt-sum(Tt[[3]]),sum(Tt[[3]])-floor(sum(Tt[[3]])/dt)*dt)
	t1 <- ifelse((sum(Tt[[3]])-t1)<dt,dt-sum(Tt[[3]])+t1,sum(Tt[[3]])-t1-floor((sum(Tt[[3]])-t1)/dt)*dt)
	if (sample[2]==ll[2]||sample[5]==ll[5]) {
		if (length(Tt[[1]])>1) {
			b <- 0
			for (i in 1:length(Tt[[1]])) {
				if (Tt[[1]][i]<t1) {
					idx <- b+1
					b<-b+1
					t1 <- t1-Tt[[1]][i]
 				} else {
 					if (t1>0) {b <- b+1}
					aa <- floor((Tt[[1]][i]-t1)/dt)
					idx <- c(b+1:(aa+1))
					b <- b+aa+1
					t1 <- dt-Tt[[1]][i]+aa*dt
				}
				if (is.element(sample[1],idx)) {
					if (i!=1) {
					sample[1] <- floor(sum(Tt[[1]][1:(i-1)])/dt)+which(idx==sample[1])
					}
					break
				}
			}
		}
		sample[1] <- sample[1]+floor((sum(Tt[[2]])+sum(Tt[[3]]))/dt)
	}
	if (sampleold[2]==ll[2]) {
		if (sampleold[3]==ll[3]&&sampleold[4]==ll[4]) {
			out <- c(sample[1],sample[1])
		} else if (sampleold[3]<ll[3]&&sampleold[4]==ll[4]) {
			out <- c(sample[3],sample[1])
		} else if (sampleold[3]==ll[3]&&sampleold[4]<ll[4]) {
			out <- c(sample[1],sample[4])
		}
	} else if (sampleold[2]==0) {
		out <- c(sample[3],sample[4])
	} else {
		if (sampleold[3]==ll[3]&&sampleold[4]==ll[4]) {
			out <- c(sample[2],sample[2])
		} else if (sampleold[3]<ll[3]&&sampleold[4]==ll[4]) {
			out <- c(sample[3],sample[2])
		} else if (sampleold[3]==ll[3]&&sampleold[4]<ll[4]) {
			out <- c(sample[2],sample[4])
		}
	}
	if (sampleold[5]==ll[5]) {
	if (sampleold[6]==ll[6]) {
		if (sampleold[7]==ll[7]) {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[1],sample[1],sample[1])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[1],sample[8],sample[1])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[1],sample[1],sample[9])
			}
		} else if (sampleold[7]==0) {
			out<-c(out,sample[1],sample[8],sample[9])
		} else {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[1],sample[7],sample[7])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[1],sample[8],sample[7])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[1],sample[7],sample[9])
			}
		}
	} else {
		if (sampleold[7]==ll[7]) {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[1],sample[1])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[8],sample[1])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[6],sample[1],sample[9])
			}
		} else if (sampleold[7]==0) {
			out<-c(out,sample[6],sample[8],sample[9])
		} else {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[7],sample[7])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[8],sample[7])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[6],sample[7],sample[9])
			}
		}
	}
	} else if (sampleold[5]==0) {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[7],sample[7])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[8],sample[7])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[6],sample[7],sample[9])
			} else {
				out<-c(out,sample[6],sample[8],sample[9])
			} 
	} else if (sampleold[5]<ll[5]&&sampleold[5]>0) {
	if (sampleold[6]==ll[6]) {
		if (sampleold[7]==ll[7]) {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[5],sample[5],sample[5])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[5],sample[8],sample[5])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[5],sample[5],sample[9])
			}
		} else if (sampleold[7]==0) {
			out<-c(out,sample[5],sample[8],sample[9])
		} else {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[5],sample[7],sample[7])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[5],sample[8],sample[7])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[5],sample[7],sample[9])
			}
		}
	} else {
		if (sampleold[7]==ll[7]) {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[5],sample[5])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[8],sample[5])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[6],sample[5],sample[9])
			}
		} else if (sampleold[7]==0) {
			out<-c(out,sample[6],sample[8],sample[9])
		} else {
			if (sampleold[8]==ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[7],sample[7])
			} else if (sampleold[8]<ll[8]&&sampleold[9]==ll[9]) {
				out<-c(out,sample[6],sample[8],sample[7])
			} else if (sampleold[8]==ll[8]&&sampleold[9]<ll[9]) {
				out<-c(out,sample[6],sample[7],sample[9])
			}
		}
	}
	}
}
out
}