in.pheno <- read.table("coxph_data.txt", header=TRUE, as.is=TRUE)
in.pheno[1:5,]
attach(in.pheno)

dose <- read.table("./test.mldose", header=FALSE, as.is=TRUE)
dose[1:5,]
table(in.pheno[,"id"]== sapply(dose[,1],
		function(x){strsplit(x,split="->")[[1]][2]})
      )

# Look at first SNP

i <- 3

library(survival)

test.coxph <- coxph(Surv(fupt_chd,chd) ~ sex + age + othercov + dose[,i],
		    x=TRUE, y=TRUE )

n <- length(test.coxph$residuals)
#  rr <- object$residuals
y <- test.coxph$y
x <- test.coxph[['x']]  # avoid matching object$xlevels
vv <- test.coxph$var

weights <- rep(1,n)

ny <- ncol(y)
status <- y[,ny,drop=TRUE]
nvar <- ncol(x)

ord <- order(y[,ny-1], -status)
newstrat <- rep(0,n)
newstrat[n] <- 1

# sort the data
x <- x[ord,]
y <- y[ord,]
score <- exp(test.coxph$linear.predictors)[ord]

method<-c("erfron")

resid <- .C("coxscore", as.integer(n),
	    as.integer(nvar),
	    as.double(y),
	    x=as.double(x),
	    as.integer(newstrat),
	    as.double(score),
	    as.double(weights[ord]),
	    as.integer(method=='efron'),
	    resid=double(n*nvar),
	    double(2*nvar))$resid

tmp.rr <-  matrix(0, n, nvar)

tmp.rr <- matrix(resid, ncol=nvar)

rr <- matrix(0, n, nvar)
rr[ord,] <- matrix(resid, ncol=nvar)
dimnames(rr) <- list(names(test.coxph$residuals),
		     names(test.coxph$coefficients))
rr <- rr %*% vv
rr <- rr * weights   ## At this point, rr is what is returned by a call to residuals(test.coxph,type="dfbeta")


## compare to call to get robust variances
test.coxph.robust <- coxph(Surv(fupt_chd,chd) ~ sex + age + othercov +
			   dose[,i],
			   x=TRUE, y=TRUE, robust=T )
robust.variances <- t(rr) %*% rr
test.coxph.robust$var
