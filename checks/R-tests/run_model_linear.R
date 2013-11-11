run.model <- function(model0.txt, model.txt, snpdata) {
    resultR <- data.frame()
    for (i in 3:dim(dose)[2]) {
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        snp      <- eval(parse(text=snpdata))

        noNA    <- which( !is.na(snp) )
        model.0 <- eval(parse(text=model0.txt))
        model   <- eval(parse(text=model.txt))

        coeff   <- summary(model)$coefficients
        if ( dim(coeff)[1] != 4 ) {
            sm <- c(NaN, NaN)
        } else {
            sm <- coeff[4, c("Estimate", "Std. Error")]
        }

        lrt <- 2 * ( logLik( model ) - logLik( model.0 ) )
        rsq <- Rsq[i-2]
        if( rsq < rsq.thresh) {
            row <- c(rsq, NaN, NaN, NaN)
        } else {
            row <- c(rsq, sm[1], sm[2], lrt)
        }
        resultR <- rbind(resultR, row)
    }
    return(resultR)
}
