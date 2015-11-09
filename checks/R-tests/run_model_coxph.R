##' This function runs the actual models for Cox regression. It is
##' called by run_models_in_R_pacox.R
##'
##'
##' @title run.model
##' @param model0.txt String containing the null model (without SNP term)
##' @param model.txt String containing the alternative model (with SNP
##' term)
##' @param snpcomponent1 String telling how the SNP term is defined
##' @param snpcomponent2 String telling how the second SNP term is
##' defined (only used in the 2 df model). By default this term is
##' constant ("1")
##' @return A data frame containing the coefficients from the
##' regression analysis and some other variables, such that this
##' output can be compared to the ProbABEL output.
##' @author L.C. Larsen
run.model <- function(model0.txt, model.txt,
                      snpcomponent1, snpcomponent2="1",
                      verbose=FALSE) {
    if (snpcomponent2 != "1") {
        ## SNP component 2 is not constant: assume we run the 2df
        ## model.
        twoDF = TRUE
    } else {
        twoDF = FALSE
    }

    resultR <- data.frame()

    for (i in 3:dim(dose)[2]) {
        if (verbose) print(paste("------- new iteration: i =", i))
        indexHom <- 3 + ( i - 3 ) * 2
        indexHet <- indexHom + 1
        snp1     <- eval(parse(text=snpcomponent1))
        snp2     <- eval(parse(text=snpcomponent2))
        snp      <- snp1 + snp2
        mu       <- rep(1., length(snp)) # Add a constant mean

        noNA    <- which( !is.na(snp) )
        model.0 <- eval(parse(text=model0.txt))

        ## Check the imputation R^2, if below threshold ProbABEL will
        ## set the coefficients to NaN and skip the rest of this
        ## iteration of the for loop.
        rsq <- Rsq[i-2]
        if (rsq < rsq.thresh ) {
            if (twoDF) {
                row <- c(rsq, NaN, NaN, NaN, NaN, NaN)
            } else {
                row <- c(rsq, NaN, NaN, NaN)
            }
            resultR <- rbind(resultR, row)
            next
        }

        ## Evaluate the model. The whole tryCatch is needed to catch
        ## problems with non-converging regression.
        model = tryCatch({
            list(
                eval(parse(text=model.txt)),
                list(message="no warnings")
                )
        }, warning = function(war) {
            return(list(
                eval(parse(text=model.txt)),
                war)
                   )
        },
            error = function(err) {
                return(list(
                    "Can't calculate model; some error occurred",
                    err)
                       )
            })

        if (verbose) print(model)

        if ( grepl("Inf", model[[2]]$message) |
            grepl("infinite", model[[2]]$message) |
            grepl("iterations", model[[2]]$message) |
            ( grepl("singular", model[[2]]$message) &
                 (regexpr("variable 1$", model[[2]]$message)[[1]] != 33)
             )
            ) {
            ## The model did not converge or some other
            ## errors/warnings occurred, fill the coefficients with
            ## NaNs
            if (twoDF) {
                smA1A2  <- c(NaN, NaN)
                smA1A1  <- c(NaN, NaN)
            } else {
                sm      <- c(NaN, NaN)
            }
            lrt <- NaN
        } else {
            ## No convergence problems, we can trust the
            ## coefficients.
            coeff <- summary(model[[1]])$coefficients
            # The SNP coefficient(s) is/are in the last row(s) of coeff
            nr <- nrow(coeff)
            if (twoDF) {
                smA1A2 <- coeff[nr - 1, c("coef", "se(coef)")]
                smA1A1 <- coeff[nr, c("coef", "se(coef)")]
            } else {
                sm     <- coeff[nr, c("coef", "se(coef)")]
            }
            lrt   <- 2 * ( model[[1]]$loglik[2] - model.0$loglik[2] )
        }

        if (twoDF) {
            row <- c(rsq, smA1A2[1], smA1A2[2], smA1A1[1], smA1A1[2], lrt)
        } else {
            row <- c(rsq, sm[1], sm[2], lrt)
        }

        resultR <- rbind(resultR, row)
    }
    return(resultR)
}
