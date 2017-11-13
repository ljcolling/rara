# These functions do things that require parameter estimates
# But

# 0. Prior specification
priorRho <- function(rho, alpha=1) {
  priorDensity <- 2^(1-2*alpha)*(1-rho^2)^(alpha-1)
  logNormalisationConstant <- -lbeta(alpha, alpha)
  result <- exp(logNormalisationConstant)*priorDensity
  return(result)
}

priorRhoPlus <- function(rho, alpha=1) {
  nonNegativeIndex <- rho >=0
  lessThanOneIndex <- rho <=1
  valueIndex <- as.logical(nonNegativeIndex*lessThanOneIndex)
  myResult <- rho*0

  myResult[valueIndex] <- 2*priorRho(rho[valueIndex], alpha)
  return(myResult)
}

# 1.0. Built-up for likelihood functions
jeffreysApproxH <- function(n, r, rho) {
  return(((1 - rho^(2))^(0.5*(n - 1)))/((1 - rho*r)^(n - 1 - 0.5)))
}



# 2.1 Two-sided main Bayes factor
# See Ly et al 2014 or 2015

# 2.2 Two-sided secondairy Bayes factor
bf10JeffreysIntegrate <- function(n, r, alpha=1) {
  # Jeffreys' test for whether a correlation is zero or not
  # Jeffreys (1961), pp. 289-292
  # This is the exact result, see EJ
  ##
  if ( any(is.na(r)) ){
    return(NaN)
  }

  # TODO: use which
  if (n > 2 && abs(r)==1) {
    return(Inf)
  }

  hyperTerm <- Re(hypergeo::hypergeo((2*n-3)/4, (2*n-1)/4, (n+2*alpha)/2, r^2))
  logTerm <- lgamma((n+2*alpha-1)/2)-lgamma((n+2*alpha)/2)-lbeta(alpha, alpha)
  myResult <- sqrt(pi)*2^(1-2*alpha)*exp(logTerm)*hyperTerm
  return(myResult)
}


# 3.0 One-sided preparation

mPlusMarginalBJeffreys <- function(n, r, alpha=1){
  # Ly et al 2014
  # This is the exact result with symmetric beta prior on rho
  # This is the contribution of one-sided test
  #
  #
  if ( any(is.na(r)) ){
    return(NaN)
  }
  if (n > 2 && r>=1) {
    return(Inf)
  } else if (n > 2 && r<=-1){
    return(0)
  }

  hyperTerm <- Re(genhypergeo(U=c(1, (2*n-1)/4, (2*n+1)/4),
                              L=c(3/2, (n+1+2*alpha)/2), z=r^2))
  logTerm <- -lbeta(alpha, alpha)
  myResult <- 2^(1-2*alpha)*r*(2*n-3)/(n+2*alpha-1)*exp(logTerm)*hyperTerm
  return(myResult)
}


bfPlus0JeffreysIntegrate <- function(n, r, alpha=1){
  # Ly et al 2014
  # This is the exact result with symmetric beta prior on rho
  #
  if ( any(is.na(r)) ){
    return(NaN)
  }
  if (n > 2 && r>=1) {
    return(Inf)
  } else if (n > 2 && r<=-1){
    return(0)
  }

  bf10 <- bf10JeffreysIntegrate(n, r, alpha)
  mPlus <- mPlusMarginalBJeffreys(n, r, alpha)

  if (is.na(bf10) || is.na(mPlus)){
    return(NA)
  }

  myResult <- bf10+mPlus
  return(myResult)
}



###############################################################################
### Compute Bayes factors for Pearson's correlation coefficient
###
### Last revised: 06-01-2015
### Author: Josine Verhagen, Alexander Ly
### www.alexander-ly.com
# See Wagenmakers, E.-J., Verhagen, A. J. & Ly, A. (2015).
# 	How to quantify the evidence for the absence of a correlation
# 		Manuscript submitted for publication.
#
###############################################################################
#
# TODOTODO: check code for robustness.

estimationPosteriorU <- function(rho, n, r, alpha=1){
  dataTerm <- (1-rho^2)^((n-1)/2)/((1-rho*r)^((2*n-3)/2))*priorRho(rho, alpha)
  hyperTerm <- Re(hypergeo(1/2, 1/2, (2*n-1)/2, 1/2+1/2*r*rho))
  myResult <- dataTerm*hyperTerm
  return(myResult)
}

estimationPosteriorNormalisationConstant <- function(n, r, alpha=1){
  # The normalisation constant for the replication Bayes factor
  integrand <- function(x){estimationPosteriorU(x, n, r, alpha)}
  myResult <- integrate(integrand, -1, 1)$value
  return(myResult)
}

estimationPosterior <- function(rho, n, r, alpha=1){
  normalisationConstant <- estimationPosteriorNormalisationConstant(n, r, alpha)
  myResult <- 1/normalisationConstant*estimationPosteriorU(rho, n, r, alpha)
  return(myResult)
}

repPrior <- function(rho, nOri, rOri){
  estimationPosterior(rho, n=nOri, r=rOri)
}

# repPriorU <- function(rho, nOri, rOri){
# 	dataTerm <- (1-rho^2)^((nOri-1)/2)/((1-rho*rOri)^((2*nOri-3)/2))
# 	hyperTerm <- Re(hypergeo(1/2, 1/2, (2*nOri-1)/2, 1/2+1/2*rOri*rho))
# 	myResult <- dataTerm*hyperTerm
# 	return(myResult)
# }
#
# repPriorNormalisationConstant <- function(nOri, rOri){
# 	# The normalisation constant for the replication Bayes factor
# 	integrand <- function(x){repPriorU(x, nOri, rOri)}
# 	myResult <- integrate(integrand, -1, 1)$value
# 	return(myResult)
# }
#
# repPrior <- function(rho, nOri, rOri){
# 	normalisationConstant <- repPriorNormalisationConstant(nOri, rOri)
# 	myResult <- 1/normalisationConstant*repPriorU(rho, nOri, rOri)
# 	return(myResult)
# }

repPosteriorU <- function(rho, nOri, rOri, nRep, rRep){
  # Unnormalised posterior for the replication Bayes factor
  dataTerm <- (1-rho^2)^((nOri+nRep-2)/2)/((1-rho*rRep)^((2*nRep-3)/2)*(1-rho*rOri)^((2*nOri-3)/2))
  hyperTerm <- Re(hypergeo(1/2, 1/2, (2*nOri-1)/2, 1/2+1/2*rOri*rho))
  myResult <- dataTerm*hyperTerm
  return(myResult)
}

repPosteriorNormalisationConstant <- function(nOri, rOri, nRep, rRep){
  # The normalisation constant for the replication Bayes factor
  integrand <- function(x){repPosteriorU(x, nOri, rOri, nRep, rRep)}
  myResult <- integrate(integrand, -1, 1)$value
  return(myResult)
}

repPosterior <- function(rho, nOri, rOri, nRep, rRep){
  normalisationConstant <- repPosteriorNormalisationConstant(nOri, rOri, nRep, rRep)
  myResult <- 1/normalisationConstant*repPosteriorU(rho, nOri, rOri, nRep, rRep)
  return(myResult)
}

repBfR0 <- function(rho=0, nOri, rOri, nRep, rRep){
  myResult <- repPrior(rho, nOri, rOri)/repPosterior(rho, nOri, rOri, nRep, rRep)
  return(myResult)
}

# Just some plotting
plotRhoPrior <- function(gammas){
  numberOfGammas <- length(gammas)

  myDomain <- seq(-1, 1, by=0.01)
  collectionOfPriors <- array(dim=c(numberOfGammas, length(myDomain)))

  for (i in 1:numberOfGammas){
    collectionOfPriors[i, ] <- priorRho(myDomain, alpha=1/gammas[i])
  }

  yLim <- max(collectionOfPriors)

  par(cex.main=1.5, mar=c(5, 6, 4, 5) + 0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
      font.lab=2, cex.axis=1.3, bty="n", las=1)
  plot(x=NULL, y=NULL, xlim=c(-1, 1), ylim=c(0, yLim),
       ylab="Prior density", xlab=expression(paste("Correlation ", rho)))

  for (i in 1:numberOfGammas){
    lines(myDomain, collectionOfPriors[i, ], col=i, lwd=2.5, lty=2)
  }

  legendText <- vector(, numberOfGammas)

  for (i in 1:numberOfGammas){
    legendText[i] <- paste(expression(gamma), "=", round(gammas[i], 2))
  }


  legend('topright', legend=legendText,
         lty=2, col=1:numberOfGammas, bty='n', cex=1.2, lwd=2)
}


plotSensitivity <- function(BF01s, myGamma, xLegend=0.52, arrowCorrect=0.3,
                            arrowStartHeight=1.25, arrowEndHeight=0.15,
                            textHeight=2.25){
  par(cex.main=1.5, mar=c(3, 5, 4, 2) + 0.1, mgp=c(2, 0.25, 0),
      cex.lab=1.5, font.lab=2, cex.axis=1.3, bty="n", las=1)

  plot(myGamma, log(BF01s), xlab=expression(gamma), ylab="", ylim=c(-log(3), log(100)),
       xlim=c(0, 1), main="Study 1", axes=FALSE, lwd=2.5, type="l")


  axis(1, at=seq(0, 1, by=0.25))
  axis(2, tcl=0, labels=c("-log(10)", "-log(3)", "log(1)", "log(3)", "log(10)",
                          "log(30)", "log(100)"), at=c(-log(10), -log(3),
                                                       log(1), log(3), log(10),
                                                       log(30), log(100)))

  abline(h=c(-log(10),-log(3), log(1),log(3),log(10),log(30),log(100)),
         col='grey', lwd=2, lty=2)
  abline(h=c(log(1)), col='darkgrey',lwd=2,lty=2)

  spreadOneIndex <- which(myGamma==1)
  points(1, log(BF01s[spreadOneIndex]), lwd=3)


  BFNumber <- round(BF01s[spreadOneIndex], 2)

  # Change BF01s[spreadOneIndex] to analyticBF01[studyNumber]
  arrows(x0=myGamma[spreadOneIndex]-arrowCorrect,
         y0=log(BF01s[spreadOneIndex]*exp(arrowStartHeight)),
         x1=myGamma[spreadOneIndex],
         y1=log(BF01s[spreadOneIndex]*exp(arrowEndHeight)),
         col="black", lwd=2, length=0.12)

  legend(xLegend, log(BF01s[spreadOneIndex]*exp(textHeight)),
         substitute(BF[0][1]==BFNumber, list(BFNumber=BFNumber)),
         bty="n", cex=1.7)
}

plotRhoPosterior <- function(n, r, alpha=1){
  numberOfGammas <- length(gammas)

  myDomain <- seq(-1, 1, by=0.01)
  collectionOfPriors <- array(dim=c(numberOfGammas, length(myDomain)))

  for (i in 1:numberOfGammas){
    collectionOfPriors[i, ] <- priorRho(myDomain, alpha=1/gammas[i])
  }

  yLim <- max(collectionOfPriors)
  yLim <- 28
  par(cex.main=1.5, mar=c(5, 6, 4, 5) + 0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
      font.lab=2, cex.axis=1.3, bty="n", las=1)
  plot(x=NULL, y=NULL, xlim=c(-1, 1), ylim=c(0, yLim),
       ylab="Prior density", xlab=expression(paste("Correlation ", rho)))

  for (i in 1:numberOfGammas){
    lines(myDomain, collectionOfPriors[i, ], col=i, lwd=2.5)
  }

  legendText <- vector(, numberOfGammas)

  for (i in 1:numberOfGammas){
    legendText[i] <- paste(expression(gamma), "=", round(gammas[i], 2))
  }


  legend('topright', legend=legendText,
         lty=1, col=1:numberOfGammas, bty='n', cex=1.2, lwd=2)
}

#  Useless
makeGammas <- function(n){
  myGamma <- sin(seq(1.5*pi, 2*pi, length=n))+1
  myGamma[1] <- myGamma[2]/10
  myGamma[n] <- 1
  return(myGamma)
}

plotRhoPriorPosterior <- function(n, r, alpha, yLim=10){
  rhoDomain <- seq(-1, 1, by=0.001)

  myGamma <- round(1/alpha, 2)
  myTitle <- substitute(paste("Prior and posterior with scale ", gamma, "=", v), list(v=myGamma))

  par(cex.main=1.5, mar=c(5, 6, 4, 5) + 0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
      font.lab=2, cex.axis=1.3, bty="n", las=1)

  plot(x=NULL, y=NULL, xlim=c(-1, 1), ylim=c(0, yLim), axes=FALSE,
       xlab=expression(paste("Correlation ", rho)), ylab="Density",
       main=myTitle)

  # Prior
  somePrior <- priorRho(rhoDomain, alpha)
  priorNull <- priorRho(0, alpha)
  lines(rhoDomain, y=somePrior, lty=2, lwd=2.5)

  # Posterior posterior(myDomain
  somePosterior <- estimationPosterior(rhoDomain, n=n, r=r, alpha=alpha)
  posteriorNull <- estimationPosterior(0, n=n, r=r, alpha=alpha)

  lines(rhoDomain, somePosterior, lwd=2.5)

  # Draw Savage Dickey points
  points(0, priorNull, pch=19, cex=2.5, col="grey")
  points(0, posteriorNull, pch=19, cex=2.5, col="grey")
  points(0, priorNull, pch=21, cex=2.5, col="black")
  points(0, posteriorNull, pch=21, cex=2.5, col="black")
  axis(1)
  axis(2)


  legendText <- c("Prior", "Posterior")

  legend("topleft", legend=legendText,
         lty=c(2, 1), bty="n", cex=1.5)
}

countNonNAEntries <- function(x, y){
  xCheck <- as.numeric(!is.na(x))
  yCheck <- as.numeric(!is.na(y))
  sum(xCheck*yCheck)
}


giveBF<-function(input){
  returnvalue = ""
  if(input <1){
    returnvalue = paste0(" = ",sprintf("%.3f",input))
  }
  if(input>1 && input < 100 ){
    returnvalue = paste0(" = ", round(input,2))
  }

  if(input>100){
    returnvalue = paste0(" = ", round(input))
  }

  if(input>1000){
    returnvalue = "> 1000"
  }
  if(input < 1/1000){
    returnvalue = "< 0.001"
  }
  return(returnvalue)
}


giveBF.noequal<-function(input){
  returnvalue = ""
  if(input <1){
    returnvalue = paste0(sprintf("%.3f",input))
  }
  if(input>1 && input < 100 ){
    returnvalue = paste0(round(input,2))
  }

  if(input>100){
    returnvalue = paste0(round(input))
  }

  if(input>1000){
    returnvalue = "> 1000"
  }
  if(input < 1/1000){
    returnvalue = "< 0.001"
  }
  return(returnvalue)
}




DoCorr<-function(x.name,y.name,data){

  x = data[,x.name]
  y = data[,y.name]

  df<-data.frame(x = x, y = y)
  df<-df[complete.cases(df),]
  x = df$x
  y = df$y

  # check if the bayesian correlation can already been run
  if(run.bayes == TRUE){



    if(file.exists(this.file)  == TRUE){
      robj = read.BayesCorrObj(this.file)
    }
    else{
      robj = BayesCorr(x,y, x.name, y.name)
      save("robj",file = this.file)

    the.file = paste0("corr.",x.name,"_",y.name,".roj")
    if(file.exists(the.file)  == TRUE){
      robj = read.BayesCorrObj(x.name,y.name)
    }
  }


  if(run.bayes == FALSE){



    cor.obj = cor.test(x,y)

    robj = list()
    robj$fit = NULL
    robj$dist = NULL
    robj$cor$r.estimate = cor.obj
    robj$cor$rho = unname(cor.obj$estimate)
    robj$cor$n = length(x)
    robj$cor$bfAlt = bf10JeffreysIntegrate(n = robj$cor$n, r = robj$cor$rho )
    robj$cor$bfAltText = giveBF(robj$cor$bfAlt)
    robj$cor$bfNull = 1/bf10JeffreysIntegrate(n = robj$cor$n, r = robj$cor$rho )
    robj$cor$bfNullText = giveBF(robj$cor$bfNull)
    robj$cor$bfEvidence = NULL
    robj$cor$bfHypothesis = NULL
    robj$cor$x.mean = mean(x)
    robj$cor$y.mean = mean(y)
    robj$cor$x.name = x.name
    robj$cor$y.name = y.name
    robj$cor$x.data = x
    robj$cor$y.data = y
    robj$cor$hdi = cor.obj$conf.int[c(1,2)]
    robj$cor$sig = cor.obj$p.value
    robj$mcmc = NULL

  }

  return(robj)
}

ParseCor<-function(corobj){

  return(paste0(paste0("*r* = ",weights::rd(corobj$cor$rho)),", ",
         paste0("95% HDI[",paste0(weights::rd(corobj$cor$hdi),collapse = "; "),"]"),", ",
         paste0("*BF*~10~ ",corobj$cor$bfAltText),", ",
         paste0("*BF*~01~ ",corobj$cor$bfNullText)))
}


# DoMeanEs.OneSample<-function(data, this.file = NULL){
#   output = list()
#   data = na.omit(data)
#
#   if(run.bayes == T){
#     if(file.exists(this.file) == T){
#       tempenv = new.env()
#       load(file = this.file, envir = tempenv)
#       best.temp = tempenv$best.temp
#     } else{
#       best.temp = BEST::BESTmcmc(data)
#       save(best.temp, file = this.file)
#     }
#
#
#     output$n = length(data)
#     output$mean = mean(best.temp$mu)
#     output$hdi = unname(HDInterval::hdi(best.temp$mu))
#     output$efsz = posteriorSummary(best.temp$mu / best.temp$sigma)$mode
#     output$efsz.ci = unname(HDInterval::hdi(best.temp$mu / best.temp$sigma))
#   }
#
#   if(run.bayes == F){
#
#
#     output$n = length(data)
#     output$mean = mean(data)
#     output$hdi = unname(t.test(data)$conf.int)[c(1,2)]
#     efsz = bootES::bootES(data,effect.type = "cohens.d")
#     output$efsz = unname(efsz$t0)
#     output$efsz.ci = unname(efsz$bounds)
#
#   }
#
#   return(output)
# }

Parse.MeanEs<-function(bestobj, sub, units){
# sub = "~diff~"
# units = "ms"
return(paste0(
  paste0("*M*",sub," = ",round(bestobj$mean,2)),", ",
  paste0("95% HDI[",paste0(round(bestobj$hdi,2),collapse = "; "),"]"),"; ",
  paste0("*d* = ",round(bestobj$efsz,2)),", ",
  paste0("95% HDI[",paste0(round(bestobj$efsz.ci,2),collapse = "; "),"]")))
}



