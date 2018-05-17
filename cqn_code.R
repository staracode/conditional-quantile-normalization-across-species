library("quantreg")
library("splines")
library("mclust")
library("nor1mix")

source("cqn_original_code/R/SQN.R")
cqn <- function(counts, x, lengths, sizeFactors = NULL, subindex = NULL,
                tau = 0.5, sqn = TRUE, lengthMethod = c("smooth", "fixed"),
                verbose = FALSE) {

#cqn(subset2[isexpr,], lengths = len_matrix[isexpr,], x = gc_matrix[isexpr,], sizeFactors = library_size, verbose = TRUE)

test = FALSE
if (test == TRUE){
	#counts = subset2[isexpr,]
	counts = cbind(sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL))
	#lengths = len_matrix[isexpr,]
	lengths = cbind(sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL))
	#x = gc_matrix[isexpr,]
	#x = cbind(seq(1,100), seq(1,100), seq(1,100), seq(1,100))
	x = cbind(sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL), 
	sample(seq(1,1000), 100, replace = TRUE, prob = NULL))
	sizeFactors = NULL
	verbose = TRUE
	lengthMethod = "smooth"
	subindex = NULL
	tau = 0.5
	verbose = TRUE
}

	#lengthMethod <- match.arg(lengthMethod)
	cl <- match.call()

	if(is.null(sizeFactors))
		sizeFactors <- colSums(counts)
  
	if(is.null(subindex))
		subindex <- which(rowMeans(counts) > 50)
  
	y <- sweep(log2(as.matrix(counts) + 1), 2, log2(sizeFactors/10^6))
	if(lengthMethod == "fixed")
		y <- sweep(y, 1, log2(lengths/1000))
	yfit <- y[subindex,, drop = FALSE]
	head(yfit)
  
	###############################################
	fixPredictor <- function(zz, varname = "B") {
		zz.fit <- zz[subindex]
		knots <- quantile(zz.fit, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) + c(+0.01, 0,0,0, -0.01)
		grid <- seq(from = min(zz.fit), to = max(zz.fit), length.out = 101)
		zz.out <- zz
		zz2 <- zz[-subindex]
		zz2[zz2 < min(zz.fit)] <- min(zz.fit)
		zz2[zz2 > max(zz.fit)] <- max(zz.fit)
		zz.out[-subindex] <- zz2
		list(fit = zz.fit, out = zz.out, knots = knots, grid = grid)
	}
	###############################################
    
	regr2 <- lapply(1:ncol(yfit), FUN = function(i) {
		x1p <- fixPredictor(x[,i])
		x1.knots <- x1p$knots
  
		x2p <- fixPredictor(log2(lengths[,i]/1000))
		x2.knots <- x2p$knots

		df.fit <- data.frame(x1 = x1p$fit, x2 = x2p$fit)
		df.predict <- data.frame(x1 = x1p$out, x2 = x2p$out)

		df.fit$y <- yfit[,i]
		fit <- rq(y ~ ns(x1, knots = x1.knots) + ns(x2, knots = x2.knots), data = df.fit, tau = tau)
		
		fitted <- predict(fit, newdata = df.predict)
		list(fitted = fitted,  coef = coef(fit))
		
	} )
	
	fitted <- do.call(cbind, lapply(regr2, function(xx) xx$fitted))
  
	#   ## This is the across sample median for the fitted values for the
	#   ## median value of the predictor amongst the subindex
	k <- lapply (seq(1, ncol(x)), function (i) order(x[,i][subindex])[length(subindex) / 2])
	k2 = unlist(k)
	offset0 <- lapply (seq(1, length(k2)), function(ii) median(fitted[subindex[k2[ii]],]) )
	offset1 = unlist(offset0)
  
	residuals <- y - fitted
	
	if(verbose) cat("SQN ")
	sqnFit <- SQN2(residuals, ctrl.id = subindex, min.q = 0.0001)
	residualsSQN <- sqnFit$yout
	if(verbose) cat(".\n")
	offset <- residualsSQN + offset1 - y
  
	rownames(offset) <- rownames(y)
	colnames(offset) <- colnames(y)
   
	glm.offset <- sweep(-offset, 2, -log2(sizeFactors / 10^6))
	glm.offset <- log(2)*glm.offset
	out <- list(counts = counts, lengths = lengths, sizeFactors = sizeFactors,
    	subindex = subindex, 
    	y = y, x = x, offset = offset, offset0 = offset1,
    	glm.offset = glm.offset)
	class(out) <- "cqn"
	out
}

