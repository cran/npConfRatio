npconfratio <- function(x,y, conf.level=0.95, alternative=c("two.sided","less", "greater"), equiv.test = TRUE, equiv.range=c(0.8, 1.25),
pest = TRUE)
{
	alternative <- match.arg(alternative)
	if (is.null(x) | length(x) < 1)
		stop("not enough x observations")
	if (is.null(y) | length(y) < 1)
		stop("not enough y observations")
	DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
		conf.level < 0 || conf.level > 1))
			stop("conf.level must be a single number between 0 and 1")
	alpha <- 1 - conf.level
	PVALUE <- !is.null(equiv.range)
	if (length(equiv.range) < 2 & alternative =="two.sided")
		stop("no two.sided equivalence range defined")

	m <- length(x)
	n <- length(y)
	N <- m + n
	EXACT <- ((m < 50) && (n < 50))
	r <- rank(c(x,y))
	TIES <- (length(r) != length(unique(r)))
	if (!TIES | !EXACT) {
		if (length(y[y < 0]) > 0 || length(x[x <= 0]) > 0)
			stop("x <= 0 and y < 0 not allowed");
 		ratio <- sort(outer(y, x, "/"))
		if (!EXACT)
			warning("Cannot compute exact confidence intervals and p-values for sample sizes > 50")	
		CINT <- switch(alternative, two.sided={
			if (EXACT)
				qu <- qwilcox(alpha/2, m, n) 
			else
				qu <- round(qnorm(alpha/2)*sqrt(m*n*(N+1)/12) + n*m/2)	
			if (qu == 0) qu <- 1
			ql <- m*n - qu
			if (length(ratio[ratio == ratio[qu]]) > 1)
				qu <- which(ratio == min(ratio[ratio > ratio[qu]]))
			if (length(ratio[ratio == ratio[ql+1]]) > 1)
				ql <- which(ratio == max(ratio[ratio < ratio[ql + 1]]))
			uci <- ratio[qu]
			lci <- ratio[ql + 1]
			c(uci, lci)
		}, greater ={
			if (EXACT)	
				qu <- qwilcox(alpha, m, n)
			else
				qu <- round(qnorm(alpha)*sqrt(m*n*(N+1)/12) + n*m/2)
			if (length(ratio[ratio == ratio[qu]]) > 1)
				qu <- which(ratio == min(ratio[ratio > ratio[qu]]))
			if (qu == 0) uci <- ratio[1]
			else uci <- ratio[qu]
			c(uci, NA)
		}, less={
			if (EXACT)
				qu <- qwilcox(alpha, m, n)
			else
				qu <- round(qnorm(alpha)*sqrt(m*n*(N+1)/12) + n*m/2)
			if (qu == 0) qu <- 1
			ql <- m*n - qu
			if (length(ratio[ratio == ratio[ql+1]]) > 1)
				ql <- which(ratio == max(ratio[ratio < ratio[ql + 1]]))
			lci <- ratio[ql+1]
			c(NA, lci)
		})

		attr(CINT, "conf.level") <- conf.level

		# Point-Estimator

		if (pest) {
			EST <- ratio[qwilcox(0.5, m, n)]
			names(EST) <- "Point-Estimator"
		}

		# Equivalence test

		if (PVALUE) {
			PVAL <- switch(alternative, two.sided={
				r1 <- rank(c(y/equiv.range[1], x))
				if (EXACT)	
					p1 <- 1 - pwilcox(sum(r1[seq(along=y)]) - n*(n+1)/2, m, n)
				else 
					p1 <- 1 - pnorm((sum(r1[seq(along=y)]) - n*(n+1)/2 - n*m/2)/sqrt(m*n*(N+1)/12))
				r2 <- rank(c(y/equiv.range[2], x))
				if (EXACT)
					p2 <- pwilcox(sum(r2[seq(along=y)]) - n*(n+1)/2, m, n)
				else
					p2 <- pnorm((sum(r2[seq(along=y)]) - n*(n+1)/2 - n*m/2)/sqrt(m*n*(N+1)/12))
				max(c(p1, p2))
			}, greater={
				r1 <- rank(c(y/equiv.range[1], x))
				if (EXACT)
					1 - pwilcox(sum(r1[seq(along=y)]) - n*(n+1)/2, m, n)
				else
					1 - pnorm((sum(r1[seq(along=y)]) - n*(n+1)/2 - n*m/2)/sqrt(m*n*(N+1)/12))
			}, less={
				r2 <- rank(c(y/equiv.range[2], x))
				if (EXACT)
					pwilcox(sum(r2[seq(along=y)]) - n*(n+1)/2, m, n)
				else
					pnorm((sum(r2[seq(along=y)]) - n*(n+1)/2 - n*m/2)/sqrt(m*n*(N+1)/12))
			})
		}						
	} else {
		if (length(y[y < 0]) > 0 || length(x[x <= 0]) > 0)
			stop("x <= 0 and y < 0 not allowed")
		ratio <- sort(outer(y, x, "/"))
		w <- cumsum(c(n*m-1, rep(-1, m*n-1)))
		dup <- which(duplicated(ratio))
		indx <- (1:length(ratio))
		indx[dup] <- NA
		w <- w[!is.na(indx)]
		ratio <- ratio[!is.na(indx)]
		scores <- rank(c(y,x))
		if (!EXACT)
			warning("Cannot compute exact confidence intervals and p-values for sample sizes > 50")	
		CINT <- switch(alternative, two.sided={
			qu <- round(qperm(alpha/2, scores, n) - n*(n+1)/2)
			ql <- round(qperm(1-alpha/2, scores, n) - n*(n+1)/2)
			if (ql >= max(w)) uci <- 0
			else uci <- max(ratio[w > ql])
			if (qu <= min(w)) lci <- max(ratio)
			else lci <- min(ratio[w < qu])
			c(uci, lci)
		}, greater ={
			ql <- round(qperm(1-alpha, scores, n) - n*(n+1)/2)
			uci <- max(ratio[w > ql])
			c(uci, NA)
		}, less={
			qu <- round(qperm(alpha, scores, n) - n*(n+1)/2)
			lci <- min(ratio[w < qu])
			c(NA, lci)
		})

		attr(CINT, "conf.level") <- conf.level

		if (pest) {

			EST <- ratio[round(qperm(0.5, scores, n) - n*(n+1)/2)]
			names(EST) <- "Point-Estimator"
		}

		# Equivalence test

		if (PVALUE) {
			PVAL <- switch(alternative, two.sided={
				r1 <- rank(c(y/equiv.range[1], x))
				p1 <- 1 - pperm(sum(r1[seq(along=y)]), r1, n)
				r2 <- rank(c(y/equiv.range[2], x))
				p2 <- pperm(sum(r2[seq(along=y)]), r2, n)
				max(c(p1, p2))
			}, greater={
				r1 <- rank(c(y/equiv.range[1],x))
				1 - pperm(sum(r1[seq(along=y)]), r1, n)
			}, less={
				r2 <- rank(c(y/equiv.range[2],x))
				pperm(sum(r2[seq(along=y)]), r2, n)
			})
		}	
	}

	METHOD <- "Nonparametric Confidence Intervals for the Ratio of Medians"

	if (PVALUE & pest)
		RVAL <- list(p.value = PVAL, alternative = alternative, method =
		METHOD, conf.int = CINT, conf.level = conf.level, estimate = EST, data.name = DNAME)
	else {
		if (PVALUE)
			RVAL <- list(p.value = PVAL, alternative = alternative, method =
			METHOD, conf.int = CINT, conf.level = conf.level, data.name = DNAME)
		else
			RVAL <- list(alternative = alternative, method =
			METHOD, conf.int = CINT, conf.level = conf.level, estimate = EST, data.name = DNAME)
	}						
	class(RVAL) <- "htest"
	return(RVAL)
}
