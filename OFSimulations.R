rm(list=ls())
gsList <- 0.3
gxList <- 0.9
mList <- c(0.2, 0.5, 1)
sdevList <- c(0.01, 0.1, 1)
replicatesList <- 3
stepSizeList <- 0.1
library(splines)
ovules <- 100
alpha.scale <- 0.1
responseScale <- function(q) as.real(q)
parameterEstimates <- function(the.model) {
	responseScale(c(the.model$coefficients[1], the.model$coefficients[2] + the.model$coefficients[1]))
}
addEst <- function(gs, gx, m, alpha, i) {
	outputTable$gsEst[i] <- gs
	outputTable$gxEst[i] <- gx
	outputTable$mEst[i] <- m
	outputTable$alphaEst[i] <- alpha
	outputTable
}
tableDimension <- length(gsList) * length(gxList) * length(mList) * length(sdevList) * length(replicatesList) * length(stepSizeList)
outputTable <- NULL
outputTable$gs <- rep(gsList, each=tableDimension / (length(gsList)))
outputTable <- as.data.frame(outputTable)
outputTable$gx <- rep(gxList, each=tableDimension / ((length(gsList) * length(gxList))))
outputTable$m <- rep(mList, each=tableDimension / ((length(gsList) * length(gxList) * length(mList))))
conditions <- (outputTable$gx > outputTable$gs & outputTable$gs > outputTable$m)
outputTable$type[conditions] <- c("a")
conditions <- (outputTable$gx > outputTable$m & outputTable$m > outputTable$gs)
outputTable$type[conditions] <- c("b")
conditions <- (outputTable$m > outputTable$gx & outputTable$gx > outputTable$gs)
outputTable$type[conditions] <- c("c")
outputTable$sdev <- rep(sdevList, each=tableDimension / ((length(gsList) * length(gxList) * length(mList) * length(sdevList))))
outputTable$replicates <- rep(replicatesList, each=tableDimension / ((length(gsList) * length(gxList) * length(mList) * length(sdevList) * length(replicatesList))))
outputTable$stepSize <- rep(stepSizeList, each=tableDimension / ((length(gsList) * length(gxList) * length(mList) * length(sdevList) * length(replicatesList) * length(stepSizeList))))
conditions <- (outputTable$gx > outputTable$m & outputTable$m > outputTable$gs)
outputTable$alpha <- ifelse(conditions, (outputTable$m-outputTable$gs)/(outputTable$gx-outputTable$gs), NA)
m1aic <- rep(NA, dim(outputTable)[1]); m2aic <- m1aic; gsEst <- m1aic; gxEst <- m1aic; mEst <- m1aic; alphaEst <- m1aic; gsError <- m1aic; gxError <- m1aic; mError <- m1aic; alphaError <- m1aic;
outputTable <- cbind(outputTable, m1aic, m2aic, gsEst, gxEst, mEst, alphaEst)
family <- quasi(var="mu(1-mu)", link="identity")
for (i in 1:dim(outputTable)[1]) {
	crossPollen <- seq(0, 1, outputTable$stepSize[i])
	seeds <- rep((crossPollen * ovules * outputTable$gx[i] + (1 - crossPollen) * ovules * outputTable$gs[i]), outputTable$replicates[i]) + ovules * rnorm(length(crossPollen) * outputTable$replicates[i], mean=0, sd=outputTable$sdev[i])
	conditions <- (seeds > outputTable$m[i] * ovules)
	seeds[conditions] <- outputTable$m[i] * ovules + ovules * rnorm(length(seeds[conditions]), mean=0, sd=outputTable$sdev[i])
	seeds[seeds > ovules] <- ovules
	seeds[seeds < 0] <- 0
	data <- cbind(crossPollen, round(seeds, 0))
	response <- cbind(data[,2], ovules)
	model.1 <- glm(response ~ data[,1],family=family)
	outputTable$m1aic[i] <- model.1$aic
	if (summary(model.1)$coefficients[2,4] > 0.05) {
		outputTable <- addEst(NA, NA, responseScale(summary(model.1)$coefficients[1,1]), NA, i)
	} else {
		deviances <- unlist(lapply(seq(0 + alpha.scale, 1 - alpha.scale, alpha.scale), function(alpha) glm(response ~ bs(data[,1], knots=alpha, degree=1), family=family)$deviance))
		best.alpha <- which.min(deviances)*alpha.scale + alpha.scale
		model.2 <- glm(response ~ bs(data[,1], knots=best.alpha, degree=1), family=family)
		outputTable$m2aic[i] <- model.2$aic
		if (model.2$aic < (model.1$aic - 2)) {
			condition <- data[,1] <= best.alpha
			model.2.truncated <- glm(cbind(data[condition,2], ovules) ~ data[condition,1], family=family)
			params <- parameterEstimates(model.2.truncated)
			outputTable <- addEst(params[1], params[2], mean(data[data[,1]==1, 2])/ovules, best.alpha, i)
		} else {
			params <- parameterEstimates(model.1)
			outputTable <- addEst(params[1], params[2], 1, NA, i)
		}
	}
}
calcError <- function(actual, est) {
	(actual - est)/actual
}
conditions <- (outputTable$gx > outputTable$gs & outputTable$m > outputTable$gs)
outputTable$gsError <- ifelse(conditions, calcError(outputTable$gs, outputTable$gsEst), NA)
outputTable$gxError <- ifelse(conditions, calcError(outputTable$gx, outputTable$gxEst), NA)
outputTable$mError <- ifelse(conditions, calcError(outputTable$m, outputTable$mEst), NA)
outputTable$alphaError <- ifelse(conditions, calcError(outputTable$alpha, outputTable$alphaEst), NA)
write.table(outputTable, "output.txt", row.names=FALSE)