# Dynamic program for ovule fates
rm(list=ls()) # Clean up
#graphics.off()
#####
# Define constants
#####
totalOvules <- 50
gs <- 0.3; gx <- 0.9
ds <- 0.7; dx <- 0.9
m <- 0.8; kappa <- 0.7
tmax <- 8	# duration of anthesis
meanCrossPollenDeposition <- 6
probabilityOfVisit <- 0.5
selfPollenDeposition <- 5
facilitatedProportion <- 0.17 # Proportion of pollinator-delivered pollen that is self pollen
decision <- c(1, 2, 3) # Choose self, choose facilitated, choose outcross; used to store data as matrix rather than data frame
#####
# Ovule fates
#####
# Embryos
selfedEmbryosF <- function(table) {
    attach(table)
    gs * (selfedOvules)
}
crossedEmbryosF <- function(table) {
    attach(table)
    gx * crossedOvules
}
# Competition
selfResourceComp <- function(selfedEmbryos, crossedEmbryos) { # Weighted lottery
	p <- function(n, x) {
		((selfedEmbryos - x) * kappa)/((selfedEmbryos - x) * kappa + crossedEmbryos - (n - x - 1))
	}
	expectation <- 0 + p(1, 0)
	for (i in 2:(m * totalOvules)) {
		expectation <- expectation + p(i, expectation)
	}
	expectation
}
ksF <- function(table) {
    attach(table)
    ks <- ifelse((selfedEmbryos + crossedEmbryos <= m * totalOvules), 1, NA)
	conditions <- c(selfedEmbryos + crossedEmbryos > m * totalOvules)
	ks[conditions] <- selfResourceComp(selfedEmbryos[conditions], crossedEmbryos[conditions])/selfedEmbryos[conditions]
	conditions <- c(selfedEmbryos == 0)
	ks[conditions] <- 1
	ks
}
kxF <- function(table) {
    attach(table)
    kx <- ifelse((selfedEmbryos + crossedEmbryos <= m * totalOvules), 1, NA)
    conditions <- c(selfedEmbryos + crossedEmbryos > m * totalOvules)
	kx[conditions] <- (m * totalOvules - selfResourceComp(selfedEmbryos[conditions], crossedEmbryos[conditions]))/crossedEmbryos[conditions]
	conditions <- c(selfedEmbryos == 0)
	kx[conditions] <- 1
	kx
}
# Seeds
selfedSeedF <- function(table) {
    attach(table)
    selfedSeed <- ks * selfedEmbryos
    conditions <- c(crossedEmbryos == 0 & selfedEmbryos > m * totalOvules) # Ensures that pure selfers can't support > m*totalOvules seeds
    selfedSeed[conditions] <- m * totalOvules
    selfedSeed
}
crossedSeedF <- function(table) {
    attach(table)
    crossedSeed <- kx * crossedEmbryos
    conditions <- c(selfedEmbryos == 0 & crossedEmbryos > m * totalOvules) # Ensures that pure crossers can't support > m*totalOvules seeds
    crossedSeed[conditions] <- m  * totalOvules
    crossedSeed
}
selfedSeedlingF <- function(table) {
    attach(table)
    ds * selfedSeed
}
crossedSeedlingF <- function(table) {
    attach(table)
    dx * crossedSeed
}
# Fitness
fitnessF <- function(table) {
    attach(table)
    (2 * selfedSeedling) + crossedSeedling
}
#####
# Ovule fates table: calculate all values & extract as necessary
#####
fillOvuleFatesTable <- function(table) {
	# Ovule values are given
	# Embryo values
	table$selfedEmbryos <- selfedEmbryosF(table)
	table$crossedEmbryos <- crossedEmbryosF(table)
	# Competition
	table$ks <- ksF(table)
    table$kx <- kxF(table)
	table$selfedSeed <- selfedSeedF(table)
	table$crossedSeed <- crossedSeedF(table)
	table$seed <- table$selfedSeed + table$crossedSeed
	# Seedlings
	table$selfedSeedling <- selfedSeedlingF(table)
	table$crossedSeedling <- crossedSeedlingF(table)
	# Fitness
	table$fitness <- fitnessF(table)
	table
}
#####
# Create ovule fates table
#####
sequence <- seq(0, totalOvules, 1)
tableDimension <- length(sequence)^2 # Start with all possible combinations
ovuleFatesTable <- NULL
ovuleFatesTable$selfedOvules <- rep(sequence, each=tableDimension / length(sequence))
ovuleFatesTable <- as.data.frame(ovuleFatesTable)
ovuleFatesTable$crossedOvules <- rep(sequence, each=tableDimension / length(sequence)^2)
ovuleFatesTable <- ovuleFatesTable[ovuleFatesTable$crossedOvules+ovuleFatesTable$selfedOvules<=totalOvules,] # Remove impossible combinations
ovuleFatesTable <- fillOvuleFatesTable(ovuleFatesTable)
#####
# Fitness calculations
#####
fertilizationValueSelfing <- function(table) { # If the plant chooses to self, it realises this number of fertilizations
	attach(table)
	newSelf <- selfedOvules + selfPollenDeposition
	conditions <- c(newSelf + crossedOvules > totalOvules)
	newSelf[conditions] <- totalOvules - crossedOvules[conditions]
	newSelf
}
fertilizationValueCrossing <- function(table, mode) { # If the plant chooses not to self, a pollinator may arrive
	ifelse(mode==1, facilitatedProportion <- 0, facilitatedProportion) # Mode 1 means pure crossing
	table$newOvules <- meanCrossPollenDeposition * probabilityOfVisit
	attach(table)
	conditions <- c(newOvules + selfedOvules + crossedOvules > totalOvules)
	table$newOvules[conditions] <- totalOvules - (selfedOvules[conditions] + crossedOvules[conditions])
	newSelf <- selfedOvules + facilitatedProportion * table$newOvules
	newCross <- crossedOvules + (1 - facilitatedProportion) * table$newOvules
	round(cbind(newSelf, newCross))
}
extractFitness <- function(table) { # Input self and cross vectors
	fitnessVector <- NULL
	for (i in 1:dim(table)[1]) { # For each self (table[ , 1]) and cross (table[ , 2]) ovule combination, lookup the fitness value
		fitnessVector[i] <- fitnessTable[fitnessTable[ , 1]==table[ , 1][i] & fitnessTable[ , 2]==table[ , 2][i], 3]
	}
	fitnessVector
}
chooseStrategy <- function(table) {
	table$choice <- 0
	attach(table)
	conditions <- (chooseSelfFitness > chooseFacilitatedFitness & chooseSelfFitness > chooseCrossFitness)
	table$choice[conditions] <- decision[1]
	conditions <- (chooseFacilitatedFitness > chooseSelfFitness & chooseFacilitatedFitness > chooseCrossFitness)
	table$choice[conditions] <- decision[2]
	conditions <- (chooseCrossFitness > chooseFacilitatedFitness & chooseCrossFitness > chooseSelfFitness)
	table$choice[conditions] <- decision[3]
	table$choice
}
theLoop <- function() {
	stateTable <- ovuleFatesTable[ , c(1,2)] # Start with self and cross ovule combinations
	# Autonomous selfing
	stateTable$chooseSelfNewSelf <- fertilizationValueSelfing(stateTable)
	stateTable$chooseSelfNewCross <- stateTable$crossedOvules
	# Facilitated selfing
	facilitatedOvules <- fertilizationValueCrossing(stateTable, 0)
	stateTable$chooseFacilitatedNewSelf <- facilitatedOvules[ , 1]
	stateTable$chooseFacilitatedNewCross <- facilitatedOvules[ , 2]
	# Crossing
	crossOvules <- fertilizationValueCrossing(stateTable, 1)
	stateTable$chooseCrossNewSelf <- stateTable$selfedOvules
	stateTable$chooseCrossNewCross <- facilitatedOvules[ , 2]
	# Calculate the fitness of each strategy
	stateTable$chooseSelfFitness <- extractFitness(stateTable[ , c(3,4)])
	stateTable$chooseFacilitatedFitness <- extractFitness(stateTable[, c(5,6)])
	stateTable$chooseCrossFitness <- extractFitness(stateTable[ , c(7,8)])
	# Choose the best strategy
	stateTable$optimalStrategy <- chooseStrategy(stateTable[ , c(9:11)])
	# Store the best strategy for each ovule combination
	optimalTable <- NULL
	optimalTable <- stateTable[ , c(1, 2)]
	optimalTable$optimalStrategy <- stateTable$optimalStrategy
	optimalTable$fitness <- ovuleFatesTable[ , 12]
	for (i in 1:length(decision)) {
		conditions <- decision[i]
		optimalTable$fitness[optimalTable$optimalStrategy == conditions] <- stateTable$chooseSelfFitness[stateTable$optimalStrategy == conditions]
	}
	optimalTable$newSelf <- optimalTable[ , 1]
	optimalTable$newCross <- optimalTable[ , 2]
	# choose self
	conditions <- decision[1]
	optimalTable$newSelf[optimalTable$optimalStrategy == conditions] <- stateTable$chooseSelfNewSelf[stateTable$optimalStrategy == conditions]
	optimalTable$newCross[optimalTable$optimalStrategy == conditions] <- stateTable$chooseSelfNewCross[stateTable$optimalStrategy == conditions]
	# choose facilitated
	conditions <- decision[2]
	optimalTable$newSelf[optimalTable$optimalStrategy == conditions] <- stateTable$chooseFacilitatedNewSelf[stateTable$optimalStrategy == conditions]
	optimalTable$newCross[optimalTable$optimalStrategy == conditions] <- stateTable$chooseFacilitatedNewCross[stateTable$optimalStrategy == conditions]
	# choose cross
	conditions <- decision[3]
	optimalTable$newSelf[optimalTable$optimalStrategy == conditions] <- stateTable$chooseCrossNewSelf[stateTable$optimalStrategy == conditions]
	optimalTable$newCross[optimalTable$optimalStrategy == conditions] <- stateTable$chooseCrossNewCross[stateTable$optimalStrategy == conditions]
	#optimalTable$time <- t
	#print(optimalTable)
	optimalTable
}
begin <- Sys.time()
anthesisArray <- NULL
#t=T-1
fitnessTable <- ovuleFatesTable[ , c(1, 2, 12)] # Fitness for t=tmax
t <- tmax
optimalTable <- theLoop()
anthesisArray <- optimalTable
fitnessTable <- optimalTable[ , c(1, 2, 4)]
for (t in (tmax-1):1) {
	#print(t)
	loopTable <- theLoop()
	anthesisArray <- cbind(loopTable, anthesisArray)
	fitnessTable <- loopTable[ , c(1, 2, 4)] # Recreate the fitness table with current expectations
}
print(Sys.time() - begin)
anthesisArray <- as.matrix(anthesisArray)
dim(anthesisArray) <- c(dim(optimalTable)[1], dim(optimalTable)[2], tmax)
plotColour <- c("black", "red", "blue", "green") # Choose none, self, facilitated, outcross
#quartz(height=8, width=8)
#windows(height=6, width=6)
#par(mfrow=c(2, 2), oma=c(0, 0, 2, 0))
decisionPlot <- function(time) {
	data <- anthesisArray[ , , time][ , 1:3]
	plot(data[, 1], data[, 2], pch=21, bg=plotColour[data[, 3]+1], col=plotColour[data[, 3]+1], xlab="Self-fertilized ovules", ylab="Cross-fertilized ovules", main=paste("Time: ", time))
}
#decisionPlot(1)
#decisionPlot(round(tmax/2))
#decisionPlot(tmax)
decisionTable <- NULL
decisionTable <- data.frame(decision=numeric(tmax), selfedOvules=numeric(tmax), crossedOvules=numeric(tmax))
decisionTable[1, ] <- anthesisArray[ , , 1][1, c(3, 5, 6)] # Start at the origin
ifelse(decisionTable[1, 1]==decision[2], decisionTable[1, ] <- round(c(decision[2], facilitatedProportion*rpois(1, meanCrossPollenDeposition), (1-facilitatedProportion)*rpois(1, meanCrossPollenDeposition))), 0)
for (i in 2:tmax) { # Make the optimal decision at each time step
	decisionTable[i, ] <- anthesisArray[ , , i][anthesisArray[ , , i][, 1]==decisionTable[i-1, 2] & anthesisArray[ , , i][, 2]==decisionTable[i-1, 3]][c(3, 5, 6)]
	if (decisionTable[i, 1]==decision[2]) {
		currentOvules <- decisionTable[i - 1, 2] + decisionTable[i - 1, 3]
		newOvules <- rpois(1, meanCrossPollenDeposition)
		ifelse(currentOvules + newOvules > totalOvules, newOvules <- totalOvules - currentOvules, 0)
		decisionTable[i, ] <- round(c(decision[2], decisionTable[i-1, 2] + facilitatedProportion * newOvules, decisionTable[i-1, 3] + (1 - facilitatedProportion) * newOvules), 0)		
	}
}
plot(c(0, decisionTable$selfedOvules), c(0, decisionTable$crossedOvules), xlim=range(0:totalOvules) , ylim=range(0:totalOvules), pch=21, bg="black", type="o", xlab="Self-fertilized ovules", ylab="Cross-fertilized ovules")
points(ovuleFatesTable[ovuleFatesTable$fitness==max(ovuleFatesTable$fitness),][1:2], pch=21, bg=plotColour[2], col=plotColour[2]) # Add the optimal ovule combination
title(main=paste("self: ", selfPollenDeposition, "cross: ", meanCrossPollenDeposition, "facilitated: ", facilitatedProportion, "visit: ", probabilityOfVisit ), outer=FALSE)
