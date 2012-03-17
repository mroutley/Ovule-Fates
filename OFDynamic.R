# Define constants
decision <- c(1, 2, 3) # Choose self, choose facilitated, choose outcross; used to store data as matrix rather than data frame
totalOvules <- 50
gs <- 0.3; gx <- 0.9
ds <- 0.7; dx <- 0.9
m <- 0.8; kappa <- 0.7
tmax <- 10	# duration of anthesis
iterations <- 4 # How many trajectories to consider

constructFatesTable <- function() { # Generates the ovule fates table
    
    # Start with ovule combinations
    sequence <- seq(0, totalOvules, 1)
    tableDimension <- length(sequence)^2 # Start with all possible combinations
    fatesTable <- NULL
    fatesTable$selfedOvules <- rep(sequence, each=tableDimension / length(sequence))
    fatesTable <- as.data.frame(fatesTable)
    fatesTable$outcrossedOvules <- rep(sequence, each=tableDimension / length(sequence)^2)
    fatesTable <- fatesTable[fatesTable$outcrossedOvules + fatesTable$selfedOvules <= totalOvules, ] # Remove impossible combinations
    attach(fatesTable)
    
    # Embryos
    fatesTable$selfedEmbryos <- gs * selfedOvules
    fatesTable$outcrossedEmbryos <- gx * outcrossedOvules
    attach(fatesTable)
    
    # Competition
    selfResourceComp <- function(selfedEmbryos, outcrossedEmbryos) { # Weighted lottery
    	p <- function(n, x) {
    		((selfedEmbryos - x) * kappa)/((selfedEmbryos - x) * kappa + outcrossedEmbryos - (n - x - 1))
    	}
    	expectation <- 0 + p(1, 0)
    	for (i in 2:(m * totalOvules)) {
    		expectation <- expectation + p(i, expectation)
    	}
    	expectation
    }
    
    fatesTable$ks <- ifelse((selfedEmbryos + outcrossedEmbryos <= m * totalOvules), 1, NA)
    conditions <- c(selfedEmbryos + outcrossedEmbryos > m * totalOvules)
    fatesTable$ks[conditions] <- selfResourceComp(selfedEmbryos[conditions], outcrossedEmbryos[conditions])/selfedEmbryos[conditions]
    fatesTable$ks[outcrossedEmbryos == 0] <- 1
    fatesTable$ks[selfedEmbryos == 0] <- 1
    
    fatesTable$kx <- ifelse((selfedEmbryos + outcrossedEmbryos <= m * totalOvules), 1, NA)
    #conditions <- c(selfedEmbryos + outcrossedEmbryos > m * totalOvules)
    fatesTable$kx[conditions] <- (m * totalOvules - selfResourceComp(selfedEmbryos[conditions], outcrossedEmbryos[conditions]))/outcrossedEmbryos[conditions]
    fatesTable$kx[outcrossedEmbryos == 0] <- 1
    fatesTable$kx[selfedEmbryos == 0] <- 1
    attach(fatesTable)
    
    # Seeds
    fatesTable$selfedSeed <- ks * selfedEmbryos
    conditions <- c(outcrossedEmbryos == 0 & selfedEmbryos > m * totalOvules)
    fatesTable$selfedSeed[conditions] <- m * totalOvules
    
    fatesTable$outcrossedSeed <- kx * outcrossedEmbryos
    conditions <- c(selfedEmbryos == 0 & outcrossedEmbryos > m * totalOvules)
    fatesTable$outcrossedSeed[conditions] <- m  * totalOvules
    attach(fatesTable)
    
    fatesTable$seed <- selfedSeed + outcrossedSeed
    attach(fatesTable)
    
    # Seedlings
    fatesTable$selfedSeedling <- ds * selfedSeed
    fatesTable$outcrossedSeedling <- dx * outcrossedSeed
    attach(fatesTable)
    
    # Fitness
    fatesTable$fitness <- (2 * selfedSeedling) + outcrossedSeedling
    
    fatesTable

}

extractFitness <- function(table) { # Input self and cross vectors
	fitnessVector <- NULL
	for (i in 1:dim(table)[1]) { # For each self (table[ , 1]) and cross (table[ , 2]) ovule combination, lookup the fitness value
		fitnessVector[i] <- fitnessTable[fitnessTable[ , 1]==table[ , 1][i] & fitnessTable[ , 2]==table[ , 2][i], 3]
	}
	fitnessVector
}

timeInterval <- function() { # For each time interval, we determine the optimal strategy
    
    timeIntervalTable <- fatesTable[ , c(1,2)] # Start with self and cross ovule combinations
    
    # For autonomous pollination
    timeIntervalTable$chooseAutonomousNewSelf <- selfedOvules + selfPollenDeposition
    attach(timeIntervalTable)
    conditions <- c(chooseAutonomousNewSelf + outcrossedOvules > totalOvules)
	timeIntervalTable$chooseAutonomousNewSelf[conditions] <- totalOvules - (outcrossedOvules[conditions])
    timeIntervalTable$chooseAutonomousNewCross <- outcrossedOvules

	# For facilitated pollination
	timeIntervalTable$chooseCrossNewSelf <- selfedOvules + facilitatedProportion * meanCrossPollenDeposition * probabilityOfVisit
	timeIntervalTable$chooseCrossNewCross <- outcrossedOvules +  meanCrossPollenDeposition * probabilityOfVisit
	attach(timeIntervalTable)
	conditions <- c(chooseAutonomousNewSelf + chooseCrossNewCross > totalOvules)
	timeIntervalTable$chooseCrossNewSelf[conditions] <- selfedOvules[conditions] + facilitatedProportion * (totalOvules - (selfedOvules[conditions] + outcrossedOvules[conditions]))
	timeIntervalTable$chooseCrossNewCross[conditions] <- outcrossedOvules[conditions] + (1 - facilitatedProportion) * (totalOvules - (selfedOvules[conditions] + outcrossedOvules[conditions]))
	timeIntervalTable$chooseCrossNewSelf <- round(timeIntervalTable$chooseCrossNewSelf, 0)
	timeIntervalTable$chooseCrossNewCross <- round(timeIntervalTable$chooseCrossNewCross, 0)
	attach(timeIntervalTable)
		
	# Calculate the fitness of each strategy
    timeIntervalTable$chooseAutonomousFitness <- extractFitness(timeIntervalTable[ , c(3,4)])
    timeIntervalTable$chooseFacilitatedFitness <- extractFitness(timeIntervalTable[ , c(5,6)])
	attach(timeIntervalTable)
	    
    # Choose the best strategy
    timeIntervalTable$optimalStrategy <- ifelse(timeIntervalTable$chooseAutonomousFitness > timeIntervalTable$chooseFacilitatedFitness, decision[1], decision[2])

    # Store the best strategy for each ovule combination
    timeIntervalOptimalTable <- NULL
    timeIntervalOptimalTable <- timeIntervalTable[ , c(1, 2, 9)]
	timeIntervalOptimalTable$fitness <- timeIntervalTable$chooseFacilitatedFitness
	timeIntervalOptimalTable$fitness[timeIntervalOptimalTable$optimalStrategy == decision[1]] <- timeIntervalTable$chooseAutonomousFitness[timeIntervalTable$optimalStrategy == decision[1]]

	timeIntervalOptimalTable$newSelf <- timeIntervalTable$chooseFacilitatedNewSelf
	timeIntervalOptimalTable$newSelf[timeIntervalOptimalTable$optimalStrategy == decision[1]] <- timeIntervalTable$chooseAutonomousNewSelf[timeIntervalTable$optimalStrategy == decision[1]]

	timeIntervalOptimalTable$newCross <- timeIntervalTable$chooseFacilitatedNewCross
	timeIntervalOptimalTable$newCross[timeIntervalOptimalTable$optimalStrategy == decision[1]] <- timeIntervalTable$chooseAutonomousNewCross[timeIntervalTable$optimalStrategy == decision[1]]

    timeIntervalOptimalTable # Return optimal table for the next iteration
    
}

fatesTable <- constructFatesTable()

quartz(height=7, width=5)
#windows(height=7, width=5)
#pdf("OFTrajectory.pdf", height=7, width=5) # save figure as pdf
par(mfrow=c(3,2), mex=1)
par(mar=c(4, 4, 1, 1), oma=c(1, 1.1, 0, 0))

panelLetters <- LETTERS[1:6]
cross <- c(   5,  20,   5,  20,   5,   1)
self <-  c(  20,   5,  20,   5,   5,  20)
fac <-   c(0.17,   0,   0,   0,   0,   0)
visit <- c( 0.5, 0.5, 0.5, 0.1, 0.9, 0.1)

for (x in 1:6) {
	panelLetter <- panelLetters[x]
	meanCrossPollenDeposition <- cross[x]
	probabilityOfVisit <- visit[x]
	selfPollenDeposition <- self[x]
	facilitatedProportion <- fac[x]
	# Simulation starts here
    anthesisArray <- NULL # Stores the simulation results in a matrix
    t <- tmax
    fitnessTable <- fatesTable[ , c(1, 2, 12)] # Fitness for t=tmax
    timeIntervalOptimalTable <- timeInterval()
    anthesisArray <- timeIntervalOptimalTable
    fitnessTable <- timeIntervalOptimalTable[ , c(1, 2, 4)]
    
    for (t in (tmax-1):1) { # Cycle through the remaining time intervals
        timeIntervalTable <- timeInterval()
        anthesisArray <- cbind(timeIntervalTable, anthesisArray)
        fitnessTable <- timeIntervalTable[ , c(1, 2, 4)] # Recreate the fitness table with current expectations
    }
    
    anthesisArray <- as.matrix(anthesisArray)
    dim(anthesisArray) <- c(dim(timeIntervalOptimalTable)[1], dim(timeIntervalOptimalTable)[2], tmax)
    
    makeDecisions <- function() {    
        decisionTable <- NULL # At each time interval, make the optimum decision
        decisionTable <- data.frame(decision=numeric(tmax), selfedOvules=numeric(tmax), outcrossedOvules=numeric(tmax))
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
       	decisionTable
    }
    #####
    # This code was used to create a trajectory matrix
    # Median and quartile ranges were then extracted
    # However, this removes contingency from the trajectories
    #####
    #decisionArray <- NULL
    #decisionArray <- makeDecisions()
    #while(dim(decisionArray)[2]<3*iterations) {decisionArray <- cbind(decisionArray, makeDecisions())}
    #decisionArray <- as.matrix(decisionArray)
    #dim(decisionArray) <- c(tmax, 3, iterations)
    #dimnames(decisionArray) <- list(1:tmax, c("decision", "selfed", "crossed"), 1:iterations)
    
    #median <- NULL
    #for (i in 1:tmax) {
    #	median$selfed[i+1] <- summary(decisionArray[i,2,])[3]
    #	median$crossed[i+1] <- summary(decisionArray[i,3,])[3]
    #}
    #median <- as.data.frame(median)
    #median[1, ] <- c(0, 0)
    #lower <- NULL
    #for (i in 1:tmax) {
    #	lower$selfed[i+1] <- summary(decisionArray[i,2,])[2]
    #	lower$crossed[i+1] <- summary(decisionArray[i,3,])[2]
    #}
    #lower <- as.data.frame(lower)
    #lower[1, ] <- c(0, 0)
    #upper <- NULL
    #for (i in 1:tmax) {
    #	upper$selfed[i+1] <- summary(decisionArray[i,2,])[5]
    #	upper$crossed[i+1] <- summary(decisionArray[i,3,])[5]
    #}
    #upper <- as.data.frame(upper)
    #upper[1, ] <- c(0, 0)

    #plot(median, xlim=range(0:totalOvules) , ylim=range(0:totalOvules), pch=21, bg="black", type="o", xlab="", ylab="", las=1)
    #points(lower, pch=25)
    #lines(lower, lty="dotted")
    #points(upper, pch=24)
    #lines(upper, lty="dotted")
    #####
    # I also considered plotting ovule fertilizations through time as bar plots
    #####
    #barplot(as.matrix(t(median)))
    
    # Currently, opting for representative trajectories
    plot(rbind(c(0, 0), makeDecisions()[,2:3]), xlim=range(0:totalOvules) , ylim=range(0:totalOvules), pch=21, bg="black", type="o", xlab="", ylab="", axes=FALSE)
    axis(1)
    axis(2, las=1)
    for (i in 1:(iterations-1)) {
        points(rbind(c(0, 0), makeDecisions()[,2:3]), pch=21, bg="black", type="o")
    }
    abline(coef=c(totalOvules,-1), lty="dashed")
    points(fatesTable[fatesTable$fitness==max(fatesTable$fitness),][1:2], pch=22)
    mtext(panelLetter, side=2, las=2, at=totalOvules, line=3)
}
mtext("Self-fertilized ovules", side=1, outer=TRUE)
mtext("Cross-fertilized ovules", side=2, outer=TRUE)
#dev.off()
