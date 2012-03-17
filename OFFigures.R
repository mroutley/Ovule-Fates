#####
# Figures for the ovule fates manuscript
#####
oldPar <- par()
OFEMPIRICALCODE <- function(){}
#bitmap("OFFig1.png", type="pnggray", res=1200)
pdf("OFEmpirical.pdf") #save figure as pdf
############
# Figure 1
# S:O empirical data for hand self & crossed and open pollinated
############
# Connect to database and extract data
library(RSQLite)
driver <- dbDriver("SQLite")
connection <- dbConnect(driver, dbname="OFData.db")
seedOvule.data <- dbGetQuery(connection, "SELECT genus, species, ovuleNumber, openSeedSet, handOutcrossSeedSet, handSelfSeedSet FROM floral_data")
# Clean up database connections
dbDisconnect(connection)
rm(connection)
dbUnloadDriver(driver)
rm(driver)
attach(seedOvule.data)
# Plot data
symbols <- c(21, 22, 23) #open, outcross, self
symbolColours <- c("black", "white", "white")
lineTypes <- c("solid", "dashed", "dotdash")
lineColours <- c("grey", "grey", "grey")
lineWidth <- 3
plot(ovuleNumber, openSeedSet, log="xy", xlab="Ovule production", ylab="Seed set", axes=FALSE, pch=symbols[1], bg=symbolColours[1])
axis(2, las=1, at=c(1, 10, 100, 1000, 10000), labels=c(1, 10, 100, 1000, 10000))
axis(1, at=c(1, 10, 100, 1000, 10000), labels=c(1, 10, 100, 1000, 10000))
points(ovuleNumber, handOutcrossSeedSet, pch=symbols[2], bg=symbolColours[2])
points(ovuleNumber, handSelfSeedSet, pch=symbols[3], bg=symbolColours[3])
legend(0, 4, c("Open pollinated", "Hand crossed", "Hand selfed"), pch=symbols, pt.bg=symbolColours, bty="n", pt.cex=1.5)
abline(log10(0.5), 1, lty="dashed", col="black", lwd=2) #50:50
abline(0, 1, lty="solid", col="black", lwd=2)#1:1
open.model <- lm(log10(openSeedSet) ~ log10(ovuleNumber))
cross.model <- lm(log10(handOutcrossSeedSet) ~ log10(ovuleNumber))
self.model <- lm(log10(handSelfSeedSet) ~ log10(ovuleNumber), data=subset(seedOvule.data, handSelfSeedSet!=0))
abline(open.model, lty=lineTypes[1], col=lineColours[1], lwd=lineWidth)
abline(cross.model, lty=lineTypes[2], col=lineColours[2], lwd=lineWidth)
abline(self.model, lty=lineTypes[3], col=lineColours[3], lwd=lineWidth)
detach(seedOvule.data)
rm(symbols, symbolColours, lineTypes, lineColours, lineWidth, open.model, cross.model, self.model)
dev.off()

#####
# Construct the general table of outcomes given a pollen table
#####
createFatesTable <- function(fatesTable) {
    
    attach(fatesTable)
    
    fatesTable$priorOvules <- ifelse(priorPollen < totalOvules, priorPollen, totalOvules)
    attach(fatesTable)
    
#    fatesTable$crossOvules <- ifelse(crossPollen + autonomousPollen + facilitatedSelfPollen < totalOvules - priorOvules, crossPollen, crossPollen / (crossPollen + autonomousPollen + facilitatedSelfPollen) * (totalOvules - priorOvules))
    fatesTable$crossOvules <- ifelse(crossPollen + autonomousPollen < totalOvules - priorOvules, crossPollen, crossPollen / (crossPollen + autonomousPollen) * (totalOvules - priorOvules))
    attach(fatesTable)
    
#    fatesTable$autonomousOvules <- ifelse(crossPollen + autonomousPollen + facilitatedSelfPollen < totalOvules - priorOvules, autonomousPollen, autonomousPollen / (crossPollen + autonomousPollen + facilitatedSelfPollen) * (totalOvules - priorOvules))
    fatesTable$autonomousOvules <- ifelse(crossPollen + autonomousPollen < totalOvules - priorOvules, autonomousPollen, autonomousPollen / (crossPollen + autonomousPollen) * (totalOvules - priorOvules))
    attach(fatesTable)
    
#    fatesTable$facilitatedSelfOvules <- ifelse(crossPollen + autonomousPollen + facilitatedSelfPollen < totalOvules - priorOvules, facilitatedSelfPollen, facilitatedSelfPollen / (crossPollen + autonomousPollen + facilitatedSelfPollen) * (totalOvules - priorOvules))
#    attach(fatesTable)
    
#    fatesTable$delayedOvules <- ifelse(delayedPollen < totalOvules - (priorOvules + crossOvules + autonomousOvules + facilitatedSelfOvules), delayedPollen, totalOvules - (priorOvules + crossOvules + autonomousOvules + facilitatedSelfOvules))
    fatesTable$delayedOvules <- ifelse(delayedPollen < totalOvules - (priorOvules + crossOvules + autonomousOvules), delayedPollen, totalOvules - (priorOvules + crossOvules + autonomousOvules))
    fatesTable$delayedOvules[fatesTable$delayedOvules < 0] <- 0
    attach(fatesTable)
    
#   fatesTable$selfedOvules <- priorOvules + autonomousOvules + facilitatedSelfOvules + delayedOvules
    fatesTable$selfedOvules <- priorOvules + autonomousOvules + delayedOvules

    fatesTable$outcrossedOvules <- crossOvules
    attach(fatesTable)
    fatesTable$unfertilizedOvules <- totalOvules - selfedOvules - outcrossedOvules
    fatesTable$unfertilizedOvules[fatesTable$unfertilizedOvules < 0] <- 0
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
    
    fatesTable # Return the table

}

#####
# Construct a table of pollen values
#####
createPollenTable <- function() {

#    facilitatedSelfPollenF <- function(table) { #additive
#        if (facilitatedMode == 1) {
#            return(pollenValue)
#        }
#        else if (facilitatedMode == 2) { #fixed total
#            return((facilitatedProportion * totalOvules) - table$crossPollen)
#        }
#        else if (facilitatedMode == 3) { #proportional
#            return(facilitatedSelfProportion * table$crossPollen)
#        }
#    }
    
    priorTable <- NULL
    priorTable$priorPollen <- sequence
    priorTable <- as.data.frame(priorTable)
    priorTable$crossPollen <- pollenValue
    priorTable$autonomousPollen <- pollenValue
    #priorTable$facilitatedSelfPollen <- facilitatedSelfPollenF(priorTable)
    #priorTable$facilitatedSelfPollen[priorTable$facilitatedSelfPollen < 0] <- 0
    priorTable$delayedPollen <- pollenValue
    attach(priorTable)
    #priorTable$po <- (priorPollen + crossPollen + autonomousPollen + facilitatedSelfPollen + delayedPollen) / totalOvules
    priorTable$po <- (priorPollen + crossPollen + autonomousPollen + delayedPollen) / totalOvules
    
    crossTable <- NULL
    crossTable$priorPollen <- rep(pollenValue, length(sequence))
    crossTable <- as.data.frame(crossTable)
    crossTable$crossPollen <- sequence
#    if (facilitatedMode == 2) {
#        crossTable$crossPollen[crossTable$crossPollen > facilitatedProportion * totalOvules] <- facilitatedProportion * totalOvules
#    }
    crossTable$autonomousPollen <- pollenValue
    #crossTable$facilitatedSelfPollen <- facilitatedSelfPollenF(crossTable)
    #crossTable$facilitatedSelfPollen[crossTable$facilitatedSelfPollen < 0] <- 0
    crossTable$delayedPollen <- pollenValue
    attach(crossTable)
#    crossTable$po <- (priorPollen + crossPollen + autonomousPollen + facilitatedSelfPollen + delayedPollen) / totalOvules
    crossTable$po <- (priorPollen + crossPollen + autonomousPollen + delayedPollen) / totalOvules
    
    rbind(priorTable, crossTable)
    
}

# Define constants
totalOvules <- 100
pollenValue <- 0.1 * totalOvules
sequence <- 0:(2 * totalOvules)
#facilitatedModes <- c("additive", "fixed total", "proportional")

gsList <- c(0.3, 0.7)
gx <- 0.9
kappa <- 0.7
m <- 0.8
ds <- 0.7
dx <- 0.9
#facilitatedMode <- 1
#facilitatedProportion <- 0.8
#facilitatedSelfProportion <- 0.5

pollenTable <- createPollenTable()
gs <- gsList[1]
gsLowFatesTable <- createFatesTable(pollenTable)
gs <- gsList[2]
gsHighFatesTable <- createFatesTable(pollenTable)

OFLIMITSCODE <- function(){}
############
# Figure 3
# Generates plots of the three limits on seed production (Panels A & D)
# Generates plots of pollen deposition on seed production and maternal fitness
############
panelLabels <- LETTERS[1:6]
pointLabels <- as.character(1:10)
lineWidth <- 2
lineTypes <- c("dashed", "solid") # prior, outcross
offset <- -5
#quartz(height=9, width=6)
pdf("OFLimits.pdf", height=9, width=6) # save figure as pdf
#bitmap("OFFig3.png", height=9, width=6, type="pnggray", res=1200)
#par(oma=c(1, 0, 0, 0))
par(mfrow=c(3,2), mex=1)
par(mar=c(4, 4, 1, 1))
LIMITS <- function(){}
limitsPlot <- function (gs, gx, kappa, m) {
	plot(c(0,1.1), c(0,1.1), type="n", lwd=2,  xlab="Proportion of self-fertilized ovules", ylab="Proportion of cross-fertilized ovules", axes=FALSE, xaxs="i", yaxs="i", main="") # Maximum ovule fertilization
	axis(side=1, at=c(0, m, 1), labels=c("0", "m", "1"))
	axis(side=2, at=c(0, m, 1), labels=c("0", "m", "1"), las=1)
	#box()
	polygon(c(0,1,0), c(0,0,1), col="lightgrey", border=NA) # Pollen limitation
	fsstar <- (gx - m) / (gx - gs)
	fxstar <- 1 - fsstar
	polygon(c(0, fsstar, 0), c(m / gx, fxstar, 1), col="slategrey", border=NA) # Resource limitation
	lines(c(fsstar,1), c(fxstar, 0), lwd=3) # Ovule limitation
	lines(c(0, 1), c((m / gx), ((m - gs) / gx)), lty="dashed", lwd=2) # Resource competition
	text(0.30, 0.08, "0")
}
addPollenLines <- function() {
	addSelfLine <- c(3 * pollenValue, pollenValue, totalOvules-pollenValue, pollenValue)/totalOvules
	segments(addSelfLine[1], addSelfLine[2], addSelfLine[3], addSelfLine[4], lty=lineTypes[1])
	addCrossLine <- c(3 * pollenValue, pollenValue,  3 * pollenValue, totalOvules-3*pollenValue)/totalOvules
	segments(addCrossLine[1], addCrossLine[2], addCrossLine[3], addCrossLine[4], lty=lineTypes[1])
}
PanelA <- function(){}
limitsPlot(gsList[1], gx, kappa, m)
addPollenLines()
points(c(0.9, 0.92, 1, 0.3, 0.32, (m-gx)/(gsList[1]-gx)), c(0.1, 0.1, 0, 0.7, 0.7, 1-(m-gx)/(gsList[1]-gx))-(offset/totalOvules), pch=pointLabels[1:6])
mtext(panelLabels[1], side=2, las=2, at=1.1, line=1)
PanelD <- function(){}
limitsPlot(gsList[2], gx, kappa, m)
addPollenLines()
points(c(0.9, 0.92, 1, 0.3, 0.3, 0.3+0.05), c(0.1, 0.1, 0, (m-0.3*gsList[2])/gx, 0.7, 0.7)-(offset/totalOvules), pch=pointLabels[1:6])
mtext(panelLabels[4], side=2, las=2, at=1.1, line=1)
rm(addPollenLines)
#####
# Seed&Fitness
#####
SEEDFITNESS <- function(){}
xAxisSequence <- seq(0.5, 2.0, 0.5)
markLine <- "I"
poMarks <- c(totalOvules, totalOvules + pollenValue, totalOvules + 3 * pollenValue)/totalOvules # no pollen limitation, delayed pre-empted, all prior
drawAxes <- function(gsValue) {
    axis(1, at=xAxisSequence, labels=(xAxisSequence))
    axis(2, las=1, at=totalOvules*c(0, gsList[gsValue], m, gx), labels=c("0", "gs", "m", "gx"))
    box()
}
makeMarks <- function(i, response, table) {
    attach(table)
    if (response=="seed") {
        points(rep(poMarks[i], 2), seed[po==poMarks[i]],  pch=markLine)
        points(rep(poMarks[i], 2), seed[po==poMarks[i]]+offset,  pch=pointLabels[c(i, i+3)])
    }
    if (response=="fitness") {
        points(rep(poMarks[i], 2), fitness[po==poMarks[i]],  pch=markLine)
        points(rep(poMarks[i], 2), fitness[po==poMarks[i]]+offset,  pch=pointLabels[c(i, i+3)])
    }
    detach(table)
}

PanelB <- function(){} #Seed; low gs
plot(gsLowFatesTable$seed[1:length(sequence)] ~ gsLowFatesTable$po[1:length(sequence)], type="l", lty=lineTypes[1], xlab="Pollen-ovule ratio", ylab="Seed-ovule ratio", las=1, xlim=range(50:200)/totalOvules, ylim=range(0:100), axes=FALSE)
lines(gsLowFatesTable$seed[(length(sequence)+1):(2*(length(sequence)+1))] ~ gsLowFatesTable$po[(length(sequence)+1):(2*(length(sequence)+1))], type="l", lty=lineTypes[2])
drawAxes(1)
makeMarks(1, "seed", gsLowFatesTable)
makeMarks(2, "seed", gsLowFatesTable)
points(poMarks[3], gsLowFatesTable$seed[po==poMarks[3]][1],  pch=markLine)
points(poMarks[3], gsLowFatesTable$seed[po==poMarks[3]][1]+offset,  pch=pointLabels[3])
points(1.55, gsLowFatesTable$seed[po==1.55][2],  pch=markLine)
points(1.55, gsLowFatesTable$seed[po==1.55][2]+offset,  pch=pointLabels[6])
mtext(panelLabels[2], side=2, las=2, at=totalOvules, line=1)

PanelE <- function(){} #Seed; high gs
plot(gsHighFatesTable$seed[1:length(sequence)] ~ gsHighFatesTable$po[1:length(sequence)], type="l", lty=lineTypes[1], xlab="Pollen-ovule ratio", ylab="Seed-ovule ratio", las=1, xlim=range(50:200)/totalOvules, ylim=range(0:100), axes=FALSE)
lines(gsHighFatesTable$seed[(length(sequence)+1):(2*(length(sequence)+1))] ~ gsHighFatesTable$po[(length(sequence)+1):(2*(length(sequence)+1))], type="l", lty=lineTypes[2])
drawAxes(2)
#makeMarks(1, "seed", gsHighFatesTable)
#makeMarks(2, "seed", gsHighFatesTable)
points(poMarks[1], gsHighFatesTable$seed[po==poMarks[1]][1],  pch=markLine)
points(poMarks[1], gsHighFatesTable $seed[po==poMarks[1]][1]+offset,  pch=pointLabels[1])
points(poMarks[2], gsHighFatesTable$seed[po==poMarks[2]][1],  pch=markLine)
points(poMarks[2], gsHighFatesTable $seed[po==poMarks[2]][1]+offset,  pch=pointLabels[2])
points(poMarks[3], gsHighFatesTable$seed[po==poMarks[3]][1],  pch=markLine)
points(poMarks[3], gsHighFatesTable $seed[po==poMarks[3]][1]+offset,  pch=pointLabels[3])
marks <- c(0.96, poMarks[1], poMarks[2])
points(marks, c(gsHighFatesTable$seed[po==marks[1]][2], gsHighFatesTable$seed[po==marks[2]][2], gsHighFatesTable$seed[po==marks[3]][2]), pch=markLine)
points(marks, c(gsHighFatesTable$seed[po==marks[1]][2], gsHighFatesTable$seed[po==marks[2]][2], gsHighFatesTable$seed[po==marks[3]][2])+offset, pch=pointLabels[4:6])
mtext(panelLabels[5], side=2, las=2, at=totalOvules, line=1)

PanelC <- function(){} #Fitness; low gs
plot(gsLowFatesTable$fitness[1:length(sequence)] ~ gsLowFatesTable$po[1:length(sequence)], type="l", lty=lineTypes[1], xlab="Pollen-ovule ratio", ylab="Fitness (genomes/ovules)", las=1, xlim=range(50:200)/totalOvules, ylim=range(0:100), axes=FALSE)
lines(gsLowFatesTable$fitness[(length(sequence)+1):(2*(length(sequence)+1))] ~ gsLowFatesTable$po[(length(sequence)+1):(2*(length(sequence)+1))], type="l", lty=lineTypes[2])
drawAxes(1)
makeMarks(1, "fitness", gsLowFatesTable)
makeMarks(2, "fitness", gsLowFatesTable)
points(poMarks[3], gsLowFatesTable$fitness[po==poMarks[3]][1],  pch=markLine)
points(poMarks[3], gsLowFatesTable$fitness[po==poMarks[3]][1]+offset,  pch=pointLabels[3])
points(1.55, gsLowFatesTable$fitness[po==1.55][2],  pch=markLine)
points(1.55, gsLowFatesTable$fitness[po==1.55][2]+offset,  pch=pointLabels[6])
mtext(panelLabels[3], side=2, las=2, at=totalOvules, line=1)

PanelF <- function(){} #Fitness; high gs
plot(gsHighFatesTable$fitness[1:length(sequence)] ~ gsHighFatesTable$po[1:length(sequence)], type="l", lty=lineTypes[1], xlab="Pollen-ovule ratio", ylab="Fitness (genomes/ovules)", las=1, xlim=range(50:200)/totalOvules, ylim=range(0:100), axes=FALSE)
lines(gsHighFatesTable$fitness[(length(sequence)+1):(2*(length(sequence)+1))] ~ gsHighFatesTable$po[(length(sequence)+1):(2*(length(sequence)+1))], type="l", lty=lineTypes[2])
drawAxes(2)
points(poMarks, c(gsHighFatesTable$fitness[po==poMarks[1]][1], gsHighFatesTable$fitness[po==poMarks[2]][1], gsHighFatesTable$fitness[po==poMarks[3]][1]), pch=markLine)
points(poMarks, c(gsHighFatesTable$fitness[po==poMarks[1]][1], gsHighFatesTable$fitness[po==poMarks[2]][1], gsHighFatesTable$fitness[po==poMarks[3]][1])+offset, pch=pointLabels[1:3])
marks <- c(0.96, poMarks[1], poMarks[2])
points(marks, c(gsHighFatesTable$fitness[po==marks[1]][2], gsHighFatesTable$fitness[po==marks[2]][2], gsHighFatesTable$fitness[po==marks[3]][2]), pch=markLine)
points(marks, c(gsHighFatesTable$fitness[po==marks[1]][2], gsHighFatesTable$fitness[po==marks[2]][2], gsHighFatesTable$fitness[po==marks[3]][2])+offset, pch=pointLabels[4:6])
mtext(panelLabels[6], side=2, las=2, at=totalOvules, line=1)
dev.off()
par(oldPar)
OFMIXEDMATINGCODE <- function(){}
############
# OptimalMixedMating: Generates fitness values for optimal mating system proportions
############
OPTIMALMIXEDMATING <- function(){}
# Generate data
optimalMixedMatingTable <- function(selfingModesProportion) {
	# Calculate pollen proportions
	crossPollen <- ((m - gs)/(gx - gs)) * pollenSequence
	selfPollen <- pollenSequence - crossPollen
	totalPollen <- crossPollen + selfPollen
	priorPollen <- selfingModesProportion[1] * selfPollen
	autonomousPollen <- selfingModesProportion[2] * selfPollen
	delayedPollen <- selfingModesProportion[3] * selfPollen
	table <- as.data.frame(cbind(priorPollen, autonomousPollen, crossPollen, delayedPollen, totalPollen))
	table
}
m <- 0.6
totalOvules <- 100
pollenSequence <- seq(50, 250, 5)
gs <- gsList[1]
equalProportions <- createFatesTable(optimalMixedMatingTable(c(1/3, 1/3, 1/3)))#prior, autonomous, delayed
allPrior <- createFatesTable(optimalMixedMatingTable(c(1, 0, 0)))#prior, autonomous, delayed, m, gs, gx
allAutonomous <- createFatesTable(optimalMixedMatingTable(c(0, 1, 0)))#prior, autonomous, delayed, m, gs, gx
allDelayed <- createFatesTable(optimalMixedMatingTable(c(0, 0, 1)))#prior, autonomous, delayed, m, gs, gx
conditions  <- c(selfedEmbryos==0)
allDelayed$seed[conditions] <- m * totalOvules
allDelayed$fitness[conditions] <- dx * m * totalOvules
# Graph code
lineTypes <- c("dashed", "dotdash", "solid", "dotted", "twodash") # prior, autonomous, outcross, delayed, mixed
symbolTypes <- c(24, 23, 22, 25, 21) # prior, autonomous, outcross, delayed, mixed
lineWidth <- 1
pollenTypeNames <- c("Prior pollen", "Simultaneous pollen", "Cross pollen", "Delayed pollen", "Equal mix of selfed")
labelSequence <- seq(50, 250, 50)
panelLabels <- LETTERS[1:2]
symbolColour <-c("white")
#bitmap("OFFig5.png", height=8, width=7, type="pnggray", res=1200)
pdf("OFMixedMating.pdf", height=8, width=7)
#quartz(height=8, width=7)
#par(mfrow=c(2,1), mar=c(2, 5, 1, 1) + 0.1, mex=0.8, oma=c(2, 0, 0, 0))
par(mfrow=c(2,1), mex=1, oma=c(2, 0, 0, 0))
par(mar=c(4, 4, 1, 1))
# Seed
plot(equalProportions$totalPollen, equalProportions$seed, type="o", lty=lineTypes[3], pch=symbolTypes[5], lwd=lineWidth, xlab="", ylab="Seed-ovule ratio", las=1, axes=FALSE, ylim=range(gs*totalOvules, m*totalOvules), bg=symbolColour)
axis(1, at=labelSequence, labels=labelSequence/totalOvules)
axis(2, las=1, at=c(gs*totalOvules, m*totalOvules), labels=c("gs", "m"))
mtext(panelLabels[1], side=2, las=1, line=2, at=m*totalOvules)
lines(allPrior$totalPollen, allPrior$seed, lty=lineTypes[3], lwd=lineWidth)
lines(allAutonomous$totalPollen, allAutonomous$seed, lty=lineTypes[3], lwd=lineWidth)
lines(allDelayed$totalPollen, allDelayed$seed, lty=lineTypes[3], lwd=lineWidth)
#segments(50, (gx*50), (m/gx)*totalOvules, m*totalOvules, lty=lineTypes[3], lwd=lineWidth)
#segments((m/gx)*totalOvules, m*totalOvules, 250, m*totalOvules, lty=lineTypes[3], lwd=lineWidth)
points(allPrior$totalPollen, allPrior$seed, pch=symbolTypes[1], bg=symbolColour)
points(allAutonomous$totalPollen, allAutonomous$seed, pch=symbolTypes[2], bg=symbolColour)
points(allDelayed$totalPollen, allDelayed$seed, pch=symbolTypes[4], bg=symbolColour)
legend(180, 50, legend=pollenTypeNames[-3], pch=symbolTypes[-3], bty="n", cex=1)
# Fitness
plot(equalProportions$totalPollen, equalProportions$fitness, type="o", lty=lineTypes[3], pch=symbolTypes[5], lwd=lineWidth, xlab="", ylab="Fitness (genomes/ovule)", las=1, axes=FALSE, ylim=range(gs*totalOvules, m*totalOvules+10), bg=symbolColour)
axis(1, at=labelSequence, labels=labelSequence/totalOvules)
axis(2, las=1, at=c(gs*totalOvules, m*totalOvules, m*totalOvules+10), labels=c("gs", "m", ""))
mtext(panelLabels[2], side=2, las=1, line=2, at=m*totalOvules+10)
lines(allPrior$totalPollen, allPrior$fitness, lty=lineTypes[3], lwd=lineWidth)
lines(allAutonomous$totalPollen, allAutonomous$fitness, lty=lineTypes[3], lwd=lineWidth)
lines(allDelayed$totalPollen, allDelayed$fitness, lty=lineTypes[3], lwd=lineWidth)
#segments(50, (gx*dx*50), (m/gx)*totalOvules, dx*m*totalOvules, lty=lineTypes[3], lwd=lineWidth)
#segments((m/gx)*totalOvules, dx*m*totalOvules, 250, dx*m*totalOvules, lty=lineTypes[3], lwd=lineWidth)
points(allPrior$totalPollen, allPrior$fitness, pch=symbolTypes[1], bg=symbolColour)
points(allAutonomous$totalPollen, allAutonomous$fitness, pch=symbolTypes[2], bg=symbolColour)
points(allDelayed$totalPollen, allDelayed$fitness, pch=symbolTypes[4], bg=symbolColour)
mtext("Pollen-ovule ratio", side=1, outer=TRUE)
dev.off()
par(oldPar)
OFTERNARYCODE <- function(){}
############
# TernaryPlots: Generates ternary plots of maximum fitness for values of autonomous, prior, & outcrossing
# Presence/absence of delayed is colour coded
############
TERNARYPLOTS <- function(){}
# Define constants & construct table
library(ade4)
totalOvules <- 10
m <- 0.6
sequence <- seq(0, 1*totalOvules, 1)
tableDimension <- length(sequence)^4
outputTable <- NULL
outputTable$priorPollen <- rep(sequence, each=tableDimension / length(sequence)^1)
outputTable <- as.data.frame(outputTable)
outputTable$crossPollen <- rep(sequence, each=tableDimension / length(sequence)^2)
outputTable$autonomousPollen <- rep(sequence, each=tableDimension / length(sequence)^3)
outputTable$delayedPollen <- rep(sequence, each=tableDimension / length(sequence)^4)
outputTable$crossPollen[crossPollen==0 & autonomousPollen==0]<- 0.01
gs <- gsList[1]
gsLow <- createFatesTable(outputTable)
gsLow <- gsLow[gsLow$fitness==max(gsLow$fitness), c(5:8)]
names(gsLow) <- c("P", "X", "A", "delayedOvules")
gs <- gsList[2]
gsHigh <- createFatesTable(outputTable)
gsHigh <- gsHigh[gsHigh$fitness==max(gsHigh$fitness), c(5:8)]
names(gsHigh) <- c("P", "X", "A", "delayedOvules")
pdf("OFTernary.pdf", height=8, width=4)
#quartz(height=8, width=4)
#bitmap("OFFig4.png", height=8, width=4, type="pnggray", res=1200)
par(mfrow=c(2,1), mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0), mex=0.1)
triangle.plot(gsLow[gsLow$delayedOvules < 1, 1:3], scale=F, show.position=F, cpoint=1.5, csub = 2, sub=c("A"), possub="topleft")
triangle.plot(gsHigh[gsHigh$delayedOvules < 1, 1:3], scale=F, show.position=F, cpoint=1.5, csub = 2, sub=c("B"), possub="topleft")
dev.off()
par(oldPar)