## Weighted Brokerage Score (Peeples and Haas 2013)

# initialize necessary libraries
library(compiler)
library(sciplot)

####################################################################################
## Intitialze primary functions  ###################################################
####################################################################################

# function for calculating total brokerage score by site
# function input x = rectangular matrix of similarities among pairs of site 
# scaled so that 0 is no similarity and 1 is perfect similarity
brk.ord2 <- function(x) { 
  ord <- combn(nrow(x),3) # find every possible combination of three nodes in matrix x (triads)
  ord <- t(ord) # covert to long format
  v <- matrix(0,nrow(x),1) # create blank matrix for output
  
  # the following for loop goes over every possible triad and caculates triad brokerage scores
  for (i in 1:nrow(ord)) {
    ab <- x[ord[i,1],ord[i,2]]
    bc <- x[ord[i,2],ord[i,3]]
    ac <- x[ord[i,1],ord[i,3]]
    out1 <- min(ab,ac)-bc 
    out2 <- min(ab,bc)-ac
    out3 <- min(ac,bc)-ab
    # set negative brokerage scores to 0
    if (out1<0) out1 <- 0
    if (out2<0) out2 <- 0
    if (out3<0) out3 <- 0
    # populate output matrix
    v[ord[i,1],1] <- v[ord[i,1],1]+out1
    v[ord[i,2],1] <- v[ord[i,2],1]+out2
    v[ord[i,3],1] <- v[ord[i,3],1]+out3}
  
  v <- round(v,3) # round to 3 digits for ease of plotting/reading
  rownames(v) <- rownames(x)
  colnames(v) <- c('total')
  return(v)}

brk.total <- cmpfun(brk.ord2) # compile function for faster analysis

#######

## This function outputs the triad brokerage score for every possible triad
brk.score <- function(x) {
  ord <- combn(nrow(x),3)
  ord <- t(ord)
  v <- matrix(0,nrow(ord),3)
  
  for (i in 1:nrow(ord)) {
    ab <- x[ord[i,1],ord[i,2]]
    bc <- x[ord[i,2],ord[i,3]]
    ac <- x[ord[i,1],ord[i,3]]
    outA <- min(ab,ac)-bc
    outB <- min(bc,ab)-ac
    outC <- min(bc,ac)-ab
    if (outA<0) outA <- 0
    if (outB<0) outB <- 0
    if (outC<0) outC <- 0
    v[i,1] <- outA
    v[i,2] <- outB
    v[i,3] <- outC}
  v <- round(v,3)
  return(v)}

brk.score2 <- cmpfun(brk.score) # compiled verison of function for faster running

#######

# function for calculating closure for every possible triad
brk.cls <- function(x) {
  ord <- combn(nrow(x),3)
  ord <- t(ord)
  v <- matrix(0,nrow(ord),3)
  for (i in 1:nrow(ord)) {
    v[i,1] <- x[ord[i,2],ord[i,3]]
    v[i,2] <- x[ord[i,1],ord[i,3]]
    v[i,3] <- x[ord[i,1],ord[i,2]]}
  v <- round(v,3)
  return(v)}

brk.closure <- cmpfun(brk.cls) # compiled verison of function for faster plotting


####################################################################################
## Run Brokerage Analysis  #########################################################
####################################################################################

### Example run
# read in sample ceramic file with sites as row names and columns represneting ceramic categories for 2 consecutive periods
cer <- read.csv('AD1200.csv',header=T,row.names=1) # Ceramic data for period 1
cer2 <- read.csv('AD1250.csv',header=T,row.names=1) # Ceramic data for period 2

# create Brainard-Robinson/Manhattan distance and rescale to similarity, rounded to 3 digits for both periods
sim <- round(1-(as.matrix(dist(prop.table(as.matrix(cer),1)*100,method='manhattan'))/200),3) # BR similarity for period 1
sim2 <- round(1-(as.matrix(dist(prop.table(as.matrix(cer2),1)*100,method='manhattan'))/200),3) # BR similarity for period 2

# run brokerage script and get total brokerage score by site
brk.out <- brk.total(sim) # Total brokerage score for period 1
brk.out2 <- brk.total(sim2) # Total brokerage score for period 2

### The following chunk of code calculates closure scores for two adjacent periods as described
### in Peeples and Haas 2013

# Calculate closure score from period 1 to period 2
# The next two objects represent sim and sim2 only for the site present in both periods
sim.a <- sim[is.element(rownames(sim),intersect(rownames(sim),rownames(sim2))),is.element(rownames(sim),intersect(rownames(sim),rownames(sim2)))]
sim2.a <- sim2[is.element(rownames(sim2),intersect(rownames(sim),rownames(sim2))),is.element(rownames(sim2),intersect(rownames(sim),rownames(sim2)))]

  brkAB <- brk.score2(sim.a) # calculate raw triad brokerage score for period 1
  brkAB.rnd <- round(brkAB*10,0)/10 # created rounded verison for plotting
  simA <- brk.closure(sim.a) # Similarities among nodes in each triad for period 1
  simB <- brk.closure(sim2.a) # Similarities among nodes in each triad for period 2
  simAB <- simB-simA # Closure score (change in similarity from period 1 to period 2)

# The next line replicates the plots for displaying closure by period from Peeples and Haas 2013
# As described in the article, this comparison excludes traids where the change in similarity
# from period 1 to period 2 was less than 1% (in either direction). This code could be modified
# to remove this feature by deleting the code within the square brackets []
lineplot.CI(brkAB.rnd[which(simAB>0.01 | simAB<(-0.01))],simAB[which(simAB>0.01 | simAB<(-0.01))],xlab='Triad Brokerage Score Time 1',ylab='Change in Similarity to Time 2')
  
  
