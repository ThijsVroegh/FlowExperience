## go to my website for a more-detailed explanation
## >> https://tvroegh.netlify.app/ <<

## -----------------------------------------------------------------------------
library(qgraph)
library(igraph)
library(bnlearn)
library(Rgraphviz)
library(tidyverse)
library(gplots)
library(plyr)
library(magrittr)
library(e1071)
library(parallel)
library(corrplot)
library(CTT)
library(networktools)
library(matrixcalc)
library(Matrix)
library(reshape2)
library(MASS)
library(lavaan)
library(lavaanPlot)
library(rcausal)
library(DOT)
library(pcalg)
library(tpc) # enables 'time-ordered data' with the pc-algorithm 
library(tidySEM)
library(DiagrammeR)

## Set seed
set.seed(123)

## -----------------------------------------------------------------------------
# first, read in factor loadings from p.142
L <- as.matrix(c(0.79, 0.65, 0.57, 0.46, 0.68, 0.85, 0.58, 0.19, 0.78))
colnames(L) <- paste0("F", 1)
rownames(L) <- paste0("V", 1:9)

## -----------------------------------------------------------------------------
F <- 1
F <- as.matrix(F)
colnames(F) <- rownames(F) <- paste0("F", 1)

## -----------------------------------------------------------------------------
R <- L %*% F %*% t(L)
diag(R) <- 1
mycor <- round(as.matrix(R),2)

## -----------------------------------------------------------------------------
colnames(mycor)<-c("balance", "merging", "goals", "feedback", "concentration","control","selfawareness", "time", "autotelic") 

rownames(mycor) <- colnames(mycor)

#check pd matrix
is.positive.definite(mycor)

## -----------------------------------------------------------------------------
# means and standard deviations from Wrigley (2005), p.143
mu <- c(3.610, 3.370, 4.120, 3.780, 3.530, 3.280,3.120,3.280,3.510)
stddev <- c(.700,.740,.570,.710,.860,.750,.950,.940,.940)

covMat <- stddev %*% t(stddev) * mycor
mycov <- round(as.matrix(covMat),2)
          
#Define variable names
colnames(mycov) = c("balance", "merging", "goals", "feedback", "concentration","control","selfawareness", "time", "autotelic") 
rownames(mycov) <- colnames(mycov)

#check pd matrix
is.positive.definite(mycov)

## -----------------------------------------------------------------------------
## Plot correlation and covariance matrices
par(mfrow=c(1,2))
corrplot(mycor, method = "color", addCoef.col = "black",
         number.cex = .7, cl.pos = "n", diag=T, insig = "blank")
title(main = "Estimated correlation matrix")

corrplot(mycov, method = "color", addCoef.col = "black",
         number.cex = .7, cl.pos = "n", diag=T, insig = "blank")
title(main="Estimated covariance matrix")

## -----------------------------------------------------------------------------
HS.model <- "flow = ~ balance + merging + goals + feedback + concentration + control + selfawareness + time + autotelic"

# fit the model
fit <- cfa(HS.model, 
           sample.cov  = covMat, 
           sample.nobs = 236)

# display summary output
summary(fit, 
        standardized = TRUE,
        fit.measures=TRUE,
        rsquare = TRUE)

# output: Latent variables: 'Estimate' are unstandardized, 'Std.all' are
# the ‘completely standardized solution’ (matching those on p.142)
# get the test statistic for the original sample

#standardizedsolution(fit)

## -----------------------------------------------------------------------------
lavaanPlot::lavaanPlot(model = fit,
           node_options = list(shape = "box", fontname = "Helvetica"),
           edge_options = list(color = "grey"),
           coefs = TRUE,
           stand = TRUE,
           stars = "latent")

## -----------------------------------------------------------------------------

community_structure <- c('1','2','1','1','2','2','2','2','2') 
is.vector(community_structure)

networkflow236 <- EBICglasso(mycor, n = 236,gamma = 0.5, nlambda = 100)
is.vector(community_structure)  

bridge_centrality_236 <- bridge(networkflow236, communities = community_structure, 
                                directed = FALSE)

plot(bridge_centrality_236, include=c("Bridge Expected Influence (1-step)"), 
     theme_bw = F, zscore = TRUE)

## -----------------------------------------------------------------------------
bridge_strength_236 <- bridge_centrality_236$`Bridge Expected Influence (1-step)`

top_bridges_236 <- names(bridge_strength_236[bridge_strength_236>quantile(bridge_strength_236, probs=0.80, na.rm=TRUE)])

# Now create a new community vector where bridges are their own community
bridge_num_flow_236 <- which(names(bridge_strength_236) %in% top_bridges_236)

new_communities <- vector()

for(i in 1:length(bridge_strength_236)) {
  if(i %in% bridge_num_flow_236) {
    new_communities[i] <- "Bridge Nodes"
  } else {new_communities[i] <- community_structure[i]}
}

new_communities

# And now use that community vector to make groups in qgraph 
qgraph(networkflow236, layout="spring", 
       groups = new_communities, color=c("yellow", "turquoise", "salmon"))

## -----------------------------------------------------------------------------
tetradrunner.getAlgorithmDescription(algoId = 'fges')

## -----------------------------------------------------------------------------
dat2 <- round(mvrnorm(n = 236, mu = mu, Sigma = covMat, empirical = T),2)

# checking the means
colMeans(dat2)

# Define variable names
colnames(dat2)=c("balance", "merging", "goals", "feedback", "concentration","control","selfawareness", "time", "autotelic") 

head(dat2)

dat2 <-as.data.frame(dat2)

## -----------------------------------------------------------------------------
temporal <- list(c("goals", "feedback"),
                 c("balance"),
                 c("merging", "concentration","control","selfawareness","time", "autotelic")) 

# Define variable names
prior <- priorKnowledge(addtemporal = temporal)

## -----------------------------------------------------------------------------
tetradrunner <- tetradrunner(algoId = 'fges', #select causal algorithm
                             scoreId = 'sem-bic', #evaluation score used 
                             df = dat2, # dataset
                             dataType = 'continuous', # type of data
                             #faithfulnessAssumed = TRUE,
                             maxDegree = -1,
                             priorKnowledge = prior, # add prior knowledge 
                             numberResampling = 100, # bootstrapping
                             resamplingEnsemble = 2, # Majority
                             addOriginalDataset = TRUE,
                             verbose=TRUE)

tetradrunner$edges 
tetradrunner$nodes 

graph <- tetradrunner$graph
graph$getAttribute('BIC')

nodes <- graph$getNodes()
for(i in 0:as.integer(nodes$size()-1)){
  node <- nodes$get(i)
  cat(node$getName(),": ",node$getAttribute('BIC'),"\n")
}

graph_dot <- tetradrunner.tetradGraphToDot(tetradrunner$graph)
#dot(graph_dot)
#graph_dot

# reduce amount of decimal places in the initial output to 2
graph_dot2 <- strsplit(x = graph_dot,
                       split = "(?<=[^0-9.])(?=\\d)|(?<=\\d)(?=[^0-9.])", perl = T) %>%
    unlist %>% 
    lapply(function(x){if(!is.na(as.numeric(x)))x<-round(as.numeric(x),2)%>%         
    format(nsmall=2);x}) %>%
    paste(collapse="")

myGraph <- grViz(graph_dot2) 
myGraph

## -----------------------------------------------------------------------------
tetradrunner.getAlgorithmDescription(algoId = 'pc-all')

## -----------------------------------------------------------------------------
# define sufficient statistics
suffStat <- list(C = mycor, n = 236)

# specify tiers like we did above for FGES, based on the distinction between
# preconditions and core flow features
tiers <- c(2,3,1,1,3,3,3,3,3)

# we use code adapted from 
# https://github.com/erenaldis/GaussianStructureLearning/blob/master/application.R

V <- colnames(mycor)

alphas    <- seq(0.01, 0.1, 0.01)
scores    <- c()

detach("package:bnlearn", unload = TRUE)

for (alpha in alphas) {
  
  ## estimate CPDAG
  tpc.fit <- tpc(suffStat  = list(C = mycor, n = 236),
                 indepTest = gaussCItest, # partial correlations
                 alpha     = alpha,
                 labels    = colnames(mycor), 
                 tiers     = tiers)
  
  s_dag <- pdag2dag(tpc.fit@graph)$graph
  
  # BIC scores
  score <- new("GaussL0penObsScore", dat2)
  scores <- c(scores, score$global.score(as(s_dag, "GaussParDAG")))
  
  optimal.alpha <- alphas[which.max(scores)] # 0.04
}

## -----------------------------------------------------------------------------
data.frame (alphas = alphas,scores = scores) %>% 
ggplot(aes(alphas, scores)) +
  geom_line(color="red") +
  geom_point(size = 2, color = "red") +
  scale_x_continuous(limits=c(0.01,0.10), breaks=seq(0.01,0.1,0.01)) +
    geom_vline(xintercept = c(0, optimal.alpha), 
             linetype="longdash", color = "grey65", size=0.5) +
  labs(title ="BIC Score versus alpha") + 
  theme_minimal()

## -----------------------------------------------------------------------------
tpc.fit <- tpc(suffStat = list(C = mycor, n = 236),
              indepTest = gaussCItest,
                  alpha = optimal.alpha, 
                 labels = colnames(mycor),
                  tiers = tiers)

#show estimated CPDAG
  par(mfrow=c(1,1))
  plot(tpc.fit, main = "Estimated CPDAG")


## -----------------------------------------------------------------------------
adjmat <-as(tpc.fit, "amat")

## -----------------------------------------------------------------------------
isValidGraph(amat = adjmat, type = "dag",verbose = T) 

## -----------------------------------------------------------------------------
isValidGraph(amat = adjmat, type = "cpdag",verbose = T) 

## -----------------------------------------------------------------------------
isValidGraph(amat = adjmat, type = "pdag",verbose = T)

## -----------------------------------------------------------------------------
tpc.res <- pdag2allDags(adjmat, verbose = FALSE)
tpc.res

# Plot function (see pcald documentation)
plotAllCPs2 <- function(tpc.res) {
  require(graph)
  p <- sqrt(ncol(tpc.res$dags))
  nDags <- ceiling(sqrt(nrow(tpc.res$dags)))
  par(mfrow = c(nDags, nDags))
  for (i in 1:nrow(tpc.res$dags)) {
    
    amat <- matrix(tpc.res$dags[i,],p,p, byrow = TRUE)
    colnames(amat) <- rownames(amat) <- tpc.res$nodeNms
    plot(as(amat, "graphNEL"))
  }
}
par(mfrow = c(3,3))
plotAllCPs2(tpc.res)


## -----------------------------------------------------------------------------
# adapted from https://github.com/btaschler/scl_replicate/tree/b4ddeddf715d26254e75f4798e4a375539cee89b/R

run_ida <- function(X, sig_level, skel.method) {
  
  G_hat           <- matrix(nrow = ncol(X), ncol = ncol(X))
  dimnames(G_hat) <- list(colnames(X), colnames(X))
  
  suffStat <- list(C = cor(X), n = nrow(X))
  
  pc.fit <- tpc::tpc(suffStat, 
                     indepTest   = gaussCItest, 
                     alpha       = sig_level, 
                     p           = ncol(suffStat$C),
                     tiers       = tiers)
                      
  covX <- cov(X)
  
  for (j in 1:ncol(X)){
    
# Idafast estimates the multiset of possible total causal effects of one variable -x- on several target variables -y- from observational data
    eff <- pcalg::idaFast(j, 1:ncol(X), covX, pc.fit@graph)
    
# https://www.datasciencemadesimple.com/row-wise-minimum-row-minimum-in-dataframe-r-2/
    # rowMins returns minimum value of each row of the matrix.
    G_hat[j,] <- matrixStats::rowMins( abs(eff) )   
    }
  
  Output <- list( g_hat = G_hat )
  
  return(Output)
}

idamatrix <- run_ida(dat2, 0.05, "stable")

## -----------------------------------------------------------------------------
#weight matrix
idamatrix <- round(as.matrix(idamatrix[["g_hat"]]),2)

W <- idamatrix
W[ which(W == 0) ] <-  .Machine$double.xmin
W[ which(W == 0) ] <- 0.001


## -----------------------------------------------------------------------------
G <- bnlearn::as.bn(tpc.fit@graph)
G <- bnlearn::as.igraph(G)

G2 <-  graph_from_adjacency_matrix(as.matrix(W) * 
       as.matrix(get.adjacency(G, type="both")),
        mode="directed", weighted=TRUE)

E(G2)[ which(E(G2)$weight == .Machine$double.xmin) ]$weight <-  0.0

plot.igraph(G2, 
            edge.label=E(G2)$weight,
            edge.arrow.size = 0.5,
            layout=layout_in_circle, rescale=T)

## -----------------------------------------------------------------------------
# convert matrix W (idamatrix) to dataframe suitable for bn
# rows    -> from
# columns -> to

g <- graph.adjacency(W, weighted = TRUE)

# values in a data frame with 3 columns: from,to,value
df <- get.data.frame(g)
df$weight <- round(df$weight,3)

# bn/ cpdag
G <- bnlearn::as.bn(tpc.fit@graph)


## -----------------------------------------------------------------------------
library(bnlearn)
detach("package:tidySEM", unload=TRUE)


## -----------------------------------------------------------------------------
plot_s <- function(x,y){
  
  # convert arc to edge 
  mat <- matrix(0,length(nodes(x)),length(nodes(x)))
  rownames(mat) <- nodes(x)
  colnames(mat) <- nodes(x)
  
  for (i in 1:dim(arcs(x))[1]){
    mat[nodes(x)==arcs(x)[i,1],nodes(x)==arcs(x)[i,2]]<-1
  }
  
  # create the graphAM object from the bn object
  g1 <- graphAM(adjMat = mat, edgemode = "directed")
  
  to.score   <- apply(y[,1:2], 1, paste, collapse = "~")
  lab        <- as.character(y[,3])
  names(lab) <- to.score
  
  g1 <- layoutGraph(g1, edgeAttrs=list(label=lab))
  
  graph.par(list(nodes = list(lwd=1,
                              lty=1,
                              fontsize=22,
                              shape="ellipse"),
                 
                 edges = list(fontsize=10,
                              textCol="black")
                 )
            )
  
  nodeRenderInfo(g1) <- list(shape="ellipse")
  renderGraph(g1)
}

par(mfrow = c(1,1))
plot_s(G,df)


## -----------------------------------------------------------------------------
Bootstraps <- lapply(1:100, function(x) {
  
  bootstrapped_model <- tpc(suffStat  = list(
    C = cor(dat2[sample(1:nrow(dat2), nrow(dat2), TRUE), ]), 
    n = 236),
    
    indepTest = gaussCItest, 
    alpha     = optimal.alpha, 
    labels    = V,
    tiers     = tiers,
    verbose   = FALSE)
  
  return(bootstrapped_model)
})

# Are individual edges included (y/n)?
resBoots_temp <- lapply(Bootstraps, function(x)
                 ifelse(wgtMatrix(x) > 0, 1, 0))

# over all 100 bootstrap iterations, how often (%) is a particular edge included?
strength <- apply(simplify2array(resBoots_temp), 1:2, mean) 

# convert matrix to dataframe for plotting
str_g       <- graph.adjacency(strength, weighted=TRUE)
strength_df <- get.data.frame(str_g)
colnames(strength_df) <- c("from","to","direction probability")

par(mfrow = c(1,1))
plot_s(G,strength_df)


## -----------------------------------------------------------------------------
# PC as bn object
bn_pc <- as.bn(tpc.fit@graph)

#FGES as bn object
nodes <-  c("balance", "merging", "goals", "feedback", "concentration","control","selfawareness", "time", "autotelic") 

bn_fges <-empty.graph(nodes) 

arc.set <- matrix(c(
"autotelic","merging",
"balance","concentration",
"balance","control",
"goals","control", 
"balance","autotelic", 
"control","selfawareness",
"feedback", "goals", # undirected 
"goals", "feedback", # undirected 
"control","concentration", 
"autotelic","selfawareness",
"goals","balance", 
"control","merging",
"control","autotelic",
"feedback","balance"),
ncol = 2, byrow = TRUE,
dimnames = list(NULL, c("from", "to")))

arcs(bn_fges) <- arc.set

# plot two graphs
#par(mfrow = c(1,1))
#graphviz.plot(bn_pc)
#graphviz.plot(bn_fges)


## -----------------------------------------------------------------------------
# the true positive (tp) arcs appear in both structures
bnlearn::compare(bn_pc, bn_fges, arcs = TRUE)

# 4 edges differ between the CPDAGs of the two network
shd(bn_pc, bn_fges) 


## -----------------------------------------------------------------------------
par(mfrow = c(2,1))
graphviz.compare(bn_pc, bn_fges,
                 shape = "circle", layout = "dot",
                 main = c("PC", "FGES"),
                 sub = paste("SHD =", c("0", shd(bn_pc, bn_fges))),
                 diff.args = list(tp.lwd = 2, tp.col = "green"))
