install.packages('psych')
library(psych)
library(ggplot2)

# Getting the data ready

trackdata = read.csv('mvs.csv')
#trackdata <- trackdata[1:54,]

View(trackdata)
summary(trackdata)
datamat<- as.matrix(trackdata[,2:9])

n <- nrow(datamat); n
p <- ncol(datamat); p
summary(datamat)
dataMean <- colMeans(datamat); dataMean
sqvar <- sqrt(diag(var(datamat)));sqvar
datavar <- var(datamat); datavar
datacor <- cor(datamat); datacor
corPlot(datacor)
cent_data <- datamat-t(matrix(rep(dataMean,n),p,n))
Stddata <- cent_data/t(matrix(rep(sqvar,n),p,n))
boxplot(datamat)

uni_mat <- datamat
uni_mat[,1] <- uni_mat[,1]/60
uni_mat[,2] <- uni_mat[,2]/60
uni_mat[,3] <- uni_mat[,3]/60
colnames(uni_mat) <- c('X100.m.min','X200.m.min','X400.m.min','X800.m.min','X1500.m.min','X5000.m.min','X10000.m.min','Marathon.min.')
uni_mat
boxplot(uni_mat, main = "Uniform Unit")
uni_mat_cor <- cor(uni_mat)
corPlot(uni_mat_cor, main = "Uniform unit") 

# =======================================================
# Factor Analysis by PCA                                
# =======================================================

#=========================================================================
# Principal Component Analysis by "prcomp" method on S (covariance matrix)
#=========================================================================
# pca
dataPCAS <- prcomp(datamat,center = FALSE,scale. = FALSE); dataPCAS[[2]][,1:2]
summary(dataPCAS)
plot(dataPCAS,type="lines", main="Screeplot")
###biplot(dataPCAS)

# factor loading
dataFactLoadS <- dataPCAS[[2]][,1:2]*t(matrix(rep(dataPCAS[[1]][1:2],p),2,p))  #[\sqrt(\lambda_1)e_1:...:\sqrt(\lambda_p)e_p]
print(dataFactLoadS)  #Loadings on factors, i.e., full L matrix

# rotated factor loading
varimax(dataFactLoadS)$loadings
# communalities
commS <- cbind(dataFactLoadS,varimax(dataFactLoadS)$loadings,'communalities'=dataFactLoadS[,1]^2 + dataFactLoadS[,2]^2)
commS
# specific vairance
specvarS <- 1 - commS[,5]; specvarS
tableS <- cbind(commS,specvarS); tableS
colnames(tableS) <- c('F1','F2','rotated F1','rotated F2','communalities','specific variance')
tableS
# factor scores
fact_scoreS <- t(dataFactLoadS)%*%solve(datavar)%*%t(cent_data)
print(fact_scoreS)
plot(fact_scoreS[1,],fact_scoreS[2,],xlab = 'Factor 1',ylab='Factor 2',main='Factor scores')

# outliers
outliers <- trackdata$Country[fact_scoreS[1,]>40 & fact_scoreS[2,] < 0]
outliers



#==========================================================================
# Principal Component Analysis by "prcomp" method on R (correlation matrix)
#==========================================================================

# pca
dataPCAR <- prcomp(datamat,center = TRUE,scale. = TRUE); dataPCAR[[2]][,1:2]
summary(dataPCAR)
plot(dataPCAR,type="lines", main="Screeplot")
###biplot(StddataspeedPCA)

# factor loading
dataFactLoadR <- dataPCAR[[2]][,1:2]*t(matrix(rep(dataPCAR[[1]][1:2],p),2,p))  #[\sqrt(\lambda_1)e_1:...:\sqrt(\lambda_p)e_p]
print(dataFactLoadR)  #Loadings on factors,i.e., full L matrix
# rotated factor loading
varimax(dataFactLoadR)$loadings
# communalities
commR <- cbind(dataFactLoadR,varimax(dataFactLoadR)$loadings,'communalities'=dataFactLoadR[,1]^2 + dataFactLoadR[,2]^2)
commR

# specific vairance
specvarR <- 1 - commR[,5]; specvarR
tableR <- cbind(commR,specvarR); tableR
colnames(tableR) <- c('F1','F2','rotated F1','rotated F2','communalities','specific variance')
tableR


fac <- fa(datacor,nfactors = 2,rotate = 'varimax',fm='mle')
fa.diagram(fac)

# Factor Score by Regression L_z'R^{-1}z'
fact_scoreR <- t(dataFactLoadR)%*%solve(datacor)%*%t(Stddata)
fact_scoreR
plot(fact_scoreR[1,],fact_scoreR[2,] , xlab = 'Factor 1' , ylab='Factor 2',main='Factor scores')

# outliers
outliers <- trackdata$Country[fact_scoreR[1,] < (-15.3)]
outliers


# residual matrix on 1 factor
Res1=datacor-dataFactLoadR[,1]%*%t(dataFactLoadR[,1])-diag(dataSpecVarR[,1])
print(Res1) #R - LL' -\Psi, when one factor is considered
sum(Res1*Res1)

# residual matrix on 2 factors
Res2=datacor-dataFactLoadR[,c(1,2)]%*%t(dataFactLoadR[,c(1,2)])-diag(dataSpecVarR[,2])
print(Res2) #R - LL' -\Psi, when two factors are considered
sum(Res2*Res2)





#===================#
#       MLE         #
#===================#

factmle <- factanal(datamat,factors = 2,scores = "none")
print(factmle)
mlespvr <- 1-t(apply(factmle$loadings^2,1,cumsum))
print(diag(mlespvr[,2]))
# R-LL'-\shi
Res2M<-datacor-factmle$loadings[,c(1,2)]%*%t(factmle$loadings[,c(1,2)])-diag(mlespvr[,2])
print(Res2M)
sum(Res2M*Res2M)


