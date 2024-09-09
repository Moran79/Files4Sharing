setwd("/Users/moran/Google_Drive/Course/Loh/Research/FoLDA/codeNplots")
library(folda)
library(ggplot2)
# library(gridExtra)

### NOTE ###
#> Please first execute all functions under the Helper section
#> before running those tests below. Thanks!

# 2.1 Premature Stop -------------------------------------------------------------------

K = 20
set.seed(443)
dat <- data.frame(X2 = rep(c(1,2),each = K),
                  X1 = c(runif(K, 0, 2), runif(K, 1, 4)),
                  Y = rep(c("A", "B"), each = K))
dat$Y[dat$X1 >= 2.5] <- "C"

(pFoLDA <- ggplot(data = dat)+
    geom_point(aes(x = X1, y = X2, fill = Y, shape = Y), size = 5, stroke = 0.3)+
    scale_y_continuous(limits = c(0.5, 2.5), expand=c(0,0))+
    scale_x_continuous(limits = c(-0.25, 4.25), expand=c(0,0))+
    scale_color_manual(values=hcl.colors(3))+
    scale_fill_manual(values=hcl.colors(3))+
    scale_shape_manual(values = c(21, 22, 24))+
    theme_bw(base_size = 8)+
    labs(fill = "Class", shape = "Class", x = expression(X[1]), y = expression(X[2])))

# ggsave("FoLDA.eps", plot = pFoLDA, width = 84, height = 30, device = "eps", units = "mm")

# 2.2 Partial Lambda Distribution -----------------------------------------

fakeData <- function(N, M, J){ # all noise variables
  response <- as.factor(sample(paste0("Y", seq_len(J)), size = N, replace = TRUE))
  m <- scale(matrix(rnorm(M*N), nrow = N))
  colnames(m) <- paste0("V", seq_len(M))
  return(list(response = response, m = m))
}

getStats <- function(...){
  # return the individual Pillai / Wilks
  datList <- fakeData(N = N, M = M, J = J)

  a <- folda:::forwardSel(m = datList$m, response = datList$response, ...)
  a$forwardInfo$statDiff
}

set.seed(443)
N <- 150; M <- 1; J = 3; p = 0; nSim <- 10000
ans1 <- replicate(nSim, getStats(testStat = "Wilks", correction = FALSE, alpha = 1))
# ans1 <- replicate(nSim, getStats(method = "Wilks", correction = FALSE, alpha = 1))
(pSimF1 <- ggplot() +
    geom_density(aes(x = (N - J - p) / (J - 1) * (1 - ans1) / ans1, color = "Simulated F"), linewidth = 1) +
    geom_density(aes(x = qf(seq(nSim) / (nSim+1), df1 = J - 1, df2 = N - J - p), color = "Theoretical F"), linewidth = 1)+
    scale_color_manual(name = "", values = c("Simulated F" = hcl.colors(3)[2],
                                             "Theoretical F" = hcl.colors(3)[1]))+
    theme_bw(base_size = 10)+
    # xlim(0,6)+
    labs(x = "F value", subtitle = "One Variable"))

N <- 150; M <- 2; J = 3; p = 1; nSim <- 10000
ans2 <- replicate(nSim, getStats(testStat = "Wilks", correction = FALSE, alpha = 1))
(pSimF2 <- ggplot() +
    geom_density(aes(x = (N - J - 0) / (J - 1) * (1 - ans2[1,]) / ans2[1,], color = "Simulated F - X1"), linewidth = 1) +
    geom_density(aes(x = (N - J - 1) / (J - 1) * (1 - ans2[2,]) / ans2[2,], color = "Simulated F - X2"), linewidth = 1) +
    geom_density(aes(x = qf(seq(nSim) / (nSim + 1), df1 = J - 1, df2 = N - J - 0), color = "Theoretical F - X1"), linewidth = 1) +
    geom_density(aes(x = qf(seq(nSim) / (nSim + 1), df1 = J - 1, df2 = N - J - 1), color = "Theoretical F - X2"), linewidth = 1) +
    scale_color_manual(name = "",
                       values = c("Simulated F - X1" = hcl.colors(3)[2],
                                  "Simulated F - X2" = hcl.colors(3)[3],
                                  "Theoretical F - X1" = hcl.colors(3)[1],
                                  "Theoretical F - X2" = hcl.colors(3)[1]),
                       labels = c("Simulated F - X1" = expression("Simulated F - " ~ X^{(1)}),
                                  "Simulated F - X2" = expression("Simulated F - " ~ X^{(2)}),
                                  "Theoretical F - X1" = expression("Theoretical F - " ~ X^{(1)}),
                                  "Theoretical F - X2" = expression("Theoretical F - " ~ X^{(2)}))) +
    theme_bw(base_size = 10) +
    # xlim(0, 6) +
    labs(x = "F value", subtitle = "Two Variables"))

# pList <- list(pSimF1, pSimF2)
# postscript("partialLambda.eps", width = 6.85, height = 4.5, paper = "special", onefile = FALSE, horizontal = FALSE)
# grid.arrange(grobs = pList, ncol = 1)
# dev.off()


# 2.3 Inflated Type I Error Rate ------------------------------------------

### NOTE ###
#> The threshold have been changed to 5 to show a larger effect.
#> In fact, the threshold of 4 aligns well with p = 0.05

getTest <- function(K, signal){
  if(signal){
    idx <- c(5,1,2,3,4)
  }else idx <- 5
  dat <- cbind.data.frame(response = iris[, idx], noise = matrix(rnorm(nrow(iris)*K), ncol = K))

  #> Fixed threshold of 4 has a smaller p-value than 0.05 when selecting the 5th variable,
  #> so it is safe to first select by 0.05 then subset.
  #> pf(4, df1 = 2, df2 = 150 - 3 - 4, lower.tail = FALSE)

  resForward <- folda:::forwardSel(m = scale(as.matrix(dat[, -1, drop = FALSE])), response = dat[, 1],
                                   testStat = "Wilks", correction = FALSE, alpha = 0.05)
  partialLambda <- resForward$forwardInfo$statDiff
  partialF <- (1 - partialLambda) / partialLambda * (150 - 3 - seq_along(partialLambda) + 1) / (3 - 1)
  idxKeep <- ifelse(any(partialF <= 5), which.max(partialF <= 5), length(partialLambda))

  errorP <- any(grepl("noise", resForward$forwardInfo$var))
  error4 <- any(grepl("noise", resForward$forwardInfo$var[seq_len(idxKeep)]))
  return(c(error4, errorP))
}

getTestHelper <- function(K, n = 2000, signal = FALSE){
  ans <- replicate(n = n, expr = getTest(K = K, signal = signal))
  data.frame(M = K, errorRate = apply(ans,1,mean), Threshold = c("fixed-4", "p-value"))
}

# 1.96 * sqrt(0.05*0.95/2000) # CI half length

set.seed(443)
datSignalNnoise <- do.call(rbind, lapply(2^(0:7), function(x) getTestHelper(K = x, signal = TRUE)))
datSignalNnoise$upper <- datSignalNnoise$errorRate + 1.96 * sqrt(datSignalNnoise$errorRate * (1 - datSignalNnoise$errorRate) / 2000)
datSignalNnoise$lower <- datSignalNnoise$errorRate - 1.96 * sqrt(datSignalNnoise$errorRate * (1 - datSignalNnoise$errorRate) / 2000)

# write.csv(datSignalNnoise, "datSignalNnoise.csv", row.names = FALSE)

(pSignalNnoise <- ggplot(data = datSignalNnoise)+
    geom_ribbon(aes(x = M, ymin = lower, ymax = upper, fill = Threshold), color = "black", linewidth = 0.05) +
    geom_point(aes(x = M, y = errorRate), color = "grey60", size = 2)+
    geom_hline(yintercept = 0.05, color = hcl.colors(3)[2], linetype = "dashed")+
    geom_label(aes(x = 100, y = 0.15), label = "alpha == 0.05", parse = TRUE, color = hcl.colors(3)[2])+
    scale_fill_manual(values = hcl.colors(2)) +
    scale_x_continuous(breaks = 2^(c(0,2:7)))+
    theme_bw(base_size = 10)+
    labs(x = "Number of Noise Features", y = "Type I Error Rate", subtitle = "Scenario 1: Signals and Noises"))


set.seed(443)
datNull <- do.call(rbind, lapply(2^(0:7), function(x) getTestHelper(K = x)))
datNull$upper <- datNull$errorRate + 1.96 * sqrt(datNull$errorRate * (1 - datNull$errorRate) / 2000)
datNull$lower <- datNull$errorRate - 1.96 * sqrt(datNull$errorRate * (1 - datNull$errorRate) / 2000)

# write.csv(datNull, "datNull.csv", row.names = FALSE)

(pNull <- ggplot(data = datNull)+
    geom_ribbon(aes(x = M, ymin = lower, ymax = upper, fill = Threshold), color = "black", linewidth = 0.05) +
    geom_point(aes(x = M, y = errorRate), color = "grey60", size = 2)+
    geom_hline(yintercept = 0.05, color = hcl.colors(3)[2], linetype = "dashed")+
    geom_label(aes(x = 100, y = 0.15), label = "alpha == 0.05", parse = TRUE, color = hcl.colors(3)[2])+
    scale_fill_manual(values = hcl.colors(2)) +
    scale_x_continuous(breaks = 2^(c(0,2:7)))+
    theme_bw(base_size = 10)+
    labs(x = "Number of Noise Features", y = "Type I Error Rate", subtitle = "Scenario 2: Pure Noises"))

# pList <- list(pSignalNnoise, pNull)
# postscript("typeIwilks.eps", width = 6.85, height = 4.5, paper = "special", onefile = FALSE, horizontal = FALSE)
# grid.arrange(grobs = pList, ncol = 1)
# dev.off()


# 3.1 Improvement in Speed ------------------------------------------------

library(parallel)
num_cores <- 10

benchMarkInner <- function(M) {
  dat <- fakeData(N = 10000, M = M, J = 10)
  time1 <- proc.time()
  Original = folda(datX = as.data.frame(dat$m), response = dat$response, subsetMethod = c("all"), algorithm = "Original")
  time2 <- proc.time()
  QR = folda(datX = as.data.frame(dat$m), response = dat$response, subsetMethod = c("all"), algorithm = "QR")
  time3 <- proc.time()
  Cholesky = folda(datX = as.data.frame(dat$m), response = dat$response, subsetMethod = c("all"), algorithm = "Cho")
  time4 <- proc.time()
  return(as.numeric(c((time2 - time1)["user.self"],
                      (time3 - time2)["user.self"],
                      (time4 - time3)["user.self"])))
}

benchMarkOuter <- function(M){
  nRep <- 30
  results <- do.call(rbind, mclapply(1:nRep, function(i) benchMarkInner(M = M), mc.cores = num_cores - 1))
  data.frame(Algorithm = c("Original", "QR", "Cholesky"),
             time = apply(results,2,mean),
             timeSd = apply(results,2,sd) / sqrt(nRep),
             M = M)
}

set.seed(443)
datQR <- do.call(rbind, mclapply(2^c(1:10), function(i) benchMarkOuter(M = i), mc.cores = num_cores - 1))
datQR$upper <- datQR$time + 1.96 * datQR$timeSd
datQR$lower <- datQR$time - 1.96 * datQR$timeSd
# write.csv(datQR, "datQR.csv", row.names = FALSE)

(pQR <- ggplot(data = datQR)+
    geom_ribbon(aes(x = M, ymin = lower, ymax = upper, fill = Algorithm), color = "black", linewidth = 0.05) +
    geom_point(aes(x = M, y = time), color = "grey60", size = 1)+
    scale_color_manual(values = hcl.colors(3)) +
    scale_fill_manual(values = hcl.colors(3)) +
    scale_x_continuous(breaks = 2^(c(0,6:10)))+
    theme_bw(base_size = 10)+
    labs(x = "Number of Features", y = "Runtime (secs)"))

# ggsave("QR.eps", plot = pQR, width = 174, height = 50, device = "eps", units = "mm")


# new 4.1 -----------------------------------------------------------------


getTestNew <- function(K, signal){
  if(signal){
    idx <- c(5,1,2,3,4)
  }else idx <- 5
  dat <- cbind.data.frame(response = iris[, idx], noise = matrix(rnorm(nrow(iris)*K), ncol = K))

  errorPillai <- any(grepl("noise", folda:::forwardSel(m = scale(as.matrix(dat[,-1, drop = FALSE])), response = dat[,1], testStat = "Pillai", correction = TRUE, alpha = 0.05)$forwardInfo$var))
  errorWilksF <- any(grepl("noise", folda:::forwardSel(m = scale(as.matrix(dat[,-1, drop = FALSE])), response = dat[,1], testStat = "Wilks", correction = FALSE, alpha = 0.05)$forwardInfo$var))
  errorWilksT <- any(grepl("noise", folda:::forwardSel(m = scale(as.matrix(dat[,-1, drop = FALSE])), response = dat[,1], testStat = "Pillai", correction = TRUE, alpha = 0.05)$forwardInfo$var))
  return(c(errorPillai, errorWilksF, errorWilksT))
}

getTestHelperNew <- function(K, n = 2000, signal = FALSE){
  ans <- replicate(n = n, expr = getTestNew(K = K, signal = signal))
  data.frame(M = K, errorRate = apply(ans,1,mean), Threshold = c("Pillai", "Wilks", "Wilks.Bonferroni"))
}


set.seed(443)
datSignalNnoise3 <- do.call(rbind, mclapply(2^(0:7), function(x) getTestHelperNew(K = x, signal = TRUE), mc.cores = num_cores - 1))
datSignalNnoise3$upper <- datSignalNnoise3$errorRate + 1.96 * sqrt(datSignalNnoise3$errorRate * (1 - datSignalNnoise3$errorRate) / 2000)
datSignalNnoise3$lower <- datSignalNnoise3$errorRate - 1.96 * sqrt(datSignalNnoise3$errorRate * (1 - datSignalNnoise3$errorRate) / 2000)

# write.csv(datSignalNnoise3, "datSignalNnoise3.csv", row.names = FALSE)

(pSignalNnoise3 <- ggplot(data = datSignalNnoise3)+
    geom_ribbon(aes(x = M, ymin = lower, ymax = upper, fill = Threshold), color = "black", linewidth = 0.05) +
    geom_point(aes(x = M, y = errorRate), color = "grey60", size = 2)+
    geom_hline(yintercept = 0.05, color = "#C1121F", linetype = "dashed")+
    geom_label(aes(x = 100, y = 0.15), label = "alpha == 0.05", parse = TRUE, color = "#C1121F")+
    scale_color_manual(values = hcl.colors(3)) +
    scale_fill_manual(values = hcl.colors(3)) +
    scale_x_continuous(breaks = 2^(c(0,2:7)))+
    theme_bw(base_size = 10)+
    labs(x = "Number of Noise Features", y = "Type I Error Rate", subtitle = "Scenario 1: Signals and Noises"))


set.seed(443)
datNull3 <- do.call(rbind, mclapply(2^(0:7), function(x) getTestHelperNew(K = x, signal = FALSE), mc.cores = num_cores - 1))
datNull3$upper <- datNull3$errorRate + 1.96 * sqrt(datNull3$errorRate * (1 - datNull3$errorRate) / 2000)
datNull3$lower <- datNull3$errorRate - 1.96 * sqrt(datNull3$errorRate * (1 - datNull3$errorRate) / 2000)

# write.csv(datNull3, "datNull3.csv", row.names = FALSE)

(pNull3 <- ggplot(data = datNull3)+
    geom_ribbon(aes(x = M, ymin = lower, ymax = upper, fill = Threshold), color = "black", linewidth = 0.05) +
    geom_point(aes(x = M, y = errorRate), color = "grey60", size = 2)+
    geom_hline(yintercept = 0.05, color = "#C1121F", linetype = "dashed")+
    geom_label(aes(x = 100, y = 0.15), label = "alpha == 0.05", parse = TRUE, color = "#C1121F")+
    scale_color_manual(values = hcl.colors(3)) +
    scale_fill_manual(values = hcl.colors(3)) +
    scale_x_continuous(breaks = 2^(c(0,2:7)))+
    theme_bw(base_size = 10)+
    labs(x = "Number of Noise Features", y = "Type I Error Rate", subtitle = "Scenario 2: Pure Noises"))

# pList <- list(pSignalNnoise3, pNull3)
# postscript("typeIwilks3.eps", width = 6.85, height = 4.5, paper = "special", onefile = FALSE, horizontal = FALSE)
# grid.arrange(grobs = pList, ncol = 1)
# dev.off()


# 4.2 ---------------------------------------------------------------------

# Scenario 1
set.seed(443)
dat <- data.frame(Class = as.factor(sample(LETTERS[1:10], 2000, replace = TRUE)))
dat <- cbind(dat, model.matrix(~.-1, data = dat))

(fitWilks <- folda(datX = dat[,-1], response = dat[,1], testStat = "Wilks", alpha = 0.05))
(fitPillai <- folda(datX = dat[,-1], response = dat[,1], testStat = "Pillai", alpha = 0.05))

# Scenario 2
dat <- read.csv("NHTSA_clean.csv")
(fitWilks <- folda(datX = dat[,-1], response = dat[,1], testStat = "Wilks", alpha = 0.05))
(fitPillai <- folda(datX = dat[,-1], response = dat[,1], testStat = "Pillai", alpha = 0.05))

# dat$ENGINE[which(dat$MODELD %in% c("RX", "RX-8"))]
# table(dat$ENGINE) # ROTR

CVtestingHelper <- function(i, idx){
  # Data Preprocessing
  dat_train <- droplevels(dat[idx != i,,drop = FALSE]) # Unused level
  idxKeep <- folda:::nonConstInd(dat_train)
  dat_train <- dat_train[, idxKeep, drop = FALSE]
  dat_test <- dat[idx == i,idxKeep, drop = FALSE]

  fitWilks <- folda(datX = dat_train[,-1], response = dat_train[,1], testStat = "Wilks", alpha = 0.05)
  fitPillai <- folda(datX = dat_train[,-1], response = dat_train[,1], testStat = "Pillai", alpha = 0.05 / 1e11)
  errorWilks <- sum(dat_test$ENGINE != predict(fitWilks, dat_test))
  errorPillai <- sum(dat_test$ENGINE != predict(fitPillai, dat_test))
  return(c(errorWilks, errorPillai))
}

CVtesting <- function(){
  idx <- sample(seq_len(kFold), size = nrow(dat), replace = TRUE)
  results <- mclapply(seq_len(kFold), function(i) CVtestingHelper(i = i, idx = idx), mc.cores = num_cores - 1)
  1 - apply(do.call(rbind, results),2,sum) / nrow(dat) # Testing accuracy
}

set.seed(443)
kFold <- 9
CVtesting() # 0.3819126 0.6495570


# Helper Functions --------------------------------------------------------

library(Rcpp)
library(RcppEigen)

cppFunction('
Rcpp::List choEigen(const Eigen::MatrixXd &A) {
  Eigen::LLT<Eigen::MatrixXd> llt(A);
  Eigen::MatrixXd L = llt.matrixL();
  Eigen::MatrixXd Lt = L.transpose();
  return Rcpp::List::create(Rcpp::Named("Lt") = Lt);
}
', depends = "RcppEigen")

folda <- function(datX,
                  response,
                  subsetMethod = c("forward", "all"),
                  testStat = c("Pillai", "Wilks"),
                  correction = TRUE,
                  alpha = 0.1,
                  prior = NULL,
                  misClassCost = NULL,
                  missingMethod = c("medianFlag", "newLevel"),
                  downSampling = FALSE,
                  kSample = NULL,
                  algorithm = "Original"){
  
  # Pre-processing: Arguments ----------------------------------
  
  if (!is.data.frame(datX)) stop("datX must be a data.frame")
  response <- droplevels(as.factor(response))
  subsetMethod <- match.arg(subsetMethod, c("forward", "all"))
  
  # Pre-processing: Data Cleaning -----------------------------------------------
  
  idxTrain <- folda:::getDownSampleInd(response = response,
                                       downSampling = downSampling,
                                       kSample = kSample)
  response <- droplevels(response[idxTrain])
  priorAndMisClassCost <- folda:::checkPriorAndMisClassCost(prior = prior, misClassCost = misClassCost, response = response)
  
  imputedSummary <- missingFix(data = datX[idxTrain, , drop = FALSE], missingMethod = missingMethod)
  datX <- imputedSummary$data # this step also removes some constant columns
  if(any(dim(datX) == 0)) stop("No available observations or features, which maybe due to preprocessing steps.")
  
  # Pre-processing: Data Transformation -----------------------------------------------
  
  modelFrame <- stats::model.frame(formula = ~.-1, datX, na.action = "na.fail")
  Terms <- stats::terms(modelFrame)
  m <- scale(stats::model.matrix(Terms, modelFrame)) # constant cols would be changed to NaN in this step
  currentVarList <- seq_len(ncol(m))
  
  # Forward Selection -----------------------------------------------
  
  if(subsetMethod == "forward"){
    forwardRes <- folda:::forwardSel(m = m,
                                     response = response,
                                     testStat = testStat,
                                     alpha = alpha,
                                     correction = correction)
    
    # When no variable is selected, use the full model
    if(length(forwardRes$currentVarList) != 0){
      # modify the design matrix to make it more compact
      selectedVarRawIdx <- unique(sort(attributes(m)$assign[forwardRes$currentVarList]))
      modelFrame <- stats::model.frame(formula = ~.-1, datX[, selectedVarRawIdx, drop = FALSE], na.action = "na.fail")
      Terms <- stats::terms(modelFrame)
      m <- scale(stats::model.matrix(Terms, modelFrame))
      currentVarList <- which(colnames(m) %in% forwardRes$forwardInfo$var)
    }
  }
  
  varSD <- attr(m,"scaled:scale")[currentVarList]
  varCenter <- attr(m,"scaled:center")[currentVarList]
  m <- m[, currentVarList, drop = FALSE]
  
  # ULDA -----------------------------------------------
  
  # Step 1: SVD on the combined matrix H
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  Hb <- sqrt(tabulate(response)) * groupMeans
  Hw <- m - groupMeans[response, , drop = FALSE]
  if(algorithm != "Original" & diff(dim(m)) < 0){ # More rows than columns
    if(algorithm == "QR"){
      qrRes <- folda:::qrEigen(Hw)
      fitSVD <- folda:::svdEigen(rbind(Hb, qrRes$R))
    }else if(algorithm == "Cho"){
      Sw <- t(Hw) %*% Hw
      cho <- choEigen(Sw)
      fitSVD <- folda:::svdEigen(rbind(Hb, cho$Lt))
    }
  } else fitSVD <- folda:::svdEigen(rbind(Hb, Hw))
  
  # Step 2: SVD on the P matrix
  N <- nrow(m); J <- nlevels(response)
  rankT <- sum(fitSVD$d >= max(dim(fitSVD$u), dim(fitSVD$v)) * .Machine$double.eps * fitSVD$d[1])
  fitSVDp <- folda:::svdEigen(fitSVD$u[seq_len(J), seq_len(rankT), drop = FALSE], uFlag = FALSE)
  rankAll <- min(J - 1, sum(fitSVDp$d >= max(J, rankT) * .Machine$double.eps * fitSVDp$d[1]))
  
  # Step 3: Transform Sw into identity matrix
  unitSD <- diag(sqrt((N - J) / abs(1 - fitSVDp$d^2 + 1e-5)), nrow = rankAll) # Scale to unit var
  scalingFinal <- (fitSVD$v[, seq_len(rankT), drop = FALSE] %*% diag(1 / fitSVD$d[seq_len(rankT)], nrow = rankT) %*% fitSVDp$v[, seq_len(rankAll), drop = FALSE]) %*% unitSD
  rownames(scalingFinal) <- colnames(m)
  groupMeans <- groupMeans %*% scalingFinal
  rownames(groupMeans) <- levels(response)
  colnames(groupMeans) <- colnames(scalingFinal) <- paste("LD", seq_len(ncol(groupMeans)), sep = "")
  
  # Summary and outputs -----------------------------------------------
  
  statPillai <- sum(fitSVDp$d[seq_len(rankAll)]^2)
  p <- rankT; s <- rankAll; numF <- N - J - p + s; denF <- abs(p - J + 1) + s
  pValue <- ifelse(numF > 0, stats::pbeta(1 - statPillai / s, shape1 = numF * s / 2, shape2 = denF * s / 2), 0)
  
  res <- list(scaling = scalingFinal, groupMeans = groupMeans, prior = priorAndMisClassCost$prior,
              misClassCost = priorAndMisClassCost$misClassCost, misReference = imputedSummary$ref,
              terms = Terms, xlevels = stats::.getXlevels(Terms, modelFrame), varIdx = currentVarList,
              varSD = varSD, varCenter = varCenter, statPillai = statPillai, pValue = pValue)
  
  if(subsetMethod == "forward"){
    res$forwardInfo = forwardRes$forwardInfo
    res$stopInfo <- forwardRes$stopInfo
  }
  
  class(res) <- "ULDA"
  pred <- factor(stats::predict(res, datX), levels = levels(response))
  res$predGini <- 1 - sum(unname(table(pred) / dim(datX)[1])^2)
  res$confusionMatrix <- table(Predicted = pred, Actual = response)
  return(res)
}


