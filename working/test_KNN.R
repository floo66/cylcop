
tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
(t$toc-t$tic)/4

tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,query = cbind(c(0.3,0.5,0.7),c(0.3,0.5,0.7)),k=5,bucketSize = 10,sort=F,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
(t$toc-t$tic)/4





cop <- normalCopula(0)
n <- 1000000
samp <- rcylcop(n,cop)
times_0 <- rep(0,10)
tic()
for(i in 1:4){
test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="STD")
}
t <- toc()
times_0[1] <- (t$toc-t$tic)/4

tic()
for(i in 1:4){
test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="MIDPT")
}
t <- toc()
times_0[2] <- (t$toc-t$tic)/4
for(i in 1:4){
tic()
test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="FAIR")
}
t <- toc()
times_0[3] <- (t$toc-t$tic)/4
for(i in 1:4){
tic()
test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
times_0[4] <- (t$toc-t$tic)/4
for(i in 1:4){
tic()
test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="SL_FAIR")
}
t <- toc()
times_0[5] <- (t$toc-t$tic)/4

for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="STD")
}
t <- toc()
times_0[6] <- (t$toc-t$tic)/4

tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="MIDPT")
}
t <- toc()
times_0[7] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="FAIR")
}
t <- toc()
times_0[8] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
times_0[9] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_FAIR")
}
t <- toc()
times_0[10] <- (t$toc-t$tic)/4






cop <- normalCopula(1)
n <- 1000000
samp <- rcylcop(n,cop)
times_1 <- rep(0,10)
tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="STD")
}
t <- toc()
times_1[1] <- (t$toc-t$tic)/4

tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="MIDPT")
}
t <- toc()
times_1[2] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="FAIR")
}
t <- toc()
times_1[3] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
times_1[4] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="SL_FAIR")
}
t <- toc()
times_1[5] <- (t$toc-t$tic)/4

for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="STD")
}
t <- toc()
times_1[6] <- (t$toc-t$tic)/4

tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="MIDPT")
}
t <- toc()
times_1[7] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="FAIR")
}
t <- toc()
times_1[8] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
times_1[9] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_FAIR")
}
t <- toc()
times_1[10] <- (t$toc-t$tic)/4



cop <- normalCopula(0.5)
n <- 1000000
samp <- rcylcop(n,cop)
times_05 <- rep(0,10)
tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="STD")
}
t <- toc()
times_05[1] <- (t$toc-t$tic)/4

tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="MIDPT")
}
t <- toc()
times_05[2] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="FAIR")
}
t <- toc()
times_05[3] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
times_05[4] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=T,search="kdtree",splitRule="SL_FAIR")
}
t <- toc()
times_05[5] <- (t$toc-t$tic)/4

for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="STD")
}
t <- toc()
times_05[6] <- (t$toc-t$tic)/4

tic()
for(i in 1:4){
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="MIDPT")
}
t <- toc()
times_05[7] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="FAIR")
}
t <- toc()
times_05[8] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_MIDPT")
}
t <- toc()
times_05[9] <- (t$toc-t$tic)/4
for(i in 1:4){
  tic()
  test <- dbscan::kNN(samp,5,sort=F,search="kdtree",splitRule="SL_FAIR")
}
t <- toc()
times_05[10] <- (t$toc-t$tic)/4
