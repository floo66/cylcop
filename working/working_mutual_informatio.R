library(copula)
library(cylcop)
library(infotheo)

cop <- normalCopula(0.1)
n <- 10
samp <- rcylcop(n,cop)

samp <- rjoint(100,cop,
               marginal_1 = list(name="norm",coef=list(2,3)),
               marginal_2 = list(name="exp",coef=list(1,5)))
emp_cop <- pobs(samp, ties.method = "average")%>%as.data.frame()
infotheo::mutinformation(infotheo::discretize(samp,nbins = 5,disc = "equalwidth"))
infotheo::mutinformation(infotheo::discretize(emp_cop,nbins = 5,disc = "equalwidth"))
infotheo::mutinformation(infotheo::discretize(samp,nbins = 5,disc = "equalfreq"))
infotheo::mutinformation(infotheo::discretize(emp_cop,nbins = 5,disc = "equalfreq"))
gptentropy(samp[,1],5)




nbins <-round(max(2,nrow(emp_cop)^(1/3)))

discr_theta <- infotheo::discretize(emp_cop[,1],nbins=nbins)
discr_x <- infotheo::discretize(emp_cop[,2],nbins=nbins)
infotheo::mutinformation(discr_theta,discr_x,method="mm") / sqrt(infotheo::entropy(discr_theta) * infotheo::entropy(discr_x))
infotheo::mutinformation(discr_theta,discr_x) / min(infotheo::entropy(discr_theta),infotheo::entropy(discr_x))
infotheo::mutinformation(discr_theta,discr_x) / max(infotheo::entropy(discr_theta),infotheo::entropy(discr_x))
2*infotheo::mutinformation(discr_theta,discr_x) / (infotheo::entropy(discr_theta)+infotheo::entropy(discr_x))



#####################direct copula entropy

log_cop_samp <- log(dcylcop(samp,copula =cop),exp(1))
entrop <- -sum(log_cop_samp)/n
entrop_var <- sum((log_cop_samp-entrop)^2)/(n+1)
sqrt(entrop_var)
entrop


entknn(samp,k=3,dt=2)


####################binning
nbins <-round(max(2,nrow(samp)^(1/3)))
nbins <-ceiling(1+log(n,base=2))
discr <- infotheo::discretize(emp_cop,nbins=nbins,disc = "equalwidth")
infotheo::mutinformation(discr)
infotheo::mutinformation(discr,method="sg")
infotheo::mutinformation(discr,method="shrink")
infotheo::mutinformation(discr,method="mm")

- (infotheo::entropy(discr,method="emp") + log((1/nbins)^2))
- (infotheo::entropy(discr,method="sg") + log((1/nbins)^2))
- (infotheo::entropy(discr,method="shrink")+ log((1/nbins)^2))
- (infotheo::entropy(discr,method="mm")+ log((1/nbins)^2))


##################k-nearest neighbors, linear
copent<-function(cop_samp,k=3,dt=2){
  xc = cop_samp
  ce1 = -entknn(xc,k,dt)
  if(is.infinite(ce1)){ # log0
    N = dim(x)[1]; d = dim(x)[2];
    for(i in 1:d){
      max1 = max(abs(x[,i]))
      if(max1 > 0){x[,i] = x[,i] + max1 * 0.00000001 * runif(N)}
      else{x[,i] = runif(N)}
    }
    xc = cop_samp
    ce1 = -entknn(xc,k,dt)
  }
  ce1
}

entknn<-function(x,k=3,dt=2){
  x = as.matrix(x)
  N = dim(x)[1];  d = dim(x)[2];
  g1 = digamma(N) - digamma(k);
  if (dt == 1){	# euciledean distance
    cd = pi^(d/2) / 2^d / gamma(1+d/2);
    distx = as.matrix(dist(x));
  }
  else {	# maximum distance
    cd = 1;
    distx = as.matrix(dist(x,method = "maximum"));
  }
  logd = 0;
  for(i in 1:N){
    distx[i,] = sort(distx[i,]);
    logd = logd + log( 2 * distx[i,k+1] ) * d / N;
  }
  g1 + log(cd) + logd
}








#find nearest neighbor fast
fast_knn <- function(emp_cop,k=3,periodic=T,metric="max",buffer_width=0.5){
  N <- nrow(emp_cop)
  ranks <- emp_cop*(N+1)
  ranks <- round(ranks[order(ranks[,1]),],0)
  if(periodic){
    buffer_width <- ceiling(buffer_width*N)
    ranks <- rbind(tail(ranks,buffer_width),ranks,head(ranks,buffer_width))
    ranks[,1] <- ranks[,1]+c(rep(-N,buffer_width),rep(0,N),rep(N,buffer_width))+buffer_width
    points_2_eval <- seq_len(N) + buffer_width
  }else{
    points_2_eval <- seq_len(N)
  }

  if(metric=="max"){
    cd <- 1
    rotate_ranks_45 <- rotate(ranks,45)
    neighbors <-  RANN1::nn2(rotate_ranks_45,
                             rotate_ranks_45[points_2_eval,],
                             k=(k+1))
    dist <- round(neighbors$nn.dists[, (k+1)] / sqrt(1^2 + 1^2), 0)
  }else if(metric=="euclid"){
    cd <- pi/4
    neighbors <-  dbscan::kNN(x=ranks,
                              query = ranks[points_2_eval,],
                              k=k+1,
                              bucketSize = 10,
                              sort=F,
                              search="kdtree",
                              splitRule="SL_MIDPT")
    dist <- neighbors$dist[, (k+1)]
  }

  if(periodic){
    u_dist_nearest_border <- pmin(seq_len(N) + buffer_width - 1, N + buffer_width - seq_len(N))
    outside_buffer <- u_dist_nearest_border < dist
    if(any(outside_buffer)){
      ind_outside_buff <- which(outside_buffer)
      new_width <- ceiling(max(dist[ind_outside_buff]))
      new_ranks <- round(ranks[order(ranks[,1]),],0)
      new_ranks <- rbind(tail(new_ranks,new_width),new_ranks,head(new_ranks,new_width))
      new_ranks[,1] <- new_ranks[,1]+c(rep(-N,new_width),rep(0,N),rep(N,new_width))+new_width
      new_points_2_eval <- (seq_len(N) + new_width)[ind_outside_buff]


      if(metric=="max"){
        new_rotate_ranks_45 <- as.data.frame(rotate(new_ranks,45))
        new_neighbors <-  RANN1::nn2(new_rotate_ranks_45,
                                     new_rotate_ranks_45[new_points_2_eval,],
                                     k=(k+1))
        dist_new <- round(new_neighbors$nn.dists[,(k+1)]/sqrt(1^2+1^2),0)
      }else if(metric=="euclid"){
        new_neighbors <-  dbscan::kNN(x=new_ranks,
                                      query = new_ranks[new_points_2_eval,],
                                      k=k+1,
                                      bucketSize = 10,
                                      sort=F,
                                      search="kdtree",
                                      splitRule="SL_MIDPT")
        dist_new <- new_neighbors$dist[,(k+1)]
      }
      dist[ind_outside_buff] <- dist_new
    }
  }

  digamma(N) - digamma(k)+log(cd) +2/N*sum(log(dist)) - 2*log(N+1)
}




rotate <- function(x, angle) {
  angle <- pi * angle / 180
  r <- sqrt(x[, 1]^2 + x[, 2]^2)
  theta <- atan2(x[, 2], x[, 1])
  cbind(r * cos(theta - angle), r * sin(theta - angle))
}



floentropy <- function(x, nbins=round(nrow(x)^(1/3),0), normalize = FALSE) {
  bins <- seq(0, 1, length.out = nbins + 1)
  freq <- as.data.frame(table(cut(x[, 1], breaks = bins), cut(x[, 2], breaks = bins)))
  freq[, 1:2] <- lapply(freq[, 1:2], as.numeric)

  freq2D <- matrix(0, nrow = nbins, ncol = nbins)
  freq2D[cbind(freq[, 1], freq[, 2])] <- freq[, 3]

  p <- freq2D / nrow(x)
  ent <- p * log(p)
  ent[is.na(ent)] <- 0

  # -(-sum(ent)  + log(1 / (nbins^2)))

  mi <- -(-sum(ent)  + log(1 / ((nbins)^2))+((nbins)^2-1-15)/(2*nrow(x)))

  if (normalize) {
    mi / log(nbins)
  } else {
    mi
  }
}

