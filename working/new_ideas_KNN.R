gptEqualFreq2 <- function(x, n) {
  nx <- length(x)
  nrepl <- floor(nx / n)
  nplus <- sample(1:n, nx - nrepl * n)
  nrep <- rep(nrepl, n)
  nrep[nplus] <- nrepl + 1
  x[order(x)] <- rep(seq.int(n), nrep)
  x
}

gptentropy <- function(samp, nbins) {
  n <- length(samp)
  test2 <- gptEqualFreq2(samp, nbins)
  temp <- split(samp, test2)

  mins <- sapply(temp, min)
  maxs <- sapply(temp, max)

  limits <- c(mins[1], (maxs[-length(maxs)] + mins[-1]) / 2, tail(maxs, 1))
  widths <- diff(limits)

  p <- table(test2) / n
  -sum(p * log(p / widths))
}






x1 <- runif(100, 0, 2*pi)
x2 <- runif(100, 0,3)
DATA <- as.data.frame(samp)
nearest <- nn2(DATA,DATA,k=3)


emp_cop <- pobs(samp, ties.method = "average")%>%as.data.frame()
plot(emp_cop*101)
ranks <- emp_cop*101
ranks_u <- round(ranks[order(ranks[,1]),],0)
ranks_v <- round(ranks[order(ranks[,2]),],0)
N <- 100
log_eps <- 0

#find nearest neighbor slow
tic()
for(u_ind in 1:N){
  if(u_ind%%1000==0)cat(u_ind)
d <- 0
n_neighbor <- 0
v_ind <- ranks_u[u_ind,2]
while(n_neighbor<k){
  # cat("\n",d,"    ",n_neighbor)
  d <- d+1
  if(u_ind-d>0){
  if(abs(v_ind-ranks_u[u_ind-d,2])<=d){
    n_neighbor <- n_neighbor+1
  }
  }
  if(u_ind+d<=N){
  if(abs(v_ind-ranks_u[u_ind+d,2])<=d){
    n_neighbor <- n_neighbor+1
  }
  }
  if(v_ind+d<=N){
  if(abs(u_ind-ranks_v[v_ind+d,1])<d){
    n_neighbor <- n_neighbor+1
  }
  }
  if(v_ind-d>0){
  if(abs(u_ind-ranks_v[v_ind-d,1])<d){
    n_neighbor <- n_neighbor+1
  }
  }
}
log_eps <- log_eps + log(2*d)
}
toc()


#find nearest neighbor copent
tic()
log_eps <- 0
distx = as.matrix(dist(as.matrix(ranks_u),method = "maximum"))
for(i in 1:N){
  d <- sort(as.vector(distx[i,]))[k+1]
  log_eps <- log_eps + log(2*d)
}
toc()


#find nearest neighbor fast
fast_knn <- function(emp_cop,periodic=T,metric="max",buffer_width=0.5){
N <- nrow(emp_cop)
ranks <- emp_cop*(N+1)
if(!periodic){
  ranks_u <- round(ranks[order(ranks[,1]),],0)
  points_2_eval <- seq(1,N)
}else{
  buffer_width <- ceiling(buffer_width*N)
  ranks_u <- round(ranks[order(ranks[,1]),],0)
  ranks_u <- rbind(tail(ranks_u,buffer_width),ranks_u,head(ranks_u,buffer_width))
  ranks_u[,1] <- ranks_u[,1]+c(rep(-N,buffer_width),rep(0,N),rep(N,buffer_width))+buffer_width
  points_2_eval <- seq(1+buffer_width,N+buffer_width)
}
if(metric=="max"){
  cd <- 1
  rotate_ranks_45 <- as.data.frame(rotate(ranks_u,45))
  neighbors <-  RANN1::nn2(rotate_ranks_45,
                           rotate_ranks_45[points_2_eval,],
                           k=(k+1))
  d <- round(neighbors$nn.dists[,4]/sqrt(1^2+1^2),0)
}else if(metric=="euclid"){
  cd <- pi/4
  neighbors <-  RANN::nn2(ranks_u,
                          ranks_u[points_2_eval,],
                          k=(k+1))
  d <- neighbors$nn.dists[,4]
}

if(periodic){
  u_dist_nearest_border <- pmin(seq(1,N)+buffer_width-1,N+buffer_width-seq(1,N))
  if(any(u_dist_nearest_border<d)){
    ind_pot_outside_buff <- which(u_dist_nearest_border<d)
    new_width <- ceiling(max(d[ind_pot_outside_buff]))
    new_ranks_u <- round(ranks[order(ranks[,1]),],0)
    new_ranks_u <- rbind(tail(new_ranks_u,new_width),new_ranks_u,head(new_ranks_u,new_width))
    new_ranks_u[,1] <- new_ranks_u[,1]+c(rep(-N,new_width),rep(0,N),rep(N,new_width))+new_width
    new_points_2_eval <- seq(1+new_width,N+new_width)[ind_pot_outside_buff]


    if(metric=="max"){
      new_rotate_ranks_45 <- as.data.frame(rotate(new_ranks_u,45))
      new_neighbors <-  RANN1::nn2(new_rotate_ranks_45,
                                   new_rotate_ranks_45[new_points_2_eval,],
                                   k=(k+1))
      d_new <- round(new_neighbors$nn.dists[,4]/sqrt(1^2+1^2),0)
    }else if(metric=="euclid"){
      new_neighbors <-  RANN::nn2(new_ranks_u,
                                  new_ranks_u[new_points_2_eval,],
                                  k=(k+1))
      d_new <- new_neighbors$nn.dists[,4]
    }
    d[ind_pot_outside_buff] <- d_new
  }
}

digamma(N) - digamma(k)+log(cd) +2/N*sum(log(2*d)) - 2*log(N+1)
}



buffer_width <- 0.5
buffer_width <- ceiling(buffer_width*N)
ranks_u <- round(ranks[order(ranks[,1]),],0)
ranks_u <- rbind(tail(ranks_u,buffer_width),ranks_u,head(ranks_u,buffer_width))
ranks_u[,1] <- ranks_u[,1]+c(rep(-N,buffer_width),rep(0,N),rep(N,buffer_width))+buffer_width
rotate_ranks_45 <- as.data.frame(rotate(ranks_u,45))
neighbors <-  RANN1::nn2(rotate_ranks_45,
                         head(rotate_ranks_45[-seq(1,buffer_width),],N),
                         k=(k+1))
d <- round(neighbors$nn.dists[,4]/sqrt(1^2+1^2),0)

u_dist_nearest_border <- pmin(seq(1,N)+buffer_width-1,N+buffer_width-seq(1,N))
if(any(u_dist_nearest_border<d)){
  ind_pot_outside_buff <- which(u_dist_nearest_border<d)
  new_width <- ceiling(max(d[ind_pot_outside_buff]))
  new_ranks_u <- round(ranks[order(ranks[,1]),],0)
  new_ranks_u <- rbind(tail(new_ranks_u,new_width),new_ranks_u,head(new_ranks_u,new_width))
  new_ranks_u[,1] <- new_ranks_u[,1]+c(rep(-N,new_width),rep(0,N),rep(N,new_width))+new_width
  new_rotate_ranks_45 <- as.data.frame(rotate(new_ranks_u,45))
  new_neighbors <-  RANN1::nn2(new_rotate_ranks_45,
                           head(new_rotate_ranks_45[-seq(1,new_width),],N)[ind_pot_outside_buff,],
                           k=(k+1))
  d_new <- round(new_neighbors$nn.dists[,4]/sqrt(1^2+1^2),0)
  d[ind_pot_outside_buff] <- d_new
}
log_eps <- sum(log(2*d))
log_eps
cd <- 1


tic()
buffer_width <- 0.01
buffer_width <- ceiling(buffer_width*N)
ranks_u <- round(ranks[order(ranks[,1]),],0)
ranks_u <- rbind(tail(ranks_u,buffer_width),ranks_u,head(ranks_u,buffer_width))
ranks_u[,1] <- ranks_u[,1]+c(rep(-N,buffer_width),rep(0,N),rep(N,buffer_width))+buffer_width
neighbors <-  RANN::nn2(ranks_u,
                         head(ranks_u[-seq(1,buffer_width),],N),
                         k=(k+1))
d <- neighbors$nn.dists[,4]/(N+1)

u_dist_nearest_border <- pmin(seq(1,N)+buffer_width-1,N+buffer_width-seq(1,N))
if(any(u_dist_nearest_border<d)){
  ind_pot_outside_buff <- which(u_dist_nearest_border<d)
  new_width <- ceiling(max(d[ind_pot_outside_buff]))
  new_ranks_u <- round(ranks[order(ranks[,1]),],0)
  new_ranks_u <- rbind(tail(new_ranks_u,new_width),new_ranks_u,head(new_ranks_u,new_width))
  new_ranks_u[,1] <- new_ranks_u[,1]+c(rep(-N,new_width),rep(0,N),rep(N,new_width))+new_width
  new_neighbors <-  RANN::nn2(new_ranks_u,
                               head(new_ranks_u[-seq(1,new_width),],N)[ind_pot_outside_buff,],
                               k=(k+1))
  d_new <- new_neighbors$nn.dists[,4]/(N+1)
  d[ind_pot_outside_buff] <- d_new
}
log_eps <- sum(log(2*d))
log_eps
cd <- pi/4
toc()


digamma(N) - digamma(k)+log(cd) +2/N*log_eps - 2*log(101)


rotate <- function(x, degree) {
  degree <- pi * degree / 180
  l <- sqrt(x[,1]^2 + x[,2]^2)
  teta <- atan(x[,2] / x[,1])
  x[,1] <- l * cos(teta - degree)
  x[,2] <- l * sin(teta - degree)
  return(x)
}
