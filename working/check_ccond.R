lin_cop <- claytonCopula(4)
u <- cbind(c(0,0.5,1,0,0.5,1,0,0.5,1),c(0,0,0,1,1,1,0.5,0.5,0.5))
u <- rbind(u,cbind(runif(20),runif(20)))
cop_lst <- list( lin_cop,
                 rotCopula(lin_cop,flip=c(T,T)),
                 rotCopula(lin_cop,flip=c(T,F)),
                 rotCopula(lin_cop,flip=c(F,T)),
                 cyl_rot_combine(lin_cop,shift = F),
                 cyl_rot_combine(lin_cop,shift = T),
                 cyl_quadsec(0.1),
                 cyl_quadsec(0),
                 cyl_cubsec(0.1),
                 cyl_cubsec(-0.1),
                cyl_rect_combine(lin_cop,flip_up = T),
                cyl_rect_combine(lin_cop,low_rect = c(0.1,0.4),flip_up = T),
                cyl_rect_combine(lin_cop,low_rect = c(0.1,0.35),up_rect = c(0.7,0.9),flip_up = T),
                cyl_rect_combine(lin_cop,background = cyl_quadsec(0.1),low_rect = c(0.1,0.35),up_rect = c(0.7,0.9),flip_up = T),
                cyl_rect_combine(lin_cop,flip_up = F),
                cyl_rect_combine(lin_cop,low_rect = c(0.1,0.4),flip_up = F),
                cyl_rect_combine(lin_cop,low_rect = c(0.1,0.35),up_rect = c(0.7,0.9),flip_up = F),
                cyl_rect_combine(lin_cop,background = cyl_quadsec(0.1),low_rect = c(0.1,0.35),up_rect = c(0.7,0.9),flip_up = F))

res <- rep(list(NULL),4)
res[[1]] <- rep(list(NULL),length(cop_lst))
res[[2]] <- res[[1]]
res[[3]] <- res[[1]]
res[[4]] <- res[[1]]
res[[5]] <- res[[1]]
res[[6]] <- res[[1]]

cond_on <- 1
inverse <- F
  for (i in 1:length(cop_lst)) {
    cop <- cop_lst[[i]]
    anal <- ccylcop(u,cop,cond_on,inverse)
    num <- numerical_conditional_cop(u,cop,cond_on)
    res[[1]][[i]] <- anal - num
  }

  cond_on <- 2
  inverse <- F
  for (i in 1:length(cop_lst)) {
    cop <- cop_lst[[i]]
    anal <- ccylcop(u,cop,cond_on,inverse)
    num <- numerical_conditional_cop(u,cop,cond_on)
    res[[2]][[i]] <- anal - num
  }


   cond_on <- 1
  inverse <- T
  for (i in 1:length(cop_lst)) {
    cop <- cop_lst[[i]]
    anal <- ccylcop(u,cop,cond_on,inverse)
    num <- numerical_inv_conditional_cop(u,cop,cond_on)
    res[[3]][[i]] <- anal - num
  }

  cond_on <- 2
  inverse <- T
  for (i in 1:length(cop_lst)) {
    cop <- cop_lst[[i]]
    anal <- ccylcop(u,cop,cond_on,inverse)
    num <- numerical_inv_conditional_cop(u,cop,cond_on)
    res[[4]][[i]] <- anal - num
  }


for (i in 1:length(cop_lst)) {
  res[[5]][[i]] <- rcylcop(1000,cop_lst[[i]])
}
