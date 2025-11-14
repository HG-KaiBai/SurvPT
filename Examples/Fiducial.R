# R code for Algorithm 1 of "A unified nonparametric fiducial approach to interval-censored data, Cui et al. (2023)"
#########

fiducial <- function(l, r, nsample, grid, testgrid, mfid = 1000, mburn = 100) {
  ngrid = length(grid)
  ntestgrid = length(testgrid)
  mid=(l+r)/2
  vu <- sort(runif(nsample))
  vu[order(mid)] <- vu # initialize an u, Lines 1-3

  fid.u <- matrix(NA,mfid+mburn,nsample)
  fid.lower <- matrix(NA,mfid+mburn,ngrid)
  fid.upper <- matrix(NA,mfid+mburn,ngrid)
  fid.lower.test <- matrix(NA,mfid+mburn,ntestgrid)
  fid.upper.test <- matrix(NA,mfid+mburn,ntestgrid)


  for (j in 1:(mfid+mburn)){ # repeat for each fiducial sample
    for (i in 1:nsample){ # Gibbs sampler corresponding to Lines 6-11
      vu.pre <- vu[-i] # Line 7
      l.pre <- l[-i]
      r.pre <- r[-i]
      index1 <- which(l[i]>=r.pre)
      if(length(index1)>0){
      u.lower <- max(vu.pre[index1]) }
      else u.lower <- 0  ##### Line 8
      index2 <- which(r[i]<=l.pre)
      if(length(index2)>0){
      u.upper <- min(vu.pre[index2])
      } else u.upper <- 1 ##### Line 9
      vu[i] <- runif(1,min=u.lower, max=u.upper) # Line 10
    }

    temp <- runif(nsample)
    vu <- sort(temp)[rank(vu)]
    fid.u[j,] <- vu # Line 12, update u according to the order


    # Post-processing Lines 15-20 for grid
    fiducial.sample=sapply(grid,function(s){
      index1 <- which(s>=r)
      if(length(index1)>0){
        lower <- max(vu[index1]) }
      else lower <- 0  ##### Line 18
      index2 <- which(s<l)
      if(length(index2)>0){
        upper <- min(vu[index2])
      } else upper <- 1 ##### Line 19
      return(rbind(lower,upper))
    })

    fid.lower[j,]=fiducial.sample[1,]
    fid.upper[j,]=fiducial.sample[2,]

    # Repeat Post-processing for testgrid
    fiducial.sample.test=sapply(testgrid,function(s){
      index1 <- which(s>=r)
      if(length(index1)>0){
        lower <- max(vu[index1]) }
      else lower <- 0
      index2 <- which(s<l)
      if(length(index2)>0){
        upper <- min(vu[index2])
      } else upper <- 1
      return(rbind(lower,upper))
    })

    fid.lower.test[j,]=fiducial.sample.test[1,]
    fid.upper.test[j,]=fiducial.sample.test[2,]
  }

  ############# fiducial estimator##############
  fid.grid=cbind(fid.lower,fid.upper)[(mburn+1):(mburn+mfid),]
  linefid1=t(apply(fid.grid,1,linear_interpolation,grid=grid))
  FiducialMidLine1= apply(linefid1,1,function(x){approx(grid,x,testgrid, method = 'linear', rule = 2)$y})
  linefidi=apply(FiducialMidLine1,1,sort)
  point_li=linefidi[mfid*0.5,] # point estimator based on linear interpolation

  return(list(point_li = point_li, lower = linefidi[mfid * 0.025, ], upper = linefidi[mfid * 0.975, ]))
}



