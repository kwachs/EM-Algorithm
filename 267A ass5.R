#•Write a function that implements the EM-algorithm for GMM (thelog-likelihood
#have to increase at each iteration)•Test it on 2 data sets (real or
#simulated)•submit the function on CANVAS as .r file•Do not include your name in
#the file.•After the submission you will have to debug and grade one ofyour
#classmate homework (randomly assigned).•I’ll regrade the homework and it will
#count in your averageassignment grade.•Due on Tuesday April 24.•Review due on
#Thursday April 26.

#ass 5
#output is list of Z and parameter estimates (mean, variance) of each cluster
#bishop: pattern recognition and machine learning chapter 9
#can either randomly start mean, variance and proporation, or randomly start 
# z-- this is faster. whole thing must sum to n and each row must sum to 1
#dmvnorm for likelihood (f(xi|MUg,SIGMAg))
#product(i=1,n)sum(g=1,G)PIg*(f(xi|MUg,SIGMAg))

iris1=iris[,-5]
faithful

ass5=function(data,G){
  set.seed(21)
  pis=rep(0,G)
  mu=matrix(0,nrow = G,ncol=ncol(data))
  sigma=array(0,dim=c(ncol(data),ncol(data),G))
  z=matrix(0,nrow=nrow(data),G)
  zig=matrix(0,nrow=nrow(data),G) #shows 1's and zeros for which group e/ obs is in
  ng=rep(0,G)#size of each group
  error=.0005
  likelihood=rep(0,100)
  output=list(z=matrix(0,nrow=nrow(data),G),
                mu=matrix(0,nrow = G,ncol=ncol(data)),pi=rep(0,G),
                sigma=array(0,dim=c(ncol(data),ncol(data),G)),loops=0)
  
  likelihood[1]=-1000000
  #initialize z, sort of like a step 0 + E-step
  for(i in 1:nrow(z)){
    temp=runif(G)
     z[i,]=temp/sum(temp)
     zig[i,which.max(z[i,])]=1
  }

#zig matrix helps for making colors for each cluster
  for(x in 2:length(likelihood)){
    zig=matrix(0,nrow=nrow(data),G)
    for(i in 1:nrow(z)){
      zig[i,which.max(z[i,])]=1
    }
  ng=apply(z,2,sum)
  #step 2: M-step
  #pi and mu
  for (g in 1:G) {
    sum_mu=rep(0,ncol(data))
    sum_pi=0
  for(j in 1:nrow(z)){
    sum_mu=sum_mu+(z[j,g]*as.numeric(data[j,]))
    sum_pi=sum_pi+(z[j,g]/nrow(data))
  }
   mu[g,]=sum_mu/ng[g]
   pis[g]=sum_pi
  }
  #sigma
  for (g in 1:G) {
  sum_sigma=matrix(0,nrow=ncol(data),ncol(data))
  for(j in 1:nrow(z)){
  sum_sigma=sum_sigma+((z[j,g])*(((as.numeric(data[j,])-mu[g,]))%*%
                         t((as.numeric(data[j,])-mu[g,]))))
  }
    sigma[,,g]=(sum_sigma/ng[g]) 
  }
  #finding likelihood
  fl=rep(0,nrow(z))
  for(j in 1:nrow(z)){
  for (g in 1:G) {
  fl[j]=fl[j]+pis[g]*dmvnorm(as.numeric(data[j,]),mean=mu[g,],sigma=sigma[,,g])
  }
  }
  likelihood[x]=sum(log(fl))
  #update z
  for(j in 1:nrow(z)){
    for (g in 1:G) {
      denom=0
      for (h in 1:G) {
       denom=denom+
         (pis[h]*dmvnorm(as.numeric(data[j,]),mean=mu[h,],sigma=sigma[,,h]))
      }
  z[j,g]=(pis[g]*dmvnorm(as.numeric(data[j,]),mean=mu[g,],sigma=sigma[,,g]))/(denom)
    }
  }
  #test for likelihood increase plateau
  if(likelihood[x]-likelihood[x-1]<error){
    break
  }
  }
  #colors for the clusters!
  clusters=rep(0,nrow(zig))
  for(r in 1:length(clusters)){
    for(c in 1:ncol(zig)){
    if(zig[r,c]==1){
      clusters[r]=c
    }
    }
  }
  plot(data,col=clusters)
  output$z=z
  output$mu=mu
  output$pi=pis
  output$sigma=sigma
  output$loops=x-1
  return((output))
}
ass5(iris1,3)
ass5(faithful,2)


