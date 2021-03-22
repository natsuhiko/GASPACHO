#MPInverse=function(x){
#    x = eigen((x+t(x))/2)
#    l = 1/x[[1]]
#    l[x[[1]]<1e-7]=0
#    x[[2]]%*%t(x[[2]]*l)
#}

reducedDim = function(Ta){
    Tall = cbind(Ta[[1]],Ta[[2]])
    plot(Tall)
    grp = cutree(hclust(dist(Tall)),h=0.2)
    w = as.numeric(1/table(grp)[grp])
    Tall = t(t(Tall)%*%(diag(length(unique(grp)))[grp,]*w))
    points(Tall,col=2)
    Ta[[1]]=Tall[,1,drop=F]
    Ta[[2]]=Tall[,2:ncol(Tall),drop=F]
    Ta
}

getTa = function(Xi, M=50){
    apply(Xi,2,function(x1){
        r = range(x1)
        r = r + c(-0.1,0.1)*diff(r)
        ppoints(M)[order(runif(M))]*diff(r)+r[1]
    })
}

updateKernel=function(Data, Param){
    updateKernelMultiplicative(Data, Param)
}

updateKernelMultiplicative=function(Data, Param){
    K1 = updateKernelPeriodic(Param$Xi[[1]], Param$Ta[[1]], Param$rho[[1]])
    K = K1$K
    Knm = K1$Knm
    for(i in 2:length(Param$Ta)){
        K2 = updateKernelSE(Param$Xi[[i]], Param$Ta[[i]], Param$rho[[i]])
        K = K * K2$K
        Knm = Knm * K2$Knm
    }
    #K3 = updateKernelSE(Param$Xi[[3]], Param$Ta[[3]], Param$rho[[3]])
    list(
        K = K,
        Knm = Knm
    )
}

updateKernelAdditive=function(Data, Param){
    K1 = updateKernelPeriodic(Param$Xi[[1]], Param$Ta[[1]], Param$rho[[1]])
    K2 = updateKernelSE(Param$Xi[[2]], Param$Ta[[2]], Param$rho[[2]])
    K3 = updateKernelSE(Param$Xi[[3]], Param$Ta[[3]], Param$rho[[3]])
    list(
        K = dbind(K1$K, K2$K, K3$K),
        Knm = cbind(K1$Knm, K2$Knm, K3$Knm)
    )
}

updateKernelPeriodic=function(Xi, Ta, rho){
    K  = exp(-sin(outer(c(Ta),c(Ta),"-"))^2/rho^2)
    Knm = exp(-sin(outer(c(Xi),c(Ta),"-"))^2/rho^2)
    list(
    K = K,
    Knm = Knm
    )
}

updateKernelSE=function(Xi, Ta, rho){
    K  = getKmm(Ta, rho)
    Knm = getKnm(Xi, Ta, rho)
    list(
        K = K,
        Knm = Knm
    )
}

updateKernelSE1=function(Xi, Ta, rho){
    K  = getKmm(Ta, rho)
    Knm = getKnm(Xi, Ta, rho)
    list(
    K = K,
    Knm = Knm
    )
}

updateKernelSE2=function(Xi, Ta, rho){
    K  = getKmm(Ta, rho)
    Knm = getKnm(Xi, Ta, rho)
    list(
    K = K,
    Knm = Knm
    )
}

getKmm = function(Ta, rho){
    H = Ta%*%(t(Ta)/rho^2)
    h = diag(H)
    
    exp((t(2*H - h) - h)/2)
}




getKnm = function(X, Ta, rho){
    # C : N x M
    # X : N x Q
    # Z : M x Q
    Q = ncol(Ta)
    f = c((X^2)%*%(1/rho^2)) # N
    g = c((Ta^2)%*%(1/rho^2)) # M
    H = X%*%(t(Ta)/rho^2) # N x M
    
    exp(t(t(2*H - f) - g)/2)
}


