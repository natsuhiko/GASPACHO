updateTLB = function(
    YMat,
    Data,
    Param,
    DEBUG=F, Verbose=0){
    if(Data$Kernel=="SE"){
        updateTLBSE(YMat, Data, Param, DEBUG, Verbose)
    }else{
        updateTLBLinear(YMat, Data, Param, DEBUG, Verbose)
    }
}

updateTLBLinear = function(
YMat,
Data,
Param,
DEBUG=F, Verbose=0){
    
    if(Verbose>0)print("updateTLB")
    
    N = Data$MatrixDims$N
    M = Data$MatrixDims$M
    Q = Data$MatrixDims$Q
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    
    Knm = Param[["Knm"]]
    K = Param[["K"]]

    if(sum(Q)>0){ RK = chol(K) }
    Z = Data[["Z"]]
    
    omega2 = Param[["omega2"]]
    sigma2 = Param[["sigma2"]]
    theta = Param[["theta"]]
    Theta = rep(theta, Q)
    tZ = cbind(Z,Knm)
    
    LD = Param[["LD"]]
    if(sum(Q)>0){
        LtD= dbind(LD,sqrt(Theta)*(backsolve(RK,diag(sum(M)))))
    }else{
        LtD = LD
    }
    tlb = c(-N/2*sum(log(sigma2)), - J/2*sum(log(omega2)), - J/2*logDet(diag(K2+sum(M)) + t(LtD)%*%(t(tZ/omega2)%*%tZ)%*%LtD))
    sum(tlb)
}

updateTLBSE = function(
    YMat,
    Data,
    Param,
    DEBUG=F, Verbose=0){
    
    if(Verbose>0)print("updateTLB")
    
    N = Data$MatrixDims$N
    M = Data$MatrixDims$M
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    
    Knm = Param[["Knm"]]
    K = Param[["K"]]
    KinvKmn = Solve(K,t(Knm))
    
    RK = chol(K)
    Z = Data[["Z"]]
    
    omega2 = Param[["omega2"]]
    sigma2 = Param[["sigma2"]]
    theta = Param[["theta"]]
    Theta = rep(theta, M)
    tZ = cbind(Z,Knm)
    
    LD = Param[["LD"]]
    LtD= dbind(LD,sqrt(Theta)*(backsolve(RK,diag(sum(M)))))
    
    tlb = c(-N/2*sum(log(sigma2)), - J/2*sum(log(omega2)), - J/2*logDet(diag(K2+sum(M)) + t(LtD)%*%(t(tZ/omega2)%*%tZ)%*%LtD),
        - J/2*sum(theta * sum(1/omega2)) + J/2*sum( ((Knm*t(KinvKmn))%*%Theta)/omega2 ), # titsias
        - sum((Param$Xi[[2]]-Param$PriorXi2$mu)^2/Param$PriorXi2$var)/2 - sum(Param$Ta[[2]]^2)/2)
    sum(tlb)
}


getTLBrho = function(Yt, YMat, Data, Param, Xi1, Ta1, rho1){
    if(sum(unlist(rho1)>2000,na.rm=T)>1){return(-10^10)}
    M = Data$MatrixDims$M
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    omega2 = Param$omega2
    sigma2 = Param$sigma2
    Z = Data$Z
    Theta = rep(Param[["theta"]], M)
    
    newK = updateKernel(Data,list(Xi=Xi1, Ta=Ta1, rho=rho1))
    tLD = dbind(Param$LD, backsolve((chol(newK$K)),diag(sum(M)))*sqrt(Theta))
    tDinv  = dbind(Param$Dinv, newK$K/Theta)
    tZ = cbind(Z,newK$Knm)
    tW = cbind(Data$W,Z%*%Param$zeta)
    tA = rbind(Param[["Alpha"]],1)
    Phiinv = tDinv+t(tZ/omega2)%*%tZ
    D = cbind(YMat$YtOinvZ,as.matrix(Yt%*%(newK$Knm/omega2)))
    E = t(tZ/omega2)%*%tW
    G = (t(D)-E%*%tA)%*%(t(t(D)-E%*%tA)/sigma2)/J
    
    tlb = c(- logDet(diag(K2+sum(M))+t(tLD)%*%(t(tZ/omega2)%*%tZ)%*%tLD)
            ,sum(diag(Solve(Phiinv,G)))
            ,sum(t(newK$Knm/omega2)*Solve(newK$K/Theta,t(newK$Knm))) # - sum(1/omega2)*Param$theta
            ,- sum((Xi1[[2]]-Param$PriorXi2$mu)^2/Param$PriorXi2$var)/J - sum(Ta1[[2]]^2)/J)
            #print(tlb)
    
    sum(tlb)
}









