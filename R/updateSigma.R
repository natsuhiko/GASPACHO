updateSigma = function(
    YMat,
    Data,
    Param,
    DEBUG=F, Yt, Verbose=0){
        
    if(Verbose>0)print("updateSigma")
    
    N = Data$MatrixDims$N
        
    Knm = Param[["Knm"]]
    K = Param[["K"]]
    KinvKmn = Solve(K,t(Knm))
    Z = Data[["Z"]]
    W = Data[["W"]]
    zeta = Param[["zeta"]]
    
    omega2 = Param[["omega2"]]
    #sigma2 = Param[["sigma2"]]
    Phiinv = Param[["Phiinv"]]
    Theta = rep(Param[["theta"]], Data$MatrixDims$M)
    
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    tW = cbind(W,Z%*%zeta)
    
    
    C = cbind(YMat$YtOinvW,YMat$YtOinvZ%*%zeta)
    D = cbind(YMat$YtOinvZ,YMat$YtOinvKnm) # J x (K2+M)
    E = t(tZ/omega2)%*%tW # (K2+M) x (K1+1)
    F = Solve(Phiinv, t(D)) # (K2+M) x J
    G = Solve(Phiinv, E)
    
    sigma2 = YMat$Yt2Oinv1 - 2*colSums(t(C)*tA) + colSums(tA*((t(tW/omega2)%*%tW)%*%tA))
    sigma2 = sigma2 - colSums(t(D)*F) + 2*colSums(tA*(t(E)%*%F)) - colSums(tA*(t(E)%*%G%*%tA))
    sigma2 = (sigma2+1)/(N+1)
    if(DEBUG){
        tD = dbind(Param$Delta,Solve(K)*Theta)
        V=diag(omega2) + tZ%*%tD%*%t(tZ)
        Vinv=Solve(V)
        tY=(t(Yt)-tW%*%tA)
        S=t(tY)%*%Vinv%*%tY/N
        plot(sigma2, diag(S));abline(0,1)
    }
    
    
    list(sigma2=sigma2)
}

getNewSigma = function(
Yt,
Data,
Param, Verbose=0){
    
    if(Verbose>0)print("updateSigma")
    
    N = Data$MatrixDims$N
    
    Knm = Param[["Knm"]]
    K = Param[["K"]]
    KinvKmn = Solve(K,t(Knm))
    Z = Data[["Z"]]
    W = Data[["W"]]
    zeta = Param[["zeta"]]
    
    omega2 = Param[["omega2"]]
    Phiinv = Param[["Phiinv"]]
    Theta = rep(Param[["theta"]], Data$MatrixDims$M)
    
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    tW = cbind(W,Z%*%zeta)
    
    YtOinvZ = Yt%*%(Z/omega2)
    C = YtOinvZ%*%zeta
    D = cbind(YtOinvZ,Yt%*%(Knm/omega2)) # J x (K2+M)
    E = t(tZ/omega2)%*%tW # (K2+M) x (K1+1)
    F = Solve(Phiinv, t(D)) # (K2+M) x J
    G = Solve(Phiinv, E)
    
    sigma2 = as.numeric((Yt^2)%*%(1/Param[["omega2"]])) - 2*colSums(t(C)*tA) + colSums(tA*((t(tW/omega2)%*%tW)%*%tA))
    sigma2 = sigma2 - colSums(t(D)*F) + 2*colSums(tA*(t(E)%*%F)) - colSums(tA*(t(E)%*%G%*%tA))
    sigma2 = (sigma2+1)/(N+1)
    
    
    sigma2
}
