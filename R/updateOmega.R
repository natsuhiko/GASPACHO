updateOmega = function(
Yt,
YMat,
Data,
Param,
DEBUG=F, Verbose=0){
    
    if(Verbose>0)print("updateOmega")
    
    J = Data$MatrixDims$J
    
    Knm = Param[["Knm"]]
    K = Param[["K"]]
    KinvKmn = Solve(K,t(Knm))
    Z = Data[["Z"]]
    W = Data[["W"]]
    zeta = Param[["zeta"]]
    
    omega2 = Param[["omega2"]]
    sigma2 = Param[["sigma2"]]
    Phiinv = Param[["Phiinv"]]
    theta = Param[["theta"]]
    Theta = rep(theta, Data$MatrixDims$M)
    
    tA = rbind(Param[["Alpha"]],1)
    tASinv = t(t(tA)/sigma2)
    tZ = cbind(Z,Knm)
    tW = cbind(W,Z%*%zeta)
    
    D = cbind(YMat$YtOinvZ,YMat$YtOinvKnm)
    E = t(tZ/omega2)%*%tW
    FF= -tW + tZ%*%Solve(Phiinv, E)
    G = Solve(Phiinv, t(tZ))
    H = rbind(YMat$ASinvYt,YMat$OneSinvYt)
    
    v = 1/omega2 - colSums(G*t(tZ))/omega2^2
    
    vsv = as.numeric(t(1/sigma2)%*%(Yt^2))/omega2^2
    vsv = vsv - 2*colSums( G * as.matrix(t(D/sigma2)%*%Yt) )/omega2^2
    vsv = vsv +   colSums( G * ((t(D/sigma2)%*%D)%*%G)     )/omega2^2
    vsv = vsv + 2*colSums(t(FF) * H)/omega2^2
    vsv = vsv - 2*colSums(t(FF) * (tASinv%*%D)%*%G)/omega2^2
    vsv = vsv +   colSums(t(FF) * (tASinv%*%t(tA))%*%t(FF))/omega2^2
    vsv = vsv/J
    
    if(DEBUG){
        tD = dbind(Param$Delta,Solve(K)*Theta)
        V=diag(omega2)+tZ%*%tD%*%t(tZ)
        Vinv=Solve(V)
        Vinv2 = diag(1/omega2) - diag(1/omega2)%*%tZ%*%Solve(Phiinv)%*%t(tZ)%*%diag(1/omega2)
        plot(diag(Vinv),diag(Vinv2))
        tY=(t(Yt)-tW%*%tA)
        S=tY%*%diag(1/sigma2)%*%t(tY)/J
        plot(vsv,diag(Vinv%*%S%*%Vinv));abline(0,1)
        plot(v, diag(Vinv));abline(0,1)
    }
    
    grad = c(- v*sqrt(omega2) + vsv*sqrt(omega2) )
    hess = c(- 2*v^2*omega2 )
    
    # titsias
    if(Data$Kernel!="Linear"){
        grad = grad +   c(sum(theta)-(Knm*t(KinvKmn))%*%Theta)/sqrt(omega2)^3
        hess = hess - 3*c(sum(theta)-(Knm*t(KinvKmn))%*%Theta)/omega2^2
    }
    
    aig = 6; big=5
    grad = grad + (2*big-2*(1+aig)*omega2)/sqrt(omega2)^3/J
    hess = hess + 2*(-3*big+(1+aig)*omega2)/omega2^2/J
    
    sumto1 = function(x){x/mean(x)}
    omega22 = c((sqrt(omega2) - grad/hess/2)^2)
    #omega22[omega22/omega2>10]  = omega2[omega22/omega2>10]*10
    #omega22[omega22/omega2<0.1] = omega2[omega22/omega2<0.1]*0.1
    list(omega2 = sumto1(omega22))
}

