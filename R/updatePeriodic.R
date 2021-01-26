updatePeriodic = function(
Yt,
YMat,
Data,
Param,
DEBUG=F,
Verbose=0){
    
    if(Verbose>0)print("update Xi & Rho with Periodic kernel")
    
    J = Data$MatrixDims$J
    N = Data$MatrixDims$N
    K2= Data$MatrixDims$K2
    M = Data$MatrixDims$M
    Q = Data$MatrixDims$Q
    
    flagZK = rep(0:1,c(K2,M))
    
    Knm = Param[["Knm"]]
    K = Param[["K"]]
    Z = Data[["Z"]]
    W = Data[["W"]]
    zeta = Param[["zeta"]]
    Xi = Param[["Xi"]]
    Ta = Param[["Ta"]]
    
    rho = Param[["rho"]]
    
    omega2 = Param[["omega2"]]
    sigma2 = Param[["sigma2"]]
    theta = Param[["theta"]]
    
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    tW = cbind(W,Z%*%zeta)
    Phiinv = dbind(Param[["Dinv"]],K/theta) + t(tZ/omega2)%*%tZ
    
    
    B = Solve(Phiinv, t(tZ))
    D = cbind(YMat$YtOinvZ,YMat$YtOinvKnm)
    E = t(tZ/omega2)%*%tW
    FF= -tW + tZ%*%Solve(Phiinv, E)
    G = t(t(t(D)-E%*%tA)/sigma2)
    H = rbind(YMat$ASinvYt,YMat$OneSinvYt)
    
    KinvKmn = Solve(K,t(Knm))
    
    
    # matrix prep for dKnm and dK
    C = Solve(Phiinv, -t(tZ) + (as.matrix(t(D/sigma2)%*%Yt) - E%*%H - (G%*%D)%*%B + (G%*%t(tA))%*%t(FF))/J)[flagZK==1,]
    C = C + theta*KinvKmn
    C = t(C)/omega2
    
    AA= Solve(K) - Solve(Phiinv)[flagZK==1,flagZK==1]/theta
    AA= AA - Solve(Phiinv, t(Solve(Phiinv, G%*%t(t(D)-E%*%tA))))[flagZK==1,flagZK==1]/J/theta
    AA= AA - theta*KinvKmn%*%(t(KinvKmn)/omega2)
    
    # Xi
    gradXi = rowSums(C*Knm*(2*cos(-outer(c(Xi[[1]]),c(Ta[[1]]),"-"))*sin(-outer(c(Xi[[1]]),c(Ta[[1]]),"-"))))/rho[[1]]^2
    
    # rho
    dqK   = K  *sin(outer(c(Ta[[1]]),c(Ta[[1]]),"-"))^2/rho[[1]]^2.
    dqKnm = Knm*sin(outer(c(Xi[[1]]),c(Ta[[1]]),"-"))^2/rho[[1]]^2.
    gradRho = sum(C*dqKnm) + sum(AA*dqK)/2
    
    gk = c(gradRho,gradXi)
    if(is.null(Param$GradRhoXi)){
        ValNew  = matrix(c(2*log(rho[[1]]),c(Xi[[1]])), N+1)
        GradNew = matrix(c(gk),N+1)
        ss = matrix(gk,N+1)
    }else{
        ValNew  = cbind(c(2*log(rho[[1]]),c(Xi[[1]])), Param$ValRhoXi)
        GradNew = cbind(c(gk),Param$GradRhoXi)
        ss = matrix(-LBFGS(ValNew,GradNew),1+N)
    }
    
    tLD = dbind(Param$LD, backsolve((chol(K)),diag(M))*sqrt(theta))
    pelbonow = c(- logDet(diag(K2+M)+t(tLD)%*%(t(tZ/omega2)%*%tZ)%*%tLD), + sum(diag(Solve(Phiinv,G%*%t(t(D)-E%*%tA)/J))), + sum(t(Knm/omega2)*Solve(K,t(Knm)))*theta, - sum(1/omega2)*theta - sum(Xi[[2]]^2)/J - sum(Ta[[2]]^2)/J )
    print(pelbonow)
    pelbonow=sum(pelbonow)
    pelbo0=pelbo1=-(10^10)
    rho0 = rho
    Xi0  = Xi
    if(Verbose){cat("ss=");print(ss[1,])}
    for(r in c(50,20,10:0,-1)){
        rho1 = list(exp(log(rho[[1]]) - (c(ss[1,])/2.)/2^r), rho[[2]])
        Xi1  = list((Xi[[1]] - ss[2:(N+1),]/2^r)%%pi, Xi[[2]])
        if(Verbose>1){cat("rho1=");print(rho1[[1]]);}
        pelbo1 = getELBOrho(Yt, YMat, Data, Param, Xi1, Ta, rho1)
        if(Verbose>1){print(c(as.integer(r),pelbo0,pelbo1))}
        if(is.na(pelbo1)){break}
        if(pelbo0>pelbo1 && pelbonow>pelbo1*1.0000001){
            break
        }else{
            pelbo0=pelbo1
            rho0 = rho1
            Xi0  = Xi1
        }
    }
    if(Verbose>0){cat("ssrho=");print(ss[1,]);cat("rho=");print(rho0[[1]])}
    list(rho=rho0, Xi=Xi0, Ta=Ta, ValRhoXi = ValNew[,seq(min(25,ncol(ValNew))),drop=F], GradRhoXi = GradNew[,seq(min(25,ncol(ValNew))),drop=F])
}



