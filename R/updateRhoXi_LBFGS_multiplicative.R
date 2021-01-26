updateRhoXi_LBFGS = function(
Yt,
YMat,
Data,
Param,
DEBUG=F,
Verbose=0){
    
    if(Verbose>0)print("update Xi & Rho with LBFGS")
    
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
    M1 = nrow(Ta[[1]])
    M2 = nrow(Ta[[2]])
    DiagM2 = array(1,c(M1,1))%x%diag(M2)
    
    
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
    C = C + theta*KinvKmn # titsias
    C = t(C)/omega2
    
    AA= Solve(K) - Solve(Phiinv)[flagZK==1,flagZK==1]/theta
    AA= AA - Solve(Phiinv, t(Solve(Phiinv, G%*%t(t(D)-E%*%tA))))[flagZK==1,flagZK==1]/J/theta
    AA= AA - theta*KinvKmn%*%(t(KinvKmn)/omega2) # titsias
    

    Ta1rep = rep(c(Ta[[1]]),rep(M2,M1))
    # Xi
    gradXi = list(rowSums((C*Knm)*(2*cos(-outer(c(Xi[[1]]),Ta1rep,"-"))*sin(-outer(c(Xi[[1]]),Ta1rep,"-"))))/rho[[1]]^2,
        (- rowSums(C*Knm)*Xi[[2]] + (C*Knm)%*%Rep(Ta[[2]],M1))%*%diag(1/rho[[2]]^2) - Xi[[2]]/J)
    
    # rho
    dqK   = sin(outer(Ta1rep,Ta1rep,"-"))^2/rho[[1]]^2. # K  *
    dqKnm = sin(outer(c(Xi[[1]]),Ta1rep,"-"))^2/rho[[1]]^2. # Knm*
    gradRho = list(sum((C*Knm)*dqKnm) + sum((AA*K)*dqK)/2, rep(0,Q))
    for(q in seq(Q)){
        Ta2rep = rep(Ta[[2]][,q], M1)
        dqK   = outer(Ta2rep,Ta2rep,"-")^2/rho[[2]][q]^2/2. # K  *
        dqKnm = outer(Xi[[2]][,q],Ta2rep,"-")^2/rho[[2]][q]^2/2. # Knm*
        gradRho[[2]][q] = sum((C*Knm)*dqKnm) + sum((AA*K)*dqK)/2
    }
    
    # Ta
    DiagM2 = Rep(diag(M2),M1)
    AAA = t(DiagM2)%*%(AA*K)%*%DiagM2
    gradTa = (t((C*Knm)%*%DiagM2)%*%Xi[[2]] - colSums((C*Knm)%*%DiagM2)*Ta[[2]])%*%diag(1/rho[[2]]^2) - Ta[[2]]/J + (AAA-diag(rowSums(AAA)))%*%Ta[[2]]%*%diag(1/rho[[2]]^2)
    
    gk = rbind(t(gradRho[[2]]),gradXi[[2]],gradTa)
    if(is.null(Param$GradRhoXi)){
        print("init LBFGS")
        ValNew  = matrix(c(2*log(rho[[1]]), c(Xi[[1]]), c(rbind(t(2*log(rho[[2]])),Xi[[2]],Ta[[2]]))), 1+N+(1+N+M2)*Q)
        GradNew = matrix(c(gradRho[[1]], gradXi[[1]], c(gk)),1+N+(1+N+M2)*Q)
        ss = c(GradNew)/10000
        ss = list(ss[1:(N+1)], matrix(ss[-seq(N+1)], 1+N+M2))
    }else{
        ValNew  = cbind(c(2*log(rho[[1]]), c(Xi[[1]]), c(rbind(t(2*log(rho[[2]])),Xi[[2]],Ta[[2]]))), Param$ValRhoXi)
        GradNew = cbind(c(gradRho[[1]], gradXi[[1]], c(gk)),Param$GradRhoXi)
        #ss = matrix(-LBFGS(ValNew,GradNew),1+N+M2)
        ss = -LBFGS(ValNew,GradNew)
        ss = list(ss[1:(N+1)], matrix(ss[-seq(N+1)], 1+N+M2))
    }
    
    tLD = dbind(Param$LD, backsolve((chol(K)),diag(M))*sqrt(theta))
    pelbonow = sum(c(- logDet(diag(K2+M)+t(tLD)%*%(t(tZ/omega2)%*%tZ)%*%tLD), + sum(diag(Solve(Phiinv,G%*%t(t(D)-E%*%tA)/J))), + sum(t(Knm/omega2)*Solve(K,t(Knm)))*theta, - sum(Xi[[2]]^2)/J - sum(Ta[[2]]^2)/J)) #- sum(1/omega2)*theta 
    pelbo0=pelbo1=-(10^10)
    rho0 = rho
    Xi0  = Xi
    Ta0  = Ta
    if(Verbose>1){cat("ss=");print(c(ss[[1]][1],ss[[2]][1,]))}
    for(r in c(50,20,10:0,-1)){
        rho1 = list(exp(log(rho[[1]]) - (ss[[1]][1]/2.)/2^r), exp(log(rho[[2]]) - (ss[[2]][1,]/2.)/2^r))
        Xi1  = list(Xi[[1]] - ss[[1]][-1]/2^r,                Xi[[2]] - ss[[2]][2:(1+N),]/2^r)
        Ta1  = list(Ta[[1]],                                  Ta[[2]] - ss[[2]][(N+2):(1+N+M2),]/2^r)
        if(Verbose>1){cat("rho1=");print(unlist(rho1))}
        pelbo1 = getELBOrho(Yt, YMat, Data, Param, Xi1, Ta1, rho1)
        if(Verbose>1){print(c(as.integer(r),pelbo0,pelbo1))}
        if(is.na(pelbo1)){break}
        if(pelbo0>pelbo1 && pelbonow>pelbo1*1.0000001){
            break
        }else{
            pelbo0=pelbo1
            rho0 = rho1
            Xi0  = Xi1
            Ta0  = Ta1
        }
    }
    if(Verbose>0){cat("rho=");print(unlist(rho0))}
    list(rho=rho0, Xi=Xi0, Ta=Ta0, ValRhoXi = ValNew[,seq(min(25,ncol(ValNew))),drop=F], GradRhoXi = GradNew[,seq(min(25,ncol(ValNew))),drop=F])
}



