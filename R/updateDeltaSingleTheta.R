initDelta = function(Data, Param){
    print("initDelta")
    nh = Data$nh
    ccat = Data$ccat
    H1 = Data$H1
    Z = Data$Z
    omega2 = Param[["omega2"]]
    
    theta = Param[["theta"]]
    K = Param[["K"]]
    
    tZ = cbind(Data[["Z"]], Param[["Knm"]])
    
    LD = (getB(nh, ccat))%*%(getC(nh, ccat)*c(H1%*%Param$delta))
    Linv = forwardsolve(LD,diag(nrow(LD)))
    Dinv = t(Linv)%*%Linv
    Phiinv = dbind(Dinv, K/theta) + t(tZ/omega2)%*%tZ
    
    list(LD=LD, Delta=LD%*%t(LD), Linv=Linv, Dinv=Dinv, Phiinv=Phiinv)
}

updateDelta = function(
YMat,
Data,
Param,
DEBUG=F, BFGS=T, Verbose=0){
    
    if(Verbose>0)print("updateDelta")
    
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    M = Data$MatrixDims$M
    
    H02 = dbind(Data$H0,matrix(rep(1,M),M))
    H1 = Data$H1
    H22 = dbind(Data$H2,matrix(1,1))
    ccat = Data$ccat
    nh = Data$nh
    Bs = Data$Bs
    Cs = Data$Cs
    
    sigma2 = Param[["sigma2"]]
    omega2 = Param[["omega2"]]
    Z = Data[["Z"]]
    W = Data[["W"]]
    K = Param[["K"]]
    Knm = Param[["Knm"]]
    zeta = Param[["zeta"]]
    theta = Param[["theta"]]
    delta = Param[["delta"]]
    
    tDinv = dbind(Param[["Dinv"]],K/theta)
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    Phiinv = tDinv + t(tZ/omega2)%*%tZ
    
    # matrix prep 1
    A = tDinv%*%(diag(K2+M)-Solve(Phiinv,tDinv)); A = (A+t(A))/2
    
    # update zeta
    ty = c(YMat$OneSinvYt-W%*%(Param[["Alpha"]]%*%(sigma2)))
    zeta = c(Solve(A[1:K2,1:K2]+diag(K2)/100,t(Z/omega2)%*%ty - t(Z/omega2)%*%tZ%*%Solve(Phiinv,t(tZ))%*%(ty/omega2)))/sum(1/sigma2)
    tW = cbind(W,Z%*%zeta)
    
    # matrix prep 2
    D = cbind(YMat$YtOinvZ,YMat$YtOinvKnm)
    E = t(tZ/omega2)%*%tW
    C0= t(D)-E%*%tA
    C = Solve(Phiinv, t(D)-E%*%tA)
    B = tDinv%*%C%*%(t(C)/sigma2)%*%tDinv/J
    G = C0%*%(t(C0)/sigma2)/J
    
    modeofig = 0.01
    iga = 0.1
    igb = modeofig*(1+iga)
    L = Data$MatrixDims$L # length(nh) ncol(H0)
    P = Data$MatrixDims$P # nrow(H2)
    grad = rep(0, P+1)
    hess = array(0, c(P+1,P+1))
    for(j in 1:(L+1)){
        Ajj = (-A+B)[H02[,j]>0,H02[,j]>0,drop=F]
        if(j==(L+1)){
            grad[H22[,j]>0] = sqrt(theta)^2*sum(diag(Solve(K,Ajj))) + 2*((-1-iga)+igb/theta)/J
        }else{
            grad[H22[,j]>0] = delta[j]^2*sum(diag(Ajj)) + 2*((-1-iga)+igb/delta[j]^2)/J
        }
    }
    # titsias
    a =  - sum(1/omega2)*Data$NGPs + sum(Solve(K,t(Knm))*t(Knm/omega2))
    grad[P+1] = grad[P+1] + a*sqrt(theta)^2

    # BFGS
    if(is.null(Param$Hinvdelta)){
        Hinv = -diag(P+1)*0.01 #solve(hess)
    }else{
        Hinvtemp = diag(P+1)-Param$ssdelta%*%t(grad-Param$graddelta)/sum(Param$ssdelta*(grad-Param$graddelta))
        Hinv = Hinvtemp%*%Param$Hinvdelta%*%t(Hinvtemp) + Param$ssdelta%*%t(Param$ssdelta)/sum(Param$ssdelta*(grad-Param$graddelta))
    }
    ss = -Hinv%*%grad

    LKinv = forwardsolve(t(chol(K)),diag(M))
    tLD = dbind(Param$LD, t(LKinv)*sqrt(theta))
    ZOinvZ = t(tZ/omega2)%*%tZ
    pelbo0 = - logDet(diag(K2+M)+t(tLD)%*%(ZOinvZ)%*%tLD) + sum(diag(Solve(Phiinv,G))) + a*theta + sum((-iga-1)*log(c(delta^2,theta))-igb/c(delta^2,theta))/J
    if(Verbose>1){cat("ss=");print(c(ss))}
    for(r in c(0:10,20,25)){
        ss1 = ss/2^r
        #deltatheta = c(delta,sqrt(theta)) + ss1
        deltatheta = exp(log(c(delta,sqrt(theta))) + ss1)
        delta1 = deltatheta[1:P]
        theta1 = deltatheta[P+1]^2
        
        if(sum(deltatheta>100)==0 && sum(deltatheta<1e-10)==0){
            LD = t(getB(nh, ccat))%*%(getC(nh, ccat)*c(H1%*%delta1))
            Delta = LD%*%t(LD)
            Linv = forwardsolve(LD,diag(nrow(LD)))
            Dinv = t(Linv)%*%Linv
            tDinv  = dbind(Dinv, K/theta1)
            Phiinv = tDinv+t(tZ/omega2)%*%tZ
            tLD = dbind(LD, t(LKinv)*sqrt(theta1))
        
            pelbo1 = - logDet(diag(K2+M)+t(tLD)%*%(ZOinvZ)%*%tLD) + sum(diag(Solve(Phiinv,G))) + a*theta1 + sum((-iga-1)*log(deltatheta^2)-igb/deltatheta^2)/J
        }else{
            pelbo1 = NA
        }
        if(Verbose>1){print(c(as.integer(r),pelbo0,pelbo1))}
        if(!is.na(pelbo1)){if(sign(pelbo1)*abs(pelbo1*1.0000001)>pelbo0){break}}
    }
    
    if(Verbose>0){
        cat("delta=");print(delta1)
        cat("theta=");print(sqrt(theta1))
    }
    
    list(zeta=zeta, delta=delta1, LD = LD, Delta= Delta, Linv=Linv, Dinv=Dinv, Phiinv=Phiinv, theta=theta1, ssdelta = ss1, graddelta = grad, Hinvdelta=Hinv)
}


