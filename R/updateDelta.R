initDelta = function(Data, Param){
    print("initDelta")
    nh = Data$nh
    ccat = Data$ccat
    H1 = Data$H1
    Z = Data$Z
    omega2 = Param[["omega2"]]
    
    if(Data$MatrixDims$M>0){
        theta = Param[["theta"]]
        Theta = rep(theta, Data$MatrixDims$M)
        K = Param[["K"]]
        
        tZ = cbind(Data[["Z"]], Param[["Knm"]])
        
        LD = (getB(nh, ccat))%*%(getC(nh, ccat)*c(H1%*%Param$delta))
        Linv = forwardsolve(LD,diag(nrow(LD)))
        Dinv = t(Linv)%*%Linv
        Phiinv = dbind(Dinv, K/Theta) + t(tZ/omega2)%*%tZ
    }else{
        tZ = cbind(Data[["Z"]])
        
        LD = (getB(nh, ccat))%*%(getC(nh, ccat)*c(H1%*%Param$delta))
        Linv = forwardsolve(LD,diag(nrow(LD)))
        Dinv = t(Linv)%*%Linv
        Phiinv = Dinv + t(tZ/omega2)%*%tZ
    }
    
    list(LD=LD, Delta=LD%*%t(LD), Linv=Linv, Dinv=Dinv, Phiinv=Phiinv)
}


updateDelta = function(
YMat,
Data,
Param,
DEBUG=F, BFGS=T, Verbose=0){
    if(Data$Kernel=="SE"){
        updateDeltaSE(YMat, Data, Param, DEBUG, BFGS, Verbose)
    }else{
        updateDeltaLinear(YMat, Data, Param, DEBUG, BFGS, Verbose)
    }
}

updateDeltaLinear = function(
YMat,
Data,
Param,
DEBUG=F, BFGS=T, Verbose=0){
    
    if(Verbose>0)print("updateDelta")
    
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    M = Data$MatrixDims$M
    
    H0K = diag(length(M))[rep(seq(length(M)),M),]
    H02 = dbind(Data$H0,H0K)
    H1 = Data$H1
    H22 = dbind(Data$H2,diag(length(M)))
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
    Theta = rep(theta, M)
    delta = Param[["delta"]]
    
    tDinv = dbind(Param[["Dinv"]],K/Theta)
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    Phiinv = tDinv + t(tZ/omega2)%*%tZ
    
    # matrix prep 1
    A = tDinv%*%(diag(K2+sum(M))-Solve(Phiinv,tDinv)); A = (A+t(A))/2
    
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
    grad = colSums(((-A+B)*dbind(diag(K2),Solve(K)))%*%H02)*c(delta^2,theta) + 2*((-1-iga)+igb/c(delta^2,theta))/J
    grad = grad[1:P]
    
    # BFGS
    if(is.null(Param$Hinvdelta)){
        Hinv = -diag(P)*0.01 #solve(hess)
    }else{
        Hinvtemp = diag(P)-Param$ssdelta%*%t(grad-Param$graddelta)/sum(Param$ssdelta*(grad-Param$graddelta))
        Hinv = Hinvtemp%*%Param$Hinvdelta%*%t(Hinvtemp) + Param$ssdelta%*%t(Param$ssdelta)/sum(Param$ssdelta*(grad-Param$graddelta))
    }
    ss = -Hinv%*%grad
    
    if(length(K)>0){
        RKinv = backsolve((chol(K)),diag(sum(M))) # Kinv = RKinv%*%LKinv & K=L%*%R => I = LKinv%*% K %*%RKinv !!!!!!!!!
    }else{
        RKinv = array(0,c(0,0))
    }
    tLD = dbind(Param$LD, RKinv*sqrt(Theta))
    ZOinvZ = t(tZ/omega2)%*%tZ
    ptlb0 = - logDet(diag(K2+sum(M))+t(tLD)%*%(ZOinvZ)%*%tLD) + sum(diag(Solve(Phiinv,G))) + sum((-iga-1)*log(c(delta^2,theta))-igb/c(delta^2,theta))/J
    if(Verbose>1){cat("ss=");print(c(ss))}
    for(r in c(0:10,20,25)){
        ss1 = ss/2^r
        #deltatheta = c(delta,sqrt(theta)) + ss1
        deltatheta = exp(log(delta) + ss1)
        delta1 = deltatheta
        theta1 = 1
        Theta1 = rep(1, M)
        
        if(sum(deltatheta>100)==0 && sum(deltatheta<1e-10)==0){
            LD = t(getB(nh, ccat))%*%(getC(nh, ccat)*c(H1%*%delta1))
            Delta = LD%*%t(LD)
            Linv = forwardsolve(LD,diag(nrow(LD)))
            Dinv = t(Linv)%*%Linv
            tDinv  = dbind(Dinv, K/Theta1)
            Phiinv = tDinv+t(tZ/omega2)%*%tZ
            tLD = dbind(LD, RKinv*sqrt(Theta1))
            
            ptlb1 = - logDet(diag(K2+sum(M))+t(tLD)%*%(ZOinvZ)%*%tLD) + sum(diag(Solve(Phiinv,G))) + sum((-iga-1)*log(deltatheta^2)-igb/deltatheta^2)/J
        }else{
            ptlb1 = NA
        }
        if(Verbose>1){print(c(as.integer(r),ptlb0,ptlb1))}
        if(!is.na(ptlb1)){if(sign(ptlb1)*abs(ptlb1*1.0000001)>ptlb0){break}}
    }
    
    if(Verbose>0){
        cat("delta=");print(delta1)
        cat("theta=");print(sqrt(theta1))
    }
    
    list(zeta=zeta, delta=delta1, LD = LD, Delta= Delta, Linv=Linv, Dinv=Dinv, Phiinv=Phiinv, theta=theta1, ssdelta = ss1, graddelta = grad, Hinvdelta=Hinv)
}






updateDeltaSE = function(
YMat,
Data,
Param,
DEBUG=F, BFGS=T, Verbose=0){
    
    if(Verbose>0)print("updateDelta")
    
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    M = Data$MatrixDims$M
    
    H0K = diag(length(M))[rep(seq(length(M)),M),]
    H02 = dbind(Data$H0,H0K)
    H1 = Data$H1
    H22 = dbind(Data$H2,diag(length(M)))
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
    Theta = rep(theta, M)
    delta = Param[["delta"]]
    
    tDinv = dbind(Param[["Dinv"]],K/Theta)
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    Phiinv = tDinv + t(tZ/omega2)%*%tZ
    
    # matrix prep 1
    A = tDinv%*%(diag(K2+sum(M))-Solve(Phiinv,tDinv)); A = (A+t(A))/2
    
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
    grad = colSums(((-A+B)*dbind(diag(K2),Solve(K)))%*%H02)*c(delta^2,theta) + 2*((-1-iga)+igb/c(delta^2,theta))/J
    
    # titsias
    a =  - sum(1/omega2) + rowSums(t(H0K)%*%(Solve(K,t(Knm))*t(Knm/omega2)))
    grad = grad + c(rep(0,P), a*theta)

    # BFGS
    if(is.null(Param$Hinvdelta)){
        Hinv = -diag(P+length(M))*0.01 #solve(hess)
    }else{
        Hinvtemp = diag(P+length(M))-Param$ssdelta%*%t(grad-Param$graddelta)/sum(Param$ssdelta*(grad-Param$graddelta))
        Hinv = Hinvtemp%*%Param$Hinvdelta%*%t(Hinvtemp) + Param$ssdelta%*%t(Param$ssdelta)/sum(Param$ssdelta*(grad-Param$graddelta))
    }
    ss = -Hinv%*%grad

    RKinv = backsolve((chol(K)),diag(sum(M))) # Kinv = RKinv%*%LKinv & K=L%*%R => I = LKinv%*% K %*%RKinv !!!!!!!!!
    tLD = dbind(Param$LD, RKinv*sqrt(Theta))
    ZOinvZ = t(tZ/omega2)%*%tZ
    ptlb0 = - logDet(diag(K2+sum(M))+t(tLD)%*%(ZOinvZ)%*%tLD) + sum(diag(Solve(Phiinv,G))) + sum(a*theta) + sum((-iga-1)*log(c(delta^2,theta))-igb/c(delta^2,theta))/J
    if(Verbose>1){cat("ss=");print(c(ss))}
    for(r in c(0:10,20,25)){
        ss1 = ss/2^r
        #deltatheta = c(delta,sqrt(theta)) + ss1
        deltatheta = exp(log(c(delta,sqrt(theta))) + ss1)
        delta1 = deltatheta[1:P]
        theta1 = deltatheta[-seq(P)]^2
        Theta1 = rep(theta1, M)
        
        if(sum(deltatheta>100)==0 && sum(deltatheta<1e-10)==0){
            LD = t(getB(nh, ccat))%*%(getC(nh, ccat)*c(H1%*%delta1))
            Delta = LD%*%t(LD)
            Linv = forwardsolve(LD,diag(nrow(LD)))
            Dinv = t(Linv)%*%Linv
            tDinv  = dbind(Dinv, K/Theta1)
            Phiinv = tDinv+t(tZ/omega2)%*%tZ
            tLD = dbind(LD, RKinv*sqrt(Theta1))
        
            ptlb1 = - logDet(diag(K2+sum(M))+t(tLD)%*%(ZOinvZ)%*%tLD) + sum(diag(Solve(Phiinv,G))) + sum(a*theta1) + sum((-iga-1)*log(deltatheta^2)-igb/deltatheta^2)/J
        }else{
            ptlb1 = NA
        }
        if(Verbose>1){print(c(as.integer(r),ptlb0,ptlb1))}
        if(!is.na(ptlb1)){if(sign(ptlb1)*abs(ptlb1*1.0000001)>ptlb0){break}}
    }
    
    if(Verbose>0){
        cat("delta=");print(delta1)
        cat("theta=");print(sqrt(theta1))
    }
    
    list(zeta=zeta, delta=delta1, LD = LD, Delta= Delta, Linv=Linv, Dinv=Dinv, Phiinv=Phiinv, theta=theta1, ssdelta = ss1, graddelta = grad, Hinvdelta=Hinv)
}


