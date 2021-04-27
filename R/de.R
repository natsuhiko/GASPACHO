RunDE = function(Yt, gplvm, G=NULL, term="CC"){
    if(term=="CC"){
        if(is.null(G)){
            res = deCC(Yt, gplvm)
            rownames(res)=rownames(Yt)
        }else{
            bf=NULL
            pp=NULL
            for(i in 1:ncol(G)){
                tmp = deGbyCC(Yt, gplvm, G[,i])
                bf = cbind(bf, tmp)
            }
            res = data.frame(bf, emn(bf)[,-(ncol(bf)+1)])
            colnames(res)=paste(rep(c("bf","pp"),rep(ncol(G),2)),rep(colnames(G),2), sep="_")
            rownames(res)=rownames(Yt)
        }
    }
    res
}

deGbyCC = function(Yt, gplvm, g, deltaG=0.1){
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    N = gplvm$Data$MatrixDims$N
    M = gplvm$Data$MatrixDims$M
    
    if(length(gplvm$Param$A)>0){
        tA = rbind(gplvm$Param$A,1)
        tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
        tY = t(as.matrix(Yt)) - tW%*%tA
    }else{
        tY = t(as.matrix(Yt)) - c(gplvm$Data$Z%*%gplvm$Param$zeta)
    }
    
    K   = gplvm$Param$K
    Knm = gplvm$Param$Knm
    
    # kernel for reduced model
    Ks = updateKernelSE(gplvm$Param$Xi[[1]], reduceDim1(gplvm$Param$Ta[[1]]), gplvm$Param$rho[[1]])
    
    tZ = cbind(gplvm$Data$Z, gplvm$Param$Knm)
    tZ1 = cbind(gplvm$Data$Z, gplvm$Param$Knm, g*Ks$Knm)
    
    Dinv  = dbind(gplvm$Param$Dinv, gplvm$Param$K/theta)
    Dinv1 = dbind(gplvm$Param$Dinv, gplvm$Param$K/theta, Ks$K/theta/deltaG)
    
    Linv  = dbind(gplvm$Param$Linv, (chol(gplvm$Param$K))/sqrt(theta))
    Linv1 = dbind(gplvm$Param$Linv, (chol(gplvm$Param$K))/sqrt(theta), chol(Ks$K)/sqrt(theta*deltaG))
    
    L  = dbind(gplvm$Param$LD, backsolve(chol(gplvm$Param$K),diag(ncol(gplvm$Param$K))) * sqrt(theta))
    L1 = dbind(gplvm$Param$LD, backsolve(chol(gplvm$Param$K),diag(ncol(gplvm$Param$K))) * sqrt(theta),
                               backsolve(chol(Ks$K),diag(ncol(Ks$K))) * sqrt(theta*deltaG))
    
    Phiinv  = Dinv  + t(tZ /omega2)%*%tZ
    Phiinv1 = Dinv1 + t(tZ1/omega2)%*%tZ1
    
    ZtOinvy   = t(tZ/omega2)%*%tY
    Zt1Oinvy  = t(tZ1/omega2)%*%tY
    
    y2 = colSums(tY^2/omega2)
    s2 = (y2 - colSums(ZtOinvy *Solve(Phiinv,  ZtOinvy)))/(N)
    s21= (y2 - colSums(Zt1Oinvy*Solve(Phiinv1, Zt1Oinvy)))/(N)
    
    bf0 = -N/2*log(s2)  - logDet(diag(nrow(Dinv)) +t(L) %*%(t(tZ /omega2)%*%tZ )%*%L )/2 +
        sum(t(Knm/omega2)   *Solve(K,   t(Knm))   )*theta/2 - sum(1/omega2)*theta/2
    bf1 = -N/2*log(s21) - logDet(diag(nrow(Dinv1))+t(L1)%*%(t(tZ1/omega2)%*%tZ1)%*%L1)/2 +
        sum(t(Knm/omega2)   *Solve(K,   t(Knm))   )*theta/2 - sum(1/omega2)*theta/2 +
        sum(g^2*colSums(t(Ks$Knm/omega2)*Solve(Ks$K,t(Ks$Knm))))*theta/2 - sum(g^2/omega2)*theta/2
    
    
    bf1-bf0
}

deCC = function(Yt, gplvm){
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    N = gplvm$Data$MatrixDims$N
    M = gplvm$Data$MatrixDims$M
    
    if(length(gplvm$Param$A)>0){
        tA = rbind(gplvm$Param$A,1)
        tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
        tY = t(as.matrix(Yt)) - tW%*%tA
    }else{
        tY = t(as.matrix(Yt)) - c(gplvm$Data$Z%*%gplvm$Param$zeta)
    }
    
    K   = gplvm$Param$K
    Knm = gplvm$Param$Knm
    
    # kernel for reduced model
    Ks = updateKernelSE(gplvm$Param$Xi[[2]], gplvm$Param$Ta[[2]], gplvm$Param$rho[[2]])
    
    tZ0= cbind(gplvm$Data$Z, Ks$Knm)
    tZ = cbind(gplvm$Data$Z, gplvm$Param$Knm)
    
    Dinv0 = dbind(gplvm$Param$Dinv, Ks$K/theta)
    Dinv  = dbind(gplvm$Param$Dinv, gplvm$Param$K/theta)
    
    Linv0 = dbind(gplvm$Param$Linv, (chol(Ks$K))/sqrt(theta))
    Linv  = dbind(gplvm$Param$Linv, (chol(gplvm$Param$K))/sqrt(theta))
    
    L0 = dbind(gplvm$Param$LD, backsolve(chol(Ks$K),diag(ncol(Ks$K))) * sqrt(theta))
    L  = dbind(gplvm$Param$LD, backsolve(chol(gplvm$Param$K),diag(ncol(gplvm$Param$K))) * sqrt(theta))
    
    Phiinv0 = Dinv0 + t(tZ0/omega2)%*%tZ0
    Phiinv  = Dinv  + t(tZ /omega2)%*%tZ
    
    Zt0Oinvy = t(tZ0/omega2)%*%tY
    ZtOinvy  = t(tZ/omega2)%*%tY
    
    y2 = colSums(tY^2/omega2)
    s20= (y2 - colSums(Zt0Oinvy*Solve(Phiinv0, Zt0Oinvy)))/(N)
    s2 = (y2 - colSums(ZtOinvy *Solve(Phiinv,  ZtOinvy)))/(N)
    
    bf1 = -N/2*log(s2)  - logDet(diag(nrow(Dinv)) +t(L) %*%(t(tZ /omega2)%*%tZ )%*%L )/2 + sum(t(Knm/omega2)   *Solve(K,   t(Knm))   )*theta/2 - sum(1/omega2)*theta/2
    bf0 = -N/2*log(s20) - logDet(diag(nrow(Dinv0))+t(L0)%*%(t(tZ0/omega2)%*%tZ0)%*%L0)/2 + sum(t(Ks$Knm/omega2)*Solve(Ks$K,t(Ks$Knm)))*theta/2 - sum(1/omega2)*theta/2
    
    data.frame(bf=bf1-bf0, pp=emn(matrix(bf1-bf0,length(bf1)))[,1])
}
