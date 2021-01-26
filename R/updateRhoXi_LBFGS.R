updateRhoXi_LBFGS = function(
Yt,
YMat,
Data,
Param,
DEBUG=F,
Verbose=0){
    
    Diag=function(x){if(length(x)==1){matrix(x,1,1)}else{diag(x)}}
    if(Verbose>0)print("update Xi & Rho with LBFGS")
    
    UKP = Param$UpdateKernelParams
    
    J = Data$MatrixDims$J
    N = Data$MatrixDims$N
    K2= Data$MatrixDims$K2
    M = Data$MatrixDims$M
    Q = Data$MatrixDims$Q
    
    flagZK = rep(0:1,c(K2,sum(M)))
    
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
    Theta = rep(Param[["theta"]],M)
    
    tA = rbind(Param[["Alpha"]],1)
    tZ = cbind(Z,Knm)
    tW = cbind(W,Z%*%zeta)
    Phiinv = dbind(Param[["Dinv"]],K/Theta) + t(tZ/omega2)%*%tZ
    
    B = Solve(Phiinv, t(tZ))
    D = cbind(YMat$YtOinvZ,YMat$YtOinvKnm)
    E = t(tZ/omega2)%*%tW
    FF= -tW + tZ%*%Solve(Phiinv, E)
    G = t(t(t(D)-E%*%tA)/sigma2)
    H = rbind(YMat$ASinvYt,YMat$OneSinvYt)
    
    KinvKmn = Solve(K,t(Knm))
    
    # matrix prep for dKnm and dK
    C = Solve(Phiinv, -t(tZ) + (as.matrix(t(D/sigma2)%*%Yt) - E%*%H - (G%*%D)%*%B + (G%*%t(tA))%*%t(FF))/J)[flagZK==1,]
    C = C + Theta*KinvKmn # titsias
    C = t(C)/omega2
    
    # NOT SYMMETRIC!!!!!!!
    AA= Solve(K) - Solve(Phiinv)[flagZK==1,flagZK==1]/Theta
    AA= AA - Solve(Phiinv, t(Solve(Phiinv, G%*%t(t(D)-E%*%tA))))[flagZK==1,flagZK==1]/J/Theta
    AA= AA - KinvKmn%*%(t(KinvKmn)/omega2)*Theta # titsias
    
    gradXi = list()
    gradTa = list()
    gradRho= list()
    tC = (C*Knm)
    tAA = (AA*K)
    for(i in 1:length(Ta)){
        if(i==1){ # Periodic
            XT = outer(c(Xi[[i]]),c(Ta[[i]]),"-")
            TT = outer(c(Ta[[i]]),c(Ta[[i]]),"-")
            dqKX = -2*cos(XT)*sin(XT)/rho[[i]]^2
            dqKT = -2*cos(TT)*sin(TT)/rho[[i]]^2
            
            gradXi = c(gradXi, list(matrix(rowSums(tC*dqKX),N)))
            
            dqK   = sin(outer(c(Ta[[i]]),c(Ta[[i]]),"-"))^2/rho[[i]]^2. # K  *
            dqKnm = sin(outer(c(Xi[[i]]),c(Ta[[i]]),"-"))^2/rho[[i]]^2. # Knm*
            gradRho = c(gradRho, list(sum(tC*dqKnm) + sum(tAA*dqK)/2))
            
            gradTa = c(gradTa, list(matrix(-colSums(tC*dqKX) + rowSums(2*tAA*dqKT)/2, M)))
        }else{ # SE
            gradXi = c(gradXi, list((- rowSums(tC)*Xi[[i]] + tC%*%Ta[[i]])%*%Diag(1/rho[[i]]^2) - (Xi[[i]]-Param$PriorXi2$mu)/Param$PriorXi2$var/J))
            
            gradRho = c(gradRho, rep(0,Q[i]))
            for(q in seq(Q[i])){
                dqK   = outer(Ta[[i]][,q],Ta[[i]][,q],"-")^2/rho[[i]][q]^2/2. # K  *
                dqKnm = outer(Xi[[i]][,q],Ta[[i]][,q],"-")^2/rho[[i]][q]^2/2. # Knm*
                gradRho[[i]][q] = sum(tC*dqKnm) + sum(tAA*dqK)/2
            }
            gradTa = c(gradTa, list((t(tC)%*%Xi[[i]] - colSums(tC)*Ta[[i]])%*%Diag(1/rho[[i]]^2) - Ta[[i]]/J + (tAA-diag(rowSums(tAA)))%*%Ta[[i]]%*%Diag(1/rho[[i]]^2)))
        }
    }
    
    ValNew = NULL
    GradNew= NULL
    strend = 0
    for(i in 1:length(Ta)){
        if(UKP[i]==7){
            ValNew  = c(ValNew,  c(rbind(2*log(rho[[i]]), Xi[[i]],     Ta[[i]])))
            GradNew = c(GradNew, c(rbind(gradRho[[i]],    gradXi[[i]], gradTa[[i]])))
        }else if(UKP[i]==6){
            ValNew  = c(ValNew,  c(rbind(2*log(rho[[i]]), Xi[[i]])))
            GradNew = c(GradNew, c(rbind(gradRho[[i]],    gradXi[[i]])))
        }else if(UKP[i]==5){
            ValNew  = c(ValNew,  c(rbind(2*log(rho[[i]]), Ta[[i]])))
            GradNew = c(GradNew, c(rbind(gradRho[[i]],    gradTa[[i]])))
        }else if(UKP[i]==1){
            ValNew  = c(ValNew,  c(rbind(Ta[[i]])))
            GradNew = c(GradNew, c(rbind(gradTa[[i]])))
        }
        strend = c(strend,length(ValNew))
    }
    ValNew  = matrix(ValNew, length(ValNew))
    GradNew = matrix(GradNew,length(GradNew))
    if(is.null(Param$GradRhoXi)){
        print("init LBFGS")
        ss = -c(GradNew)/5000
    }else{
        ValNew = cbind(ValNew,Param$ValRhoXi)
        GradNew = cbind(GradNew,Param$GradRhoXi)
        ss = -LBFGS(ValNew, GradNew)
    }
    ss2 = list()
    for(i in 1:length(Ta)){
        if(UKP[i]==7){
            ss2 = c(ss2, list(matrix(ss[(strend[i]+1):strend[i+1]],1+N+M)))
        }else if(UKP[i]==6){
            ss2 = c(ss2, list(matrix(ss[(strend[i]+1):strend[i+1]],1+N)))
        }else if(UKP[i]==5){
            ss2 = c(ss2, list(matrix(ss[(strend[i]+1):strend[i+1]],1+M)))
        }else if(UKP[i]==1){
            ss2 = c(ss2, list(matrix(ss[(strend[i]+1):strend[i+1]],M)))
        }
    }
    ss = ss2
    
    tLD = dbind(Param$LD, backsolve((chol(K)),diag(sum(M)))*sqrt(Theta))
    pelbonow = sum(c(- logDet(diag(K2+sum(M))+t(tLD)%*%(t(tZ/omega2)%*%tZ)%*%tLD), + sum(diag(Solve(Phiinv,G%*%t(t(D)-E%*%tA)/J))),
                + sum(t(Knm/omega2)*Solve(K/Theta,t(Knm))), - sum((Xi[[2]]-Param$PriorXi2$mu)^2/Param$PriorXi2$var)/J - sum(Ta[[2]]^2)/J )) #- sum(1/omega2)*theta 
    pelbo0=pelbo1=-(10^10)
    rho0 = rho1 = rho
    Xi0  = Xi1  = Xi
    Ta0  = Ta1  = Ta
    if(Verbose>1){cat("ss=");print(c(ss[[1]][1,],ss[[2]][1,]))}
    if(0){ # reverse
        for(r in c(50,20,10:0,-1)){
            rho1 = list(exp(log(rho[[1]]) - (ss[[1]][1]/2.)/2^r),         exp(log(rho[[2]]) - (ss[[2]][1,]/2.)/2^r))
            Xi1  = list(Xi[[1]]           -  ss[[1]][2:(1+N)]/2^r,        Xi[[2]]           - ss[[2]][2:(1+N),]/2^r)
            Ta1  = list(Ta[[1]]           -  ss[[1]][(N+2):(1+N+M[1])]/2^r, Ta[[2]]           - ss[[2]][(N+2):(1+N+M[2]),]/2^r)
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
    }else{ # forward
        for(r in c(0:10,20,50)){
            for(i in 1:length(Ta)){
                if(UKP[i]==7){
                    rho1[[i]] = exp(log(rho[[i]]) - (ss[[i]][1,]/2.)/2^r)
                    Xi1[[i]]  = Xi[[i]]           -  ss[[i]][2:(1+N),]/2^r
                    Ta1[[i]]  = Ta[[i]]           -  ss[[i]][(N+2):(1+N+M),]/2^r
                }else if(UKP[i]==6){
                    rho1[[i]] = exp(log(rho[[i]]) - (ss[[i]][1,]/2.)/2^r)
                    Xi1[[i]]  = Xi[[i]]           -  ss[[i]][2:(1+N),]/2^r
                }else if(UKP[i]==5){
                    rho1[[i]] = exp(log(rho[[i]]) - (ss[[i]][1,]/2.)/2^r)
                    Ta1[[i]]  = Ta[[i]]           -  ss[[i]][(2):(1+M),]/2^r
                }else if(UKP[i]==1){
                    Ta1[[i]]  = Ta[[i]]           -  ss[[i]][(1):(M),]/2^r
                }
            }
            if(Verbose>1){cat("rho1=");print(unlist(rho1))}
            pelbo1 = getELBOrho(Yt, YMat, Data, Param, Xi1, Ta1, rho1)
            if(Verbose>1){print(c(as.integer(r),pelbo0,pelbo1))}
            if(is.na(pelbo1)){break}
            if(pelbonow<pelbo1*1.0000001){
                rho0 = rho1
                Xi0  = Xi1
                Ta0  = Ta1
                break
            }
        }
    }
    if(Verbose>0){cat("rho=");print(unlist(rho0))}
    list(rho=rho0, Xi=Xi0, Ta=Ta0, ValRhoXi = ValNew[,seq(min(25,ncol(ValNew))),drop=F], GradRhoXi = GradNew[,seq(min(25,ncol(ValNew))),drop=F])
}



