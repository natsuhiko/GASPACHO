updateDxCSE <-
function(
Yt,
gplvm,
donorID="donor",
d2dxc0=0.01,
DEBUG=F, BFGS=T, Verbose=0){
    
    library(Matrix)
    
    if(Verbose>0)print("updateDxC")
    
    Data = gplvm$Data
    Param= gplvm$Param
    
    J = Data$MatrixDims$J
    K2= Data$MatrixDims$K2
    
    H1 = Data$H1
    ccat = Data$ccat
    nh = Data$nh
    Bs = Data$Bs
    Cs = Data$Cs
    
    Nd = nh[donorID];
    Zd = Data$Z[,Data$H0[,names(nh)==donorID]==1]
    
    sigma2 = getNewSigma(Yt,gplvm$Data,gplvm$Param) # Param[["sigma2"]]
    omega2 = Param[["omega2"]]
    Z = Data[["Z"]]
    W = Data[["W"]]
    K = Param[["K"]]
    Knm = Param[["Knm"]]
    zeta = Param[["zeta"]]
    theta = Param[["theta"]]
    delta = Param[["delta"]]
    
    # only immune
    K1 = updateKernelSE(Param$Xi[[2]][,-c(1,3,5)],reduceDim1(Param$Ta[[2]][,-c(1,3,5)],4),Param$rho[[2]][-c(1,3,5)])
    
    tDinv = dbind(Param[["Dinv"]],K/theta,Solve(K1$K)/theta/d2dxc0)
    tA = rbind(Param[["Alpha"]][,1:nrow(Yt)],1)
    tZ = cbind(Z,Knm,Zd[,1]*(K1$Knm%*%Solve(K1$K)))
    nk=nrow(Param[["Dinv"]])+nrow(K)
    Kd1=K1$K*theta
    Kdinv=Kd1 # diag(K/th,...,K/th)^-1
    for(i in 2:Nd){
        Kdinv = dbind(Kdinv, Kd1)
        tDinv = dbind(tDinv, Solve(K1$K)/theta/d2dxc0)
        tZ = cbind(tZ, Zd[,i]*(K1$Knm%*%Solve(K1$K)))
    }
    nk=c(nk,nrow(tDinv))
    #tDinv = drop0(Matrix(tDinv))
    #tZ    = drop0(Matrix(tZ))
    #ZtOinvZ = drop0(Matrix(t(tZ/omega2)%*%tZ))
    ZtOinvZ = t(tZ/omega2)%*%tZ
    Phiinv = tDinv + ZtOinvZ
    tW = cbind(W,Z%*%zeta)
    
    # matrix prep 1
    A = tDinv%*%(diag(nrow(Phiinv))-Solve(Phiinv,tDinv)); A = (A+t(A))/2
    
    # matrix prep 2
    YtOinvZ = as.matrix(Yt%*%(Z/omega2))
    ty2Oinv = as.numeric(as.matrix(Yt^2)%*%(1/omega2)) - 2*colSums(t(YtOinvZ%*%zeta)*tA) + colSums(tA*((t(tW/omega2)%*%tW)%*%tA))
    D = as.matrix(Yt%*%(tZ/omega2))
    E = as.matrix(t(tZ/omega2)%*%tW)
    C0= t(D)-E%*%tA # tZtOinvtY
    C = Solve(Phiinv, C0)
    B = tDinv%*%C%*%(t(C)/sigma2)%*%tDinv/J
    G = C0%*%(t(C0)/sigma2)/J
   
    d2dxc = d2dxc0
    ## titsias
    tpenalty =  - sum(1/omega2) + sum(Solve(K1$K,t(K1$Knm))*t(K1$Knm/omega2))
    for(itr in 1:100){
        
        # grad & hess
        grad = sum( (-A+B)[(nk[1]+1):nk[2],(nk[1]+1):nk[2]]*Kdinv ) + tpenalty*theta
        hess = A[(nk[1]+1):nk[2],(nk[1]+1):nk[2]]%*%Kdinv
        hess = -sum(hess*t(hess))
        grad = grad*d2dxc
        hess = hess*d2dxc^2 + grad
        
        # parameter update
        
        d2dxc = exp(log(d2dxc) - grad/hess)
        print(d2dxc)
        
        # matrix update
        tDinv = dbind(Param[["Dinv"]],K/theta)
        for(i in 1:Nd){
            tDinv = dbind(tDinv, Solve(K1$K)/theta/d2dxc)
        }
        #tDinv = drop0(Matrix(tDinv))
        Phiinv = tDinv + ZtOinvZ # t(tZ/omega2)%*%tZ
        A = tDinv%*%(diag(nrow(Phiinv))-Solve(Phiinv,tDinv)); A = (A+t(A))/2
        C = Solve(Phiinv, C0)

        sigma2 = (ty2Oinv - colSums(C0*C)+1)/(gplvm$Data$MatrixDims$N+1)

        B = tDinv%*%C%*%(t(C)/sigma2)%*%tDinv/J
        G = C0%*%(t(C0)/sigma2)/J
    }
    gplvm$Param$Phiinvdxc = Phiinv
    gplvm$Param$d2dxc = d2dxc
    gplvm$Param$Zdxc = tZ
    gplvm$Param$Dinvdxc = tDinv
    gplvm$Param$sigma2DxC = sigma2
    invisible(gplvm)
    
}
