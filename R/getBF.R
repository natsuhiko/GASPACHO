

getBF = function(
yj,
gplvm,
G, # N_L x D non-scaled
did, # donor ID
d2G=0.001
){

    af = rowMeans(G)/2
    G = (G-2*af)/sqrt(2*af*(1-af))
    G = t(G)
    nu = d2G # hyperparameter of d2G s.t. sqrt(d2G) ~ N(0,nu)
    
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    N = gplvm$Data$MatrixDims$N
    L = ncol(G)

    tA = rbind(gplvm$Param$Alpha,1)[,1]
    tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
    yj = yj - tW%*%tA

    # kernel for reduced model
    Ks1 = updateKernelSE(gplvm$Param$Xi[[2]][,-c(1,3,5)],reduceDim1(gplvm$Param$Ta[[2]][,-c(1,3,5)],4),gplvm$Param$rho[[2]][-c(1,3,5)]) # immune
    M1 = nrow(Ks1$K)+1
    Ks2 = updateKernel(gplvm$Data, gplvm$Param)

    tZ    = cbind(gplvm$Data$Z,     Ks2$Knm)
    tDinv = dbind(gplvm$Param$Dinv, Ks2$K/theta)
    for(i in 1:length(unique(did))){
        tDinv = dbind(tDinv, Ks1$K/theta/gplvm$Param$d2dxc)
        tZ = cbind(tZ, Ks1$Knm*as.numeric(did==i))
    }
    tPhiinv  = tDinv  + t(tZ/omega2)%*%tZ
    Phi = Solve(tPhiinv)

    X = cbind(Ks1$Knm,1)

    ZtOinvy = t(tZ/omega2)%*%yj
    PhiZtOinvy = Solve(tPhiinv, ZtOinvy)
    XtVinvX=array(0,c(M1*M1,L))
    XtVinvy=array(0,c(M1,L))
    for(l in 1:L){
        if(l%%100==0)print(l)
        ZtOinvX = t(tZ/omega2)%*%(X*G[did,l])
        XtVinvX[,l] = c(t(X/omega2*G[did,l])%*%(X*G[did,l]) - t(ZtOinvX) %*% Phi %*% ZtOinvX) # M1 x M1 x L
        XtVinvy[,l] = c(t(X/omega2*G[did,l])%*%yj           - t(ZtOinvX) %*% PhiZtOinvy) # M1 x L
    }

    y2 = sum(yj^2/omega2)
    s2 = y2 - sum(ZtOinvy *Solve(tPhiinv,  ZtOinvy)) # yVinvy

    titsias1 = theta*(colSums(G[did,]^2/omega2)-colSums(G[did,]^2*colSums(t(Ks1$Knm/omega2)*Solve(Ks1$K,t(Ks1$Knm))))) # tr{Oinv theta Gl tKNN^(1) Gl}

    DeltaGprime = dbind(Solve(Ks1$K)*theta,matrix(1,1))
    DeltaGinv   = dbind(Ks1$K/theta,matrix(1,1))/d2G

    ## BCs
    VBc=NULL
    for(l in 1:L){
        VBc = cbind(VBc,
            c(Solve( DeltaGinv + XtVinvX[,l] ))
        )
    }
    Bc = XtVinvy
    for(l in 1:L){
        Bc[,l] = matrix(VBc[,l],M1)%*%XtVinvy[,l]
    }
    
    s2j = bfs = rep(NA,L)
    LOWERBOUND=F
    if(LOWERBOUND){
        s2j = s2 - 2*colSums(XtVinvy*Bc) + colSums(XtVinvX*Bc[rep(1:M1,M1),]*Bc[rep(1:M1,rep(M1,M1)),])
        bfs = (-N*log(s2j/N) + N*log(s2/N) - titsias1*d2G - colSums(XtVinvX*VBc))/2
    }else{ # exact
        s2j = s2 - colSums(XtVinvy*Bc)
        bfs = -N*log(s2j/N)/2 + N*log(s2/N)/2 - apply(XtVinvX, 2, function(hj){logDet(matrix(hj,M1)+DeltaGinv)})/2 - logDet(DeltaGprime*d2G)/2 - titsias1*d2G
    }
    
    SCORETEST=T
    pval=rep(NA,L)
    if(SCORETEST){
        scr = colSums((DeltaGprime%*%XtVinvy)*XtVinvy)/s2*N
        tmp = eigen((DeltaGprime+t(DeltaGprime))/2,symmetric=T)
        P = tmp[[2]]%*%diag(sqrt(tmp[[1]]))
        library(CompQuadForm)
        for(l in 1:L){
            hk=matrix(XtVinvX[,l],M1);
            lambda=eigen(t(P)%*%hk%*%P, symmetric=T)[[1]];
            pval[l] = davies(scr[l],lambda)$Qq
        }
    }
    
    # penalty
    finfo = titsias1 + 1/nu
    for(l in 1:L){
        tmp = matrix(XtVinvX[,l],M1)%*%DeltaGprime
        finfo[l] = finfo[l] + 2*sum(tmp*t(tmp))*d2G
    }
    cbind(bfs - d2G/nu/2 - log(nu)/2 - log(finfo)/2, pval)
}
