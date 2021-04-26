reduceDimAll <-
function(Ta){
    Tall = cbind(Ta[[1]],Ta[[2]])
    grp = cutree(hclust(dist(Tall)),h=2)
    w = as.numeric(1/table(grp)[grp])
    Tall = t(t(Tall)%*%(diag(length(unique(grp)))[grp,]*w))
    Ta=list(Tall[,1,drop=F],Tall[,-1])
    Ta
}

reduceDim1 <-
    function(Ta){
        Tall = Ta;#cbind(Ta[[1]],Ta[[2]])
        grp = cutree(hclust(dist(Tall)),h=2)
        w = as.numeric(1/table(grp)[grp])
        Tall = t(t(Tall)%*%(diag(length(unique(grp)))[grp,]*w))
        Ta=Tall #Tall[,2:ncol(Tall),drop=F]
        Ta
    }

getSpatialDE = function(
Yt,
gplvm,
W0,
B0=NULL,
XB=NULL,
ItrMax=100
){
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    N = gplvm$Data$MatrixDims$N
    M = gplvm$Data$MatrixDims$M
    
    tA = rbind(gplvm$Param$A,1)
    tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
    tY = t(as.matrix(Yt)) - tW%*%tA
    
    K   = gplvm$Param$K
    Knm = gplvm$Param$Knm
    
    # kernel for reduced model
    Param1=gplvm$Param
    Param1$Xi[[2]] =Param1$Xi[[2]][,c(1,3,5),drop=F]
    Param1$Ta[[2]] =Param1$Ta[[2]][,c(1,3,5),drop=F]
    Param1$Ta = reduceDimAll(Param1$Ta)
    Param1$rho[[2]]=Param1$rho[[2]][c(1,3,5)]
    Ks = updateKernel(gplvm$Data, Param1)
    
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
    
    # immune clustering
    Ks2 = updateKernelSE(gplvm$Param$Xi[[2]][,-c(1,3,5)], reduceDim1(gplvm$Param$Ta[[2]][,-c(1,3,5)]), gplvm$Param$rho[[2]][-c(1,3,5)])
    #mvnprop = mvnMix2(tY%*%diag(1/sqrt(gplvm$Param$sigma2)), t(Solve(Ks2$K,t(Ks2$Knm))), gplvm$Param$omega2, tZ0, Phiinv0, Solve(Ks2$K), W0)
    X = t(Solve(Ks2$K,t(Ks2$Knm)))
    if(is.null(B0)){
        B0 = NULL
        for(k in 1:ncol(W0)){
            yy = rowMeans(tY[,W0[,k]>0.5])
            B0 = cbind(B0, coef(lm(yy~X-1)))
            gplot(U$layout,,X%*%B0[,k],sym=F,col=magma(100))
        }
    }
    if(!is.null(XB)){
        B0 = as.matrix(coef(lm(XB~X-1)))
        cat("dim(B0)=");print(dim(B0));print(dim(X))
    }
    #mvnprop = mvnMix2(tY%*%diag(1/sqrt(gplvm$Param$sigma2)), X, B0, gplvm$Param$omega2, tZ0, Dinv0, theta*Ks2$K, W0)
    mvnprop = mvnMix2(t(t(tY)/sqrt(gplvm$Param$sigma2)), X, B0, gplvm$Param$omega2, tZ0, Dinv0, theta*Ks2$K, W0, itrmax=ItrMax)
    
    list(cbind(bf1-bf0,emn(matrix(bf1-bf0,length(bf1)))), mvnprop)
}





getBeta = function(
y,
j,
gplvm,
g,
deltaG = 1.0,
rhoscale=1.0
){
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    Phiinv0 = gplvm$Param$Phiinv
    N = gplvm$Data$MatrixDims$N
    M = gplvm$Data$MatrixDims$M
    
    
    tA = rbind(gplvm$Param$A,1)[,j,drop=F]
    tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
    yt = y - tW%*%tA
    
    # kernel for g
    Ks = updateKernelSE(gplvm$Param$Xi[[2]][,-c(1,3,5)],reduceDim1(gplvm$Param$Ta[[2]][,-c(1,3,5)]),gplvm$Param$rho[[2]][-c(1,3,5)]*rhoscale)
    #Ks = updateKernelSE(gplvm$Param$Xi[[2]][,-1],gplvm$Param$Ta[[2]][,-1],gplvm$Param$rho[[2]][-1]*rhoscale)
    tZ0= cbind(gplvm$Data$Z, gplvm$Param$Knm)
    tZ = cbind(gplvm$Data$Z, gplvm$Param$Knm, g*Ks$Knm,scale(g))
    
    #tZ = cbind(gplvm$Data$Z, gplvm$Param$Knm, scale(g))
    
    Dinv0 = dbind(gplvm$Param$Dinv, gplvm$Param$K/gplvm$Param$theta)
    Dinv = dbind(gplvm$Param$Dinv, gplvm$Param$K/gplvm$Param$theta, Ks$K/gplvm$Param$theta/deltaG, 1/deltaG)
    
    Linv = dbind(gplvm$Param$Linv, (chol(gplvm$Param$K))/sqrt(gplvm$Param$theta), (chol(Ks$K))/sqrt(gplvm$Param$theta*deltaG), 1/sqrt(deltaG))
    Linv0 = dbind(gplvm$Param$Linv, (chol(gplvm$Param$K))/sqrt(gplvm$Param$theta))
    
    L  = dbind(gplvm$Param$LD, backsolve(chol(gplvm$Param$K),diag(ncol(gplvm$Param$K))) * sqrt(gplvm$Param$theta), backsolve(chol(Ks$K),diag(ncol(Ks$K)))*sqrt(gplvm$Param$theta*deltaG), 1/sqrt(deltaG))
    L0 = dbind(gplvm$Param$LD, backsolve(chol(gplvm$Param$K),diag(ncol(gplvm$Param$K))) * sqrt(gplvm$Param$theta))
    
    Phiinv = Dinv + t(tZ/omega2)%*%tZ
    
    ZtOinvy = c(t(tZ/omega2)%*%yt)
    Zt0Oinvy = c(t(tZ0/omega2)%*%yt)
    beta = Solve(Phiinv,ZtOinvy)
    s20= (sum(yt^2/omega2) - sum(Zt0Oinvy*Solve(Phiinv0,Zt0Oinvy)))/(N)
    s2 = (sum(yt^2/omega2) - sum(ZtOinvy*beta))/(N)
    
    bf = -N/2*log(s2) - logDet(diag(nrow(Dinv))+t(L)%*%(t(tZ/omega2)%*%tZ)%*%L)/2 - (-N/2*log(s20) - logDet(diag(nrow(Linv0))+t(L0)%*%(t(tZ0/omega2)%*%tZ0)%*%L0)/2) + sum(t(Ks$Knm/omega2)*Solve(Ks$K,t(Ks$Knm)))*theta*deltaG/2 - sum(1/omega2)*theta*deltaG/2
    
    betas = c(cbind(Ks$Knm,1/sd(g))%*%rev(rev(beta)[1:(ncol(Ks$Knm)+1)]))
    attr(betas,"bf") = bf
    betas
}


GetBFAll = function(
y,
j,
gplvm,
G,
deltaG = 1.0,
rhoscale=1.0,
LF=10
){
    bfs = NULL
    cat("|",sep="")
    for(i in 1:(((ncol(G)-1)%/%100)+1)){
        cat("=",sep="")
        bfs = c(bfs, getBFAll(y,j,gplvm,G[,((i-1)*100+1):min(i*100,ncol(G)),drop=F],deltaG=deltaG,rhoscale=rhoscale,LF=LF))
    }
    cat("|100%",sep="\n")
    bfs
}

getBFAll = function(
y,
j,
gplvm,
G,
deltaG = 1.0,
rhoscale=1.0,
LF=10
){
    N = gplvm$Data$MatrixDims$N
    NL = ncol(G)
    
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    Phiinv = gplvm$Param$Phiinv

    tA = rbind(gplvm$Param$A,1)[,j,drop=F]
    tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
    tZ = cbind(gplvm$Data$Z,gplvm$Param$Knm)
    
    yt = y - tW%*%tA
    
    # extra columns of Z
    reduceDim1 <-
    function(Ta){
        Tall = Ta;#cbind(Ta[[1]],Ta[[2]])
        grp = cutree(hclust(dist(Tall)),h=2)
        w = as.numeric(1/table(grp)[grp])
        Tall = t(t(Tall)%*%(diag(length(unique(grp)))[grp,]*w))
        Ta=Tall #Tall[,2:ncol(Tall),drop=F]
        Ta
    }
    
    Ks = updateKernelSE(gplvm$Param$Xi[[2]][,-c(1,3,5)],reduceDim1(gplvm$Param$Ta[[2]][,-c(1,3,5)]),gplvm$Param$rho[[2]][-c(1,3,5)]*rhoscale)
    #U = Ks$Knm
    U = t(Solve(Ks$K,t(Ks$Knm)))
    Gk = apply(G, 2, scale)
    K2 = 1
    if(LF>0){
        K2 = ncol(U)+1
        Gk = ( t(rep(1,NL))%x%cbind(U,1) )*( G%x%t(rep(1,K2)) ) # N x (K2 * NL)
        Gk[,1:NL*K2] = apply(G, 2, scale)
    }
    
    
    
    
    # temp vars
    ZtOinvy = t(tZ/omega2)%*%yt
    ZtOinvGk= t(tZ/omega2)%*%Gk # ncol(Z) x (K2 * NL)
    PhiZtOinvy = Solve(Phiinv,ZtOinvy)
    a = sum(yt^2/omega2)
    bk = matrix(c(t(yt/omega2)%*%Gk),K2) # K2 x NL
    Ck = NULL # t(Gk/omega2)%*%Gk  # K2^2 x NL
    for(i in 1:NL){
        Ck=cbind(Ck,c(t(Gk[,((i-1)*K2+1):(i*K2)]/omega2)%*%Gk[,((i-1)*K2+1):(i*K2)]))
    }
    d  = sum(PhiZtOinvy*ZtOinvy)
    ek = matrix(c(t(ZtOinvGk)%*%PhiZtOinvy),K2) # K2 x NL
    Fk = NULL # t(ZtOinvGk)%*%solve(A, ZtOinvGk) # K2^2 x NL
    for(i in 1:NL){
        Fk=cbind(Fk,c(t(ZtOinvGk[,((i-1)*K2+1):(i*K2)])%*%Solve(Phiinv,ZtOinvGk[,((i-1)*K2+1):(i*K2)])))
    }
    bfs=NULL
    for(deltaG in c(0.01,0.1,1,10)){
        if(LF>0){
            #DeltaG    = deltaG*dbind(theta*Solve(Ks$K),matrix(1,1))
            #DeltaGinv = dbind(Ks$K/theta,matrix(1,1)) / deltaG
            
            DeltaG    = deltaG*dbind(theta*Ks$K,matrix(1,1))
            DeltaGinv = dbind(Solve(Ks$K)/theta,matrix(1,1)) / deltaG
        }else{
            DeltaG = matrix(1,1,1)*deltaG
            DeltaGinv = matrix(1,1,1)/deltaG
        }
        Hk = Ck - Fk + c(DeltaGinv)  # K2^2 x NL
        Hkinv=NULL
        for(i in 1:NL){Hkinv=cbind(Hkinv, c(Solve(matrix(Hk[,i],K2))))}
        
        
        # new sigma2
        N = length(y)
        s2nul = (a-d)/N
        s2new = (a-d-colSums(Hkinv*((ek-bk)%x%rep(1,K2))*(rep(1,K2)%x%(ek-bk))))/N
        
        # BF
        #-N*log(s2nul)/2
        bfs = cbind(bfs,
            -N*log(s2new)/2 - apply(Hk, 2, function(Hk1){logDet(matrix(Hk1,K2))})/2 - logDet(DeltaG)/2 + N*log(s2nul)/2 +
                sum(t(Ks$Knm/omega2)*Solve(Ks$K,t(Ks$Knm)))*theta*deltaG/2 - sum(1/omega2)*theta*deltaG/2
            )
    }
    mbfs=apply(bfs,1,max)
    log(apply(exp(bfs-mbfs),1,mean,na.rm=T))+mbfs
}
















