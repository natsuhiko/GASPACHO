
getResid=function(Yt,gplvm,Kernel="cc"){
    tZ = cbind(gplvm$Data$Z, gplvm$Param$Knm)
    tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
    tA = rbind(gplvm$Param$Alpha,1)
    Phiinv = gplvm$Param$Phiinv
    D = Yt%*%(tZ/gplvm$Param$omega2)
    E = t(tZ/gplvm$Param$omega2)%*%tW
    C0= t(D)-E%*%tA
    C = solve(Phiinv, C0)
    pred=tZ%*%C
    resi=t(Yt)-pred-tW%*%tA
    if(Kernel=="cc"){
        newK = updateKernelPeriodic(gplvm$Param$Xi[[1]], gplvm$Param$Ta[[1]], gplvm$Param$rho[[1]])
        pred = newK$Knm%*%C[-seq(gplvm$Data$MatrixDims$K2),]
    }else{
        #newK = updateKernelSE(gplvm$Param$Xi[[2]][,-c(1,3,5)],reduceDim1(gplvm$Param$Ta[[2]][,-c(1,3,5)]),gplvm$Param$rho[[2]][-c(1,3,5)])
        newK = updateKernelSE(gplvm$Param$Xi[[2]][,-c(1,3,5)], gplvm$Param$Ta[[2]][,-c(1,3,5)], gplvm$Param$rho[[2]][-c(1,3,5)])
        #newK = updateKernelSE(gplvm$Param$Xi[[2]][,-1], gplvm$Param$Ta[[2]][,-1], gplvm$Param$rho[[2]][-1])
        print(dim(newK$Knm))
        pred = newK$Knm%*%C[-seq(gplvm$Data$MatrixDims$K2),]
    }
    list(residual=resi+pred,prediction=pred)
}


getResidualOld = function(Yt, gplvm){
    
    YMat = updateYMat(Yt, gplvm$Data, gplvm$Param)
    
    M = gplvm$Data$MatrixDims$M
    K2= gplvm$Data$MatrixDims$K2
    flag = rep(0:1,c(K2,M))
    
    tZ = cbind(gplvm$Data$Z,gplvm$Param$Knm)
    tW = cbind(gplvm$Data$W, gplvm$Data$Z%*%gplvm$Param$zeta)
    tA = rbind(gplvm$Param$Alpha,1)
    
    D = cbind(YMat$YtOinvZ,YMat$YtOinvKnm)
    E = t(tZ/gplvm$Param$omega2)%*%tW
    B = Solve(gplvm$Param$Phiinv,t(D)-E%*%tA)
    list(predicted=tZ[,flag==1] %*% B[flag==1,], residuals=as.matrix(t(Yt)) - tW%*%tA - tZ[,flag==0] %*% B[flag==0,])
}


