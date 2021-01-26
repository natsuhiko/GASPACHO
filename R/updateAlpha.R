updateAlpha = function(
YMat,
Data,
Param,
DEBUG=F,
Verbose=0){
    
    if(Verbose>0){print("update Alpha")}
    W = Data$W
    Z = Data$Z
    Knm = Param[["Knm"]]
    omega2 = Param[["omega2"]]
    
    tZ = cbind(Z,Knm)
    
    zeta = Param$zeta
    
    
    WtVinv = t(W/omega2) - (t(W/omega2)%*%tZ)%*%Solve(Param$Phiinv)%*%(t(tZ/omega2))
    WtVinvW = WtVinv%*%W
    Atmp = solve(WtVinvW, WtVinv)
    A = solve(WtVinvW, t(YMat$YtOinvW)-(t(W/omega2)%*%tZ)%*%Solve(Param$Phiinv)%*%t(cbind(YMat$YtOinvZ,YMat$YtOinvKnm))) - c(solve(WtVinvW, WtVinv%*%(Z%*%zeta)))
    
    
    
    list(Alpha = A)
}

