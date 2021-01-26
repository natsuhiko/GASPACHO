GPMix = function(Y, W, xi, ti, lam, th, K=2){
    M = length(ti)
    N = length(xi)
    J = ncol(Y)
    K=function(xi, ti, lam, th){th*exp(-outer(xi,ti,"-")^2/lam)}
    KMM = K(ti, ti)
    KNM = K(xi, ti)
    tmp = eigen(KMM,T)
    P = tmp[[2]]
    L = sqrt(tmp[[1]])
    
    for(itr in 1:100){
        W1=W*0
        s2new = 0
        yhat = NULL
        for(j in 1:J){for(k in 1:K){
            yk = c(Y[,j],rep(0,M))
            Xk = rbind(sqrt(W[,k]/sigma2) * th * KNM %*% P %*% diag(1/L)/sqrt(th), diag(M))
            Xk.qr = qr(Xk)
            titsias = W[,k]*th*(1 - rowSums((KNM%*%P%*%diag(1/L))^2))/sigma2
            Qk = qr.Q(Xk.qr)
            Rk = qr.R(Xk.qr)
            yhat = c(yhat, (Qk[1:N,]%*%(t(Qk)%*%yk))/W[,k]*sqrt(th))
            eps = c(yk - Qk%*%(t(Qk)%*%yk))[1:N]
            W1[,k] = W1[,k] - eps^2/2 - rowSums(Qk[1:N,]^2)/2 - titsias/2
            s2new = s2new + sum(eps^2) + sum(Qk[1:N,]^2) + sum(titsias)
        }}
        W = exp(W1/W)/rowSums(exp(W1/W))
        sigma2 = s2new*sigma2/N/J
        
    }
    
}
