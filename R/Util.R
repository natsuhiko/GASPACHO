# H = rbind(c(H11),c(H12),c(H13),...,c(H22),c(H33),c(H44),...)
# ( H11 H12 H13 H14 ... )  %*% c(t(X))
# ( H21 H22  0   0  ... )
# ( H31  0  H33  0  ... )
# ( H41  0   0  H44 ... )
# ( ... ... ... ... ... )
BlockProduct = function(H, X, debug=F){
    Q = ncol(X)
    N = nrow(X)-1
    # H : (2N+1) x Q^2
    # X : (N+1) x Q
    Xone = X%x%t(rep(1,Q))
    oneX = t(rep(1,Q))%x%X
    B = diag(Q)[rep(1:Q,rep(Q,Q)),]
    res = c(
        rowSums(matrix(colSums((H[1:(N+1),]*Xone)),Q)),
        c(t((t(t(H[2:(N+1),])*oneX[1,])+H[(2+N):(2*N+1),]*oneX[-1,])%*%B))
    )
    if(debug){
        HH = diag((N+1)*Q)
        for(i in 1:(N+1)){
            for(j in 1:(N+1)){
                if(i==1 && j==1){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = matrix(H[1,],Q)
                }else if(i==1 && j>1){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = matrix(H[j,],Q)
                }else if(j==1 && i>1){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = t(matrix(H[i,],Q))
                }else if(i==j){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = t(matrix(H[N+i,],Q))
                }
            }
        }
        image(HH)
        plot(res, HH%*%c(t(X)))
    }
    t(matrix(res,Q))
}

BlockProduct2=function(S,Y){
    Q = ncol(S)
    rbind(t(matrix(c(S[1,]%*%t(c(t(Y)))),Q^2)), (t(rep(1,Q))%x%S[-1,])*(Y[-1,]%x%t(rep(1,Q))))
}

BFGS = function(H,S,Y,G,debug=F){
    print(dim(H))
    print(dim(S))
    print(dim(Y))
    print(dim(G))
    sy = sum(S*Y)
    Hy = BlockProduct(H,Y)
    yHy = sum(Hy*Y)
    Hnew = H + (sy+yHy)/sy^2*BlockProduct2(S,S) - BlockProduct2(Hy,S)/sy - BlockProduct2(S,Hy)/sy
    
    if(debug){
        Q = ncol(S)
        N = nrow(S)-1
        sy=sum(S*Y)
        y=c(t(Y))
        s=c(t(S))
        HH = diag((N+1)*Q)
        for(i in 1:(N+1)){
            for(j in 1:(N+1)){
                if(i==1 && j==1){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = matrix(H[1,],Q)
                }else if(i==1 && j>1){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = matrix(H[j,],Q)
                }else if(j==1 && i>1){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = t(matrix(H[i,],Q))
                }else if(i==j){
                    HH[((i-1)*Q+1):(i*Q),((j-1)*Q+1):(j*Q)] = t(matrix(H[N+i,],Q))
                }
            }
        }
        HHnew=(diag(length(s))-s%*%t(y)/sy)%*%HH%*%(diag(length(s))-y%*%t(s)/sy)+s%*%t(s)/sy
        HHnew[HH==0]=0
    
        plot(c(t(BlockProduct(Hnew,G))), c(HHnew%*%c(t(G))))
    }
    
    list(H=Hnew, ss=BlockProduct(Hnew,G))
}

LBFGS_old=function(X, G){
    print(dim(X))
    print(dim(G))
    M = ncol(X)
    S = X[,1:(M-1),drop=F]-X[,2:M,drop=F]
    Y = G[,1:(M-1),drop=F]-G[,2:M,drop=F]
    q = G[,1]
    a = rep(0, M-1)
    for(i in 1:(M-1)){
        a[i] = sum(S[,i]*q)/sum(S[,i]*Y[,i])
        if(!is.na(a[i])){q = q - a[i]*Y[,i]}else{print("ss too small")}
    }
    gammak = sum(S[,1]*Y[,1])/sum(Y[,1]^2)
    z = gammak*q
    for(i in (M-1):1){
        bi = sum(Y[,i]*z)/sum(S[,i]*Y[,i])
        if(!is.na(a[i]+bi)){z = z + S[,i]*(a[i]-bi)}
    }
    -z
}

LBFGS=function(X, G){
    print(dim(X))
    print(dim(G))
    M = ncol(X)
    S = X[,1:(M-1),drop=F]-X[,2:M,drop=F]
    Y = G[,1:(M-1),drop=F]-G[,2:M,drop=F]
    q = G[,1]
    a = rep(0, M-1)
    flag = rep(T,M-1)
    for(i in 1:(M-1)){
        a[i] = sum(S[,i]*q)/sum(S[,i]*Y[,i])
        if(!is.na(a[i])){q = q - a[i]*Y[,i]}else{flag[i]=F; print(i); print("ss too small")}
    }
    effi = rev(seq(M-1)[flag])
    gammak = sum(S[,min(effi)]*Y[,min(effi)])/sum(Y[,min(effi)]^2)
    print(gammak)
    z = gammak*q
    for(i in effi){
        bi = sum(Y[,i]*z)/sum(S[,i]*Y[,i])
        if(!is.na(a[i]+bi)){z = z + S[,i]*(a[i]-bi)}
    }
    -z
}

Rep = function(X,n){
    res = NULL
    for(i in 1:n){
        res = rbind(res,X)
    }
    res
}

cprod=function(a,b){
    res = NULL
    for(i in 1:ncol(a)){
        res = cbind(res, a[,i]*b)
    }
    res
}
Cbind=function(lis){do.call(cbind,lis)}
dbind=function(...){
    dims=NULL
    if(length(list(...))==1){
        lis = list(...)[[1]]
    }else{
        lis = list(...)
    }
    N = length(lis)
    for(i in 1:length(lis)){
        x = dim(lis[[i]])
        if(is.null(x)){
            dims = rbind(dims, c(length(lis[[i]]),1))
        }else{
            dims = rbind(dims, x)
        }
    }
    A = array(0, apply(dims,2,sum))
    cdims = apply(rbind(c(0,0),dims),2,cumsum)
    for(i in 1:N){
        A[(cdims[i,1]+1):cdims[i+1,1], (cdims[i,2]+1):cdims[i+1,2]]= lis[[i]]
    }
    A
}

dbind0=function(a,b){
    d = array(0,dim(a)+dim(b))
    d[1:nrow(a),1:ncol(a)]=a
    d[(nrow(a)+1):(nrow(a)+nrow(b)), (ncol(a)+1):(ncol(a)+ncol(b))] = b
    d
}

Solve = function(x,y=diag(nrow(x))){
	r  = chol((x+t(x))/2)
    if(is.character(r)){
        return(y*NA)
    }
	ra = forwardsolve(t(r),y)
	backsolve(r, ra)
}
logDet=function(x){
    x=(x+t(x))/2;
    eval = try(eigen(x,T)[[1]]);
    if(is.character(eval)){
        return(NA)
    }
    sum(log(eval))
}

logDetChol=function(x){
    x=(x+t(x))/2;
    2*sum(log(diag(chol(x))))
}


getB <-
function(nh,ccat){
	id=NULL; 
	J=length(nh)
	for(i in 1:J){
		if(ccat[i]<3){
			id=c(id,nh[i])
		}else{
			id=c(id,sum(nh[i]:1))
		}
	};
	B=array(0,c(sum(id),sum(nh)))
	nh=cumsum(c(1,nh))
	id=cumsum(c(1,id))
	for(i in 1:J){
		if(ccat[i]<3){
			B[id[i]:(id[i+1]-1), nh[i]:(nh[i+1]-1)] = diag(nh[i+1]-nh[i])
		}else{
			B[id[i]:(id[i+1]-1), nh[i]:(nh[i+1]-1)] = getB1(nh[i+1]-nh[i])
		}
	}
	B
}

getB1 <-
function(n){
	B=diag(n); 
	if(n>1){
		for(k in (n-1):1){
			B=rbind(B,cbind(matrix(rep(0,(n-k)*k),k),diag(k)))
		}
	}; 
	B
}

getC <-
function(nh,ccat){
	id=NULL; 
	J=length(nh)
	for(i in 1:J){
		if(ccat[i]<3){
			id=c(id,nh[i])
		}else{
			id=c(id,sum(nh[i]:1))
		}
	};
	B=array(0,c(sum(id),sum(nh)))
	nh=cumsum(c(1,nh))
	id=cumsum(c(1,id))
	for(i in 1:J){
		if(ccat[i]<3){
			B[id[i]:(id[i+1]-1), nh[i]:(nh[i+1]-1)] = diag(nh[i+1]-nh[i])
		}else{
			B[id[i]:(id[i+1]-1), nh[i]:(nh[i+1]-1)] = getC1(nh[i+1]-nh[i])
		}
	}
	B
}

getC1 <-
function(n){
	diag(n)[rep(1:n,n:1),,drop=F]
}

getH0 <-
function(nh){
	diag(length(nh))[rep(seq(length(nh)),nh),]
}

getH1 <-
function(nh,ccat){
	id=NULL; 
	for(i in 1:length(nh)){
		if(ccat[i]==1){
			id=c(id,nh[i])
		}else if(ccat[i]==2){
			id=c(id,rep(1,nh[i]))
		}else{
			id=c(id,rep(1,sum(nh[i]:1)))
		}
	}; 
	diag(length(id))[rep(1:length(id),id),]
}

getH2 <-
function(nh,ccat){
	nh[ccat==1]=1; 
	nh[ccat==2]=nh[ccat==2] 
	nh[ccat==3]=sapply(nh[ccat==3],function(x1){sum(x1:1)}); 
	diag(length(nh))[rep(seq(length(nh)),nh),]
}

Collapse <-
function(x,flag){
	x1=x[cumsum(flag)==0]
	x2=x[rev(cumsum(rev(flag)))==0]
	res = sum(x[flag])
	names(res) = paste(names(x[flag]),collapse="/")
	return(c(x1, res, x2))
}


emn <-
function(X){
    N=ncol(X)
    geta = apply(X,1,max)-10
    X = X-geta
    print(range(-geta))
    #x[x==Inf&!is.na(x)]=max(x[x<Inf&!is.na(x)],na.rm=T)
    p=c(0.9,rep(0.1/N,N))
    lkhd=0
    lkhd0=-10000000000
    for(itr in 1:1000){
        d=c(exp(-geta)*p[1]+exp(X)%*%p[-1])
        z=cbind(exp(-geta)*p[1],t(t(exp(X))*p[-1]))/d
        p=c(apply(z,2,mean,na.rm=T))
        p=p/sum(p)
        lkhd=mean(log(d),na.rm=T)
        if(lkhd>lkhd0 && abs(lkhd-lkhd0)<1e-7){cat("Converged.\n");return(z[,c(2:(N+1),1)])}else{print(c(lkhd0,lkhd,p)); lkhd0=lkhd}
    }
}


mvnMix = function(Mu, Mu0=NULL, Sinv){
    P = nrow(Mu)
    if(is.null(Mu0)){
        Mu0 = matrix(rnorm(P*10,0,sd(c(Mu))),P)/1000+rowMeans(Mu)
    }
    p=rep(1,ncol(Mu0))
    lkhd0=lkhd=-1e10
    for(itr in 1:1000){
        f = colSums(Mu*(Sinv%*%Mu))
        g = colSums(Mu0*(Sinv%*%Mu0))
        H = t(Mu)%*%Sinv%*%Mu0
        Z = - (outer(f,g,"+")-2*H)/2
        Z = exp(Z-max(Z)+10)%*%diag(p)
        
        d = rowSums(Z)
        lkhd = sum(log(d))
        if(abs(lkhd-lkhd0)<1e-7){break}else{lkhd0=lkhd;print(c(lkhd,lkhd0))}
        Z = Z/d
        p = colSums(Z)
        p=p/sum(p)
        Mu0 = Mu%*%Z%*%diag(1/colSums(Z))
    }
    list(Mu0,p)
}

mvnMix = function(Y, X, omega2, Z, Phiinv, K, W){
    XSinvX = t(X/omega2)%*%X - (t(X/omega2)%*%Z)%*%Solve(Phiinv, t(Z)%*%(X/omega2))
    YSinvX = t(Y/omega2)%*%X - (t(Y/omega2)%*%Z)%*%Solve(Phiinv, t(Z)%*%(X/omega2))
    for(itr in 1:20){
        # M step
        p=colMeans(W); p=p/sum(p)
        Yk = Y%*%W
        B = NULL
        for(k in 1:ncol(W)){
            B = cbind(B, Solve(sum(W[,k])*XSinvX+K, t(X/omega2)%*%Yk[,k]-(t(X/omega2)%*%Z)%*%Solve(Phiinv, t(Z)%*%(Yk[,k]/omega2))))
        }
        
        # E step
        W = - outer(colSums(Y^2/omega2)+colSums((t(Z/omega2)%*%Y)*Solve(Phiinv, t(Z)%*%(Y/omega2))), colSums(B*(XSinvX%*%B)), "+") + 2*YSinvX%*%B
        W = exp(W-apply(W,1,max))%*%diag(p)
        W = W/rowSums(W)
        
        print(p)
    }
    list(W=W, XB=X%*%B, p=p)
}

# Z    : cbind(Z, Knm)
# Dinv : dbind(Dinv, theta*Kinv)
# X    : Knm(imm) Kinv(imm)
# K    : theta*K(imm)
mvnMix2 = function(Y, X, B, omega2, Z, Dinv, K, W, itrmax=100){
    print(dim(W))
    KK=ncol(W)
    #W = cbind(W/2,0.5)
    y2 = colSums(Y^2/omega2)
    XOinvY = t(X/omega2)%*%Y
    ZOinvY = t(Z/omega2)%*%Y
    ZOinvX = t(Z/omega2)%*%X
    XOinvX = t(X/omega2)%*%X
    Grad = NULL
    Vals = NULL
    lkhd_all=NULL
    for(itr in 1:itrmax){
        print(itr)
        # M step
        p=colMeans(W)+1e-7; p=p/sum(p)
        print(p)
        grad = NULL
        for(k in 1:KK){
            Zk = cbind(Z, X%*%B[,k])
            ZkOinvY = rbind(ZOinvY, t((X%*%B[,k])/omega2)%*%Y)
            ZkOinvX = t(Zk)%*%(X/omega2)
            Phiinvk = dbind(Dinv,1) + t(Zk/omega2)%*%Zk
            XVinvX = XOinvX - t(ZkOinvX)%*%Solve(Phiinvk, ZkOinvX)
            XVinvY = XOinvY - t(ZkOinvX)%*%Solve(Phiinvk, ZkOinvY)
            
            grad = cbind(grad,-sum(W[,k])*XVinvX%*%B[,k] + (XVinvY*W[,k])%*%t(XVinvY)%*%B[,k] - Solve(K,B[,k]))
        }
        # L-BFGS
        if(is.null(Grad)){
            Vals = matrix(c(B),   length(c(B)))
            Grad = matrix(c(grad),length(c(grad)))
            ss = -grad/10000000
            print(range(c(ss)))
        }else{
            Vals = cbind(c(B),Vals)
            Grad = cbind(c(grad),Grad)
            Vals=Vals[,1:min(20,ncol(Vals))]
            Grad=Grad[,1:min(20,ncol(Grad))]
            ss = matrix(-LBFGS(Vals, Grad), nrow(B))
        }
        # lkhd
        lkhd = 0
        for(k in 1:KK){
            Zk = cbind(Z, X%*%B[,k])
            ZkOinvY = rbind(ZOinvY, t((X%*%B[,k])/omega2)%*%Y)
            Phiinvk = dbind(Dinv, 1) + t(Zk/omega2)%*%Zk
            Rk = chol(Solve(dbind(Dinv, 1)))
            
            lkhd = lkhd + sum((-logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2)*W[,k])
            lkhd = lkhd - sum(Solve(K,B[,k])*B[,k])/2 + log(p[k])*sum(W[,k])
        }
        # line search
        for(r in c(0:10,20,50)){
            B1 = B-ss/2^r
            lkhd1 = 0
            for(k in 1:KK){
                Zk = cbind(Z, X%*%B1[,k])
                ZkOinvY = rbind(ZOinvY, t((X%*%B1[,k])/omega2)%*%Y)
                Phiinvk = dbind(Dinv, 1) + t(Zk/omega2)%*%Zk
                Rk = chol(Solve(dbind(Dinv, 1)))
                
                lkhd1 = lkhd1 + sum((-logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2)*W[,k])
                lkhd1 = lkhd1 - sum(Solve(K,B1[,k])*B1[,k])/2 + log(p[k])*sum(W[,k])
            }
            
            if(lkhd1*ifelse(sign(lkhd1)>0,1.0000001,0.9999999)>lkhd || r==50){
                print(c(r,lkhd,lkhd1))
                B = B1
                lkhd_all=c(lkhd_all,lkhd1)
                break;
            }else{
                print(c(r,lkhd,lkhd1))
            }
        }
        if(itr>1)plot(lkhd_all[-1])
        
        # E step
        W = NULL
        penaltyAll = 0
        for(k in 1:KK){
            Zk = cbind(Z, X%*%B[,k])
            ZkOinvY = rbind(ZOinvY, t((X%*%B[,k])/omega2)%*%Y)
            ZkOinvX = t(Zk/omega2)%*%X
            Phiinvk = dbind(Dinv, 1) + t(Zk/omega2)%*%Zk
            Rk = chol(Solve(dbind(Dinv, 1)))
            
            # hessian for laplace approximation
            XVkinvX = t(X/omega2)%*%X - t(ZkOinvX)%*%Solve(Phiinvk, ZkOinvX)
            penaltyk = -logDet(sum(t(XVkinvX*B[,k])*B[,k])*XVkinvX + (XVkinvX%*%B[,k])%*%t(XVkinvX%*%B[,k]))/2 + log(2*pi)/2*nrow(B)
            #penaltyAll = penaltyAll + penaltyk
            W = cbind(W, -logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2 )
        }
        #Zk = Z
        #ZkOinvY = ZOinvY
        #Phiinvk = Dinv + t(Zk/omega2)%*%Zk
        #Rk = chol(Solve(Dinv))
        #W = cbind(W+penaltyAll/KK, -logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2)
        
        W = exp(W-apply(W,1,max))%*%diag(p)
        W = W/rowSums(W)
        
        
        if(is.na(p[1])){break}
    }
    
    Phiinv = Dinv + t(Z/omega2)%*%Z
    BXVinvXB = t(B)%*%t(X/omega2)%*%(X%*%B) - t(B)%*%t(ZOinvX)%*%Solve(Phiinv, ZOinvX%*%B)
    Del = t((t(B)%*%XOinvY - t(B)%*%t(X/omega2)%*%Z%*%Solve(Phiinv, ZOinvY)) / diag(BXVinvXB))
    
    list(W=W, X=X, B=B, p=p, Delta=Del)
}






mvnMix2 = function(Y, X, B, omega2, Z, Dinv, K, W, itrmax=100){
    print(dim(W))
    KK=ncol(W)
    W = cbind(W*0.9,0.1)
    y2 = colSums(Y^2/omega2)
    XOinvY = t(X/omega2)%*%Y
    ZOinvY = t(Z/omega2)%*%Y
    ZOinvX = t(Z/omega2)%*%X
    XOinvX = t(X/omega2)%*%X
    Phiinv = Dinv + t(Z/omega2)%*%Z
    R = chol(Solve(Dinv))
    W0 = -logDet(diag(nrow(Phiinv))+R%*%t(Z/omega2)%*%Z%*%t(R))/2 - y2/2 + colSums(ZOinvY*Solve(Phiinv, ZOinvY))/2
    
    Grad = NULL
    Vals = NULL
    lkhd_all=NULL
    for(itr in 1:itrmax){
        print(itr)
        # M step
        p=colMeans(W)+1e-7; p=p/sum(p)
        print(p)
        grad = NULL
        for(k in 1:KK){
            Zk = cbind(Z, X%*%B[,k])
            ZkOinvY = rbind(ZOinvY, t((X%*%B[,k])/omega2)%*%Y)
            ZkOinvX = t(Zk)%*%(X/omega2)
            Phiinvk = dbind(Dinv,1) + t(Zk/omega2)%*%Zk
            XVinvX = XOinvX - t(ZkOinvX)%*%Solve(Phiinvk, ZkOinvX)
            XVinvY = XOinvY - t(ZkOinvX)%*%Solve(Phiinvk, ZkOinvY)
            
            grad = cbind(grad,-sum(W[,k])*XVinvX%*%B[,k] + (XVinvY*W[,k])%*%t(XVinvY)%*%B[,k] - Solve(K,B[,k]))
        }
        # L-BFGS
        if(is.null(Grad)){
            Vals = matrix(c(B),   length(c(B)))
            Grad = matrix(c(grad),length(c(grad)))
            ss = -grad/10000000
            print(range(c(ss)))
        }else{
            Vals = cbind(c(B),Vals)
            Grad = cbind(c(grad),Grad)
            Vals=Vals[,1:min(20,ncol(Vals))]
            Grad=Grad[,1:min(20,ncol(Grad))]
            ss = matrix(-LBFGS(Vals, Grad), nrow(B))
        }
        # lkhd
        lkhd = sum((-logDet(diag(nrow(Phiinv))+R%*%t(Z/omega2)%*%Z%*%t(R))/2 - y2/2 + colSums(ZOinvY*Solve(Phiinv, ZOinvY))/2)*W[,KK+1])
        for(k in 1:KK){
            Zk = cbind(Z, X%*%B[,k])
            ZkOinvY = rbind(ZOinvY, t((X%*%B[,k])/omega2)%*%Y)
            Phiinvk = dbind(Dinv, 1) + t(Zk/omega2)%*%Zk
            Rk = chol(Solve(dbind(Dinv, 1)))
            
            lkhd = lkhd + sum((-logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2)*W[,k])
            lkhd = lkhd - sum(Solve(K,B[,k])*B[,k])/2 + log(p[k])*sum(W[,k])
        }
        # line search
        for(r in c(0:10,20,50)){
            B1 = B-ss/2^r
            lkhd1 = sum((-logDet(diag(nrow(Phiinv))+R%*%t(Z/omega2)%*%Z%*%t(R))/2 - y2/2 + colSums(ZOinvY*Solve(Phiinv, ZOinvY))/2)*W[,KK+1])
            for(k in 1:KK){
                Zk = cbind(Z, X%*%B1[,k])
                ZkOinvY = rbind(ZOinvY, t((X%*%B1[,k])/omega2)%*%Y)
                Phiinvk = dbind(Dinv, 1) + t(Zk/omega2)%*%Zk
                Rk = chol(Solve(dbind(Dinv, 1)))
                
                lkhd1 = lkhd1 + sum((-logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2)*W[,k])
                lkhd1 = lkhd1 - sum(Solve(K,B1[,k])*B1[,k])/2 + log(p[k])*sum(W[,k])
            }
            
            if(lkhd1*ifelse(sign(lkhd1)>0,1.0000001,0.9999999)>lkhd || r==50){
                print(c(r,lkhd,lkhd1))
                B = B1
                lkhd_all=c(lkhd_all,lkhd1)
                break;
            }else{
                print(c(r,lkhd,lkhd1))
            }
        }
        if(itr>1)plot(lkhd_all[-1])
        
        # E step
        W = NULL
        penaltyAll = 0
        for(k in 1:KK){
            Zk = cbind(Z, X%*%B[,k])
            ZkOinvY = rbind(ZOinvY, t((X%*%B[,k])/omega2)%*%Y)
            ZkOinvX = t(Zk/omega2)%*%X
            Phiinvk = dbind(Dinv, 1) + t(Zk/omega2)%*%Zk
            Rk = chol(Solve(dbind(Dinv, 1)))
            
            # hessian for laplace approximation
            XVkinvX = t(X/omega2)%*%X - t(ZkOinvX)%*%Solve(Phiinvk, ZkOinvX)
            penaltyk = -logDet(sum(t(XVkinvX*B[,k])*B[,k])*XVkinvX + (XVkinvX%*%B[,k])%*%t(XVkinvX%*%B[,k]))/2 + log(2*pi)/2*nrow(B)
            #penaltyAll = penaltyAll + penaltyk
            W = cbind(W, -logDet(diag(nrow(Phiinvk))+Rk%*%t(Zk/omega2)%*%Zk%*%t(Rk))/2 - y2/2 + colSums(ZkOinvY*Solve(Phiinvk, ZkOinvY))/2 + penaltyk)
        }
        W = cbind(W, W0)
        
        W = exp(W-apply(W,1,max))%*%diag(p)
        W = W/rowSums(W)
        
        if(is.na(p[1])){break}
    }
    
    
    BXVinvXB = t(B)%*%t(X/omega2)%*%(X%*%B) - t(B)%*%t(ZOinvX)%*%Solve(Phiinv, ZOinvX%*%B)
    Del = t((t(B)%*%XOinvY - t(B)%*%t(X/omega2)%*%Z%*%Solve(Phiinv, ZOinvY)) / diag(BXVinvXB))
    
    list(W=W, X=X, B=B, p=p, Delta=Del)
}











