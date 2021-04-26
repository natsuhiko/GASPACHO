GPLVM <-
function(Yt, X, Fixed=NULL, ccat=NULL, delta0=NULL, sigma0=NULL, zeta0=NULL, omega0=NULL, Beta0=NULL, Alpha0=NULL, Q0=0, ITRMAX=20, forced=T, Intercept=T, Xi0, Ta0, rho0=NULL, Nu0=NULL, theta0=NULL, DEBUG=T, Kernel="SE", Verbose=2, Plot=F){

    library(Matrix)
	if(nrow(X)!=ncol(Yt)){warning("Input TPM matrix Y is not compatible with covariate matrix X."); return()}
	if(!is.data.frame(X)){warning("Covariate matrix X is not a data frame."); return()}
	
	nh = 1
	Z = rep(1,ncol(Yt))
	isnum = 0
    N=ncol(Yt)
	for(i in seq(ncol(X))){
		if(is.character(X[,i]) || is.factor(X[,i])){
			isnum = c(isnum, 0)
			Z1 = array(0,c(N,length(table(X[,i]))))
			Z1[!is.na(X[,i]),] = model.matrix(~0+X[,i])
		}else{
			isnum = c(isnum, 1)
			Z1 = matrix(scale(as.numeric(X[,i])),N)
			Z1[is.na(Z1)]=0
		}
		Z  = cbind(Z, Z1)
		nh = c(nh, ncol(Z1))
	}
    if(Intercept==F){
        nh=nh[-1]; Z = Z[,-1]; isnum=isnum[-1]
        names(nh)=c(names(X))
    }else{
        names(nh)=c("Intercept", names(X))
    }
	print(dim(Z))
	
	
	isnum = cumsum(1-isnum)*isnum
	for(i in seq(max(isnum))){
		if(sum(isnum==i)>1){
			if(!forced){
				mflag = readline(paste("Do you want to merge ", paste(names(nh)[isnum==i],collapse=","), " factors? [N/y]: ",sep=""))
			}else{
				mflag="N"
			}
			if(sum(mflag%in%c("y","Y"))){
				nh    = Collapse(nh, isnum==i)
				isnum = Collapse(isnum, isnum==i)
			}
		}
	}
	
    K = length(nh)
	H = diag(K)[rep(seq(K),nh),]
    
    
    # fixed effect
    Wfix=array(0,c(ncol(Yt),0))
    if(!is.null(Fixed)){
        H1=H[,names(nh)%in%Fixed,drop=F]
        H = H[ apply(H1,1,sum)==0,!names(nh)%in%Fixed,drop=F]
        Wfix= Z[,apply(H1,1,sum)==1, drop=F][,-1,drop=F]
        Z = Z[,apply(H1,1,sum)==0, drop=F]
        nh = nh[!names(nh)%in%Fixed]
    }
    
    X=Z
    if(is.null(ccat)){ccat=rep(1,length(nh))}
    
    
    # Parameters init
    print("Parameter Init")
    
    H0=getH0(nh) # covariate stratification for Z
    H1=getH1(nh,ccat) # parameter expansion for delta
    H2=getH2(nh,ccat) # parameter stratification for delta
    Bs=Cs=NULL
    for(i in 1:length(nh)){ if(ccat[i]>2){ Bs=c(Bs,list(getB1(nh[i]))) }else{ Bs=c(Bs,list(diag(nh[i]))) } }
    for(i in 1:length(nh)){ if(ccat[i]>2){ Cs=c(Cs,list(getC1(nh[i]))) }else{ Cs=c(Cs,list(diag(nh[i]))) } }
    
    # LINEAR
    if(Kernel=="Linear"){
        Data = list(Kernel=Kernel,W=Wfix,Z=Z,nh=nh,H0=H0,H1=H1,H2=H2,ccat=ccat,Bs=Bs,Cs=Cs,
        MatrixDims=list(J=nrow(Yt),N=ncol(Yt),M=Q0,Q=Q0,K1=ifelse(is.null(Wfix),0,ncol(Wfix)),K2=ncol(Z),L=ncol(H1),P=nrow(H2)))
        Param = list(Xi=array(0,c(N,Q0)),Ta=NA)
    }else{
        Data = list(Kernel=Kernel,W=Wfix,Z=Z,nh=nh,H0=H0,H1=H1,H2=H2,ccat=ccat,Bs=Bs,Cs=Cs,
                MatrixDims=list(J=nrow(Yt),N=ncol(Yt),M=unique(unlist(lapply(Ta0,nrow))),Q=unlist(lapply(Ta0,ncol)),K1=ifelse(is.null(Wfix),0,ncol(Wfix)),K2=ncol(Z),L=ncol(H1),P=nrow(H2)))
        Param = list(Xi=Xi0,Ta=Ta0)
    
        Param = c(Param, list(UpdateKernelParams=rep(7,length(Ta0))))
        
        if(is.null(rho0)){  Param = c(Param, list(rho=rep(1,ncol(Ta)))) }else{ Param = c(Param, list(rho = rho0))   }
        
        PriorXi2=list(mu=array(0,c(Data$MatrixDims$N,Data$MatrixDims$Q[2])), var=array(1,c(Data$MatrixDims$N,Data$MatrixDims$Q[2])))
        Param=c(Param, PriorXi2=list(PriorXi2))
    }
    
    
    if(!is.null(Fixed)){if(is.null(Alpha0)){Param = c(Param, list(Alpha=array(0,c(ncol(Wfix),nrow(Yt)))))}else{Param=c(Param,list(Alpha=Alpha0))}}else{ Param = c(Param, list(Alpha = array(0,c(0,Data$MatrixDims$J))))   }
    if(is.null(zeta0)){ Param = c(Param, list(zeta=rep(0,ncol(Z))))                  }else{ Param = c(Param, list(zeta = zeta0)) }
    
    if(is.null(theta0)){Param = c(Param, list(theta=1))                           }else{ Param = c(Param, list(theta=theta0)) }
    if(is.null(delta0)){Param = c(Param, list(delta = c(1., rep(.1,ncol(H0)-1))))   }else{ Param = c(Param, list(delta = delta0)) }
    
    if(is.null(sigma0)){Param = c(Param, list(sigma2 = rowMeans(Yt^2) - rowMeans(Yt)^2)) }else{ Param = c(Param, list(sigma2 = sigma0)); if(sum(is.na(sigma0))>0){warning("NA in sigma");return()} }
    if(is.null(omega0)){Param = c(Param, list(omega2 = rep(1,ncol(Yt)))) }else{       Param = c(Param, list(omega2 = omega0)); if(sum(is.na(omega0))>0){warning("NA in omega");return()} }
    
    
    
    #load("Param_init_randomTa.Rbin")
    #load("../Param_init.Rbin") # cc based on marker genes
    #load("Param_init.Rbin")
    #load("Param_init_pc.Rbin")
    #load("Param_init_cc.Rbin")
    #load("Param_init_interaction_btw_batch_cc.Rbin")
    #load("../Param_init_cc_based_on_PCs.Rbin")
    #load("Param_init_Ta_from_Xi.Rbin")
    #print("Param_init_Ta_from_Xi.Rbin")
    #load("Param_init_DE_result.Rbin")
    #print("Param_init_DE_result.Rbin")
    #load("Param400.Rbin")

    #load("Param1260.Rbin")
    #load("PriorXi2.Rbin")
    #Param=c(Param, PriorXi2=list(PriorXi2))
    #print(c(Data$MatrixDims$N,Data$MatrixDims$Q[2]))
    
    #print(names(Param))
    #print(dim(Param$PriorXi2$var))
    #print("Param1260 with PriorXi2")
    #Param$delta=c(Param$delta,0.3)
    #Param$zeta = c(Param$zeta,0,0)
    
    
    Param=c(Param, updateKernel(Data, Param))
    Param=c(Param, initDelta(Data, Param))
    
    # random effect setting
    print(data.frame(names(nh),nh,ccat))
    
    # matrix prep
    print("Matrix Prep");
    #load("Param4750.Rbin");cat("Param4750.Rbin")
    YMat = updateYMat(Yt, Data, Param)
    
    
    # Debug
    if(DEBUG){
        cat("dim(Yt)=");print(dim(Yt))
        cat("dim(Z)=");print(dim(Z))
        cat("dim(X)=");print(dim(X))
        cat("dim(H0)=");print(dim(H0))
        cat("dim(H1)=");print(dim(H1))
        cat("dim(H2)=");print(dim(H2))
        cat("dim(K)=");print(dim(Param$K))
        cat("length(zeta=");print(length(Param$zeta))
        cat("length(delta)=");print(length(Param$delta))
        cat("dim(Phiinv)=");print(dim(Param$Phiinv))
        cat("dim(Xi)=");print(dim(Param$Xi))
        print(names(Param))
        print(unlist(Data$MatrixDims))
    }
    
    ##################
    #    ITERATION   #
    ##################
    print("Iteration starts")
    tlb0 = updateTLB(YMat,Data,Param,F)
    tlb.all=NULL
    for(itr in 1:ITRMAX){
        cat(c(paste("[",itr,"]",sep=""),""),sep="\n");
        
        # intermediate parameter saving off
        if(Plot){if(itr%%20==0){
            save(Param,file=paste("Param",itr,".Rbin",sep=""))
            #png("pairs.png",width=2000,height=2000,res=150);pairs(rbind(Param$Xi,Param$Ta),col=rep(1:2,c(22188,50)));dev.off()
            #hoge=Param$Xi
            #save(hoge,file="hoge.Rbin")
        }}
        
        # Rho & Xi
        if(sum(Data$MatrixDims$Q)>0){
            Param[c("rho","Xi","Ta","ValRhoXi","GradRhoXi")] = updateRhoXiTa(Yt, YMat, Data, Param, F, Verbose=Verbose)
            Param[c("K","Knm")] = updateKernel(Data,Param)
            YMat[["YtOinvKnm"]] = as.matrix(Yt%*%(Param$Knm/Param$omega2))
            tDinv = dbind(Param$Dinv,Param$K/rep(Param$theta,Data$MatrixDims$M))
            tZ = cbind(Data$Z,Param$Knm)
            Param[["Phiinv"]] = tDinv + t(tZ/Param$omega2)%*%tZ
            hist(Param$Xi[[1]]%%pi,breaks=100);abline(v=Param$Ta[[1]],col=2:4)
            #plot(Param$Ta[[3]],xlim=c(0,17),ylim=c(0,25))
        }
        
        # Delta
        Param[c("zeta", "delta", "LD", "Delta", "Linv", "Dinv", "Phiinv", "theta", "ssdelta", "graddelta", "Hinvdelta")] = updateDelta(YMat, Data, Param, F, Verbose=Verbose)
        
        # Alpha
        if(Data$MatrixDims$K1>0){
            Param["Alpha"] = updateAlpha(YMat, Data, Param, F, Verbose=Verbose)
            YMat[["ASinvYt"]] = as.matrix(t(t(Param$Alpha)/Param$sigma2)%*%Yt)
        }
        
        # Omega
        #if(itr%%20==0){
            Param["omega2"]=updateOmega(Yt,YMat, Data, Param,F, Verbose=Verbose)
            YMat[c("YtOinvW","YtOinvZ","YtOinvKnm","Yt2Oinv1")] = updateYMatOmega(Yt, Data, Param)
            tDinv = dbind(Param$Dinv,Param$K/rep(Param$theta,Data$MatrixDims$M))
            tZ = cbind(Data$Z,Param$Knm)
            Param[["Phiinv"]] = tDinv + t(tZ/Param$omega2)%*%tZ
            #}
        
        # Sigma
        Param["sigma2"]=updateSigma(YMat, Data, Param,F,Yt, Verbose=Verbose)
        YMat[c("ASinvYt","BSinvYt","USinvYt","OneSinvYt")] = updateYMatSigma(Yt, Data, Param)
        #plot(sigma0, Param$sigma2);abline(0,1)
        #plot(omega0,Param$omega2);abline(0,1)
        #hist(Param$omega2)
        
        
        # Titsias lower bound
        tlb = updateTLB(YMat,Data,Param,F)
        cat("TLB=");print(tlb)
        if(is.na(tlb)){break}
        tlb.all=c(tlb.all,tlb)
        if(Plot){png("tlb.png",width=3000,height=1000,res=150);par(mfcol=c(1,3));plot(tlb.all);plot(rev(rev(tlb.all)[1:100]));plot(rev(rev(tlb.all)[1:20]));dev.off()}
    }
    list(Data=Data,Param=Param,tlb=tlb.all)
}



updateYMat=function(Yt, Data, Param){
    flag=rep(1:4,c(Data$MatrixDims$K1,Data$MatrixDims$K2,sum(Data$MatrixDims$M),1))
    flag2=rep(1:2,c(Data$MatrixDims$K1,1))
    tmp1 = as.matrix(Yt%*%(cbind(Data$W,Data$Z,Param$Knm,1)/Param$omega2))
    tmp2 = as.matrix(t(cbind(t(Param$Alpha),rep(1,Data$MatrixDims$J))/Param$sigma2)%*%Yt)
    print(dim(tmp2))
    print(length(flag))
    list(
        YtOinvW = tmp1[,flag==1,drop=F],
        YtOinvZ = tmp1[,flag==2,drop=F],
        YtOinvKnm = tmp1[,flag==3,drop=F],
        Yt2Oinv1 = as.numeric((Yt^2)%*%(1/Param[["omega2"]])),
        ASinvYt = tmp2[flag2==1,,drop=F],
        OneSinvYt = as.numeric(tmp2[flag2==2,])
    )
}

updateYMatSigma=function(Yt, Data, Param){
    flag=rep(1:2,c(Data$MatrixDims$K1,1))
    tmp2 = as.matrix(t(cbind(t(Param$Alpha),rep(1,Data$MatrixDims$J))/Param$sigma2)%*%Yt)
    list(
        ASinvYt = tmp2[flag==1,,drop=F],
        OneSinvYt = as.numeric(tmp2[flag==2,])
    )[c("ASinvYt","BSinvYt","USinvYt","OneSinvYt")]
}

updateYMatOmega=function(Yt, Data, Param){
    flag=rep(1:4,c(Data$MatrixDims$K1,Data$MatrixDims$K2,sum(Data$MatrixDims$M),1))
    batch=list(runif(nrow(Yt))<0.1, runif(ncol(Yt))<0.1)
    tmp1 = as.matrix(Yt%*%(cbind(Data$W,Data$Z,Param$Knm,1)/Param$omega2))
    
    list(
    YtOinvW = tmp1[,flag==1,drop=F],
    YtOinvZ = tmp1[,flag==2,drop=F],
    YtOinvKnm = tmp1[,flag==3,drop=F],
    Yt2Oinv1 = as.numeric((Yt^2)%*%(1/Param[["omega2"]]))
    )[c("YtOinvW","YtOinvZ","YtOinvKnm","Yt2Oinv1")]
    
}
