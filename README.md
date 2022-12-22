# GASPACHO

GAuSsian Processes for Association mapping leveraging Cell HeterOgeneity.

## Fitting GP-LVM
The example log normalised CPM (count per million) data is available [here](https://drive.google.com/file/d/1voKdmSMBHW_UET2TKdezMX4szREhEMgN/view?usp=sharing). Please dowonload it and put in the `Data` directory.

	% R
	metadata = readRDS("Data/metadata.RDS")
	cpm = readRDS("Data/log_cpm_4999_22188.RDS")
	init_param = readRDS("Data/init_param.RDS")
	gplvm = GPLVM(cpm, metadata,
		Xi     = init_param$Xi, # latent variables
		Ta     = init_param$Ta, # inducing variables
		delta0 = init_param$delta, # variance parameters for metadata
		Alpha0 = init_param$Alpha, # fixed effect (NULL)
		sigma0 = init_param$sigma2, # residual variance for genes
		omega0 = init_param$omega2, # residual variance for cells
		zeta0  = init_param$zeta, # grand mean
		rho0   = init_param$rho, # length parameters for Kernels
		theta0 = init_param$theta, # variance parameter for GP
		ITRMAX = 1000)

## Mapping dynamic eQTLs
In order to map eQTL, GASPACHO requires to set the variance parameter of donor x context interaction effect. It is available from the function `updateDxCSE` which returns the gplvm object with the estimate located at `gplvm$Param$d2dxc`. The function `getBF` then computes both the Bayes factors (BFs) and the Score based statistics for a given genotype dosage matrix (the third argument).  Note that, the column number of the genotype matrix has to be compatible with the donor ID provided by `did`. The example below computes the OAS1 eQTL BFs and Score statistics. The 3,329th row of the `cpm` matrix is the normalised OAS1 expression and `metadata$donor` indicates which cell belngs to the donor or origin.

	% R
	# Estimating the Donor by Context interaction effect
	gplvm = updateDxCSE(as.matrix(cpm), gplvm)
	
	# Computing eQTL Bayes factors
	G = readRDS("Data/G_OAS1.RDS")
	bfs = getBF(yj = as.numeric(cpm[3329,]), gplvm, G = as.matrix(G[,10:77]), 
	            did = as.numeric(as.factor(metadata$donor)))
