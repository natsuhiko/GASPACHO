# GASPACHO

GAuSsian Processes for Association mapping leveraging Cell HeterOgeneity.

# Fitting GP-LVM
The example log normalised CPM (count per million) data is available [here](https://drive.google.com/file/d/1voKdmSMBHW_UET2TKdezMX4szREhEMgN/view?usp=sharing).

	% R
	> metadata = readRDS("Data/metadata.RDS")
	> cpm = readRDS("Data/log_cpm_4999_22188.RDS")
	> init_param = readRDS("init_param.RDS")
	> 

