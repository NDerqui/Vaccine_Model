
data {
	int <lower = 0, upper = 1>	IncludeIntercept; // Boolean for intercept
	int <lower = 0, upper = 1>	IncludeScaling; 	// Boolean for scaling
	int <lower = 0, upper = 1>	DoKnots;        	// Boolean for spline
	int <lower = 0, upper = 1>	Quadratic;        // Boolean for quadratic (If 1, quadratic spline, if 0, linear)
	int <lower = 0, upper = 1>  DoVariants;       // Boolean for Variants model
	int <lower = 0, upper = 1>  DoAge;            // Boolean for age model
	
	int<lower = 1>	NumDatapoints; 	// Number of Rt values (all timepoints x regions) 
	int<lower = 1>	NumDoses; 		  // Number of parameters: Vax doses
	int<lower = 1>  NumVar;         // Number of parameters: SARS-CoV-2 variants
	int<lower = 1>  NumGroup;       // Number of parameters: age groups
	int<lower = 1>	NumLTLAs;		    // Number of regions / LTLAs
	int<lower = 1>	NumTimepoints;	// Number of timepoints (weeks)
	int<lower = 1>	NumKnots;       // Number of knots
	int<lower = 1>	NumPointsLine;  // Number of line-lambdas (knots x regions)
	int<lower = 1>  NumTrendPar;    // Number of NationalTrend pars (Knots or Free)
	
	int LTLAs[NumDatapoints];  // vector giving LTLA number for each Rt-region combination.
	int Groups[NumDatapoints]; // vector giving age group number for each Rt-region-age combi
	    
	vector[NumKnots]		Knots;	        // Sequence of knots
	vector[NumDatapoints]   Timepoints; // Sequence of timepoints
	
	vector[NumDatapoints] 			RtVals; 	    // y: Rt values across all time points and regions (expressed as giant vector) 
	matrix[NumDatapoints, NumDoses]	VaxProp;	// x predictor: Binary design matrix. Each row is a region and date combination.
	                                          // Each column has the proportion of vax at that time/LTLA with 1, 2, 3 doses
	matrix[NumDatapoints, NumVar] VarProp;    // x predictor: Each col has the proportion of variants 1, 2, etc. up to NumVar
	vector[NumDatapoints]       AgeProp;      // x predictor: Each age group proportion in each LTLA
}

transformed data{
	int<lower = 1> 	IntDim = 1;			// internal dimension of matrix factors - number of latent factors.
}

parameters {
	real<lower = 0> 				sigma_nc; 
	real<lower = 0> 				phi_nc;
	real<lower = 0> 				phi2_nc;
	real<lower = 0> 				phi3_nc;
	real<lower = 0> 				phi4_nc;
	
	matrix[NumLTLAs,IntDim] 		gamma_nc; 	  	// Prev alpha, indexed by: i) LTLA; ii) factor 
	vector[NumLTLAs] 				intercept_nc;     	// indexed by: i) LTLA
	
  matrix[NumTrendPar,IntDim] 		lambda_raw_nc; // indexed by: i) Knots/Timepoints, ii) factor
  
	matrix<lower = 0>[NumDoses*NumGroup, NumVar] 	VaxEffect_nc;	
	matrix<lower = 0>[1, NumVar] VarAdvantage_nc;
}

transformed parameters{
	real sigma 	= 0;
	real phi 	= 0;
	real phi2 	= 0;
	real phi3 	= 0;
	real phi4 	= 0;
	
	// allocate
	vector[NumLTLAs] 			intercept 	= rep_vector(0, NumLTLAs); 
	matrix[NumLTLAs,IntDim] 	gamma 		= rep_matrix(0, NumLTLAs, IntDim); 

	matrix[NumTrendPar,IntDim] lambda_raw 		= rep_matrix(0, NumTrendPar, IntDim); // has NumKnots/NumTimepoints rows (not NumDatapoints)
	matrix[NumPointsLine,IntDim] lambda 	= rep_matrix(0, NumPointsLine, IntDim);  // has NumKnots*NumLTLAs rows (not NumTimepoints)
	
	matrix[NumKnots-1, IntDim] origin; 	// intercept
	matrix[NumKnots-1, IntDim] slope;  	// slope
	matrix[NumKnots-2, IntDim] a; 		// Coeff a for quadratic eq
	matrix[NumKnots-2, IntDim] b; 		// Coeff b for quadratic eq
	matrix[NumKnots-2, IntDim] c; 		// Coeff c for quadratic eq
	
	matrix[NumDatapoints, IntDim]	NationalTrend	= rep_matrix(0, NumDatapoints, IntDim); // To calculate lambda from line or free
	vector[NumDatapoints] RegionalTrends 		= rep_vector(0, NumDatapoints);
	
	matrix[NumDatapoints, NumVar] VacEffects_Regional 	= rep_matrix(0, NumDatapoints, NumVar);
	matrix[NumDoses*NumGroup, NumVar] 			VaxEffect 	= rep_matrix(0, NumDoses*NumGroup, NumVar);
	
	matrix[1, NumVar] VarAdvantage = rep_matrix(1, 1, NumVar);
	
	vector[NumDatapoints] LogPredictions 		= rep_vector(0, NumDatapoints);
	
	// initialize - get centered parameter values from their non-centered equivalents.
	phi  		= phi_nc		* 2.0;
	phi2 		= phi2_nc		* 0.5;
	phi3 		= phi3_nc		* 0.5;
	phi4 		= phi4_nc		* 0.5;
	
	sigma 		= sigma_nc		* 0.5;
	gamma 		= gamma_nc 		* phi2;
	intercept 	= intercept_nc	* phi3;
	VaxEffect 	= VaxEffect_nc	* phi;
	
	if(DoVariants) {
	 VarAdvantage[,2:NumVar] 	= VarAdvantage_nc[,2:NumVar]	* phi;
	}
	
	if (DoKnots) {  
	 
		// Lambda for the SPLINE
		lambda_raw 	= lambda_raw_nc 	* phi4;
		{  
			// initialize lambda matrix to have same values for every LTLA
			int ind = 0; // initialize index
			
			for (i in 1:NumLTLAs){
				for(j in 1:IntDim)
					lambda[(ind + 1):(ind + NumKnots), j] = lambda_raw[1:NumKnots, j];
		    
				ind = ind + NumKnots; // update index
			}
		}
		
		// spline to fit the lambda: Quadratic or Linear
		if (Quadratic) {
		
			for (i in 1:(NumKnots-2))
				for (j in 1:IntDim) {
					a[i,j] = (lambda[i+2, j] - lambda[i+1, j] - ((Knots[i+2] - Knots[i+1]) * (lambda[i, j] - lambda[i+1, j]) / (Knots[i] - Knots[i+1]))) / ((Knots[i+2]*Knots[i+2]) - (Knots[i+1]*Knots[i+1]) - (Knots[i] + Knots[i+1])*(Knots[i+2] - Knots[i+1]));
					b[i,j] = ((lambda[i, j] - lambda[i+1, j]) / (Knots[i] - Knots[i+1])) - (Knots[i] + Knots[i+1])*a[i,j];
					c[i,j] = lambda[i, j] - (b[i,j] * Knots[i]) - (a[i,j] * (Knots[i]*Knots[i]));
				}
		    
		} else {
		
			for (i in 1:(NumKnots-1))
				for (j in 1:IntDim) {
					slope[i,j] = (lambda[i+1, j] - lambda[i, j])/(Knots[i+1] - Knots[i]);
					origin[i,j] = lambda[i, j] - slope[i, j]*Knots[i]; 
				}
		}
		
		// Calculate lambda from the spline ONLY IF we are doing knots
		if (Quadratic) {
		
		for (i in 1:NumDatapoints) 
			for (j in 1:IntDim)
		  		if (LTLAs[i] == 1) {
			        
			        if (Timepoints[i] >= Knots[NumKnots-1] && Timepoints[i] < Knots[NumKnots]) {
			              NationalTrend[i,j] = a[NumKnots-2,1]*(Timepoints[i] * Timepoints[i]) + b[NumKnots-2,1]*Timepoints[i] + c[NumKnots-2,1]; 
			        
			        } else for (k in 1:(NumKnots-2))
			            if (Timepoints[i] >= Knots[k] && Timepoints[i] < Knots[k+1])
			              NationalTrend[i,j] = a[k,j]*(Timepoints[i] * Timepoints[i]) + b[k,j]*Timepoints[i] + c[k,j];
			
		  		} else NationalTrend[i,j] = NationalTrend[(i - NumTimepoints),j]; // i.e. make equal to previous LTLA's NationalTrend 
		  
		} else { // i.e. doing knots, linear spline 
   
			for (i in 1:NumDatapoints)
				for (j in 1:IntDim)
					if (LTLAs[i] == 1) {     
						for (k in 1:(NumKnots-1)) 
							if (Timepoints[i] >= Knots[k] && Timepoints[i] < Knots[k+1])
								NationalTrend[i,j] = origin[k,j] + slope[k,j] * Timepoints[i];
   
					} else NationalTrend[i,j] = NationalTrend[(i - NumTimepoints),j]; // i.e. make equal to previous LTLA's NationalTrend 
		}

	} else { // i.e. step function (each week a free parameter)
 
		// Initialise lambda if we are not doing knots: free parameters allowed
		lambda_raw 	= lambda_raw_nc 	* phi4;
		{  
			// initialize lambda matrix to have same values for every LTLA 
			int ind_par = 0; // initialize index
			for (i in 1:NumLTLAs){
		    
				for(j in 1:IntDim)
					NationalTrend[(ind_par + 1):(ind_par + NumTimepoints), j] = lambda_raw[1:NumTimepoints, j];
				ind_par = ind_par + NumTimepoints; // update index
			}
		}
 	}

 // Calculate regional trend from national trend
	for (i in 1:NumDatapoints){
		for(j in 1:IntDim){
			if (IncludeIntercept) {
			
				if (IncludeScaling) {
					RegionalTrends[i] += NationalTrend[i,j] * gamma[LTLAs[i],j] + intercept[LTLAs[i]]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
				} else {
					RegionalTrends[i] += NationalTrend[i,j] + intercept[LTLAs[i]]; }
			
			} else {
			
				if (IncludeScaling) {
					RegionalTrends[i] += NationalTrend[i,j] * gamma[LTLAs[i],j]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
				} else {
					RegionalTrends[i] += NationalTrend[i,j]; }
			}
		}
	}
	
	
	// CONTINUE WITH LOG PRED MODEL WHETHER LAMBDA WAS FIT IN THE LINE OR NOT
	// VacEffects_Regional[1:NumDatapoints] = VaxProp * VaxEffect; // x * -beta in manuscript
	for (k in 1:NumGroup)
	  for (i in 1:NumDatapoints)
	    for(l in 1:NumVar)
	    {
		    VacEffects_Regional[i, l] = 0; // initialize to zero for each timepoint
		     if(Groups[i] == k) {
		        for (j in 1:NumDoses) {
			        VacEffects_Regional[i, l] += VaxProp[i,j] * VaxEffect[j, l] * AgeProp[i];
		      }
		    }
	    }
	
	// final (logged) regional Rt predictions are regional trends minus regional vaccine effects
	// LogPredictions[1:NumDatapoints] = RegionalTrends[1:NumDatapoints] - VacEffects_Regional[1:NumDatapoints];
	for (i in 1:NumDatapoints)
	  for (j in 1:NumVar)
		  LogPredictions[i] = VarProp[i, j]*(RegionalTrends[i]*VarAdvantage[1, j] - VacEffects_Regional[i, j]);
	
}

model {
	
	for (i in 1:NumLTLAs)
		for(j in 1:IntDim)
			gamma_nc[i,j] ~ std_normal();
	
	for (i in 1:NumTrendPar)
		for(j in 1:IntDim)
			lambda_raw_nc[i,j] ~ std_normal(); 
	
	intercept_nc 	~ std_normal();
	
	phi_nc 			~ std_normal();
	phi2_nc 		~ std_normal();
	phi3_nc 		~ std_normal();
	phi4_nc 		~ std_normal();
	sigma_nc 		~ std_normal();
	
	for (i in 1:NumDoses*NumGroup)
	  for (j in 1:NumVar)
		  VaxEffect_nc[i, j] ~ std_normal();
		  
	for (i in 2:NumVar)
		VarAdvantage_nc[,i] ~ std_normal();
		
	RtVals 			~ normal(LogPredictions, sigma);
}

generated quantities {
	vector[NumDatapoints] log_lik;
	
	for (i in 1:NumDatapoints) {
		log_lik[i] = normal_lpdf(RtVals[i] | LogPredictions[i], sigma); // NOTES: Log of normal function: real normal_lpdf(reals y | reals mu, reals sigma)
	}
}
