
data {
  int <lower = 0, upper = 1>					IncludeIntercept; 	// Boolean for intercept
  int <lower = 0, upper = 1>					IncludeScaling; 	// Boolean for scaling

  int<lower = 1> 							NumDatapoints; 	// Number of Rt values (all timepoints x regions) 
  int<lower = 1> 							NumDoses; 			// Number of parameters: Vax doses
  int<lower = 1> 							NumLTLAs;				// Number of regions / LTLAs
  int<lower = 1> 							NumTimepoints;	// Number of timepoints (weeks)
  int<lower = 1>              NumKnots;       // Number of knots
  
  vector[NumDatapoints] 				RtVals; 				
      // y: Rt values across all time points and regions (expressed as giant vector) 
  matrix[NumDatapoints,NumDoses] 	VaxProp; 			
      // x predictor: Binary design matrix. Each row is a region and date combination.
      // Each column has the proportion of vax at that time/LTLA with 1, 2, 3 doses
  int 									LTLAs[NumDatapoints]; 	
      // vector giving LTLA number for each Rt-region combination.
      
  vector[NumKnots]            Knots;
      // Sequence of knots
}

transformed data{
  int<lower = 1> 	IntDim = 1;			// internal dimension of matrix factors - number of latent factors.
  int<lower = 1>              NumLines = NumKnots - 1;
}

parameters {
  matrix[NumLTLAs,IntDim] 			gamma_nc; 	  	// Prev alpha, indexed by: i) LTLA; ii) factor 
  vector[NumLTLAs] 					intercept_nc;     	// indexed by: i) LTLA
  matrix[NumTimepoints,IntDim] 		lambda_raw_nc; // indexed by: i) DataPoints, ii) factor
  vector<lower = 0>[NumDoses] 	VaxEffect_nc;
  real<lower = 0> 		phi_nc;
  real<lower = 0> 		sigma_nc;
  real<lower = 0> 		phi2_nc;
  real<lower = 0> 		phi3_nc;
  
  matrix[NumLines, IntDim] origin; //intercept
  matrix[NumLines, IntDim] slope;  //slope
}

transformed parameters{

  // allocate
  vector[NumLTLAs] 			intercept 	= rep_vector(0, NumLTLAs); 
  matrix[NumLTLAs,IntDim] 	gamma 		= rep_matrix(0, NumLTLAs, IntDim); 
  vector[NumDoses] 	VaxEffect = rep_vector(0, NumDoses);
  real phi 		= 0;
  real sigma 	= 0;
  real phi2 	= 0;
  real phi3 	= 0;
  vector[NumDatapoints] random_effects 		= rep_vector(0, NumDatapoints);
  vector[NumDatapoints] fixed_effects 		= rep_vector(0, NumDatapoints);
  vector[NumDatapoints] LogPredictions 		= rep_vector(0, NumDatapoints);
  matrix[NumTimepoints,IntDim] lambda_raw 	= rep_matrix(0, NumTimepoints, IntDim); // has NumTimepoints rows (not NumDatapoints)
  matrix[NumDatapoints,IntDim] lambda 		= rep_matrix(0, NumDatapoints, IntDim); // has NumDatapoints rows (not NumTimepoints)
  
  matrix[NumKnots,IntDim] lambda_y;
  matrix[NumLines,IntDim] lambda_change;
  
  // initialize
  phi  			= phi_nc			* 2.0;
  phi2 			= phi2_nc			* 0.5;
  phi3 			= phi3_nc			* 0.5;
  sigma 		= sigma_nc			* 0.5;
  VaxEffect 	= VaxEffect_nc 	* phi;
  gamma 		= gamma_nc 			* phi2;
  intercept 	= intercept 		* phi3;
  lambda_raw 	= lambda_raw_nc 	* phi3;
  {  
  	// initialize lambda matrix to have same values for every LTLA
    int ind = 0; // initialize index
    for (i in 1:NumLTLAs){
      for(j in 1:IntDim){
        lambda[(ind + 1):(ind + NumTimepoints), j] = lambda_raw[1:NumTimepoints, j];
      }
      ind = ind + NumTimepoints; // update index
    }
  }
  
  
  fixed_effects[1:NumDatapoints] = VaxProp * - VaxEffect; // x * -beta in manuscript
  
  for (i in 1:NumDatapoints){
    for(j in 1:IntDim){
      
      if (IncludeIntercept) {
       
        if (IncludeScaling) {
          random_effects[i] += lambda[i,j] * gamma[LTLAs[i],j] + intercept[LTLAs[i]]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
        } else {
          random_effects[i] += lambda[i,j] + intercept[LTLAs[i]]; }
        
        } else {
          
        if (IncludeScaling) {
          random_effects[i] += lambda[i,j] * gamma[LTLAs[i],j]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
	      } else {
          random_effects[i] += lambda[i,j]; }
        }

        }
      }
  LogPredictions[1:NumDatapoints] = fixed_effects[1:NumDatapoints] + random_effects[1:NumDatapoints];

  lambda_y = lambda[1:NumKnots];

  for (i in 2:(NumKnots)){
    for (j in 1:NumLines) {
      lambda_change[j] = lambda_y[i] - lambda_y[i-1]; 
      lambda_change[j] = origin[j] + slope[j]*(Knots[i] - Knots[i-1]);
    }
  }
}

model {
  VaxEffect_nc ~ std_normal();
  for (i in 1:NumLTLAs){
    for(j in 1:IntDim){
      gamma_nc[i,j] ~ std_normal();
    }
  }
  for (i in 1:NumTimepoints){
    for(j in 1:IntDim){
      lambda_raw_nc[i,j] ~ std_normal();
    }
  }  
  intercept_nc 	~ std_normal();
  phi_nc 		~ std_normal();
  phi2_nc 		~ std_normal();
  phi3_nc 		~ std_normal();
  sigma_nc 		~ std_normal();
  RtVals 		~ normal(LogPredictions, sigma);
  
}

generated quantities {
  vector[NumDatapoints] log_lik;
  for (i in 1:NumDatapoints) {
   log_lik[i] = normal_lpdf(RtVals[i] | LogPredictions[i], sigma); // NOTES: Log of normal function: real normal_lpdf(reals y | reals mu, reals sigma)
  }
}
