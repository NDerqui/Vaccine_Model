
data {
  int <lower = 0, upper = 1>					IncludeIntercept; 	// Boolean for intercept
  int <lower = 0, upper = 1>					IncludeScaling; 	// Boolean for scaling

  int<lower = 1> 							NumDatapoints; 	// Number of Rt values (all timepoints x regions) 
  int<lower = 1> 							NumDoses; 			// Number of parameters: Vax doses
  int<lower = 1> 							NumLTLAs;				// Number of regions / LTLAs
  int<lower = 1> 							NumTimepoints;	// Number of timepoints (weeks)
  int<lower = 1>              NumKnots;       // Number of knots
  int<lower = 1>              NumPointsLine;  // Number of lambdas (knots x regions)
  
  vector[NumDatapoints] 				RtVals; 				
      // y: Rt values across all time points and regions (expressed as giant vector) 
  matrix[NumDatapoints,NumDoses] 	VaxProp; 			
      // x predictor: Binary design matrix. Each row is a region and date combination.
      // Each column has the proportion of vax at that time/LTLA with 1, 2, 3 doses
  int 									LTLAs[NumDatapoints]; 	
      // vector giving LTLA number for each Rt-region combination.
      
  vector[NumKnots]            Knots;
      // Sequence of knots
  vector[NumDatapoints]       Timepoints;
      // Sequence of timepoints
}

transformed data{
  int<lower = 1> 	IntDim = 1;			// internal dimension of matrix factors - number of latent factors.
}

parameters {
  matrix[NumLTLAs,IntDim] 			gamma_nc; 	  	// Prev alpha, indexed by: i) LTLA; ii) factor 
  vector[NumLTLAs] 					intercept_nc;     	// indexed by: i) LTLA
  matrix[NumKnots,IntDim] 		lambda_raw_nc; // indexed by: i) Knots, ii) factor
  matrix[NumTimepoints,IntDim] 		lambda_raw_nc_par; // indexed by: i) DataPoints, ii) factor
  vector<lower = 0>[NumDoses] 	VaxEffect_nc;
  real<lower = 0> 		phi_nc;
  real<lower = 0> 		sigma_nc;
  real<lower = 0> 		phi2_nc;
  real<lower = 0> 		phi3_nc;
  
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

  matrix[NumKnots,IntDim] lambda_raw 	= rep_matrix(0, NumKnots, IntDim); // has NumKnots rows (not NumDatapoints)
  matrix[NumPointsLine,IntDim] lambda 		= rep_matrix(0, NumPointsLine, IntDim); // has NumKnots*NumLTLAs rows (not NumTimepoints)
  matrix[NumKnots-1, IntDim] origin; //intercept
  matrix[NumKnots-1, IntDim] slope;  //slope

  matrix[NumTimepoints,IntDim] lambda_raw_par 	= rep_matrix(0, NumTimepoints, IntDim); // has NumTimepoints rows (not NumDatapoints)
  matrix[NumDatapoints, IntDim] lambda_parameters 		= rep_matrix(0, NumDatapoints, IntDim); // To calculate lambda from line

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
        lambda[(ind + 1):(ind + NumKnots), j] = lambda_raw[1:NumKnots, j];
      }
      ind = ind + NumKnots; // update index
    }
  }
  
  // line to fit the lambda
  for (i in 1:(NumKnots-1)){
    slope[i] = (lambda[i+1] - lambda[i])/(Knots[i+1] - Knots[i]);
    origin[i] = lambda[i] - slope[i]*Knots[i];
  }

  // initialise lambda from line
  lambda_raw_par 	= lambda_raw_nc_par 	* phi3;
  {  
  	int ind_2 = 0; // initialize index
    for (i in 1:NumLTLAs){
      for(j in 1:IntDim){
        lambda_parameters[(ind_2 + 1):(ind_2 + NumTimepoints), j] = lambda_raw_par[1:NumTimepoints, j];
      }
      ind_2 = ind_2 + NumTimepoints; // update index
    }
  }

  for (i in 1:NumDatapoints){
    if(Timepoints[i] < 6) {
      lambda_parameters[i] = origin[1] + slope[1]*Timepoints[i];
    } else {
      if (Timepoints[i] < 12) {
        lambda_parameters[i] = origin[2] + slope[2]*Timepoints[i];
      } else {
        if (Timepoints[i] < 16) {
          lambda_parameters[i] = origin[3] + slope[3]*Timepoints[i];
        } else {
          if (Timepoints[i] < 25) {
            lambda_parameters[i] = origin[4] + slope[4]*Timepoints[i];
          } else {
            if (Timepoints[i] < 42) {
              lambda_parameters[i] = origin[5] + slope[5]*Timepoints[i];
            }
          }
        }
      }
    }
  }
  
  fixed_effects[1:NumDatapoints] = VaxProp * - VaxEffect; // x * -beta in manuscript
  
  for (i in 1:NumDatapoints){
    for(j in 1:IntDim){
      
      if (IncludeIntercept) {
       
        if (IncludeScaling) {
          random_effects[i] += lambda_parameters[i,j] * gamma[LTLAs[i],j] + intercept[LTLAs[i]]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
        } else {
          random_effects[i] += lambda_parameters[i,j] + intercept[LTLAs[i]]; }
        
        } else {
          
        if (IncludeScaling) {
          random_effects[i] += lambda_parameters[i,j] * gamma[LTLAs[i],j]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
	      } else {
          random_effects[i] += lambda_parameters[i,j]; }
        }

        }
      }
  LogPredictions[1:NumDatapoints] = fixed_effects[1:NumDatapoints] + random_effects[1:NumDatapoints];

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
      lambda_raw_nc_par[i,j] ~ std_normal();
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
