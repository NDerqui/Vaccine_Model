
data {
  int <lower = 0, upper = 1>		IncludeIntercept; 	// Boolean for intercept
  int <lower = 0, upper = 1>		IncludeScaling; 	// Boolean for scaling
  int <lower = 0, upper = 1>		DoKnots;        	// Boolean for spline
  int <lower = 0, upper = 1>		Quadratic;        	// Boolean for quadratic (If 1, quadratic spline, if 0, linear)

  int<lower = 1> 		NumDatapoints; 	// Number of Rt values (all timepoints x regions) 
  int<lower = 1> 		NumDoses; 		// Number of parameters: Vax doses
  int<lower = 1> 		NumLTLAs;		// Number of regions / LTLAs
  int<lower = 1> 		NumTimepoints;	// Number of timepoints (weeks)
  int<lower = 1>        NumKnots;       // Number of knots
  int<lower = 1>        NumPointsLine;  // Number of line-lambdas (knots x regions)
  
  vector[NumDatapoints] 			RtVals; 	// y: Rt values across all time points and regions (expressed as giant vector) 
  matrix[NumDatapoints,NumDoses] 	VaxProp;	// x predictor: Binary design matrix. Each row is a region and date combination.
      // Each column has the proportion of vax at that time/LTLA with 1, 2, 3 doses
  int 								LTLAs[NumDatapoints];  // vector giving LTLA number for each Rt-region combination.
      
  vector[NumKnots]            Knots;	  // Sequence of knots
  vector[NumDatapoints]       Timepoints; // Sequence of timepoints
}

transformed data{
  int<lower = 1> 	IntDim = 1;			// internal dimension of matrix factors - number of latent factors.
}

parameters {
  matrix[NumLTLAs,IntDim] 		gamma_nc; 	  		// Prev alpha, indexed by: i) LTLA; ii) factor 
  vector[NumLTLAs] 				intercept_nc;     	// indexed by: i) LTLA

  //// looks wrong to have two things here - your parameters are only the knots - not two versions of the knots.
  matrix[NumKnots,IntDim] 		lambda_raw_nc; 		// indexed by: i) Knots, ii) factor
  matrix[NumTimepoints,IntDim] 	lambda_raw_nc_par; 	// indexed by: i) Timepoint, ii) factor

  vector<lower = 0>[NumDoses] 	VaxEffect_nc;
  real<lower = 0> 		sigma_nc; 
  real<lower = 0> 		phi_nc;
  real<lower = 0> 		phi2_nc;
  real<lower = 0> 		phi3_nc;
  real<lower = 0> 		phi4_nc;
}

transformed parameters{

  // allocate
  vector[NumLTLAs] 			intercept 	= rep_vector(0, NumLTLAs); 
  matrix[NumLTLAs,IntDim] 	gamma 		= rep_matrix(0, NumLTLAs, IntDim); 
  vector[NumDoses] 			VaxEffect 	= rep_vector(0, NumDoses);
  real sigma 	= 0;
  real phi 		= 0;
  real phi2 	= 0;
  real phi3 	= 0;
  real phi4 	= 0;
  vector[NumDatapoints] RegionalTrends 		= rep_vector(0, NumDatapoints);
  vector[NumDatapoints] fixed_effects 		= rep_vector(0, NumDatapoints);
  vector[NumDatapoints] LogPredictions 		= rep_vector(0, NumDatapoints);

  matrix[NumKnots,IntDim] lambda_raw 		= rep_matrix(0, NumKnots, IntDim); 		// has NumKnots rows (not NumDatapoints)
  matrix[NumPointsLine,IntDim] lambda 		= rep_matrix(0, NumPointsLine, IntDim); // has NumKnots*NumLTLAs rows (not NumTimepoints)
  matrix[NumKnots-1, IntDim] origin; 						//intercept
  matrix[NumKnots-1, IntDim] slope;  						//slope
  matrix[NumKnots-2, IntDim] a; 							//Coeff a for quadratic eq
  matrix[NumKnots-2, IntDim] b; 							//Coeff b for quadratic eq
  matrix[NumKnots-2, IntDim] c; 							//Coeff c for quadratic eq


  matrix[NumTimepoints,IntDim] 	lambda_raw_par 		= rep_matrix(0, NumTimepoints, IntDim); // has NumTimepoint rows (not NumDatapoints)
  matrix[NumDatapoints, IntDim] lambda_parameters 	= rep_matrix(0, NumDatapoints, IntDim); // To calculate lambda from line

  // initialize - get centered parameter values from their non-centered equivalents.
  phi  			= phi_nc			* 2.0;
  phi2 			= phi2_nc			* 0.5;
  phi3 			= phi3_nc			* 0.5;
  phi4 			= phi4_nc			* 0.5;
  sigma 		= sigma_nc			* 0.5;
  VaxEffect 	= VaxEffect_nc 	* phi;
  gamma 		= gamma_nc 		* phi2;
  intercept 	= intercept_nc	* phi3;
  
  //Lambda for the SPLINE
  lambda_raw 	= lambda_raw_nc 	* phi4;
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
  
 if (DoKnots) {  
  // spline to fit the lambda
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
  
  }

  //Calculate lambda from the spline ONLY IF we are doing knots
 if (DoKnots) {
 
   if (Quadratic) {

     // very inefficient, in that you're calculating x and y coordinates for spline for every LTLA, when before vaccination they are all the same. 
     // Make the calculation once for the spline, then populate for every LTLA. See (and check) below 
     for (i in 1:NumDatapoints) 
     	for (j in 1:IntDim)
       		if (LTLAs[i] == 1) {
 		        
 		        if (Timepoints[i] >= Knots[NumKnots-1] && Timepoints[i] < Knots[NumKnots]) {
		              lambda_parameters[i,j] = a[NumKnots-2,1]*(Timepoints[i] * Timepoints[i]) + b[NumKnots-2,1]*Timepoints[i] + c[NumKnots-2,1]; 
		        
		        } else for (k in 1:(NumKnots-2))
		            if (Timepoints[i] >= Knots[k] && Timepoints[i] < Knots[k+1])
		              lambda_parameters[i,j] = a[k,j]*(Timepoints[i] * Timepoints[i]) + b[k,j]*Timepoints[i] + c[k,j];
     	
       		} else lambda_parameters[i,j] = lambda_parameters[(i - NumTimepoints),j]; // i.e. make equal to previous LTLA's lambda_parameters 
       
   } else { // i.e. doing knots, linear spline 
   
     // very inefficient, in that you're calculating x and y coordinates for spline for every LTLA, when before vaccination they are all the same. 
     // Make the calculation once for the spline, then populate for every LTLA. See (and check) below 
     for (i in 1:NumDatapoints)
      for (j in 1:IntDim)
		if (LTLAs[i] == 1) {     
        for (k in 1:(NumKnots-1)) 
            if (Timepoints[i] >= Knots[k] && Timepoints[i] < Knots[k+1])
              lambda_parameters[i,j] = origin[k,j] + slope[k,j] * Timepoints[i];
        
        } else lambda_parameters[i,j] = lambda_parameters[(i - NumTimepoints),j]; // i.e. make equal to previous LTLA's lambda_parameters 
   }
        
 } else { // i.e. step function (each week a free parameter)
 
   // Initialise lambda if we are not doing knots: free parameters allowed
  lambda_raw_par 	= lambda_raw_nc_par 	* phi4;
  {  
  	// initialize lambda matrix to have same values for every LTLA 
    int ind_par = 0; // initialize index
    for (i in 1:NumLTLAs){
      
      for(j in 1:IntDim)
        lambda_parameters[(ind_par + 1):(ind_par + NumTimepoints), j] = lambda_raw_par[1:NumTimepoints, j];
     
      ind_par = ind_par + NumTimepoints; // update index
    }
  }
 }

  // CONTINUE WITH LOG PRED MODEL WHETHER LAMBDA WAS FIT IN THE LINE OR NOT
  fixed_effects[1:NumDatapoints] = VaxProp * - VaxEffect; // x * -beta in manuscript
  
  for (i in 1:NumDatapoints){
    for(j in 1:IntDim){
      
      if (IncludeIntercept) {
       
        if (IncludeScaling) {
          RegionalTrends[i] += lambda_parameters[i,j] * gamma[LTLAs[i],j] + intercept[LTLAs[i]]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
        } else {
          RegionalTrends[i] += lambda_parameters[i,j] + intercept[LTLAs[i]]; }
        
        } else {
          
        if (IncludeScaling) {
          RegionalTrends[i] += lambda_parameters[i,j] * gamma[LTLAs[i],j]; // lambda * gamma^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda
	      } else {
          RegionalTrends[i] += lambda_parameters[i,j]; }
        }

        }
      }
  LogPredictions[1:NumDatapoints] = fixed_effects[1:NumDatapoints] + RegionalTrends[1:NumDatapoints];

}

model {
  
  VaxEffect_nc ~ std_normal();
  
  for (i in 1:NumLTLAs){
    for(j in 1:IntDim){
      gamma_nc[i,j] ~ std_normal();
    }
  }
  
  for (i in 1:NumKnots){
   for(j in 1:IntDim){
    lambda_raw_nc[i,j] ~ std_normal();
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
  phi4_nc 		~ std_normal();
  sigma_nc 		~ std_normal();
  RtVals 		~ normal(LogPredictions, sigma);
  
}

generated quantities {
  vector[NumDatapoints] log_lik;
  for (i in 1:NumDatapoints) {
   log_lik[i] = normal_lpdf(RtVals[i] | LogPredictions[i], sigma); // NOTES: Log of normal function: real normal_lpdf(reals y | reals mu, reals sigma)
  }
}
