
data {
	int <lower = 0, upper = 1>	IncludeIntercept; // Boolean for intercept
	int <lower = 0, upper = 1>	IncludeScaling; 	// Boolean for scaling
	int <lower = 0, upper = 1>	DoKnots;        	// Boolean for spline
	int <lower = 0, upper = 1>	Quadratic;        // Boolean for quadratic (If 1, quadratic spline, if 0, linear)
	int <lower = 0, upper = 1>  DoVariants;       // Boolean for Variants model
	int <lower = 0, upper = 1>  DoVaxVariants;    // Boolean for specific variant's VEs
	int <lower = 0, upper = 1>  DoAge;            // Boolean for age model
	int <lower = 0, upper = 1>  DoVaxAge;         // Boolean for specific age's VEs
	
	int<lower = 1>	NumDatapoints; 	// Number of Rt values (all timepoints x regions) 
	int<lower = 1>	NumDoses; 		  // Number of parameters: Vax doses
	int<lower = 1>  NumVar;         // Number of parameters: SARS-CoV-2 variants
	int<lower = 1>  NumVaxVar;      // Number of parameters: variants considered for the VEs
	int<lower = 1>  NumGroup;       // Number of parameters: age groups
	int<lower = 1>  NumVaxGroup;    // Number of parameters: agr groups considered for the VEs
	int<lower = 1>	NumLTLAs;		    // Number of regions / LTLAs
	int<lower = 1>	NumTimepoints;	// Number of timepoints (weeks)
	int<lower = 1>	NumKnots;       // Number of knots
	int<lower = 1>	NumPointsLine;  // Number of line-lambdas (knots x regions)
	int<lower = 1>  NumTrendPar;    // Number of NationalTrend pars (Knots or Free)
	
	int LTLAs[NumDatapoints];  // vector giving LTLA number for each Rt-region combination.
	int Groups[NumDatapoints]; // vector giving age group number for each Rt-region-age combi
	    
	vector[NumKnots]		Knots;	        // Sequence of knots
	vector[NumDatapoints]   Timepoints; // Sequence of timepoints
	
	vector <lower = 0> [NumDatapoints] 			                RtVals;   // y: Rt values across all time points and regions (expressed as giant vector) 
	// matrix <lower = 0, upper = 1> [NumDatapoints, NumDoses]	VaxProp;	// x predictor: Binary design matrix. Each row is a region and date combination.
	// instead of the above, try this.
	real<lower = 0, upper = 1> VaxProp[NumDatapoints, NumDoses, NumGroup];
	                                                                  // Each column has the proportion of vax at that time/LTLA with 1, 2, 3 doses
	matrix <lower = 0, upper = 1> [NumDatapoints, NumVar]   VarProp;  // x predictor: Each col has the proportion of variants 1, 2, etc. up to NumVar
	//matrix <lower = 0, upper = 1> [NumDatapoints, 1]        AgeProp;  // x predictor: Each age group proportion in each LTLA
	// Two changes here - first, make AgeProp a vector, not a matrix.
	// Second, it only needs a few numbers in it, namely the number of age groups, not the number of data points (i.e. LTLAs x Timepoints)
	vector <lower = 0, upper = 1> [NumGroup] AgeProp;  // x predictor: Each age group proportion in each LTLA
}


parameters {
	real<lower = 0> 				sigma_nc; 
	real<lower = 0> 				phi_nc;
	real<lower = 0> 				phi2_nc;
	real<lower = 0> 				phi3_nc;
	real<lower = 0> 				phi4_nc;
	
	vector <lower = 0> [NumLTLAs] 	RegionalScale_nc; 	  	// Prev alpha, indexed by: i) LTLA; 
	vector[NumLTLAs] 				intercept_nc;     	// indexed by: i) LTLA
	
	vector <lower = 0> [NumTrendPar] 		lambda_raw_nc; // indexed by: i) Knots/Timepoints, 
  
  //matrix <lower = 0, upper = 1> [NumDoses, NumVaxVar*NumVaxGroup] 	VaxEffect_nc;
  // instead of the above, try this.
	real<lower = 0, upper = 1> VaxEffect_nc[NumDoses, NumVaxVar, NumVaxGroup];


	vector <lower = 1> [NumVar-1] VarAdvantage_nc; /// move this to data while debugging using Knock Whittles values.
}

transformed parameters{
  
	real sigma 	= 0;
	real phi 	  = 0;
	real phi2 	= 0;
	real phi3 	= 0;
	real phi4 	= 0;
	real FinalRtperVariantTimeRegion 	= 0;

	// allocate
	vector[NumLTLAs] 	intercept 	  = rep_vector(0, NumLTLAs); 
	vector[NumLTLAs] 	RegionalScale = rep_vector(1, NumLTLAs); 

	vector <lower = 0> [NumTrendPar] lambda_raw = rep_vector(0, NumTrendPar); 	// has NumKnots/NumTimepoints rows (not NumDatapoints)
	vector <lower = 0> [NumPointsLine] lambda 	= rep_vector(0, NumPointsLine);	// has NumKnots*NumLTLAs rows (not NumTimepoints)
	
	vector[NumKnots-1] origin; 	// intercept
	vector[NumKnots-1] slope;  	// slope
	vector[NumKnots-2] a; 		// Coeff a for quadratic eq
	vector[NumKnots-2] b; 		// Coeff b for quadratic eq
	vector[NumKnots-2] c; 		// Coeff c for quadratic eq
	
	vector [NumDatapoints] NationalTrend	= rep_vector(0, NumDatapoints); // To calculate lambda from line or free
	vector [NumDatapoints] RegionalTrends 	= rep_vector(0, NumDatapoints);
	
	real<lower = 0, upper = 1> VaxEffect[NumDoses, NumVaxVar, NumVaxGroup];
	vector <lower = 1> [NumVar] VarAdvantage;
	
	vector [NumDatapoints] LogPredictions 		= rep_vector(0, NumDatapoints);
	
	// VarAdvantage for base variant set to 1
	VarAdvantage[1] = 1;
	
	// initialize - get centered parameter values from their non-centered equivalents.
	phi  		= phi_nc		* 2.0;
	phi2 		= phi2_nc		* 0.5;
	phi3 		= phi3_nc		* 0.5;
	phi4 		= phi4_nc		* 0.5;
	
	sigma 		= sigma_nc		* 0.5;
	RegionalScale 		= RegionalScale_nc 		* phi2;
	intercept 	= intercept_nc	* phi3;

  // VaxEffect 	= VaxEffect_nc	* phi; /// Danny wants to get rid of this. // reckon fine to do it by entire array, but can do indicies version below if that doesn't work.
  for (VaxVariant in 1:NumVaxVar)
    for (VaxGroup in 1:NumVaxGroup)
      for (Dose in 1:NumDoses)
        VaxEffect[Dose, VaxVariant, VaxGroup] = VaxEffect_nc[Dose, VaxVariant, VaxGroup] * phi; 
  


//	if (NumVaxVar == NumVar && NumGroup == NumVaxGroup) {
//	  
//	  VaxEffect 	= VaxEffect_nc	* phi; /// Danny wants to get rid of this
//	}
//	   else { //// NumVaxVar != NumVar OR NumGroup != NumVaxGroup (or both)
//	     
//	     if (NumGroup == NumVaxGroup) { // i.e. NumVaxVar
//	       
//	       for (Variant in 1:NumVar) // i.e. NumVaxVar < NumVar (should be NumVaxVar = 1 and NumVar = 1,2,3, or 4 depending on which variants we're modelling)
//          for (Dose in 1:NumDoses)  
//	            VaxEffect[Dose, Variant] = VaxEffect_nc[Dose, 1] * phi;
//	     }
//	       else if (NumVaxVar == NumVar){
//	         
//	         for (Group in 1:NumGroup) // i.e. NumVaxGroup < NumGroup (should be NumVaxVar = 1 and NumVar = 1,2,3, or 4 depending on which variants we're modelling)
//	          for (Dose in 1:NumDoses)  
//	            VaxEffect[Dose, Group] = VaxEffect_nc[Dose, 1] * phi;
//	         
//	       } else { // i.e NumVaxVar != NumVar AND NumGroup != NumVaxGroup
//	       
//	       for (Variant in 1:NumVar) // i.e. NumVaxVar < NumVar (should be NumVaxVar = 1 and NumVar = 1,2,3, or 4 depending on which variants we're modelling)
//	          for (Dose in 1:NumDoses)  
//	         for (Group in 1:NumGroup) // i.e. NumVaxGroup < NumGroup (should be NumVaxVar = 1 and NumVar = 1,2,3, or 4 depending on which variants we're modelling)
//	       	            VaxEffect[Dose, Group] = VaxEffect_nc[Dose, 1] * phi; /// WRONG!!!
//
//	         
//	       }
//	     
//	   }
	
	
	if(DoVariants) {
	 VarAdvantage[2:NumVar] 	= VarAdvantage_nc	* phi;
	}
	
	// Initialise lambda if we are not doing knots: free parameters allowed
	lambda_raw 	= lambda_raw_nc 	* phi4;
	{  
			// initialize lambda matrix to have same values for every LTLA 
			int ind_par = 0; // initialize index
			for (i in 1:NumLTLAs){
		    NationalTrend[(ind_par + 1):(ind_par + NumTimepoints)] = lambda_raw[1:NumTimepoints];
				ind_par = ind_par + NumTimepoints; // update index
			}
	}
	
	for (i in 1:NumDatapoints)
	RegionalTrends[i] += NationalTrend[i] * RegionalScale[LTLAs[i]]; // lambda * RegionalScale^T in manuscript. Note use of intercept makes this line akin to IntDim (B) = 2 with column vector of 1s for one column of lambda

	
	// Danny changes - consolodating loops below
	{ // Bracket to compile
	int VaxVariantIndex;  // set to 1 by default
	int VaxAgeIndex    ;  // set to 1 by default
	
	real ProductOfDoses_PerVariantPerAgeGroup; 
	real WeightedSumEfficacyOverAgeGroups_PerVariant;   
	
	for (TimeRegion in 1:NumDatapoints)
	{
	  LogPredictions[TimeRegion] = 0; // initialize final Rt to zero for each time point and region.
	  
	  for (Variant in 1:NumVar) // sum over variants
	  {
	      // choose appropriate index for variant vaccine 
	      if (NumVaxVar == NumVar) VaxVariantIndex = Variant; else VaxVariantIndex = 1; // don't need the else as is specified by default above, but you get idea.
        
        // (re-)initialize to 0
        WeightedSumEfficacyOverAgeGroups_PerVariant = 0;   

	      FinalRtperVariantTimeRegion = VarProp[TimeRegion, Variant] * VarAdvantage[Variant] * RegionalTrends[TimeRegion]; 
	      
	      for (Group in 1:NumGroup)
	      {
  	      // choose appropriate index for age vaccine 
	        if (NumVaxGroup == NumGroup) VaxAgeIndex = Group; else VaxAgeIndex = 1; // don't need the else as is specified by default above, but you get idea.
	        
	        // calculate product (i.e. effect) of all doses (for this variant and age group)
	        ProductOfDoses_PerVariantPerAgeGroup = 1; // (re-)initialize.
  	      for (Dose in 1:NumDoses)
  	      {
  	        ProductOfDoses_PerVariantPerAgeGroup                   *= 
  	            (1 - (VaxProp[TimeRegion,Dose, Group]              *     // note has Group, not index
  	            VaxEffect[Dose, VaxVariantIndex, VaxAgeIndex]));         // note this has indices, not Group/Variant.
  	      }
  	      // add to age group sum
  	      WeightedSumEfficacyOverAgeGroups_PerVariant += AgeProp[Group] * ProductOfDoses_PerVariantPerAgeGroup;
	      }
	      
	      FinalRtperVariantTimeRegion *= WeightedSumEfficacyOverAgeGroups_PerVariant; 

	      LogPredictions[TimeRegion] += FinalRtperVariantTimeRegion; 
	  }
	}
	
	} // Bracket to compile
}

model {
	
	for (i in 1:NumLTLAs)
		RegionalScale_nc[i] ~ std_normal();
	
	for (i in 1:NumTrendPar)
		lambda_raw_nc[i] ~ std_normal(); 
	
	intercept_nc 	~ std_normal();
	
	phi_nc 			~ std_normal();
	phi2_nc 		~ std_normal();
	phi3_nc 		~ std_normal();
	phi4_nc 		~ std_normal();
	sigma_nc 		~ std_normal();
	
	for (i in 1:NumDoses)
	  for (j in 1:NumVaxVar)
	    for(k in 1:NumVaxGroup)
		    VaxEffect_nc[i, j, k] ~ std_normal();
		  
	for (i in 1:(NumVar-1))
	  VarAdvantage_nc[i] ~ std_normal();
		
	RtVals 			~ normal(LogPredictions, sigma);
}

generated quantities {
	vector[NumDatapoints] log_lik;
	
	for (i in 1:NumDatapoints) {
		log_lik[i] = normal_lpdf(RtVals[i] | LogPredictions[i], sigma); // NOTES: Log of normal function: real normal_lpdf(reals y | reals mu, reals sigma)
	}
}
