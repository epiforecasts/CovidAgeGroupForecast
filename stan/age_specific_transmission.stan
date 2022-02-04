
data {
  int<lower=0> T; //number of days in time-series
  int<lower=0> W; //number of weeks
  int<lower=0> A; //age groups
  
  real inf_mu[T,A];                    // infections parameter matrix
  real inf_sd[T,A];                    // infections parameter matrix
  real anb_mu[T,A];                    // antibodies parameter matrix
  real anb_sd[T,A];                    // antibodies parameter matrix
  matrix[A,A] contact_matrices [W];      // contact matricies
  int day_to_week_converter[T];      // id vector to determine which cm to use for each day
}

transformed data{
  matrix [T,A] infections_fixed;
  infections_fixed = to_matrix(inf_mu);
}

parameters {
  // incidence of infection and ab prevalence modelled as parameters drawn from 
  // independent normal distributions per time step.
  real <lower=0> infections[T,A]; 
  real <lower=0, upper=1> antibodies[T,A];
  
  //parameters to fit
  real<lower=0,upper=1> inf_rate[A];       // age specific rate of infection conditional on contact with 100% susceptibility
  real<lower=0,upper=1> susceptibility[A]; // age specific susceptibility 
  real<lower=0,upper=1> ab_protection; // protection offered by antibodies
  real<lower=0> sigma; // uncertainty in estimates
}




model {
  

  vector[A] full_susceptibility;                // container for ab and inherent susceptibiluty
  matrix[A,A] transmissibility_correction;      // container for conversion of CM to NGM
  vector[A] next_gen;


  inf_rate ~ beta(5.0, 1.0);       // age specific rate of infection on contact
  susceptibility ~ beta(5.0, 1.0);   // age specific inherent susceptibility 
  ab_protection ~ beta(5.0, 1.0);                  // protectiveness of antibodies
  
  sigma ~ gamma(0.05, 0.01);
  // sample priors for infections and antibodies 
  for(t in 1:T){
    for(a in 1:A){
      infections[t,a] ~ normal(inf_mu[t,a], inf_sd[t,a])T[0,];
      antibodies[t,a] ~ normal(anb_mu[t,a], anb_sd[t,a])T[0,1]; 
      
    }
    
  }
  

  

  
  for(t in 1:(T-3)){
    // combine inherent and ab susceptibility
    full_susceptibility = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t,:]) * ab_protection));
    
    // construct matrix to correct contact matrix to NGM
    transmissibility_correction = full_susceptibility * to_row_vector(inf_rate);
    
    // estimate next generation of infections
    next_gen = (diag_matrix(full_susceptibility) * contact_matrices[day_to_week_converter[t]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections_fixed[t,:]);
    
    // fit model 
    for(a in 1:A){
      //infections_fixed[t+3, a] ~ normal(next_gen[a], inf_sd[t,a]);
      next_gen[a] ~ normal(inf_mu[t+3,a], inf_sd[t+3, a]);
    }
    
  }
  
  
}


generated quantities{
  vector[A] full_susceptibility[T];                // container for ab and inherent susceptibiluty
  matrix[A,A] transmissibility_correction[T];      // container for conversion of CM to NGM
  vector[A] next_gens[T];
  
  
  for(t in 4:T){
    // combine inherent and ab susceptibility
    full_susceptibility[t-3] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t-3,:]) * ab_protection));
    
    // construct matrix to correct contact matrix to NGM
    transmissibility_correction[t-3] = full_susceptibility[t-3] * rep_row_vector(1,7);//to_row_vector(inf_rate);
    
    // estimate next generation of infections
    next_gens[t] = (diag_matrix(full_susceptibility[t-3]) * contact_matrices[day_to_week_converter[t-3]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections_fixed[t-3,:]);
    //next_gens[t] = (contact_matrices[day_to_week_converter[t-3]] .* transmissibility_correction[t-3])  * to_vector(infections_fixed[t-3,:]);
  }
  
}
