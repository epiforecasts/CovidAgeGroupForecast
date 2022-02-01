
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

parameters {
  // incidence of infection and ab prevalence modelled as parameters drawn from 
  // independent normal distributions per time step.
  real <lower=0> infections[T,A]; 
  real <lower=0> antibodies[T,A];
  
  //parameters to fit
  real<lower=0> inf_rate[A];       // age specific rate of infection conditional on contact with 100% susceptibility
  real<lower=0> susceptibility[A]; // age specific susceptibility 
  real<lower=0,upper=1> ab_protection; // protection offered by antibodies
  real<lower=0> sigma; // uncertainty in estimates
}


model {
  
  
  vector[A] full_susceptibility;                // container for ab and inherent susceptibiluty
  matrix[A,A] transmissibility_correction;      // container for conversion of CM to NGM
  vector[A] next_gen;

  inf_rate ~ gamma(rep_array(0.5, A), 0.2);
  susceptibility ~ gamma(rep_array(1.0, A), 0.2);
  ab_protection ~ beta(5.0, 1.0);
  
  sigma ~ gamma(0.05, 0.01);
  
  for(t in 1:T){
    for(a in 1:A){
      infections[t,a] ~ normal(inf_mu[t,a], inf_sd[t,a])T[0,];
      antibodies[t,a] ~ normal(anb_mu[t,a], anb_sd[t,a])T[0,]; 
      
    }
    
  }
  

  
  for(t in 4:T){
    
    full_susceptibility = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t,:])));
    transmissibility_correction = full_susceptibility * to_row_vector(inf_rate);
    
    next_gen = (contact_matrices[day_to_week_converter[t]] .* transmissibility_correction)  * to_vector(infections[t,:]);
    for(a in 1:A){
      infections[t-3, a] ~ normal(next_gen[a], sigma);
    }
    
  }
  
  
}
