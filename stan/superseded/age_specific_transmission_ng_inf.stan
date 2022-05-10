
data {
  int<lower=0> T; //number of days in time-series
  int<lower=0> W; //number of weeks
  int<lower=0> A; //age groups
  
  real inf_mu[T,A];                    // infections parameter matrix
  real inf_sd[T,A];                    // infections parameter matrix
  real inf_mu_ng[T,A];
  real anb_mu[T,A];                    // antibodies parameter matrix
  real anb_sd[T,A];                    // antibodies parameter matrix
  matrix[A,A] contact_matrices [W];      // contact matricies
  int day_to_week_converter[T];      // id vector to determine which cm to use for each day
}

transformed data{
  matrix [T,A] infections_fixed;
  matrix [T,A] antibodies_fixed;
  infections_fixed = to_matrix(inf_mu);
  antibodies_fixed = to_matrix(anb_mu);
}

parameters {
  // incidence of infection and ab prevalence modelled as parameters drawn from 
  // independent normal distributions per time step.
  real <lower=0> infections[T,A]; 
  real <lower=0, upper=1> antibodies[T,A];
  
  //parameters to fit
  real <lower=0, upper=1> inf_rate_hyper_mu;
  real <lower=0> inf_rate_hyper_sd; 
  real <lower=0, upper=1> suscept_hyper_mu;
  real <lower=0> suscept_hyper_sd; 
  real inf_prime[A];
  real sus_prime[A];
  //real<lower=0,upper=1> inf_rate[A];       // age specific rate of infection conditional on contact with 100% susceptibility
  //real<lower=0,upper=1> susceptibility[A]; // age specific susceptibility 
  real<lower=0,upper=1> ab_protection; // protection offered by antibodies
 // real<lower=0> sigma; // uncertainty in estimates
}

transformed parameters{
  real inf_rate[A];
  real susceptibility[A];
  for(a in 1:A){
    inf_rate[a] = inf_rate_hyper_mu + inf_rate_hyper_sd * inf_prime[a];
    susceptibility[a] = suscept_hyper_mu + suscept_hyper_sd * sus_prime[a];
  }
  
  
}


model {
  

  vector[A] full_susceptibility;                // container for ab and inherent susceptibiluty
  vector[A] next_gen[T];


  //inf_rate ~ beta(5.0, 1.0);       // age specific rate of infection on contact
  //susceptibility ~ beta(5.0, 1.0);   // age specific inherent susceptibility 
  ab_protection ~ beta(5.0, 1.0);                  // protectiveness of antibodies
  
  inf_rate_hyper_mu ~ normal(0.25, 0.01)T[0,1];
  inf_rate_hyper_sd ~ normal(0.05,0.001)T[0,];
  suscept_hyper_mu ~ normal(0.5, 0.01)T[0,1];
  suscept_hyper_sd ~ normal(0.1, 0.001)T[0,];
  
  for(a in 1:A){
    inf_prime[a] ~ normal(0,0.1);
    sus_prime[a] ~ normal(0,0.1);
    }
  
  
  //sigma ~ gamma(0.05, 0.01);
  // sample priors for infections and antibodies 
  for(t in 1:T){
    for(a in 1:A){
      infections[t,a] ~ normal(inf_mu[t,a], inf_sd[t,a])T[0,];
      antibodies[t,a] ~ normal(anb_mu[t,a], anb_sd[t,a])T[0,1]; 
      
    }
    
  }
  

  

  

  
  for(t in 1:(T-5)){
    // combine inherent and ab susceptibility
    full_susceptibility = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies_fixed[t])));
    
    // estimate next generation of infections
    next_gen[t] = (diag_matrix(full_susceptibility) * contact_matrices[day_to_week_converter[t]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections_fixed[t]);
    
    // fit model 
    for(a in 1:A){
      //infections_fixed[t+3, a] ~ normal(next_gen[a], inf_sd[t,a]);
      next_gen[t][a] ~ normal(inf_mu_ng[t+5,a], inf_sd[t+5, a]);
    }
    
  }
  
  
}


generated quantities{
  vector[A] full_susceptibility[T];                // container for ab and inherent susceptibiluty
  vector[A] next_gens[T];
  
  
  for(t in 6:T){
    // combine inherent and ab susceptibility
    full_susceptibility[t-5] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies_fixed[t-5])));

    // estimate next generation of infections
    next_gens[t] = (diag_matrix(full_susceptibility[t-5]) * contact_matrices[day_to_week_converter[t-5]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections_fixed[t-3]);
    //next_gens[t] = (contact_matrices[day_to_week_converter[t-3]] .* transmissibility_correction[t-3])  * to_vector(infections_fixed[t-3,:]);
  }
  
}
