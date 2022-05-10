
data {
  int<lower=0> T; //number of days in time-series
  int<lower=0> W; //number of weeks
  int<lower=0> A; //age groups

  real inf_mu[T,A];                    // infections parameter matrix
  real inf_sd[T,A];                    // infections parameter matrix
  real anb_mu[T,A];                    // antibodies parameter matrix
  real anb_sd[T,A];                    // antibodies parameter matrix
  matrix[A,A] contact_matrices_mu [W];      // contact matricies means
  matrix[A,A] contact_matrices_sd [W];      // contact matricies means
  real population[A];
  
  int day_to_week_converter[T];      // id vector to determine which cm to use for each day
  
  int<lower=0> smax;
  real<lower=0> w_g[smax]; 
  int<lower=0> horizon;
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
  matrix<lower=0>[A,A] contact_matrices[W];
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
  matrix[A,A] contact_matrices_sym[W];
  matrix[A,A] contact_matrices_aug[W];
  
  for(a in 1:A){
    inf_rate[a] = inf_rate_hyper_mu + inf_rate_hyper_sd * inf_prime[a];
    susceptibility[a] = suscept_hyper_mu + suscept_hyper_sd * sus_prime[a];
  }
  
  
  
  
  for( w in 1:W){
      for(ai in 1:A){
        for(aj in 1:A){
          if (ai <= aj) contact_matrices_aug[w, ai, aj] = contact_matrices[w, ai, aj] * population[ai];
          else contact_matrices_aug[w, ai, aj] = contact_matrices[w, aj, ai] * population[ai];
    }}}
  
  
}
    




model {
  

  vector[A] full_susceptibility;                // container for ab and inherent susceptibiluty
  vector[A] next_gen[T];
  matrix[smax,A] next_gen_smax[T]; 


  //inf_rate ~ beta(5.0, 1.0);       // age specific rate of infection on contact
  //susceptibility ~ beta(5.0, 1.0);   // age specific inherent susceptibility 
  ab_protection ~ beta(5.0, 1.0);                  // protectiveness of antibodies
  
  inf_rate_hyper_mu ~ normal(0.25, 0.01)T[0,1];
  inf_rate_hyper_sd ~ normal(0.05,0.001)T[0,];
  suscept_hyper_mu ~ normal(0.5, 0.01)T[0,1];
  suscept_hyper_sd ~ normal(0.1, 0.001)T[0,];
  
  for(w in 1:W){
    for(ai in 1:A){
      for(aj in 1:A){
    
    contact_matrices[w, ai, aj] ~ normal(contact_matrices_mu[w, ai, aj], contact_matrices_sd[w, ai, aj])T[0,];
  }}}
  
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
  

  

  

  
  for(t in 21:T){
    // combine inherent and ab susceptibility
    full_susceptibility = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t])*ab_protection));
    

    
    // estimate next generation of infections based based on previous 20 days with each day weighted by the distribution of generation interval w_g
    for(s in 1:smax){
          next_gen_smax[t][s] = to_row_vector(w_g[s] * (diag_matrix(full_susceptibility) * contact_matrices_aug[day_to_week_converter[t-smax + s]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections[t-smax + s]));
    }
    
    next_gen[t] = to_vector(rep_row_vector(1,smax) *  to_matrix(next_gen_smax[t])) ;
    // next_gen[t] = (diag_matrix(full_susceptibility) * contact_matrices_aug[day_to_week_converter[t]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections[t]);
    
    // fit model 
    for(a in 1:A){
      //infections_fixed[t+3, a] ~ normal(next_gen[a], inf_sd[t,a]);
      next_gen[t][a] ~ normal(inf_mu[t,a], inf_sd[t, a]);
    }
    
  }
  
  
}


generated quantities{
  vector[A] full_susceptibility[T];                // container for ab and inherent susceptibiluty
  vector[A] next_gens[T];
  matrix[smax,A] next_gens_smax[T];
  vector[A] forecast_gens[T + horizon]; 
  matrix[smax,A] forecast_gens_smax[T + horizon]; 

  
  
  
  for (t in 1:T){
    for(a in 1:A){
      forecast_gens[t][a] = infections[t][a];
    }
  }
  
  
  
  for(t in 1:T){
    // combine inherent and ab susceptibility
    full_susceptibility[t] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t])*ab_protection ));
    
    // construct matrix to correct contact matrix to NGM
    if(t>20){
    
    for(s in 1:smax){
        next_gens_smax[t][s] = to_row_vector(w_g[s] * (diag_matrix(full_susceptibility[t-smax + s]) * contact_matrices_aug[day_to_week_converter[t-smax + s]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections[t-smax + s]));
    }
    
    next_gens[t] = to_vector(rep_row_vector(1,smax) *  to_matrix(next_gens_smax[t])) ;
      
}
    // estimate next generation of infections
    //next_gens[t] = (diag_matrix(full_susceptibility[t-5]) * contact_matrices_aug[day_to_week_converter[t-5]] * diag_matrix(to_vector(inf_rate)))  * to_vector(infections[t-3]);
    //next_gens[t] = (contact_matrices[day_to_week_converter[t-3]] .* transmissibility_correction[t-3])  * to_vector(infections_fixed[t-3,:]);
  }
  
  // forecast cases
  for(f in 1:horizon){
    
      //full_susceptibility[T+f] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[T])*ab_protection ));
      for(s in 1:smax){
        forecast_gens_smax[T+f][s] = to_row_vector(w_g[s] * diag_matrix(full_susceptibility[min(T+f-smax + s, T)]) * contact_matrices_aug[day_to_week_converter[min(T+f-smax + s, T)]] * diag_matrix(to_vector(inf_rate))  * to_vector(forecast_gens[min(T+f-smax + s, T)]));
      }
     forecast_gens[T+f] = to_vector(rep_row_vector(1,smax) *  to_matrix(forecast_gens_smax[T+f])) ;
   
  }

  
}
