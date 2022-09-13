data {
  
  int<lower=0> T;                           //number of days in time-series
  int<lower=0> W;                           //number of weeks
  int<lower=0> A;                           //age groups

  real inf_mu[T,A];                         // infections parameter matrix
  real inf_sd[T,A];                         // infections parameter matrix
  real anb_mu[T,A];                         // antibodies parameter matrix
  real anb_sd[T,A];                         // antibodies parameter matrix
  
  matrix[A,A] contact_matrices_mu [W];      // contact matricies means
  matrix[A,A] contact_matrices_sd [W];      // contact matricies means
  real mean_contacts_mu[W,A];
  real mean_contacts_sd[W,A];
  real mean_contacts_mu_mat[W];
  real mean_contacts_sd_mat[W];
  
  int day_to_week_converter[T];             // id vector to determine which cm to use for each day
  
  real population[A];
  
  int<lower=0> smax;                        // maximum generation interval
  int<lower=0> horizon;                     // number of days to forecast
  int<lower=0> contact_option;              // 1-full matrices, 2-age-group means, 3-overall means, 4-no contact data.
  int<lower=0> sigma_option;                // 1-fit sd, 2-fit (pseudo) COV [COV normalised by the mean input infections as opposed to the fit value]

}


parameters {
  // incidence of infection and ab prevalence modelled as parameters drawn from 
     // independent normal distributions per time step.
  real <lower=0> infections[T,A]; 
  real <lower=0, upper=1> antibodies[T,A];
  
  //parameters to fit
  real <lower=0, upper=1> inf_rate_hyper_mu;  // hyper prior mean for inf rate
  real <lower=0> inf_rate_hyper_sd;           // hyper prior sd for inf rate
  real <lower=0, upper=1> suscept_hyper_mu;   // hyper prior mean for susceptibility
  real <lower=0> suscept_hyper_sd;            // hyper prior sd for susceptibility
  matrix<lower=0>[A,A] contact_matrices[W];   // contact matrices - Q pof
  real inf_prime[A];
  real sus_prime[A];
  //real<lower=0,upper=1> inf_rate[A];        // age specific rate of infection conditional on contact with 100% susceptibility
  //real<lower=0,upper=1> susceptibility[A];  // age specific susceptibility 
  real<lower=0,upper=1> ab_protection;        // protection offered by antibodies

  real<lower=0> w_mu; 
  real<lower=0> w_sig;
  
  real<lower=0> sigma_inf;                    // uncertainty in model
  
  matrix<lower=0>[A,A] sigma_cm;
  real<lower=0> sigma_mca[A];
  real<lower=0> sigma_mc;

}

transformed parameters{
  real inf_rate[A];                          // vector of relative infectiousness by age
  real susceptibility[A];                    // vector of relative susceptibility by age
  matrix[A,A] contact_matrices_aug[W];       // population corrected contact matrix parameter
  real w_g[smax];
  
  
  real<lower=0> combined_sigma[T,A];
  real<lower=0> combined_sigma_cm[contact_option == 1 ? W : 0, contact_option == 1 ? A : 0,  contact_option == 1 ? A : 0];
  real<lower=0> combined_sigma_mca[contact_option == 2 ? W : 0, contact_option == 2 ? A : 0];
  real<lower=0> combined_sigma_mc[contact_option == 3 ? W : 0];
  
  // inf and susc vectors are calculated from the hyper parameters and unique offset parameter under NCP framework
  for (a in 1:A) {
    inf_rate[a] = exp(inf_rate_hyper_mu + inf_rate_hyper_sd * inf_prime[a]);
    susceptibility[a] = exp(suscept_hyper_mu + suscept_hyper_sd * sus_prime[a]);
  }
  
  // calculate contact matrix weighted by the population distribution to maintain reciprocity
  for (w in 1:W) {
    for (ai in 1:A) {
      for (aj in 1:A) {
        if (ai <= aj) contact_matrices_aug[w, ai, aj] = contact_matrices[w, ai, aj] * population[ai];
        else contact_matrices_aug[w, ai, aj] = contact_matrices[w, aj, ai] * population[ai];
      }
    }
  }
    
  if (contact_option == 5) {
    for (w in 1:W) {
      for (ai in 1:A) {
        for (aj in 1:A) {
          if (ai <= aj) contact_matrices_aug[w, ai, aj] = (contact_matrices[w, ai, aj] * population[ai]) .* diag_matrix(rep_vector(1,A))[ai, aj];
          else contact_matrices_aug[w, ai, aj] = (contact_matrices[w, aj, ai] * population[ai]) .*  diag_matrix(rep_vector(1,A))[ai, aj];
        }
      }
    }
  }
  
  for (s in 1:smax) {
    w_g[s] = lognormal_cdf(s+1, w_mu, w_sig) - lognormal_cdf(s, w_mu, w_sig);
  }
  for (s in 1:(smax-1)) {
    w_g[s] = w_g[s]/sum(w_g);
  }
  
  
  if (sigma_option == 1) {
    for (t in 1:T) {
      for (a in 1:A) {
        combined_sigma[t,a] = sqrt(sigma_inf^2 + inf_sd[t,a]^2);
      }
    }
  }
    
  if (sigma_option == 2) {
    for (t in 1:T) {
      for (a in 1:A) {
        combined_sigma[t,a] = sqrt((sigma_inf * inf_mu[t,a])^2 + inf_sd[t,a]^2);
      }
    }
  }
  
  for (w in 1:W) {
    if (contact_option == 1) {
      for (a in 1:A) {
        for (b in 1:A) {
          combined_sigma_cm[w,a,b] = sqrt(sigma_cm[a,b]^2 + contact_matrices_sd[w,a,b]^2);
        }
      }
    }
    
    if (contact_option == 2) {
      for (a in 1:A) {
        combined_sigma_mca[w,a] = sqrt(sigma_mca[a]^2 + mean_contacts_sd[w,a]^2);
      }
    }
    
    if (contact_option == 3) {
      combined_sigma_mc[w] = sqrt(sigma_mc^2 + mean_contacts_sd_mat[w]^2);
    }
  }

}
    
model {
  
  vector[A] full_susceptibility[T];                // container for antibody and inherent susceptibiluty combined 
  vector[A] next_gens[T];                          // container for timeseries of infections generated for fit
  matrix[smax,A] next_gens_smax[T];                // container for timeseries of infections generated for fit as caused by infections originating at each date in the past
  row_vector[A] mean_conts[W];                     // mean contacts per age group
  real mat_mean[W];                                // mean contacts overall
  
  row_vector[A] ones = rep_row_vector(1,A);        // a row vector of ones 'A' long 
  
  ab_protection ~ beta(5.0, 1.0);                  // protectiveness of antibodies
  inf_rate_hyper_mu ~ normal(0.25, 0.05)T[0,1];    // priors for inf and susc hyper parameters 
  inf_rate_hyper_sd ~ normal(0.05,0.01)T[0,];
  suscept_hyper_mu ~ normal(0.5, 0.1)T[0,1];
  suscept_hyper_sd ~ normal(0.1, 0.02)T[0,];
  
  w_mu ~ normal(5.0/7.0, 1.0/7.0)T[0,];
  w_sig ~ normal(1.7/7.0, 0.17/7.0)T[0,];
  
  sigma_inf ~ normal(0.005, 0.0025) T[0,];
  
  for (ai in 1:A) {
    for (aj in 1:A) {
      sigma_cm[ai,aj] ~ normal(0.005, 0.0025) T[0,];
    }
  }
  
  for (a in 1:A) {
    sigma_mca[a] ~ normal(0.005, 0.0025) T[0,]; 
  }
  
  sigma_mc ~ normal(0.005, 0.0025) T[0,];

  if (contact_option == 1) {
    // prior for symetric contact matrix using full data
    for (w in 1:W) {
      for (ai in 1:A) {
        for (aj in 1:A) {
          contact_matrices[w, ai, aj] ~ gamma(2,2);
          contact_matrices_mu[w, ai, aj] ~ normal(contact_matrices[w, ai, aj], combined_sigma_cm[w,ai,aj])T[0,];
        }
      }
    }
  } else if (contact_option == 2) {
    // prior for symetric contact matrix using only age band means
    for (w in 1:W) {
      for (ai in 1:A) {
        for (aj in 1:A) {
          contact_matrices[w, ai, aj] ~ gamma(2,2);
        }
      }
      mean_conts[w] = ones * contact_matrices_aug[w];
      mean_contacts_mu[w] ~ normal(to_array_1d(mean_conts[w]), combined_sigma_mca[w]);
    }
  } else if (contact_option == 3) {
    // prior for symetric contact matrix using overall means
    for (w in 1:W) {
      for (ai in 1:A) {
        for (aj in 1:A) {
          contact_matrices[w, ai, aj] ~ gamma(2,2);
        }
      }
              
      mean_conts[w] = ones * contact_matrices_aug[w];
      mat_mean[w] = mean(to_array_1d(mean_conts[w]));
      mean_contacts_mu_mat[w] ~ normal(mat_mean[w], combined_sigma_mc[w]);
    }
  } else if (contact_option == 4) {
    // prior for symetric contact matrix using no contact data
    for (w in 1:W) {
      for (ai in 1:A) {
        for (aj in 1:A) {
          contact_matrices[w, ai, aj] ~ gamma(2,2);
        }
      }
    }
  } else if (contact_option == 5) {
    // prior for symetric contact matrix using no contact data
    for (w in 1:W) {
      for (ai in 1:A) {
        for (aj in 1:A) {
          contact_matrices[w, ai, aj] ~ gamma(2,2);
        }
      }
    }
  }
  
  // prior for age-specific susc and inf offset in NCP framework
  for (a in 1:A) {
    inf_prime[a] ~ normal(0,1);
    sus_prime[a] ~ normal(0,1);
  }
  
  // sample priors for infections and antibodies
  for (t in 1:T) {
    for (a in 1:A) {
      infections[t,a] ~ normal(inf_mu[t,a], inf_sd[t,a])T[0,];
      antibodies[t,a] ~ normal(anb_mu[t,a], anb_sd[t,a])T[0,1]; 
    }
  }

  for (t in 1:T) {
    // combine inherent and ab susceptibility
    full_susceptibility[t] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t])*ab_protection ));
    
    // construct matrix to correct contact matrix to NGM
    if (t>smax) {
      // calculate contribution for infections at t from dates between t-smax and t based on weights w_g
      for (s in 1:smax) {
        next_gens_smax[t][s] = to_row_vector(w_g[s] * (diag_matrix(full_susceptibility[t]) * contact_matrices_aug[day_to_week_converter[t]] * diag_matrix(to_vector(inf_rate)))  * to_vector(inf_mu[t-s]));
      }
    
      // sum over all dates
      next_gens[t] = to_vector(rep_row_vector(1,smax) *  to_matrix(next_gens_smax[t])) ;
    
      // sampling statement for calculated infections vs modelled infections
      for (a in 1:A) {
        inf_mu[t,a] ~ normal(next_gens[t][a], combined_sigma[t, a]);
      }
    }
  }

}
  
generated quantities {

  vector[A] full_susceptibility[T];                // container for ab and inherent susceptibiluty
  vector[A] next_gens[T];
  matrix[smax,A] next_gens_smax[T];
  vector[A] forecast_gens[T + horizon]; 
  matrix[smax,A] forecast_gens_smax[T + horizon]; 

  for (t in 1:T) {
    for (a in 1:A) {
      forecast_gens[t][a] = infections[t][a];
    }
  }
  
  for (t in 1:T) {
    // combine inherent and ab susceptibility
    full_susceptibility[t] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[t])*ab_protection ));
    
    // construct matrix to correct contact matrix to NGM
    if (t>smax) {
      for (s in 1:smax) {
        next_gens_smax[t][s] = to_row_vector(w_g[s] * (diag_matrix(full_susceptibility[t]) * contact_matrices_aug[day_to_week_converter[t]] * diag_matrix(to_vector(inf_rate)))  * to_vector(inf_mu[t-s]));
      }
      next_gens[t] = to_vector(rep_row_vector(1,smax) *  to_matrix(next_gens_smax[t])) ;
    }
  }
  
  // forecast cases
  for (f in 1:horizon) {
    //full_susceptibility[T+f] = to_vector(susceptibility) .* (1.0 - (to_vector(antibodies[T])*ab_protection ));
    for (s in 1:smax) {
      forecast_gens_smax[T+f][s] = to_row_vector(w_g[s] * diag_matrix(full_susceptibility[T]) * contact_matrices_aug[day_to_week_converter[T]] * diag_matrix(to_vector(inf_rate))  * to_vector(forecast_gens[T+f-s]));
    }
    forecast_gens[T+f] = to_vector(normal_rng(rep_row_vector(1,smax) *  to_matrix(forecast_gens_smax[T+f]), combined_sigma[T]));
  }

}
