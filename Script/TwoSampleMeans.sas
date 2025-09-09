/* 
    refmean        Specifies the reference mean 
    meandiff       Specifies the mean difference 
    stddev         Specifies the common standard deviation 
    npergroup      Specifies the common sample size per group 
    sides          Specifies the direction of the statistical test 
    lambda         Specifies the promising threshold at the final analysis 
    interim        Specifies the analysis time points (0=baseline, 1=final analysis)
    margin         Specifies the equivalence or noninferiority or superiority margin 
    prior0_eta     Specifies the mean parameter of the normal prior for group 0  
    prior0_tau     Specifies the stddev parameter of the normal prior for group 0  
    prior1_eta     Specifies the mean parameter of the normal prior for group 1  
    prior1_tau     Specifies the stddev parameter of the normal prior for group 1  
    gamma_L        Specifies the futility threshold at interim analyses 
    gamma_U        Specifies the promising threshold at interim analyses 
    sims           Specifies the number of simulations 
    nmc            Specifies the number of Monte Carlo samples 
    seed           Specifies the seed for random number generation 
*/
%macro m_Bayes_TwoSampleMeans(
    refmean, 
    meandiff, 
    stddev, 
    npergroup, 
    sides, 
    lambda,
	
    interim  = %quote({0.0 1.0}),
    margin      =    0,
    prior0_eta  =  0.0,
    prior0_tau  =  100,
    prior1_eta  =  0.0,
    prior1_tau  =  100,
    gamma_L  =   .,
    gamma_U  =   .,

    sims   =  1000,
    nmc    =  1000,
    seed   =  2025
);
    %* time to start;
    data start;
        start = time();
    run;
    
    
	proc iml;
    	reset noprint nolog;
    	
        * set seed;
        call streaminit(&seed.);
        call   randseed(&seed.);

		*****  posterior probability *****;

	    start f_posterior_prob(
	    	npergroup,
            prior0_eta,
            prior0_tau,
            prior1_eta,
            prior1_tau,
	    	y0_mean,
	    	y1_mean,
            stddev,
	    	lambda,
	    	sides,
            margin
		);

	        * posterior probability;
			if npergroup = 0 then do;
                eta0 = prior0_eta;
                tau0 = prior0_tau;
                eta1 = prior1_eta;
                tau1 = prior1_tau;
			end;
			else do;
                eta0 =  ((prior0_tau**2)        *    y0_mean) / ((stddev**2)/npergroup + prior0_tau**2) 
                      + ((stddev**2)/npergroup  * prior0_eta) / ((stddev**2)/npergroup + prior0_tau**2);
                tau0 =  ((1/(prior0_tau**2) + npergroup/(stddev**2))**(-1))**(1/2);
                eta1 =  ((prior1_tau**2)        *    y1_mean) / ((stddev**2)/npergroup + prior1_tau**2) 
                      + ((stddev**2)/npergroup  * prior1_eta) / ((stddev**2)/npergroup + prior1_tau**2);
                tau1 =  ((1/(prior1_tau**2) + npergroup/(stddev**2))**(-1))**(1/2);
            end;

	        if sides = "U" then do;
		        posterior_prob  = 1 - cdf("normal", margin, eta1 - eta0, sqrt(tau1**2 + tau0**2));
		    end;
		    if sides = "L" then do;
		        posterior_prob  =     cdf("normal", margin, eta1 - eta0, sqrt(tau1**2 + tau0**2));
		    end;

			return posterior_prob;
		finish;

		*****  predictive probability *****;

	    start f_predictive_prob(
	    	npergroup,
	    	ncurrent,
            prior0_eta,
            prior0_tau,
            prior1_eta,
            prior1_tau,
	    	y0_mean_obs,
	    	y1_mean_obs,
            stddev,
	    	lambda,
	    	sides,
            margin
		);
            
	        * posterior probability;
			if ncurrent = 0 then do;
				eta0 = prior0_eta;
                tau0 = prior0_tau;
				eta1 = prior1_eta;
                tau1 = prior1_tau;
			end;
			else do;
                eta0 =  ((prior0_tau**2)        * y0_mean_obs) / ((stddev**2)/npergroup + prior0_tau**2) 
                      + ((stddev**2)/npergroup  *  prior0_eta) / ((stddev**2)/npergroup + prior0_tau**2);
                tau0 =  ((1/(prior0_tau**2) + npergroup/(stddev**2))**(-1))**(1/2);
                eta1 =  ((prior1_tau**2)        * y1_mean_obs) / ((stddev**2)/npergroup + prior1_tau**2) 
                      + ((stddev**2)/npergroup  *  prior1_eta) / ((stddev**2)/npergroup + prior1_tau**2);
                tau1 =  ((1/(prior1_tau**2) + npergroup/(stddev**2))**(-1))**(1/2);
            end;
            
            if (npergroup-ncurrent) > 0 then do;
			
				* predictive response;
                y0_mean_pred = J(&sims., &nmc., .);
                y1_mean_pred = J(&sims., &nmc., .);
                do k = 1 to &sims. by 100; * simulate with size 100 to prevent memory error;
                    * number of rows;
                    nrow = min(k+100-1, &sims.) -k +1;
                    * extract rows from eta;
                    _eta0 = eta0[min(k:(k+nrow-1), nrow(eta0)),];
                    _eta1 = eta1[min(k:(k+nrow-1), nrow(eta1)),];
                    
                    * arm 0;
                    * sampling from the posterior distribution (parameter);
                    epsilon = J(nrow, &nmc., .);
                    call randgen(epsilon, "normal", 0, tau0); * mean=0;
                    * mean=0 -> mean=eta -> long; 
                    mu0_sample_long = shape( epsilon + _eta0 , nrow*&nmc.);
                    * sampling from the posterior distribution (outcome);
                    epsilon = J(nrow*&nmc., npergroup-ncurrent, .);
                    call randgen(epsilon, "normal", 0, stddev); * mean=0;
                    * mean=0 -> mean=mu; 
                    mu0_sample = mu0_sample_long + epsilon;
                    * rowMeans -> wide; 
                    y0_mean_pred[k:(k+nrow-1),] = shape( mu0_sample[,:] , nrow);

                    * arm 1;
                    * sampling from the posterior distribution (parameter);
                    epsilon = J(nrow, &nmc., .);
                    call randgen(epsilon, "normal", 0, tau1); * mean=0;
                    * mean=0 -> mean=eta -> long; 
                    mu1_sample_long = shape( epsilon + _eta1 , nrow*&nmc.);
                    * sampling from the posterior distribution (outcome);
                    epsilon = J(nrow*&nmc., npergroup-ncurrent, .);
                    call randgen(epsilon, "normal", 0, stddev); * mean=0;
                    * mean=0 -> mean=mu; 
                    mu1_sample = mu1_sample_long + epsilon;
                    * rowMeans -> wide; 
                    y1_mean_pred[k:(k+nrow-1),] = shape( mu1_sample[,:] , nrow);
                end;
                
                posterior_prob = f_posterior_prob( npergroup, prior0_eta, prior0_tau, prior1_eta, prior1_tau, 
                                                  (ncurrent*y0_mean_obs + (npergroup-ncurrent)*y0_mean_pred) / npergroup, 
                                                  (ncurrent*y1_mean_obs + (npergroup-ncurrent)*y1_mean_pred) / npergroup, 
                                                   stddev, lambda, sides, margin);
			end;
            else do;
                posterior_prob = f_posterior_prob( npergroup, prior0_eta, prior0_tau, prior1_eta, prior1_tau, 
                                                   y0_mean_obs, y1_mean_obs, 
                                                   stddev, lambda, sides, margin);
            end;
            predictive_prob = (posterior_prob >= lambda)[,:]; *rowMeans;

            return predictive_prob;
        finish;

		*****  execution *****;
		
		results = J(ncol(&npergroup.), 8, .);
		create Results from results[colname=({lambda gamma_U gamma_L npergroup power ESS interim se})];
		do i = 1 to ncol(&npergroup.);
		
            n       = J(     1, ncol(&interim.), .);               n[,1] = 0;
            y0_mean = J(&sims., ncol(&interim.), .);         y0_mean[,1] = 0;
            y1_mean = J(&sims., ncol(&interim.), .);         y1_mean[,1] = 0;
            posterior_prob  = J(&sims., ncol(&interim.), .);
            predictive_prob = J(&sims., ncol(&interim.), .);
			promising       = J(&sims., ncol(&interim.), 0);
			futility        = J(&sims., ncol(&interim.), 0);
			ongoing         = J(&sims., ncol(&interim.), 1);
		    
		    do j = 1 to ncol(&interim.);
			
				* if already stopped then 0;
				ongoing[,j]   = (promising[,1:j][,+] + futility[,1:j][,+] = 0);
		    
				%* after enrollment;
				if j > 1 then do;
                    %* number of subjects;
                    n[,j] = floor( &npergroup.[i]*(&interim.[j]-&interim.[j-1]) );
                    %* response;
                    y0     = J(&sims., n[,j], .);
                    call randgen(y0, "normal", &refmean.,            &stddev.);
                    y0_mean[,j] = y0[,:]; *RowMeans;
                    y1     = J(&sims., n[,j], .);
                    call randgen(y1, "normal", &refmean.+&meandiff., &stddev.);
                    y1_mean[,j] = y1[,:]; *RowMeans;
                end;

                posterior_prob[,j] = f_posterior_prob( n[,1:j][,+], &prior0_eta., &prior0_tau., &prior1_eta., &prior1_tau.,
                                                      (n[,1:j]#y0_mean[,1:j])[,+] / max(n[,1:j][,+], 1), 
                                                      (n[,1:j]#y1_mean[,1:j])[,+] / max(n[,1:j][,+], 1),  /*weighted sum*/
                                                      &stddev., &lambda., &sides., &margin.);
				
                predictive_prob[,j] = f_predictive_prob(&npergroup.[i],
                                                         n[,1:j][,+], &prior0_eta., &prior0_tau., &prior1_eta., &prior1_tau.,
                                                        (n[,1:j]#y0_mean[,1:j])[,+] / max(n[,1:j][,+], 1), 
                                                        (n[,1:j]#y1_mean[,1:j])[,+] / max(n[,1:j][,+], 1),  /*weighted sum*/
                                                        &stddev., &lambda., &sides., &margin.);
				
				if 1 < j & j < ncol(&interim.) then do;
					* stop early with promising result;
                    if &gamma_U. = . then promising[,j] = (posterior_prob[,j]  > &lambda. );
                    if &gamma_U.^= . then promising[,j] = (predictive_prob[,j] > &gamma_U.);
					* stop early with futility  result;
                    futility[,j]  = (predictive_prob[,j] < &gamma_L.);
				end;
				else if j = ncol(&interim.) then do;
					* termination with promising result;
                    promising[,j] = (posterior_prob[,j]  > &lambda. );
				end;

			end;
			
			met = (promising[,+] > 0);
			/*print n y0_mean y1_mean posterior_prob predictive_prob met;*/

			* samplesize & power;
			npergroup = &npergroup.[i];
			interim = ncol(&interim.)-2;
			power = mean(met);
			se    =  std(met) / sqrt(&sims.);
			ESS     = (n ## ( ongoing=1 ))[,+][:,];
			results = &lambda. || &gamma_U. || &gamma_L. || npergroup || power || ESS || interim || se;
			
			append from results;
        end;
    quit;
	
	proc print data = work.RESULTS; 
        format se  8.3;
        format ESS 8.1;
    run;

    %* time to stop;
    data stop;
        set start;
        stop = (time()-start)/60;
        put "runtime(min) : " stop;
    run;

%mend;
