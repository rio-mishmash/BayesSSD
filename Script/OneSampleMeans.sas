
%macro m_Bayes_OneSampleMeans(
    nullmean,
    mean,
    stddev,
    ntotal,
    interim,
    sides,
    lambda,

    margin     =    0,
    prior_eta  =  0.0, 
    prior_tau  =  100, 
    gamma_L    =    0,
    gamma_U    =    1,

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
	    	ntotal,
            prior_eta,
            prior_tau,
	    	y_mean,
            stddev,
            nullmean,
	    	lambda,
	    	sides,
            margin
		);

	        * posterior probability;
			if ntotal = 0 then do;
                eta = prior_eta;
                tau = prior_tau;
			end;
			else do;
                eta = ((prior_tau**2)      *    y_mean) / ((stddev**2)/ntotal + prior_tau**2) 
                    + ((stddev**2)/ntotal  * prior_eta) / ((stddev**2)/ntotal + prior_tau**2);
                tau =  ((1/(prior_tau**2) + ntotal/(stddev**2))**(-1))**(1/2);
            end;

	        if sides = "U" then do;
		        posterior_prob  = 1 - cdf("normal", nullmean + margin, eta, tau);
		    end;
		    if sides = "L" then do;
		        posterior_prob  =     cdf("normal", nullmean + margin, eta, tau);
		    end;

			return posterior_prob;
		finish;

		*****  predictive probability *****;
	
	    start f_predictive_prob(
	    	ntotal,
	    	ncurrent,
            prior_eta,
            prior_tau,
	    	y_mean_obs,
            stddev,
            nullmean,
	    	lambda,
	    	sides,
            margin
		);

	        * posterior probability;
			if ncurrent = 0 then do;
				eta = prior_eta;
                tau = prior_tau;
			end;
			else do;
                eta = ((prior_tau**2)        * y_mean_obs) / ((stddev**2)/ncurrent + prior_tau**2) 
                    + ((stddev**2)/ncurrent  *  prior_eta) / ((stddev**2)/ncurrent + prior_tau**2);
                tau =  ((1/(prior_tau**2) + ncurrent/(stddev**2))**(-1))**(1/2);
            end;
            
            if (ntotal-ncurrent) > 0 then do;
			
				* predictive response;			
                y_mean_pred = J(&sims., &nmc., .);
                do k = 1 to &sims. by 100; * simulate with size 100 to prevent memory error;
                    * number of rows;
                    nrow = min(k+100-1, &sims.) -k +1;
                    * extract rows from eta;
                    _eta = eta[min(k:(k+nrow-1), nrow(eta)),];
                    
                    * sampling from the posterior distribution (parameter);
                    epsilon = J(nrow, &nmc., .);
                    call randgen(epsilon, "normal", 0, tau); * mean=0;
                    * mean=0 -> mean=eta -> long; 
                    mu_sample_long = shape( epsilon + _eta , nrow*&nmc.);
                    * sampling from the posterior distribution (outcome);
                    epsilon = J(nrow*&nmc., ntotal-ncurrent, .);
                    call randgen(epsilon, "normal", 0, stddev); * mean=0;
                    * mean=0 -> mean=mu; 
                    mu_sample = mu_sample_long + epsilon;
                    * rowMeans -> wide; 
                    y_mean_pred[k:(k+nrow-1),] = shape( mu_sample[,:] , nrow);
                end;
                
                posterior_prob = f_posterior_prob( ntotal, prior_eta, prior_tau, 
                                                  (ncurrent*y_mean_obs
                                                   + (ntotal-ncurrent)*y_mean_pred) / ntotal, 
                                                   stddev, nullmean, lambda, sides, margin);
			end;
			else do;
                posterior_prob = f_posterior_prob( ntotal, prior_eta, prior_tau, 
                                                   y_mean_obs, 
                                                   stddev, nullmean, lambda, sides, margin);
			end;
			predictive_prob = (posterior_prob >= lambda)[,:]; *rowMeans;
			
			return predictive_prob; 
        finish;

		*****  execution *****;
				
		results = J(ncol(&ntotal.), 7, .);
		create Results from results[colname=({lambda gamma_L samplesize power ESS interim se})];
		do i = 1 to ncol(&ntotal.);
		    
            n      = J(     1, ncol(&interim.), .);               n[,1] = 0;
            y_mean = J(&sims., ncol(&interim.), .);          y_mean[,1] = 0;
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
                    n[,j] = floor( &ntotal.[i]*(&interim.[j]-&interim.[j-1]) );
                    %* response;
                    y     = J(&sims., n[,j], .);
                    call randgen(y, "normal", &mean., &stddev.);
                    y_mean[,j] = y[,:]; *RowMeans;
                end;

                posterior_prob[,j] = f_posterior_prob(n[,1:j][,+], &prior_eta., &prior_tau.,
                                                      (n[,1:j]#y_mean[,1:j])[,+] / max(n[,1:j][,+],1),  /*weighted sum*/
                                                      &stddev., &nullmean., &lambda., &sides., &margin.);
				
                predictive_prob[,j] = f_predictive_prob(&ntotal.[i],
                                                        n[,1:j][,+], &prior_eta., &prior_tau.,
                                                        (n[,1:j]#y_mean[,1:j])[,+] / max(n[,1:j][,+],1), 
                                                        &stddev., &nullmean., &lambda., &sides., &margin.);
				
				if j > 1 then do;
					* stopping with promising result;
                    promising[,j] = (posterior_prob[,j]  > &lambda. );
					* stopping with futility  result;
                    futility[,j]  = (predictive_prob[,j] < &gamma_L.);
				end;

			end;
			
			met = (promising[,+] > 0);
			/* print n y_mean posterior_prob[format=8.3] promising futility; */

			* samplesize & power;
			samplesize = &ntotal.[i];
			interim = ncol(&interim.)-2;
			power = mean(met);
			se    =  std(met) / sqrt(&sims.);
			ESS     = (n ## ( ongoing=1 ))[,+][:,];
			results = &lambda. || &gamma_L. || samplesize || power || ESS || interim || se;
			
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
