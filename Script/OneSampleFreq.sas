
%macro m_Bayes_OneSampleFreq(
    nullproportion,	/* Specifies the null proportion */
    proportion,		/* Specifies the binomial proportion */
    ntotal,			/* Specifies the sample size */
    sides,			/* Specifies the direction of the statistical test */
    lambda,			/* Specifies the promising threshold at the final analysis */
	
    interim  =  {0, 1},	/* Specifies the analysis time points (0=baseline, 1=final analysis)  */
    margin   =   0, /* Specifies the equivalence or noninferiority or superiority margin */
    prior_a  =   1, /* Specifies the first  parameter of the beta prior */
    prior_b  =   1, /* Specifies the second parameter of the beta prior */
    gamma_L  =   ., /* Specifies the futility threshold at interim analyses */
    gamma_U  =   ., /* Specifies the promising threshold at interim analyses */

    sims   =  1000, /* Specifies the number of simulations */
    nmc    =  1000, /* Specifies the number of Monte Carlo samples */
    seed   =  2025  /* Specifies the seed for random number generation */
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
			prior_a,
			prior_b,
			y_sum,
			nullproportion,
			lambda,
			sides,
			margin
		);
			* posterior probability;
			if ntotal = 0 then do;
				a = prior_a;
				b = prior_b;
			end;
			else do;
				a = prior_a + y_sum;
				b = prior_b + ntotal - y_sum;
			end;

			if sides = "U" then do;
				posterior_prob  = 1 - cdf("beta", nullproportion + margin, a, b);
			end;
			if sides = "L" then do;
				posterior_prob  =     cdf("beta", nullproportion + margin, a, b);
			end;
	        
			return posterior_prob;
		finish;
			
		*****  predictive probability *****;
	
		start f_predictive_prob(
			ntotal,
			ncurrent,
			prior_a,
			prior_b,
			y_sum_obs,
			nullproportion,
			lambda,
			sides,
			margin
		);
			* posterior distribution;
			if ncurrent = 0 then do;
				a = prior_a;
				b = prior_b;
			end;
			else do;
				a = prior_a + y_sum_obs;
				b = prior_b + ntotal - y_sum_obs;
			end;
			
			if (ntotal-ncurrent) > 0 then do;
			
				* predictive response;
                y_sum_pred = J(&sims., &nmc., .);
                do k = 1 to &sims. by 100; * simulate with size 100 to prevent memory error;
                    * number of rows;
                    nrow = min(k+100-1, &sims.) -k +1;
					
					* sampling from the posterior distribution;
					p_sample = J(nrow, &nmc., .);
					call randgen(p_sample, "beta",  a, b);
					y_sum_pred[k:(k+nrow-1),]  = rand("binomial", p_sample, ntotal-ncurrent);
                end;
				
				posterior_prob = f_posterior_prob( ntotal, prior_a, prior_b, 
										     	   y_sum_obs + y_sum_pred,
											       nullproportion, lambda, sides, margin);
			end;
			else do;
				posterior_prob = f_posterior_prob( ntotal, prior_a, prior_b, 
										     	   y_sum_obs,
											       nullproportion, lambda, sides, margin);
			end;
			predictive_prob = (posterior_prob >= lambda)[,:]; *rowMeans;
			
			return predictive_prob;
		finish;

		*****  execution *****;
		
		results = J(ncol(&ntotal.), 8, .);
		create Results from results[colname=({lambda gamma_U gamma_L ntotal power ESS interim se})];
		do i = 1 to ncol(&ntotal.);
		    
			%* init;
			n     = J(     1, ncol(&interim.), .);                n[,1] = 0;
			y_sum = J(&sims., ncol(&interim.), .);		      y_sum[,1] = 0;
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
					%* response probability;
					_p = J(&sims., 1, .);
					do k = 1 to &sims.;
						_p[k,] = &proportion.;
					end;
					%* response;
					y_sum[,j] = rand("binomial", _p, n[,j]);
				end;
					
				posterior_prob[,j] = f_posterior_prob(  n[,1:j][,+], &prior_a., &prior_b., 
														y_sum[,1:j][,+], &nullproportion., 
														&lambda., &sides., &margin.);

				predictive_prob[,j] = f_predictive_prob(&ntotal.[i], 
														n[,1:j][,+], &prior_a., &prior_b.,
														y_sum[,1:j][,+], &nullproportion., 
														&lambda., &sides., &margin.);

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
			
			met = promising[,+] > 0;
			/* print n y_sum posterior_prob[format=8.3] predictive_prob[format=8.3] promising; */
			
			* ntotal & power;
			ntotal  = &ntotal.[i];
			interim = ncol(&interim.)-2;
			power = mean(met);
			se    =  std(met) / sqrt(&sims.);
			ESS     = (n ## ( ongoing=1 ))[,+][:,];
			results = &lambda. || &gamma_U. || &gamma_L. || ntotal || power || ESS || interim || se;

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
