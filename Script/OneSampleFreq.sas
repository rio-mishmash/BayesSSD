
proc datasets lib=work memtype=data kill nolist; quit;

%macro m_Bayes_OneSampleFreq(
    nullproportion,
    proportion,
    ntotal,
    interim,
    sides,
    lambda,

    margin   =   0,
    prior_a  =   1, 
    prior_b  =   1, 
    gamma_L  =-999,
    gamma_U  =+999,

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
			
			* predictive response;
			if (ntotal-ncurrent) > 0 then do;
			
				* sampling from the posterior distribution;
				p_sample = J(&sims., &nmc., .);
				call randgen(p_sample, "beta",  a, b);

				y_sum_pred  = rand("binomial", p_sample, ntotal-ncurrent);
				
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
		create Results from results[colname=({lambda gamma_L gamma_U ntotal power ESS interim se})];
		do i = 1 to ncol(&ntotal.);
		    
			%* init;
			n     = J(     1, ncol(&interim.), .);              n[,1] = 0;
			y_sum = J(&sims., ncol(&interim.), .);		      y_sum[,1] = 0;
			posterior_prob  = J(&sims., ncol(&interim.), .);
			predictive_prob = J(&sims., ncol(&interim.), .);
			promising       = J(&sims., ncol(&interim.), .);
			futility        = J(&sims., ncol(&interim.), .);
			
			do j = 1 to ncol(&interim.);
		    
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

				if j >= 2 then do;
					* early stopping with promising result;
                    promising[,j] = ((posterior_prob[,j] >= &lambda.) | (predictive_prob[,j] >= &gamma_U.)) + promising[,1:j-1][,+];
					* early stopping with futility result;
                    futility[,j]  = (                                   (predictive_prob[,j] <= &gamma_L.)) +  futility[,1:j-1][,+];
				end;

			end;
			
			met = promising[,+] > 0;
			/* print n y_sum posterior_prob predictive_prob met; */
			
			* ntotal & power;
			ntotal  = &ntotal.[i];
			interim = ncol(&interim.)-2;
			power = mean(met);
			se    =  std(met) / sqrt(&sims.);
			ESS     = (n ## ( promising + futility <= 1))[,+][:,];
			results = &lambda. || &gamma_L. || &gamma_U. || ntotal || power || ESS || interim || se;

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
