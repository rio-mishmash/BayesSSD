 /* 
	refproportion    Specifies the reference proportion 
	proportiondiff   Specifies the proportion difference 
	npergroup        Specifies the common sample size per group 
	sides            Specifies the direction of the statistical test 
	lambda           Specifies the promising threshold at the final analysis 
	interim          Specifies the analysis time points (0=baseline, 1=final analysis)
	margin           Specifies the equivalence or noninferiority or superiority margin 
	prior0_a         Specifies the first parameter of the beta prior for the group 0 
	prior0_b         Specifies the second parameter of the beta prior for the group 0 
	prior1_a         Specifies the first parameter of the beta prior for the group 1 
	prior1_b         Specifies the second parameter of the beta prior for the group 1 
	gamma_L          Specifies the futility threshold at interim analyses 
	gamma_U          Specifies the promising threshold at interim analyses 
	sims             Specifies the number of simulations 
	nmc              Specifies the number of Monte Carlo samples 
	seed             Specifies the seed for random number generation 
*/
%macro m_Bayes_TwoSampleFreq(
    refproportion, 
	proportiondiff, 
	npergroup, 
	sides, 
	lambda,
    
    interim  = %quote({0.0 1.0}),
    margin    =    0,
    prior0_a  =    1,
    prior0_b  =    1,
    prior1_a  =    1,
    prior1_b  =    1,
    gamma_L   =    .,
    gamma_U   =    .,

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
			prior0_a,
			prior0_b,
			prior1_a,
			prior1_b,
			y0_sum,
			y1_sum,
			lambda,
			sides,
			margin
		);

			* posterior probability;
			if npergroup = 0 then do;
				a0 = prior0_a;
				b0 = prior0_b;
				a1 = prior1_a;
				b1 = prior1_b;
			end;
			else do;
				a0 = prior0_a + y0_sum;
				b0 = prior0_b + npergroup - y0_sum;
				a1 = prior1_a + y1_sum;
				b1 = prior1_b + npergroup - y1_sum;
			end;
			p0 = a0 / (a0 + b0); * posterior mean;
			p1 = a1 / (a1 + b1); * posterior mean;

			* https://doi.org/10.1002/pst.1571;
			if sides = "U" then do;
				posterior_prob  = 1 - cdf("normal", margin, p1 - p0, 
											sqrt((p1#(1-p1))/(a1+b1+1) +(p0#(1-p0))/(a0+b0+1)));
			end;
			if sides = "L" then do;
				posterior_prob  =     cdf("normal", margin, p1 - p0, 
											sqrt((p1#(1-p1))/(a1+b1+1) +(p0#(1-p0))/(a0+b0+1)));
			end;
	        
			return posterior_prob;
		finish;

		*****  predictive probability *****;
	
		start f_predictive_prob(
			npergroup,
			ncurrent,
			prior0_a,
			prior0_b,
			prior1_a,
			prior1_b,
			y0_sum_obs,
			y1_sum_obs,
			lambda,
			sides,
			margin
		);

			* posterior distribution;
			if ncurrent = 0 then do;
				a0 = prior0_a;
				b0 = prior0_b;
				a1 = prior1_a;
				b1 = prior1_b;
			end;
			else do;
				a0 = prior0_a + y0_sum_obs;
				b0 = prior0_b + ncurrent - y0_sum_obs;
				a1 = prior1_a + y1_sum_obs;
				b1 = prior1_b + ncurrent - y1_sum_obs;
			end;

			if (npergroup-ncurrent) > 0 then do;			
			
				* predictive response;
                y0_sum_pred = J(&sims., &nmc., .);
                y1_sum_pred = J(&sims., &nmc., .);
                do k = 1 to &sims. by 100; * simulate with size 100 to prevent memory error;
                    * number of rows;
                    nrow = min(k+100-1, &sims.) -k +1;

					* arm 0;
					* sampling from the posterior distribution;
					p0_sample = J(nrow, &nmc., .);
					call randgen(p0_sample, "beta",  a0, b0);
					y0_sum_pred[k:(k+nrow-1),]  = rand("binomial", p0_sample, npergroup-ncurrent);

					* arm 1;
					* sampling from the posterior distribution;
					p1_sample = J(nrow, &nmc., .);
					call randgen(p1_sample, "beta",  a1, b1);
					y1_sum_pred[k:(k+nrow-1),]  = rand("binomial", p1_sample, npergroup-ncurrent);
                end;

				posterior_prob = f_posterior_prob(npergroup, 
												prior0_a, prior0_b, prior1_a, prior1_b, 
												y0_sum_obs + y0_sum_pred, y1_sum_obs + y1_sum_pred,
												lambda, sides, margin);				
			end;
			else do;
				posterior_prob = f_posterior_prob(npergroup, 
												prior0_a, prior0_b, prior1_a, prior1_b, 
												y0_sum_obs, y1_sum_obs,
												lambda, sides, margin);				
			end;
			predictive_prob = (posterior_prob >= lambda)[,:]; *rowMeans;
			
			return predictive_prob;
		finish;

		*****  execution *****;
		
		results = J(ncol(&npergroup.), 8, .);
		create Results from results[colname=({lambda gamma_U gamma_L npergroup power ESS interim se})];
		do i = 1 to ncol(&npergroup.);
		    
			n      = J(     1, ncol(&interim.), .);               n[,1] = 0;
			y0_sum = J(&sims., ncol(&interim.), .);          y0_sum[,1] = 0;
			y1_sum = J(&sims., ncol(&interim.), .);          y1_sum[,1] = 0;
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
					%* response probability;
					_p0 = J(&sims., 1, .);
					_p1 = J(&sims., 1, .);
					do k = 1 to &sims.;
						_p0[k,] = &refproportion.;
						_p1[k,] = _p0[k,]+&proportiondiff.;
					end;
					%* response;
					y0_sum[,j] = rand("binomial", _p0, n[,j]);
					y1_sum[,j] = rand("binomial", _p1, n[,j]);
				end;
									
				posterior_prob[,j] = f_posterior_prob(n[,1:j][,+], 
													&prior0_a., &prior0_b., &prior1_a., &prior1_b.,
													y0_sum[,1:j][,+], y1_sum[,1:j][,+],
													&lambda., &sides., &margin.);

				predictive_prob[,j] = f_predictive_prob(&npergroup.[i], n[,1:j][,+], 
														&prior0_a., &prior0_b., &prior1_a., &prior1_b.,
														y0_sum[,1:j][,+], y1_sum[,1:j][,+],
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
			/* print n y0_sum y1_sum posterior_prob predictive_prob met; */
			
			* npergroup & power;
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
