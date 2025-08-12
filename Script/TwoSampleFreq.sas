
proc datasets lib=work memtype=data kill nolist; quit;

%macro m_Bayes_TwoSampleFreq(
    refproportion,
    proportiondiff,
    npergroup,
    interim,
    sides,
    lambda,

    margin    =    0,
    prior0_a  =    1, 
    prior0_b  =    1, 
    prior1_a  =    1, 
    prior1_b  =    1, 
    gamma_L   = -999,
    gamma_U   = +999,

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

			* predictive response;
			if (npergroup-ncurrent) > 0 then do;
			
				* sampling from the posterior distribution;
				p0_sample = J(&sims., &nmc., .);
				call randgen(p0_sample, "beta",  a0, b0);
				p1_sample = J(&sims., &nmc., .);
				call randgen(p1_sample, "beta",  a1, b1);

				y0_sum_pred  = rand("binomial", p0_sample, npergroup-ncurrent);
				y1_sum_pred  = rand("binomial", p1_sample, npergroup-ncurrent);

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
		create Results from results[colname=({lambda gamma_L gamma_U npergroup power ESS interim se})];
		do i = 1 to ncol(&npergroup.);
		    
			n      = J(     1, ncol(&interim.), .);               n[,1] = 0;
			y0_sum = J(&sims., ncol(&interim.), .);          y0_sum[,1] = 0;
			y1_sum = J(&sims., ncol(&interim.), .);          y1_sum[,1] = 0;
			posterior_prob  = J(&sims., ncol(&interim.), .);
			predictive_prob = J(&sims., ncol(&interim.), .);
			promising       = J(&sims., ncol(&interim.), .);
			futility        = J(&sims., ncol(&interim.), .);
		
			do j = 1 to ncol(&interim.);
		    
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

				if j >= 2 then do;
					* early stopping with promising result;
                    promising[,j] = ((posterior_prob[,j] >= &lambda.) | (predictive_prob[,j] >= &gamma_U.)) + promising[,1:j-1][,+];
					* early stopping with futility result;
                    futility[,j]  = (                                   (predictive_prob[,j] <= &gamma_L.)) +  futility[,1:j-1][,+];
				end;

			end;
			
			met = promising[,+] > 0;
			/* print n y0_sum y1_sum posterior_prob predictive_prob met; */
			
			* npergroup & power;
			npergroup = &npergroup.[i];
			interim = ncol(&interim.)-2;
			power = mean(met);
			se    =  std(met) / sqrt(&sims.);
			ESS     = (n ## ( promising + futility <= 1))[,+][:,];
			results = &lambda. || &gamma_L. || &gamma_U. || npergroup || power || ESS || interim || se;

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
