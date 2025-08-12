/*  ====================================================

    [Project] SASÉÜÅ[ÉUÅ[ëçâÔ2025
    [Name]    Ryo Mishima

==================================================== */

* SAS Enterprise Guide & Visual Studio Code;
%let prg=&_SASPROGRAMFILE.;
* SAS 9.4;
%*let prg=%sysget(SAS_EXECFILEPATH);
%put &=prg.;

%let dir=%substr(&prg., 1, %length(&prg.)-%length(%scan(&prg.,-1,'\'))-1);
%put &=dir.;


/* ------------------------------------------------- */

%include "&dir.\OneSampleFreq.sas"  / source2;
%include "&dir.\TwoSampleFreq.sas"  / source2;
%include "&dir.\OneSampleMeans.sas" / source2;
%include "&dir.\TwoSampleMeans.sas" / source2;

/* ------------------------------------------------- */

* number of simulations;
%global sims;
%let    sims =  500;
%global nmc;
%let    nmc  =  200;

/* ------------------------------------------------- */


title1 "One Sample Freq";

title2 "Type I error";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.30,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);
title4 "with interim analysis (0.7N)";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.30,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 0.7 1.0},
    lambda         = 0.93,
    gamma_L        = 0.10,
    gamma_U        = 0.80,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "optimistic prior";

title4 "without interim analysis";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.30,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior_a        = 3,
    prior_b        = 2,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "pessimistic prior";

title4 "without interim analysis";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.30,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior_a        = 1,
    prior_b        = 4,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title2 "Power";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.40,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);
title4 "with interim analysis (0.7N)";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.40,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 0.7 1.0},
    lambda         = 0.93,
    gamma_L        = 0.10,
    gamma_U        = 0.80,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "optimistic prior";

title4 "without interim analysis";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.40,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior_a        = 3,
    prior_b        = 2,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "pessimistic prior";

title4 "without interim analysis";
%m_Bayes_OneSampleFreq(
    nullproportion = 0.30,
    proportion     = 0.40,
    ntotal         = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior_a        = 1,
    prior_b        = 4,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

/* ------------------------------------------------- */

title1 "One Sample Means";

title2 "Type I error";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     = 0.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 1.0},
    lambda   = 0.90,
    sides    = "L",
    sims     = &sims.,
    nmc      = &nmc.
);

title4 "with interim analysis (0.7N)";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     = 0.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 0.7 1.0},
    lambda   = 0.93,
    sides    = "L",
    gamma_L  = 0.10,
    gamma_U  = 1.00,
    sims     = &sims.,
    nmc      = &nmc.
);
title4;

title3 "optimistic prior";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     = 0.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 1.0},
    lambda   = 0.90,
    prior_eta=-2.0,
    prior_tau= 2.0,
    sides    = "L",
    sims     = &sims.,
    nmc      = &nmc.
);

title3 "pessimistic prior";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     = 0.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 1.0},
    lambda   = 0.90,
    prior_eta=+1.0,
    prior_tau= 2.0,
    sides    = "L",
    sims     = &sims.,
    nmc      = &nmc.
);

title2 "Power";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     =-1.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 1.0},
    lambda   = 0.90,
    sides    = "L",
    sims     = &sims.,
    nmc      = &nmc.
);

title4 "with interim analysis (0.7N)";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     =-1.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 0.7 1.0},
    lambda   = 0.93,
    sides    = "L",
    gamma_L  = 0.20,
    gamma_U  = 0.80,
    sims     = &sims.,
    nmc      = &nmc.
);
title4;

title3 "optimistic prior";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     =-1.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 1.0},
    lambda   = 0.90,
    prior_eta=-2.0,
    prior_tau= 2.0,
    sides    = "L",
    sims     = &sims.,
    nmc      = &nmc.
);

title3 "pessimistic prior";
%m_Bayes_OneSampleMeans(
    nullmean = 0.0,
    mean     =-1.0,
    stddev   = 5.0,
    ntotal   = {20 40 60 80 100 150 200 300 500},
    interim  = {0.0 1.0},
    lambda   = 0.90,
    prior_eta=+1.0,
    prior_tau= 2.0,
    sides    = "L",
    sims     = &sims.,
    nmc      = &nmc.
);

/* ------------------------------------------------- */

title1 "Two Sample Freq";

title2 "Type I error";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);
title4 "with interim analysis (0.5N 0.7N)";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 0.5 0.7 1.0},
    lambda         = 0.95,
    gamma_L        = 0.10,
    gamma_U        = 0.80,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "optimistic prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior1_a       = 3,
    prior1_b       = 2,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "pessimistic prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior1_a       = 1,
    prior1_b       = 4,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title2 "Power";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleFreq(
    refproportion = 0.30,
    proportiondiff= 0.10,
    npergroup     = {20 40 60 80 100 150 200 300 500},
    interim       = {0.0 1.0},
    lambda        = 0.90,
    sides         = "U",
    sims           = &sims.,
    nmc            = &nmc.
);
title4 "with interim analysis (0.5N 0.7N)";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0.10,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 0.5 0.7 1.0},
    lambda         = 0.95,
    gamma_L        = 0.10,
    gamma_U        = 0.80,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "optimistic prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0.10,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior1_a       = 3,
    prior1_b       = 2,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

title3 "pessimistic prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleFreq(
    refproportion  = 0.30,
    proportiondiff = 0.10,
    npergroup      = {20 40 60 80 100 150 200 300 500},
    interim        = {0.0 1.0},
    lambda         = 0.90,
    prior1_a       = 1,
    prior1_b       = 4,
    sides          = "U",
    sims           = &sims.,
    nmc            = &nmc.
);

/* ------------------------------------------------- */

title1 "Two Sample Means";

title2 "Type I error";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = 0.0,
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 1.0},
    lambda     = 0.90,
    sides      = "L",
    sims       = &sims.,
    nmc        = &nmc.
);

title4 "with interim analysis (0.5N, 0.7N)";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = 0.0,
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 0.5 0.7 1.0},
    lambda     = 0.95,
    sides      = "L",
    gamma_L    = 0.10,
    gamma_U    = 1.00,
    sims       = &sims.,
    nmc        = &nmc.
);
title4;

title3 "optimistic prior";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = 0.0,
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 1.0},
    lambda     = 0.90,
    prior1_eta =-2.0,
    prior1_tau = 2.0,
    sides      = "L",
    sims       = &sims.,
    nmc        = &nmc.
);

title3 "pessimistic prior";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = 0.0,
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 1.0},
    lambda     = 0.90,
    prior1_eta =+1.0,
    prior1_tau = 2.0,
    sides      = "L",
    sims       = &sims.,
    nmc        = &nmc.
);

title2 "Power";

title3 "non-informative prior";

title4 "without interim analysis";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = rand("normal",-1.0,0.2),
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 1.0},
    lambda     = 0.90,
    sides      = "L",
    sims       = &sims.,
    nmc        = &nmc.
);

title4 "with interim analysis (0.5N 0.7N)";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = rand("normal",-1.0,0.2),
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 0.5 0.7 1.0},
    lambda     = 0.95,
    sides      = "L",
    gamma_L    = 0.20,
    gamma_U    = 0.80,
    sims       = &sims.,
    nmc        = &nmc.
);
title4;

title3 "optimistic prior";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = rand("normal",-1.0,0.2),
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 1.0},
    lambda     = 0.90,
    prior1_eta =-2.0,
    prior1_tau = 2.0,
    sides      = "L",
    sims       = &sims.,
    nmc        = &nmc.
);

title3 "pessimistic prior";
%m_Bayes_TwoSampleMeans(
    refmean    = 0.0,
    meandiff   = rand("normal",-1.0,0.2),
    stddev     = 5.0,
    npergroup  = {20 40 60 80 100 150 200 300 500},
    interim    = {0.0 1.0},
    lambda     = 0.90,
    prior1_eta =+1.0,
    prior1_tau = 2.0,
    sides      = "L",
    sims       = &sims.,
    nmc        = &nmc.
);

/* ------------------------------------------------- */
title4;
title3;
title2;
title1;


