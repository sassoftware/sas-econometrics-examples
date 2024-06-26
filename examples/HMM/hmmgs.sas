/*--------------------------------------------------------------

                    SAS Sample Library

        Name: hmmgs.sas
 Description: Example program from SAS Econometrics User's Guide,
              The HMM Procedure
       Title: Getting Started Example for PROC HMM
     Product: SAS Econometrics Software
        Keys: Hidden Markov Model
        PROC: HMM
       Notes:

--------------------------------------------------------------*/
ods graphics on / antialiasmax=10000;

%let pi1 = 0.5;
%let a11 = 0.001;
%let a22 = 0.001;
%let mu1_1 = 1;
%let mu1_2 = 1;
%let sigma1_11 = 1;
%let sigma1_21 = 0;
%let sigma1_22 = 1;
%let mu2_1 = -1;
%let mu2_2 = -1;
%let sigma2_11 = 1;
%let sigma2_21 = 0;
%let sigma2_22 = 1;
%let T = 10000;
%let seed = 1234;

data DGPone;
   retain cd1_11 cd1_21 cd1_22 cd2_11 cd2_21 cd2_22;
   do t = 1 to &T.;
      if(t=1) then do;
         /* initial probability distribution */
         p = &pi1.;
         /* Cholesky decomposition of sigma1 */
         cd1_11 = sqrt(&sigma1_11.);
         cd1_21 = &sigma1_21./sqrt(&sigma1_11.);
         cd1_22 = sqrt(&sigma1_22.-&sigma1_21.*&sigma1_21./&sigma1_11.);
         /* Cholesky decomposition of sigma2 */
         cd2_11 = sqrt(&sigma2_11.);
         cd2_21 = &sigma2_21./sqrt(&sigma2_11.);
         cd2_22 = sqrt(&sigma2_22.-&sigma2_21.*&sigma2_21./&sigma2_11.);
      end;
      else do;
         /* transition probability matrix */
         if(lags=1) then p = &a11.;
         else p = 1-&a22.;
      end;
      u = uniform(&seed.);
      if(u<=p) then s=1;
      else s = 2;
      e1 = normal(&seed.); e2 = normal(&seed.);
      if(s=1) then do;
         /* (x,y) ~ N(mu1, Sigma1) at state 1 */
         x = &mu1_1. + cd1_11*e1;
         y = &mu1_2. + cd1_21*e1+cd1_22*e2;
      end;
      else do;
         /* (x,y) ~ N(mu2, Sigma2) at state 2 */
         x = &mu2_1. + cd2_11*e1;
         y = &mu2_2. + cd2_21*e1+cd2_22*e2;
      end;
      output;
      lags = s;
   end;
run;

data One;
   set DGPone;
   keep t x y;
run;

proc sgplot data=One;
   series x=t y=x;
run;

proc sgplot data=One;
   series x=t y=y;
run;

proc sgplot data=One;
   scatter y=y x=x;
run;

/* classify the data points by the line y = -x */
data Cluster;
   set One;
   if(y>-x) then state = 1;
   else state = 2;
   keep t state;
run;

data clusterCheck;
   merge DGPone(in=a) Cluster(in=b);
   by t;
   if(s=state) then correct = 1;
   else correct = 0;
   if(a=b);
   keep correct;
run;

data clusterAccuracy;
   set clusterCheck;
   retain count 0 correctCount 0;
   count = count + 1;
   correctCount = correctCount + correct;
   if(count>&T.-0.5) then do;
      accuracy = correctCount / count;
      if(accuracy<0.5) then accuracy = 1 - accuracy;
      output;
   end;
   keep accuracy;
run;

proc print data = clusterAccuracy noobs; run;

   data mylib.One; set One; run;

proc hmm data=mylib.One labelSwitch=(sort=desc(mu));
   id time=t;
   model x y / type=gaussian nstate=2;
   optimize printLevel=3 printIterFreq=1;
   decode out=mylib.oneDecode;
run;

data gHmmCheck;
   merge DGPone(in=a) mylib.oneDecode(in=b);
   by t;
   if(s=state) then do; correct=1; ds = s; end;
   else do; correct=0; ds = -s; end;
   if(a=b);
   keep t x y s correct ds;
run;

%let qFlipState = -1;
data gHmmAccuracy;
   set gHmmCheck;
   retain count 0 correctCount 0;
   count = count + 1;
   correctCount = correctCount + correct;
   if(count>&T.-0.5) then do;
      accuracy = correctCount / count;
      if(accuracy<0.5) then call symputx("qFlipState",1,'G');
      else call symputx("qFlipState",0,'G');
      if(accuracy<0.5) then accuracy = 1 - accuracy;
      output;
   end;
   keep accuracy;
run;

data gHmmCheck;
   set gHmmCheck;
   if(&qFlipState.=1) then do;
      ds = -ds;
      correct = 1 - correct;
   end;
run;

proc print data = gHmmAccuracy noobs; run;

proc sgplot data=gHmmCheck;
   scatter y=y x=x / group=ds;
run;

%let nSections = 10000;
%let T = 100;
%let seed = 1234;
%let nStates = 6;
%let mu1 = 0;
%let sigma1 = 0.0625;
%let mu2 = 0;
%let sigma2 = 0.25;
%let mu3 = 0;
%let sigma3 = 1;
%let mu4 = 0;
%let sigma4 = 4;
%let mu5 = 0;
%let sigma5 = 16;
%let mu6 = 0;
%let sigma6 = 64;
%let selfTransProb = 0.95;

%macro createCstsTable(tableName);
   data &tableName.;
      do sec = 1 to &nSections.;
         state = ceil(uniform(&seed.)*&nStates.);
         do t = 1 to &T.;
            %do i = 1 %to &nStates.;
               if(state=&i.) then do;
                  y = &&mu&i.. + sqrt(&&sigma&i..)*normal(&seed.);
                  output;
               end;
            %end;
            u = uniform(&seed.);
            if(u>&selfTransProb.) then do;
               u = (u-&selfTransProb.)/(1-&selfTransProb.)*(&nStates.-1);
               state=state + ceil(u);
               if(state>&nStates.) then state = state - &nStates.;
            end;
         end;
      end;
   run;
%mend;

%createCstsTable(cstsDGP);

data mylib.cstsDGP;
   set cstsDGP;
run;

proc hmm data=mylib.cstsDGP
   outstat=mylib.cstsStat;
   id section=sec time=t;
   model y  / method=ml type=gaussian nstate=6;
   optimize ALGORITHM=activeset printlevel=3 printIterFreq=1;
   estimate out=mylib.cstsEst;
   evaluate out=mylib.cstsEval;
   decode out=mylib.cstsDecode;
   filter out=mylib.cstsFilter;
   smooth out=mylib.cstsSmooth;
   forecast out=mylib.cstsForcast;
   score outmodel=mylib.cstsModel;
run;

%let a11 = 0.98;
%let a22 = 0.95;
%let c1_1 = 0.6;
%let c1_2 = 0.4;
%let c2_1 = 0.7;
%let c2_2 = 0.3;
%let mu1_1_1 = 1;
%let mu1_1_2 = 1;
%let mu1_2_1 = -3;
%let mu1_2_2 = -3;
%let sigma1_1_11 = 1;
%let sigma1_1_21 = 0;
%let sigma1_1_22 = 1;
%let sigma1_2_11 = 1.69;
%let sigma1_2_21 = 1;
%let sigma1_2_22 = 1.69;
%let mu2_1_1 = -1;
%let mu2_1_2 = -1;
%let mu2_2_1 = 3;
%let mu2_2_2 = 3;
%let sigma2_1_11 = 1;
%let sigma2_1_21 = 0;
%let sigma2_1_22 = 1;
%let sigma2_2_11 = 1.69;
%let sigma2_2_21 = 1;
%let sigma2_2_22 = 1.69;
%let seed = 1234;
%let T = 200000;

data gmDGP;
   retain cd1_1_11 cd1_1_21 cd1_1_22
          cd1_2_11 cd1_2_21 cd1_2_22
          cd2_1_11 cd2_1_21 cd2_1_22
          cd2_2_11 cd2_2_21 cd2_2_22;
   do t = 1 to &T.;
      if(t=1) then do;
         * initial probability distribution;
         p = (1-&a22.)/((1-&a11.)+(1-&a22.));;
         * Cholesky decomposition;
         cd1_1_11 = sqrt(&sigma1_1_11.);
         cd1_1_21 = &sigma1_1_21./cd1_1_11;
         cd1_1_22 = sqrt(&sigma1_1_22.-cd1_1_21*cd1_1_21);
         cd1_2_11 = sqrt(&sigma1_2_11.);
         cd1_2_21 = &sigma1_2_21./cd1_2_11;
         cd1_2_22 = sqrt(&sigma1_2_22.-cd1_2_21*cd1_2_21);
         cd2_1_11 = sqrt(&sigma2_1_11.);
         cd2_1_21 = &sigma2_1_21./cd2_1_11;
         cd2_1_22 = sqrt(&sigma2_1_22.-cd2_1_21*cd2_1_21);
         cd2_2_11 = sqrt(&sigma2_2_11.);
         cd2_2_21 = &sigma2_2_21./cd2_2_11;
         cd2_2_22 = sqrt(&sigma2_2_22.-cd2_2_21*cd2_2_21);
      end;
      else do;
         * transition probability matrix;
         if(lags=1) then p = &a11.;
         else p = 1-&a22.;
      end;
      u = uniform(&seed.);
      if(u<=p) then s=1;
      else s = 2;
      e1 = normal(&seed.); e2 = normal(&seed.);
      if(s=1) then do;
         * choose component;
         u = uniform(&seed.);
         if(u<=&c1_1.) then do;
            * (x,y) ~ N(mu, Sigma) at state 1, component 1;
            x = &mu1_1_1. + cd1_1_11*e1;
            y = &mu1_1_2. + cd1_1_21*e1+cd1_1_22*e2;
         end;
         else do;
            * (x,y) ~ N(mu, Sigma) at state 1, component 2;
            x = &mu1_2_1. + cd1_2_11*e1;
            y = &mu1_2_2. + cd1_2_21*e1+cd1_2_22*e2;
         end;
      end;
      else do;
         * choose component;
         u = uniform(&seed.);
         if(u<=&c2_1.) then do;
            * (x,y) ~ N(mu, Sigma) at state 2, component 1;
            x = &mu2_1_1. + cd2_1_11*e1;
            y = &mu2_1_2. + cd2_1_21*e1+cd2_1_22*e2;
         end;
         else do;
            * (x,y) ~ N(mu, Sigma) at state 2, component 2;
            x = &mu2_2_1. + cd2_2_11*e1;
            y = &mu2_2_2. + cd2_2_21*e1+cd2_2_22*e2;
         end;
      end;
      output;
      lags = s;
   end;
run;

data gmTrain;
   set gmDGP(where=(t<=&T./2));
   keep t x y;
run;

data gmTest;
   set gmDGP(where=(t>&T./2));
   keep t x y;
run;

proc sgplot data=gmTrain;
   series x=t y=x;
run;

proc sgplot data=gmTrain;
   series x=t y=y;
run;

proc sgplot data=gmTrain;
   scatter y=y x=x;
run;

proc sgplot data=gmDGP(where=(t<=&T./2 and s=1));
   scatter y=y x=x;
run;

proc sgplot data=gmDGP(where=(t<=&T./2 and s=2));
   scatter y=y x=x;
run;

   data mylib.gmTrain; set gmTrain; run;
   data mylib.gmTest; set gmTest; run;

proc hmm data=mylib.gmTrain labelSwitch=(sort=none);
   id time=t;
   model x y / type=gaussianMixture nstate=2 ncomponent=2;
   optimize printLevel=3 printIterFreq=1;
   score outmodel=mylib.gmModel;
run;

proc hmm data=mylib.gmTest;
   score inmodel=mylib.gmModel;
   decode out=mylib.gmDecode;
run;

data gmHmmCheck;
   merge gmDGP(in=a where=(t>&T./2)) mylib.gmDecode(in=b);
   by t;
   if(s=state) then do; correct=1; ds = s; end;
   else do; correct=0; ds = -s; end;
   if(a=b);
   keep t x y s correct ds;
run;

%let qFlipState = -1;
data gmHmmAccuracy;
   set gmHmmCheck;
   retain count 0 correctCount 0;
   count = count + 1;
   correctCount = correctCount + correct;
   if(count>&T./2-0.5) then do;
      accuracy = correctCount / count;
      if(accuracy<0.5) then call symputx("qFlipState",1,'G');
      else call symputx("qFlipState",0,'G');
      if(accuracy<0.5) then accuracy = 1 - accuracy;
      output;
   end;
   keep accuracy;
run;

data gmHmmCheck;
   set gmHmmCheck;
   if(&qFlipState.=1) then do;
      ds = -ds;
      correct = 1 - correct;
   end;
run;

proc print data = gmHmmAccuracy noobs; run;

%let a11 = 0.98;
%let a22 = 0.95;
%let c1_1 = 0.6;
%let c1_2 = 0.4;
%let c2_1 = 0.7;
%let c2_2 = 0.3;
%let mu1_1_1 = 1;
%let mu1_1_2 = 1;
%let mu1_2_1 = -3;
%let mu1_2_2 = -3;
%let sigma1_1_11 = 1;
%let sigma1_1_21 = 0;
%let sigma1_1_22 = 1;
%let sigma1_2_11 = 1.69;
%let sigma1_2_21 = 1;
%let sigma1_2_22 = 1.69;
%let mu2_1_1 = -1;
%let mu2_1_2 = -1;
%let mu2_2_1 = 3;
%let mu2_2_2 = 3;
%let sigma2_1_11 = 1;
%let sigma2_1_21 = 0;
%let sigma2_1_22 = 1;
%let sigma2_2_11 = 1.69;
%let sigma2_2_21 = 1;
%let sigma2_2_22 = 1.69;
%let seed = 1234;
%let T = 200000;
%let nSections = 200;

data gmDGP;
   retain cd1_1_11 cd1_1_21 cd1_1_22
          cd1_2_11 cd1_2_21 cd1_2_22
          cd2_1_11 cd2_1_21 cd2_1_22
          cd2_2_11 cd2_2_21 cd2_2_22;
   do t = 1 to &T.;
      sec = ceil(t*&nSections./&T.);
      if(t=1) then do;
         * initial probability distribution;
         p = (1-&a22.)/((1-&a11.)+(1-&a22.));;
         * Cholesky decomposition;
         cd1_1_11 = sqrt(&sigma1_1_11.);
         cd1_1_21 = &sigma1_1_21./cd1_1_11;
         cd1_1_22 = sqrt(&sigma1_1_22.-cd1_1_21*cd1_1_21);
         cd1_2_11 = sqrt(&sigma1_2_11.);
         cd1_2_21 = &sigma1_2_21./cd1_2_11;
         cd1_2_22 = sqrt(&sigma1_2_22.-cd1_2_21*cd1_2_21);
         cd2_1_11 = sqrt(&sigma2_1_11.);
         cd2_1_21 = &sigma2_1_21./cd2_1_11;
         cd2_1_22 = sqrt(&sigma2_1_22.-cd2_1_21*cd2_1_21);
         cd2_2_11 = sqrt(&sigma2_2_11.);
         cd2_2_21 = &sigma2_2_21./cd2_2_11;
         cd2_2_22 = sqrt(&sigma2_2_22.-cd2_2_21*cd2_2_21);
      end;
      else do;
         * transition probability matrix;
         if(lags=1) then p = &a11.;
         else p = 1-&a22.;
      end;
      u = uniform(&seed.);
      if(u<=p) then s=1;
      else s = 2;
      e1 = normal(&seed.); e2 = normal(&seed.);
      if(s=1) then do;
         * choose component;
         u = uniform(&seed.);
         if(u<=&c1_1.) then do;
            * (x,y) ~ N(mu, Sigma) at state 1, component 1;
            x = &mu1_1_1. + cd1_1_11*e1;
            y = &mu1_1_2. + cd1_1_21*e1+cd1_1_22*e2;
         end;
         else do;
            * (x,y) ~ N(mu, Sigma) at state 1, component 2;
            x = &mu1_2_1. + cd1_2_11*e1;
            y = &mu1_2_2. + cd1_2_21*e1+cd1_2_22*e2;
         end;
      end;
      else do;
         * choose component;
         u = uniform(&seed.);
         if(u<=&c2_1.) then do;
            * (x,y) ~ N(mu, Sigma) at state 2, component 1;
            x = &mu2_1_1. + cd2_1_11*e1;
            y = &mu2_1_2. + cd2_1_21*e1+cd2_1_22*e2;
         end;
         else do;
            * (x,y) ~ N(mu, Sigma) at state 2, component 2;
            x = &mu2_2_1. + cd2_2_11*e1;
            y = &mu2_2_2. + cd2_2_21*e1+cd2_2_22*e2;
         end;
      end;
      output;
      lags = s;
   end;
run;

data gmTrain;
   set gmDGP(where=(sec<=&nSections./2));
   keep sec t x y;
run;

data gmTest;
   set gmDGP(where=(sec>&nSections./2));
   keep sec t x y;
run;

data mylib.gmTrain; set gmTrain; run;
data mylib.gmTest; set gmTest; run;

proc hmm data=mylib.gmTrain labelSwitch=(sort=none);
   id time=t section=sec;
   model x y / type=gaussianMixture nstate=2 ncomponent=2;
   optimize printLevel=3 printIterFreq=1;
   score outmodel=mylib.gmModel;
   filter out=mylib.gmFilter;
   evaluate out=mylib.gmEval;
run;

proc hmm data=mylib.gmTest;
   score inmodel=mylib.gmModel;
   decode out=mylib.gmDecode;
   filter out=mylib.gmFilter2;
   evaluate out=mylib.gmEval2;
run;

data gmHmmCheck;
   merge gmDGP(in=a where=(sec>&nSections./2)) mylib.gmDecode(in=b);
   by t;
   if(s=state) then do; correct=1; ds = s; end;
   else do; correct=0; ds = -s; end;
   if(a=b);
   keep sec t x y s correct ds;
run;

%let qFlipState = -1;
data gmHmmAccuracy;
   set gmHmmCheck;
   retain count 0 correctCount 0;
   count = count + 1;
   correctCount = correctCount + correct;
   if(count>&T./2-0.5) then do;
      accuracy = correctCount / count;
      if(accuracy<0.5) then call symputx("qFlipState",1,'G');
      else call symputx("qFlipState",0,'G');
      if(accuracy<0.5) then accuracy = 1 - accuracy;
      output;
   end;
   keep accuracy;
run;

data gmHmmCheck;
   set gmHmmCheck;
   if(&qFlipState.=1) then do;
      ds = -ds;
      correct = 1 - correct;
   end;
run;

proc print data = gmHmmAccuracy noobs; run;


%let x_lb = -8;
%let x_ub = 0;
%let pi1 = 0.5;
%let a11 = 0.95;
%let a22 = 0.95;
%let const1_1 = 0;
%let xl1_0_1_1 = 1;
%let cov1_11 = 2.56;
%let const2_1 = 4;
%let xl2_0_1_1 = 1.5;
%let cov2_11 = 4;
%let T = 4000;
%let seed = 1234;

data rsregdgp;
   retain cd1_11 cd2_11;
   do t = 1 to &T.;
      if(t=1) then do;
         /* initial probability distribution */
         p = &pi1.;
         /* Cholesky decomposition of COV1 */
         cd1_11 = sqrt(&cov1_11.);
         /* Cholesky decomposition of COV2 */
         cd2_11 = sqrt(&cov2_11.);
      end;
      else do;
         /* transition probability matrix */
         if(lags=1) then p = &a11.;
         else p = 1-&a22.;
      end;
      u = uniform(&seed.);
      if(u<=p) then s=1;
      else s = 2;
      x = &x_lb. + (&x_ub. - &x_lb.)*uniform(&seed.);
      e = normal(&seed.);
      if(s=1) then do;
         /* y ~ N(beta1_0+beta1_1*x, Sigma1) at state 1 */
         y = &const1_1. + &xl1_0_1_1 * x + cd1_11*e;
      end;
      else do;
         /* y ~ N(beta2_0+beta2_1*x, Sigma2) at state 2 */
         y = &const2_1. + &xl2_0_1_1 * x + cd2_11*e;
      end;
      output;
      lags = s;
   end;
run;

data rsreg;
   set rsregdgp;
   keep t y x;
run;


proc sgplot data=rsreg;
   scatter y=y x=x;
run;

data mylib.rsreg; set rsreg; run;
proc hmm data=mylib.rsreg;
   id time=t;
   model y = x / type=reg nstate=2 method=ml;
   optimize algorithm=activeset maxiter=256 printLevel=3 printIterFreq=1;
   learn out=mylib.mylearn;
   filter out=mylib.myfilter;
   smooth out=mylib.mysmooth;
   decode out=mylib.mydecode;
   evaluate out=mylib.myeval;
run;

data rsregCheck;
   merge rsregdgp(in=a) mylib.mydecode(in=b);
   by t;
   if(s=state) then do; correct=1; ds = s; end;
   else do; correct=0; ds = -s; end;
   if(a=b);
   keep t x y s correct ds;
run;

%let qFlipState = -1;
data rsregAccuracy;
   set rsregCheck;
   retain count 0 correctCount 0;
   count = count + 1;
   correctCount = correctCount + correct;
   if(count>&T.-0.5) then do;
      accuracy = correctCount / count;
      if(accuracy<0.5) then call symputx("qFlipState",1,'G');
      else call symputx("qFlipState",0,'G');
      if(accuracy<0.5) then accuracy = 1 - accuracy;
      output;
   end;
   keep accuracy;
run;

data rsregCheck;
   set rsregCheck;
   if(&qFlipState.=1) then do;
      ds = -ds;
      correct = 1 - correct;
   end;
run;

proc print data = rsregAccuracy noobs; run;

proc sgplot data=rsregCheck;
   scatter y=y x=x / group=ds;
run;

%let pi1 = 0.5;
%let a11 = 0.95;
%let a22 = 0.95;
%let ar1_1_1_1 = 0.8;
%let cov1_11 = 2.56;
%let ar2_1_1_1 = -0.7;
%let cov2_11 = 4;
%let T = 1000;
%let seed = 1234;

data rsardgp;
   retain cd1_11 cd2_11;
   ylag = 0;
   do t = 1 to &T.;
      if(t=1) then do;
         /* initial probability distribution */
         p = &pi1.;
         /* Cholesky decomposition of COV1 */
         cd1_11 = sqrt(&cov1_11.);
         /* Cholesky decomposition of COV2 */
         cd2_11 = sqrt(&cov2_11.);
      end;
      else do;
         /* transition probability matrix */
         if(lags=1) then p = &a11.;
         else p = 1-&a22.;
      end;
      u = uniform(&seed.);
      if(u<=p) then s=1;
      else s = 2;
      e = normal(&seed.);
      if(s=1) then do;
         /* y ~ N(beta1*ylag, Sigma1) at state 1 */
         y = &ar1_1_1_1 * ylag + cd1_11*e;
      end;
      else do;
         /* y ~ N(beta2*ylag, Sigma2) at state 2 */
         y = &ar2_1_1_1 * ylag + cd2_11*e;
      end;
      output;
      lags = s;
      ylag = y;
   end;
run;

data rsar;
   set rsardgp;
   keep t y;
run;

proc sgplot data=rsar;
   series y=y x=t;
run;

data mylib.rsar; set rsar; run;
proc hmm data=mylib.rsar;
   id time=t;
   model y / type=ar noint ylag=1 nstate=2 method=ml;
   optimize algorithm=activeset printLevel=3 printIterFreq=1;
   learn out=mylib.mylearn;
   filter out=mylib.myfilter;
   smooth out=mylib.mysmooth;
   decode out=mylib.mydecode;
   evaluate out=mylib.myeval;
run;

data rsarCheck;
   merge rsardgp(in=a) mylib.mydecode(in=b);
   by t;
   if(s=state) then do; correct=1; ds = s; end;
   else do; correct=0; ds = -s; end;
   if(a=b);
   keep t y s correct ds;
run;

%let qFlipState = -1;
data rsarAccuracy;
   set rsarCheck;
   retain count 0 correctCount 0;
   count = count + 1;
   correctCount = correctCount + correct;
   if(count>&T.-0.5) then do;
      accuracy = correctCount / count;
      if(accuracy<0.5) then call symputx("qFlipState",1,'G');
      else call symputx("qFlipState",0,'G');
      if(accuracy<0.5) then accuracy = 1 - accuracy;
      output;
   end;
   keep accuracy;
run;

data rsarCheck;
   set rsarCheck;
   if(&qFlipState.=1) then do;
      ds = -ds;
      correct = 1 - correct;
   end;
run;

proc print data = rsarAccuracy noobs; run;

proc sgplot data=rsarCheck;
   scatter y=y x=t / group=ds;
run;

data mylib.finite1;
   input t  y ;
datalines;
1   6
2   3
3   6
4   5
5   1
6   2
7   1
8   6
9   4
10  6
11  1
12  3
13  2
14  6
15  4
16  2
17  6
18  6
19  2
20  2
21  2
22  6
23  6
24  2
25  1
26  3
27  2
28  2
29  1
30  1
31  1
32  1
33  5
34  4
35  5
36  6
37  4
38  1
39  2
40  3
41  6
42  1
43  5
44  3
45  2
46  4
47  5
48  6
49  4
50  4
51  5
52  6
53  4
54  1
55  3
56  1
57  3
58  2
59  1
60  1
61  3
62  2
63  1
64  1
65  3
66  2
67  2
68  3
69  4
70  4
71  2
72  5
73  6
74  1
75  2
76  2
77  2
78  1
79  2
80  6
81  6
82  6
83  4
84  4
85  5
86  5
87  4
88  6
89  5
90  5
91  6
92  4
93  4
94  5
95  6
96  6
97  2
98  2
99  3
100 3
101 2
102 5
103 6
104 5
105 6
106 4
107 1
108 5
109 3
110 5
111 6
112 2
113 4
114 4
115 4
116 6
117 5
118 6
119 6
120 6
121 2
122 2
123 5
124 5
125 1
126 3
127 3
128 4
129 4
130 6
131 5
132 5
133 5
134 6
135 5
136 6
137 6
138 5
139 5
140 6
141 1
142 3
143 1
144 4
145 6
146 2
147 1
148 3
149 3
150 1
151 1
152 2
153 3
154 4
155 4
156 6
157 6
158 6
159 5
160 1
161 6
162 4
163 5
164 5
165 6
166 3
167 2
168 2
169 6
170 5
171 3
172 3
173 1
174 1
175 1
176 3
177 4
178 5
179 1
180 1
181 6
182 4
183 4
184 1
185 3
186 4
187 6
188 4
189 4
190 4
191 4
192 5
193 6
194 4
195 6
196 6
197 5
198 4
199 4
200 6
201 4
202 1
203 2
204 2
205 4
206 4
207 4
208 5
209 5
210 4
211 5
212 5
213 5
214 4
215 6
216 5
217 4
218 1
219 2
220 1
221 1
222 2
223 6
224 5
225 4
226 5
227 3
228 5
229 6
230 6
231 6
232 1
233 1
234 6
235 2
236 2
237 1
238 1
239 4
240 2
241 2
242 6
243 5
244 5
245 4
246 6
247 2
248 5
249 3
250 3
251 2
252 1
253 1
254 2
255 3
256 5
257 2
258 2
259 2
260 1
261 1
262 3
263 6
264 5
265 5
266 6
267 5
268 5
269 5
270 5
271 5
272 3
273 2
274 2
275 1
276 3
277 4
278 4
279 3
280 1
281 3
282 5
283 4
284 5
285 3
286 6
287 3
288 3
289 1
290 2
291 1
292 2
293 6
294 5
295 6
296 6
297 6
298 4
299 4
300 3
301 4
302 6
303 3
304 2
305 2
306 2
307 1
308 5
309 3
310 3
311 6
312 5
313 5
314 6
315 4
316 6
317 4
318 5
319 6
320 5
321 6
322 1
323 3
324 3
325 2
326 3
327 5
328 4
329 1
330 5
331 6
332 6
333 5
334 2
335 3
336 6
337 5
338 5
339 6
340 4
341 6
342 5
343 4
344 6
345 4
346 4
347 3
348 5
349 4
350 1
351 4
352 4
353 5
354 1
355 3
356 1
357 1
358 4
359 4
360 6
361 4
362 6
363 4
364 2
365 3
366 3
367 3
368 3
369 2
370 5
371 5
372 6
373 4
374 1
375 3
376 3
377 3
378 5
379 2
380 2
381 2
382 3
383 3
384 3
385 3
386 1
387 4
388 4
389 3
390 2
391 3
392 2
393 6
394 4
395 5
396 4
397 4
398 5
399 5
400 3
401 6
402 2
403 4
404 2
405 5
406 3
407 6
408 4
409 4
410 1
411 5
412 5
413 4
414 4
415 6
416 1
417 2
418 3
419 2
420 2
421 3
422 3
423 1
424 2
425 1
426 5
427 6
428 4
429 2
430 1
431 1
432 1
433 1
434 2
435 5
436 1
437 5
438 3
439 3
440 1
441 1
442 5
443 3
444 2
445 3
446 4
447 4
448 6
449 4
450 4
451 5
452 1
453 1
454 6
455 4
456 4
457 3
458 2
459 5
460 3
461 3
462 1
463 2
464 5
465 4
466 5
467 6
468 3
469 6
470 4
471 3
472 2
473 1
474 3
475 3
476 3
477 5
478 4
479 5
480 6
481 1
482 3
483 1
484 2
485 3
486 2
487 3
488 1
489 3
490 2
491 1
492 1
493 1
494 5
495 1
496 1
497 2
498 1
499 4
500 4
501 3
502 3
503 1
504 1
505 2
506 2
507 3
508 1
509 2
510 3
511 4
512 6
513 6
514 2
515 6
516 4
517 2
518 1
519 6
520 5
521 2
522 3
523 5
524 5
525 6
526 4
527 3
528 3
529 2
530 3
531 5
532 4
533 3
534 1
535 1
536 6
537 1
538 2
539 5
540 3
541 3
542 4
543 5
544 4
545 6
546 6
547 4
548 6
549 4
550 3
551 3
552 2
553 3
554 2
555 3
556 2
557 3
558 3
559 2
560 2
561 1
562 3
563 6
564 4
565 5
566 5
567 5
568 3
569 6
570 6
571 4
572 6
573 5
574 5
575 4
576 6
577 5
578 2
579 2
580 1
581 5
582 6
583 5
584 6
585 2
586 1
587 1
588 1
589 4
590 4
591 5
592 4
593 1
594 3
595 3
596 6
597 1
598 1
599 6
600 2
601 4
602 3
603 3
604 1
605 2
606 1
607 2
608 5
609 4
610 6
611 4
612 2
613 4
614 3
615 1
616 3
617 5
618 6
619 4
620 6
621 3
622 1
623 4
624 4
625 2
626 4
627 1
628 2
629 6
630 6
631 3
632 2
633 2
634 3
635 5
636 5
637 5
638 3
639 5
640 5
641 2
642 3
643 2
644 1
645 1
646 3
647 2
648 5
649 3
650 3
651 3
652 4
653 2
654 2
655 5
656 5
657 4
658 5
659 6
660 5
661 4
662 1
663 6
664 4
665 5
666 6
667 1
668 5
669 6
670 6
671 4
672 4
673 4
674 4
675 6
676 6
677 3
678 5
679 2
680 6
681 6
682 3
683 5
684 6
685 6
686 1
687 2
688 2
689 1
690 6
691 1
692 3
693 2
694 3
695 2
696 3
697 1
698 3
699 4
700 5
701 6
702 4
703 3
704 5
705 5
706 5
707 5
708 5
709 6
710 3
711 2
712 3
713 4
714 5
715 5
716 3
717 2
718 2
719 1
720 4
721 4
722 6
723 6
724 5
725 4
726 6
727 6
728 6
729 6
730 6
731 6
732 1
733 3
734 5
735 4
736 2
737 3
738 1
739 4
740 3
741 1
742 1
743 2
744 1
745 3
746 1
747 6
748 4
749 2
750 6
751 5
752 4
753 5
754 5
755 6
756 5
757 4
758 2
759 2
760 1
761 1
762 6
763 4
764 5
765 6
766 5
767 5
768 6
769 1
770 1
771 3
772 2
773 5
774 4
775 3
776 4
777 6
778 4
779 5
780 5
781 5
782 4
783 6
784 6
785 4
786 6
787 5
788 4
789 3
790 6
791 5
792 4
793 1
794 3
795 1
796 1
797 1
798 4
799 4
800 6
801 5
802 2
803 2
804 2
805 6
806 4
807 6
808 6
809 5
810 5
811 6
812 6
813 5
814 3
815 1
816 6
817 4
818 6
819 5
820 4
821 3
822 3
823 1
824 3
825 3
826 3
827 2
828 2
829 3
830 3
831 1
832 6
833 5
834 4
835 6
836 1
837 1
838 2
839 6
840 6
841 6
842 5
843 4
844 1
845 3
846 3
847 4
848 4
849 6
850 4
851 3
852 6
853 5
854 1
855 5
856 5
857 3
858 2
859 2
860 5
861 2
862 4
863 1
864 6
865 4
866 2
867 1
868 3
869 3
870 1
871 2
872 5
873 5
874 4
875 1
876 3
877 6
878 5
879 5
880 5
881 5
882 3
883 1
884 4
885 6
886 6
887 3
888 2
889 1
890 1
891 1
892 1
893 4
894 1
895 2
896 1
897 1
898 6
899 4
900 4
901 2
902 1
903 6
904 2
905 3
906 3
907 3
908 2
909 4
910 5
911 1
912 3
913 3
914 1
915 5
916 4
917 2
918 2
919 6
920 5
921 5
922 4
923 5
924 2
925 1
926 3
927 1
928 3
929 1
930 3
931 3
932 1
933 2
934 3
935 2
936 1
937 5
938 4
939 5
940 4
941 1
942 2
943 5
944 6
945 4
946 1
947 1
948 1
949 5
950 3
951 1
952 6
953 2
954 6
955 2
956 4
957 6
958 4
959 4
960 5
961 2
962 1
963 1
964 3
965 2
966 2
967 3
968 2
969 4
970 6
971 2
972 3
973 6
974 6
975 4
976 2
977 1
978 5
979 5
980 4
981 3
982 1
983 1
984 3
985 2
986 2
987 5
988 6
989 6
990 2
991 1
992 1
993 2
994 3
995 6
996 6
997 5
998 5
999 2
1000    6
;


proc hmm data=mylib.finite1;
   id time=t;
   model y / type=finite nstate=1 nCategory=6;
   estimate outall=mylib.myEstimate;
run;

proc iml;
  use mylib.myEstimate;
  read all var {Estimate cov_2 cov_3 cov_4 cov_5 cov_6} into m;
  close mylib.myEstimate;

cov = m[2:6,2:6];
invcov = inv(cov);
v = m[2:6,1];
c = 1/6 // 1/6 // 1/6 // 1/6 // 1/6 ;
v1 = v-c;

stat = v1`*invcov *v1;
df = 5;
pvalue = 1- cdf('CHISQ',stat,df);
varname={"stat" "df" "pvalue"};
create table var varname;
append;
close table;
quit;

proc print data=table noobs label;
   label stat="Stat" df = "DF" pvalue = "Pr > ChiSq";
run;

proc hmm data=mylib.finite1;
   id time=t;
   model y / type=finite nstate=2 nCategory=6;
   score outmodel=mylib.myScore;
run;

data mylib.finite2;
   input t  y ;
datalines;
1   6
2   4
3   6
4   1
5   1
6   1
7   1
8   4
9   3
10  6
11  5
12  1
13  5
14  5
15  5
16  4
17  6
18  4
19  6
20  6
21  5
22  3
23  1
24  2
25  3
26  2
27  2
28  6
29  6
30  4
31  2
32  3
33  3
34  6
35  6
36  2
37  1
38  3
39  5
40  1
41  1
42  3
43  1
44  5
45  6
46  6
47  6
48  2
49  1
50  3
51  2
52  1
53  2
54  1
55  2
56  3
57  2
58  2
59  1
60  1
61  2
62  2
63  2
64  3
65  3
66  1
67  3
68  2
69  2
70  1
71  6
72  6
73  2
74  3
75  4
76  3
77  2
78  5
79  6
80  5
81  3
82  5
83  3
84  3
85  4
86  6
87  6
88  4
89  5
90  6
91  4
92  4
93  6
94  6
95  5
96  2
97  1
98  6
99  4
100 6
;


proc hmm data=mylib.finite2;
   id time=t;
   score inmodel=mylib.myScore;
   forecast out=mylib.myforecast online;
run;
data mylib.forecast;
   set mylib.myforecast;
   if (Category1+Category2+Category3>0.5) then bet = 0;
   else bet = 1;
run;
data mylib.forecast;
   set mylib.forecast;
   t=t+step; *so t equals the time when forecast is on;
 run;
data mylib.finite2;
   set mylib.finite2;
 if(y<=3) then big=0;
 else big=1;
run;
data fnHmmCheck;
   merge mylib.finite2(in=a) mylib.forecast(in=b);
   by t;
run;
data fnHmmCheck;
   set fnHmmCheck;
   *for first round bet according to inital probability;
   if(t=1) then bet=1;
   if(big=bet) then gain=1;
   else gain=-1;
   keep t y gain;
run;
data totalgain;
   set fnHmmCheck;
   Totalgain + gain;
   if (t=100) then do
      Zscore = Totalgain/10;
      pvalue = 2*(1- cdf('NORMAL',Zscore));
      output totalgain;
   end;
   keep Totalgain Zscore pvalue;
run;
proc print data = totalgain noobs label;
   label totalgain = "Total gain" Zscore="Z-score" pvalue="Pr > |z|";
run;

%let pi1 = 0.5;*initial probabilty in state 1;
%let a11 = 0.75;
%let a22 = 0.75;

%macro simData(name, T, seed);
data mylib.&name.(keep = t y);
  p = &pi1.;
  call streaminit(&seed);
  do t = 1 to &T.;
   if(t=1) then
      /* initial probability distribution */
     p = &pi1.;
   else
      /* transition probability matrix */
     if(lagx=1) then p = &a11.;
     else p = 1-&a22.;
     u = rand("Uniform");
     if(u<=p) then x=1;
     *Using the random value u to control regimes, set x=number of regimes;
     else x = 2;
     if(x=1) then
     do;
      /* at state 1 */
     y = rand("Table", 1/3, 1/3, 1/3, 0 , 0 , 0);
   end;
   else do;
      /* at state 2 */
     y = rand("Table", 0 , 0 , 0 , 1/3, 1/3, 1/3);
   end;
   output;
   lagx=x;
end;
run;
%mend;

%simData(finite1,1000, 3);
%simData(finite2,100, 4);


%let pi1 = 0.7;
%let a11 = 0.9;
%let a22 = 0.8;
%let lambda1 = 5;
%let lambda2 = 10;

%macro dgpPsHMM(table,T,nSections,seed);
data &table;
   call streaminit(&seed.);
   secLen = ceil(&T./&nSections.);
   do t = 1 to &T.;
      sec = floor((t-1)/secLen)+1;
      if(t=1) then do;

   /* initial probability distribution */
      p = &pi1.;
   end;
   else do;
   /* transition probability matrix */
   if(lags=1) then p = &a11.;
   else p = 1-&a22.;
   end;
   u = rand("Uniform");
   if(u<=p) then s=1;
   else s = 2;
   if(s=1) then do;
   /* at state 1 */
      y = rand("Poisson",&lambda1.);
   end;
   else do;
   /* at state 2 */
      y = rand("Poisson",&lambda2.);
   end;
   output;
   lags = s;
   end;
run;
proc sort data=&table.; by sec t; run;
%mend;

%let T = 1000 ;
%let nSections = 2;
%let seed = 1234;
%dgpPsHMM(hmmpsdt,&T.,&nSections.,&seed.);

data mylib.hmmpsdt;
   set hmmpsdt;
run;

proc hmm data=mylib.hmmpsdt;
   id time=t;
   model y / type=Poisson nstate=2;
run;

   %let a11 = 0.98;
   %let a22 = 0.98;
   %let c1_1 = 0.3;
   %let c1_2 = 0.4;
   %let c1_3 = 0.3;
   %let c2_1 = 0.3;
   %let c2_2 = 0.4;
   %let c2_3 = 0.3;
   %let mu1_1_1 = -1;
   %let mu1_2_1 = -6;
   %let mu1_3_1 = -11;
   %let sigma1_1_11 = 1;
   %let sigma1_2_11 = 2;
   %let sigma1_3_11 = 3;
   %let mu2_1_1 = 6;
   %let mu2_2_1 = 12;
   %let mu2_3_1 = 20;
   %let sigma2_1_11 = 4;
   %let sigma2_2_11 = 5;
   %let sigma2_3_11 = 6;
   %let seed = 1234;
   %let T = 400000;
   %let nSections = 200;

   data labelSwitchDGP;
      do t = 1 to &T.;
         sec = ceil(t*&nSections./&T.);
         if(t=1) then do;
            * initial probability distribution;
            p = (1-&a22.)/((1-&a11.)+(1-&a22.));;
         end;
         else do;
            * transition probability matrix;
            if(lags=1) then p = &a11.;
            else p = 1-&a22.;
         end;
         u = uniform(&seed.);
         if(u<=p) then s=1;
         else s = 2;
         e1 = normal(&seed.);
         if(s=1) then do;
            * choose component;
            u = uniform(&seed.);
            if(u<=&c1_1.) then do;
               * x ~ N(mu, Sigma) at state 1, component 1;
               x = &mu1_1_1. + sqrt(&sigma1_1_11.)*e1;
            end;
            else do;
               if(u<=&c1_1.+&c1_2.) then do;
                  * x ~ N(mu, Sigma) at state 1, component 2;
                  x = &mu1_2_1. + sqrt(&sigma1_2_11.)*e1;
               end;
               else do;
                  * x ~ N(mu, Sigma) at state 1, component 3;
                  x = &mu1_3_1. + sqrt(&sigma1_3_11.)*e1;
               end;
            end;
         end;
         else do;
            * choose component;
            u = uniform(&seed.);
            if(u<=&c2_1.) then do;
               * x ~ N(mu, Sigma) at state 2, component 1;
               x = &mu2_1_1. + sqrt(&sigma2_1_11.)*e1;
            end;
            else do;
               if(u<=&c1_1.+&c1_2.) then do;
                  * x ~ N(mu, Sigma) at state 2, component 2;
                  x = &mu2_2_1. + sqrt(&sigma2_2_11.)*e1;
               end;
               else do;
                  * x ~ N(mu, Sigma) at state 2, component 3;
                  x = &mu2_3_1. + sqrt(&sigma2_3_11.)*e1;
               end;
            end;
         end;
         output;
         lags = s;
      end;
   run;

proc sgplot data=labelSwitchDGP;
   histogram x;
run;

   data mylib.labelSwitchDGP; set labelSwitchDGP; run;

proc hmm data=mylib.labelSwitchDGP;
   id time=t section=sec;
   model x / type=gaussian nstate=6;
   optimize printLevel=3 printIterFreq=1 algorithm=interiorpoint;
   initial mu={20, 12, 6, -1, -6, -11};
run;

proc hmm data=mylib.labelSwitchDGP labelSwitch=(out=mylib.olabel1);
   id time=t section=sec;
   model x / type=gaussian nstate=6;
   optimize printLevel=3 printIterFreq=1 algorithm=interiorpoint;
   initial mu={20, 12, 6, -1, -6, -11};
run;

data olabel1; set mylib.olabel1; run;
proc sort data=olabel1; by oldStateLabel; run;
proc print data=olabel1 noobs label; run;

proc hmm data=mylib.labelSwitchDGP labelSwitch=(sort=asc(mu) out=mylib.olabel2);
   id time=t section=sec;
   model x / type=gaussian nstate=6;
   optimize printLevel=3 printIterFreq=1 algorithm=interiorpoint;
   initial mu={20, 12, 6, -1, -6, -11};
run;

proc hmm data=mylib.labelSwitchDGP labelSwitch=(out=mylib.olabel3);
   id time=t section=sec;
   model x / type=gaussianMixture nstate=2 ncomponent=3;
   optimize printLevel=3 printIterFreq=1 algorithm=interiorpoint;
   initial mu={20, 12, 6, -1, -6, -11};
   score outmodel=mylib.omodel;
run;

data olabel3; set mylib.olabel3; run;
proc sort data=olabel3; by oldStateLabel oldComponentLabel; run;
proc print data=olabel3 noobs label; run;

proc hmm data=mylib.labelSwitchDGP
         labelSwitch=(out=mylib.olabel4);;
   score inmodel=mylib.omodel;
run;

data olabel4; set mylib.olabel4; run;
proc sort data=olabel4; by oldStateLabel oldComponentLabel; run;
proc print data=olabel4 noobs label; run;

proc hmm data=mylib.labelSwitchDGP
         labelSwitch=(reorder=(2 1) reorderComponent=(2 1 3 2 3 1)
                      out=mylib.olabel5);
   score inmodel=mylib.omodel;
run;

data olabel5; set mylib.olabel5; run;
proc sort data=olabel5; by oldStateLabel oldComponentLabel; run;
proc print data=olabel5 noobs label; run;

%let pi1 = 0.5;
%let a11 = 0.2;
%let a22 = 0.2;
%let mu1_1 = 1;
%let mu1_2 = 1;
%let sigma1_11 = 1;
%let sigma1_21 = 0;
%let sigma1_22 = 1;
%let mu2_1 = -1;
%let mu2_2 = -1;
%let sigma2_11 = 1;
%let sigma2_21 = 0;
%let sigma2_22 = 1;

%macro dgpGHMM(table,T,nSections,seed);
   data &table;
      retain cd1_11 cd1_21 cd1_22 cd2_11 cd2_21 cd2_22;
      secLen = ceil(&T./&nSections.);
      do t = 1 to &T.;
         sec = floor((t-1)/secLen)+1;
         if(t=1) then do;
            /* initial probability distribution */
            p = &pi1.;
            /* Cholesky decomposition of sigma1 */
            cd1_11 = sqrt(&sigma1_11.);
            cd1_21 = &sigma1_21./sqrt(&sigma1_11.);
            cd1_22 = sqrt(&sigma1_22.-&sigma1_21.*&sigma1_21./&sigma1_11.);
            /* Cholesky decomposition of sigma2 */
            cd2_11 = sqrt(&sigma2_11.);
            cd2_21 = &sigma2_21./sqrt(&sigma2_11.);
            cd2_22 = sqrt(&sigma2_22.-&sigma2_21.*&sigma2_21./&sigma2_11.);
         end;
         else do;
            /* transition probability matrix */
            if(lags=1) then p = &a11.;
            else p = 1-&a22.;
         end;
         u = uniform(&seed.);
         if(u<=p) then s=1;
         else s = 2;
         e1 = normal(&seed.); e2 = normal(&seed.);
         if(s=1) then do;
            /* (x,y) ~ N(mu1, Sigma1) at state 1 */
            x = &mu1_1. + cd1_11*e1;
            y = &mu1_2. + cd1_21*e1+cd1_22*e2;
         end;
         else do;
            /* (x,y) ~ N(mu2, Sigma2) at state 2 */
            x = &mu2_1. + cd2_11*e1;
            y = &mu2_2. + cd2_21*e1+cd2_22*e2;
         end;
         output;
         lags = s;
      end;
   run;
   proc sort data=&table.; by sec t; run;
%mend;

* prepare the training data;
%let T = 10000;
%let nSections = 10;
%let seed = 1234;
%dgpGHMM(gHmmTrain,&T.,&nSections.,&seed.);
data mylib.gHmmTrain; set gHmmTrain; run;

* estimate the model and save it in an analytic store table;
proc hmm data=mylib.gHmmTrain labelSwitch=(sort=desc(mu));
   id time=t section=sec;
   model x y / type=gaussian nstate=2;
   store out=mylib.ghmmModel;
run;

* produce a local analytic store file;
proc astore;
    download rstore=mylib.ghmmModel store="ls_ghmm";
run;

* get new data;
%let T = 40;
%let nSections = 2;
%let seed = 12345;
%dgpGHMM(gHmmTest,&T.,&nSections.,&seed.);

* score new data;
proc astore;
    setoption decode 1;
    setoption decode_window 8;
    setoption evaluate 1;
    setoption filter 1;
    setoption forecast 1;
    setoption forecast_alpha 0.8;
    setoption forecast_cov 1;
    setoption forecast_lead 4;
    setoption smooth 1;
    setoption smooth_window 12;
    describe store="ls_ghmm";
    score data=gHmmTest store="ls_ghmm" out=ghmmOut;
run;

* get new data;
%let T = 40;
%let nSections = 2;
%let seed = 12345;
%dgpGHMM(gHmmTest,&T.,&nSections.,&seed.);
data mylib.gHmmTest; set gHmmTest; run;

* upload the local analytic store;
proc astore;
    upload store="ls_ghmm" rstore=mylib.rs_ghmm;
run;

* score new data;
proc astore;
    setoption decode 1;
    setoption decode_window 8;
    setoption evaluate 1;
    setoption filter 1;
    setoption forecast 1;
    setoption forecast_alpha 0.8;
    setoption forecast_cov 1;
    setoption forecast_lead 4;
    setoption smooth 1;
    setoption smooth_window 12;
    describe rstore=mylib.rs_ghmm;
    score data=mylib.gHmmTest rstore=mylib.rs_ghmm out=mylib.ghmmCasOut;
run;

data one;
   call streaminit(12345);
   do sec = 1 to 10000;
      t = 1;
      u = rand('UNIFORM');
      if(u<=0.75) then do;
         s = 1;
         x = 2 + rand('NORMAL');
         y = 2 + rand('NORMAL');
      end;
      else do;
         s = 2;
         x = -2 + rand('NORMAL');
         y = -2 + rand('NORMAL');
      end;
      output;
   end;
run;

data mylib.one; set one; run;

proc sgplot data=one;
   scatter x=x y=y / group=s;
run;

proc hmm data=mylib.one labelSwitch=(sort=desc(mu));
   id time=t section=sec;
   model x y / nstate=2 type=gaussian estispv;
   decode out=mylib.mydecode;
   store mylib.mygmm;
run;

data mydecode; set mylib.mydecode; run;
proc sort data=mydecode; by sec; run;

data combine;
   merge mydecode(in=a) one(in=b);
   by sec;
run;

data checkAccuracy;
   set combine;
   retain accuracy 0;
   if(state=s) then accuracy+1;
   if(_N_=10000) then do;
      accuracy = accuracy/_N_;
      output;
   end;
   keep accuracy;
run;

proc print data=checkAccuracy noobs; run;

proc sgplot data=combine;
   scatter x=x y=y / group=state;
run;

data two;
   call streaminit(12346);
   do sec = 1 to 10000;
      t = 1;
      u = rand('UNIFORM');
      if(u<=0.75) then do;
         s = 1;
         x = 2 + rand('NORMAL');
         y = 2 + rand('NORMAL');
      end;
      else do;
         s = 2;
         x = -2 + rand('NORMAL');
         y = -2 + rand('NORMAL');
      end;
      output;
   end;
run;

proc astore;
   download rstore=mylib.mygmm store='ls_gmm';
   setoption decode 1;
   setoption decode_window 1;
   setoption evaluate 0;
   describe store='ls_gmm';
   score data=two store='ls_gmm' out=gmmout;
run;

data combine2;
   merge gmmout(in=a) two(in=b);
   by sec;
run;

data checkAccuracy2;
   set combine2;
   retain accuracy 0;
   if(DecodedState_BackStep0=s) then accuracy+1;
   if(_N_=10000) then do;
      accuracy = accuracy/_N_;
      output;
   end;
   keep accuracy;
run;

proc print data=checkAccuracy2 noobs; run;

