/*--------------------------------------------------------------

                    SAS Sample Library

        Name: hmmex01.sas
 Description: Example program from SAS Econometrics User's Guide,
              The HMM Procedure
       Title: Discovering the Hidden Market States
     Product: SAS Econometrics Software
        Keys: Hidden Markov Model
        PROC: HMM
       Notes:
--------------------------------------------------------------*/

ods graphics on;

%let ds = vwmi;
%let cutDate =  '31DEC2000'd;
data &ds. &ds.In;
   set &ds.O;
   retain cumReturn 0;
   return = (log(price)-log(lag(price)))*100;
   if(return~=.) then cumReturn + return;
   if(mod(_N_,5)=1 and _N_>1) then do;
      returnw = cumReturn;
      w + 1;
      output &ds.;
      if (date<=&cutDate.) then output &ds.In;
      cumReturn = 0;
   end;
   keep w date returnw price;
   label w="Week Index" date="Date" returnw="Weekly Return"
         price="Market Index";
run;

proc sgplot data=&ds.;
   series x=date y=returnw;
   series x=date y=price / y2axis;
   refline 0 / axis=y lineattrs=(color=gray thickness=1);
   refline &cutDate. / axis=x
                       lineattrs=(color=gray thickness=2 pattern=dash);
   yaxis max=20 min=-20;
run;

   data mylib.&ds.In; set &ds.In; run;
   data mylib.&ds.; set &ds.; run;

***************************************************************************;
* ANALYSIS WITH A SIMPLE MODEL                                             ;
***************************************************************************;

proc hmm data=mylib.&ds.In;
   id time=date;
   model returnw / type=ar nstate=2 ylag=0;
   smooth out=mylib.&ds.inSmooth_k2To2_p0To0;
   decode out=mylib.&ds.inDecode_k2To2_p0To0;
   evaluate out=mylib.&ds.inEval_k2To2_p0To0;
run;

*print sample mean and standard error for each state
 according to the weighted sample,
 where the smoothed probabilities are the weight;
data &ds.inSmooth_k2to2_p0to0;
   set mylib.&ds.inSmooth_k2to2_p0to0;
run;
proc sort data=&ds.inSmooth_k2to2_p0to0; by date; run;
data &ds.inWeighted;
   merge &ds.inSmooth_k2to2_p0to0(in=a) &ds.in(in=b);
   by date;
   if(a=b);
run;
%macro printMeanAndStdErrOfAState(k);
   data &ds.inMeanStdErr;
      set &ds.inWeighted end=eof;
      array cWeights{&k.} _TEMPORARY_;
      array cMeans[&k.] mean1-mean&k.;
      array cVars[&k.] sd1-sd&k.;
      array weights[&k.] state1-state&k.;
      retain cWeights cMeans cVars;
      if(_N_=1) then do;
         do i = 1 to &k.;
            cWeights[i] = weights[i];
            cMeans[i] = weights[i]*returnw;
            cVars[i] = weights[i]*returnw*returnw;
         end;
      end;
      else do;
         do i = 1 to &k.;
            cWeights[i] = cWeights[i] + weights[i];
            cMeans[i] = cMeans[i] + weights[i]*returnw;
            cVars[i] = cVars[i] + weights[i]*returnw*returnw;
         end;
      end;
      if eof then do;
         do i = 1 to &k.;
            if(cWeights[i]>0) then do;
               cMeans[i] = cMeans[i] / cWeights[i];
               cVars[i] = cVars[i] / cWeights[i] - cMeans[i]*cMeans[i];
               if(cVars[i]>=0) then cVars[i] = sqrt(cVars[i]);
            end;
            else do;
               cMeans[i] = .;
               cVars[i] = .;
            end;
         end;
         %do i = 1 %to &k.;
            call symputx("mu&i.",cMeans[&i.],'G');
            call symputx("sigma&i.",cVars[&i.],'G');
            state = &i.;  mean = mean&i.; stdErr = sd&i.;
            output;
         %end;
      end;
      label state="State" mean="Sample Mean"
            stdErr="Sample Standard Error";
      keep state mean stdErr;
   run;
   proc print data=&ds.inMeanStdErr noobs label; run;
%mend;
%printMeanAndStdErrOfAState(2);

%macro plotState(myData,myColumn,ks,ke);
   proc sgplot data=&myData.;
      refline 0 / axis=x lineattrs=(thickness=3);
      histogram &myColumn. / nbins=200;
      %do i = &ks. %to &ke.;
         density &myColumn. / type=normal(mu=&&mu&i. sigma=&&sigma&i.)
            name="state&i."
            legendlabel="State &i.";
      %end;
      keylegend %do i = &ks. %to &ke.; "state&i." %end;;
      xaxis min=-20 max=20 ;
      yaxis min=0 max=15;
   run;
%mend;
%plotState(&ds.In,returnw,1,2);

* plot decoded results;
data &ds.InDecode_k2to2_p0to0; set mylib.&ds.InDecode_k2to2_p0to0; run;
proc sort data=&ds.InDecode_k2to2_p0to0; by date; run;
data decodeMaketStates;
   set &ds.In;
   set &ds.InDecode_k2to2_p0to0;
   array returnws returnw1-returnw2;
   array prices price1-price2;
   if(state~=.) then do;
      if(state=1) then do;
         returnws[1] = returnw; prices[1] = price;
         marketState="Bull Market";
      end;
      if(state=2) then do;
         returnws[2] = returnw; prices[2] = price;
         marketState="Bear Market";
      end;
   end;
   label marketState="Market State";
run;
proc sgplot data=decodeMaketStates;
   styleattrs datacontrastcolors=(green red);
   scatter x=date y=marketState / group=marketState;
   yaxis label="Market State" offsetmin=0.1 offsetmax=0.1 reverse;
run;

proc sgplot data=decodeMaketStates;
   series x=date y=returnw1 / break lineattrs=(color=green)
                              legendlabel="Bull Market";
   series x=date y=returnw2 / break lineattrs=(color=red)
                              legendlabel="Bear Market";
   refline 0 / axis=y lineattrs=(color=gray thickness=1);
   yaxis max=20 min=-20 label="Weekly Return";
run;

proc sgplot data=decodeMaketStates;
   series x=date y=price1 / break lineattrs=(color=green)
                            legendlabel="Bull Market";
   series x=date y=price2 / break lineattrs=(color=red)
                            legendlabel="Bear Market";
   yaxis label="Market Index";
run;

* plot conditional log likelihood;
%macro plotEvaluation(inds, outds, qConditional);
   data &outds.;
      set &inds.;
      format date MMDDYY10.;
   run;
   %if (&qConditional.=1) %then %do;
      proc sort data=&outds.;
         by date;
      run;
      data &outds.;
         set &outds.;
         if(_N_=1) then
            condLogLikelihood = logLikelihood;
         else
            condLogLikelihood = logLikelihood - lag(logLikelihood);
      run;
      proc sgplot data=&outds.;
         series x=date y=condLogLikelihood;
         yaxis label='Conditional Log Likelihood' min=-15 max=0;
      run;
   %end;
   %else %do;
      proc sgplot data=&outds.;
         series x=date y=logLikelihood;
         yaxis label='Log Likelihood';
      run;
   %end;
%mend plotEvaluation;

%plotEvaluation(mylib.&ds.inEval_k2to2_p0to0, &ds.inEval_k2to2_p0to0, 1);

***************************************************************************;
* MODEL SELECTION                                                          ;
***************************************************************************;

%macro estimateRSAR(myds, inEstDs, kStart, kEnd, pStart, pEnd, method,
                    maxiter, qMultiStart);
   proc hmm data=&myds. labelSwitch=(sort=none)
      outstat=&myds.&method.Stat_k&kStart.To&kEnd._p&pStart.To&pEnd.;
      id time=date;
      model returnw / type=ar nstate=&kStart.:&kEnd. ylag=&pStart.:&pEnd.
         method=&method.;
      optimize printLevel=3 printIterFreq=1 algorithm=interiorpoint
         maxiter=&maxiter. Multistart=&qMultiStart.;
      score outmodel=&myds&method.Model_k&kStart.To&kEnd._p&pStart.To&pEnd.;
      learn out=&myds.&method.Learn_k&kStart.To&kEnd._p&pStart.To&pEnd.
         %if %length(&inEstDs.)>0 %then %do; in=&inEstDs. %end;
         ;
      display / excludeall;
   run;
%mend;

* (1) MAP + MULTISTART for AR(0);
* be aware that the following macro might take tens of hours to finish;
* uncomment it to run;
* %estimateRSAR(myds=mylib.&ds.In, inEstDs=,
   kStart=2, kEnd=9, pStart=0, pEnd=0,
   method=MAP, maxiter=128, qMultiStart=1);

* (2) ML + initial values estimated from MAP's AR(0) for AR(0)
      and initial values estimated from AR(p-1) for AR(p), p>0;
* be aware that the following macro might take tens of minutes to finish;
* uncomment it to run;
* %estimateRSAR(myds=mylib.&ds.In, inEstDs=mylib.&ds.InMAPLearn_k2to9_p0to0,
   kStart=2, kEnd=9, pStart=0, pEnd=5,
   method=ML, maxiter=128, qMultiStart=0);

* print fit statistics in matrix form;
data &ds.InMLStat_k2to9_p0to5; set mylib.&ds.InMLStat_k2to9_p0to5; run;
proc sort data=&ds.InMLStat_k2to9_p0to5; by modelIndex; run;
%macro printFitStat(fitStat);
   data &ds.In&fitStat.;
      set &ds.InMLStat_k2to9_p0to5 end=eof;
      array &fitStat.s{8,6} _TEMPORARY_;
      array yLags{6} yLag0-yLag5;
      &fitStat.s[nState-1,yLag+1] = &fitStat.;
      if eof then do;
         do i = 1 to 8;
            nState = i+1;
            do j = 1 to 6; yLags(j) = &fitStat.s[i,j]; end;
            output;
         end;
      end;
      keep nState yLag0-yLag5;
   run;
   proc print data=&ds.In&fitStat. label noobs
      style(header)={textalign=center};
      var nState yLag0-yLag5;
      label nState='k' yLag0='p = 0' yLag1='p = 1' yLag2='p = 2'
         yLag3='p = 3' yLag4='p = 4' yLag5='p = 5';
   run;
%mend;
%printFitStat(AIC);

* Uncomment the following macros if you want to see other fit statistics;
* %printFitStat(logLikelihood);
* %printFitStat(AICC);
* %printFitStat(BIC);
* %printFitStat(HQC);

***************************************************************************;
* THE SELECTED MODEL AND IN-SAMPLE INFERENCE                               ;
***************************************************************************;

* select eight-state RS-AR(4) model according to AIC or AICC;
* adjust data to be the same sample size for AR(5);
proc hmm data=mylib.&ds.In(where=(date>='30JAN1926'd));
   id time=date;
   model returnw / type=ar nstate=8 ylag=4;
   optimize printLevel=3 algorithm=interiorpoint maxiter=0;
   score outmodel=mylib.&ds.InMLModel_k8To8_p4To4;
   learn out=mylib.&ds.InMLLearn_k8To8_p4To4
         in=mylib.&ds.InMLLearn_k2to9_p0to5;
   smooth out=mylib.&ds.inSmooth_k8to8_p4to4;
   decode out=mylib.&ds.inDecode_k8to8_p4to4;
   evaluate out=mylib.&ds.inEval_k8to8_p4to4;
   display tpm;
run;

%plotEvaluation(mylib.&ds.inEval_k8to8_p4to4, &ds.inEval_k8to8_p4to4, 1);

* print observation parameters;
data &ds.InObsParm;
   set mylib.&ds.InMLLearn_k8To8_p4To4(where=(type='EST' and index>=17));
   if(index=17) then parameter='Constant';
   if(index=19) then parameter='AR(1)';
   if(index=21) then parameter='AR(2)';
   if(index=23) then parameter='AR(3)';
   if(index=25) then parameter='AR(4)';
   if(index=27) then parameter='Variance';
run;
proc sort data=&ds.InObsParm; by index; run;
proc print data=&ds.InObsParm noobs label;
   var parameter state1-state8;
   label state1='State 1' state2='State 2' state3='State 3'
         state4='State 4' state5='State 5' state6='State 6'
         state7='State 7' state8='State 8';
run;

data &ds.inSmooth_k8to8_p4to4;
   set mylib.&ds.inSmooth_k8to8_p4to4(where=(date>='25FEB1926'd));
run;
proc sort data=&ds.inSmooth_k8to8_p4to4; by date; run;
data &ds.inTruncated; set &ds.in(where=(date>='25FEB1926'd)); run;
data &ds.inWeighted;
   merge &ds.inSmooth_k8to8_p4to4(in=a) &ds.inTruncated(in=b);
   by date;
   if(a=b);
run;
%printMeanAndStdErrOfAState(8);

%plotState(&ds.In,returnw,1,8);

%plotState(&ds.In,returnw,1,3);

%plotState(&ds.In,returnw,4,6);

%plotState(&ds.In,returnw,7,8);

* plot decoded results;
data &ds.InDecode_k8to8_p4to4; set mylib.&ds.InDecode_k8to8_p4to4; run;
proc sort data=&ds.InDecode_k8to8_p4to4; by date; run;
data decodeMaketStates;
   set &ds.In(where=(date>='30JAN1926'd));
   set &ds.InDecode_k8to8_p4to4;
   array returnws returnw1-returnw3;
   array prices price1-price3;
   if(state~=.) then do;
      if(state=1 | state=2 | state=3) then do;
         returnws[1] = returnw; prices[1] = price;
         marketState="Bull Market      ";
      end;
      if(state=4 | state=5 | state=6) then do;
         returnws[2] = returnw; prices[2] = price;
         marketState="Transition Market";
      end;
      if(state=7 | state=8) then do;
         returnws[3] = returnw; prices[3] = price;
         marketState="Bear Market      ";
      end;
   end;
   label marketState="Market State";
run;
proc sgplot data=decodeMaketStates;
   styleattrs datacontrastcolors=(red yellow green);
   scatter x=date y=marketState / group=marketState nomissinggroup;
   yaxis label="Market State" offsetmin=0.1 offsetmax=0.1;
   keylegend / sortorder=reverseauto;
run;

proc sgplot data=decodeMaketStates;
   series x=date y=returnw1 / break lineattrs=(color=green)
                              legendlabel="Bull Market";
   series x=date y=returnw2 / break lineattrs=(color=yellow)
                              legendlabel="Transition Market";
   series x=date y=returnw3 / break lineattrs=(color=red)
                              legendlabel="Bear Market";
   refline 0 / axis=y lineattrs=(color=gray thickness=1);
   yaxis max=20 min=-20 label="Weekly Return";
run;

proc sgplot data=decodeMaketStates;
   series x=date y=price1 / break lineattrs=(color=green)
                            legendlabel="Bull Market";
   series x=date y=price2 / break lineattrs=(color=yellow)
                            legendlabel="Transition Market";
   series x=date y=price3 / break lineattrs=(color=red)
                            legendlabel="Bear Market";
   yaxis label="Market Index";
run;

* print Predictive Performance of VaR Forecasts;
%let alpha1 = 0.80;
%let alpha2 = 0.90;
%let alpha3 = 0.98;
%macro forecastVaR(myLib, dsData, dsModel, timeId, varReturn,
   k, p, iStart, iEnd, oosStart, oosEnd,
   dsForecastPrefix);
   %do i = &iStart. %to &iEnd.;
      proc hmm data=&myLib..&dsData.;
         score inmodel=&myLib..&dsModel.;
         forecast out=&myLib..&dsForecastPrefix.&i. alpha=&&alpha&i. online;
      run;
      data &dsForecastPrefix.&i.;
         set &myLib..&dsForecastPrefix.&i.;
      run;
      proc sort data=&dsForecastPrefix.&i.; by &timeId.; run;
      data &dsForecastPrefix.&i.;
         set &dsForecastPrefix.&i. (FIRSTOBS=&oosStart. OBS=&oosEnd.
            keep=&timeId. &varReturn._Q1
            rename=(&varReturn._Q1=&varReturn._Q1_&k.&p.&i.));
         time=_N_; nState = &k; yLag = &p;
      run;
   %end;
%mend forecastVaR;

%macro evaluateVaR(dsData,dsForecastPrefix,varReturn,k,p,iStart,iEnd,n,dsOut);
   data &dsData.out;
      set &dsData.;
      if (date>&cutDate.) then output;
   run;
   data forecastData;
      set &dsData.out;
      time=_N_;
   run;
   data forecastData;
      merge &dsForecastPrefix.: forecastData;
      by time;
   run;
   data &dsOut.;
      set forecastData;
      retain %do i = &iStart. %to &iEnd.; cviol&k.&p.&i. 0 %end; ;
      retain %do i = &iStart. %to &iEnd.; c&varReturn.&k.&p.&i. 0 %end; ;
      %do i = &iStart. %to &iEnd.;
         if &varReturn. le &varReturn._Q1_&k.&p.&i. then do;
            cviol&k.&p.&i.=cviol&k.&p.&i.+1;
         end;
         c&varReturn.&k.&p.&i. = c&varReturn.&k.&p.&i.
                                  + &varReturn._Q1_&k.&p.&i.;
      %end;
      if _N_ = &n. then do;
         %do i = &iStart. %to &iEnd.;
            Norminal = (1-&&alpha&i.)/2;
            Viol = cviol&k.&p.&i. / &n.;
            LR = 2*(cviol&k.&p.&i.*log(Viol)+(&n.-cviol&k.&p.&i.)*log(1-Viol)
               -(cviol&k.&p.&i.*log(Norminal)+(&n.-cviol&k.&p.&i.)
                *log(1-Norminal)));
            pValue = 1 - cdf("CHISQUARE", LR, 1);
            meanVaR = c&varReturn.&k.&p.&i. / &n.;
            output;
         %end;
      end;
      label Norminal='Target Prob.' Viol='Violation Ratio' LR='LR Stat.'
            pValue='Pr > ChiSq' meanVaR='Avg. of VaR' nState='Number of States'
            yLag='AR Lag';
      keep nState yLag Norminal Viol LR pValue meanVaR;
   run;
   proc print data=&dsOut. noobs label;  format Viol LR pValue meanVaR 6.4; run;
%mend evaluateVaR;

%macro VaR(k,p);
   %forecastVaR(myLib=mylib, dsData=&ds.,
      dsModel=&ds.inMLModel_k&k.To&k._p&p.To&p.,
      timeId=date, varReturn=returnw,
      k=&k., p=&p., iStart=1, iEnd=3,
      oosStart=%eval(3999-(5-&p.)), oosEnd=%eval(4854-1-(5-&p.)),
      dsForecastPrefix=&ds.Forecastk&k._p&p.);
   %evaluateVaR(dsData=&ds.,dsForecastPrefix=&ds.Forecastk&k._p&p.,
      varReturn=returnw,k=&k.,p=&p.,iStart=1,iEnd=3,n=855,
      dsOut=VaR_outputk&k._p&p.);
%mend VaR;

%VaR(8,4);

* print in-sample and out-of-sample average weekly log likelihood;
proc hmm data=mylib.&ds.
         outstat=mylib.&ds.Stat_k2to9_p0to5;
   score inmodel=mylib.&ds.inmlmodel_k2to9_p0to5;
   evaluate out=mylib.&ds.Eval_k2to9_p0to5;
run;
data inLL;
   set mylib.&ds.Eval_k2to9_p0to5(where=(date='26DEC2000'd)
      rename=(logLikelihood=isLL));
run;
data fullLL;
   set mylib.&ds.Eval_k2to9_p0to5(where=(date='22DEC2017'd)
      rename=(logLikelihood=fullLL));
run;
data oosLL;
   merge inLL(in=a) fullLL(in=b);
   by modelIndex;
   oosLL = (fullLL - isLL)/855;
   if(a=b);
   drop date isLL fullLL;
run;

data inLL;
   set mylib.&ds.InMLStat_k2to9_p0to5;
   isLL = logLikelihood / (3999-5);
   keep modelIndex nState yLag isLL;
run;
proc sort data=inLL; by modelIndex; run;

data avgLL;
   merge inLL(in=a) oosLL(in=b);
   by modelIndex;
   if(a=b and (modelIndex=1 or modelIndex=41 or modelIndex=48));
   drop modelIndex;
run;

proc print data=avgLL noobs label;
   label nState='k' yLag='p' isLL='In-Sample' oosLL='Out-Of-Sample';
run;

proc hmm data=mylib.&ds.;
   score inmodel=mylib.&ds.InMLModel_k8To8_p4To4;
   forecast out=mylib.&ds.Forecast_k8to8_p4to4 online;
   display / excludeall;
run;

%let forecastStartDate = '26DEC2000'd;
%let forecastEndDate = '22DEC2017'd;
data ForecastProb;
   set mylib.&ds.Forecast_k8to8_p4to4
          (where=(date>=&forecastStartDate. and date<=&forecastEndDate.));
run;
proc sort data=ForecastProb; by date; run;
data ForecastProb;
   set ForecastProb(keep=state1 - state8) end=last;
   set &ds.(where=(date>=&forecastStartDate. and date<=&forecastEndDate.)
            keep=date price);
   * probability distribution of market states;
   bullProb = state1+state2+state3;
   transitionProb = state4+state5+state6;
   bearProb = state7+state8;
   bearTransitionProb = bearProb + transitionProb;
run;
* plot forecasted prob. dist. of market states and trading strategies;
%macro plotForecastAndTradingStrategies(ds, bearThreshold, bullThreshold,
                                        oneMinusBullThreshold);
   proc sgplot data=&ds.;
      band x=date lower=0 upper=bearProb /
         fillattrs=(color=red) name="bear" legendlabel="Bear Market";
      band x=date lower=bearProb upper=bearTransitionProb /
         fillattrs=(color=yellow)
         name="transition" legendlabel="Transition Market";
      band x=date lower=bearTransitionProb upper=1 /
         fillattrs=(color=green) name="bull" legendlabel="Bull Market";
      series x=date y=price / y2axis lineattrs=(color=blue)
         name="index" legendlabel="Market Index";
      %if %length(&bearThreshold.)>0 %then %do;
         refline &bearThreshold. / axis=y lineattrs=(color=red thickness=3)
            label="Bear Threshold=&bearThreshold." labelloc=inside
            labelattrs=(weight=bold);
      %end;
      %if %length(&bearThreshold.)>0 %then %do;
         refline &oneMinusBullThreshold. / axis=y
            lineattrs=(color=green thickness=3)
            label="Bull Threshold=&bullThreshold." labelloc=inside
            labelattrs=(weight=bold);
      %end;
      xaxis values=("01JAN2001"d to "31DEC2017"d by year) display=(nolabel)
         tickvalueformat=year4. offsetmin=0.03 offsetmax=0.15;
      yaxis label="Forecasted Prob. Dist. of Market States";
      y2axis label="Market Index";
      keylegend "bull" "transition" "bear" "index"/
         location=outside position=bottom;
   run;
%mend;

%plotForecastAndTradingStrategies(ForecastProb);

%let bearThreshold=0.30;
%let bullThreshold=0.10;
%let oneMinusBullThreshold=0.90;
%plotForecastAndTradingStrategies(ForecastProb, &bearThreshold.,
                                  &bullThreshold., &oneMinusBullThreshold.);

%let strategy1=RiskFree;
%let strategy2=Market;
%let strategy3=Bear;
%let strategy4=Bull;
%let strategy1Label=Risk Free;
%let strategy2Label=Market;
%let strategy3Label=Bear;
%let strategy4Label=Bull;
* assume that 1 year = 252 trading days = 50.4 weeks;
* weekly risk free rate = exp(ln(1+APY)/50.4)-1;
%let riskFreeRate = 0.0002661; * APY=1.35%;
%let initialMoney=1;

%macro port_val_cal(stockPosition, stockPrice, riskFreeRate, strategy,
                    initialMoney, rowNo);
   if (&rowNo. = 1) then do;
      &strategy.CashValue = &initialMoney.-&initialMoney.*&stockPosition.;
      &strategy. = &initialMoney.;
      &strategy.StockAmount = &strategy.*&stockPosition./&stockPrice.;
      &strategy.CurrentCashValue = &strategy.CashValue;
  end;
  else do;
      &strategy.CurrentCashValue = &strategy.CashValue*(1+&riskFreeRate.);
      &strategy. = &strategy.CurrentCashValue
                  + &strategy.StockAmount*&stockPrice.;
      &strategy.StockAmount = &strategy.*&stockPosition./&stockPrice.;
      &strategy.CashValue = &strategy.*(1-&stockPosition.);
  end;
  &strategy.ExcessReturn = &strategy./lag(&strategy.)-1-&riskFreeRate.;
  &strategy.ExcessReturn2 = &strategy.ExcessReturn**2;
%mend port_val_cal;

%macro port_val_gen(ForecastProb, bearThreshold, bullThreadhold, riskFreeRate,
                    initialMoney, startDate, endDate);
   data PortValPath;
      set &ForecastProb. end=last;
      array stockPos[4] %do i = 1 %to 4; &&strategy&i..StockPos %end; ;
      retain
         %do i = 1 %to 4;
            &&strategy&i..CashValue 0 &&strategy&i..StockAmount 0
            &&strategy&i..CurrentCashValue 0 &&strategy&i..CumExcessReturn 0
            &&strategy&i..CumExcessReturn2 0
         %end; ;
      * stock position for each strategy;
      stockPos[1] = 0;
      stockPos[2] = 1;
      if(bearProb>=&bearThreshold.) then stockPos[3]=0;
      else stockPos[3]=1;
      if(bullProb>=&bullThreshold.) then stockPos[4]=1;
      else stockPos[4]=0;
      * wealth for each strategy;
      %do i = 1 %to 4;
         %port_val_cal(stockPosition=stockPos[&i.], stockPrice=price,
                       riskFreeRate=&riskFreeRate., strategy=&&strategy&i..,
                       initialMoney=&initialMoney., rowNo=_N_);
         &&strategy&i..CumExcessReturn+&&strategy&i..ExcessReturn;
         &&strategy&i..CumExcessReturn2+&&strategy&i..ExcessReturn2;
      %end;
      * Sharpe ratio and portfolio final wealth;
      if last then do;
         &strategy1.LastPorVal = &strategy1.;
         %do i = 2 %to 4;
            &&strategy&i..LastPorVal = &&strategy&i..;
            &&strategy&i..ExcessReturnMean =
               &&strategy&i..CumExcessReturn/(_N_-1);
            &&strategy&i..ExcessReturnStd =
               sqrt(&&strategy&i..CumExcessReturn2/(_N_-1)
                  - &&strategy&i..ExcessReturnMean**2);
            &&strategy&i..SharpeRatio = round(&&strategy&i..ExcessReturnMean
               /&&strategy&i..ExcessReturnStd, 0.0001);
            call symputx("&&strategy&i..SharpeRatio",
               &&strategy&i..SharpeRatio, 'G');
         %end;
      end;
   run;
%mend;

* plot portfolio values of different trading strategies;
%let color1=black;
%let color2=blue;
%let color3=red;
%let color4=green;
%macro perf_comp_plot_gen(nStrategies);
   proc sgplot data=PortValPath cycleattrs;
      %do i = 1 %to &nStrategies.;
         %let strategy=&&strategy&i..;
         series x=date y=&&strategy&i.. / name="&strategy."
            lineattrs=(color=&&color&i..)
            %if (&i.=1) %then %do;
               legendlabel="&&strategy&i.Label.";
            %end;
            %else %do;
               legendlabel="&&strategy&i.Label. (SR=&&&strategy.SharpeRatio.)";
            %end;
         scatter x=date y=&&strategy&i..LastPorVal / datalabelpos=topright
            markerattrs=(color=&&color&i..)
            datalabel=&&strategy&i..LastPorVal;
      %end;
      xaxis values=("01JAN2001"d to "31DEC2017"d by year) display=(nolabel)
         tickvalueformat=year4. offsetmin=0.03 offsetmax=0.15;
      yaxis max=2.6 min=0.5 label="Wealth Curve";
      keylegend %do i = 1 %to &nStrategies.; "&&strategy&i.." %end; /
         location=outside position=bottom;
   run;
%mend perf_comp_plot_gen;

%port_val_gen(ForecastProb=ForecastProb,
              bearThreshold=&bearThreshold., bullThreadhold=&bullThreshold.,
              riskFreeRate=&riskFreeRate., initialMoney=&initialMoney.,
              startDate=&forecastStartDate., endDate=&forecastEndDate.);
%perf_comp_plot_gen(nStrategies=4);

