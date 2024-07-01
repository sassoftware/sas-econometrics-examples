/*--------------------------------------------------------------

                    SAS Sample Library

        Name: hmmex02.sas
 Description: Example program from SAS Econometrics User's Guide,
              The HMM Procedure
       Title: Clustering Time Series
     Product: SAS Econometrics Software
        Keys: Hidden Markov Model
        PROC: HMM
       Notes:

--------------------------------------------------------------*/

ods graphics on;

title "Clustering Time Series";

* create a data table dgpTable
  consisting of nSections sections and T observations for each section;
* then copy the data table dgpTable to the data table trainingTable
  and the data table libTrainingTable;
%macro createClustersTable(nSections,T,
                           dgpTable,trainingTable,libTrainingTable);
   %let seed = 12345;
   %let nClustersDGP = 2;
   %let c1portion = 0.5;

   %let c1pi1 = 0.5;
   %let c1a11 = 0.6;
   %let c1a22 = 0.6;
   %let c1mu1 = 0;
   %let c1sigma1 = 1;
   %let c1mu2 = 3;
   %let c1sigma2 = 1;

   %let c2pi1 = 0.5;
   %let c2a11 = 0.4;
   %let c2a22 = 0.4;
   %let c2mu1 = 0;
   %let c2sigma1 = 1;
   %let c2mu2 = 3;
   %let c2sigma2 = 1;
   data &dgpTable.;
      do sec = 1 to &nSections.;
         * which cluster;
         * (1) randomly assign clusters;
         * u = uniform(&seed.);
         * if(u<=&c1portion.) then cluster = 1;
         * else cluster = 2;
         * (2) all time series for one cluster are in one block;
         if(sec/&nSections.<=&c1portion.) then cluster = 1;
         else cluster = 2;
         %do i = 1 %to &nClustersDGP.;
            if(cluster=&i.) then do;
               do t = 1 to &T.;
                  if(t=1) then p = &&c&i.pi1;
                  else do;
                     if(lagState=1) then p = &&c&i.a11;
                     else p = 1-&&c&i.a22;
                  end;
                  u = uniform(&seed.);
                  if(u<=p) then state = 1;
                  else state = 2;
                  if(state=1) then
                     y = &&c&i.mu1. + sqrt(&&c&i.sigma1.)*normal(&seed.);
                  else
                     y = &&c&i.mu2. + sqrt(&&c&i.sigma2.)*normal(&seed.);
                  output;
                  lagState = state;
               end;
            end;
         %end;
      end;
   run;

   * trainingTable is sorted by sec and t;
   data &trainingTable.;
      set &dgpTable.;
      keep sec t y;
   run;

   data &libTrainingTable.;
      set &trainingTable.;
   run;
%mend;

%macro estSection(libTrainingTable,casSectionTable,selectedSec,i,T,
                  evalTable,libModelTable);
   data &casSectionTable.&selectedSec.;
      set &libTrainingTable.;
      if(sec=&selectedSec.);
   run;

   proc hmm data=&casSectionTable.&selectedSec.;
      id section=sec time=t;
      model y  / method=ml type=gaussian nstate=&k.;
      optimize ALGORITHM=interiorpoint printlevel=3 printIterFreq=1;
      score outmodel=&libModelTable.&i.;
   run;

   proc hmm data=&libTrainingTable.;
      evaluate out=mylib.tempEval;
      score inmodel=&libModelTable.&i.;
   run;

   data &evalTable.&i.;
      set mylib.tempEval;
      logLikelihood&i. = logLikelihood;
      if(t=&T.) then output;
      keep sec logLikelihood&i.;
   run;

   proc sort data=&evalTable.&i.; by sec; run;
%mend;

%macro assignSeriesToClusters(evalTable,nClusters);
   data &evalTable.;
      merge &evalTable.1 - &evalTable.&nClusters.;
      by sec;
      retain sumLogLikelihood 0;
      array arr[*] logLikelihood1 - logLikelihood&nClusters.;
      sumLogLikelihood = sumLogLikelihood + max(of arr[*]);
      pCluster = whichn(max(of arr[*]), of arr[*]);
      keep sec pCluster sumLogLikelihood;
   run;
%mend;

* initialize clusters;
* each section has T observations;
%macro ctsInitClusters(libTrainingTable,nClusters,nSections,T,
                       evalTable,libModelTable);
   * as a start, use first section or randomly select a section;
   %let selectedSec = 1;
   %let i = 0;
   %estSection(&libTrainingTable.,mylib.tempSection,
               &selectedSec.,&i.,&T.,
               &evalTable.,&libModelTable.);
   proc sort data=&evalTable.&i.; by descending logLikelihood&i.; run;

   %do i = 1 %to &nClusters.;
      data _null_;
         set &evalTable.0;
         pos = floor(1+(&nSections.-1)*(&i.-1)/(&nClusters.-1)+0.5);
         if(_N_=pos) then
            call symput('selectedSec',trim(left(put(sec,12.))));
      run;
      %estSection(&libTrainingTable.,mylib.tempSection,
                  &selectedSec.,&i.,&T.,
                  &evalTable.,&libModelTable.);
   %end;

   %assignSeriesToClusters(&evalTable.,&nClusters.);
%mend;

* estimate HMM for cluster i;
* each section has T observations;
%macro estCluster(libClusteredTable,libTrainingTable,i,T,
                  evalTable,libModelTable);
   data &libClusteredTable.&i.;
      set &libClusteredTable.;
      if(pCluster=&i.);
   run;

   proc hmm data=&libClusteredTable.&i.;
      id section=sec time=t;
      model y  / method=ml type=gaussian nstate=&k.;
      optimize ALGORITHM=interiorpoint printlevel=3 printIterFreq=1;
      score outmodel=&libModelTable.&i.;
   run;

   proc hmm data=&libTrainingTable.;
      evaluate out=mylib.tempEval;
      score inmodel=&libModelTable.&i.;
   run;

   data &evalTable.&i.;
      set mylib.tempEval;
      logLikelihood&i. = logLikelihood;
      if(t=&T.) then output;
      keep sec logLikelihood&i.;
   run;

   proc sort data=&evalTable.&i.; by sec; run;
%mend;

* moves time series between clusters based on their likelihoods;
* each section has T observations;
* (1) &libModelTable.&i. has the model estimates for cluster i,
* and (2) &evalTable. has the clustering result;
%macro ctsClustering(trainingTable,libTrainingTable,nClusters,T,
                     libClusteredTable,evalTable,libModelTable);
   * save last clustering result;
   data prev&evalTable.; set &evalTable.; run;

   * clustering data;
   data &libClusteredTable.;
      merge &trainingTable.(in=a) &evalTable.(in=b);
      by sec;
   run;

   %do i = 1 %to &nClusters.;
      %estCluster(&libClusteredTable.,&libTrainingTable.,&i.,&T.,
                  &evalTable.,&libModelTable.);
   %end;

   %assignSeriesToClusters(&evalTable.,&nClusters.);

   * compare current and previous clustering results;
   proc compare data=&evalTable. compare=prev&evalTable.; run;
   %let converged=%eval(&sysinfo.=0);
   %let nIter=%eval(&nIter.+1);
   %let done=%eval(&converged. or &nIter.>=&maxIter.);

   %put nIter     = &nIter.;
   %put converged = &converged.;
   %put done      = &done.;
%mend;

* cluster time series;
* each section has T observations;
%macro ctsLoopClustering(trainingTable,libTrainingTable,nClusters,T,
                         libClusteredTable,evalTable,libModelTable);
   %let done = 0;     * terminate clustering;
   %do %while(&done.=0);
         %ctsClustering(&trainingTable.,&libTrainingTable.,&nClusters.,&T.,
                     &libClusteredTable.,&evalTable.,&libModelTable.);
   %end;
%mend;

* compare the clustering results with the true clusters in the DGP;
* it is applicable only when the DGP is available;
%macro evaluateClusteringResults(dgpTable,evalTable,nSections,T);
   data true&evalTable.;
      set &dgpTable.;
      if(t=&T.);
      keep sec cluster;
   run;

   data ctsCheck;
      merge true&evalTable.(in=a) &evalTable.(in=b);
      by sec;
      if(cluster=pCluster) then do; correct=1; end;
      else do; correct=0; end;
      if(a=b);
   run;

   data ctsAccuracy;
      set ctsCheck;
      retain count 0 correctCount 0;
      count = count + 1;
      correctCount = correctCount + correct;
      if(count>&nSections.-0.5) then do;
         accuracy = correctCount / count;
         if(accuracy<0.5) then accuracy = 1 - accuracy;
         nIterations = &nIter.;
         converged = &converged.;
         output;
      end;
      keep nIterations converged accuracy;
   run;

   proc print data = ctsAccuracy noobs;
      var nIterations converged accuracy;
      format accuracy 6.4;
   run;
%mend;

%let N = 10000;
%let T = 200;
%createClustersTable(&N.,&T.,ctsDGP,ctsTrain,mylib.ctsTrain);
%let nClusters = 2; * number of clusters;
%let k=2;           * number of states in the HMM for each cluster;
%ctsInitClusters(mylib.ctsTrain,&nClusters.,&N.,&T.,
    evalCluster,mylib.modelCluster);
%let maxIter = 100;* maximum number of iterations;
%let nIter = 0;    * number of iterations;
%let converged = 0;* converge status;
%ctsLoopClustering(ctsTrain,mylib.ctsTrain,&nClusters.,&T.,
    mylib.ctsTrainClustered,evalCluster,mylib.modelCluster);
%evaluateClusteringResults(ctsDGP,evalCluster,&N.,&T.);

* print parameter estimates for cluster 1;
proc hmm data=mylib.ctsTrain;
   score inmodel=mylib.modelCluster1;
run;

* print parameter estimates for cluster 2;
proc hmm data=mylib.ctsTrain;
   score inmodel=mylib.modelCluster2;
run;

%let N = 10000;
%let T = 400;
%createClustersTable(&N.,&T.,ctsDGP,ctsTrain,mylib.ctsTrain);
%let nClusters = 2; * number of clusters;
%let k=2;           * number of states in the HMM for each cluster;
%ctsInitClusters(mylib.ctsTrain,&nClusters.,&N.,&T.,
    evalCluster,mylib.modelCluster);
%let maxIter = 100;* maximum number of iterations;
%let nIter = 0;    * number of iterations;
%let converged = 0;* converge status;
%ctsLoopClustering(ctsTrain,mylib.ctsTrain,&nClusters.,&T.,
    mylib.ctsTrainClustered,evalCluster,mylib.modelCluster);
%evaluateClusteringResults(ctsDGP,evalCluster,&N.,&T.);

* print parameter estimates for cluster 1;
proc hmm data=mylib.ctsTrain;
   score inmodel=mylib.modelCluster1;
run;

* print parameter estimates for cluster 2;
proc hmm data=mylib.ctsTrain;
   score inmodel=mylib.modelCluster2;
run;


%let N = 10000;
%let T = 20;
%createClustersTable(&N.,&T.,ctsDGP,ctsTrain,mylib.ctsTrain);
%let nClusters = 2; * number of clusters;
%let k=2;           * number of states in the HMM for each cluster;
%ctsInitClusters(mylib.ctsTrain,&nClusters.,&N.,&T.,
    evalCluster,mylib.modelCluster);
%let maxIter = 100;* maximum number of iterations;
%let nIter = 0;    * number of iterations;
%let converged = 0;* converge status;
%ctsLoopClustering(ctsTrain,mylib.ctsTrain,&nClusters.,&T.,
    mylib.ctsTrainClustered,evalCluster,mylib.modelCluster);
%evaluateClusteringResults(ctsDGP,evalCluster,&N.,&T.);

* print parameter estimates for cluster 1;
proc hmm data=mylib.ctsTrain;
   score inmodel=mylib.modelCluster1;
run;

* print parameter estimates for cluster 2;
proc hmm data=mylib.ctsTrain;
   score inmodel=mylib.modelCluster2;
run;

