/*--------------------------------------------------------------

                    SAS Sample Library

        Name: hmmex03.sas
 Description: Example program from SAS/ETS User's Guide,
              The HMM Procedure
       Title: Analysis of the Business Cycle
     Product: SAS/ETS Software
        Keys: Hidden Markov Model
        PROC: HMM
       Notes:

--------------------------------------------------------------*/

ods graphics on;

title 'Analysis of the Business Cycle';
data gnpHamilton;
   input date: anydtdte. gnp;
   format date yyq6.;
   dgnp = 100*(log(gnp)-log(lag(gnp)));
   if cmiss(of _all_) then delete;
   label date="Date" gnp="GNP" dgnp="Growth Rate of GNP";
datalines;
1951q1  1286.6
1951q2  1320.4
1951q3  1349.8
1951q4  1356.0
1952q1  1369.2
1952q2  1365.9
1952q3  1378.2
1952q4  1406.8
1953q1  1431.4
1953q2  1444.9
1953q3  1438.2
1953q4  1426.6
1954q1  1406.8
1954q2  1401.2
1954q3  1418.0
1954q4  1438.3
1955q1  1469.6
1955q2  1485.7
1955q3  1505.5
1955q4  1518.7
1956q1  1515.7
1956q2  1522.6
1956q3  1523.7
1956q4  1540.6
1957q1  1553.3
1957q2  1552.4
1957q3  1561.5
1957q4  1537.3
1958q1  1506.1
1958q2  1514.2
1958q3  1550.0
1958q4  1568.7
1959q1  1606.4
1959q2  1634.0
1959q3  1629.5
1959q4  1634.4
1960q1  1671.6
1960q2  1666.8
1960q3  1668.4
1960q4  1654.1
1961q1  1671.3
1961q2  1692.1
1961q3  1716.3
1961q4  1754.9
1962q1  1777.9
1962q2  1796.4
1962q3  1813.1
1962q4  1810.1
1963q1  1834.6
1963q2  1860.0
1963q3  1892.5
1963q4  1906.1
1964q1  1948.7
1964q2  1965.4
1964q3  1985.2
1964q4  1993.7
1965q1  2036.9
1965q2  2066.4
1965q3  2099.3
1965q4  2147.6
1966q1  2190.1
1966q2  2195.8
1966q3  2218.3
1966q4  2229.2
1967q1  2241.8
1967q2  2255.2
1967q3  2287.7
1967q4  2300.6
1968q1  2327.3
1968q2  2366.9
1968q3  2385.3
1968q4  2383.0
1969q1  2416.5
1969q2  2419.8
1969q3  2433.2
1969q4  2423.5
1970q1  2408.6
1970q2  2406.5
1970q3  2435.8
1970q4  2413.8
1971q1  2478.6
1971q2  2478.4
1971q3  2491.1
1971q4  2491.0
1972q1  2545.6
1972q2  2595.1
1972q3  2622.1
1972q4  2671.3
1973q1  2734.0
1973q2  2741.0
1973q3  2738.3
1973q4  2762.8
1974q1  2747.4
1974q2  2755.2
1974q3  2719.3
1974q4  2695.4
1975q1  2642.7
1975q2  2669.6
1975q3  2714.9
1975q4  2752.7
1976q1  2804.4
1976q2  2816.9
1976q3  2828.6
1976q4  2856.8
1977q1  2896.0
1977q2  2942.7
1977q3  3001.8
1977q4  2994.1
1978q1  3020.5
1978q2  3115.9
1978q3  3142.6
1978q4  3181.6
1979q1  3181.7
1979q2  3178.7
1979q3  3207.4
1979q4  3201.3
1980q1  3233.4
1980q2  3157.0
1980q3  3159.1
1980q4  3199.2
1981q1  3261.1
1981q2  3250.2
1981q3  3264.6
1981q4  3219.0
1982q1  3170.4
1982q2  3179.9
1982q3  3154.5
1982q4  3159.3
1983q1  3190.6
1983q2  3259.3
1983q3  3303.4
1983q4  3357.2
1984q1  3449.4
1984q2  3492.6
1984q3  3510.4
1984q4  3515.6
;

proc sgplot data=gnpHamilton;
   series x=date y=dgnp;
   series x=date y=gnp / y2axis;
   refline 0 / axis=y lineattrs=(color=gray thickness=1);
run;

data mylib.gnpHamilton; set gnpHamilton; run;
proc hmm data=mylib.gnpHamilton labelswitch=(sort=desc(const));
   id time=date;
   model dgnp / armean=adjusted type=ar ylag=4 method=ml nstate=2
                stateIndependent=(ar cov);
   optimize maxiter=0;
   initial tpm={0.9049 0.0951, 0.2450 0.7550},
      ar={0.014 -0.058 -0.247 -0.213, 0.014 -0.058 -0.247 -0.213},
      const={1.1643, -0.3577},
      cov={0.5914, 0.5914};
   filter out=mylib.filter;
   decode out=mylib.decode;
   smooth out=mylib.smooth;
run;

data filter; set mylib.filter; run;
proc sort data=filter; by date; run;
proc sgplot data=filter;
   series x=date y=state2;
   yaxis label="Filtered Prob. of Recession";
run;

data plot2state;
   merge mylib.decode mylib.gnpHamilton;
   if state=2 then do; recession=1; boom=0; end;
   else if state=1 then do; recession=0; boom=1; end;
   rec = recession*(-3);
   boo = boom*3;
   by date;
   if cmiss(of _all_) then delete;
run;

proc sgplot data=plot2state;
   band x=date upper=0 lower=rec / legendlabel='Recession' type=step
                                   fillattrs=(color=red);
   band x=date upper=boo lower=0 / legendlabel='Boom' type=step
                                   fillattrs=(color=green);
   series x=date y=dgnp / lineattrs=(thickness=2);
run;

data mylib.nber;
   input date: anydtdte. recessionNBER;
   format date yyq6.;
   label date="Date" recessionNBER="Recession Reported from NBER";
datalines;
1951q1  0
1951q2  0
1951q3  0
1951q4  0
1952q1  0
1952q2  0
1952q3  0
1952q4  0
1953q1  0
1953q2  0
1953q3  1
1953q4  1
1954q1  1
1954q2  1
1954q3  0
1954q4  0
1955q1  0
1955q2  0
1955q3  0
1955q4  0
1956q1  0
1956q2  0
1956q3  0
1956q4  0
1957q1  0
1957q2  0
1957q3  1
1957q4  1
1958q1  1
1958q2  1
1958q3  0
1958q4  0
1959q1  0
1959q2  0
1959q3  0
1959q4  0
1960q1  0
1960q2  1
1960q3  1
1960q4  1
1961q1  1
1961q2  0
1961q3  0
1961q4  0
1962q1  0
1962q2  0
1962q3  0
1962q4  0
1963q1  0
1963q2  0
1963q3  0
1963q4  0
1964q1  0
1964q2  0
1964q3  0
1964q4  0
1965q1  0
1965q2  0
1965q3  0
1965q4  0
1966q1  0
1966q2  0
1966q3  0
1966q4  0
1967q1  0
1967q2  0
1967q3  0
1967q4  0
1968q1  0
1968q2  0
1968q3  0
1968q4  0
1969q1  0
1969q2  0
1969q3  0
1969q4  1
1970q1  1
1970q2  1
1970q3  1
1970q4  1
1971q1  0
1971q2  0
1971q3  0
1971q4  0
1972q1  0
1972q2  0
1972q3  0
1972q4  0
1973q1  0
1973q2  0
1973q3  0
1973q4  1
1974q1  1
1974q2  1
1974q3  1
1974q4  1
1975q1  1
1975q2  0
1975q3  0
1975q4  0
1976q1  0
1976q2  0
1976q3  0
1976q4  0
1977q1  0
1977q2  0
1977q3  0
1977q4  0
1978q1  0
1978q2  0
1978q3  0
1978q4  0
1979q1  0
1979q2  0
1979q3  0
1979q4  0
1980q1  1
1980q2  1
1980q3  1
1980q4  0
1981q1  0
1981q2  0
1981q3  1
1981q4  1
1982q1  1
1982q2  1
1982q3  1
1982q4  1
1983q1  0
1983q2  0
1983q3  0
1983q4  0
1984q1  0
1984q2  0
1984q3  0
1984q4  0
;

data gnpsmooth;
   merge mylib.smooth mylib.nber mylib.gnpHamilton;
   if state2>=0.5 then recession=1; else recession=0;
   by date;
   gnp2=gnp;
   if cmiss(of _all_) then delete;
   label recession="Recession Found by the Model" gnp2="GNP";
run;
proc transpose data=gnpsmooth  out=gnpsmooth_t;
   by date;
   var recessionNBER recession;
run;
proc transpose data=gnpsmooth  out=gnpsmooth_t2;
   by date;
   var gnp gnp2;
run;
data gnpplot;
   merge gnpsmooth_t(rename=(_name_=name1 _label_=label1 col1=recession))
         gnpsmooth_t2(rename=(_name_=name2 _label_=label2 col1=gnp));
   by date;
   rec=recession*3520;
run;
proc sgpanel data=gnpplot;
   panelby label1 / layout=rowlattice uniscale=column novarname;
   band x=date upper=rec lower=0 / legendlabel='Recession' type=step;
   series x=date y=gnp / legendlabel='GNP';
   rowaxis display=(nolabel) values=(1000 to 3520 by 500);
run;

proc hmm data=mylib.gnpHamilton labelswitch=(sort=desc(const));
   id time=date;
   model dgnp / armean=adjusted type=ar ylag=4 method=ml nstate=2
             stateIndependent=(ar cov);
   optimize  printIterFreq=1 printlevel=3 ;
   filter out=mylib.filter;
   decode out=mylib.decode;
   smooth out=mylib.smooth;
run;


data filter; set mylib.filter; run;
proc sort data=filter; by date; run;
proc sgplot data=filter;
   series x=date y=state2;
   yaxis label="Filtered Prob. of Recession";
run;

data plot2state;
   merge mylib.decode mylib.gnpHamilton;
   if state=2 then do; recession=1; boom=0; end;
   else if state=1 then do; recession=0; boom=1; end;
   rec = recession*(-3);
   boo = boom*3;
   by date;
   if cmiss(of _all_) then delete;
run;

proc sgplot data=plot2state;
   band x=date upper=0 lower=rec / legendlabel='Recession' type=step
                                   fillattrs=(color=red);
   band x=date upper=boo lower=0 / legendlabel='Boom' type=step
                                   fillattrs=(color=green);
   series x=date y=dgnp / lineattrs=(thickness=2);
run;

data gnpsmooth;
   merge mylib.smooth mylib.nber mylib.gnpHamilton;
   if state2>=0.5 then recession=1; else recession=0;
   by date;
   gnp2=gnp;
   if cmiss(of _all_) then delete;
   label recession="Recession Found by the Model" gnp2="GNP";
run;
proc transpose data=gnpsmooth  out=gnpsmooth_t;
   by date;
   var recessionNBER recession;
run;
proc transpose data=gnpsmooth  out=gnpsmooth_t2;
   by date;
   var gnp gnp2;
run;
data gnpplot;
   merge gnpsmooth_t(rename=(_name_=name1 _label_=label1 col1=recession))
         gnpsmooth_t2(rename=(_name_=name2 _label_=label2 col1=gnp));
   by date;
   rec=recession*3520;
run;
proc sgpanel data=gnpplot;
   panelby label1 / layout=rowlattice uniscale=column novarname;
   band x=date upper=rec lower=0 / legendlabel='Recession' type=step;
   series x=date y=gnp / legendlabel='GNP';
   rowaxis display=(nolabel) values=(1000 to 3520 by 500);
run;

