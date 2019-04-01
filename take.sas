data melanoma_file;
set "D:\KUMC COURSE\spring 2019\survival analysis\Exam\melanoma.sas7bdat";
run;


/* ===============================================================
                      Ques No 1
==================================================================*/

/*Part 1(a)*/

PROC LIFETEST data = melanoma_file;
TIME SURV*DEAD(0);
STRATA TRT;
RUN;

/*Part 1(b)*/

PROC LIFETEST data = melanoma_file;
TIME REMDUR*RECUR(0);
STRATA TRT;
RUN;

/*Part 1(c)*/



PROC LIFETEST data = melanoma_file;
TIME SURV*DEAD(0);
*STRATA TRT;
TEST TRT;
RUN;

PROC LIFETEST data = melanoma_file;
TIME SURV*DEAD(0);
*STRATA TRT;
TEST SEX;
RUN;

PROC LIFETEST data = melanoma_file;
TIME SURV*DEAD(0);
*STRATA TRT;
TEST AGE;
RUN;

PROC LIFETEST data = melanoma_file;
TIME SURV*DEAD(0);
*STRATA TRT;
TEST TRT SEX AGE;
RUN;

/*Part 1(d)*/


PROC LIFETEST data = melanoma_file PLOTS = (LS LLS h);
TIME SURV*DEAD(0);
*STRATA TRT;
RUN;

/* ===============================================================
                      Ques No 2
==================================================================*/

/*Ques No 2 (a)*/

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = EXPONENTIAL;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = weibull;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = lognormal;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = logistic;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = llogistic;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = gamma;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE/ DIST = normal;
RUN;
title "NLMIXED: Gompertz distribution";
proc nlmixed data=melanoma_file;
parms log_gamma -20 to 0 by 5;
gamma = exp(log_gamma);
linp = b0 + b1*trt + b2*sex +b3*age;
alpha = exp(-linp);
G_t = exp((alpha/gamma)*(1 - exp(gamma*surv)));
g = alpha*exp(gamma*surv)*G_t;
ll = (dead = 0)*log(g) + /* ll for observed failures */
(dead = 1)*log(G_t); /* ll for censored failures */
model ll ~ general(ll);
estimate "gamma" exp(log_gamma);
run;
/**/
/*Nested model checking for log normal and generalizd gamma*/

data lltest;
test = 2*( -19.39491312 + 21.0152007);
p_value = 1-PROBCHI(test,1);
run;

proc print;run;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE / DIST = lognormal;
RUN;

PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE TRT*SEX/ DIST = lognormal;
RUN;
PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE TRT*AGE/ DIST = lognormal;
RUN;
PROC LIFEREG DATA = melanoma_file;
CLASS TRT SEX;
MODEL SURV*DEAD(0) = TRT SEX AGE TRT*SEX*AGE/ DIST = lognormal;
RUN;

/*Ques No 2 (b)*/
TITLE;
DATA ONE;
C=1; AGE = 40; TRT = 0; OUTPUT;
C=1; AGE = 40; TRT = 1; OUTPUT;
RUN;

DATA MELANOMA_ADD;
SET melanoma_file ONE;
RUN;

PROC LIFEREG DATA = MELANOMA_ADD;
CLASS TRT ;
MODEL SURV*DEAD(0) = TRT / DIST = LNORMAL;
OUTPUT OUT = SURV_EST QUANTILES = 0.2 TO 0.98 BY 0.02 PREDICTED = PRED CONTROL = C;
RUN;

DATA SURV_EST_2;
SET SURV_EST;
SDF = 100*(1-_PROB_);
RUN;


PROC GPLOT DATA = SURV_EST_2;
PLOT SDF*PRED = TRT ;
SYMBOL1 I = SPLINE COLOR = RED L=1;
SYMBOL2 I = SPLINE COLOR = BLUE L=1;
RUN;




/* ===============================================================
                      Ques No 3
==================================================================*/

/* Mostly calculation   */



/* ===============================================================
                      Ques No 4
==================================================================*/

/*Ques No 4 (a)*/
data addict;
set "D:\KUMC COURSE\spring 2019\survival analysis\Exam\addicts_1.sas7bdat";
run;

/*Ques No 4 (c)*/

/*Only including dose*/
PROC PHREG DATA = addict;
CLASS CLINIC;
MODEL LENGTH*STATUS(0)  = CLINIC DOSE/RISKLIMITS;
RUN;

/*Only including prison*/

PROC PHREG DATA = addict;
CLASS CLINIC PRISON;
MODEL LENGTH*STATUS(0)  = CLINIC PRISON/RISKLIMITS;
RUN;

PROC PHREG DATA = addict;
CLASS CLINIC PRISON;
MODEL LENGTH*STATUS(0)  = CLINIC PRISON DOSE/RISKLIMITS;
RUN;

/*Ques No 4 (d)*/

data addict_mod;
set addict;
if dose < 60 then dose_cat = 1;
Else if 60 <= dose and dose < 80 then dose_cat = 2;
Else if dose >= 80 then dose_cat = 3;
run; 

PROC PHREG DATA = addict_mod;
CLASS CLINIC PRISON dose_cat;
MODEL LENGTH*STATUS(0) = CLINIC PRISON dose_cat/RISKLIMITS;
RUN;



PROC PHREG DATA = addict_mod;
CLASS CLINIC PRISON dose_cat;
MODEL LENGTH*STATUS(0) = CLINIC PRISON dose_cat/RISKLIMITS;
output out=Outp xbeta=Xb resmart=Mart resdev=Dev;
RUN;


proc sgplot data=Outp;
      yaxis grid;
      refline 0 / axis=y;
      scatter y=Mart x=Xb;
      run;
   proc sgplot data=Outp;
      yaxis grid;
      refline 0 / axis=y;
      scatter y=Dev x=Xb;
      run;
