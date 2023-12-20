/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
* 
* Project: HTE Analysis - limiting cohort to severe sepsis patients only
*
* Author: Sarah Seelye (modified Jennifer Cano's code)
*
* Description: Running robust Poisson models/creating figures
* 
/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/;


libname h 'filepath';
libname new 'filepath';
libname rr 'filepath';

*import data set;
proc import datafile="filepath"
out=h.happi
dbms=dta replace;
run;

%let Tx=Tx0;

data sepsis2;
set h.happi(drop=aod_ind where=(sum(SEPTIC_SHOCK,SEVERE_SEPSIS) ge 1));
AOD=sum(AOD_LUNG,AOD_KIDNEY,AOD_LIVER,AOD_HEME,AOD_LACTATE,pressor_in_72hr);
if aod ge 3 then aod=3;
survival=1-mort30_ed;
cancer=0;
if cancer_nonmet eq 1 then cancer=1;
if cancer_met eq 1 then cancer=2;
Tx=0;
if .<TIME_TO_ABX_min<150 then Tx=1;
*alternative treatment;
if .<TIME_TO_ABX_min<=180 then Tx0=1; *including 2-3 time group;
if TIME_TO_ABX_min>180 then Tx0=0;
if .<age <40 then ageg=1;
else if 40<=age<50 then ageg=3;
else if 50<=age<60 then ageg=4;
else if 60<=age<70 then ageg=5;
else if 70<=age<80 then ageg=6;
else if 80<=age<90 then ageg=7;
else if age ge 90 then ageg=8; 
if aod ge 1 then aod_ind = 1;
	else aod_ind =0;
run;

*/*/**/*/*/ Create New Primary Cohort */*/*/*/;

*calculate 12 hours by using the time_to_abx_min var
*remove facilities with <15 sepsis hospitalizations
*restrict to severe sepsis only (exclude shock);
data rr.happi_hte_noshock_cohort;
set sepsis2;
if time_to_abx_min <=720 then abx_in_12hr =1;
else abx_in_12hr =0;
if (hospid not in ("523A5", "542", "598A0", "657A0")) and (abx_in_12hr =1) and (severe_sepsis = 1);
run;

proc sql;
select count(distinct uniqid)
from rr.happi_hte_noshock_cohort;
quit;

proc sql;
select count(distinct uniqid), data
from rr.happi_hte_noshock_cohort
group by data;
quit;

/*/*/*/*Create Tables */*/*/*/;

proc means data=rr.happi_hte_noshock_cohort n median p25 p75 maxdec=1;
class tx0;
var TIME_TO_ABX_HR age;
ways 0 1;
run;

*p values;
proc npar1way data = rr.happi_hte_noshock_cohort wilcoxon;
class tx0;
var TIME_TO_ABX_HR age ;
output out=ps wilcoxon;
run;

*cat vars;
proc freq data=rr.happi_hte_noshock_cohort;
tables tx0  / nocum;
run;
proc tabulate data=rr.happi_hte_noshock_cohort /*missing order=freq*/
format =8.1;
class  tx0 sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
	coag cancer_met cancer_nonmet;

table sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
	coag cancer_met cancer_nonmet
all ='Column Total',
all(colpctn) tx0*(n colpctn) / nocellmerge;
run;

proc freq data=rr.happi_hte_noshock_cohort;
table cancer*CANCER_MET*CANCER_NONMET/list;
table Tx0 Tx*Tx0/missing;
table aod ageg sirs_temp;
table AOD_LUNG*AOD_KIDNEY*AOD_LIVER*AOD_HEME*AOD_LACTATE*pressor_in_72hr/list;
run;

*p values;

proc freq data=rr.happi_hte_noshock_cohort;
tables tx0*(male sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
		coag cancer_met cancer_nonmet)
	/ chisq exact;
run;

proc freq data=rr.happi_hte_noshock_cohort;
table aod liver coag sirs_temp;
run;

proc sql;
select count(uniqid)
from rr.happi_hte_noshock_cohort
where CANCER_MET = 1 and CANCER_NONMET = 1;
quit;


*7/21/23 lump 2 anemia vars to single condition;
proc print data= rr.happi_hte_noshock_cohort (obs=100);
var cancer:;
run;


proc freq data=rr.happi_hte_noshock_cohort;
tables cancer cancer*tx0 cancer_met*cancer_nonmet;
run;


data cancer;
set rr.happi_hte_noshock_cohort;
if cancer_met = 1 then cancer_nonmet = 0;
	else cancer_nonmet = cancer_nonmet;
	run;

proc freq data=cancer;
tables cancer cancer_nonmet cancer*tx0 cancer_met*cancer_nonmet;
run;
proc freq data=cancer ;
tables tx0*cancer_nonmet / chisq exact;
run;

data new;
set rr.happi_hte_noshock_cohort;
if anemia_cbl = 1 or anemia_def = 1 then anemia_new =  1;
	else anemia_new = 0;
run;

proc print data= new (obs=100);
var anemia:;
run;

********************************************************************************************************************;
*Set up data sets for Poisson models;

*create output KP + VA data set used in Poisson models;
proc logistic data=new;
class cancer(ref='0');
model mort30_ed(event='1')= SEPTIC_SHOCK SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr HTN
MALE AGE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA CANCER /*CANCER*/ RA COAG OBESITY WTLOSS FEN anemia_new
output out=rr.sepsis_risk_VA_KP_RR_noshock pred=risk;
run;

*KP + VA data set for new results of causal forest in R;
data rr.sepsis_risk_VA_KP_R_RR_noshock;
set rr.sepsis_risk_VA_KP_RR_noshock;
if cancer_met = 1 then cancer_nonmet = 0;
	else cancer_nonmet = cancer_nonmet;
keep Tx Tx0 AOD AGEG
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr /*AOD_IND*/ HTN
MALE AGE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
run;


/*Model with no Interaction*/

proc glimmix data=rr.sepsis_risk_VA_KP_RR_noshock;
class Tx0(ref='0') cancer uniqid;
model mort30_ed(event='1')=Tx0 
SIRS_PULSE SIRS_RR SIRS_WBC SIRS_TEMP AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr /*AOD_IND*/ MALE AGE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA RA COAG CANCER OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//dist=poisson link=log;
lsmeans Tx0/pdiff cl;
ods output lsmeans=lsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=diffs;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

proc print data=lsm;
run;
proc print data=diffs;
run;

*LSM;
data lsm_v2;
set lsm;
where Tx0 = 0;
	BL_RR = exp(Estimate);
	BL_LCL = exp(Lower);
	BL_UCL = exp(Upper);
keep BL_RR BL_LCL BL_UCL;
run;

*DIFFSA;
data diffsa;
set diffs;
keep RR RR_LCL RR_UCL;
RR = exp(Estimate);
RR_LCL = exp(Lower);
RR_UCL = exp(Upper);
run;

*Merge;
data merge_RR;
merge lsm_v2 diffsa;
run;


*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*these are the parameters;
/*%macro RDCI(BR, BR_LOW, BR_HI, RR, RR_LOW, RR_HI);*/

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set merge_RR;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(tx0, 1, 0.1081065628,0.1061404559,0.1101090891,0.9250450558,0.9048080174,0.9457347181);

data tx0;
set tx0;
drop level;
run;

*merge all tables;
data merge_rr_ad;
merge merge_RR tx0;
run;

 *AR and RR columns;
data tab2;
set merge_rr_ad;
length AMR $ 25;
/*AMR = compress((abs(round(RD*100, .01))||'('||abs(round(RDU*100, .01))||','||abs(round(RDL*100, .01))||')'));*/
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep AMR RelRis;
run;


/*/*/*/*/* Sum of AODs */*/*/*/*/;
%let Tx=Tx0;

proc glimmix data=rr.sepsis_risk_VA_KP_RR_noshock;
class AOD(ref='1') &Tx. uniqid cancer;
model mort30_ed(event='1')=&Tx. AOD*&Tx. AOD
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
MALE AGE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA cancer RA COAG OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//dist=poisson link=log;
lsmeans AOD*&Tx./cl ilink diff ;
/*slice AOD*&Tx./sliceby=AOD exp means cl;*/
lsmeans AOD/cl ilink;
lsmeans &Tx./pdiff cl ilink;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
ods output lsmeans=sumofaodlsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=diffs;
run;

*Create table for plotting;
*LSM;
proc sort data=sumofaodlsm ;
by effect aod &tx. ;
run;
data sumofaodlsm_v2;
set sumofaodlsm;
where Tx0 = 0 and effect ne 'Tx0';
keep AOD Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL
	   AOD = Level;
run;

*DIFFSA;
data diffsa;
set diffs;
if aod eq _aod;
run;
proc sort data=diffsa;
by aod ;
run;
data sumofaoddiffsa;
set diffsa;
where effect ne 'Tx0';
keep AOD RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
rename AOD = Level;
run;

*Merge;
proc sql;
create table sumaod_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from sumofaodlsm_v2 a
left join sumofaoddiffsa b
	on a.level = b.level;
quit;

*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*these are the parameters;
/*%macro RDCI(BR, BR_LOW, BR_HI, RR, RR_LOW, RR_HI);*/

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set sumaod_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(aod1, 1, 0.0946236373,0.0925890574,0.0967029257,0.9382652147,0.9119437311,0.9653464168);
%RDCI(aod2, 2, 0.1408653337,0.1367831199,0.1450693787,0.9312988334,0.8959315239,0.9680622836);
%RDCI(aod3, 3, 0.2275115537,0.2174903483,0.2379945018,0.8092616385,0.7585719617,0.8633385264);

*merge all tables;
data sumaod_plot_rr_ad;
merge sumaod_plot_rr aod1 aod2 aod3;
by level;
run;

proc freq data=rr.sepsis_risk_VA_KP_RR_noshock;
tables aod;
run;

 *AR and RR columns;
data tab2;
set sumaod_plot_rr_ad;
length AMR $ 25;
/*AMR = compress((abs(round(RD*100, .01))||'('||abs(round(RDU*100, .01))||','||abs(round(RDL*100, .01))||')'));*/
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep AMR RelRis;
run;



/*/*/*/*/*/*/* Cancer */*/*/*/*/*/*/;

options mprint symbolgen;
*Generic macro;
%let Tx=Tx0;
%let cov=cancer;	
%let adj=sirs_temp;
%let ref=(ref='0');
proc glimmix data=rr.sepsis_risk_VA_KP_RR_noshock;
class &Tx. &cov.&ref. uniqid;
model mort30_ed(event='1')=&Tx. &cov.*&Tx. &cov. AGE &adj.
SIRS_PULSE SIRS_RR SIRS_WBC AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr /*AOD_IND*/ MALE /*AGE*/ HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA RA COAG OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//dist=poisson link=log;
lsmeans &cov.*&Tx./cl ilink diff ;
slice &cov.*&Tx./sliceby=&cov. exp means cl;
lsmeans &cov./cl ilink diff;
lsmeans &Tx./pdiff cl;
ods output lsmeans=cancerlsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=cancerdiffs;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

proc print data=cancerlsm;
run;
proc print data= cancerdiffs;
run;

*Create table for plotting;
*LSM;
proc sort data=cancerlsm ;
by effect cancer &tx. ;
run;
data cancerlsm_v2;
set cancerlsm;
where Tx0 = 0 and effect ne 'Tx0';
keep cancer Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL
	   cancer = Level;
run;

*DIFFSA;
data cancerdiffsa;
set cancerdiffs;
if cancer eq _cancer;
run;
proc sort data=cancerdiffsa;
by cancer ;
run;
data cancerdiffsa;
set cancerdiffsa;
where effect ne 'Tx0';
keep cancer RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
rename cancer = Level;
run;

*Merge;
proc sql;
create table cancer_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from cancerlsm_v2 a
left join cancerdiffsa b
	on a.level = b.level;
quit;

*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*these are the parameters;
/*%macro RDCI(BR, BR_LOW, BR_HI, RR, RR_LOW, RR_HI);*/

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set cancer_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(cancer0, 0, 0.0717119302,0.0702668811, 0.073186697,0.9789272237,0.9525439138,1.0060412917);
%RDCI(cancer1, 1, 0.091063752,0.0877561424,0.0944960284,0.8885405314,0.8406823121,0.9391232153);
%RDCI(cancer2, 2, 0.2042467804,0.1976865369,0.2110247261,0.8007865151,0.7637967453,0.8395676556);

*merge all tables;
data cancer_plot_rr_ad;
merge cancer_plot_rr cancer0 cancer1 cancer2;
by level;
run;

proc freq data=rr.sepsis_risk_VA_KP_RR_noshock;
tables cancer;
run;

 *AR and RR columns;
data tab2;
set cancer_plot_rr_ad;
length AMR $ 25;
/*AMR = compress((abs(round(RD*100, .01))||'('||abs(round(RDU*100, .01))||','||abs(round(RDL*100, .01))||')'));*/
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;

/*/*/*/*/*/*/* SIRS temp */*/*/*/*/*/*/;

options mprint symbolgen;
*Generic macro;
%let Tx=Tx0;
%let cov=sirs_temp;	
%let adj=cancer;
%let ref=(ref='0');
proc glimmix data=rr.sepsis_risk_VA_KP_RR_noshock;
class &Tx. &cov.&ref. uniqid cancer;
model mort30_ed(event='1')=&Tx. &cov.*&Tx. &cov. AGE &adj.
SIRS_PULSE SIRS_RR SIRS_WBC AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr /*AOD_IND*/ MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD/*AH*/
LYMPHOMA RA COAG OBESITY WTLOSS FEN anemia_new 
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//dist=poisson link=log;
lsmeans &cov.*&Tx./cl ilink diff ;
slice &cov.*&Tx./sliceby=&cov. exp means cl;
lsmeans &cov./cl ilink diff;
lsmeans &Tx./pdiff cl;
ods output lsmeans=sirstemplsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=sirstempdiffs;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

proc print data=sirstemplsm;
run;
proc print data=sirstempdiffs;
run;

*Create table for plotting;
*LSM;
proc sort data=sirstemplsm ;
by effect sirs_temp &tx. ;
run;
data sirstemplsm_v2;
set sirstemplsm;
where Tx0 = 0 and effect ne 'Tx0';
keep sirs_temp Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL
	   sirs_temp = Level;
run;

*DIFFSA;
data sirstempdiffsa;
set sirstempdiffs;
if sirs_temp eq _sirs_temp;
run;
proc sort data=sirstempdiffsa;
by sirs_temp ;
run;
data sirstempdiffsa;
set sirstempdiffsa;
where effect ne 'Tx0';
keep sirs_temp RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
rename sirs_temp = Level;
run;

*Merge;
proc sql;
create table sirs_temp_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from sirstemplsm_v2 a
left join sirstempdiffsa b
	on a.level = b.level;
quit;

*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*these are the parameters;
/*%macro RDCI(BR, BR_LOW, BR_HI, RR, RR_LOW, RR_HI);*/

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set sirs_temp_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(sirs_temp0, 0, 0.114829512,0.1121991568,0.1175215324,0.9448061444,0.9157801089,0.9747521724);
%RDCI(sirs_temp1, 1, 0.1035172283,0.1010819559,0.1060111714,0.9068976681,0.8798643897, 0.934761527);

*merge all tables;
data sirs_temp_plot_rr_ad;
merge sirs_temp_plot_rr sirs_temp0 sirs_temp1 ;
by level;
run;

proc freq data=rr.sepsis_risk_VA_KP_RR_noshock;
tables sirs_temp;
run;

 *AR and RR columns;
data tab2;
set sirs_temp_plot_rr_ad;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;
