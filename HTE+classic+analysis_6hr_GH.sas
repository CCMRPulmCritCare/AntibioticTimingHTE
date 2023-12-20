/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
* 
* Project: HTE Analysis - limiting cohort to ABX in 6 hours
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
*NOTE: The import data set has 1560126 observations and 143 variables.;

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
*restrict to abx in 6 hr;

data rr.happi_hte_6hr_cohort;
set sepsis2;
if time_to_abx_min <=360 then abx_in_6hr =1;
else abx_in_6hr =0;
if (hospid not in ("523A5", "542", "598A0", "657A0")) and (abx_in_6hr =1) and (septic_shock = 1 or severe_sepsis = 1);
run;
/*NOTE: The data set has 218111 observations and 150 variables.*/

proc sql;
select count(distinct uniqid)
from rr.happi_hte_6hr_cohort;
quit;

proc sql;
select count(distinct uniqid), data
from rr.happi_hte_6hr_cohort
group by data;
quit;

/*/*/*/*Create Tables */*/*/*/;
proc means data=rr.happi_hte_6hr_cohort n median p25 p75 maxdec=1;
class tx0;
var TIME_TO_ABX_HR age;
ways 0 1;
run;

proc npar1way data = rr.happi_hte_6hr_cohort wilcoxon;
class tx0;
var TIME_TO_ABX_HR age ;
output out=ps wilcoxon;
run;

*cat vars;
proc freq data=rr.happi_hte_6hr_cohort;
tables tx0  / nocum;
run;
proc tabulate data=rr.happi_hte_6hr_cohort /*missing order=freq*/
format =8.1;
class  tx0 sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
	coag cancer_met cancer_nonmet;

table sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
	coag cancer_met cancer_nonmet
all ='Column Total',
all(colpctn) tx0*(n colpctn) / nocellmerge;
run;

proc freq data=rr.happi_hte_6hr_cohort;
table cancer*CANCER_MET*CANCER_NONMET/list;
table Tx0 Tx*Tx0/missing;
table aod ageg sirs_temp;
table AOD_LUNG*AOD_KIDNEY*AOD_LIVER*AOD_HEME*AOD_LACTATE*pressor_in_72hr/list;
run;

proc freq data=rr.happi_hte_6hr_cohort;
tables tx0*(male sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
		coag cancer_met cancer_nonmet)
	/ chisq exact;
run;

proc freq data=rr.happi_hte_6hr_cohort;
table aod liver coag sirs_temp;
run;

proc sql;
select count(uniqid)
from rr.happi_hte_6hr_cohort
where CANCER_MET = 1 and CANCER_NONMET = 1;
quit;

*lump 2 anemia vars to single condition;
proc print data= rr.happi_hte_6hr_cohort (obs=100);
var cancer:;
run;
*want to get counts of the categorial cancer var instead of the individual dichotomous ones for tables;
proc freq data=rr.happi_hte_6hr_cohort;
tables cancer cancer*tx0 cancer_met*cancer_nonmet;
run;
*get p value - will have to create new dich var ;
data cancer;
set rr.happi_hte_6hr_cohort;
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
set rr.happi_hte_6hr_cohort;
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
output out=rr.sepsis_risk_VA_KP_cancer_v4_RR pred=risk;
run;

*KP + VA data set for causal forest in R;
data new.causalforest_ptid_R;
set new.sepsis_risk_VA_KP_cancer_v4;
if cancer_met = 1 then cancer_nonmet = 0;
	else cancer_nonmet = cancer_nonmet;
keep patientid uniqid Tx Tx0 AOD AGEG
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr HTN
MALE AGE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD
LYMPHOMA CANCER_MET CANCER_NONMET RA COAG OBESITY WTLOSS FEN anemia_new
mort30_ed survival SEPTIC_SHOCK;
run;

/*Model with no Interaction*/
proc glimmix data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set merge_RR;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(tx0, 1, 0.1282146902,0.1258312734,0.1306432521,0.9370062221, 0.918067935,0.9563351761);

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

proc glimmix data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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
%RDCI(aod1, 1, 0.0966480879, 0.094195839,0.0991641774,0.9788431485,0.9488686687,1.0097645132);
%RDCI(aod2, 2, 0.1537130371,0.1488011907,0.1587870208,0.9738121178, 0.936614699,1.0124868227);
%RDCI(aod3, 3, 0.3200466889,0.3101085248,0.3303033451,0.8788959675,0.8475117297,0.9114423961);

*merge all tables;
data sumaod_plot_rr_ad;
merge sumaod_plot_rr aod1 aod2 aod3;
by level;
run;

proc freq data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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


/*/*/*/* Individual AODs */*/*/*/*/;
%let Tx = Tx0;

proc glimmix data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
class AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr &Tx. uniqid cancer;
model mort30_ed(event='1')=&Tx.
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr /*AOD_IND*/
&Tx.*AOD_LUNG &Tx.*AOD_KIDNEY &Tx.*AOD_LIVER &Tx.*AOD_HEME &Tx.*AOD_LACTATE &Tx.*pressor_in_72hr 
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
MALE AGE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA cancer RA COAG OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//dist=poisson link=log ;
slice AOD_LUNG*&Tx./sliceby=AOD_LUNG exp means cl;
estimate 'AOD_LUNG*Tx' AOD_LUNG*&Tx. 1 -1 -1 1/ilink cl ;
slice AOD_KIDNEY*&Tx./sliceby=AOD_KIDNEY exp means cl;
estimate 'AOD_KIDNEY*Tx' AOD_KIDNEY*&Tx. 1 -1 -1 1/ilink cl ;
slice AOD_LIVER*&Tx./sliceby=AOD_LIVER exp means cl;
estimate 'AOD_LIVER*Tx' AOD_LIVER*&Tx. 1 -1 -1 1/ilink cl ;
slice AOD_HEME*&Tx./sliceby=AOD_HEME exp means cl;
estimate 'AOD_HEME*Tx' AOD_HEME*&Tx. 1 -1 -1 1/ilink cl ;
slice AOD_LACTATE*&Tx./sliceby=AOD_LACTATE exp means cl;
estimate 'AOD_LACTATE*Tx' AOD_LACTATE*&Tx. 1 -1 -1 1/ilink cl ;
slice pressor_in_72hr*&Tx./sliceby=pressor_in_72hr exp means cl;
estimate 'pressor_in_72hr*Tx' pressor_in_72hr*&Tx. 1 -1 -1 1/ilink cl ;
lsmeans AOD_LUNG*&Tx./pdiff cl ilink;
lsmeans AOD_KIDNEY*&Tx./pdiff cl ilink;
lsmeans AOD_LIVER*&Tx./pdiff cl ilink;
lsmeans AOD_HEME*&Tx./pdiff cl ilink;
lsmeans AOD_LACTATE*&Tx./pdiff cl ilink;
lsmeans pressor_in_72hr*&Tx./pdiff cl ilink;
lsmeans &Tx./pdiff cl ;
ods output lsmeans=indaodlsm 
slices=slices slicediffs=slicediff 
slicetests=slicetest diffs=indaoddiffs estimates=est;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
/*store glim; *to get rate ratios;*/
run;

proc print data=indaodlsm;
where tx0 = 0;
run;
proc print data=indaoddiffs;
run;

*Create table for plotting;
*LSM;
data indaodlsm_v2;
length Level $20;
set indaodlsm;
where Tx0 = 0 and effect ne 'Tx0';
array var [6] aod_lung aod_kidney aod_liver aod_heme aod_lactate pressor_in_72hr;
array yesaod [6] $20 _temporary_ ('Lung' 'Kidney' 'Liver' 'Hematologic' 'Lactate' 'Shock');
array noaod [6] $20 _temporary_ ('No Lung' 'No Kidney' 'No Liver' 'No Hematologic' 'No Lactate' 'No Shock');
	do i = 1 to 6;
		if var[i] = 1 then Level = yesaod[i];
			else if var[i] = 0 then Level = noaod[i];
	end;
keep Level Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL;
run;

*DIFFSA;
data indaoddiffsa;
length Level $20;
set indaoddiffs;
array var [6] aod_lung aod_kidney aod_liver aod_heme aod_lactate pressor_in_72hr;
array var2 [6] _aod_lung _aod_kidney _aod_liver _aod_heme _aod_lactate _pressor_in_72hr;
array yesaod [6] $20 _temporary_ ('Lung' 'Kidney' 'Liver' 'Hematologic' 'Lactate' 'Shock');
array noaod [6] $20 _temporary_ ('No Lung' 'No Kidney' 'No Liver' 'No Hematologic' 'No Lactate' 'No Shock');
	do i = 1 to 6;
		if var[i] = var2[i];
			if var[i] = 1 then Level = yesaod[i];
				else if var[i] = 0 then Level = noaod[i];
	end;
where effect ne 'Tx0';
keep Level RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
run;

*Merge;
proc sql;
create table indaod_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from indaodlsm_v2 a
left join indaoddiffsa b
	on a.level = b.level;
quit;

*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set indaod_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
keep Level macstring;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(aod1, 'Hematologic', 0.30637041, 0.294245115,0.3189953658,0.8243741862,0.7847803642,0.8659655998);
%RDCI(aod2, 'Kidney',  0.296167997,0.2870838361,0.3055396071, 0.853242983,0.8213631929,0.8863601319);
%RDCI(aod3, 'Lactate', 0.2998928276,0.2908719333,0.3091934895,0.8579760854,0.8275254097,0.8895472628);
%RDCI(aod4, 'Liver',  0.2988605266, 0.285872449,0.3124386931,0.8175959178,0.7716327677,0.8662969131);
%RDCI(aod5, 'Lung', 0.3695426199, 0.354278757,0.3854641161,0.8647095808,0.8215576095,0.9101280915);
%RDCI(aod6, 'No Hematologic', 0.2183200791,0.2118351067, 0.225003578,0.9065884831,0.8731620347,0.9412945652);
%RDCI(aod7, 'No Kidney', 0.2258407823,0.2179024522,0.2340683111,0.8759147838,0.8379968729,0.9155484147);
%RDCI(aod8, 'No Lactate', 0.2230357181,0.2150219175,0.2313481905, 0.871082721, 0.830875208,0.9132359463);
%RDCI(aod9, 'No Liver', 0.2238061109, 0.217866493,0.2299076586,0.9141045433,0.8856271048,0.9434976769);
%RDCI(aod10, 'No Lung', 0.180998912,0.1753943382,0.1867825752,0.8642995979,0.8312126726,0.8987035684);
%RDCI(aod11, 'No Shock', 0.1930724798,0.1864688626,0.1999099577,0.8971724474,0.8591553519,0.9368717759);
%RDCI(aod12, 'Shock', 0.3464336929, 0.333536716,0.3598293615,0.8330261871,0.7957269704,0.8720737818);

*merge all tables;
data indaod_plot_rr_ad;
merge indaod_plot_rr aod1 aod2 aod3
 aod4 aod5 aod6 aod7 aod8 aod9 aod10 aod11 aod12;
by level;
run;
proc print data=indaod_plot_rr_ad;
run;

proc format;
invalue order
'Lactate' = 1 
'No Lactate' = 2 
'Kidney' = 3
'No Kidney' = 4
'Hematologic' = 5
'No Hematologic' = 6 
'Shock' = 7 
'No Shock' = 8
'Liver' = 9
'No Liver' = 10 
'Respiratory' = 11 
'No Respiratory' = 12;
quit;

ods output onewayfreqs= overall (drop=cumfrequency cumpercent F_:);
proc freq data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
tables aod_lung aod_kidney aod_liver aod_heme aod_lactate pressor_in_72hr;
run;
ods output close;
proc sort data=overall;
by descending aod_lung descending aod_kidney descending aod_liver 
descending aod_heme descending aod_lactate descending pressor_in_72hr;
run;
data overall2;
set overall;
length tot $ 40;
tot = (strip(Frequency)||' ('||strip(round(Percent, .1))||')');
drop table frequency percent;
run;

*AR and RR columns;
data tab2;
set indaod_plot_rr_ad;
neworder = input(level, order.);
run;
proc sort data=tab2;
by neworder;
run;
data tab2;
set tab2;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep Level AMR RelRis;
run;

/*/*/*/*/*/*/* Cancer */*/*/*/*/*/*/;
options mprint symbolgen;
*Generic macro;
%let Tx=Tx0;
%let cov=cancer;	
%let adj=sirs_temp;
%let ref=(ref='0');
proc glimmix data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set cancer_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(cancer0, 0, 0.0910174191,0.0891468565,0.0929272316,0.9828228179,0.9586302165,1.0076259592);
%RDCI(cancer1, 1, 0.1112416795,0.1068622198,0.1158006194,0.8962177876,0.8507525689,0.9441127211);
%RDCI(cancer2, 2, 0.2248366847,0.2168049186, 0.233165996, 0.819075805,0.7822943693,0.8575866076);

*merge all tables;
data cancer_plot_rr_ad;
merge cancer_plot_rr cancer0 cancer1 cancer2;
by level;
run;

proc freq data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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
proc glimmix data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set sirs_temp_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(sirs_temp0, 0, 0.1350570255,0.1317489694,0.1384481429,0.9617420677,0.9333678195, 0.990978889);
%RDCI(sirs_temp1, 1, 0.1240788847,0.1211166414,0.1271135777, 0.916870296,0.8921565084,0.9422686848);

*merge all tables;
data sirs_temp_plot_rr_ad;
merge sirs_temp_plot_rr sirs_temp0 sirs_temp1 ;
by level;
run;

*Overall column;
proc freq data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
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


/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/* AGE */*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/;

proc summary data=rr.sepsis_risk_VA_KP_cancer_v4_rr nway;
class ageg tx0;
var mort30_ed;
output out=mr mean=;
run;

proc transpose data=mr out=mrt;
by ageg;
var mort30_ed;
id tx0;
run;

proc glimmix data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
class &Tx. CANCER ageg uniqid;
model mort30_ed(event='1')=&Tx. ageg*&Tx. ageg 
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD/*AH*/
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//
dist=poisson link=log solution;
lsmeans ageg*&Tx./cl ilink diff ;
slice ageg*&Tx./sliceby=ageg exp means cl;
lsmeans ageg/cl ilink diff;
lsmeans &Tx./pdiff cl;
ods output lsmeans=ageglsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=agegdiffs;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

proc print data=ageglsm;
where tx0 = 0;
run;
proc print data=agegdiffs;
run;

*Create table for plotting;
*LSM;
proc sort data=ageglsm ;
by effect ageg &tx. ;
run;
data ageglsm_v2;
set ageglsm;
where Tx0 = 0 and effect ne 'Tx0';
keep ageg Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL
	   ageg = Level;
run;

*DIFFSA;
data agegdiffsa;
set agegdiffs;
if ageg eq _ageg;
run;
proc sort data=agegdiffsa;
by ageg ;
run;
data agegdiffsa;
set agegdiffsa;
where effect ne 'Tx0';
keep ageg RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
rename ageg = Level;
run;

*Merge;
proc sql;
create table ageg_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from ageglsm_v2 a
left join agegdiffsa b
	on a.level = b.level;
quit;

*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set ageg_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(ageg1, 1, 0.0505295489,0.0434057995,0.0588224463,0.9021082871, 0.748307467,1.0875200337);
%RDCI(ageg3, 3, 0.0725664789,0.0650306536,0.0809755641,0.8802919291,0.7666294069,1.0108063604);
%RDCI(ageg4, 4, 0.093393177,0.0883991386,0.0986693497,0.8309841112,0.7748812811,0.8911488894);
%RDCI(ageg5, 5, 0.1130642378,0.1094821062,0.1167635727,0.9051900786,0.8692824799,0.9425809185);
%RDCI(ageg6, 6, 0.1396286227, 0.135111405,0.1442968657, 0.946106963,0.9087565721,0.9849924752);
%RDCI(ageg7, 7, 0.2025032988, 0.195962747,0.2092621514,0.9742925083,0.9373881632,1.0126497528);
%RDCI(ageg8, 8, 0.2923150645,0.2784702082, 0.306848253,0.9809147859,0.9267316219,1.0382658738);		

*merge all tables;
data ageg_plot_rr_ad;
merge ageg_plot_rr ageg1 ageg3 ageg4 ageg5 ageg6 ageg7 ageg8;
by level;
run;

*Overall column;
proc freq data=rr.sepsis_risk_VA_KP_cancer_v4_rr;
tables ageg;
run;

 *AR and RR columns;
data tab2;
set ageg_plot_rr_ad;
length AMR $ 25;
/*AMR = compress((abs(round(RD*100, .01))||'('||abs(round(RDU*100, .01))||','||abs(round(RDL*100, .01))||')'));*/
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;


/*/*/*/*/*/*/*/*/*/*/*/**Unadjusted continuous age HTE*/*/*/*/*/*/*/*/*/*/*/;

data sepsis_risk0_v3;
set rr.sepsis_risk_VA_KP_cancer_v4_rr;
if age ge 99 then age=99;
run;

*Absolute Risk;
proc summary data=sepsis_risk0_v3;
var SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA RA COAG OBESITY WTLOSS FEN anemia_new 
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*/;
output out=base mean=;
run;

%let Tx = Tx0;
data base0;
set base;
do age=18 to 105;
do &Tx.=0 to 1;
cancer=1;
output;
end;
end;
drop _:;
run;

data sepsis_riskp;
set sepsis_risk0_v3 base0;
*sage=sqrt(age);
keep &Tx. age  mort30_ed
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*/ uniqid;
run;

*AGE CONTINUOUS;
proc glimmix data=sepsis_riskp;
class &Tx. CANCER uniqid;
*effect sage=spline(age);*/KNOTMETHOD=list(50 60 70 80 90) basis=tpf);
effect sage=poly(age/degree=2);
model mort30_ed(event='1')=&Tx. sage*&Tx. sage 
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD /*AH*/
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new
/*ANEMIA_CBL ANEMIA_DEF ETOH DRUG PSYCHOSES DEPRESSION*//
dist=poisson link=log solution;
output out=preds predicted(blup ilink)=pmort
stderr(blup ilink)=stde predicted(blup noilink)=xbeta
stderr(blup noilink)=xbetase;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

data predsa;
set preds(where=(mort30_ed eq .));
run;

proc transpose data=predsa out=estimate;
where mort30_ed eq .;
by age;
var pmort;
id Tx0;
run;

proc transpose data=predsa out=error;
where mort30_ed eq .;
by age;
var stde;
id Tx0;
run;

data morts;
merge estimate(rename=('0'n=control '1'n=case))
error(rename=('0'n=controlse '1'n=casese));
by age;
diff=case-control;
diff_se=sqrt(controlse**2+casese**2);
lower=diff-1.96*diff_se;
upper=diff+1.96*diff_se;
run;
