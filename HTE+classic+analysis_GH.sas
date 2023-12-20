/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
* 
* Project: HTE Analysis
*
* Author: Jennifer Cano (modified Patricia Kipnis's code)
*
* Description: Running robust Poisson models/creating figures
* 
/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/;

libname h 'filepath';
libname new 'filepath';

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
*remove facilities with <15 sepsis hospitalizations;

data h.happi_hte_cohort;
set sepsis2;
if time_to_abx_min <=720 then abx_in_12hr =1;
else abx_in_12hr =0;
if (hospid not in ("523A5", "542", "598A0", "657A0")) and (abx_in_12hr =1) and (septic_shock = 1 or severe_sepsis = 1);
run;
/*NOTE: The data set WORK.SEPSIS2 has 273255 observations and 150 variables.*/

*Fill in eTable 1;
proc sql;
select count(distinct uniqid), data
from sepsis2
where (hospid not in ("523A5", "542", "598A0", "657A0"))
group by data;
quit;
proc sql;
select count(distinct uniqid)
from h.happi_hte_cohort;
quit;
proc sql;
select count(distinct uniqid), data
from h.happi_hte_cohort
group by data;
quit;

/*/*/*/*Create Tables */*/*/*/;

*Table 1c - new cohort;
*cont vars;
proc means data=h.happi_hte_cohort n median p25 p75 maxdec=1;
class tx0;
var TIME_TO_ABX_HR age;
ways 0 1;
run;

*p values;
proc npar1way data = h.happi_hte_cohort wilcoxon;
class tx0;
var TIME_TO_ABX_HR age ;
output out=ps wilcoxon;
run;

*cat vars;
proc freq data=h.happi_hte_cohort;
tables tx0  / nocum;
run;

proc tabulate data=h.happi_hte_cohort 
format =8.1;
class  tx0 sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
	coag cancer_met cancer_nonmet;

table sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
	coag cancer_met cancer_nonmet
all ='Column Total',
all(colpctn) tx0*(n colpctn) / nocellmerge;
run;

proc freq data=h.happi_hte_cohort;
table cancer*CANCER_MET*CANCER_NONMET/list;
table Tx0 Tx*Tx0/missing;
table aod ageg sirs_temp;
table AOD_LUNG*AOD_KIDNEY*AOD_LIVER*AOD_HEME*AOD_LACTATE*pressor_in_72hr/list;
run;

*p values;
proc freq data=h.happi_hte_cohort;
tables tx0*(male sirs_wbc sirs_temp sirs_pulse sirs_rr aod: pressor_in_72hr renal liver chf pulm htn obesity pvd cardic_arrhym
		coag cancer_met cancer_nonmet)
	/ chisq exact;
run;

proc freq data=h.happi_hte_cohort;
table aod liver coag sirs_temp;
run;

proc sql;
select count(uniqid)
from h.happi_hte_cohort
where CANCER_MET = 1 and CANCER_NONMET = 1;
quit;

*lump 2 anemia vars to single condition;
proc print data= h.happi_hte_cohort (obs=100);
var cancer:;
run;
*want to get counts of the categorial cancer var instead of the individual dichotomous ones for tables;
proc freq data=h.happi_hte_cohort;
tables cancer cancer*tx0 cancer_met*cancer_nonmet;
run;
*get p value - will have to create new dich var ;
data cancer;
set h.happi_hte_cohort;
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
set h.happi_hte_cohort;
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
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new;
output out=new.sepsis_risk_VA_KP_cancer_v4 pred=risk;
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


/*Model with no Interaction for Table 2*/

proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class Tx0(ref='0') cancer uniqid;
model mort30_ed(event='1')=Tx0 
SIRS_PULSE SIRS_RR SIRS_WBC SIRS_TEMP AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr MALE AGE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD
LYMPHOMA RA COAG CANCER OBESITY WTLOSS FEN anemia_new /dist=poisson link=log;
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
%RDCI(tx0, 1, 0.1310723854, 0.129093687,0.1330814125,0.9101124278,0.8938755953,0.9266441947);

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
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep AMR RelRis;
run;

* crude mortality;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables mort30_ed;
run;


/*/*/*/*/* Sum of AODs */*/*/*/*/;

%let Tx=Tx0;

proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class AOD(ref='1') &Tx. uniqid cancer;
model mort30_ed(event='1')=&Tx. AOD*&Tx. AOD
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
MALE AGE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA cancer RA COAG OBESITY WTLOSS FEN anemia_new /dist=poisson link=log;
lsmeans AOD*&Tx./cl ilink diff ;
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

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set sumaod_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(aod1, 1, 0.0999659904, 0.097963868,0.1020090308,0.9480270774,0.9225426468,0.9742154931);
%RDCI(aod2, 2, 0.1588451445, 0.154809341,0.1629861594, 0.943446951,0.9119380444,0.9760445402);
%RDCI(aod3, 3, 0.3285693968,0.3203701776,0.3369784581, 0.855255249,0.8286442978,0.8827207801);

*merge all tables;
data sumaod_plot_rr_ad;
merge sumaod_plot_rr aod1 aod2 aod3;
by level;
run;

*get numbers for Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables aod;
run;

 *AR and RR columns;
data tab2;
set sumaod_plot_rr_ad;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep AMR RelRis;
run;

*crude mortality;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables aod*mort30_ed;
run;

*create plots;
proc format;
value numaod
1 = "1"
2 = "2"
3 = "3+";
quit;

*AR;
ods graphics on;
proc sgplot data=sumaod_plot_rr_ad;
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12 color=deypk)
							errorbarattrs=(color=black) ;
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Number of Dysfunctional Organs" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.2 offsetmin=0.2 
			type=discrete;
format level numaod.;
refline 0;
run;

*RR;
ods graphics on ;
proc sgplot data=sumaod_plot_rr_ad;
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) ;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Number of Dysfunctional Organs" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.2 offsetmin=.2 
			type=discrete;
format level numaod.;
refline 1;
run;



/*/*/*/*/* Individual AODs */*/*/*/*/;
%let Tx = Tx0;

proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr &Tx. uniqid cancer;
model mort30_ed(event='1')=&Tx.
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
&Tx.*AOD_LUNG &Tx.*AOD_KIDNEY &Tx.*AOD_LIVER &Tx.*AOD_HEME &Tx.*AOD_LACTATE &Tx.*pressor_in_72hr 
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
MALE AGE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA cancer RA COAG OBESITY WTLOSS FEN anemia_new /dist=poisson link=log ;
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
%RDCI(aod1, 'Hematologic', 0.3192418526,0.3093927733,0.3294044633,0.7822872458,0.7495567304,0.8164469881);
%RDCI(aod2, 'Kidney',  0.305744549,0.2983805931,0.3132902454,0.8204625456,0.7936110908,0.8482225067);
%RDCI(aod3, 'Lactate', 0.3084468051,0.3010687469,0.3160056715,0.8325022407,0.8065926948,0.8592440587);
%RDCI(aod4, 'Liver',  0.3121375139,0.3017252318, 0.322909115,0.7757646295,0.7373830608,0.8161439995);
%RDCI(aod5, 'Lung', 0.3778847269,0.3655649091, 0.390619732,0.8422831301,0.8056246836,0.8806096508);
%RDCI(aod6, 'No Hematologic',  0.2245699223,0.2193377061,0.2299269509,0.8826611888,0.8540265287,0.9122559405);
%RDCI(aod7, 'No Kidney', 0.2344837161,0.2280206609,0.2411299612,0.8415918484,0.8096440976,0.8748002257);
%RDCI(aod8, 'No Lactate', 0.2324294395,0.2259872884,0.2390552352,0.8294206989,0.7957713575,0.8644929091);
%RDCI(aod9, 'No Liver', 0.2296811977,0.2248224442,0.2346449562,0.8900825896,0.8658478155,0.9149956865);
%RDCI(aod10, 'No Lung', 0.1897195439,0.1851266003,0.1944264372,0.8197891726,0.7922616407,0.8482731625);
%RDCI(aod11, 'No Shock', 0.2001793521,0.1948428713,0.2056619918,0.8589684446,0.8271561896,0.8920041923);
%RDCI(aod12, 'Shock', 0.3581394248,0.3475456053,0.3690561632,0.8038649088,0.7723263936,0.8366913224);

*merge all tables;
data indaod_plot_rr_ad;
merge indaod_plot_rr aod1 aod2 aod3
 aod4 aod5 aod6 aod7 aod8 aod9 aod10 aod11 aod12;
by level;
run;
proc print data=indaod_plot_rr_ad;
run;

*get numbers for Table 2 in 'TablesFigures_TTA_2023.06.19_JC.docx';
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


*get numbers for Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
ods output onewayfreqs= overall (drop=cumfrequency cumpercent F_:);
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
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

*crude mortality;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables aod_lung*mort30_ed aod_kidney*mort30_ed aod_liver*mort30_ed aod_heme*mort30_ed aod_lactate*mort30_ed pressor_in_72hr*mort30_ed;
run;

*create plots;
data indaod_plot_rr_ad;
set indaod_plot_rr_ad;
if level in ('Lactate' 'No Lactate' 'Hematologic' 'No Hematologic' 'Liver' 'No Liver') then bl = 0;
	else if level in ('Kidney' 'No Kidney' 'Shock' 'No Shock' 'Lung' 'No Lung') then bl = 1;
if level = 'Lung' then level = 'Respiratory';
if level = 'No Lung' then level  = 'No Respiratory';
run;


*AR;
ods graphics on / width=8in height=6in;
proc sgplot data=indaod_plot_rr_ad noautolegend;
block x= level block = bl / filltype= alternate  fillattrs=(color=white) altfillattrs=(color=cxcacfd2) 
									nooutline novalues transparency=0.75;
styleattrs datacontrastcolors= (cx993329 cx993329 cx993329 cx993329 cx993329 cx003178 cx003178 cx003178 cx003178 cx003178 cx003178 cx993329);
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12) group = level 
						errorbarattrs=(color=black);
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Organ" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.08 offsetmin=.08 type=discrete
	values=('Lactate' 'No Lactate' 'Kidney' 'No Kidney' 'Hematologic' 'No Hematologic' 'Shock' 'No Shock'
			'Liver' 'No Liver' 'Respiratory' 'No Respiratory');
refline 0;
format rd rdu rdl percentn5.2;
run;


*RR;
ods graphics on / width=8in height=6in;
proc sgplot data=indaod_plot_rr_ad noautolegend;
block x= level block = bl / filltype= alternate  fillattrs=(color=white) altfillattrs=(color=cxcacfd2) 
								nooutline novalues transparency=0.75;
styleattrs datacontrastcolors= (cx993329 cx993329 cx993329 cx993329 cx993329 cx003178 cx003178 cx003178 cx003178 cx003178 cx003178 cx993329);
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) group= level;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Organ" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.08 offsetmin=0.08 
	values=('Lactate' 'No Lactate' 'Kidney' 'No Kidney' 'Hematologic' 'No Hematologic' 'Shock' 'No Shock'
			'Liver' 'No Liver' 'Respiratory' 'No Respiratory');
refline 1;
run;


/*/*/*/*/*/*/* Cancer */*/*/*/*/*/*/;
options mprint symbolgen;
%let Tx=Tx0;
%let cov=cancer;	
%let adj=sirs_temp;
%let ref=(ref='0');
proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class &Tx. &cov.&ref. uniqid;
model mort30_ed(event='1')=&Tx. &cov.*&Tx. &cov. AGE &adj.
SIRS_PULSE SIRS_RR SIRS_WBC AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA RA COAG OBESITY WTLOSS FEN anemia_new /dist=poisson link=log;
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
%RDCI(cancer0, 0, 0.0926420642,0.0911305192,0.0941786805,0.9567934049,0.9361521259,0.9778898049);
%RDCI(cancer1, 1,  0.1120394593,0.1085765563,0.1156128069,0.8788947833, 0.839400908,0.9202468484);
%RDCI(cancer2, 2, 0.2308130692,0.2243354951,  0.23747768, 0.783306759,0.7523309253,0.8155579654);

*merge all tables;
data cancer_plot_rr_ad;
merge cancer_plot_rr cancer0 cancer1 cancer2;
by level;
run;

*get numbers to Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables cancer;
run;

 *AR and RR columns;
data tab2;
set cancer_plot_rr_ad;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;

*crude mortality;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables cancer*mort30_ed;
run;

*create plots;
proc format;
value numcancer
0 = "No Cancer"
1 = "Non-Metastatic Cancer"
2 = "Metastatic Cancer";
quit;

*AR;
ods graphics on ;
proc sgplot data=cancer_plot_rr_ad;
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12 color=deypk)
							errorbarattrs=(color=black) ;
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Cancer Status" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.2 offsetmin=0.2 
			type=discrete fitpolicy=split;
format level numcancer.;
refline 0;
run;

*RR;
ods graphics on ;
proc sgplot data=cancer_plot_rr_ad;
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) ;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Cancer Status" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.2 offsetmin=.2 
			type=discrete fitpolicy=split;
format level numcancer.;
refline 1;
run;

/*/*/*/*/*/*/* SIRS temp */*/*/*/*/*/*/;
options mprint symbolgen;
%let Tx=Tx0;
%let cov=sirs_temp;	
%let adj=cancer;
%let ref=(ref='0');
proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class &Tx. &cov.&ref. uniqid cancer;
model mort30_ed(event='1')=&Tx. &cov.*&Tx. &cov. AGE &adj.
SIRS_PULSE SIRS_RR SIRS_WBC AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD
LYMPHOMA RA COAG OBESITY WTLOSS FEN anemia_new /dist=poisson link=log;
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
%RDCI(sirs_temp0, 0, 0.1367495403,0.1340867956,0.1394651628,0.9429146194,0.9183375832,0.9681493992);
%RDCI(sirs_temp1, 1, 0.1274496956,0.1249996772,0.1299477349,0.8844350645,0.8636052209,0.9057673164);

*merge all tables;
data sirs_temp_plot_rr_ad;
merge sirs_temp_plot_rr sirs_temp0 sirs_temp1 ;
by level;
run;

*get numbers to Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
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

*crude mortality;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables sirs_temp*mort30_ed;
run;

*create plots;
proc format;
value numsirs_temp
0 = "Normal Temp"
1 = "Abnormal SIRS Temp";
quit;

*AR;
ods graphics on ;
proc sgplot data=sirs_temp_plot_rr_ad;
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12 color=deypk)
							errorbarattrs=(color=black) ;
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Body Temperature" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.33 offsetmin=0.33 
			type=discrete fitpolicy=split;
format level numsirs_temp.;
refline 0;
run;

*RR;
ods graphics on ;
proc sgplot data=sirs_temp_plot_rr_ad;
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) ;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Body Temperature" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.33 offsetmin=.33 
			type=discrete fitpolicy=split;
format level numsirs_temp.;
refline 1;
run;


/*/*/*/*/*/*/* Liver */*/*/*/*/*/*/;

*freqs;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
table liver coag/list nocum;
run;

options mprint symbolgen;
%let Tx=Tx0;
%let cov=liver;	
%let adj=coag;
%let ref=(ref='0');
proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class &Tx. &cov.&ref. uniqid cancer;
model mort30_ed(event='1')=&Tx. &cov.*&Tx. &cov. AGE &adj.
SIRS_PULSE SIRS_RR SIRS_WBC SIRS_TEMP AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL PUD cancer 
LYMPHOMA RA  OBESITY WTLOSS FEN anemia_new /dist=poisson link=log;
lsmeans &cov.*&Tx./cl ilink diff ;
slice &cov.*&Tx./sliceby=&cov. exp means cl;
lsmeans &cov./cl ilink diff;
lsmeans &Tx./pdiff cl;
ods output lsmeans=liverlsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=liverdiffs;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

proc print data=liverlsm;
run;
proc print data=liverdiffs;
run;

*Create table for plotting;
*LSM;
proc sort data=liverlsm ;
by effect liver &tx. ;
run;
data liverlsm_v2;
set liverlsm;
where Tx0 = 0 and effect ne 'Tx0';
keep liver Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL
	   liver = Level;
run;

*DIFFSA;
data liverdiffsa;
set liverdiffs;
if liver eq _liver;
run;
proc sort data=liverdiffsa;
by liver ;
run;
data liverdiffsa;
set liverdiffsa;
where effect ne 'Tx0';
keep liver RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
rename liver = Level;
run;

*Merge;
proc sql;
create table liver_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from liverlsm_v2 a
left join liverdiffsa b
	on a.level = b.level;
quit;

*Add LAD and UAD columns from Excel Macro;
*First run Patricia's macro;
%include 'filepath';

*get RD, RDL, and RDU;

*create macstring var to then copy and paste into the macro parameters below;
data rdmacro;
set liver_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(liver0, 0, 0.1235346865,0.1214323737,0.1256733959,0.9538574217,0.9342640897,0.9738616639);
%RDCI(liver1, 1, 0.1599155656, 0.155878288,0.1640574094, 0.797382409,0.7705369164,0.8251631982);

*merge all tables;
data liver_plot_rr_ad;
merge liver_plot_rr liver0 liver1 ;
by level;
run;

*get numbers to Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables liver;
run;

 *AR and RR columns;
data tab2;
set liver_plot_rr_ad;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;

*create plots;
proc format;
value numliver
0 = "No Liver Disease"
1 = "Liver Disease";
quit;

*AR;
ods graphics on ;
proc sgplot data=liver_plot_rr_ad;
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12 color=deypk)
							errorbarattrs=(color=black) ;
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Disease" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.33 offsetmin=0.33 
			type=discrete fitpolicy=split;
format level numliver.;
refline 0;
run;

*RR;
ods graphics on ;
proc sgplot data=liver_plot_rr_ad;
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) ;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Disease" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.33 offsetmin=.33 
			type=discrete fitpolicy=split;
format level numliver.;
refline 1;
run;


/*/*/*/*/*/*/* Coag */*/*/*/*/*/*/;

*freqs;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
table liver coag/list nocum;
run;

options mprint symbolgen;
%let Tx=Tx0;
%let cov=coag;	
%let adj=liver;
%let ref=(ref='0');
proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class &Tx. &cov.&ref. uniqid cancer;
model mort30_ed(event='1')=&Tx. &cov.*&Tx. &cov. AGE &adj.
SIRS_PULSE SIRS_RR SIRS_WBC SIRS_TEMP AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME 
AOD_LACTATE pressor_in_72hr MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL PUD cancer 
LYMPHOMA RA  OBESITY WTLOSS FEN anemia_new /dist=poisson link=log;
lsmeans &cov.*&Tx./cl ilink diff ;
slice &cov.*&Tx./sliceby=&cov. exp means cl;
lsmeans &cov./cl ilink diff;
lsmeans &Tx./pdiff cl;
ods output lsmeans=coaglsm 
slices=slices slicediffs=slicediff slicetests=slicetest diffs=coagdiffs;
random _residual_ / subject=uniqid type=cs; *use of robust poisson;
run;

proc print data=coaglsm;
run;
proc print data=coagdiffs;
run;

*Create table for plotting;
*LSM;
proc sort data=coaglsm ;
by effect coag &tx. ;
run;
data coaglsm_v2;
set coaglsm;
where Tx0 = 0 and effect ne 'Tx0';
keep coag Mu LowerMu UpperMu;
rename Mu = BL_RR 
	   LowerMu = BL_LCL
	   UpperMu = BL_UCL
	   coag = Level;
run;

*DIFFSA;
data coagdiffsa;
set coagdiffs;
if coag eq _coag;
run;
proc sort data=coagdiffsa;
by coag ;
run;
data coagdiffsa;
set coagdiffsa;
where effect ne 'Tx0';
keep coag RR RR_LCL RR_UCL;
RR = exp(-Estimate);
RR_LCL = exp(-Upper);
RR_UCL = exp(-Lower);
rename coag = Level;
run;

*Merge;
proc sql;
create table coag_plot_rr as
select a.Level, a.BL_RR, a.BL_LCL, a.BL_UCL, b.RR, b.RR_LCL, b.RR_UCL
from coaglsm_v2 a
left join coagdiffsa b
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
set coag_plot_rr;
macstring= BL_RR || ',' || BL_LCL || ',' || BL_UCL || ',' || RR || ',' || RR_LCL || ',' || RR_UCL ;
run;

*run macro to create data sets for each level and then merge all data sets together;
%RDCI(coag0, 0, 0.1249478292,0.1228543999,0.1270769303,0.9429002032,0.9231449925,0.9630781735);
%RDCI(coag1, 1, 0.1545676684,   0.1505078,0.1587370497,0.8358849693,0.8093278755,0.8633135012);

*merge all tables;
data coag_plot_rr_ad;
merge coag_plot_rr coag0 coag1 ;
by level;
run;

*get numbers to Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables coag;
run;

 *AR and RR columns;
data tab2;
set coag_plot_rr_ad;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;

*create plots;
proc format;
value numcoag
0 = "No Pre-Existing Coagulopathy"
1 = "Pre-Existing Coagulopathy";
quit;

*AR;
ods graphics on;
proc sgplot data=coag_plot_rr_ad;
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12 color=deypk)
							errorbarattrs=(color=black) ;
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Coagulopathy History" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.4 offsetmin=0.4 
			type=discrete fitpolicy=split;
format level numcoag.;
refline 0;
run;

*RR;
ods graphics on;
proc sgplot data=coag_plot_rr_ad;
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) ;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Coagulopathy History" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.4 offsetmin=.4 
			type=discrete fitpolicy=split;
format level numcoag.;
refline 1;
run;


/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/* AGE */*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/;

proc summary data=new.sepsis_risk_VA_KP_cancer_v4 nway;
class ageg tx0;
var mort30_ed;
output out=mr mean=;
run;

proc transpose data=mr out=mrt;
by ageg;
var mort30_ed;
id tx0;
run;

proc glimmix data=new.sepsis_risk_VA_KP_cancer_v4;
class &Tx. CANCER ageg uniqid;
model mort30_ed(event='1')=&Tx. ageg*&Tx. ageg 
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new /
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
%RDCI(ageg1, 1, 0.0584221806,0.0519993681,0.0656383204,0.8292082407,0.7064197856,0.9733395361);
%RDCI(ageg3, 3, 0.0775313152,0.0712905805,0.0843183602,0.8454123861,0.7502274626,0.9526738732);
%RDCI(ageg4, 4, 0.0968760871,0.0928985297,0.1010239482,0.8049137787,0.7575098313,0.8552842015);
%RDCI(ageg5, 5, 0.1167196218,0.1138505785,0.1196609653,0.8688992859,0.8387327841, 0.900150779);
%RDCI(ageg6, 6,  0.142568252,0.1389393051,0.1462919832, 0.916259523,0.8846755742,0.9489710555);
%RDCI(ageg7, 7, 0.2035691073,0.1982370931,0.2090445376,0.9584734754,0.9264665891,0.9915861121);
%RDCI(ageg8, 8, 0.2939156929,0.2824660017,0.3058294947,0.9667172079,0.9195515052,1.0163021373);


*merge all tables;
data ageg_plot_rr_ad;
merge ageg_plot_rr ageg1 ageg3 ageg4 ageg5 ageg6 ageg7 ageg8;
by level;
run;

*get numbers to Table 2 of TablesFigures_TTA_2023.06.19_JC saved in 'filepath';
*Overall column;
proc freq data=new.sepsis_risk_VA_KP_cancer_v4;
tables ageg;
run;

 *AR and RR columns;
data tab2;
set ageg_plot_rr_ad;
length AMR $ 25;
AMR = (abs(round(RD*100, .01))||' ('||strip(abs(round(RDU*100, .01)))||', '||strip(abs(round(RDL*100, .01))||')'));
RelRis = (abs(round(RR, .01))||' ('||strip(abs(round(RR_LCL, .01)))||', '||strip(abs(round(RR_UCL, .01))||')'));
keep level AMR RelRis;
run;

*create plots;
proc format;
value numageg
1 = "<40"
3 = "40-49"
4 = "50-59"
5 = "60-69"
6 = "70-79"
7 = "80-89"
8 = "90+";
quit;

*AR;
ods graphics on / width=8in height=6in;
proc sgplot data=ageg_plot_rr_ad;
scatter x=level y=rd / yerrorupper=rdu yerrorlower= rdl  markerattrs=(symbol=circlefilled size=12 color=deypk)
							errorbarattrs=(color=black) ;
yaxis label="Absolute Mortality Difference" values=(-.07 to .01 by .01) grid labelattrs=(size=16) valueattrs=(size=16)
			valuesformat=percentn12.1;
xaxis label="Age" labelattrs=(size=16) valueattrs=(size=16) offsetmax=0.08 offsetmin=0.08 
			type=discrete fitpolicy=split;
format level numageg.;
refline 0;
run;


*RR;
ods graphics on / width=8in height=6in;
proc sgplot data=ageg_plot_rr_ad;
scatter x=level y=rr / yerrorupper=rr_ucl yerrorlower= rr_lcl  markerattrs=(symbol=circlefilled size=12) errorbarattrs=(color=black) ;
yaxis label="Relative Risk" values=(0.6 to 1.1 by 0.1) grid labelattrs=(size=16) valueattrs=(size=16);
xaxis label="Age" labelattrs=(size=16) valueattrs=(size=16) offsetmax=.08 offsetmin=.08 
			type=discrete fitpolicy=split;
format level numageg.;
refline 1;
run;



/*/*/*/*/*/*/*/*/*/*/*/**Unadjusted continuous age HTE*/*/*/*/*/*/*/*/*/*/*/;

data sepsis_risk0_v3;
set new.sepsis_risk_VA_KP_cancer_v4;
if age ge 99 then age=99;
run;

*Absolute Risk;
proc summary data=sepsis_risk0_v3;
var SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE HTN CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA RA COAG OBESITY WTLOSS FEN anemia_new ;
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
keep &Tx. age  mort30_ed
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new uniqid;
run;

*AGE CONTINUOUS;
proc glimmix data=sepsis_riskp;
class &Tx. CANCER uniqid;
effect sage=poly(age/degree=2);
model mort30_ed(event='1')=&Tx. sage*&Tx. sage 
SIRS_TEMP SIRS_PULSE SIRS_RR SIRS_WBC
AOD_LUNG AOD_KIDNEY AOD_LIVER AOD_HEME AOD_LACTATE pressor_in_72hr 
MALE CHF CARDIC_ARRHYM VALVULAR_D2 PULM_CIRC PVD PARALYSIS
NEURO PULM DM_UNCOMP DM_COMP HYPOTHYROID RENAL LIVER PUD 
LYMPHOMA CANCER RA COAG OBESITY WTLOSS FEN anemia_new /
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

*AR Figure;
goptions reset=all ftext="Arial"  ;
ods graphics on / width=11 in height=8 in  height=1.2 outputfmt=png imagename='age cont AR_v2';
symbol1 i=sm10 color=black v=none h=.8 line=1 width=2;
symbol2 i=sm10 color=grey v=none line=16 h=.8;
symbol3 i=sm10 color=grey v=none line=16 h=.8;
legend1 label=none ;
axis1 label=(h=2 "Age") value=(h=2) order=(10 to 100 by 10) ;
axis2 label=(h=2 a=90 r=0 "Absolute Mortality Difference") 
value=(h=2) minor=none order=(-.07 to .01 by .01);
proc gplot data=morts;
   plot (diff lower upper)*age/ overlay
   haxis=axis1 
   vaxis=axis2 nolegend vref=0 cvref=red wvref=2;
   format age 6.0 diff lower upper percentn6.1;
run;


**********************************;

*Adjusted RR;
proc transpose data=predsa out=estimate;
where mort30_ed eq .;
by age;
var xbeta;
id Tx0;
run;

proc transpose data=predsa out=error;
where mort30_ed eq .;
by age;
var xbetase;
id Tx0;
run;

data morts;
merge estimate(rename=('0'n=control '1'n=case))
error(rename=('0'n=controlse '1'n=casese));
by age;
diff=case-control;
diff_se=sqrt(controlse**2+casese**2);
low=diff-1.96*diff_se;
up=diff+1.96*diff_se;
ratio=exp(diff);
lower=exp(low);
upper=exp(up);
run;

*RR Figure;
goptions reset=all ftext="Arial/bold" ;
ods graphics on / width=11 in height=8 in  height=1.2 outputfmt=png imagename='age cont RR_v2';
symbol1 i=sm10 color=black v=none h=.8 line=1 width=2;
symbol2 i=sm10 color=grey v=none line=16 h=.8;
symbol3 i=sm10 color=grey v=none line=16 h=.8;
legend1 label=none ;
axis1 label=(h=2 "Age") value=(h=2) order=(10 to 100 by 10);
axis2 label=(h=2 a=90 r=0 "Relative Risk") 
value=(h=2) minor=none order=(.6 to 1.1 by .1);
proc gplot data=morts;
   plot (ratio lower upper)*age/ overlay
   haxis=axis1
   vaxis=axis2 nolegend vref=1 cvref=red wvref=2;
   format age 6.0 ratio 6.1 ;
run;
