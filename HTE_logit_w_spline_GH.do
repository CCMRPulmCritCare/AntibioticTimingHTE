version 17
set more off
cap log close
clear all
set linesize 80

cd ""

local c_date = c(current_date)
*local c_time = c(current_time)
*local c_time_date = "`c_date'"+"_" +"`c_time'"
*local time_string = subinstr("`c_time_date'", ":", "_", .)
*local time_string = subinstr("`time_string'", " ", "_", .)
local date = subinstr("`c_date'", " ", "", .)
*display "`time_string'"

log using "Results\HTE R&R\Logs\HTE_RR_`date'.log", replace

***************************************************
* HTE - Logistic regression w/ spline TTA 
* Author: Sarah Seelye
* Date created: 11 Dec 2023
* Date updated: 20 Dec 2023
****************************************************

*-------------------------------------------------------------
* Logistic regression w/ spline variable for TTA
*-------------------------------------------------------------

use Data\HTE_RR\sepsis_risk_VA_KP_cancer_v4, clear

mkspline sp_hr = time_to_abx_hr, cubic nknots(3)
 
* metastatic cancer vs non-metastatic cancer
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney pressor_in_72hr aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new   							
		  
logit mort30_ed i.cancer_met##(c.sp_hr?) `covar' , or

preserve 

replace male=1
replace age=69
replace sirs_temp=1
replace sirs_rr=1 
replace sirs_pulse=1 
replace sirs_wbc=	1						
replace	aod_lactate= 1
replace aod_kidney= 1
replace pressor_in_72hr= 0
replace aod_liver= 0
replace aod_heme= 0
replace aod_lung=0	
replace	pulm= 0
replace chf= 0
replace dm_uncomp =0
replace dm_comp	=0
replace liver =0
replace neuro =	0				
replace renal =0
replace htn =1
replace cardic_arrhym =0
replace valvular_d2 =0
replace pulm_cir =0
replace pvd =0
replace paralysis =		0
replace pud =0
replace hypothyroid=0
replace lymphoma =0
replace ra =0
replace coag =0
replace obesity =0
replace wtloss =0
replace fen		=	1
replace	anemia_new =0

predict pr_mort, pr 
predict xb, xb 
predict error, stdp 

gen lb = xb - invnormal(0.975)*error
gen ub = xb + invnormal(0.975)*error
gen plb = invlogit(lb)
gen pub = invlogit(ub)

separate pr_mort, by(cancer_met) veryshortlabel
separate plb, by(cancer_met) veryshortlabel
separate pub, by(cancer_met) veryshortlabel 

twoway line pr_mort? time_to_abx_hr, sort 			///
						xla(0(3)12) xtick(0(3)12) 	///
						yla(0(0.1)0.3) ytick(0(0.1)0.3) ///
						xtitle("Time-to-antibiotics (hrs)") ///
						ytitle("Probability 30-day mortality") ///
						legend(label(1 "No metastatic" "cancer") ///
							   label(2 "Metastatic" "cancer") ///
							   symysize (6) size(small))	
						   
							   
sort time_to_abx_hr	
graph twoway (rarea plb1 pub1 time_to_abx_hr, fcolor(gs10) lwidth(none)) ///
				(line pr_mort1 time_to_abx_hr,			///
						xla(0(3)12) xtick(0(3)12) 	///
						yla(0(0.1)0.3) ytick(0(0.1)0.3) ///
						xtitle("Time-to-antibiotics (hrs)") ///
						ytitle("Probability 30-day mortality") ///
						lcolor(gs0))  	///
			 || (rarea plb0 pub0 time_to_abx_hr, fcolor(gs10) lwidth(none)) ///
					(line pr_mort0 time_to_abx_hr, lcolor(gs0) lpattern(dash))	///
			 || , legend(order(4 "No metastatic" "cancer" 2 "Metastatic" "cancer") ///
					 symysize (6) size(small))					

*graph save "Graph" "P:\ORD_Prescott_202208055D\HAPPI\HAPPI_KP\Jen\Figures\RR\Figure1a_cancer.gph"
*graph export "P:\ORD_Prescott_202208055D\HAPPI\HAPPI_KP\Jen\Figures\RR\Figure1a_cancer.png", as(png) name("Graph")					 
					 
restore


* Shock vs no shock
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney  aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new cancer_met cancer_nonmet  							
			   
logit mort30_ed i.pressor_in_72hr##(c.sp_hr?) `covar' , or

preserve 

replace male=1
replace age=69
replace sirs_temp=1
replace sirs_rr=1 
replace sirs_pulse=1 
replace sirs_wbc=	1						
replace	aod_lactate= 1
replace aod_kidney= 1
replace aod_liver= 0
replace aod_heme= 0
replace aod_lung=0	
replace	pulm= 0
replace chf= 0
replace dm_uncomp =0
replace dm_comp	=0
replace liver =0
replace neuro =	0				
replace renal =0
replace htn =1
replace cardic_arrhym =0
replace valvular_d2 =0
replace pulm_cir =0
replace pvd =0
replace paralysis =		0
replace pud =0
replace hypothyroid=0
replace lymphoma =0
replace ra =0
replace coag =0
replace obesity =0
replace wtloss =0
replace fen		=	1
replace	anemia_new =0
replace cancer_met =0
replace cancer_nonmet=0 

predict pr_mort, pr 
predict xb, xb 
predict error, stdp 

gen lb = xb - invnormal(0.975)*error
gen ub = xb + invnormal(0.975)*error
gen plb = invlogit(lb)
gen pub = invlogit(ub)

separate pr_mort, by(pressor_in_72hr) veryshortlabel
separate plb, by(pressor_in_72hr) veryshortlabel
separate pub, by(pressor_in_72hr) veryshortlabel 

twoway line pr_mort? time_to_abx_hr, sort ///
						xla(0(3)12) xtick(0(3)12) 	///
						yla(0(0.1)0.3) ytick(0(0.1)0.3) ///
						xtitle("Time-to-antibiotics (hrs)") ///
						ytitle("Probability 30-day mortality") ///
						legend(label(1 "No shock") ///
							   label(2 "Shock") ///
							   size(small))	

sort time_to_abx_hr	
graph twoway (rarea plb1 pub1 time_to_abx_hr, fcolor(gs10) lwidth(none)) ///
				(line pr_mort1 time_to_abx_hr,			///
						xla(0(3)12) xtick(0(3)12) 	///
						yla(0(0.1)0.3) ytick(0(0.1)0.3) ///
						xtitle("Time-to-antibiotics (hrs)") ///
						ytitle("Probability 30-day mortality") ///
						lcolor(gs0))  	///
			 || (rarea plb0 pub0 time_to_abx_hr, fcolor(gs10) lwidth(none)) ///
					(line pr_mort0 time_to_abx_hr, lcolor(gs0) lpattern(dash))	///
			 || , legend(order(4 "No shock" 2 "Shock") ///
					 symysize (6) size(small))								   

*graph save "Graph" "P:\ORD_Prescott_202208055D\HAPPI\HAPPI_KP\Jen\Figures\RR\Figure1b_shock.gph"
*graph export "P:\ORD_Prescott_202208055D\HAPPI\HAPPI_KP\Jen\Figures\RR\Figure1b_shock.png", as(png) name("Graph")					 
					 
restore

			  
*-------------------------------------------------------------
* spline variable for time to abx in minutes 
*-------------------------------------------------------------

foreach var in male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney pressor_in_72hr aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new   {
	sum `var'			
			}

mkspline sp_min = time_to_abx_min, cubic nknots(3)

 
* metastatic cancer vs non-metastatic cancer
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney pressor_in_72hr aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new   							
		  
logit mort30_ed i.cancer_met##(c.sp_min?) `covar' , or

preserve 

replace male=1
replace age=69
replace sirs_temp=1
replace sirs_rr=1 
replace sirs_pulse=1 
replace sirs_wbc=	1						
replace	aod_lactate= 1
replace aod_kidney= 1
replace pressor_in_72hr= 0
replace aod_liver= 0
replace aod_heme= 0
replace aod_lung=0	
replace	pulm= 0
replace chf= 0
replace dm_uncomp =0
replace dm_comp	=0
replace liver =0
replace neuro =	0				
replace renal =0
replace htn =1
replace cardic_arrhym =0
replace valvular_d2 =0
replace pulm_cir =0
replace pvd =0
replace paralysis =		0
replace pud =0
replace hypothyroid=0
replace lymphoma =0
replace ra =0
replace coag =0
replace obesity =0
replace wtloss =0
replace fen		=	1
replace	anemia_new =0

predict pr_mort, pr 

separate pr_mort, by(cancer_met) veryshortlabel

twoway line pr_mort? time_to_abx_min, sort 

restore


* Shock vs no shock
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney  aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new cancer_met cancer_nonmet  							
			   
logit mort30_ed i.pressor_in_72hr##(c.sp_min?) `covar' , or

preserve 

replace male=1
replace age=69
replace sirs_temp=1
replace sirs_rr=1 
replace sirs_pulse=1 
replace sirs_wbc=	1						
replace	aod_lactate= 1
replace aod_kidney= 1
replace aod_liver= 0
replace aod_heme= 0
replace aod_lung=0	
replace	pulm= 0
replace chf= 0
replace dm_uncomp =0
replace dm_comp	=0
replace liver =0
replace neuro =	0				
replace renal =0
replace htn =1
replace cardic_arrhym =0
replace valvular_d2 =0
replace pulm_cir =0
replace pvd =0
replace paralysis =		0
replace pud =0
replace hypothyroid=0
replace lymphoma =0
replace ra =0
replace coag =0
replace obesity =0
replace wtloss =0
replace fen		=	1
replace	anemia_new =0
replace cancer_met =0
replace cancer_nonmet=0 

predict pr_mort, pr 

separate pr_mort, by(pressor_in_72hr) veryshortlabel

twoway line pr_mort? time_to_abx_min, sort 

restore


*----------------------------------------
* spline diagnostics - using nknot(3)	
*----------------------------------------

* cancer
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney pressor_in_72hr aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new   							
		  
logit mort30_ed i.cancer_met##(c.sp_hr?) `covar' , or
est store spline
testparm sp_hr1 sp_hr2
estat ic

local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney pressor_in_72hr aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new   	

logit mort30_ed i.cancer_met##c.time_to_abx_hr `covar' , or
est store nospline
estat ic

lrtest spline nospline


* shock
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new cancer_met cancer_nonmet 							
		  
logit mort30_ed i.pressor_in_72hr##(c.sp_hr?) `covar' , or
est store spline2
testparm sp_hr1 sp_hr2
estat ic

local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new cancer_met cancer_nonmet   	

logit mort30_ed i.pressor_in_72hr##c.time_to_abx_hr `covar' , or
est store nospline2
estat ic

lrtest spline2 nospline2

*---------------------------
* Sensitivity Analysis 
*---------------------------

* examine relationship when TTA is categorical

sum time_to_abx_hr, de
version 16: table Tx0, c(min time_to_abx_hr max time_to_abx_hr)

* create a categorical variable for TTA
gen tta_hr_cat = .
replace tta_hr_cat = 0 if time_to_abx_hr<1
replace tta_hr_cat = 1 if time_to_abx_hr>=1 & time_to_abx_hr<2
replace tta_hr_cat = 2 if time_to_abx_hr>=2 & time_to_abx_hr<3
replace tta_hr_cat = 3 if time_to_abx_hr>=3 & time_to_abx_hr<4
replace tta_hr_cat = 4 if time_to_abx_hr>=4 & time_to_abx_hr<5
replace tta_hr_cat = 5 if time_to_abx_hr>=5 & time_to_abx_hr<6
replace tta_hr_cat = 6 if time_to_abx_hr>=6 & time_to_abx_hr<7
replace tta_hr_cat = 7 if time_to_abx_hr>=7 & time_to_abx_hr<8
replace tta_hr_cat = 8 if time_to_abx_hr>=8 & time_to_abx_hr<9
replace tta_hr_cat = 9 if time_to_abx_hr>=9 & time_to_abx_hr<10
replace tta_hr_cat = 10 if time_to_abx_hr>=10 & time_to_abx_hr<11
replace tta_hr_cat = 11 if time_to_abx_hr>=11 & time_to_abx_hr<=12

tab tta_hr_cat

* Metastatic Cancer vs no metastatic cancer (=no cancer + non-metastatic cancer)
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney pressor_in_72hr aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new   							
			   
logit mort30_ed i.cancer_met##i.tta_hr_cat `covar' , or
margins i.tta_hr_cat, over(cancer_met)
marginsplot , xtitle("Time-to-antibiotics (hrs)")  ///
			  ytitle("Probability of 30-day mortality") ylabel(0(0.1)0.4) ///
			  title("") ///
			  plot(, label("No metastatic cancer" "Metastatic cancer")) 			  


* Shock vs no shock
local covar male age														///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			aod_lactate aod_kidney  aod_liver aod_heme aod_lung	///
			pulm chf dm_uncomp dm_comp	liver neuro 						///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_new cancer_met cancer_nonmet  							
			   
logit mort30_ed i.pressor_in_72hr##i.tta_hr_cat `covar' , or
margins i.tta_hr_cat, over(pressor_in_72hr)
marginsplot , xtitle("Time-to-antibiotics (hrs)")  ///
			  ytitle("Probability of 30-day mortality") ylabel(0(0.1)0.4) ///
			  title("") ///
			  plot(, label("No shock" "Shock"))

log close