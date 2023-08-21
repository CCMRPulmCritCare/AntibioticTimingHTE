/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/
*
* Project: HTE Analysis
*
* Author: Patricia Kipnis
*
* Description: Translation of Excel Sheet to get RD, RDL, and RDU from:
*			   Newcombe, Robert G. "MOVER-R confidence intervals for ratios and products of two independently estimated quantities." 
*		   	   Statistical methods in medical research 25.5 (2016): 1774-1778.
*
/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/**/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/;

%macro RDCI(Data, Level, BR, BR_LOW, BR_HI, RR, RR_LOW, RR_HI);
data &Data;
Level = &Level;
t1=&RR.-1;
t2=1/&BR.;
RD=t1/t2;
L1=&RR_LOW.-1;
U1=&RR_HI.-1;
U2=1/&BR_LOW.; 
L2=1/&BR_HI.;
r1=(t1*t2-sqrt((t1*t2)**2-(L1*U2)*(2*t1-L1)*(2*t2-U2)))/(U2*(2*t2-U2));
r2=(t1*t2+sqrt((t1*t2)**2-(U1*L2)*(2*t1-U1)*(2*t2-L2)))/(L2*(2*t2-L2));
r3=(t1*t2-sqrt((t1*t2)**2-(L1*L2)*(2*t1-L1)*(2*t2-L2)))/(L2*(2*t2-L2));
r4=(t1*t2+sqrt((t1*t2)**2-(U1*U2)*(2*t1-U1)*(2*t2-U2)))/(U2*(2*t2-U2));
if U2 eq 2*t2 then do;
  r1=(L1*(1-(L1/(2*t1))))/t2;
  r4=(U1*(1-(U1/(2*t1))))/t2;
end;
RDL=r3;
if L1>0 then RDL=r1;
RDU=r4;
if U1>=0 then RDU=r2;
if &RR.<0 then do;
put "RR must be positive";
RD=.;
RL=.;
RU=.;
end;
if &BR.<0 or &BR.>1 then do;
put "BR must be between 0 and 1";
RD=.;
RL=.;
RU=.;
end;
if &RR_LOW.<0 or &RR_HI.<0 then do;
put "RR confidence levels must be positive";
RDL=.;
RDU=.;
end;
if &BR_LOW.<0 or &BR_HI.<0 then do;
put "BR confidence levels must be positive";
RDL=.;
RDU=.;
end;
if &BR_LOW.>=&BR. or &BR.>=&BR_HI. then do;
put "BR confidence levels are not well defined";
RDL=.;
RDU=.;
end;
if &RR_LOW.>=&RR. or &RR.>=&RR_HI. then do;
put "RR confidence levels are not well defined";
RDL=.;
RDU=.;
end;

keep Level RD RDL RDU;
run;

%mend;
