// C++ script derived from the original O'Hara-Rudy dynamic model
// including a Markov model for hERG gating similar with that described in
// Li et al. 2017 https://doi.org/10.1161/CIRCEP.116.004628
// but with simplified pharmacodynamic component including only blocking/unblocking rates
// used to compute several proarrhythmogenic risk predictors, particularly Qnet
// as described in Dutta et al. 2017 https://doi.org/10.3389/fphys.2017.00616
// pharmacology data for a 12-compounds set (CiPA training set)
// and a supplementary 16-compounds set (CiPA validation set)
// and for chloroquine and hydroxychloroquine using experimental data from:
// Thomet U, Amuzescu B, Knott T, Mann SA, Mubagwa K, Radu BM (2021)
// Assessment of Proarrhythmogenic Risk for Chloroquine and Hydroxychloroquine Using the CiPA Concept
// European Journal of Pharmacology (in review)



// Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy,
//                            Washington University in St. Louis.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the names of the copyright holders nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//

// C++ Implementation of the O'Hara-Rudy dynamic (ORd) model for the
// undiseased human ventricular action potential and calcium transient
//
// The ORd model is described in the article "Simulation of the Undiseased
// Human Cardiac Ventricular Action Potential: Model Formulation and
// Experimental Validation"
// by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
//
// The article and supplemental materails are freely available in the
// Open Access jounal PLoS Computational Biology
// Link to Article:
// http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
//
// Email: tom.ohara@gmail.com / rudy@wustl.edu
// Web: http://rudylab.wustl.edu
//

#include <math.h>
#include <iostream>
using namespace std;
#include <fstream>
using std::ifstream;
#include <string>


void revpots();//compute reversal potentials
void RGC();//compute rates, gates, and currents
void stimulus();//determine the value for the periodic stimulus
void voltage();//calculate the new membrane voltage
void dVdt_APD_CaD();//calculate voltage derivative and APD90
void FBC();//calculate fluxes, buffers, and concentrations

const double CL=2000;//pacing cycle length - orig CL=1000
const int beats=40;// Number of Beats
const double ft=beats*CL;//orig 3*CL final time - orig ft=1000*CL
const int skip=200;//orig skip=20 number of timesetps to skip in sampling of data in output file - orig skip=10
const double safetime=25.0;//time from the beginning of each beat during which dt is fixed to small values - orig safetime=25.0
const double beatssave=6;//orig 3  number of beats to save in the output - orig betassave=2

//const double amp=-80.0;//orig amp=-80  stimulus amplitude in uA/uF
const double start=500;//orig 500 start time of the stimulus, relative to each beat - orig start=0
const double duration=0.5;//duration of the stimulus in ms

const int celltype=0;  //endo = 0, epi = 1, M = 2

//initial values for state variables, there are 41 of them
double v=-87.5; //orig v=-87.5;
double nai=7;
double nass=nai;
double ki=145;
double kss=ki;
double cai=1.0e-4;
double cass=cai;
double cansr=1.2;
double cajsr=cansr;
double m=0;
double hf=1;
double hs=1;
double j=1;
double hsp=1;
double jp=1;
double mL=0;
double hL=1;
double hLp=1;
double a=0;
double iF=1;
double iS=1;
double ap=0;
double iFp=1;
double iSp=1;
double d=0;
double ff=1;
double fs=1;
double fcaf=1;
double fcas=1;
double jca=1;
double nca=0;
double ffp=1;
double fcafp=1;
double xrf=0;
double xrs=0;
double xs1=0;
double xs2=0;
double xk1=1;
double Jrelnp=0;
double Jrelp=0;
double CaMKt=0;
/* Hyperpolarization-activated (Funny) Current*/
	void comp_statesIf();  // Function that will compute state probabilities of If
	double ifna;  // If Na Current (uA/uF)
	double ifk;   // If K Current (uA/uF)
	double iftot; // If Total Current (uA/uF)
	//double scgf;  // Scaling Factor for Max conductance of If
	double gf=0.05;    // Max conductance of If (mS/uF)
	double gfna=0.3833;  // Na conductance of If (proportion)
	double gfk=0.6167;   // K conductance of If (proportion)
	double ach=0.0;  // ACh concentration (mM)
	double dvy;   // ACh-induced Voltage Shift in If V1/2
	double af;    // If alpha rate (ms^-1)
	double bf;    // If beta rate (ms^-1)
	const double kf = 0.92; //If kappa rate constant
	const double lf = 10;   //If lambda rate constant
	const double fff = 10.9; //If f rate constant
	double c1f=0.682;   // If probability of state C1
	double c2f=0.231;   // If probability of state C2
	double o1f=0.068;   // If probability of state O1
	double o2f=0.019;   // If probability of state O2

//constants
//double const nao=144.5;//extracellular sodium in mM
//double const cao=1.8;//extracellular calcium in mM
//double const ko=2.5;//extracellular potassium in mM
double const nao=140.0;//extracellular sodium in mM
double const cao=1.8;//extracellular calcium in mM
double const ko=5.4;//extracellular potassium in mM

//buffer paramaters
double const BSRmax=0.047;
double const KmBSR=0.00087;
double const BSLmax=1.124;
double const KmBSL=0.0087;
double const cmdnmax=0.05;
double const kmcmdn=0.00238;
double const trpnmax=0.07;
double const kmtrpn=0.0005;
double const csqnmax=10.0;
double const kmcsqn=0.8;

//CaMK paramaters
double const aCaMK=0.05;
double const bCaMK=0.00068;
double const CaMKo=0.05;
double const KmCaM=0.0015;
double const KmCaMK=0.15;

//physical constants
double const R=8314.0;
double const T=310.0;//default T=310.0  PT=310.0 (physiological temp. 37'C)  RT=298.0 (room temp. 25'C)
double const F=96485.0;

//cell geometry
//double const L=0.01;  //orig definition
double const rad=0.0011;  //orig definition
//double const vcell=1000*3.14*rad*rad*L;  //orig definition
//double const Ageo=2*3.14*rad*rad+2*3.14*rad*L;
//double const Acap=2*Ageo;  //orig definition, resulting in Acap=0.000153358 cm^2
double const Acap=153.358*0.000001;  //Cm (pF), orig Cm=153.358 pF
double const amp=-80;//amp=-2000/(Acap/0.000001); //scaled stimulus amplitude equivalent to 2000 pA total
double const rgc=2;//1.2;  //capacitive/geometric area ratio, orig rgc=2 for mature ventricular cm, for hiPSC-cm use rgc=1.2
double const Ageo=Acap/rgc;
double const L=(Ageo-2*3.14*rad*rad)/(2*3.14*rad);
double const vcell = 1000*3.14*rad*rad*L;
double const rvmyo=0.68;  //orig rvmyo=0.68
double const rvmito=0.26;  //orig rvmito=0.26
double const rvsr=0.06;  //orig rvsr=0.06
double const rvnsr=0.0552;  //orig rvnsr=0.0552
double const rvjsr=0.0048;  //orig rvjsr=0.0048
double const rvss=0.02;  //orig rvss=0.02
double const vmyo=rvmyo*vcell;  //orig vmyo=0.68*vcell
double const vmito=rvmito*vcell;  //orig vmito=0.26*vcell
double const vsr=rvsr*vcell;  //orig vsr=0.06*vcell
double const vnsr=rvnsr*vcell;  //orig vnsr=0.0552*vcell
double const vjsr=rvjsr*vcell;  //orig vjsr=0.0048*vcell
double const vss=rvss*vcell;  //orig vss=0.02*vcell

//introduce varaibles for reversal potentials, currents, fluxes, and CaMK
double ENa,EK,EKs;
double INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist,If; //orig without If
double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
double CaMKa,CaMKb;

//variables and parameters for IKr dynamic (Markov) model (from Li,Dutta,Colatsky 2017)
void comp_statesIKr();  // Function that will compute state probabilities of IKr
double IKrB1=0.00004631;
double IKrB2=-0.004226;
double IKrA2=0.000004986;
double IKrA1=0.0264;
double IKrA3=0.001214;
double IKrB3=0.008516;
double IKrA4=0.00001854;
double IKrB4=-0.04641;
double IKrA6=0.03776;
double IKrB6=-0.00000006304;
double IKrA5=0.000000007057;
double IKrB5=0.000000009502;
//drug binding&unbinding parameters
double IKrKf=0.05;
double IKrKf2=IKrKf;
double IKrKt=0.000035;
double IKrKt2=0.000035;
double IKrVhalf=1;
double IKrKu=0.01;
double IKrnH=1;
double IKrhalfmax=1;
double IKrunnamed=1;
//double IKrtrap=1;
double IKrA11=0.0007868;
double IKrB11=0.00000001535;
double IKrq11=4.942;
double IKrA21=0.000005455;
double IKrB21=-0.1688;
double IKrq21=4.156;
double IKrq1=4.843;
double IKrq2=4.23;
double IKrq31=4.22;
double IKrq41=1.459;
double IKrA31=0.005509;
double IKrB31=0.000000007771;
double IKrA41=0.001416;
double IKrB41=-0.02877;
double IKrq3=4.962;
double IKrq4=3.769;
double IKrq51=5;
double IKrq61=5.568;
double IKrA51=0.4492;
double IKrB51=0.008595;
double IKrA61=0.01241;
double IKrB61=0.1725;
double IKrq53=2.412;
double IKrq63=5.682;
double IKrA53=0.149;
double IKrB53=0.004668;
double IKrA63=0.008978;
double IKrB63=-0.02215;
double IKrq52=4.663;
double IKrq62=5;
double IKrA52=0.3181;
double IKrB52=0.00000003613;
double IKrA62=0.3226;
double IKrB62=-0.0006575;

double IKrIC1=1;
double IKrIC2=0;
double IKrC1=0;
double IKrC2=0;
double IKrO=0;
double IKrIO=0;
//double IKrIObound=0;
//double IKrObound=0;
double IKrCbound=0;
//double dIKrIC1dt;
//double dIKrIC2dt;
//double dIKrC1dt;
//double dIKrC2dt;
//double dIKrOdt;
//double dIKrIOdt;
//double dIKrIObounddt;
//double dIKrObounddt;
//double dIKrCbounddt;
double ReactionFlux1;
double ReactionFlux2;
double ReactionFlux3;
double ReactionFlux4;
double ReactionFlux5;
double ReactionFlux6;
double ReactionFlux7;
double ReactionFlux8;
double ReactionFlux9;
//double ReactionFlux10;
//double ReactionFlux11;
double IKrkb; //new hERG1 blocking rate from onset-of-block kinetics analysis
double IKrkub; //new hERG1 unblocking rate from onset-of-block kinetics analysis


//introduce APD, timing, and counting parameters
int APD50_flag=0;
int APD90_flag=0;
int CaD50_flag=0;
int CaD90_flag=0;
int vdotmax_flag=0;
int caidotmax_flag=0;
int caimax_flag=0;
double t_vdot_max[beats+1];
double t_caidot_max[beats+1];
double vo=v;
double caio=cai;
double dt=0.005;
double t0=0;
double t=0;
double dto;
double vdot_old;
double vdot=0;
double caidot_old;
double caidot=0;
int p=1;
int n=0;
int count1=1;

//introduce proarrhythmogenesis prediction metrics
double Iqinward;
double Inet;
double cqInward[beats];
double Qnet[beats];
double caimax[beats];
double caidot_max[beats];
double cairest[beats];
double CaD50[beats];
double CaD90[beats];
double APD50[beats];
double APD90[beats];
double vmax[beats];
double vdot_max[beats];
double vrest[beats];
int EAD_flag[beats+1];
double APD;

//scaling factors for ion conductance densities
double scgna=1.0; //0.00375;  //orig scgna=1.0
double scgnal=scgna*2.661;  //orig scgnal=scgna //scgnal=2.661 in Dutta et al. 2017
double scgcal=1.007; //0.774;  //orig scgcal=1.0 //1.007 in Dutta et al. 2017
double scgkr=1.013;  //orig scgkr=1.0 //1.013 in Dutta et al. 2017
double scgks=1.87;  //orig scgks=1.0 //1.87 in Dutta et al. 2017
double scgki=1.698;  //orig scgki=1.0 //1.698 in Dutta et al 2017
double scgtof=1.0;  //orig scgtof=1.0
double scgtos=1.0;  //orig scgtos=1.0
double scgf=0.0; //0.124;  //orig scgf=1.0
double scgncx=1.0;  //orig scgncx=1.0
double scgnak=1.0;  //orig scgnak=1.0
double scgpca=1.0;  //orig scgpca=1.0
double scgrel=1.0;  //orig scgrel=1.0 [total Ca2+ release via ryanodine receptors from jsr to myoplasm]
double scgup=1.0;  //orig scgup=1.0 [total Ca2+ uptake via SERCA pumps from myoplasm to nsr]
double scgkb=1.0;  //orig scgkb=1.0
double scgnab=1.0;  //orig scgnab=1.0
double scgcab=1.0;  //orig scgcab=1.0

//value holders for state varaibles in the case that the increase in dt was too aggressive, so a smaller one can be taken
double nai0,nass0,ki0,kss0,cai0,cass0,cansr0,cajsr0,m0,hf0,hs0,jO,hsp0,jp0,mL0,hL0,hLp0,a0,iF0,iS0,ap0,iFp0,iSp0,d0,ff0,fs0,fcaf0,fcas0,jca0,nca0,ffp0,fcafp0,xrf0,xrs0,xs10,xs20,xk10,Jrelnp0,Jrelp0,CaMKt0;


//pharmacodynamic data
int nCmax;
string cpd;
string IKrmode;
double Cmax;
double IKrIC50;
double IKrnHill;
double INaLIC50;
double INaLnHill;
double ICaLIC50;
double ICaLnHill;
double INaIC50;
double INanHill;
double ItoIC50;
double ItonHill;
double IK1IC50;
double IK1nHill;
double IKsIC50;
double IKsnHill;
double actshift;


int main()
{
string cpd;
cpd="chloroquine";
if (cpd=="CTRL")
	{Cmax=0;IKrIC50=1;IKrnHill=1;INaLIC50=1;INaLnHill=1;ICaLIC50=1;ICaLnHill=1;INaIC50=1;INanHill=1;ItoIC50=1;ItonHill=1;IK1IC50=1;IK1nHill=1;IKsIC50=1;IKsnHill=1;IKrKf=0;IKrKu=0;IKrhalfmax=100000000;IKrnH=1;IKrVhalf=0;actshift=0;}
//CiPA training set
//values from Li et al. 2017 https://doi.org/10.1161/CIRCEP.116.004628
//and Dutta et al. 2017 https://doi.org/10.3389/fphys.2017.00616
//low-risk
if (cpd=="mexiletine")
	{Cmax=4129;IKrIC50=28880;IKrnHill=0.9;INaLIC50=8956.8;INaLnHill=1.4;ICaLIC50=38243.6;ICaLnHill=1;INaIC50=100000000;INanHill=1;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=9.996;IKrKu=0.09967;IKrhalfmax=2308000;IKrnH=1.304;IKrVhalf=-86.26;actshift=0;}
if (cpd=="diltiazem")
	{Cmax=122;IKrIC50=13150;IKrnHill=0.9;INaLIC50=21868.5;INaLnHill=0.7;ICaLIC50=112.1;ICaLnHill=0.7;INaIC50=110859;INanHill=0.7;ItoIC50=2820000000;ItonHill=0.2;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=251;IKrKu=0.2816;IKrhalfmax=1000000;IKrnH=0.9485;IKrVhalf=-90.89;actshift=0;}
if (cpd=="ranolazine")
	{Cmax=1948.2;IKrIC50=8270;IKrnHill=0.9;INaLIC50=7884.5;INaLnHill=0.9;ICaLIC50=100000000;ICaLnHill=1;INaIC50=68774;INanHill=1.4;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=36155020;IKsnHill=0.5;IKrKf=55.84;IKrKu=0.01929;IKrhalfmax=147200;IKrnH=0.95;IKrVhalf=-94.87;actshift=0;}
if (cpd=="verapamil")
	{Cmax=81;IKrIC50=288;IKrnHill=1;INaLIC50=7028;INaLnHill=1;ICaLIC50=201.8;ICaLnHill=1.1;INaIC50=100000000;INanHill=1;ItoIC50=13429.2;ItonHill=0.8;IK1IC50=349000000;IK1nHill=0.3;IKsIC50=100000000;IKsnHill=1;IKrKf=46460;IKrKu=0.0007927;IKrhalfmax=9184000;IKrnH=1.043;IKrVhalf=-100;actshift=0;}
//intermediate-risk
if (cpd=="chlorpromazine")
	{Cmax=38;IKrIC50=929.2;IKrnHill=0.8;INaLIC50=4559.6;INaLnHill=0.9;ICaLIC50=8191.9;ICaLnHill=0.8;INaIC50=4535.6;INanHill=2;ItoIC50=17616711;ItonHill=0.4;IK1IC50=9269.9;IK1nHill=0.7;IKsIC50=100000000;IKsnHill=1;IKrKf=206000;IKrKu=0.03866;IKrhalfmax=56770000;IKrnH=0.8871;IKrVhalf=-14.57;actshift=0;}
if (cpd=="terfenadine")
	{Cmax=4;IKrIC50=23;IKrnHill=0.6;INaLIC50=20056;INaLnHill=0.6;ICaLIC50=700.4;ICaLnHill=0.7;INaIC50=4803.2;INanHill=1;ItoIC50=239960.8;ItonHill=0.3;IK1IC50=100000000;IK1nHill=1;IKsIC50=399754;IKsnHill=0.5;IKrKf=9884;IKrKu=0.0000818;IKrhalfmax=41380;IKrnH=0.65;IKrVhalf=-77.49;actshift=0;}
if (cpd=="ondansetron")
	{Cmax=139;IKrIC50=1320;IKrnHill=0.9;INaLIC50=19180.8;INaLnHill=1;ICaLIC50=22551.4;ICaLnHill=0.8;INaIC50=57666.4;INanHill=1;ItoIC50=1023378;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=569807;IKsnHill=0.7;IKrKf=33540;IKrKu=0.02325;IKrhalfmax=9950000;IKrnH=0.8874;IKrVhalf=-82.11;actshift=0;}
if (cpd=="cisapride")
	{Cmax=2.6;IKrIC50=10.1;IKrnHill=0.7;INaLIC50=100000000;INaLnHill=1;ICaLIC50=9258076;ICaLnHill=0.4;INaIC50=100000000;INanHill=1;ItoIC50=219112.4;ItonHill=0.2;IK1IC50=29498;IK1nHill=0.5;IKsIC50=81192862;IKsnHill=0.3;IKrKf=9.997;IKrKu=0.0004161;IKrhalfmax=42.06;IKrnH=0.9728;IKrVhalf=-199.5;actshift=0;}
//high-risk
if (cpd=="sotalol")
	{Cmax=14690;IKrIC50=110600;IKrnHill=0.8;INaLIC50=100000000;INaLnHill=1;ICaLIC50=7061527;ICaLnHill=0.9;INaIC50=1140000000;INanHill=0.5;ItoIC50=43143455;ItonHill=0.7;IK1IC50=3050260;IK1nHill=1.2;IKsIC50=4221856;IKsnHill=1.2;IKrKf=2403;IKrKu=0.01985;IKrhalfmax=9619000;IKrnH=0.7516;IKrVhalf=-55;actshift=0;}
if (cpd=="bepridil")
	{Cmax=33;IKrIC50=50;IKrnHill=0.9;INaLIC50=1813.9;INaLnHill=1.4;ICaLIC50=2808.1;ICaLnHill=0.6;INaIC50=2929.3;INanHill=1.2;ItoIC50=8594;ItonHill=3.5;IK1IC50=100000000;IK1nHill=1;IKsIC50=28628.3;IKsnHill=1.4;IKrKf=37350000;IKrKu=0.0001765;IKrhalfmax=1000000000;IKrnH=0.9365;IKrVhalf=-54.93;actshift=0;}
if (cpd=="dofetilide")
	{Cmax=2;IKrIC50=4.9;IKrnHill=0.9;INaLIC50=753160.4;INaLnHill=0.3;ICaLIC50=260.3;ICaLnHill=1.2;INaIC50=380.5;INanHill=0.9;ItoIC50=18.8;ItonHill=0.8;IK1IC50=394.3;IK1nHill=0.8;IKsIC50=100000000;IKsnHill=1;IKrKf=100000000;IKrKu=0.0000179;IKrhalfmax=548300000;IKrnH=0.9999;IKrVhalf=-1.147;actshift=0;}
if (cpd=="quinidine")
	{Cmax=3237;IKrIC50=992;IKrnHill=0.8;INaLIC50=9417;INaLnHill=1.3;ICaLIC50=51592.3;ICaLnHill=0.6;INaIC50=12329;INanHill=1.5;ItoIC50=3487.4;ItonHill=1.3;IK1IC50=39589919;IK1nHill=0.4;IKsIC50=4898.9;IKsnHill=1.4;IKrKf=5770;IKrKu=0.01;IKrhalfmax=1000000;IKrnH=0.8311;IKrVhalf=-64.87;actshift=0;}
//CiPA validation set
//described in Han et al. 2020 https://doi.org/10.1016/j.vascn.2020.106890
//and Li et al. 2020 https://doi.org/10.1002/cpt.1647
//alternate data from Yamazaki et al 2018 https://doi.org/10.1016/j.jphs.2018.02.005
//general IC50 & nHill data from: https://github.com/FDA/CiPA/blob/Lab_Specific_Validation_Calibration_2020/chantest_Hill_fitting/results/drug_name/optimal.csv
//hERG kinetics data from: https://github.com/FDA/CiPA/blob/Lab_Specific_Validation_Calibration_2020/hERG_fitting/results/drug_name/pars.txt
//Cmax data from https://github.com/FDA/CiPA/blob/Lab_Specific_Validation_Calibration_2020/chantest_AP_simulation/data/newCiPA.csv
if (cpd=="astemizole")
	{Cmax=0.26;IKrIC50=27.01;IKrnHill=1.504;INaLIC50=786.5;INaLnHill=1.98;ICaLIC50=972;ICaLnHill=1.98;INaIC50=1862;INanHill=1.98;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=77530;IKrKu=0.00003073;IKrhalfmax=169600;IKrnH=0.6762;IKrVhalf=-1.001;actshift=0;}
//Yamazaki 2018 fETPC 0.51
if (cpd=="azimilide")
	{Cmax=70;IKrIC50=90.79;IKrnHill=0.505;INaLIC50=3700;INaLnHill=1.352;ICaLIC50=9587;ICaLnHill=1.98;INaIC50=18970;INanHill=1.98;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=78990;IKrKu=0.007555;IKrhalfmax=1296000;IKrnH=0.5785;IKrVhalf=-10.16;actshift=0;}
//Yamazaki 2018 fETPC 72
if (cpd=="clarithromycin")
	{Cmax=1206;IKrIC50=21180;IKrnHill=0.8009;INaLIC50=269100;INaLnHill=1.335;ICaLIC50=128900;ICaLnHill=1.98;INaIC50=243300;INanHill=1.98;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=68.65;IKrKu=0.01517;IKrhalfmax=56880;IKrnH=0.7455;IKrVhalf=-161.7;actshift=0;}
if (cpd=="clozapine")
	{Cmax=71;IKrIC50=980.6;IKrnHill=1.98;INaLIC50=2235;INaLnHill=1.347;ICaLIC50=4477;ICaLnHill=1.98;INaIC50=9338;INanHill=1.889;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=6.331;IKrKu=0.02552;IKrhalfmax=33760;IKrnH=1.314;IKrVhalf=-8.325;actshift=0;}
//Yamazaki 2018 fETPC 320
if (cpd=="disopyramide")
	{Cmax=742;IKrIC50=906.2;IKrnHill=0.9325;INaLIC50=20050;INaLnHill=0.7218;ICaLIC50=110000;ICaLnHill=0.8718;INaIC50=73110;INanHill=0.9003;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=1.756;IKrKu=0.0364;IKrhalfmax=1560000;IKrnH=1.999;IKrVhalf=-72.32;actshift=0;}
if (cpd=="domperidone")
	{Cmax=19;IKrIC50=20.31;IKrnHill=1.378;INaLIC50=669;INaLnHill=0.8924;ICaLIC50=14120;ICaLnHill=1.291;INaIC50=6659;INanHill=1.425;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=677100;IKrKu=0.0002305;IKrhalfmax=3576000;IKrnH=0.6137;IKrVhalf=-24.19;actshift=0;}
//Yamazaki 2018 fETPC 27
if (cpd=="droperidol")
	{Cmax=6.33;IKrIC50=30.97;IKrnHill=0.505;INaLIC50=530.8;INaLnHill=1.017;ICaLIC50=7427;ICaLnHill=0.9114;INaIC50=8063;INanHill=1.067;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=15.85;IKrKu=0.001267;IKrhalfmax=90.09;IKrnH=0.5489;IKrVhalf=-78.78;actshift=0;}
if (cpd=="ibutilide")
	{Cmax=100;IKrIC50=2.02;IKrnHill=1.98;INaLIC50=811;INaLnHill=0.9446;ICaLIC50=31200;ICaLnHill=0.9494;INaIC50=6957;INanHill=0.6987;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=14.45;IKrKu=0.00006104;IKrhalfmax=35.07;IKrnH=0.902;IKrVhalf=-9.817;actshift=0;}
//Yamazaki 2018 fETPC 150
if (cpd=="loratadine")
	{Cmax=0.45;IKrIC50=915.7;IKrnHill=1.98;INaLIC50=2187;INaLnHill=1.087;ICaLIC50=2652;ICaLnHill=1.522;INaIC50=18180;INanHill=1.98;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=98960;IKrKu=0.01055;IKrhalfmax=106500000;IKrnH=0.8403;IKrVhalf=-1;actshift=0;}
//Yamazaki 2018 fETPC 0.6
if (cpd=="metoprolol")
	{Cmax=1800;IKrIC50=13550;IKrnHill=0.9451;INaLIC50=24350;INaLnHill=0.7559;ICaLIC50=391000;ICaLnHill=0.8924;INaIC50=148800;INanHill=1.004;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=82210;IKrKu=0.9785;IKrhalfmax=53800000;IKrnH=0.6873;IKrVhalf=-86.19;actshift=0;}
if (cpd=="nifedipine")
	{Cmax=7.7;IKrIC50=28770;IKrnHill=1.126;INaLIC50=781.5;INaLnHill=0.7852;ICaLIC50=33.89;ICaLnHill=0.8526;INaIC50=41720;INanHill=1.538;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=1.029;IKrKu=0.3565;IKrhalfmax=130200;IKrnH=0.9613;IKrVhalf=-23.05;actshift=0;}
if (cpd=="nitrendipine")
	{Cmax=3.02;IKrIC50=10840;IKrnHill=1.404;INaLIC50=356.2;INaLnHill=0.8922;ICaLIC50=527.1;ICaLnHill=0.6373;INaIC50=10410;INanHill=1.211;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=1.847;IKrKu=0.991;IKrhalfmax=15820000;IKrnH=1.561;IKrVhalf=-68.23;actshift=0;}
if (cpd=="pimozide")
	{Cmax=0.431;IKrIC50=47.02;IKrnHill=1.865;INaLIC50=6624;INaLnHill=0.5151;ICaLIC50=443.3;ICaLnHill=1.151;INaIC50=7943000000;INanHill=1.235;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=12.74;IKrKu=0.00004288;IKrhalfmax=5.314;IKrnH=0.57;IKrVhalf=-71.47;actshift=0;}
//Yamazaki 2018 fETPC 0.43
if (cpd=="risperidone")
	{Cmax=1.81;IKrIC50=78.25;IKrnHill=1.124;INaLIC50=3301;INaLnHill=1.98;ICaLIC50=8503;ICaLnHill=1.98;INaIC50=7943000000;INanHill=1.389;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=3.726;IKrKu=0.0009118;IKrhalfmax=93.16;IKrnH=0.8593;IKrVhalf=-72.15;actshift=0;}
//Yamazaki 2018 fETPC 1.5
if (cpd=="tamoxifen")
	{Cmax=21;IKrIC50=2115;IKrnHill=1.98;INaLIC50=4515;INaLnHill=1.98;ICaLIC50=3944;ICaLnHill=1.98;INaIC50=9227;INanHill=1.98;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=4.429;IKrKu=0.006673;IKrhalfmax=144900;IKrnH=1.715;IKrVhalf=-1;actshift=0;}
if (cpd=="vandetanib")
	{Cmax=255.4;IKrIC50=458.8;IKrnHill=0.8127;INaLIC50=4654;INaLnHill=1.812;ICaLIC50=3829;ICaLnHill=1.98;INaIC50=6335;INanHill=1.98;ItoIC50=100000000;ItonHill=1;IK1IC50=100000000;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=25.67;IKrKu=0.01765;IKrhalfmax=797.7;IKrnH=0.7128;IKrVhalf=-41.57;actshift=0;}

//new
if (cpd=="chloroquine")
	{Cmax=410;IKrIC50=1820;IKrnHill=1.36;INaLIC50=30400;INaLnHill=0.85;ICaLIC50=30700;ICaLnHill=1;INaIC50=159000;INanHill=0.98;ItoIC50=4600000;ItonHill=1;IK1IC50=5860;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=1;IKrKu=1;IKrhalfmax=1;IKrnH=1;IKrVhalf=0; IKrkb=0.000366; IKrkub=0.000108;actshift=0.0052;}
if (cpd=="hydroxychloroquine")
	{Cmax=495;IKrIC50=3420;IKrnHill=1.02;INaLIC50=64900;INaLnHill=0.89;ICaLIC50=90000;ICaLnHill=1;INaIC50=96200;INanHill=0.81;ItoIC50=4600000;ItonHill=1;IK1IC50=29280;IK1nHill=1;IKsIC50=100000000;IKsnHill=1;IKrKf=1;IKrKu=1;IKrhalfmax=1;IKrnH=1;IKrVhalf=0; IKrkb=0.000213; IKrkub=0.000566;actshift=0.001;}
//switch to include effects of chloroquine/hydroxychloroquine on hERG voltage-dependent activation
string IKrmode; //possible values: "vshift", "novshift"
IKrmode="vshift";

//establish the output file "output.txt"
FILE*output;
output = fopen("output.txt","w");
//fprintf(output,"%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\n","t","v","nai","nass","ki","kss","cai","cass","cansr","cajsr","Jrel","CaMKt","Jup","Jtr","Jdiff","JdiffNa","JdiffK","Jleak","INa","INaL","Ito","ICaL","ICaNa","ICaK","IKr","IKs","IK1","If","INaCa_i","INaCa_ss","INaCa","INaK","IKb","INab","IpCa","ICab","Ist","dt","APD");
fprintf(output,"%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\n","t","v","Qnet","cqInward","Inet","Iqinward","cai","cass","cansr","cajsr","Jrel","n Cmax");
fclose(output);

FILE*fmaxs;
fmaxs = fopen("fmaxs","w");
fprintf(fmaxs,"%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\t%-07s\%-07s\n","n Cmax","Qnet","cqInward","caimax","cairest","Catri","CaD50","CaD90","APDtri","APD50","APD90","vmax","vdot_max","vrest","EAD_flag");
fclose(fmaxs);

for (nCmax=1; nCmax<=25; nCmax++){
cout<<cpd<<" n Cmax = "<<nCmax<<endl;//output runtime progress to the screen
v=-87.5; //orig v=-87.5;
nai=7;
nass=nai;
ki=145;
kss=ki;
cai=1.0e-4;
cass=cai;
cansr=1.2;
cajsr=cansr;
APD50_flag=0;
APD90_flag=0;
CaD50_flag=0;
CaD90_flag=0;
vdotmax_flag=0;
caidotmax_flag=0;
caimax_flag=0;
vo=v;
caio=cai;
dt=0.005;
t0=0;
t=0;
vdot=0;
caidot=0;
p=1;
n=0;
count1=1;

while (t<=ft)
	{
	//rules for dynamic dt choice, and model integration, comment to use fixed time steps
/*	if ((t>=(start+n*CL-2) && t<(start+duration+n*CL)) || (n>=1 && t<(start+duration+(n-1)*CL+safetime)) || (APD_flag==1 && v<0.7*vrest))
		{
		dt=0.005;
		t=t+dt;
		revpots();
		RGC();
		stimulus();
		vo=v;
		voltage();
		dVdt_APD();
		FBC();
		}
	else if (fabs(v-vo)<0.2)
		{
		dt=fabs(0.8/vdot);
		if (dt>1.0)
			{
			dt=1.0;
			}
		t=t+dt;
		revpots();
		RGC();
		stimulus();
		vo=v;
		voltage();
		dVdt_APD();
		FBC();
		}
	else if (fabs(v-vo)>0.8)
		{
		nai0=nai;
		nass0=nass;
		ki0=ki;
		cai0=cai;
		cass0=cass;
		cansr0=cansr;
		cajsr0=cajsr;
		m0=m;
		hf0=hf;
		hs0=hs;
		jO=j;
		hsp0=hsp;
		jp0=jp;
		mL0=mL;
		hL0=hL;
		hLp0=hLp;
		a0=a;
		iF0=iF;
		iS0=iS;
		ap0=ap;
		iFp0=iFp;
		iSp0=iSp;
		d0=d;
		ff0=ff;
		fs0=fs;
		fcaf0=fcaf;
		fcas0=fcas;
		jca0=jca;
		nca0=nca;
		ffp0=ffp;
		fcafp=fcafp;
		xrf0=xrf;
		xrs0=xrs;
		xs10=xs1;
		xs20=xs2;
		xk10=xk1;
		Jrelnp0=Jrelnp;
		Jrelp0=Jrelp;
		CaMKt0=CaMKt;

		t0=t;
		dto=dt;
		dt=fabs(0.2/vdot);
		t=t+dt;
		revpots();
		RGC();
		stimulus();
		vo=v;
		voltage();
		dVdt_APD();
		FBC();
		while (fabs(v-vo)>0.8)
			{
			v=vo;
			nai=nai0;
			nass=nass0;
			ki=ki0;
			cai=cai0;
			cass=cass0;
			cansr=cansr0;
			cajsr=cajsr0;
			m=m0;
			hf=hf0;
			hs=hs0;
			j=jO;
			hsp=hsp0;
			jp=jp0;
			mL=mL0;
			hL=hL0;
			hLp=hLp0;
			a=a0;
			iF=iF0;
			iS=iS0;
			ap=ap0;
			iFp=iFp0;
			iSp=iSp0;
			d=d0;
			ff=ff0;
			fs=fs0;
			fcaf=fcaf0;
			fcas=fcas0;
			jca=jca0;
			nca=nca0;
			ffp=ffp0;
			fcafp=fcafp0;
			xrf=xrf0;
			xrs=xrs0;
			xs1=xs10;
			xs2=xs20;
			xk1=xk10;
			Jrelnp=Jrelnp0;
			Jrelp=Jrelp0;
			CaMKt=CaMKt0;
			if (p==1)
				{
				dt=dto-0.01;
				p=0;
				}
			else
				{
				dt=dt-0.01;
				}
			if (dt<=0)
				{
				dt=1e-6;
				}
			t=t0+dt;
			revpots();
			RGC();
			stimulus();
			voltage();
			dVdt_APD();
			FBC();
			}
		p=1;
		}
	else
		{
		t=t+dt;
		revpots();
		RGC();
		stimulus();
		vo=v;
		voltage();
		dVdt_APD();
		FBC();
		}*/

	//uncomment below, and comment above to use a fixed dt
	t=t+dt; //fixed time step
	revpots();
	RGC();
	stimulus();
	vo=v;
	voltage();
	caio=cai;
	FBC();
	dVdt_APD_CaD();



//	if (count1%500000==0)
//		{
//		cout<<t/ft*100<<"% complete"<<endl;//output runtime progress to the screen
//		}

	if ((nCmax==1||(nCmax)%5==0) && count1%skip==0 && t>=ft-beatssave*CL)//save results ot output file when the sampling interval and time are correct
		{
//		fprintf(output,"%-018e\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\n",
//		t-(ft-beatssave*CL),v,nai,nass,ki,kss,cai,cass,cansr,cajsr,Jrel,CaMKt,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak,INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,If,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist,dt,APD);
		output = fopen("output.txt","a");//open the output file for appending new data
		fprintf(output,"%-018e\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%i\n",
		t-(ft-beatssave*CL),v,Qnet[n],cqInward[n],Inet,Iqinward,cai,cass,cansr,cajsr,Jrel,nCmax);
		fclose(output);//close the output file
		}

	count1++;//increase the loop counter

	}//end of numeric integration



fmaxs = fopen("fmaxs","a");
	for (n = beats-10; n <= beats; n++)
		{
		fprintf(fmaxs,"%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%i\n",nCmax,Qnet[n],cqInward[n],caimax[n],cairest[n],CaD90[n]-CaD50[n],CaD50[n],CaD90[n],APD90[n]-APD50[n],APD50[n],APD90[n],vmax[n],vdot_max[n],vrest[n],EAD_flag[n]);
		}
fclose(fmaxs);

}//end of nCmax loops

return 0;
} //end of main loop

void revpots()
{
ENa=(R*T/F)*log(nao/nai);
EK=(R*T/F)*log(ko/ki);
EKs=(R*T/F)*log((ko+0.01833*nao)/(ki+0.01833*nai));
}

void RGC()
{
CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
CaMKa=CaMKb+CaMKt;
double vffrt=v*F*F/(R*T);
double vfrt=v*F/(R*T);

double mss=1.0/(1.0+exp((-(v+39.57))/9.871));
double tm=pow(1.98,-(T-310)/10)/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955)); //orig tm=1.0/...
m=mss-(mss-m)*exp(-dt/tm);
double hss=1.0/(1+exp((v+82.90)/6.086));
double thf=pow(1.19,-(T-310)/10)/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));  //orig thf=1.0/...
double ths=pow(1.19,-(T-310)/10)/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));  //orig ths=1.0/...
double Ahf=0.99;
double Ahs=1.0-Ahf;
hf=hss-(hss-hf)*exp(-dt/thf);
hs=hss-(hss-hs)*exp(-dt/ths);
double h=Ahf*hf+Ahs*hs;
double jss=hss;
double tj=pow(1.19,-(T-310)/10)*(2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45)));  //orig tj=...
j=jss-(jss-j)*exp(-dt/tj);
double hssp=1.0/(1+exp((v+89.1)/6.086));
double thsp=3.0*ths;
hsp=hssp-(hssp-hsp)*exp(-dt/thsp);
double hp=Ahf*hf+Ahs*hsp;
double tjp=1.46*tj;
jp=jss-(jss-jp)*exp(-dt/tjp);
double GNa=75*scgna;  //orig GNa=75
double fINap=(1.0/(1.0+KmCaMK/CaMKa));
INa=(1/(1+pow(((nCmax*Cmax)/INaIC50),INanHill)))*GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);

double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
double tmL=tm;
mL=mLss-(mLss-mL)*exp(-dt/tmL);
double hLss=1.0/(1.0+exp((v+87.61)/7.488));
double thL=pow(1.19,-(T-310)/10)*200.0;  //orig thL=200.0
hL=hLss-(hLss-hL)*exp(-dt/thL);
double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
double thLp=3.0*thL;
hLp=hLssp-(hLssp-hLp)*exp(-dt/thLp);
double GNaL=0.0075*scgnal;  //orig GNaL=0.0075
if (celltype==1)
{
GNaL*=0.6;
}
double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
INaL=(1/(1+pow(((nCmax*Cmax)/INaLIC50),INaLnHill)))*GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
double ta=pow(1.19,(T-310)/10)*1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));  //orig ta=1.0515/...
a=ass-(ass-a)*exp(-dt/ta);
double iss=1.0/(1.0+exp((v+43.94)/5.711));
double delta_epi;
if (celltype==1)
{
delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
}
else
{
delta_epi=1.0;
}
double tiF=pow(1.72,-(T-310)/10)*(4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59)));  //orig tiF=4.562+...
double tiS=pow(1.72,-(T-310)/10)*(23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079)));  //orig tiS=23.62+...
tiF*=delta_epi;
tiS*=delta_epi;
double AiF=1.0/(1.0+exp((v-213.6)/151.2));
double AiS=1.0-AiF;
iF=iss-(iss-iF)*exp(-dt/tiF);
iS=iss-(iss-iS)*exp(-dt/tiS);
double i=AiF*iF+AiS*iS;
double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
ap=assp-(assp-ap)*exp(-dt/ta);
double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
double tiFp=dti_develop*dti_recover*tiF;
double tiSp=dti_develop*dti_recover*tiS;
iFp=iss-(iss-iFp)*exp(-dt/tiFp);
iSp=iss-(iss-iSp)*exp(-dt/tiSp);
double ip=AiF*iFp+AiS*iSp;
double Gto=0.02*scgtof;  //orig Gto=0.02
if (celltype==1)
{
Gto*=4.0;
}
if (celltype==2)
{
Gto*=4.0;
}
double fItop=(1.0/(1.0+KmCaMK/CaMKa));
Ito=(1/(1+pow(((nCmax*Cmax)/ItoIC50),ItonHill)))*Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
double td=pow(3.18,-(T-310)/10)*(0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0))));  //orig td=0.6+...
d=dss-(dss-d)*exp(-dt/td);
double fss=1.0/(1.0+exp((v+19.58)/3.696));
double tff=pow(3.18,-(T-310)/10)*(7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0)));  //orig tff=7.0+...
double tfs=pow(3.18,-(T-310)/10)*(1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0)));  //orig tfs=1000.0+...
double Aff=0.6;
double Afs=1.0-Aff;
ff=fss-(fss-ff)*exp(-dt/tff);
fs=fss-(fss-fs)*exp(-dt/tfs);
double f=Aff*ff+Afs*fs;
double fcass=fss;
double tfcaf=pow(3.18,-(T-310)/10)*(7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0)));  //orig tfcaf=7.0+...
double tfcas=pow(3.18,-(T-310)/10)*(100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0)));  //orig tfcas=100.0+...
double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
double Afcas=1.0-Afcaf;
fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
double fca=Afcaf*fcaf+Afcas*fcas;
double tjca=pow(3.18,-(T-310)/10)*75.0;  //orig tjca=75.0;
jca=fcass-(fcass-jca)*exp(-dt/tjca);
double tffp=2.5*tff;
ffp=fss-(fss-ffp)*exp(-dt/tffp);
double fp=Aff*ffp+Afs*fs;
double tfcafp=2.5*tfcaf;
fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
double fcap=Afcaf*fcafp+Afcas*fcas;
double Kmn=0.002;
double k2n=1000.0;
double km2n=jca*1.0;
double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));
nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
double PhiCaL=4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
double PhiCaNa=1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
double PhiCaK=1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
double zca=2.0;
double PCa=0.0001*scgcal;  //orig PCa=0.0001
if (celltype==1)
{
PCa*=1.2;
}
if (celltype==2)
{
PCa*=2.5;
}
double PCap=1.1*PCa;
double PCaNa=0.00125*PCa;
double PCaK=3.574e-4*PCa;
double PCaNap=0.00125*PCap;
double PCaKp=3.574e-4*PCap;
double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
ICaL=(1/(1+pow(((nCmax*Cmax)/ICaLIC50),ICaLnHill)))*(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaNa=(1/(1+pow(((nCmax*Cmax)/ICaLIC50),ICaLnHill)))*(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaK=(1/(1+pow(((nCmax*Cmax)/ICaLIC50),ICaLnHill)))*(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);

//IKr classical (Hodgkin-Huxley) model
/*double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
double txrf=pow(3.95,(310-T)/10)*(12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38)));  //orig txrf=12.98+...
double txrs=pow(3.95,(310-T)/10)*(1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94)));  //orig txrs=1.865+...
double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
double Axrs=1.0-Axrf;
xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
double xr=Axrf*xrf+Axrs*xrs;
double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
double GKr=0.046*scgkr;  //orig GKr=0.046
if (celltype==1)
{
GKr*=1.3;
}
if (celltype==2)
{
GKr*=0.8;
}
IKr=(1/(1+pow(((nCmax*Cmax)/IKrIC50),IKrnHill)))*GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);*/

//IKr dyn. (Markov) model (from Li,Dutta,Colatsky 2017 & Dutta et al.2017)
    if (t==0) // Steady State Values
        {
            int iindex = 0;
            int ilength = 30000;
            for (iindex = 0; iindex < ilength; iindex++)
            	{comp_statesIKr();}
		}
	else
	{
	comp_statesIKr();
	double GKr=0.0418*scgkr;  //orig GKr=0.046 in ORD2011
//		if (celltype==1)
//			GKr*=1.3;
//		if (celltype==2)
//			GKr*=0.8;
	IKr=GKr*sqrt(ko/5.4)*IKrO*(v-EK);
	}

double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
double txs1=pow(1.51,(310-T)/10)*(817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0)));  //orig txs1=817.3+...
xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
double xs2ss=xs1ss;
double txs2=pow(1.51,(310-T)/10)*1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));  //orig txs2=1.0/...
xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
double KsCa=1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
double GKs=0.0034*scgks;  //orig GKs=0.0034
if (celltype==1)
{
GKs*=1.4;
}
IKs=(1/(1+pow(((nCmax*Cmax)/IKsIC50),IKsnHill)))*GKs*KsCa*xs1*xs2*(v-EKs);

double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
double txk1=pow(1.26,(310-T)/10)*122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));  //orig txk1=122.2/...
xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
double GK1=0.1908*scgki;  //orig GK1=0.1908
if (celltype==1)
{
GK1*=1.2;
}
if (celltype==2)
{
GK1*=1.3;
}
IK1=(1/(1+pow(((nCmax*Cmax)/IK1IC50),IK1nHill)))*GK1*sqrt(ko)*rk1*xk1*(v-EK);

double kna1=15.0;
double kna2=5.0;
double kna3=88.12;
double kasymm=12.5;
double wna=6.0e4;
double wca=6.0e4;
double wnaca=5.0e3;
double kcaon=1.5e6;
double kcaoff=5.0e3;
double qna=0.5224;
double qca=0.1670;
double hca=exp((qca*v*F)/(R*T));
double hna=exp((qna*v*F)/(R*T));
double h1=1+nai/kna3*(1+hna);
double h2=(nai*hna)/(kna3*h1);
double h3=1.0/h1;
double h4=1.0+nai/kna1*(1+nai/kna2);
double h5=nai*nai/(h4*kna1*kna2);
double h6=1.0/h4;
double h7=1.0+nao/kna3*(1.0+1.0/hna);
double h8=nao/(kna3*hna*h7);
double h9=1.0/h7;
double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
double h11=nao*nao/(h10*kna1*kna2);
double h12=1.0/h10;
double k1=h12*cao*kcaon;
double k2=kcaoff;
double k3p=h9*wca;
double k3pp=h8*wnaca;
double k3=k3p+k3pp;
double k4p=h3*wca/hca;
double k4pp=h2*wnaca;
double k4=k4p+k4pp;
double k5=kcaoff;
double k6=h6*cai*kcaon;
double k7=h5*h2*wna;
double k8=h8*h11*wna;
double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
double E1=x1/(x1+x2+x3+x4);
double E2=x2/(x1+x2+x3+x4);
double E3=x3/(x1+x2+x3+x4);
double E4=x4/(x1+x2+x3+x4);
double KmCaAct=150.0e-6;
double allo=1.0/(1.0+pow(KmCaAct/cai,2.0));
double zna=1.0;
double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
double JncxCa=E2*k2-E1*k1;
double Gncx=0.0008*scgncx;  //orig Gncx=0.0008
if (celltype==1)
{
Gncx*=1.1;
}
if (celltype==2)
{
Gncx*=1.4;
}
INaCa_i=pow(2.5,(T-310)/10)*0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);  //Q10 2.5 for NCX1 (Elias CL, Xue XH, Marshall CR, Omelchenko A, Hryshko LV, Tibbits GF. Temperature dependence of cloned mammalian and salmonid cardiac Na(+)/Ca(2+) exchanger isoforms. Am J Physiol Cell Physiol. 2001 Sep;281(3):C993-C1000)

h1=1+nass/kna3*(1+hna);
h2=(nass*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+nass/kna1*(1+nass/kna2);
h5=nass*nass/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+nao/kna3*(1.0+1.0/hna);
h8=nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
h11=nao*nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*cass*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0/(1.0+pow(KmCaAct/cass,2.0));
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
INaCa_ss=pow(2.5,(T-310)/10)*0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);  //Q10 2.5 for NCX1 (Elias CL, Xue XH, Marshall CR, Omelchenko A, Hryshko LV, Tibbits GF. Temperature dependence of cloned mammalian and salmonid cardiac Na(+)/Ca(2+) exchanger isoforms. Am J Physiol Cell Physiol. 2001 Sep;281(3):C993-C1000)

INaCa=INaCa_i+INaCa_ss;

double k1p=949.5;
double k1m=182.4;
double k2p=687.2;
double k2m=39.4;
k3p=1899.0;
double k3m=79300.0;
k4p=639.0;
double k4m=40.0;
double Knai0=9.073;
double Knao0=27.78;
double delta=-0.1550;
double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
double Kki=0.5;
double Kko=0.3582;
double MgADP=0.05;
double MgATP=9.8;
double Kmgatp=1.698e-7;
double H=1.0e-7;
double eP=4.2;
double Khp=1.698e-7;
double Knap=224.0;
double Kxkur=292.0;
double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
double b1=k1m*MgADP;
double a2=k2p;
double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
double zk=1.0;
double JnakNa=3.0*(E1*a3-E2*b3);
double JnakK=2.0*(E4*b1-E3*a1);
double Pnak=30*scgnak;  //orig Pnak=30
if (celltype==1)
{
    Pnak*=0.9;
}
if (celltype==2)
{
    Pnak*=0.7;
}
INaK=Pnak*(zna*JnakNa+zk*JnakK);

double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
double GKb=0.003*scgkb; //orig GKb=0.003
if (celltype==1)
{
    GKb*=0.6;
}
IKb=GKb*xkb*(v-EK);

double PNab=3.75e-10*scgnab; //orig PNab=3.75e-10
INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);

double PCab=2.5e-8*scgcab;  //orig PCab=2.5e-8
ICab=PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);

double GpCa=0.0005*scgpca;  //orig GpCa=0.0005
IpCa=GpCa*cai/(0.0005+cai);

dvy = -10*pow(ach/(ach+0.0000025),4);
    	af = pow(3.0,(T-298)/10)*0.18*6*exp(-F*(v+65-dvy)/(R*T));//orig. (v + 64.0) //orig ...*5*exp(...   Q10 3 from Pena F et al. 2006
    	bf = pow(3.0,(T-298)/10)*0.18*exp(F*(v+65-dvy)/(R*T));//orig ...*5*exp(...    Q10 3 from Pena F et al. 2006

        if (t==0) // Steady State Values
        {
            int iindex = 0;
            int ilength = 30000;
            for (iindex = 0; iindex < ilength; iindex++)
            	{comp_statesIf();}
	}
	else
	{
	comp_statesIf();
	ifk = scgf*gf*gfk*(o1f+o2f)*(v-EK);
	ifna = scgf*gf*gfna*(o1f+o2f)*(v-ENa);
	If = ifna+ifk;
	}


}

void comp_statesIKr()
{
IKrKf2=IKrKf;

if ((cpd=="chloroquine" || cpd=="hydroxychloroquine") && IKrmode=="novshift")
	{
	    actshift=0;
	}
ReactionFlux1=IKrA11*exp(IKrB11*(v+actshift*Cmax*nCmax))*IKrIC1*pow(IKrq11,(T-293)/10)-IKrA21*exp(IKrB21*(v+actshift*Cmax*nCmax))*IKrIC2*pow(IKrq21,(T-293)/10);
//ReactionFlux1=IKrA11*exp(IKrB11*v)*IKrIC1*pow(IKrq11,(T-293)/10)-IKrA21*exp(IKrB21*v)*IKrIC2*pow(IKrq21,(T-293)/10);
ReactionFlux2=IKrA1*exp(IKrB1*(v+actshift*Cmax*nCmax))*IKrC1*pow(IKrq1,(T-293)/10)-IKrA2*exp(IKrB2*(v+actshift*Cmax*nCmax))*IKrC2*pow(IKrq2,(T-293)/10);
//ReactionFlux2=IKrA1*exp(IKrB1*v)*IKrC1*pow(IKrq1,(T-293)/10)-IKrA2*exp(IKrB2*v)*IKrC2*pow(IKrq2,(T-293)/10);
ReactionFlux3=IKrA31*exp(IKrB31*(v+actshift*Cmax*nCmax))*IKrC2*pow(IKrq31,(T-293)/10)-IKrA41*exp(IKrB41*(v+actshift*Cmax*nCmax))*IKrO*pow(IKrq41,(T-293)/10);
//ReactionFlux3=IKrA31*exp(IKrB31*v)*IKrC2*pow(IKrq31,(T-293)/10)-IKrA41*exp(IKrB41*v)*IKrO*pow(IKrq41,(T-293)/10);
ReactionFlux4=IKrA3*exp(IKrB3*(v+actshift*Cmax*nCmax))*IKrIC2*pow(IKrq3,(T-293)/10)-IKrA4*exp(IKrB4*(v+actshift*Cmax*nCmax))*IKrIO*pow(IKrq4,(T-293)/10);
//ReactionFlux4=IKrA3*exp(IKrB3*v)*IKrIC2*pow(IKrq3,(T-293)/10)-IKrA4*exp(IKrB4*v)*IKrIO*pow(IKrq4,(T-293)/10);
ReactionFlux5=IKrA51*exp(IKrB51*v)*IKrC1*pow(IKrq51,(T-293)/10)-IKrA61*exp(IKrB61*v)*IKrIC1*pow(IKrq61,(T-293)/10);
ReactionFlux6=IKrA52*exp(IKrB52*v)*IKrC2*pow(IKrq52,(T-293)/10)-IKrA62*exp(IKrB62*v)*IKrIC2*pow(IKrq62,(T-293)/10);
ReactionFlux7=IKrA53*exp(IKrB53*v)*IKrO*pow(IKrq53,(T-293)/10)-IKrA63*exp(IKrB63*v)*IKrIO*pow(IKrq63,(T-293)/10);
if (cpd=="chloroquine" || cpd=="hydroxychloroquine")
	{
	    ReactionFlux8=IKrkb*Cmax*nCmax*IKrIO-IKrkub*IKrCbound*IKrA53*exp(IKrB53*v)*pow(IKrq53,(T-293)/10)/(IKrA63*exp(IKrB63*v)*pow(IKrq63,(T-293)/10));
		ReactionFlux9=IKrkb*Cmax*nCmax*IKrO-IKrkub*IKrCbound;
	}
else
	{
		ReactionFlux8=IKrKf2*IKrKt2*IKrIO*pow(nCmax*Cmax,IKrnH)/(pow(nCmax*Cmax,IKrnH)+IKrhalfmax)-IKrKt2*IKrCbound/(1+exp(-(v-IKrVhalf)/6.789))*IKrA53*exp(IKrB53*v)*pow(IKrq53,(T-293)/10)/(IKrA63*exp(IKrB63*v)*pow(IKrq63,(T-293)/10));
		ReactionFlux9=IKrKf*IKrKt*IKrO*pow(nCmax*Cmax,IKrnH)/(pow(nCmax*Cmax,IKrnH)+IKrhalfmax)-IKrKt*IKrCbound/(1+exp(-(v-IKrVhalf)/6.789));
	}
//ReactionFlux10=IKrKt*IKrCbound/(1+exp(-(v-IKrVhalf)/6.789))-IKrKt*IKrObound;
//ReactionFlux11=IKrKt2*IKrCbound/(1+exp(-(v-IKrVhalf)/6.789))-IKrKt2*IKrIObound;

IKrIC1+=dt*(1/IKrunnamed)*(-ReactionFlux1+ReactionFlux5);
if (IKrIC1>=1)
	IKrIC1=1;
if (IKrIC1<=0)
	IKrIC1=0;
IKrIC2+=dt*(1/IKrunnamed)*(ReactionFlux1-ReactionFlux4+ReactionFlux6);
if (IKrIC2>=1)
	IKrIC2=1;
if (IKrIC2<=0)
	IKrIC2=0;
IKrC1+=dt*(1/IKrunnamed)*(-ReactionFlux2-ReactionFlux5);
if (IKrC1>=1)
	IKrC1=1;
if (IKrC1<=0)
	IKrC1=0;
IKrC2+=dt*(1/IKrunnamed)*(ReactionFlux2-ReactionFlux3-ReactionFlux6);
if (IKrC2>=1)
	IKrC2=1;
if (IKrC2<=0)
	IKrC2=0;
IKrO+=dt*(1/IKrunnamed)*(ReactionFlux3-ReactionFlux7-ReactionFlux9);
if (IKrO>=1)
	IKrO=1;
if (IKrO<=0)
	IKrO=0;
IKrIO+=dt*(1/IKrunnamed)*(ReactionFlux4+ReactionFlux7-ReactionFlux8);
if (IKrIO>=1)
	IKrIO=1;
if (IKrIO<=0)
	IKrIO=0;
//IKrIObound+=dt*(1/IKrunnamed)*(ReactionFlux8+ReactionFlux11);
//IKrObound+=dt*(1/IKrunnamed)*(ReactionFlux9+ReactionFlux10);
//IKrCbound+=dt*(1/IKrunnamed)*(-ReactionFlux10-ReactionFlux11);
IKrCbound+=dt*(1/IKrunnamed)*(ReactionFlux8+ReactionFlux9);
IKrCbound+=dt*(1/IKrunnamed)*(ReactionFlux8+ReactionFlux9);
if (IKrCbound>=1)
	IKrCbound=1;
if (IKrCbound<=0)
	IKrCbound=0;
}

void comp_statesIf()
{
    	double newc1f = c1f+(c2f*lf+o1f*bf-c1f*(kf+af))*dt*0.001;
    	double newc2f = c2f+(c1f*kf+o2f*bf/fff-c2f*(af*fff+lf))*dt*0.001;
    	double newo1f = o1f+(c1f*af+o2f*lf/fff-o1f*(bf+kf*fff))*dt*0.001;
    	double newo2f = o2f+(c2f*af*fff+o1f*kf*fff-o2f*(lf+bf)/fff)*dt*0.001;

    	c1f = newc1f;
    	c2f = newc2f;
    	o1f = newo1f;
    	o2f = newo2f;
}

void FBC()
{
double CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
CaMKa=CaMKb+CaMKt;
CaMKt+=dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

JdiffNa=(nass-nai)/2.0;
JdiffK=(kss-ki)/2.0;
Jdiff=(cass-cai)/0.2;

double bt=4.75;
double a_rel=0.5*bt;
double Jrel_inf=a_rel*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
if (celltype==2)
{
Jrel_inf*=1.7;
}
double tau_rel=bt/(1.0+0.0123/cajsr);
if (tau_rel<0.005)
{
tau_rel=0.005;
}
Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
double btp=1.25*bt;
double a_relp=0.5*btp;
double Jrel_infp=a_relp*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
if (celltype==2)
{
Jrel_infp*=1.7;
}
double tau_relp=btp/(1.0+0.0123/cajsr);
if (tau_relp<0.005)
{
tau_relp=0.005;
}
Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
Jrel=scgrel*((1.0-fJrelp)*Jrelnp+fJrelp*Jrelp);  //orig Jrel=(1.0-...

double Jupnp=0.004375*cai/(cai+0.00092);
double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
if (celltype==1)
{
Jupnp*=1.3;
Jupp*=1.3;
}
double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
Jleak=0.0039375*cansr/15.0;
Jup=scgup*((1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak);  //orig Jup=(1.0-...

Jtr=(cansr-cajsr)/100.0;

//nai=6.0;
nai+=dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab+gfna*If)*Acap/(F*vmyo)+JdiffNa*vss/vmyo); //modified by addition of If
nass+=dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

//ki=125.0;
ki+=dt*(-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK+gfk*If)*Acap/(F*vmyo)+JdiffK*vss/vmyo); //modified by addition of If
kss+=dt*(-(ICaK)*Acap/(F*vss)-JdiffK);

double Bcai;
if (celltype==1)
{
Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
}
else
{
Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn+cai,2.0));
}
cai+=dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));

double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+cass,2.0)+BSLmax*KmBSL/pow(KmBSL+cass,2.0));
cass+=dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));

cansr+=dt*(Jup-Jtr*vjsr/vnsr);

double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+cajsr,2.0));
cajsr+=dt*(Bcajsr*(Jtr-Jrel));
}

void voltage()
{
v+=-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist+If);  //orig without If
Iqinward=ICaL+INaL;
cqInward[n]+=dt*Iqinward;
Inet=ICaL+INaL+IKr+IKs+IK1+Ito;
Qnet[n]+=dt*Inet;
}

void stimulus()
{
if ((t>(start+n*CL) && t<(start+duration+n*CL-dt)))
	{
	if (Ist==0)
		{
		vrest[n+1]=v;
		cairest[n+1]=cai;
		Qnet[n+1]=0;
		cqInward[n+1]=0;
		vdotmax_flag=1;
		caidotmax_flag=1;
		caimax_flag=1;
		EAD_flag[n]=0;
		}
	Ist=amp;
	}
else if (t>(start+duration+n*CL-dt))
	{
	Ist=0.0;
	n=n+1;
	}
}

void dVdt_APD_CaD()
{
vdot_old=vdot;
vdot=(v-vo)/dt;
caidot_old=caidot;
caidot=(cai-caio)/dt;
if (vdotmax_flag==1 && APD50_flag==0 && APD90_flag==0 && v>-40 && vdot<vdot_old)
	{
	vdot_max[n]=vdot_old;
	t_vdot_max[n]=t-dt;
	vdotmax_flag=0;
	APD50_flag=1;
	APD90_flag=1;
	}
if (caidotmax_flag==1 && CaD50_flag==0 && CaD90_flag==0 && v>-40 && caidot<caidot_old)
	{
	caidot_max[n]=caidot_old;
	t_caidot_max[n]=t-dt;
	caidotmax_flag=0;
	CaD50_flag=1;
	CaD90_flag=1;
	}
if (vdotmax_flag==0 && APD50_flag==1 && v>vo)
	vmax[n]=v;
if (caidotmax_flag==0 && caimax_flag==1 && CaD50_flag==1 && caidot_old>0 && caidot<=0)
	{caimax[n]=cai;
	caimax_flag=0;}
if	(vdotmax_flag==0 && APD50_flag==1 && v<(vmax[n]-0.5*(vmax[n]-vrest[n])))
	{
	APD50[n]=t-t_vdot_max[n];
	APD50_flag=0;
	}
if	(vdotmax_flag==0 && APD50_flag==0 && APD90_flag==1 && v<(vmax[n]-0.9*(vmax[n]-vrest[n])))
	{
	APD90[n]=t-t_vdot_max[n];
	APD90_flag=0;
	}
if	(caidotmax_flag==0 && caimax_flag==0 && CaD50_flag==1 && cai<(caimax[n]-0.5*(caimax[n]-cairest[n])))
	{
	CaD50[n]=t-t_caidot_max[n];
	CaD50_flag=0;
	}
if	(caidotmax_flag==0 && caimax_flag==0 && CaD50_flag==0 && CaD90_flag==1 && cai<(caimax[n]-0.9*(caimax[n]-cairest[n])))
	{
	CaD90[n]=t-t_caidot_max[n];
	CaD90_flag=0;
	}
if  (t>(start+duration+(n-1)*CL+20) && v>vrest[n] && (vdot_old<=0 && vdot>0) && EAD_flag[n]==0)
	EAD_flag[n]=1;
}
