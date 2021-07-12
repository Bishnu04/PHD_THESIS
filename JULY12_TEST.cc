// July 12, 2021
// Author Bishnu Pandey
// Al data involved in tune
extern double calcf2t_th(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_ph(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_mom(double* P, 
			  double xf, double xpf,
			  double yf, double ypf,double);

const double  XFPm=-0.7,  XpFPm=-0.15; // m is the mean from the old definition
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; // tm = target offset.. MOmm is the momentum offset
const double  XFPr=1.3,   XpFPr=0.27; // r is the scaling factor or range
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; // tr is the target range
const double  PLm = 25.4, PLr=0.7; // m is the offset and PLr is the path laegth range
const double  Ztm = -0.15,Ztr=0.35; //Ztm  z position at target  point offset
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double tune(double* pa, int j);

const int nmax = 3000; // was 3000 before Dec5
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;

const int npeak = 2;
double Lambda_width[npeak] = {2.37, 2.3}; // +/- from mean position
double Lambda_cent[npeak] ={1115.75,1192.6}; // current peak location
double Lambda_real[npeak] ={1115.683,1192.642}; // Mev// nominal

double p10[nmax],p11[nmax],p12[nmax];
double p13[nmax],p14[nmax],p15[nmax];
double p16[nmax],p17[nmax],p18[nmax],p19[nmax];
double phir[nmax];
double phil[nmax];
// ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
const int nmax_2 = 3000;   // 2400 before Dec 5
double x_2[nmax_2], y_2[nmax_2]; 
double xp_2[nmax_2], yp_2[nmax_2];
double z_recon_2[nmax_2];
int foil_flag_2[nmax_2];

const int npeak_2 = 1;
double Lambda_width_2[npeak_2] = {2.38}; //6.5,8.8}
double Lambda_cent_2[npeak_2] ={1115.71};
double Lambda_real_2[npeak_2] ={1115.683}; // Mev

double p10_2[nmax_2],p11_2[nmax_2],p12_2[nmax_2];
double p13_2[nmax_2],p14_2[nmax_2],p15_2[nmax_2];
double p16_2[nmax_2];
double phir_2[nmax_2];
double phil_2[nmax_2];
int ntune_event_2 = 0;
//===============================================Jan 07 2020==========================================
const int nmax_4 = 3000;   
double x_4[nmax_4], y_4[nmax_4]; 
double xp_4[nmax_4], yp_4[nmax_4];
double z_recon_4[nmax_4];
int foil_flag_4[nmax_4];

const int npeak_4 = 3; // JAn 07
double Lambda_width_4[npeak_4] ={1.63,2.38,2.7};// +/- from mean position
double Lambda_cent_4[npeak_4] ={-13.24,-2.125,6.49};//current peak location
double Lambda_real_4[npeak_4] ={-13.24,-2.125,6.49};// nominal position

double p10_4[nmax_4],p11_4[nmax_4],p12_4[nmax_4];
double p13_4[nmax_4],p14_4[nmax_4],p15_4[nmax_4];
double p16_4[nmax_4];
double phir_4[nmax_4];
double phil_4[nmax_4];
int ntune_event_4 = 0;
//))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
// ((((((((((((((((((((((((((((((((   t3   (((((((((((((((((((((((((((((((((((((((((((((((((((   // Jan_02
const int nmax_3 =3000;
double x_3[nmax_3], y_3[nmax_3];
double xp_3[nmax_3], yp_3[nmax_3];
double z_recon_3[nmax_3];
int foil_flag_3[nmax_3];

const int npeak_3 = 3;
double Lambda_width_3[npeak_3] ={1.63,2.38,2.7};// +/- from mean position
double Lambda_cent_3[npeak_3] ={-13.24,-2.125,6.49};//current peak location
double Lambda_real_3[npeak_3] ={-13.24,-2.125,6.49};// nominal position

double p10_3[nmax_3],p11_3[nmax_3],p12_3[nmax_3];
double p13_3[nmax_3],p14_3[nmax_3],p15_3[nmax_3];
double p16_3[nmax_3];
double phir_3[nmax_3];
double phil_3[nmax_3];
int ntune_event_3 = 0;

const double m_Al =25.1267; // Al target mass // BY Dr. Tang on Dec 19 2019
const double m_T = 2.808921; // for tritium target by Gogami Tritium target mass
//))))))))))))))))))))))))))))))))   t3  ))))))))))))))))))))))))))))))))))))))))))))))))))))))
//========================================
const int Total_Par = 126;
double thetaL_opt[nmax];
double phiL_opt[nmax];
double thetaR_opt[nmax];
double phiR_opt[nmax];
double momL_opt[nmax];
double momR_opt[nmax];
const int Mom_Par = 252;
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
const double me = 0.000511;
const double mk = 0.493677;
const double mp = 0.938272;
const double mL = 1.115683;

extern double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt);
void  JULY12_TEST(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ========================================  
  TChain * t1 = new TChain("T");  
  TChain * t2 = new TChain("T"); 
  TChain * t3 = new TChain("T");
 
  t1->Add("./Rootfiles/DEC17_Rootfiles/DEC17_H149_542.root");// replayed on Dec 17, 2019 replay by ole
  t2->Add("./Rootfiles/DEC17_Rootfiles/DEC17_HT_552_716.root");
  t3->Add("./Rootfiles/DEC17_Rootfiles/DEC23_T221_830.root");  
   
  double ent = t1->GetEntries();
  double ent_2 = t2->GetEntries();
  double ent_3 = t3->GetEntries(); 

  // ent = 50;
  // ent_2= 50;
  //  ent_3 =50;
  
  cout<<"entry in the t1=="<<ent<<endl;
  cout<<"entry in the t2=="<<ent_2<<endl;
  cout<<"entry in the t3=="<<ent_3<<endl;   
  
  const int max = 100;
  Double_t trig5[max]; 
  double momL[max];
  double momR[max]; 
  
  double lvz[max],rvz[max];// raster corrected 
  double th1[max], ph1[max];// RHRS angle 
  double th2[max], ph2[max];  
  double delta_pep[max];     // target straggling
  double pep_real[max]; 
  double delta_pk[max];
  double pk_real[max];
  double par_ep[3];
  double par_k[3];
  double mm; 
  double hallap;

  double mmT_T; //sept4, 2020.. H/T considering as Tritium target
  
  double l_th_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double l_y_fp[max];

  double r_th_fp[max];
  double r_ph_fp[max];
  double r_x_fp[max];
  double r_y_fp[max];
  const int n = 16; 
  double ctime;  
 
  double z_av[nmax];
  double z_av_1[nmax];
  double a1, a2;
  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  Double_t trig5_2[max]; // JUly 01, 2019 
  double momL_2[max];
  double momR_2[max]; 
  
  double lvz_2[max],rvz_2[max];// raster corrected 
  double th1_2[max], ph1_2[max];// RHRS angle 
  double th2_2[max], ph2_2[max];  
  double delta_pep_2[max];     // target straggling
  double pep_real_2[max]; 
  double delta_pk_2[max];
  double pk_real_2[max];
  double par_ep_2[3];
  double par_k_2[3];
  double mm_2;
  double mm_4;
  double mm_ht;
  double hallap_2;

  double par_ht_ep_2[3]; 
  double par_ht_k_2[3]; 
  double pep_ht_real_2[max]; 
  double pk_ht_real_2[max]; 
  double delta_ht_pep_2[max];  
  double delta_ht_pk_2[max]; 
  double hallap_ht_2; 
  double z_av_ht_2[nmax];    
  
  double l_th_fp_2[max];
  double l_ph_fp_2[max];
  double l_x_fp_2[max];
  double l_y_fp_2[max];

  double r_th_fp_2[max];
  double r_ph_fp_2[max];
  double r_x_fp_2[max];
  double r_y_fp_2[max];
  double ctime_2; 
 
  double z_av_2[nmax];
  double z_av_1_2[nmax];
  double a1_2, a2_2;
  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t3 variables starts here    ((((((((((((((((((((((((((((((((((((((((((((((((((( Jan_02
  Double_t trig5_3[max]; // JUly 01, 2019 
  double momL_3[max];
  double momR_3[max]; 
  double lvz_3[max],rvz_3[max];// raster corrected 
  double th1_3[max], ph1_3[max];// RHRS angle 
  double th2_3[max], ph2_3[max];  
  double delta_pep_3[max];     // target straggling
  double pep_real_3[max]; 
  double delta_pk_3[max];
  double pk_real_3[max];
  double par_ep_3[3];
  double par_k_3[3];

  double par_tt_ep_3[3];//// May 01, 2021 , tt means T/T 
  double par_tt_k_3[3];//// May 01, 2021 , tt means T/T
  double hallap_tt_3;//// May 01, 2021 , tt means T/T
  double z_av_tt_3[nmax];//// May 01, 2021 , tt means T/T
  double delta_tt_pep_3[max];//// May 01, 2021 , tt means T/T
  double delta_tt_pk_3[max];//// May 01, 2021 , tt means T/T
  double pep_tt_real_3[max];//// May 01, 2021 , tt means T/T
  double pk_tt_real_3[max];//// May 01, 2021 , tt means T/T   
  
  double mm_3;
  double mm_t;  
  double mm_Al;
  double mm_Al1;
  double a1_3, a2_3;
  double mm_h;
  double hallap_3;
  double l_th_fp_3[max];
  double l_ph_fp_3[max];
  double l_x_fp_3[max];
  double l_y_fp_3[max];
  double r_th_fp_3[max];
  double r_ph_fp_3[max];
  double r_x_fp_3[max];
  double r_y_fp_3[max];
  double ctime_3;   
  double z_av_3[nmax];
  double z_av_1_3[nmax];
  //))))))))))))))))))))))))))))))))   t3 variables up to here   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t1 branch address  (((((((((((((((((((((((((((((((((((((((((((((((((((
  t1->SetBranchAddress("HALLA_p", &hallap);  
  t1->SetBranchAddress("DR.T5", &trig5);  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);   
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp); 
  t1->SetBranchAddress("coin_time",  &ctime); 
  t1->SetBranchAddress("ztR_wRC",  &rvz);
  t1->SetBranchAddress("ztL_wRC",  &lvz);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  t2->SetBranchAddress("HALLA_p", &hallap_2);  
  t2->SetBranchAddress("DR.T5", &trig5_2);  
  t2->SetBranchAddress("L.tr.x",   &l_x_fp_2);
  t2->SetBranchAddress("L.tr.y",   &l_y_fp_2);
  t2->SetBranchAddress("L.tr.th",  &l_th_fp_2);
  t2->SetBranchAddress("L.tr.ph",  &l_ph_fp_2);   
  t2->SetBranchAddress("R.tr.x",   &r_x_fp_2);
  t2->SetBranchAddress("R.tr.y",   &r_y_fp_2);
  t2->SetBranchAddress("R.tr.th",  &r_th_fp_2);
  t2->SetBranchAddress("R.tr.ph",  &r_ph_fp_2); 
  t2->SetBranchAddress("coin_time",  &ctime_2); 
  t2->SetBranchAddress("ztR_wRC",  &rvz_2);
  t2->SetBranchAddress("ztL_wRC",  &lvz_2);
  t2->SetBranchAddress("R.a1.asum_c", &a1_2);
  t2->SetBranchAddress("R.a2.asum_c", &a2_2);
  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t3   ((((((((((((((((((((((((((((((((((((((((((((((((((( Jan_02
  t3->SetBranchAddress("HALLA_p", &hallap_3);  
  t3->SetBranchAddress("DR.T5", &trig5_3);   
  t3->SetBranchAddress("L.tr.x",   &l_x_fp_3);
  t3->SetBranchAddress("L.tr.y",   &l_y_fp_3);
  t3->SetBranchAddress("L.tr.th",  &l_th_fp_3);
  t3->SetBranchAddress("L.tr.ph",  &l_ph_fp_3);
  t3->SetBranchAddress("R.tr.x",   &r_x_fp_3);
  t3->SetBranchAddress("R.tr.y",   &r_y_fp_3);
  t3->SetBranchAddress("R.tr.th",  &r_th_fp_3);
  t3->SetBranchAddress("R.tr.ph",  &r_ph_fp_3); 
  t3->SetBranchAddress("coin_time",  &ctime_3); 
  t3->SetBranchAddress("ztR_wRC",  &rvz_3);
  t3->SetBranchAddress("ztL_wRC",  &lvz_3);
  t3->SetBranchAddress("R.a1.asum_c", &a1_3);
  t3->SetBranchAddress("R.a2.asum_c", &a2_3);
  //))))))))))))))))))))))))))))))))   t3 barnch address up to here    ))))))))))))))))))))))))))))))))))))))))))))))))))))))  
  TFile* fnew = new TFile("./output_root/paper_prep.root","recreate"); 
  TTree* tnew = new TTree("tree","For z calibration (LHRS)");
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D");
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D");
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D");
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D");
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D");  
  tnew->Branch("L.tr.tg_th_TH2", &th2, "L.tr.tg_th_TH2[100]/D");
  tnew->Branch("L.tr.tg_ph_PH2", &ph2, "L.tr.tg_ph_PH2[100]/D");  
 
  double XFP, XpFP;
  double YFP, YpFP;
  double R_XFP, R_XpFP; 
  double R_YFP, R_YpFP;
  // ((((((((((((((((((((((((((((((((((((((((((( for t2 ((((((((((
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double R_XFP_2, R_XpFP_2; 
  double R_YFP_2, R_YpFP_2;
  //)))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((((( for t3 (((((((((( Jan_02
  double XFP_3, XpFP_3;
  double YFP_3, YpFP_3;
  double R_XFP_3, R_XpFP_3; 
  double R_YFP_3, R_YpFP_3; 
  //)))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((((( for t2 ((((((((((
  double XFP_4, XpFP_4;
  double YFP_4, YpFP_4;  // jan 07 2020
  double R_XFP_4, R_XpFP_4; 
  double R_YFP_4, R_YpFP_4;
  // ===============or LHRS  theta information input==========   
  ntune_event = 0;
  for(int i=0 ; i<Total_Par; i++){
    thetaL_opt[i] = -2222.0;
  }
  
  char name_Angle_L[500]; 
  // sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");//theta_3rd_LHRS_Opt_7.dat
  sprintf(name_Angle_L,"./matrices/theta_L4th_4th_6.dat");// optimized on OCT 23, 2019 with SS data
  ifstream Angle_L(name_Angle_L);
  double Theta_L[Total_Par];    
  for(int i =0; i<Total_Par;i++){
    double par1 =0.0;
    int p1 =0;
    Angle_L>>par1>>p1>>p1>>p1>>p1>>p1;
    Theta_L[i]=par1;
    thetaL_opt[i] = Theta_L[i];
  }
  Angle_L.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // for LHRS  phi information input
  //--------------------------------------------------------  
  ntune_event = 0;
  for(int i =0;i<Total_Par;i++){
    phiL_opt[i] = -2222.0;
  }
  char name_angle_phi[500];
  //  sprintf(name_angle_phi,"./matrices/phi_LHRS_3rd_Opt_9.dat");
  sprintf(name_angle_phi,"./matrices/phi_L4th_5th_5.dat");// optimized on OCT 23, 2019  
  ifstream angle_phi(name_angle_phi);
  double PHI_L[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par2 =0.0;
    int p2 =0;
    angle_phi>>par2>>p2>>p2>>p2>>p2>>p2;
    PHI_L[i] = par2;
    phiL_opt[i]= PHI_L[i];
  }
  angle_phi.close();
  // LHRS momentum information========================July 20, 2019
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momL_opt[i] = -2222.0;
  }
  char name_Mom_lhrs[500];  
  //  sprintf(name_Mom_lhrs,"./MOM_MATRICES/LMOM5_1st_0.dat"); /// Previously optimized without Al Data May 02, 2021
  sprintf(name_Mom_lhrs,"./MOM_MATRICES/LMOM5_July8_1st_2.dat");// tuned with Al data (no shift is included)
  ifstream Mom_lhrs(name_Mom_lhrs);
  double mom_L[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par5 = 0.0;
    int p5 =0;
    Mom_lhrs>>par5>>p5>>p5>>p5>>p5>>p5;
    mom_L[i]= par5;
    momL_opt[i] = mom_L[i];
  }
  Mom_lhrs.close();
  // up to here thelhrs momentum matrix======================  
  // =======RHRS theta input information
  ntune_event =0;
  for(int i =0;i<Total_Par;i++){
    thetaR_opt[i] = -2222.0;
  }
  char name_Angle_R[500]; 
  sprintf(name_Angle_R,"./All_Matrices/xpt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  ifstream Angle_R(name_Angle_R);
  double Theta_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par3 =0.0;
    int p3 = 0;
    Angle_R>>par3>>p3>>p3>>p3>>p3>>p3;
    Theta_R[i]=par3;
    thetaR_opt[i] = Theta_R[i];
  }
  Angle_R.close();
  //====================================================
  //=======RHRS phi input information===============
  ntune_event = 0;
  for(int i =0;i<Total_Par;i++){
    phiR_opt[i] = -2222.0;
  }
  char name_phi_Rhrs[500]; 
  sprintf(name_phi_Rhrs,"./All_Matrices/ypt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  ifstream phi_Rhrs(name_phi_Rhrs);
  double PHI_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par4 =0.0;
    int p4 =0;
    phi_Rhrs>>par4>>p4>>p4>>p4>>p4>>p4;
    PHI_R[i] = par4;
    phiR_opt[i]= PHI_R[i];
  }
  phi_Rhrs.close();
  //==================================================
  // =====RHRS momentum recon==========================6 
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momR_opt[i] = -2222.0;
  }
  char name_Mom_rhrs[500];   
  // sprintf(name_Mom_rhrs,"./MOM_MATRICES/RMOM5_2nd_0.dat"); // Optimized without Al data May 02, 2021
  sprintf(name_Mom_rhrs,"./MOM_MATRICES/RMOM5_July8_1st_0.dat");// tuned with Al data (no shift is included)  
  ifstream Mom_rhrs(name_Mom_rhrs);
  double mom_R[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par6 = 0.0;
    int p6 =0;
    Mom_rhrs>>par6>>p6>>p6>>p6>>p6>>p6;
    mom_R[i]= par6;
    momR_opt[i] = mom_R[i];
  }
  Mom_rhrs.close();
  // =====RHRS momentum recon up to here=============== 
  gStyle->SetTickLength(0.055,"X");
  TH1F *h = new TH1F("h"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); // to plot lambda and sigma (H/H data)  
  TH1F *h_2 = new TH1F("h_2"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); //Lambda (H/T data)  
   
  TH1F *hh = new TH1F("hh","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",208,-100,212);
  TH1F *ht = new TH1F("ht","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",208,-100,212);
  TH1F *htt = new TH1F("htt","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",208,-100,212);

  TH1F *hb = new TH1F("hb","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",208,-100,212);
  TH1F *hb1 = new TH1F("hb1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",208,-100,212);
  TH1F *hb2 = new TH1F("hb2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",208,-100,212);

  TH1F *hal_15 = new TH1F("hal_15","Al Spectrum,(H/T + T/T data)  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",151,-100,202);
  TH1F *hal_15_1 = new TH1F("hal_15_1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",151,-100,202);
  TH1F *hal_15_2 = new TH1F("hal_15_2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/2.0 MeV ",151,-100,202);

  TH1F *hbal_15 = new TH1F("hbal_15","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",151,-100,202);
  TH1F *hbal_15_1 = new TH1F("hbal_15_1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",151,-100,202);
  TH1F *hbal_15_2 = new TH1F("hbal_15_2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/2.0 MeV ",151,-100,202);  

  TH1F *h_al = new TH1F("h_al","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",240,-100,200);
  TH1F *h_al1 = new TH1F("h_al1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",240,-100,200);
  TH1F *h_al2 = new TH1F("h_al2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",240,-100,200);

  TH1F *hb_al = new TH1F("hb_al","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",240,-100,200);
  TH1F *hb_al1 = new TH1F("hb_al1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",240,-100,200);
  TH1F *hb_al2 = new TH1F("hb_al2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/MeV ",240,-100,200);

  TH1F *hal_20 = new TH1F("hal_20","Al Spectrum (H/T+T/T data)  ; -B_{#Lambda}(MeV);Counts/2 MeV ",155,-100,210);
  TH1F *hal_20_1 = new TH1F("hal_20_1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/2 MeV ",155,-100,210);
  TH1F *hal_20_2 = new TH1F("hal_20_2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/2 MeV ",155,-100,210);

  TH1F *hbal_20 = new TH1F("hbal_20","Al Spectrum, H/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",155,-100,210);
  TH1F *hbal_20_1 = new TH1F("hbal_20_1","Al Spectrum, T/T data ; -B_{#Lambda}(MeV);Counts/2 MeV ",155,-100,210);
  TH1F *hbal_20_2 = new TH1F("hbal_20_2","Al Spectrum(H/T + T/T data); -B_{#Lambda}(MeV);Counts/2 MeV ",155,-100,210);  
  
  TH1F *h20 = new TH1F("h20","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",244,-100.5,204.5);
  TH1F *h21 = new TH1F("h21","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",244,-100.5,204.5 ); 
  TH1F *HT = new TH1F("HT","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",244,-100.5,204.5);
 
  TH1F *h20_b = new TH1F("h20_b","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",244,-100.5,204.5);
  TH1F *h21_b = new TH1F("h21_b","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",244,-100.5,204.5); 
  TH1F *HT_b = new TH1F("HT_b","Al Spectrum( H/T+ TT data); -B_{#Lambda}(MeV);Counts/ 1.25 MeV ",244,-100.5,204.5);

  TH1F *hd = new TH1F("hd","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",202,-100.5,202.5);//need here
  TH1F *he = new TH1F("he","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/1.5MeV ",202,-100.5,202.5); 
  TH1F *hf = new TH1F("hf","Al Spectrum(H/T+ TT); -B_{#Lambda}(MeV);Counts/1.5MeV ",202,-100.5,202.5);

  TH1F *hb5 = new TH1F("hb5","Al Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",202,-100.5,202.5);
  TH1F *hb6 = new TH1F("hb6","Al Spectrum,H/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",202,-100.5,202.5); 
  TH1F *hb7 = new TH1F("hb7","Al Spectrum(H/T+ TT); -B_{#Lambda}(MeV);Counts/ 1.5 MeV ",202,-100.5,202.5);
 
  TH1F *h50 = new TH1F("h50","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-100.5,149.5);
  TH1F *h50b = new TH1F("h50b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-100.5,149.5);
  
  TH1F *h25_2 = new TH1F("h25_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",102,-99.5,154.5);  
  TH1F *h25_2b = new TH1F("h25_2b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",102,-99.5,154.5);  

  TH1F *h30_1 = new TH1F("h30_1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.3 MeV ",110,-100,153);
  TH1F *h30_2 = new TH1F("h30_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.3 MeV ",110,-101,152);
  TH1F *h30_3 = new TH1F("h30_3","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.3 MeV ",110,-99,154);
 
  TH1F *h35_1 = new TH1F("h35_1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.35 MeV ",108,-100,154);
  TH1F *h35_2 = new TH1F("h35_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.35 MeV ",108,-104,150);
  TH1F *h35_3 = new TH1F("h35_3","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.35 MeV ",108,-99,155);
  TH1F *h35_4 = new TH1F("h35_4","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.35 MeV ",108,-98,156);
 
  TH1F *h54 = new TH1F("h54","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",125,-100,150); 
  TH1F *hb54 = new TH1F("hb54","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",125,-100,150);

  TH1F *h1_2 = new TH1F("h1_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",128,-104,152);  // C/1.52 in really
  TH1F *h1_2b = new TH1F("h1_2b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",128,-104,152); // C/1.52 in really

  TH1F *h1_qu = new TH1F("h1_qu","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.6 MeV ",160,-100,156);
  TH1F *h1_qub = new TH1F("h1_qub","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.6 MeV ",160,-100,156);

  TH1F *h1_q1 = new TH1F("h1_q1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-102,150);
  TH1F *h1_qb1 = new TH1F("h1_qb1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-102,150);
  /////////////////////++++++++++ Dec01, 2020&&&&&&&&&&&&&&&&&&##################@@@@@@@@@@@@@@@@@!!!!!!!!!
  TH1F *h1_0 = new TH1F("h1_0","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",174,-110,151);

  TH1F *h1_01 = new TH1F("h1_01","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-100,150);
  TH1F *h1_01b = new TH1F("h1_01b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-100,150);

  TH1F *h1_may = new TH1F("h1_may","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",173,-103.5,156.5);
  TH1F *h1b_may = new TH1F("h1b_may","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",173,-103.5,156.5);
  
  TH1F *h1_02 = new TH1F("h1_02","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",100,-100,150);
  TH1F *h1_02b = new TH1F("h1_02b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",100,-100,150);
  
  TH1F *h1_03 = new TH1F("h1_03","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",127,-102.5,151.5);
  TH1F *h1_03b = new TH1F("h1_03b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",127,-102.5,151.5);
  
  TH1F *h1_04 = new TH1F("h1_04","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",102,-105,150);
  TH1F *h1_05 = new TH1F("h1_05","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",102,-100,155);
  TH1F *h1_05b = new TH1F("h1_05b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.5 MeV ",102,-100,155);
   
  TH1F *h150_1 = new TH1F("h150_1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",130,-100,160);
  TH1F *h150_1b = new TH1F("h150_1b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",130,-100,160);

  TH1F *h150_2 = new TH1F("h150_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",126,-101,151);  
  TH1F *h150_2b = new TH1F("h150_2b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",126,-101,151);

  TH1F *h150_3 = new TH1F("h150_3","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",127,-102,152);
  
  TH1F *h150_4 = new TH1F("h150_4","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",130,-105,155);
  TH1F *h150_5 = new TH1F("h150_5","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",125,-100,150);
  TH1F *h150_6 = new TH1F("h150_6","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",132,-104,160);

  TH1F *h125_1 = new TH1F("h125_1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",112,-100,152);
  TH1F *h125_2 = new TH1F("h125_2","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",112,-104,148);
  
  TH1F *h125_3 = new TH1F("h125_3","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",112,-102,150);
  TH1F *h125_3b = new TH1F("h125_3b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",112,-102,150);
 
  TH1F *h125_4 = new TH1F("h125_4","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",112,-105,147);
  TH1F *h125_4b = new TH1F("h125_4b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",112,-105,147);

  TH1F *h125_5 = new TH1F("h125_5","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",116,-103,158);
  TH1F *h125_5b = new TH1F("h125_5b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2.25 MeV ",116,-103,158);

  TH1F *h125_6 = new TH1F("h125_6","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100,152);
  TH1F *h125_6b = new TH1F("h125_6b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100,152);
  /////////////////////++++++++++ Dec01, 2020&&&&&&&&&&&&&&&&&&##################@@@@@@@@@@@@@@@@@!!!!!!!!!
  TH1F *h_75_1 = new TH1F("h_75_1","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ 0.76 MeV ",350,-100,166); 
  TH1F *h_75_3 = new TH1F("h_75_3","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ 0.76 MeV ",350,-103,163); 
  TH1F *h_75_1b = new TH1F("h_75_1b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ 0.76 MeV ",350,-100,166); 
  TH1F *h_75_3b = new TH1F("h_75_3b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/ 0.76 MeV ",350,-103,163);  

  // H contamination
  TH1F *h_h = new TH1F("h_h","  ;-B_{#Lambda}(MeV);Counts / 1.5 MeV ",134,-101,100); // H in  T/T data or to see H contamination
  TH1F *h_hbg = new TH1F("h_hbg","  ;-B_{#Lambda}(MeV);Counts / 1.5 MeV ",134,-101,100);// bg in H contamination

  TH1F *h_hc0 = new TH1F("h_hc0"," H Contamination;-B_{#Lambda}(MeV);Counts / 1.5 MeV ",120,-90,90);
  TH1F *h_hc0b = new TH1F("h_hc0b","H Contamination;-B_{#Lambda}(MeV);Counts / 1.5 MeV ",120,-90,90);
  
  TH1F *h_hc1 = new TH1F("h_hc1"," H Contamination;-B_{#Lambda}(MeV);Counts / 1.5 MeV ",128,-96,96);
  TH1F *h_hc1b = new TH1F("h_hc1b","H contamination;-B_{#Lambda}(MeV);Counts / 1.5 MeV ",128,-96,96); 
  
  TH1F *H_T = new TH1F("H_T","H in  T/T data  ;-B_{#Lambda}(MeV);Counts/ 1.52MeV ",150,-100,128);
  TH1F *H_TB = new TH1F("H_TB","H in  T/T data  ;-B_{#Lambda}(MeV);Counts/ 1.52MeV ",150,-100,128);
    
  TH1F *h53_t = new TH1F("h53_t","nnL Spectrum, T/T data(-100.5, 151.5) ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100.5,151.5);
  TH1F *h53_tb = new TH1F("h53_tb","nnL Spectrum, T/T data(-100.5, 151.5) ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100.5,151.5);

  gStyle->SetOptStat(111111);  
  TH1F *h_1 = new TH1F("h_1"," ",120,1000, 1300);
  
  h->GetXaxis()->SetTitle("Missing Mass (MeV/c^{2})"); // HH data for lambda and sigma
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleFont(62);//  32 gives the times italic bold
  h->GetXaxis()->SetTitleOffset(1.05);  
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetLabelSize(0.06);
 
  h->GetYaxis()->SetTitle("Counts / 0.75 MeV");
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetTitleFont(62);//  32 gives the times italic bold
  h->GetYaxis()->SetTitleOffset(0.80);
  h->GetYaxis()->SetLabelSize(0.06);

  h_h->GetXaxis()->SetTitle("-B_{#Lambda} (MeV)"); // Tritium target for H kinematics to see the H contamination
  h_h->GetXaxis()->SetTitleSize(0.07);
  h_h->GetXaxis()->SetTitleFont(62);//  32 gives the times italic bold
  h_h->GetXaxis()->SetTitleOffset(0.85);  
  h_h->GetXaxis()->CenterTitle();
  h_h->GetXaxis()->SetLabelSize(0.06);
 
  h_h->GetYaxis()->SetTitle("Counts / 1.5 MeV");
  h_h->GetYaxis()->CenterTitle();
  h_h->GetYaxis()->SetTitleSize(0.07);
  h_h->GetYaxis()->SetTitleFont(62);//  32 gives the times italic bold
  h_h->GetYaxis()->SetTitleOffset(0.70);
  h_h->GetYaxis()->SetLabelSize(0.06);
  
  TH1F *h6 = new TH1F("h6",";RHRS reconstructed Momentum;Counts/ 14.4 mev",250,1.7,2.0); 
  char tempc[500];
  // ======================================================  
  bool rtrig = false; 
  for(int i=0; i<nmax; i++){
    x[i]    = -2222.0; 
    y[i]    = -2222.0; 
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    z_av[i] = -2222.0;
    z_av_1[i] = -2222.0;
    phir[i] = -2222.0;
    phil[i] = -2222.0;
    z_recon[i] = -2222.0;  /// Jan 04, 2019
    foil_flag[i] = -1;
  }
  // ((((((((((((((((((((((((((((((((((((((((((((
  bool rtrig_2 = false; 
  for(int i=0; i<nmax_2; i++){
    x_2[i]    = -2222.0; 
    y_2[i]    = -2222.0; 
    xp_2[i]   = -2222.0;
    yp_2[i]   = -2222.0;
    z_av_2[i] = -2222.0;
    z_av_1_2[i] = -2222.0;
    phir_2[i] = -2222.0;
    phil_2[i] = -2222.0;
    z_recon_2[i] = -2222.0; ///Jan 04, 2019
    foil_flag_2[i] = -1;

    x_4[i]    = -2222.0; 
    y_4[i]    = -2222.0; 
    xp_4[i]   = -2222.0;
    yp_4[i]   = -2222.0; 
   
    phir_4[i] = -2222.0;
    phil_4[i] = -2222.0;
    z_recon_4[i] = -2222.0; ///Jan 04, 2019
    foil_flag_4[i] = -1;
  }  
  // +++++++++++++++++++++++++ for t1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (int i=0; i< ent; i++){
    for(int j=0; j<max; j++){
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0; 
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      th1[j] = -2222.0;
      th2[j] = -2222.0;
      ph1[j] =-2222.0;
      ph2[j] =-2222.0;    
     
      delta_pep[j]= -2222.0;
      pep_real[j] =-2222.0;
      delta_pk[j]= -2222.0;
      pk_real[j] = -2222.0;
     
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;      
      
      trig5[j] = 0.0;
      rtrig = false;
    }   
    trig5[0] = 0.0;
    rtrig = false;   
    t1->GetEntry(i);   
   
    if(trig5[0]>1.0) rtrig = true; //JUly 01, 2019
    else rtrig = false;

    z_av[0] = (lvz[0] + rvz[0])/2.0;
    z_av_1[0] =  z_av[0];
   
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    R_XFP   = r_x_fp[0];
    R_XpFP  = r_th_fp[0];
    R_YFP   = r_y_fp[0];
    R_YpFP  = r_ph_fp[0];
    
    if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.040 & fabs(z_av_1[0])<0.10){ // to plot lambda and sigma      
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      
      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;

      z_av[0] =(z_av[0]- Ztm)/Ztr;
      
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av_1[0]);
      th2[0] = th2[0]*Xptr + Xptm; 
     
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP, z_av_1[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;          
     
      momL[0] =  calcf2t_mom(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]); 
      momL[0] = momL[0]*Momr + Momm;      

      par_ep[1] = -th2[0]; // right handed system
      par_ep[2] = -ph2[0]; // right handed system
      double holiang;
      
      // Target struggling LHRS step #7
      if(z_av_1[0]<8.0e-2){	
	
	holiang =  par_ep[2] + hrs_ang;
	holiang=-holiang;
	delta_pep[0] = -1.35758 * sin(-4.59571 * holiang) + 2.09093;

      } 
      else{
	holiang =  par_ep[2] + hrs_ang;
	holiang=-holiang;	
	delta_pep[0] = 6.23409e-3 * holiang + 4.03363e-1;
      }       
      pep_real[0] = momL[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV      
      par_ep[0] = pep_real[0]; 
      // RHRS angle and momentum calculation      
      th1[0] = calcf2t_th(Theta_R, R_XFP, R_XpFP, R_YFP, R_YpFP,   z_av[0]);
      th1[0] = th1[0]*Xptr + Xptm;
     
      ph1[0] = calcf2t_ph(PHI_R, R_XFP, R_XpFP, R_YFP, R_YpFP, z_av[0]);
      ph1[0] = ph1[0]*Yptr + Yptm;     
      
      momR[0] =  calcf2t_mom(mom_R, R_XFP, R_XpFP, R_YFP, R_YpFP,  z_av[0]);
      momR[0] = momR[0]*Momr+Momm;     
      
      par_k[1] = -th1[0]; /// DEc4 2019 2 lines
      par_k[2] = -ph1[0];
      double holiang1;
 
      // target struggling step #11	
      if(z_av_1[0]<8.0e-2){ 
	holiang1=  par_k[2] - hrs_ang;
	delta_pk[0] =-1.31749 * sin(-4.61513* holiang1) + 2.03687;
	
      } 
      else{
 	holiang1=  par_k[2] - hrs_ang;
	delta_pk[0] = 3.158e-2 * holiang1 + 4.05819e-1;	
      }
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; // kaon momentum at the reaction point      
      par_k[0] = pk_real[0];
      
      // missing mass calculation==============================      
      hallap = hallap - 0.1843 ;// must be -ve
      hallap = hallap/1000.0; // MeV-->GeV
     
      mm = CalcMM(hallap, par_ep, par_k, mp);     
      mm = (mm)*1000.; // MeV--->GeV
      
      h->Fill(mm);      

      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;      
      
      R_XFP  = R_XFP*XFPr +XFPm ; 
      R_XpFP = R_XpFP*XpFPr+XpFPm;
      R_YFP  = R_YFP*YFPr+ YFPm;
      R_YpFP = R_YpFP*YpFPr +YpFPm; 
      
      z_av[0] =z_av[0]*Ztr + Ztm;
    
      tnew->Fill();
      
      bool lambdaflag=false;  
      int peak_with_hit= -1; 
      for(int j=0; j<npeak; j++){
	if(Lambda_cent[j]-Lambda_width[j]<mm
	   &&mm < Lambda_cent[j]+Lambda_width[j]){	 
	  
	  lambdaflag=true;
	  peak_with_hit=j; 
	  h_1 ->Fill(mm);
	  h_1 ->SetLineColor(j+2);	  
	}
	else lambdaflag=false;  
	if(ntune_event<nmax && lambdaflag==true){
	  foil_flag[ntune_event] = peak_with_hit;	  
	  p10[ntune_event]  = par_ep[0];
	  p11[ntune_event]  = par_ep[1];
	  p12[ntune_event]  = par_ep[2];
	  p13[ntune_event]  = par_k[0];
	  p14[ntune_event]  = par_k[1];
	  p15[ntune_event]  = par_k[2];
	  p16[ntune_event]  = hallap;
	    
	  x[ntune_event]  = R_XFP; ////RHRS  open  these lines only when RHRS mom matrix is tuning
	  y[ntune_event]  = R_YFP;
	  xp[ntune_event] = R_XpFP;
	  yp[ntune_event] = R_YpFP;

	  // x[ntune_event]  = XFP; ////LHRS open these line only when LHRS mom matrix is tuning
	  // y[ntune_event]  = YFP;
	  // xp[ntune_event] = XpFP;
	  // yp[ntune_event] = YpFP;

	  z_recon[ntune_event] = z_av_1[0];
	  phir[ntune_event] =ph1[0];
	  phil[ntune_event] =ph2[0];	 
	  ntune_event++;	    
	}	  
      }//int j	
    }    
  }  
  tnew->Write();
  // ((((((((((((((((((((((((((((((((((((((((( t2 ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  for (int i=0 ; i< ent_2 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_2[j]  = -2222.0;
      l_th_fp_2[j] = -2222.0; 
      l_y_fp_2[j]  = -2222.0;
      l_ph_fp_2[j] = -2222.0;
      th1_2[j] = -2222.0;
      th2_2[j] = -2222.0;
      ph1_2[j] =-2222.0;
      ph2_2[j] =-2222.0;    
     
      delta_pep_2[j]= -2222.0;
      pep_real_2[j] =-2222.0;
      delta_pk_2[j]= -2222.0;
      pk_real_2[j] = -2222.0;

      delta_ht_pep_2[j]= -2222.0; //May 01, 2021
      pep_ht_real_2[j] =-2222.0;
      delta_ht_pk_2[j]= -2222.0;
      pk_ht_real_2[j] = -2222.0;      
     
      r_x_fp_2[j]  = -2222.0;
      r_th_fp_2[j] = -2222.0;
      r_y_fp_2[j]  = -2222.0;
      r_ph_fp_2[j] = -2222.0;     
      
      trig5_2[j] = 0.0;
      rtrig_2 = false;
    }   
    trig5_2[0] = 0.0;
    rtrig_2 = false;   
    t2->GetEntry(i); 
   
    if(trig5_2[0]>1.0) rtrig_2 = true; //JUly 01, 2019
    else rtrig_2 = false;

    z_av_2[0] = (lvz_2[0] + rvz_2[0])/2.0;
    z_av_1_2[0] =  z_av_2[0];
    z_av_ht_2[0] =  z_av_2[0];     
   
    XFP_2   = l_x_fp_2[0];
    XpFP_2  = l_th_fp_2[0];
    YFP_2   = l_y_fp_2[0];
    YpFP_2  = l_ph_fp_2[0];
    
    R_XFP_2   = r_x_fp_2[0];
    R_XpFP_2  = r_th_fp_2[0];
    R_YFP_2   = r_y_fp_2[0];
    R_YpFP_2  = r_ph_fp_2[0];    
       
    XFP_2  = (XFP_2-XFPm)/XFPr;
    XpFP_2 = (XpFP_2-XpFPm)/XpFPr;
    YFP_2  = (YFP_2-YFPm)/YFPr;
    YpFP_2 = (YpFP_2-YpFPm)/YpFPr;
    
    R_XFP_2  = (R_XFP_2-XFPm)/XFPr; 
    R_XpFP_2 = (R_XpFP_2-XpFPm)/XpFPr;
    R_YFP_2  = (R_YFP_2-YFPm)/YFPr;
    R_YpFP_2 = (R_YpFP_2-YpFPm)/YpFPr;
    /// ================= for Al event in tune ==== jan 07, 2020
    XFP_4  = XFP_2; 
    XpFP_4 = XpFP_2;
    YFP_4  = YFP_2; 
    YpFP_4 = YpFP_2;

    R_XFP_4  =R_XFP_2;
    R_XpFP_4 = R_XpFP_2;
    R_YFP_4  =R_YFP_2;
    R_YpFP_4 =R_YpFP_2;
    
    z_av_2[0] =(z_av_2[0]- Ztm)/Ztr;
    
    th2_2[0] = calcf2t_th(Theta_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_1_2[0]);
    th2_2[0] = th2_2[0]*Xptr + Xptm; 
    
    ph2_2[0] = calcf2t_ph(PHI_L, XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_1_2[0] );
    ph2_2[0] = ph2_2[0]*Yptr + Yptm;       
    
    momL_2[0] =  calcf2t_mom(mom_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2[0]); 
    momL_2[0] =momL_2[0]*Momr + Momm; 
    
    momL_2[0] = momL_2[0]*2.21784/2.10; //  Original for H/T and tritium only     
    
    par_ep_2[1] = -th2_2[0]; // Dec4, 2019
    par_ep_2[2] = -ph2_2[0];    
    ///// From here is the energy struggling for the H data in the T kinematics======LHRS===**************************
    par_ht_ep_2[1] = -th2_2[0]; //May 01, 2021
    par_ht_ep_2[2] = -ph2_2[0];

    double holiang2_ht;
    // Target struggling LHRS step #7
    if( z_av_ht_2[0]<8.0e-2){
      holiang2_ht = par_ht_ep_2[2] + hrs_ang;	
      holiang2_ht = - holiang2_ht;
      
      delta_ht_pep_2[0] = -1.35758*sin(-4.59571* holiang2_ht) + 2.09093;      
    } 
    else{
      holiang2_ht = par_ht_ep_2[2] + hrs_ang;
      holiang2_ht = - holiang2_ht;
      delta_ht_pep_2[0] = 6.23409e-3* holiang2_ht + 4.03363e-1;      
    }     
    pep_ht_real_2[0] = momL_2[0] + delta_ht_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
    par_ht_ep_2[0] = pep_ht_real_2[0] ;     
    ///// Upto here is the energy struggling for the H data in the T kinematics======LHRS==**************************
    ///// From here is the energy struggling for the Al data in the H/T kinematics==========******* LHRS ***************
    double holiang2;
    // Target struggling LHRS step #7
    if( z_av_1_2[0]<8.0e-2){
      holiang2 = par_ep_2[2] + hrs_ang;	
      holiang2 = - holiang2;
      
      delta_pep_2[0] = -1.35758*sin(-4.59571* holiang2) + 2.09093;
      delta_pep_2[0] = delta_pep_2[0] + 0.063; //May 07 2021***************************           
    } 
    else{
      holiang2 = par_ep_2[2] + hrs_ang;
      holiang2 = - holiang2;
      // delta_pep_2[0] = 6.23409e-3* holiang2 + 4.03363e-1;
      delta_pep_2[0] =  0.3004; //May 07 2021************************************
    }     
    pep_real_2[0] = momL_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
    par_ep_2[0] = pep_real_2[0] ; 
    ///// upto here is the energy struggling for the Al data in the T kinematics=====LHRS =====**************************    
    // RHRS angle and momentum calculation      
    th1_2[0] = calcf2t_th(Theta_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,   z_av_2[0]);
    th1_2[0] = th1_2[0]*Xptr + Xptm;
    
    ph1_2[0] = calcf2t_ph(PHI_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2, z_av_2[0]);
    ph1_2[0] = ph1_2[0]*Yptr + Yptm;    
    
    momR_2[0] =  calcf2t_mom(mom_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,  z_av_2[0]);
    momR_2[0] = momR_2[0]*Momr+Momm;    
    
    par_k_2[1] = -th1_2[0]; // Dec2, 2019
    par_k_2[2] = -ph1_2[0];
    
    ///// From here is the energy struggling for the H data in the T kinematics  for RHRS=====**************************
    par_ht_k_2[1] = -th1_2[0]; // May 01, 2021
    par_ht_k_2[2] = -ph1_2[0];

    double holiang3_ht;
    // target struggling step #11	
    if(z_av_ht_2[0]<8.0e-2){ 
      holiang3_ht = par_ht_k_2[2] - hrs_ang;
      holiang3_ht = holiang3_ht;
      delta_ht_pk_2[0] =-1.31749*sin(-4.61513*holiang3_ht) + 2.03687;
    } 
    else{
      holiang3_ht = par_ht_k_2[2] - hrs_ang;
      holiang3_ht = holiang3_ht;
      delta_ht_pk_2[0] = 3.158e-2*holiang3_ht + 4.05819e-1;      
    }
    pk_ht_real_2[0] = momR_2[0] + delta_ht_pk_2[0]/1000.0; // kaon momentum at the reaction point     
    par_ht_k_2[0] = pk_ht_real_2[0];
    
    hallap_ht_2 = hallap_2-0.1843 ;// must be -ve  for the H data in the tritium kinematics
    hallap_ht_2 = hallap_ht_2/1000.0; // MeV-->GeV
    
    // missing mass calculation==============================   
    mm_2 = CalcMM(hallap_ht_2, par_ht_ep_2, par_ht_k_2, mp);     
    mm_2 = (mm_2)*1000.; // MeV--->GeV

    mmT_T = CalcMM(hallap_ht_2, par_ht_ep_2, par_ht_k_2, m_T);// to analyze H/T  and considering H mass as tritium mass  
    mmT_T = (mmT_T)*1000.; // MeV--->GeV  // m_T is the tritium target
    mmT_T = mmT_T - 2994.814;
    ///// Upto here is the energy struggling for the H data in the T kinematics  for RHRS=====**************************
    ///// From here is the energy struggling for the Al data in the H/T kinematics  for RHRS=====**************************
    double holiang3;
    // target struggling step #11	
    if(z_av_1_2[0]<8.0e-2){ 
      holiang3 = par_k_2[2] - hrs_ang;
      holiang3 = holiang3;
      delta_pk_2[0] =-1.31749*sin(-4.61513*holiang3) + 2.03687; 
      delta_pk_2[0] =delta_pk_2[0] + 0.0627;//May 07 2021************************************      
    } 
    else{
      holiang3 = par_k_2[2] - hrs_ang;
      holiang3 = holiang3;
      //  delta_pk_2[0] = 3.158e-2*holiang3 + 4.05819e-1;      
      delta_pk_2[0] = 0.2962;//May 07 2021************************************ 
    }
    pk_real_2[0] = momR_2[0] + delta_pk_2[0]/1000.0; // kaon momentum at the reaction point     
    par_k_2[0] = pk_real_2[0];
    // missing mass calculation==============================    
    if(z_av_1_2[0]<8.0e-2){  //May 07 2021************************************
      hallap_2 = hallap_2 - 0.1175;
    }
    else{
      hallap_2 = hallap_2 - 0.2257;//May 07 2021************************************
    }    
    hallap_2 = hallap_2/1000.0; // MeV-->GeV
   
    mm_4 = CalcMM(hallap_2, par_ep_2, par_k_2, m_Al); 
    mm_4 = (mm_4)*1000.0;
    mm_ht = mm_4 - 25.3123*1000;   
    //// MAy 17, 2021.... Now adjusting the Al entrance and Exit window.          
    bool HTflag = false;   
    if(rtrig_2==true &&  fabs(ctime_2)<1.0  && fabs(lvz_2[0]-rvz_2[0])<0.04){ // Oct 28 2020 closed the ct 0.04 is original
      HTflag = true;
    }
    else  HTflag = false;    
       
    if(HTflag == true  && fabs(z_av_1_2[0])<0.10){
      
      h_2->Fill(mm_2);
               
      XFP_2 = XFP_2 * XFPr + XFPm;
      XpFP_2 = XpFP_2 * XpFPr + XpFPm;
      YFP_2 = YFP_2 * YFPr + YFPm;
      YpFP_2 = YpFP_2 * YpFPr + YpFPm;
      
      
      R_XFP_2  = R_XFP_2*XFPr +XFPm ; 
      R_XpFP_2 = R_XpFP_2*XpFPr+XpFPm;
      R_YFP_2  = R_YFP_2*YFPr+ YFPm;
      R_YpFP_2 = R_YpFP_2*YpFPr +YpFPm; 
      
      z_av_2[0] =z_av_2[0]*Ztr + Ztm;
      
      tnew->Fill();
      
      bool lambdaflag_2=false;  
      int peak_with_hit_2= -1; 
      for(int j=0; j<npeak_2; j++){
	if(Lambda_cent_2[j]-Lambda_width_2[j]<mm_2
	   &&mm_2 < Lambda_cent_2[j]+Lambda_width_2[j]){	 
	  
	  lambdaflag_2=true;
	  peak_with_hit_2=j; 	  
	}
	else lambdaflag_2=false;  

	if(ntune_event_2<nmax_2 && lambdaflag_2==true){
	  foil_flag_2[ntune_event_2] = peak_with_hit_2;	  
	  p10_2[ntune_event_2]  = par_ht_ep_2[0];  /// for H data in T kinematics
	  p11_2[ntune_event_2]  = par_ht_ep_2[1]; // ht added on may 05, 2021
	  p12_2[ntune_event_2]  = par_ht_ep_2[2];
	  p13_2[ntune_event_2]  = par_ht_k_2[0];
	  p14_2[ntune_event_2]  = par_ht_k_2[1];
	  p15_2[ntune_event_2]  = par_ht_k_2[2];
	  p16_2[ntune_event_2]  = hallap_ht_2; // may 05
	    
	  x_2[ntune_event_2]  = R_XFP_2; ////RHRS open these line only when RHRS mom matrix is tuning
	  y_2[ntune_event_2]  = R_YFP_2;
	  xp_2[ntune_event_2] = R_XpFP_2;
	  yp_2[ntune_event_2] = R_YpFP_2;

	  // x_2[ntune_event_2]  = XFP_2; ////LHRS open these line only when LHRS mom matrix is tuning
	  // y_2[ntune_event_2]  = YFP_2;
	  // xp_2[ntune_event_2] = XpFP_2;
	  // yp_2[ntune_event_2] = YpFP_2;

	  z_recon_2[ntune_event_2] = z_av_1_2[0];
	  phir_2[ntune_event_2] =ph1_2[0];
	  phil_2[ntune_event_2] =ph2_2[0];	 
	  ntune_event_2++;	    
	}	  
      }//int j	
    }
    /// ==================== to include the Al/HT evenys for tune========== jan 07 2020
    bool htflag1 = false;   
    if(rtrig_2==true &&  fabs(ctime_2)<1.0  && fabs(lvz_2[0]-rvz_2[0])<0.040){ 
      htflag1 = true;
    }
    else  htflag1 = false;     
    bool HTallflag = false;
    if(a1_2<120.0 && a2_2>1650.0 && a2_2<6800.0 &&  //for entrance and exit aluminum caps  
       ((z_av_1_2[0] > -0.14 && z_av_1_2[0]< -0.11) || (z_av_1_2[0] > 0.11 && z_av_1_2[0]< 0.14))){
      HTallflag = true;
    }
    else HTallflag = false;
    ///// AL background analysis
    bool hbgflag = false;
    if(rtrig_2==true &&((ctime_2>-49.39 && ctime_2 < -9.06)||(ctime_2>13.18 && ctime_2 < 48.6))&& fabs(lvz_2[0]-rvz_2[0])<0.040){ 
      hbgflag = true;
    }
    else  hbgflag = false;

    if(HTallflag == true && hbgflag == true)
      {
	hb-> Fill(mm_ht);
	hbal_15_1-> Fill(mm_ht); // c/0.15
	
	hb5-> Fill(mm_ht);
	hb_al1-> Fill(mm_ht); // c/0.25
	hbal_20_1-> Fill(mm_ht); // c/0.2
      }
    ///// Real spectrum
    if(htflag1 == true  && HTallflag == true/* && (mm_ht> -5.0 && mm_ht <-1.0) ||(mm_ht> 19.5 && mm_ht <32.4))*/){ //plot spectrum
      
      h21->Fill(mm_ht);      
      hh->Fill(mm_ht);
      hal_15_1->Fill(mm_ht); // NOV 30, 2020// counts /0.15 MeV
      
      he->Fill(mm_ht); 
      h_al1->Fill(mm_ht); // NOV 16, 2020// counts /0.25 MeV
      hal_20_1->Fill(mm_ht); // NOV 30, 2020// counts /0.2 MeV     

      XFP_4 = XFP_4 * XFPr + XFPm;
      XpFP_4 = XpFP_4 * XpFPr + XpFPm;
      YFP_4 = YFP_4 * YFPr + YFPm;
      YpFP_4 = YpFP_4 * YpFPr + YpFPm;      
      
      R_XFP_4  = R_XFP_4*XFPr +XFPm ; 
      R_XpFP_4 = R_XpFP_4*XpFPr+XpFPm;
      R_YFP_4  = R_YFP_4*YFPr+ YFPm;
      R_YpFP_4 = R_YpFP_4*YpFPr +YpFPm;    
      
      tnew->Fill(); 
      bool lambdaflag_4=false;  
      int peak_with_hit_4= -1; 
      for(int j=0; j<npeak_4; j++){
	if(Lambda_cent_4[j]-Lambda_width_4[j]<mm_ht
	   &&mm_ht < Lambda_cent_4[j]+Lambda_width_4[j]){	 
	  
	  lambdaflag_4=true;
	  peak_with_hit_4=j; 
	}
	else lambdaflag_4=false;  

	if(ntune_event_4<nmax_4 && lambdaflag_4==true){
	  foil_flag_4[ntune_event_4] = peak_with_hit_4;
	  
	  p10_4[ntune_event_4]  = par_ep_2[0]; // right side should be _2 // for Al data
	  p11_4[ntune_event_4]  = par_ep_2[1];
	  p12_4[ntune_event_4]  = par_ep_2[2];
	  p13_4[ntune_event_4]  = par_k_2[0];
	  p14_4[ntune_event_4]  = par_k_2[1];
	  p15_4[ntune_event_4]  = par_k_2[2];
	  p16_4[ntune_event_4]  = hallap_2;
	    
	  x_4[ntune_event_4]  = R_XFP_4; ////RHRS open these lines only when RHRS mom matrix is tuning
	  y_4[ntune_event_4]  = R_YFP_4;
	  xp_4[ntune_event_4] = R_XpFP_4;
	  yp_4[ntune_event_4] = R_YpFP_4;

	  // x_4[ntune_event_4]  = XFP_4; ////LHRS open these lines only when LHRS mom matrix is tuning
	  // y_4[ntune_event_4]  = YFP_4;
	  // xp_4[ntune_event_4] = XpFP_4;
	  // yp_4[ntune_event_4] = YpFP_4;

	  z_recon_4[ntune_event_4] = z_av_1_2[0];
	  phir_4[ntune_event_4] =ph1_2[0];
	  phil_4[ntune_event_4] =ph2_2[0];
	  ntune_event_4++;	   
	}
      }      
    }//  j=0; j< npeak_4
  }  
  tnew->Write();  
  // ))))))))))))))))))))))))))))))))))))))))) t2 ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((( t3 ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  bool rtrig_3 = false; 
  for(int i=0; i<nmax_3; i++){
    x_3[i]    = -2222.0; 
    y_3[i]    = -2222.0; 
    xp_3[i]   = -2222.0;
    yp_3[i]   = -2222.0;
    z_av_3[i] = -2222.0;
    z_av_1_3[i] = -2222.0;
    phir_3[i] = -2222.0;
    phil_3[i] = -2222.0;
    z_recon_3[i] = -2222.0;  //// Jan 04, 2019  
    
    foil_flag_3[i] = -1;
  }  
  // ///////////////////////////////////////////////////////////////////////////////////
  for (int i=0 ; i< ent_3 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_3[j]  = -2222.0;
      l_th_fp_3[j] = -2222.0; 
      l_y_fp_3[j]  = -2222.0;
      l_ph_fp_3[j] = -2222.0;
      th1_3[j] = -2222.0;
      th2_3[j] = -2222.0;
      ph1_3[j] =-2222.0;
      ph2_3[j] =-2222.0;    
      
      delta_pep_3[j]= -2222.0;
      pep_real_3[j] =-2222.0;
      delta_pk_3[j]= -2222.0;
      pk_real_3[j] = -2222.0;

      delta_tt_pep_3[j]= -2222.0; // May 01, 2021
      pep_tt_real_3[j] =-2222.0;
      delta_tt_pk_3[j]= -2222.0;
      pk_tt_real_3[j] = -2222.0;      
      
      r_x_fp_3[j]  = -2222.0;
      r_th_fp_3[j] = -2222.0;
      r_y_fp_3[j]  = -2222.0;
      r_ph_fp_3[j] = -2222.0;        
      
      trig5_3[j] = 0.0;
      rtrig_3 = false;
    }
    
    trig5_3[0] = 0.0;
    rtrig_3 = false;
    
    t3->GetEntry(i);

    if(trig5_3[0]>1.0) rtrig_3 = true; //JUly 01, 2019
    else rtrig_3 = false;
    
    z_av_3[0] = (lvz_3[0] + rvz_3[0])/2.0;
    z_av_1_3[0] =  z_av_3[0];
    z_av_tt_3[0] =  z_av_3[0];
    
    XFP_3   = l_x_fp_3[0];
    XpFP_3  = l_th_fp_3[0];
    YFP_3   = l_y_fp_3[0];
    YpFP_3  = l_ph_fp_3[0];
    
    R_XFP_3   = r_x_fp_3[0];
    R_XpFP_3  = r_th_fp_3[0];
    R_YFP_3   = r_y_fp_3[0];
    R_YpFP_3  = r_ph_fp_3[0];    
   
    XFP_3  = (XFP_3-XFPm)/XFPr;
    XpFP_3 = (XpFP_3-XpFPm)/XpFPr;
    YFP_3  = (YFP_3-YFPm)/YFPr;
    YpFP_3 = (YpFP_3-YpFPm)/YpFPr;
      
    R_XFP_3  = (R_XFP_3-XFPm)/XFPr; 
    R_XpFP_3 = (R_XpFP_3-XpFPm)/XpFPr;
    R_YFP_3  = (R_YFP_3-YFPm)/YFPr;
    R_YpFP_3 = (R_YpFP_3-YpFPm)/YpFPr;
      
    z_av_3[0] =(z_av_3[0]- Ztm)/Ztr;
      
    th2_3[0] = calcf2t_th(Theta_L, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_1_3[0]);
    th2_3[0] = th2_3[0]*Xptr + Xptm;     
    ph2_3[0] = calcf2t_ph(PHI_L, XFP_3, XpFP_3, YFP_3, YpFP_3, z_av_1_3[0] );
    ph2_3[0] = ph2_3[0]*Yptr + Yptm; 
     
    momL_3[0] =  calcf2t_mom(mom_L, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_3[0]); 
    momL_3[0] =momL_3[0]*Momr + Momm;     
    momL_3[0] = momL_3[0]*2.21784/2.10; //  Original for H/T and tritium only     

    par_ep_3[1] = -th2_3[0]; // Dec4, 2019
    par_ep_3[2] = -ph2_3[0];
    /// May 01, 2021, ************************* from here "seperate target struggling for NNL"  ***********========================
    /// making the different variable to store the variable to calculate the NNL spectrum.====LHRS =========================
    par_tt_ep_3[1] = -th2_3[0]; // May 01, 2021
    par_tt_ep_3[2] = -ph2_3[0];
    
    double holiang5_tt;
    // Target struggling LHRS step #7...... MAy 01, 2021
    if( z_av_tt_3[0]<8.0e-2){
      holiang5_tt = par_tt_ep_3[2] + hrs_ang;	
      holiang5_tt = - holiang5_tt;
      
      delta_tt_pep_3[0] = -1.35758*sin(-4.59571* holiang5_tt) + 2.09093;      
    } 
    else{
      holiang5_tt = par_tt_ep_3[2] + hrs_ang;
      holiang5_tt = - holiang5_tt;
      delta_tt_pep_3[0] = 6.23409e-3* holiang5_tt + 4.03363e-1;      
    }     
    pep_tt_real_3[0] = momL_3[0] + delta_tt_pep_3[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
    par_tt_ep_3[0] = pep_tt_real_3[0] ;        
    /// May 01, 2021, ************************* upto here "seperate target struggling for NNL"  ***********======================    
    ///// The following are the target energy struggling for the Al (T/T) target====LHRS =======================
    double holiang5;
    if( z_av_1_3[0]<8.0e-2){
      holiang5 = par_ep_3[2] + hrs_ang;	
      holiang5 = - holiang5;
      
      delta_pep_3[0] = -1.35758*sin(-4.59571* holiang5) + 2.09093;
      delta_pep_3[0] =  delta_pep_3[0] + 0.0524; //May 03, 2021********************************************      
    } 
    else{
      holiang5 = par_ep_3[2] + hrs_ang;
      holiang5 = - holiang5;
      // delta_pep_3[0] = 6.23409e-3* holiang5 + 4.03363e-1;
      delta_pep_3[0] =  0.3027;	//May 07, 2021********************************************
    }       
    pep_real_3[0] = momL_3[0] + delta_pep_3[0]/1000.0; //LHRS  momentum at the reaction point in GeV      
    par_ep_3[0] = pep_real_3[0] ; 
    ///// up to here is the target energy struggling for the Al target===LHRS ========================= 
    // RHRS angle and momentum calculation      
    th1_3[0] = calcf2t_th(Theta_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3,   z_av_3[0]);
    th1_3[0] = th1_3[0]*Xptr + Xptm;
    
    ph1_3[0] = calcf2t_ph(PHI_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3, z_av_3[0]);
    ph1_3[0] = ph1_3[0]*Yptr + Yptm;   
           
    momR_3[0] =  calcf2t_mom(mom_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3,  z_av_3[0]);
    momR_3[0] = momR_3[0]*Momr+Momm;      

    par_k_3[1] = -th1_3[0]; // Dec2, 2019
    par_k_3[2] = -ph1_3[0];
    /// May 01, 2021, ************************* from here "seperate target struggling for NNL"  ***********=====================
    /// making the different variable to store the variable to calculate the NNL spectrum.==RHRS ================================
    par_tt_k_3[1] = -th1_3[0]; // May 01, 2021
    par_tt_k_3[2] = -ph1_3[0];

    double holiang6_tt;
    // target struggling step #11	
    if(z_av_tt_3[0]<8.0e-2){ 
      holiang6_tt = par_tt_k_3[2] - hrs_ang;
      holiang6_tt = holiang6_tt;
      delta_tt_pk_3[0] =-1.31749*sin(-4.61513*holiang6_tt) + 2.03687;	
    } 
    else{
      holiang6_tt = par_tt_k_3[2] - hrs_ang;
      holiang6_tt = holiang6_tt;
      delta_tt_pk_3[0] = 3.158e-2*holiang6_tt + 4.05819e-1;
    }
    pk_tt_real_3[0] = momR_3[0] + delta_tt_pk_3[0]/1000.0; // kaon momentum at the reaction point      
    par_tt_k_3[0] = pk_tt_real_3[0];
    /// May 01, 2021, ************************* upto here "seperate target struggling for NNL"  ***********===============
    //// Now MM calculation for the Lnn spectrum===========================================================

    hallap_tt_3 = hallap_3-0.1843 ;// must be -ve
    hallap_tt_3 = hallap_tt_3/1000.0; // MeV-->GeV

    mm_h = CalcMM(hallap_tt_3, par_tt_ep_3, par_tt_k_3, mp); //// to see hydrogen in tritium data
    mm_h = (mm_h)*1000.;
    mm_h =  mm_h -1115.683;

    mm_3 = CalcMM(hallap_tt_3, par_tt_ep_3, par_tt_k_3, m_T); 
    mm_3 = (mm_3)*1000.; // GeV--->MeV
    mm_t = mm_3 -2994.814; // for tritium target only By TOSHI when consider the tritium mass (recoil mass)
    //// up to here is  MM calculation for the Lnn spectrum===========================================================
    
    //// The following is target eenrgy struggling and 27_Mg_L calculation===RHRS===Al target(T/T)=========
    double holiang6;
    // target struggling step #11	
    if(z_av_1_3[0]<8.0e-2){ 
      holiang6 = par_k_3[2] - hrs_ang;
      holiang6 = holiang6;
      delta_pk_3[0] =-1.31749*sin(-4.61513*holiang6) + 2.03687;
      delta_pk_3[0] =	delta_pk_3[0] + 0.0512;//May 07, 2021********************************************
    } 
    else{
      holiang6 = par_k_3[2] - hrs_ang;
      holiang6 = holiang6;
      //  delta_pk_3[0] = 3.158e-2*holiang6 + 4.05819e-1;
      delta_pk_3[0] =  0.2993;//May 07, 2021********************************************
    }
    pk_real_3[0] = momR_3[0] + delta_pk_3[0]/1000.0; // kaon momentum at the reaction point       
     
    par_k_3[0] = pk_real_3[0];

    // missing mass calculation====Al data ==========================
    // hallap_3 = hallap_3-0.1843 ;// must be -ve
    if(z_av_1_3[0]<8.0e-2){
      hallap_3 = hallap_3 - 0.1065;//May 07, 2021********************************************
    }
    else{
      hallap_3 = hallap_3 - 0.2318;//May 07, 2021********************************************
    }    
    hallap_3 = hallap_3/1000.0; // MeV-->GeV
   
    mm_Al = CalcMM(hallap_3, par_ep_3, par_k_3, m_Al);
    mm_Al = (mm_Al)*1000.0; 
    mm_Al1 = mm_Al -25.3123*1000; // for Al  Al kinematics only. when consider Al as target  
    /////+++++++++++++++++++++++++++++++++++===========================- aerogel hist for thesis sept 04, 2020++++++++++++++++  
    //// ========================= from here is tritium data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    bool Tritium_flag = false;
    if(rtrig_3==true && fabs(lvz_3[0]-rvz_3[0])< 0.053 &&   
       a1_3<160.0 && a2_3>1685.0 && a2_3<8000.0 && fabs(z_av_1_3[0])<0.10){//chnged  Jan 11 2020  
      Tritium_flag = true;
    }
    else Tritium_flag = false;
    //// for nnl real analysis
    if(Tritium_flag == true && fabs(ctime_3)<1.0){  //to plot real nnL spectrum 
      
      h50->Fill(mm_t); // For nnL spectrum      
      h25_2->Fill(mm_t);   
     
      // H contamination +++++++++++++++++++++++++++++
      h_h->Fill(mm_h); // H Contamination with 1.5 MeV/bin
      h_hc0->Fill(mm_h);
      h_hc1->Fill(mm_h);       
      h1_2->Fill(mm_t); // nnL spectrum
      ///&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$# DEC 01, 2020%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      h150_1->Fill(mm_t);
      h150_2->Fill(mm_t);
      h150_3->Fill(mm_t);
      h150_4->Fill(mm_t);
      h150_5->Fill(mm_t);
      h150_6->Fill(mm_t);
 
      h125_1->Fill(mm_t);
      h125_2->Fill(mm_t);
      h125_3->Fill(mm_t);
      h125_4->Fill(mm_t);
      h125_5->Fill(mm_t);
      h125_6->Fill(mm_t);
      h1_qu->Fill(mm_t);
      h1_q1->Fill(mm_t);
     
      h1_0->Fill(mm_t);
      h1_01->Fill(mm_t);
      h1_02->Fill(mm_t);
      h1_03->Fill(mm_t);
      h1_04->Fill(mm_t);
      h1_05->Fill(mm_t);
      h1_may->Fill(mm_t);
    }    
    //// for nnl backkground analysis
    bool bg_nnlflag = false;
    if((ctime_3>-49.39 && ctime_3 < -9.06)||(ctime_3> 13.18 && ctime_3 < 48.6)){//chnged  Jan 11 2020  
      bg_nnlflag = true;
    }
    else bg_nnlflag = false;
   
    if(Tritium_flag == true && bg_nnlflag ==true){
         
      h25_2b->Fill(mm_t);     
      h50b->Fill(mm_t);
      // background for H contamination ++++++++++++++++++++++++
      h1_2b->Fill(mm_t);
      h_hbg->Fill(mm_h); // bg in H contamination
      h_hc0b->Fill(mm_h);
      h_hc1b->Fill(mm_h);    

      h1_01b->Fill(mm_t);
      h1_03b->Fill(mm_t);
      h1_05b->Fill(mm_t);
      h125_6b->Fill(mm_t);
      h125_4b->Fill(mm_t);
      h125_5b->Fill(mm_t);
      h125_3b->Fill(mm_t);
      h1b_may->Fill(mm_t);
      
      h150_1b->Fill(mm_t);
      h150_2b->Fill(mm_t);
      h1_qub->Fill(mm_t);
      h1_qb1->Fill(mm_t);
     
      h1_02b->Fill(mm_t);     
    }
    ////// ========================= upt o here is tritium data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    // //// ========================= from here is Al data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    ////// for tune a1 nd a2 = 120 1650, 6800& 0.04 and to see the spectrum a1 , a2 = 160, 1585, 8000 & 0.5   
    bool TTallflag = false;
    if(rtrig_3==true && fabs(lvz_3[0]-rvz_3[0])<0.040 && a1_3<120.0
       && a2_3>1650.0 && a2_3<6800.0 &&
       ((z_av_1_3[0]>-0.14 && z_av_1_3[0]<-0.11) || (z_av_1_3[0]>0.11  && z_av_1_3[0]<0.14))){//for entrance and exit aluminum caps  
      TTallflag = true;
    }
    else TTallflag = false;
    ///// for Al background analysis ======================================================================
    bool bg1_flag = false;
    if((ctime_3> -49.39 && ctime_3 < -9.06) ||(ctime_3> 13.18 && ctime_3 < 48.6)){
      bg1_flag = true;
    }
    else bg1_flag = false;   

    if(TTallflag == true && bg1_flag  == true)
      {
	hb1->Fill(mm_Al1);
	hb6->Fill(mm_Al1);

	h21_b->Fill(mm_Al1);
	hb_al2->Fill(mm_Al1); // count/0.25
	hbal_20_2->Fill(mm_Al1);// count/0.2
	hbal_15_2->Fill(mm_Al1);// count/0.15
      } 
    /// real spectrum AL 
    if(TTallflag == true && fabs(ctime_3) < 1.0 /*&&((mm_Al1 > -5.3 &&mm_Al1<-1.1)||(mm_Al1 >19.5 &&mm_Al1<32.4))*/){//plot Al spectrum
     
      h20->Fill(mm_Al1);
      ht->Fill(mm_Al1);
      hal_15_2->Fill(mm_Al1); // NOV 30 2020  counts/0.15 MeV 
      
      hd->Fill(mm_Al1);
      h_al2->Fill(mm_Al1); // NOV 16 2020  counts/0.25 MeV
      hal_20_2->Fill(mm_Al1); // NOV 30 2020  counts/0.2 MeV     
      
      XFP_3 = XFP_3 * XFPr + XFPm;
      XpFP_3 = XpFP_3 * XpFPr + XpFPm;
      YFP_3 = YFP_3 * YFPr + YFPm;
      YpFP_3 = YpFP_3 * YpFPr + YpFPm;      
      
      R_XFP_3  = R_XFP_3*XFPr +XFPm ; 
      R_XpFP_3 = R_XpFP_3*XpFPr+XpFPm;
      R_YFP_3  = R_YFP_3*YFPr+ YFPm;
      R_YpFP_3 = R_YpFP_3*YpFPr +YpFPm; 
      
      z_av_3[0] =z_av_3[0]*Ztr + Ztm;
      
      tnew->Fill();
         
      bool lambdaflag_3=false;  // need adjustment for tune Jan_02
      int peak_with_hit_3= -1; 
      for(int j=0; j<npeak_3; j++){
  	if(Lambda_cent_3[j]-Lambda_width_3[j]<mm_Al1  // mm_3 need to be adjusted for event selection
  	   &&mm_Al1 < Lambda_cent_3[j]+Lambda_width_3[j]){	 
	  
  	  lambdaflag_3=true;
  	  peak_with_hit_3=j; 	  
  	}
  	else lambdaflag_3=false;  
	
  	if(ntune_event_3<nmax_3 && lambdaflag_3==true){
  	  foil_flag_3[ntune_event_3] = peak_with_hit_3;
	  
  	  p10_3[ntune_event_3]  = par_ep_3[0];
  	  p11_3[ntune_event_3]  = par_ep_3[1];
  	  p12_3[ntune_event_3]  = par_ep_3[2];
  	  p13_3[ntune_event_3]  = par_k_3[0];
  	  p14_3[ntune_event_3]  = par_k_3[1];
  	  p15_3[ntune_event_3]  = par_k_3[2];
  	  p16_3[ntune_event_3]  = hallap_3;
	    
  	  x_3[ntune_event_3]  = R_XFP_3; ////RHRS open these lines only when RHRS momentum matrix  is tuning
  	  y_3[ntune_event_3]  = R_YFP_3;
  	  xp_3[ntune_event_3] = R_XpFP_3;
  	  yp_3[ntune_event_3] = R_YpFP_3;

  	  // x_3[ntune_event_3]  = XFP_3; ////LHRS open these lines only when LHRS momentum matrix  is tuning
  	  // y_3[ntune_event_3]  = YFP_3;
  	  // xp_3[ntune_event_3] = XpFP_3;
  	  // yp_3[ntune_event_3] = YpFP_3;

  	  z_recon_3[ntune_event_3] = z_av_1_3[0];
  	  phir_3[ntune_event_3] =ph1_3[0];
  	  phil_3[ntune_event_3] =ph2_3[0];
	 
  	  ntune_event_3++;	    
  	}
	  
      }//int j	
    }    
  }
  tnew->Write();  
  // ))))))))))))))))))))))))))))))))))))))))) t3 )))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))  
  // // HH data
 
  TF1 *f1 = new TF1("f1","gaus",1112.65,1118.74);////1112.74,1118.64 --->With Al and 1112.38,1118.93
  TF1 *f2 = new TF1("f2","gaus",1189.95,1195.15);//1189.75,1195.15 --->With Al and with out al 1189.71,1195.42
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h->Draw();
  
  h->Fit("f1","MR+");  /// these 3 lines are closed on June 09, 2020 to remove the Gaussian fit
  h->Fit("f2","MR+"); /// during the paper is already spread to Hall A collaboration
  f1->Draw("same"); 
    
  TLatex l;
  l.SetTextSize(0.06);
  l.SetTextFont(62);// time bold italic
  l.DrawLatex(1126,200,Form("#Lambda"));  
  l.DrawLatex(1126,179,Form("#color[2]{Mean = %.6g MeV}",f1->GetParameter(1)));  // closed on June 09, 2020
  l.DrawLatex(1126,158,Form("#color[2]{    #sigma = %.4g MeV}",f1->GetParameter(2))); // while getting parameter from fit
  l.DrawLatex(1164,116,Form("#Sigma^{0}")); 
  l.DrawLatex(1164,100,Form("#color[2]{Mean = %.6g MeV}",f2->GetParameter(1))); // closed on June 09, 2020
  l.DrawLatex(1164,79,Form("#color[2]{    #sigma = %.4g MeV}",f2->GetParameter(2))); // while getting parameter from fit
  ///// H in tritium data  H contamination ++++++++++++++++++++++++++++++++++++++++++++++++=
  h_hbg->Scale(1.0/38.0);
  TCanvas *c_h = new TCanvas("c_h","c_h", 600,600);
  c_h->cd();
  h_h->Draw();
  h_hbg->Draw("E2 same");
  h_hbg->SetFillStyle(3002);
  h_hbg->SetMarkerStyle(28);
  h_hbg->SetMarkerColor(kGreen);  
  
  h_hc0b->Scale(1.0/38.0);
  TCanvas *ch_0 = new TCanvas("ch_0","ch_0", 600,600);
  ch_0->cd();
  h_hc0->Draw();
  h_hc0b->Draw("E2 same");
  h_hc0b->SetFillStyle(3002);
  h_hc0b->SetMarkerStyle(28);
  h_hc0b->SetMarkerColor(kGreen);
  
  h_hc1b->Scale(1.0/38.0);
  TCanvas *ch_1 = new TCanvas("ch_1","ch_1", 600,600);
  ch_1->cd();
  h_hc1->Draw();
  h_hc1b->Draw("E2 same");
  h_hc1b->SetFillStyle(3002);
  h_hc1b->SetMarkerStyle(28);
  h_hc1b->SetMarkerColor(kGreen);
  ////// up to here is tritium data for H kinematics to see the H contamination in the T gas(((((()()()()()()()()()()()()()(()()()(
  
  ////  For H data with T kinematics
  TF1 *f1_2 = new TF1("f1_2","gaus",1112.52,1118.8);//1112.7,1118.76  june 15, 2021
  TCanvas* c2_2 = new TCanvas("c2_2","c2_2",600,600);
  c2_2->cd();
  h_2->Draw();  
  h_2->Fit("f1_2","MR+");  
  TLatex l2;
  l2.SetTextSize(0.06);
  l2.SetTextFont(62);
  l2.DrawLatex(1140,80,Form("#Lambda"));  // when use the fitting 
  l2.DrawLatex(1140,60,Form("#color[2]{#sigma = %.6g}",f1_2->GetParameter(2)));
  l2.DrawLatex(1140,70,Form("#color[2]{mean = %.6g}",f1_2->GetParameter(1)));  
  ////0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000------------888888888  
  TF1 *ft2 = new TF1("ft2","gaus",-4.5686,-2.1805);
  TF1 *f3 = new TF1("f3","gaus",19.4155,22.883);
  TF1 *f4 = new TF1("f4","gaus",29.4463,31.8078);

  HT_b ->Add(h20_b,h21_b,1.0,1.0);
  HT_b->Scale(1.0/38.0);
  TCanvas* cT = new TCanvas("cT","cT",600,600);
  cT->cd();
  HT->Add(h20,h21,1.0,1.0);
  HT->Draw();
  
  HT_b->Draw("E2 same"); // background
  HT_b->SetFillStyle(3002);
  HT_b->SetMarkerStyle(28);
  HT_b->SetMarkerColor(kGreen); 
  ////// counts /0.5   =========================================================================================
  hb2->Add(hb,hb1,1.0,1.0); 
  hb2->Scale(1.0/38.0);  
  TF1 *f10 = new TF1("f10","gaus",-4.6245,-3.0746);
  TF1 *f30 = new TF1("f30","gaus",29.183,31.3214);
 
  TCanvas* c11 = new TCanvas("c11","c11",600,600);
  c11->cd();
  htt->Add(hh,ht,1.0,1.0);
  htt->Draw();

  hb2->Draw("E2 same"); // background
  hb2->SetFillStyle(3002);
  hb2->SetMarkerStyle(28);
  hb2->SetMarkerColor(kGreen); 
  // /// counts / 0.25 MeV Al spectrum +++++++++++++++++++++++++ 
  hb_al->Add(hb_al1,hb_al2,1.0,1.0);
  hb_al->Scale(1.0/38.0);
  TCanvas* c_al = new TCanvas("c_al","c_al",600,600);
  c_al->cd();
  h_al->Add(h_al1,h_al2,1.0,1.0);
  h_al->Draw(); 
  
  hb_al->Draw("E2 same"); // background
  hb_al->SetFillStyle(3002);
  hb_al->SetMarkerStyle(28);
  hb_al->SetMarkerColor(kGreen); 

  hbal_20->Add(hbal_20_1,hbal_20_2,1.0,1.0);//NOV 30, 2020
  hbal_20->Scale(1.0/38.0);
  TCanvas* cal_20 = new TCanvas("cal_20","cal_20",600,600);
  cal_20->cd();
  hal_20->Add(hal_20_1,hal_20_2,1.0,1.0);
  hal_20->Draw();
  
  hbal_20->Draw("E2 same"); // background
  hbal_20->SetFillStyle(3002);
  hbal_20->SetMarkerStyle(28);
  hbal_20->SetMarkerColor(kGreen);
  ////0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 
  hbal_15->Add(hbal_15_1,hbal_15_2,1.0,1.0);  //NOV 30, 2020
  hbal_15->Scale(1.0/38.0);
  TCanvas* cal_15 = new TCanvas("cal_15","cal_15",600,600);
  cal_15->cd();
  hal_15->Add(hal_15_1,hal_15_2,1.0,1.0);
  hal_15->Draw();

  hbal_15->Draw("E2 same"); // background
  hbal_15->SetFillStyle(3002);
  hbal_15->SetMarkerStyle(28);
  hbal_15->SetMarkerColor(kGreen); 
 
  TF1 *f6 = new TF1("f6","gaus",-4.93,-2.582);
  TF1 *f5 = new TF1("f5","gaus",4.499,6.8395);  
  TF1 *f7 = new TF1("f7","gaus",19.4075,20.6485);
  TF1 *f8 = new TF1("f8","gaus",29.9542,31.1404);
  hb7->Add(hb5,hb6,1.0,1.0);
  hb7->Scale(1.0/38.0);  
  TCanvas* cf = new TCanvas("cf","cf",600,600);
  cf->cd();
  hf->Add(hd,he,1.0,1.0);
  hf->Draw();
  f5->SetLineWidth(1);
  f6->SetLineWidth(1);
  f7->SetLineWidth(1);
  f8->SetLineWidth(1);

  hb7->Draw("E2 same"); // background
  hb7->SetFillStyle(3002);
  hb7->SetMarkerStyle(28);
  hb7->SetMarkerColor(kGreen);    
  ////000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000---------------------8888888
  h50b ->Scale(1.0/38.0);
  TH1F * h50b1 = (TH1F*)h50b->Clone(); 
  TCanvas* c50 = new TCanvas("c50","c50",600,600);
  c50->cd();
  h50->Draw();
  h50b1->Draw("E2 same");
  h50b1->SetFillStyle(3002);
  h50b1->SetMarkerStyle(28);
  h50b1->SetMarkerColor(kGreen);  
  
  h1_qb1 ->Scale(1.0/38.0);
  TCanvas* c_q1 = new TCanvas("c_q1","c_q1",600,600);
  c_q1->cd();
  h1_q1->Draw(); 
  h1_qb1->Draw("E2 same");
  h1_qb1->SetFillStyle(3002);
  h1_qb1->SetMarkerStyle(28);
  h1_qb1->SetMarkerColor(kGreen);  
  //////// up to here is the quasi frre shape calculation /////////////////////////////////////////////////////////// 
  h1_2b->Scale(1.0/38.0);
  TCanvas *c1_2 = new TCanvas("c1_2","c1_2", 600,600);
  c1_2->cd();
  h1_2->Draw();    
  h1_2b->Draw("E2 same");
  
  h1_2b->SetFillStyle(3002);
  h1_2b->SetMarkerStyle(28);
  h1_2b->SetMarkerColor(kGreen);
  //////////&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%$$$$$$$$$$$ Dec 01, 2020 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  h1_01b->Scale(1.0/38.0);
  TCanvas *c1_01 = new TCanvas("c1_01","c1_01", 600,600);
  c1_01->cd();
  h1_01->Draw();
  h1_01b->Draw("E2 same");
  h1_01b->SetFillStyle(3002);
  h1_01b->SetMarkerStyle(28);
  h1_01b->SetMarkerColor(kGreen);

  h1b_may->Scale(1.0/38.0);
  TCanvas *c1_may = new TCanvas("c1_may","c1_may", 600,600);
  c1_may->cd();
  h1_may->Draw();
  h1b_may->Draw("E2 same");
  h1b_may->SetFillStyle(3002);
  h1b_may->SetMarkerStyle(28);
  h1b_may->SetMarkerColor(kGreen);
  
  
  h1_02b->Scale(1.0/38.0);
  TCanvas *c1_02 = new TCanvas("c1_02","c1_02", 600,600);
  c1_02->cd();
  h1_02->Draw();
  h1_02b->Draw("E2 same");
  h1_02b->SetFillStyle(3002);
  h1_02b->SetMarkerStyle(28);
  h1_02b->SetMarkerColor(kGreen);  

  h1_03b->Scale(1.0/38.0);
  TCanvas *c1_03 = new TCanvas("c1_03","c1_03", 600,600);
  c1_03->cd();
  h1_03->Draw();
  h1_03b->Draw("E2 same");
  h1_03b->SetFillStyle(3002);
  h1_03b->SetMarkerStyle(28);
  h1_03b->SetMarkerColor(kGreen);
 
  h125_6b->Scale(1.0/38.0);
  TCanvas *c125_6 = new TCanvas("c125_6","c125_6", 600,600);
  c125_6->cd();
  h125_6->Draw();
  h125_6b->Draw("E2 same");
  h125_6b->SetFillStyle(3002);
  h125_6b->SetMarkerStyle(28);
  h125_6b->SetMarkerColor(kGreen); 
  //////////&&&&&&&&&&&&&&&&&&&&&&%%%%%%%%%%$$$$$$$$$$$ Dec 01, 2020 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  h25_2b->Scale(1.0/38.0);
  TCanvas *c25_2 = new TCanvas("c25_2","c25_2", 600,600);
  c25_2->cd();
  h25_2->Draw();    
  h25_2b->Draw("E2 same");
  
  h25_2b->SetFillStyle(3002);
  h25_2b->SetMarkerStyle(28);
  h25_2b->SetMarkerColor(kGreen);  
  
  TFile* f_new = new TFile("./output_root/July12_12pm.root","recreate"); // also open the  f_new->Close(); aslo
  
  h->Write();
  h_2->Write();  
  hf->Write();
  hb7->Write();
  hal_20->Write();
  hbal_20->Write();
  HT->Write();
  HT_b->Write();
  hal_15->Write();
  hbal_15->Write();

  htt->Write();
  hb2->Write();
  
  h_al->Write();
  hb_al->Write();  

  h125_6->Write();
  h125_6b->Write();
  
  h1_q1->Write();
  h1_qb1->Write();
  
  h25_2->Write();
  h25_2b->Write();

  h1_03->Write();
  h1_03b->Write(); 

  h1_2->Write();  
  h1_2b->Write();

  h1_01->Write();
  h1_01b->Write();

  h1_may->Write();
  h1b_may->Write(); 
  
  f_new->Close(); 
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  const int nite =0;
  double temp[nite]; 
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  for(int i=0 ; i<nite ; i++){
    x[i] = i+1;
    // temp[i] = tune(momL_opt,i); /// open when LHRS mom matrix is tuned
    temp[i] = tune(momR_opt,i);  /// open when RHRS mom matrix is tuned
    
    sprintf(tempc, "./MOM_MATRICES/RMOM5_July8_1st_%d.dat",i); // output matrix
    ofstream * ofs = new ofstream(tempc); 
    int nppp = 0;
    const int nn = 5;  // 5 for 5th order 4 for 4th order jan 31
    for(int i=0; i<nn+1; i++){
      for(int e=0; e<nn+1; e++){
	for(int d=0; d<nn+1; d++){ 
	  for(int c=0; c<nn+1; c++){
	    for(int b=0; b<nn+1; b++){
	      for(int a=0; a<nn+1; a++){  
		if(a+b+c+d+e==i){
		  *ofs <<momR_opt[nppp] // MomR_opt[] for RHRS tune and  MomL_opt[] for LHRS tune
		       << " " << a 
		       << " " << b
		       << " " << c
		       << " " << d
		       << " " << e << endl;
		  nppp++; 
		  
		}
	      }
	    }
	  }
	}
      }
    }
    ofs->close();
    ofs->clear();
    
    cout << temp[i]<<endl; 
  }      
  if(nite>0){
    TGraph * gr = new TGraph(nite,x,temp);  
    TCanvas * c4 = new TCanvas("c4","",600,600); 
    gr->Draw("*la"); 
  }
} //end of  main function

//////////////////////////////////////////////////
double calcf2t_th(double* P, double xf, double xpf, 
		  double yf, double ypf,double zt)
//////////////////////////////////////////////////
{  // -----4th order -----   
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n =   
  return Y;
}
// ////////////////////////////////////////////////
//////////////////////////////////////////////////
double calcf2t_ph(double* P, double xf, double xpf, 
		  double yf, double ypf, double zt)
//////////////////////////////////////////////////
{  // -----4th order -----   
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
  
  return Y; 
}
//////////////////////////////////////////////////
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
//////////////////////////////////////////////////
{  // -----5th order -----   
  const int nMatT=5;  
  const int nXf=5;
  const int nXpf=5;
  const int nYf=5;
  const int nYpf=5;
  const int nZt=5;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
  
  return Y; 
}
// missing mass function definition====================
double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt){ // Dec 3,2019
  
  double pe = ee;
  double Ee = sqrt(me*me + pe*pe);
  Ee = Ee - 0.000111; // GeV /// 0.000163;
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = pvec_ep[0];
  double xpep = pvec_ep[1];
  double ypep = pvec_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  vec_ep.RotateX(hrs_ang);
  double Eep = sqrt(pep*pep + me*me);
  
 
  double pk  = pvec_k[0];
  double xpk = pvec_k[1];
  double ypk = pvec_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
  vec_k.RotateX(-hrs_ang);
  double Ek = sqrt(pk*pk + mk*mk);  
 
  double missingE2 = 0.0, missingP2 = 0.0, missingM2 = 0.0;
  missingE2 = pow(Ee + mt - Ek - Eep, 2.0);
  missingP2 = (vec_e - vec_ep - vec_k) * (vec_e - vec_ep - vec_k);
  missingM2 = missingE2 - missingP2;
  
  double MissingMass = 0.0;
  MissingMass = sqrt(missingM2);

  return MissingMass;  
}
//############### up to hear missing mass #####################
// #############################################################
double tune(double* pa, int j) // tune fun defn
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam =Mom_Par; // for momentum tune jan 31
 
  TMinuit* minuit = new TMinuit(allparam); 
  minuit->SetFCN(fcn); // very imp function setying for chi square 
    
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1); 
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
 
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    //  step[i] = 1.0e-3; /// Original 
    step[i] = 2.0*0.5;  // for rough matrix, when matrix is far from reality
    
    LLim[i] = pa[i] -10; // pa[i]*0.8; // KI
    ULim[i] = pa[i] + 10; //pa[i]*0.8; // KI
    // LLim[i] = pa[i] - pa[i]*0.8; // KI
    // ULim[i] = pa[i] + pa[i]*0.8; // KI
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  //  arglist[0] = 2.0; // was active before
  arglist[0] = 1.0;  // KI
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  // ~~~~ Migrad + Simplex  ~~~~ one of the way to get optimized parameter
  arglist[0] = 20000;
  arglist[1] = 0.01; // To make more presise
  minuit -> mnexcm("MINImize",arglist,2,ierflg); 
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat); 
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){  
     
    // / minuit -> GetParameter(i,momL_opt[i],er); // open this line only when LHRS momentum matrix is tuned
    minuit -> GetParameter(i,momR_opt[i],er); // open this line only when RHRS momentum matrix is tuned   
  }  
  return chi2; 
}
// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################  
{
  double chi2 = 0.0;
  double chi_12 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  const double sigma = 0.0045; 
  double ref_mm = 0.0; 
  double residual = 0.0;

  double par_ep[3];
  double par_k[3];
  double halla_p;
  double momr[100];
  double moml[100];
  double z_av;
  double z_av_sc;
  double rvz;
  double ph1;
  double th1;
  double ph2;
  double th2;
  double  delta_pk[100];
  double delta_pep[100];
  double pk_real[100];
  double pep_real[100];
  double MM;
  double THL;
  double PHL;
  double THR;
  double PHR;
  // (((((((((((((((((((((((((((((((((((((((( t2 (((((((((((((((((((((((((((((((((
 
  double chi_22 = 0.0;
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double ref_mm_2 = 0.0; 
  double residual_2 = 0.0;

  double par_ep_2[3];
  double par_k_2[3];
  double halla_p_2;
  double momr_2[100];
  double moml_2[100];
  double z_av_2;
  double z_av_sc_2;
  double ph1_2;
  double th1_2;
  double ph2_2;
  double th2_2;
  double delta_pk_2[100];
  double delta_pep_2[100];
  double pk_real_2[100];
  double pep_real_2[100];
  double MM_2;
  double THL_2;
  double PHL_2;
  double THR_2;
  double PHR_2;
  //)))))))))))))))))))))))))))))))))))))))))))))   t2  ))))))))))))))))))))))))))))  
  // (((((((((((((((((((((((((((((((((((((((( t3 (((((((((((((((((((((((((((((((((
  double chi_33 = 0.0;
  double XFP_3, XpFP_3;
  double YFP_3, YpFP_3;
  double ref_mm_3 = 0.0; 
  double residual_3 = 0.0;

  double par_ep_3[3];
  double par_k_3[3];
  double halla_p_3;
  double momr_3[100];
  double moml_3[100];
  double z_av_3;
  double z_av_sc_3;

  double ph1_3;
  double th1_3;
  double ph2_3;
  double th2_3;
  double delta_pk_3[100];
  double delta_pep_3[100];
  double pk_real_3[100];
  double pep_real_3[100];
  double MM_3;
  double THL_3;
  double PHL_3;
  double THR_3;
  double PHR_3;
  //)))))))))))))))))))))))))))))))))))))))))))))   t3  ))))))))))))))))))))))))))))
  // (((((((((((((((((((((((((((((((((((((((( to include the Al/HT data for tune Jan 07, 2020 (((((((((((((((((((((((((((((((((
  double chi_44 = 0.0;
  double XFP_4, XpFP_4;
  double YFP_4, YpFP_4;
  double ref_mm_4 = 0.0; 
  double residual_4 = 0.0;

  double par_ep_4[3];
  double par_k_4[3];
  double halla_p_4;
  double momr_4[100];
  double moml_4[100];
  double z_av_4;
  double z_av_sc_4;

  double ph1_4;
  double th1_4;
  double ph2_4;
  double th2_4;
  double delta_pk_4[100];
  double delta_pep_4[100];
  double pk_real_4[100];
  double pep_real_4[100];
  double MM_4;
  double THL_4;
  double PHL_4;
  double THR_4;
  double PHR_4;
  //)))))))))))))))))))))))))))))))))))))))))))))   up here is _4  ))))))))))))))))))))))))))))
  for(int i=0; i<ntune_event; i++){ 
    residual = 0.0;
    ref_mm = 0.0; 
    ref_mm  = Lambda_real[foil_flag[i]];    
    ref_mm = ref_mm/1000.0;
    
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    z_av = z_recon[i]; 
    ph1 = phir[i];  // open when calibrate the Momentum 
    ph1= -ph1;
    ph2 = phil[i];    
    ph2 = -ph2; 

    XFP   =(XFP -XFPm)/XFPr;  
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    z_av_sc = (z_av - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml[0] =  calcf2t_mom(param, XFP, XpFP, YFP, YpFP,  z_av_sc);
    moml[0] = moml[0]*Momr+Momm; 
    double hang;
    if( z_av<8.0e-2){
      hang = ph2 + hrs_ang;
      hang = -hang;
      delta_pep[0] = -1.35758*sin(-4.59571* hang) + 2.09093;      
    } 
    else{
      hang = ph2 + hrs_ang;
      hang = -hang;
      delta_pep[0] = 6.23409e-3*hang + 4.03363e-1;
    } 
    
    pep_real[0] = moml[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV  
    
    //  par_ep[0] = pep_real[0];// When LHRS momentum  tuned keep this line open otherwise comment it
    par_ep[0] = p10[i];// comment this line if the previous line is open or not commented
    par_ep[1] = p11[i]; 
    par_ep[2] = p12[i];   
     
    //  for the RHRS momentum tunning
    momr[0] =  calcf2t_mom(param, XFP, XpFP, YFP, YpFP,  z_av_sc);
    momr[0] = momr[0]*Momr+Momm;
    double hang1;    

    if(z_av<8.0e-2){
      hang1= ph1 - hrs_ang;
      delta_pk[0] =-1.31749*sin(-4.61513* hang1) + 2.03687;      
    } 
    else{
      hang1= ph1 - hrs_ang;
      delta_pk[0] = 3.158e-2*hang1 + 4.05819e-1; 
    }
    pk_real[0] = momr[0] + delta_pk[0]/1000.0; 

    par_k[0] = pk_real[0];// open this line only when RHRS momentum matrix tuned and the next line closed
    //  par_k[0] = p13[i];   // if previous line is open (not commented), comment this line 
    par_k[1] = p14[i]; // jan 31    
    par_k[2] = p15[i];

    halla_p = p16[i];
    MM = CalcMM(halla_p, par_ep, par_k, mp);    
    residual = MM-ref_mm;
   
    //   chi_12 = chi_12 + pow(residual,2.0);
    ///////  if need to use the sigma statistical weigh
    if(foil_flag[i] ==0)
      {chi_12 = chi_12 +6*pow(residual,2.0);}
    else
      {chi_12 = chi_12 +10*pow(residual,2.0);}    
  }
  // (((((((((((((((((((((((((((((((((((((((( t2 ((((((((((((((((((((((((((((((((( 
  for(int i=0; i<ntune_event_2; i++){ 
    residual_2 = 0.0;
    ref_mm_2 = 0.0; 
    ref_mm_2  = Lambda_real_2[foil_flag_2[i]];    
    ref_mm_2 = ref_mm_2/1000.0;
    
    XFP_2   = x_2[i];
    XpFP_2  = xp_2[i];
    YFP_2   = y_2[i];
    YpFP_2  = yp_2[i];
    z_av_2 = z_recon_2[i]; 
    ph1_2 = phir_2[i];  // open when calibrate the Momentum 
    ph1_2 = - ph1_2;
    ph2_2 = phil_2[i];    
    ph2_2 = - ph2_2;

    XFP_2   =(XFP_2 -XFPm)/XFPr;  
    XpFP_2  =(XpFP_2-XpFPm)/XpFPr;
    YFP_2   =(YFP_2 -YFPm)/YFPr;
    YpFP_2  =(YpFP_2-YpFPm)/YpFPr;
    z_av_sc_2 = (z_av_2 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_2[0] =  calcf2t_mom(param, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_sc_2);
    moml_2[0] = moml_2[0]*Momr+Momm;

    moml_2[0] = moml_2[0]*2.21784/2.1; //for H/T and tritium only
    double hang2;   
    if( z_av_2<8.0e-2){

      hang2 =  ph2_2 + hrs_ang;
      hang2 = -hang2;
      delta_pep_2[0] = -1.35758*sin(-4.59571* hang2) + 2.09093;
     
    } 
    else{
      hang2 =  ph2_2 + hrs_ang;
      hang2 = -hang2;
      delta_pep_2[0] = 6.23409e-3*hang2 + 4.03363e-1;
    }    
    pep_real_2[0] = moml_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV     
       
    //  par_ep_2[0] = pep_real_2[0];// Wwhen LHRS momentum  tuned, uncoment this line
    par_ep_2[0] = p10_2[i];
    par_ep_2[1] = p11_2[i];
    par_ep_2[2] = p12_2[i];   
   
    //  for the RHRS momentum tunning
    momr_2[0] =  calcf2t_mom(param, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_sc_2);
    momr_2[0] = momr_2[0]*Momr+Momm;
    double hang3;    
    if(z_av_2<8.0e-2){

      hang3 = ph1_2 - hrs_ang;
      delta_pk_2[0] =-1.31749*sin(-4.61513* hang3) + 2.03687;     
    } 
    else{
      hang3 = ph1_2 - hrs_ang;
      delta_pk_2[0] = 3.158e-2*hang3 + 4.05819e-1; 
    }
    pk_real_2[0] = momr_2[0] + delta_pk_2[0]/1000.0;     
    
    par_k_2[0] = pk_real_2[0];// when RHRS matrix tuned
    //  par_k_2[0] = p13_2[i];    
    par_k_2[1] = p14_2[i];   // jan 31 
    par_k_2[2] = p15_2[i]; 
       
    halla_p_2 = p16_2[i];
    MM_2 = CalcMM(halla_p_2, par_ep_2, par_k_2, mp);    
    residual_2 = MM_2-ref_mm_2;
   
    chi_22 = chi_22 +25*pow(residual_2,2.0);
    // chi_22 = chi_22 +pow(residual_2,2.0);
  }
  //)))))))))))))))))))))))))))))))))))))))))))))   t2  ))))))))))))))))))))))))))))
  // (((((((((((((((((((((((((((((((((((((((( t3  T/T data for Al ((((((((((((((((((((((((((((((
  for(int i=0; i<ntune_event_3; i++){ 
    residual_3 = 0.0;
    ref_mm_3 = 0.0; 
    ref_mm_3  = Lambda_real_3[foil_flag_3[i]];    
    ref_mm_3 = ref_mm_3/1000.0;
    
    XFP_3   = x_3[i];
    XpFP_3  = xp_3[i];
    YFP_3   = y_3[i];
    YpFP_3  = yp_3[i];
    z_av_3 = z_recon_3[i]; 
    ph1_3 = phir_3[i];  // open when calibrate the Momentum 
    ph1_3 = - ph1_3;
    ph2_3 = phil_3[i];    
    ph2_3 = - ph2_3;

    XFP_3   =(XFP_3 -XFPm)/XFPr;  
    XpFP_3  =(XpFP_3-XpFPm)/XpFPr;
    YFP_3   =(YFP_3 -YFPm)/YFPr;
    YpFP_3  =(YpFP_3-YpFPm)/YpFPr;
    z_av_sc_3 = (z_av_3 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_3[0] =  calcf2t_mom(param, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_sc_3);
    moml_3[0] = moml_3[0]*Momr+Momm;

    moml_3[0] = moml_3[0]*2.21784/2.1; //for H/T and tritium only
    double hang9;   
    if( z_av_3<8.0e-2){

      hang9 =  ph2_3 + hrs_ang;
      hang9 = -hang9;
      delta_pep_3[0] = -1.35758*sin(-4.59571* hang9) + 2.09093;
      delta_pep_3[0] = delta_pep_3[0] + 0.0524; /// may 07, 2021************
    } 
    else{
      hang9 =  ph2_3 + hrs_ang;
      hang9 = -hang9;
      //  delta_pep_3[0] = 6.23409e-3*hang9 + 4.03363e-1;
      delta_pep_3[0] = 0.3027;/// may 07, 2021************
    }    
    pep_real_3[0] = moml_3[0] + delta_pep_3[0]/1000.0; //LHRS  momentum at the reaction point in GeV 

    
    //  par_ep_3[0] = pep_real_3[0];// Wwhen LHRS momentum  tuned, uncoment this line
    par_ep_3[0] = p10_3[i];
    par_ep_3[1] = p11_3[i];
    par_ep_3[2] = p12_3[i];   
   
    //  for the RHRS momentum tunning
    momr_3[0] =  calcf2t_mom(param, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_sc_3);
    momr_3[0] = momr_3[0]*Momr+Momm;
    double hang10;    
    if(z_av_3<8.0e-2){
      hang10 = ph1_3 - hrs_ang;
      delta_pk_3[0] =-1.31749*sin(-4.61513* hang10) + 2.03687;
      delta_pk_3[0] = delta_pk_3[0] + 0.0512; /// may 07, 2021************
    } 
    else{
      hang10 = ph1_3 - hrs_ang;
      // delta_pk_3[0] = 3.158e-2*hang10 + 4.05819e-1;
      delta_pk_3[0] =   0.2993;  /// may 07, 2021************
    }
    pk_real_3[0] = momr_3[0] + delta_pk_3[0]/1000.0;   

    par_k_3[0] = pk_real_3[0];// when RHRS matrix tuned
    // par_k_3[0] = p13_3[i];    
    par_k_3[1] = p14_3[i];  // jan 31  
    par_k_3[2] = p15_3[i];   
    
    halla_p_3 = p16_3[i];
    MM_3 = CalcMM(halla_p_3, par_ep_3, par_k_3, m_Al);  //need to adjust  
    MM_3 = MM_3 -25.3123;   
    
    residual_3 = MM_3-ref_mm_3;    
    
    if(foil_flag_3[i]==0)
      {chi_33 = chi_33 +25*pow(residual_3,2.0);}
    
    else if(foil_flag_3[i]== 1)      
      {chi_33 = chi_33 +25*pow(residual_3,2.0);}
    
    else       
      {chi_33 = chi_33 +15*pow(residual_3,2.0);}    
    
    
    // chi_33 = chi_33 +9*pow(residual_3,2.0);/// colsed on May 24, 2021
  }
  // ((((((((((((((((((((((((4(((((((((((((((( _4 H/T data for Al  Jan 07, 2020((((((((((((((((((((((((((((((
  for(int i=0; i<ntune_event_4; i++){ 
    residual_4 = 0.0;
    ref_mm_4 = 0.0; 
    ref_mm_4  = Lambda_real_4[foil_flag_4[i]];    
    ref_mm_4 = ref_mm_4/1000.0;
    
    XFP_4   = x_4[i];
    XpFP_4  = xp_4[i];
    YFP_4   = y_4[i];
    YpFP_4  = yp_4[i];
    z_av_4 = z_recon_4[i]; 
    ph1_4 = phir_4[i];  // open when calibrate the Momentum 
    ph1_4 = - ph1_4;
    ph2_4 = phil_4[i];    
    ph2_4 = - ph2_4;

    XFP_4   =(XFP_4 -XFPm)/XFPr;  
    XpFP_4  =(XpFP_4-XpFPm)/XpFPr;
    YFP_4   =(YFP_4 -YFPm)/YFPr;
    YpFP_4  =(YpFP_4-YpFPm)/YpFPr;
    z_av_sc_4 = (z_av_4 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_4[0] =  calcf2t_mom(param, XFP_4, XpFP_4, YFP_4, YpFP_4,  z_av_sc_4);
    moml_4[0] = moml_4[0]*Momr+Momm;

    moml_4[0] = moml_4[0]*2.21784/2.1; //for H/T and tritium only
    double hang9;   
    if(z_av_4<8.0e-2){
      hang9 =  ph2_4 + hrs_ang;
      hang9 = -hang9;
      delta_pep_4[0] = -1.35758*sin(-4.59571* hang9) + 2.09093;
      delta_pep_4[0] = delta_pep_4[0] + 0.0637;/// may 07, 2021************
    } 
    else{
      hang9 =  ph2_4 + hrs_ang;
      hang9 = -hang9;
      // delta_pep_4[0] = 6.23409e-3*hang9 + 4.03363e-1;
      delta_pep_4[0] =   0.3004;/// may 07, 2021************
    }    
    pep_real_4[0] = moml_4[0] + delta_pep_4[0]/1000.0; //LHRS  momentum at the reaction point in GeV 
    
    //  par_ep_4[0] = pep_real_4[0];// Wwhen LHRS momentum  tuned, uncomment this line
    par_ep_4[0] = p10_4[i];
    par_ep_4[1] = p11_4[i];
    par_ep_4[2] = p12_4[i];   
   
    ////  for the RHRS momentum tunning
    momr_4[0] =  calcf2t_mom(param, XFP_4, XpFP_4, YFP_4, YpFP_4,  z_av_sc_4);
    momr_4[0] = momr_4[0]*Momr+Momm;
    double hang10;    
    if(z_av_4<8.0e-2){
      hang10 = ph1_4 - hrs_ang;
      delta_pk_4[0] =-1.31749*sin(-4.61513* hang10) + 2.03687;
      delta_pk_4[0] = delta_pk_4[0] + 0.0627;/// may 07, 2021************
    } 
    else{
      hang10 = ph1_4 - hrs_ang;
      // delta_pk_4[0] = 3.158e-2*hang10 + 4.05819e-1;
      delta_pk_4[0] = 0.2962;/// may 07, 2021************
    }
    pk_real_4[0] = momr_4[0] + delta_pk_4[0]/1000.0;     

    par_k_4[0] = pk_real_4[0];// when RHRS matrix tuned un coment this line
    //  par_k_4[0] = p13_4[i];    
    par_k_4[1] = p14_4[i];  // jan 31  
    par_k_4[2] = p15_4[i];   

    halla_p_4 = p16_4[i];
    MM_4 = CalcMM(halla_p_4, par_ep_4, par_k_4, m_Al);  //need to adjust  
    MM_4 = MM_4 -25.3123;
     
    residual_4 = MM_4-ref_mm_4;    
    
    if(foil_flag_4[i]==0)
      {chi_44 = chi_44 +25*pow(residual_4,2.0);}
    
    else if(foil_flag_4[i]==1)
      {chi_44 = chi_44 +25*pow(residual_4,2.0);}
    
    else
      {chi_44 = chi_44 +15*pow(residual_4,2.0);}    
    
    // chi_44 = chi_44 +9*pow(residual_4,2.0);/// colsed on May 24, 2021
  }
  //)))))))))))))))))))))))))))))))))))))))))))))   _4 Jan 07, 2020  ))))))))))))))))))))))))))))
  //)))))))))))))))))))))))))))))))))))))))))))))  Jan 12, 2020  ))))))))))))))))))))))))))))
  chi2 = chi_12 +chi_22 + chi_33 + chi_44;
  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  fval = chi2;
}
