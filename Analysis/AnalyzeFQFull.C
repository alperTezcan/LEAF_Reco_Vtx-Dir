#include <iostream>
#include <stdio.h>     
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TMath.h>
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"
//#define ONLYGOODEVENTS
using namespace std;
double TankSize = 1100/2;//in cm, the maximal size of the plots for vertices etc...
double TankRadius = 742/2;//in cm, the maximal size of the plots for vertices etc...
double TankHalfHeight = 1042/2;//in cm, the maximal size of the plots for vertices etc...
double TotalNHits = 1e4;
bool containedOnly = true;
bool stopReconstructed = true;
bool dwallCut = false; double dwallCutValue = 30;//To remove interactions in PMTs
int MRTuneCut=false;

double scalar_product(double vrec[3], double vtrue[3], bool divide_by_norm = false){
  double prod_scal = vrec[0]*vtrue[0] + vrec[1]*vtrue[1] + vrec[2]*vtrue[2];
  double norm_rec = TMath::Sqrt( pow(vrec[0],2) + pow(vrec[1],2) + pow(vrec[2],2));
  double norm_true = TMath::Sqrt( pow(vtrue[0],2) + pow(vtrue[1],2) + pow(vtrue[2],2));
  if(divide_by_norm) return prod_scal/(norm_rec*norm_true);
  else return prod_scal;
}

double vector_difference(double vrec[3], double vtrue[3]){
  return TMath::Sqrt(pow(vrec[0]-vtrue[0],2)+pow(vrec[1]-vtrue[1],2)+pow(vrec[2]-vtrue[2],2));
}

bool outOfTank(double particlePos[3]){
    bool offTank = false;
    double Radius = TMath::Sqrt(pow(particlePos[0],2)+pow(particlePos[1],2));
    double HalfHeight = TMath::Abs(particlePos[2]);    
    if(Radius > TankRadius || HalfHeight > TankHalfHeight) offTank = true;

    return offTank;
}

double FindToWall(double vtx[3], double dir[3]){
  double stepSize = 10.;//in cm
  double trackPos[3]={0.};
  for(int i=0;i<3;i++) trackPos[i] = vtx[i];
  double trackRadius = TMath::Sqrt(pow(trackPos[0],2) + pow(trackPos[1],2));
  double towall=0;
  
  while( ( (trackRadius) < TankRadius) && (TMath::Abs(trackPos[2]) < TankHalfHeight) ){
    //cout << (trackRadius) << endl;
    for(int i=0;i<3;i++){
      trackPos[i] += dir[i]*stepSize;
    }
    trackRadius = TMath::Sqrt(pow(trackPos[0],2) + pow(trackPos[1],2));
    towall += stepSize;
  }

  return towall;
}

int DeterminePID(int lastdg){
  int q_tenthousands = (int (lastdg%10000));
  int div_thousands = (int (q_tenthousands/1000));
  int q_thousands = q_tenthousands%1000;
  int div_hundreds = (int (q_thousands/100));
  int q_hundreds = q_thousands%100;
  int div_tens = (int (q_hundreds/10));
  int q_tens = q_hundreds%10;
  //int q_units = lastdigits%1;
  //int pid = std::max(div_thousands-1.,0.)*pow(3,3) + std::max(div_hundreds-1.,0.)*pow(3,2) + std::max(div_tens-1.,0.)*pow(3,1) + std::max(q_tens-1.,0.);
  int pid = div_hundreds*pow(4,2) + div_tens*pow(4,1) + q_tens;
  std::cout<<"Last digits = "<<lastdg<<", "<<div_thousands<<" "<<div_hundreds<<" "<<div_tens<<" "<<q_tens<<std::endl; 
  return pid;
}
// Simple example of reading a generated Root file
int main(int argc, char **argv){

  char * Suffix="ola";
  int ParticleType = 0;  //int ParticleType=2;//Particle type: 0->gamma, 1->electron, 2->muon, 3->pion and then see the code.
  bool isHK=false;
  bool verbose=true;
  bool Input = false;
  char * filename = new char[256];
  char * outfilename = new char[256];
  int c=-1;
  while( (c = getopt(argc,argv,"s:p:vhf:o:")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
      switch(c){
      case 'f':
	Input=true;
	sprintf(filename,"%s",optarg);
	break;
      case 'o':
	sprintf(outfilename,"%s",optarg);
	break;
      case 's':
	Suffix = optarg;
	break;
      case 'p':
	ParticleType = atoi(optarg);
	break;
      case 'h':
	isHK = true;
	break;
      case 'v':
	verbose = true;
	break;
      }
  }
  if(isHK){
    TotalNHits = 1e5;
    TankSize = 7100/2;//in cm, the maximal size of the plots for vertices etc...
    //TankRadius = 7080/2;//in cm, the maximal size of the plots for vertices etc...
    //TankHalfHeight = 5480/2;//in cm, the maximal size of the plots for vertices etc...
    //nevt = std::min( 1e2 , ((double) WSTreeHits->GetEntries()) );
  }

  TFile * fOutput;
  TFile * fInput;
  //int c=-1;
  int pdgID[5] = {22, 11, 13, 211, 111 };
  double Mass[5]={0, 0.511, 105.658, 139.570, 134.9766};
  int Nfiles = 1;//2000;//333;//2*3332;//2000;//2900;//2000;//1000;//200;
  int maxEntries = 20000;
  TChain * FQTree = new TChain();
  TChain * WCSimTree; 
  TChain * GeoTree;// = (TTree*)file->Get("wcsimGeoT");
  int nevt = 0;
  bool hybrid=true;
  char inputFolder[512]="/disk01/usr5/bquilain/fiTQun/analyze/output_full";
  cout << "Folder is : " << inputFolder << endl;
  cout << "Suffix is : " << Suffix << ", particle type = " <<ParticleType<<", pdg = "<<pdgID[ParticleType]<<endl;

  if(Input){
    cout<<"File used = "<<filename<<", output file name = "<<outfilename<<endl;
    fInput = new TFile(filename,"read");
    FQTree = (TChain*) fInput->Get("fiTQun");
    WCSimTree = (TChain*) fInput->Get("wcsimT");
    GeoTree = (TChain*) fInput->Get("wcsimGeoT");
    nevt = std::min(2e4 , (double) FQTree->GetEntries());
    //nevt = FQTree->GetEntries();
    //nevt = maxEntries;//FQTree->GetEntries();
  }
  else{
    FQTree = new TChain("fiTQun","fiTQun");
    int nevtPerFile=5;
    ///disk01/usr5/bquilain/fiTQun/analyze/output/fiTQun_%s_%d.root
    
    for(int i=0;i<Nfiles;i++){
      cout<<"File open "<<i<<endl;
      //FQTree->Add(Form("%s/fiTQun_event%d_%s.root",inputFolder,i,Suffix));
      FQTree->Add(Form("%s/fiTQun_%s_%d.root",inputFolder,Suffix,i));
      //FQTree->Add(Form("%s/fiTQun_%s_full.root",inputFolder,Suffix,i));
      if(i==0) nevtPerFile=FQTree->GetEntries();
    }
    nevt = TMath::Min(Nfiles*nevtPerFile,maxEntries);//FQTree->GetEntries();
    cout<<"Number of events = "<<nevt<<endl;
    WCSimTree = new TChain("wcsimT","wcsimT");
    GeoTree = new TChain("wcsimGeoT","wcsimGeoT");
    //GeoTree->Add(Form("%s/fiTQun_event%d_%s.root",inputFolder,0,Suffix));
    //WCSimTree->Add(Form("%s/fiTQun_event%d_%s.root",inputFolder,0,Suffix));
    GeoTree->Add(Form("%s/fiTQun_%s_%d.root",inputFolder,Suffix,0));
    WCSimTree->Add(Form("%s/fiTQun_%s_%d.root",inputFolder,Suffix,0));
  }

  cout<<"Chaining is over, we will analyze "<<nevt<<" events"<<endl;


    //////////////////////////////////////
  char* wcsimdirenv;
  //wcsimdirenv = getenv ("WCSIMDIR");
  //if(wcsimdirenv !=  NULL){
  //gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  //}else{
  //gSystem->Load("../libWCSimRoot.so");
  //}
  // Get the number of events
  int nevent = WCSimTree->GetEntries();
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the WCSimTree in

  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the WCSimTree
  TBranch *branch = WCSimTree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  // Force deletion to prevent memory leak 
  WCSimTree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
  TBranch *branch2;
  if(hybrid){
    branch2 = WCSimTree->GetBranch("wcsimrootevent2");
    branch2->SetAddress(&wcsimrootsuperevent2);
  // Force deletion to prevent memory leak 
    WCSimTree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
  }



  // Geometry WCSimTree - only need 1 "event"
  WCSimRootGeom *geo = 0; 
  GeoTree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "GeoTree has " << GeoTree->GetEntries() << " entries" << std::endl;
  if (GeoTree->GetEntries() == 0) {
      exit(9);
  }
  GeoTree->GetEntry(0);
  TankRadius = geo->GetWCCylRadius();//in cm, the maximal size of the plots for vertices etc...
  TankHalfHeight = geo->GetWCCylLength()/2;//in cm, the maximal size of the plots for vertices etc...
  

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;

  int num_trig=0;


  
  //TCanvas * MultiPlots[7];
  //Dimensions of 7: 0->gamma, 1->electron, 2->muon, 3->pion and then see the code.
  //fq1rmom[7]: momentum of the particle assuming there is only one ring and assuming the PID hypothesis [i].
  //fq1rpos[7]: initial position of the particle assuming ...
  //fq1rdir[7]: initial direction of the particle assuming ...
  cout<<"Number of events: "<<nevt<<endl;
  const int nPID = 10;
  const int nsevt = 10;
  int evtId;//used to match WCSim event and fiTQun ones
  int fqnse;
  int fqnhitpmt[nsevt];
  float fqtotq[nsevt];
  int fq1rpcflg[nsevt][nPID];
  float fq1rmom[nsevt][nPID];
  float fq1rpos[nsevt][nPID][3];
  float fq1rdir[nsevt][nPID][3];
  float fq1rnll[nsevt][nPID];
  float fq1rtotmu[nsevt][nPID];
  float fqpi0nll[2];
  float fqpi0mass[2];
  float fqpi0totmu[2];
  float fqpi0mom1[2];
  float fqpi0mom2[2];
  float fqpi0momtot[2];
  float fqpi0pos[2][3];
  
  const int max_nfit=100;
  int fqnmrfit;//Number of MR fit results that are available
  int fqmrifit[max_nfit];//Fit ID of each MR fit result
  int fqmrnring[max_nfit];//Number of rings for this fit [1-6] fqmrpcflg[fqnmrfit]/I <0 if MINUIT did not converge during the fit
  int fqmrpcflg[max_nfit];//Number of rings for this fit [1-6] fqmrpcflg[fqnmrfit]/I <0 if MINUIT did not converge during the fit
  float fqmrnll[max_nfit];//Number of rings for this fit [1-6] fqmrpcflg[fqnmrfit]/I <0 if MINUIT did not converge during the fit
  float fqmrdconv[max_nfit][6];
  float fqmrmom[max_nfit][6];
  float fqmrpos[max_nfit][6][3];
  float fqmrdir[max_nfit][6][3];
  int fqmrpid[max_nfit][6];
    
  FQTree -> SetBranchAddress("nevt",&evtId);

  FQTree -> SetBranchAddress("fqnse",&fqnse);
  FQTree -> SetBranchAddress("fqtotq",fqtotq);
  FQTree -> SetBranchAddress("fqnhitpmt",fqnhitpmt);
  FQTree -> SetBranchAddress("fq1rpcflg",fq1rpcflg);
  FQTree -> SetBranchAddress("fq1rmom",fq1rmom);
  FQTree -> SetBranchAddress("fq1rpos",fq1rpos);
  FQTree -> SetBranchAddress("fq1rdir",fq1rdir);
  FQTree -> SetBranchAddress("fq1rnll",fq1rnll);
  FQTree -> SetBranchAddress("fq1rtotmu",fq1rtotmu);
  FQTree -> SetBranchAddress("fqpi0nll",fqpi0nll);
  FQTree -> SetBranchAddress("fqpi0mass",fqpi0mass);

  FQTree -> SetBranchAddress("fqpi0totmu",fqpi0totmu);
  FQTree -> SetBranchAddress("fqpi0mom1",fqpi0mom1);
  FQTree -> SetBranchAddress("fqpi0mom2",fqpi0mom2);
  FQTree -> SetBranchAddress("fqpi0momtot",fqpi0momtot);
  FQTree -> SetBranchAddress("fqpi0pos",fqpi0pos);

  FQTree -> SetBranchAddress("fqnmrfit",&fqnmrfit);
  FQTree -> SetBranchAddress("fqmrifit",fqmrifit);
  FQTree -> SetBranchAddress("fqmrnring",fqmrnring);
  FQTree -> SetBranchAddress("fqmrpcflg",fqmrpcflg);
  FQTree -> SetBranchAddress("fqmrnll",fqmrnll);
  FQTree -> SetBranchAddress("fqmrdconv",fqmrdconv);
  FQTree -> SetBranchAddress("fqmrmom",fqmrmom);
  FQTree -> SetBranchAddress("fqmrpos",fqmrpos);
  FQTree -> SetBranchAddress("fqmrdir",fqmrdir);
  FQTree -> SetBranchAddress("fqmrpid",fqmrpid);

  int distVtxNBins=50; double distVtxMin=0; double distVtxMax=500;
  int costhNBins=100; double costhMin=-1; double costhMax=1;
  int dwallNBins=70; double dwallMin=0; double dwallMax=3500;
  int towallNBins=200; double towallMin=0; double towallMax=10000;
  int MomentumNBins=500; double MomentumMin=0; double MomentumMax=1000;
  int LikelihoodNBins=750; double LikelihoodMin=-2500; double LikelihoodMax=5000;
  int Likelihood_eNBins=900; double Likelihood_eMin=0; double Likelihood_eMax=9e4;
  int Likelihood_pi0NBins=500; double Likelihood_pi0Min=0; double Likelihood_pi0Max=3e3;
  int TotalChargeNBins=1000; double TotalChargeMin=0; double TotalChargeMax=5e3;
  int Mass_pi0NBins=300; double Mass_pi0Min=0; double Mass_pi0Max=6e2;

  int NRingMin=0;int NRingMax=10;
  int Likelihood_1RMRNBins = 1500; double Likelihood_1RMRMin=-100;double Likelihood_1RMRMax=2900;
  int Likelihood_1RMR_simpleNBins = 300; double Likelihood_1RMR_simpleMin=-1500;double Likelihood_1RMR_simpleMax=1500;
  int PID_1RMRNBins = 64;double PID_1RMRMin = 0; double PID_1RMRMax = 64;
  int ConvL_1RMRNBins = 350;double ConvL_1RMRMin = -3500.; double ConvL_1RMRMax = 3500.;
  
  TH1D * EMom = new TH1D("EMom","Rec momentum assuming electron hypothesis",100,0,1000);
  TH1D * EPos = new TH1D("EPos","Rec position assuming electron hypothesis",100,0,1);
  TH1D * EDir = new TH1D("EDir","Rec direction assuming electron hypothesis",100,0,1);
  TH2D * NLL  = new TH2D("NLL","NLL for mu hypothesis X e hypothesis",500,0,5e4,500,0,5e4);

  //TH1D * Momentum = new TH1D("Momentum","",50,300,700);
  TH1D * Momentum = new TH1D("Momentum","",MomentumNBins,MomentumMin,MomentumMax);
  Momentum->GetXaxis()->SetTitle("Momentum (MeV/c)");
  Momentum->SetLineWidth(2);

  TH2D * MomentumXdWall = new TH2D("MomentumXdWall","",dwallNBins,dwallMin,dwallMax,MomentumNBins,MomentumMin,MomentumMax);
  MomentumXdWall->GetXaxis()->SetTitle("dWall (cm)");
  MomentumXdWall->GetYaxis()->SetTitle("Momentum (MeV/c)");

  TH2D * MomentumXtoWall = new TH2D("MomentumXtoWall","",towallNBins,towallMin,towallMax,MomentumNBins,MomentumMin,MomentumMax);
  MomentumXtoWall->GetXaxis()->SetTitle("toWall (cm)");
  MomentumXtoWall->GetYaxis()->SetTitle("Momentum (MeV/c)");

  TH1D * Likelihood = new TH1D("Likelihood","",LikelihoodNBins,LikelihoodMin,LikelihoodMax);
  Likelihood->GetXaxis()->SetTitle("Likelihood");
  Likelihood->SetLineWidth(2);

  TH1D * Likelihood_e = new TH1D("Likelihood_e","",Likelihood_eNBins,Likelihood_eMin,Likelihood_eMax);
  Likelihood_e->GetXaxis()->SetTitle("Likelihood e");
  Likelihood_e->SetLineWidth(2);

  TH2D * Likelihood_pi0 = new TH2D("Likelihood_pi0","",Mass_pi0NBins,Mass_pi0Min,Mass_pi0Max,Likelihood_pi0NBins,Likelihood_pi0Min,Likelihood_pi0Max);
  Likelihood_pi0->GetXaxis()->SetTitle("Likelihood #pi^{0}");
  Likelihood_pi0->SetLineWidth(2);
  
  TH2D * LikelihoodXdWall = new TH2D("LikelihoodXdWall","",dwallNBins,dwallMin,dwallMax,LikelihoodNBins,LikelihoodMin,LikelihoodMax);
  LikelihoodXdWall->GetXaxis()->SetTitle("LikelihoodXdWall");
  LikelihoodXdWall->SetLineWidth(2);

  TH2D * LikelihoodXtoWall = new TH2D("LikelihoodXtoWall","",towallNBins,towallMin,towallMax,LikelihoodNBins,LikelihoodMin,LikelihoodMax);
  LikelihoodXtoWall->GetXaxis()->SetTitle("LikelihoodXtoWall");
  LikelihoodXtoWall->SetLineWidth(2);
  
  TH2D * LikelihoodXMomentum = new TH2D("LikelihoodXMomentum","",MomentumNBins,MomentumMin,MomentumMax,LikelihoodNBins,LikelihoodMin,LikelihoodMax);
  LikelihoodXMomentum->GetXaxis()->SetTitle("LikelihoodXMomentum");
  LikelihoodXMomentum->SetLineWidth(2);

  TH1D * TotalCharge = new TH1D("TotalCharge","",TotalChargeNBins,TotalChargeMin,TotalChargeMax);
  TotalCharge->GetXaxis()->SetTitle("TotalCharge");
  TotalCharge->SetLineWidth(2);

  TH2D * TotalChargeXdWall = new TH2D("TotalChargeXdWall","",dwallNBins,dwallMin,dwallMax,TotalChargeNBins,TotalChargeMin,TotalChargeMax);
  TotalChargeXdWall->GetXaxis()->SetTitle("TotalChargeXdWall");
  TotalChargeXdWall->SetLineWidth(2);

  TH2D * TotalChargeXtoWall = new TH2D("TotalChargeXtoWall","",towallNBins,towallMin,towallMax,TotalChargeNBins,TotalChargeMin,TotalChargeMax);
  TotalChargeXtoWall->GetXaxis()->SetTitle("TotalChargeXtoWall");
  TotalChargeXtoWall->SetLineWidth(2);
  

  TH3D * Likelihood_pi0XdWall = new TH3D("Likelihood_pi0XdWall","",dwallNBins,dwallMin,dwallMax,Mass_pi0NBins,Mass_pi0Min,Mass_pi0Max,Likelihood_pi0NBins,Likelihood_pi0Min,Likelihood_pi0Max);
  TH3D * Likelihood_pi0XtoWall = new TH3D("Likelihood_pi0XtoWall","",towallNBins,towallMin,towallMax,Mass_pi0NBins,Mass_pi0Min,Mass_pi0Max,Likelihood_pi0NBins,Likelihood_pi0Min,Likelihood_pi0Max);

  TH1D * NRing = new TH1D("NRing","",NRingMax,NRingMin,NRingMax);
  NRing->GetXaxis()->SetTitle("NRing");
  NRing->SetLineWidth(2);

  TH2D * NRingXdWall = new TH2D("NRingXdWall","",dwallNBins,dwallMin,dwallMax,NRingMax,NRingMin,NRingMax);
  NRingXdWall->GetXaxis()->SetTitle("NRingXdWall");
  NRingXdWall->SetLineWidth(2);

  TH2D * NRingXtoWall = new TH2D("NRingXtoWall","",towallNBins,towallMin,towallMax,NRingMax,NRingMin,NRingMax);
  NRingXtoWall->GetXaxis()->SetTitle("NRingXtoWall");
  NRingXtoWall->SetLineWidth(2);

  TH1D * Likelihood_1RMR = new TH1D("Likelihood_1RMR","",Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMR->GetXaxis()->SetTitle("Likelihood");
  Likelihood_1RMR->SetLineWidth(2);

  TH1D * Likelihood_1RMR_e = new TH1D("Likelihood_1RMR_e","",Likelihood_eNBins,Likelihood_eMin,Likelihood_eMax);
  Likelihood_1RMR_e->GetXaxis()->SetTitle("Likelihood");
  Likelihood_1RMR_e->SetLineWidth(2);

  TH1D * Likelihood_1RMR_ee = new TH1D("Likelihood_1RMR_ee","",Likelihood_eNBins,Likelihood_eMin,Likelihood_eMax);
  Likelihood_1RMR_ee->GetXaxis()->SetTitle("Likelihood");
  Likelihood_1RMR_ee->SetLineWidth(2);

  TH2D * NRingXLikelihood_1RMR = new TH2D("NRingXLikelihood_1RMR","",Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax,NRingMax,NRingMin,NRingMax);
  NRingXLikelihood_1RMR->GetXaxis()->SetTitle("Likelihood");
  NRingXLikelihood_1RMR->SetLineWidth(2);

  TH2D * DistVtxRecXTrue = new TH2D("DistVtxRecXTrue","",distVtxNBins,distVtxMin,distVtxMax,distVtxNBins,distVtxMin,distVtxMax);
  DistVtxRecXTrue->GetXaxis()->SetTitle("true dist conv (cm)");
  DistVtxRecXTrue->GetYaxis()->SetTitle("rec dist conv (cm)");
  DistVtxRecXTrue->SetLineWidth(2);

  TH2D * Likelihood_1RMRXdistVtx = new TH2D("Likelihood_1RMRXdistVtx","",distVtxNBins,distVtxMin,distVtxMax,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXdistVtx->GetXaxis()->SetTitle("true dist conv (cm)");
  Likelihood_1RMRXdistVtx->GetYaxis()->SetTitle("NLL MR");
  Likelihood_1RMRXdistVtx->SetLineWidth(2);

  TH2D * Likelihood_1RMRXcosth = new TH2D("Likelihood_1RMRXcosth","",costhNBins,costhMin,costhMax,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXcosth->GetXaxis()->SetTitle("cos #theta");
  Likelihood_1RMRXcosth->GetYaxis()->SetTitle("NLL MR");
  Likelihood_1RMRXcosth->SetLineWidth(2);

  TH2D * Likelihood_1RMRXdWall = new TH2D("Likelihood_1RMRXdWall","",dwallNBins,dwallMin,dwallMax,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXdWall->GetXaxis()->SetTitle("LikelihoodXdWall");
  Likelihood_1RMRXdWall->SetLineWidth(2);

  TH2D * Likelihood_1RMRXtoWall = new TH2D("Likelihood_1RMRXtoWall","",towallNBins,towallMin,towallMax,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXtoWall->GetXaxis()->SetTitle("LikelihoodXtoWall");
  Likelihood_1RMRXtoWall->SetLineWidth(2);

  TH1D * PID_1RMR = new TH1D("PID_1RMR","",PID_1RMRNBins,PID_1RMRMin,PID_1RMRMax);
  PID_1RMR->GetXaxis()->SetTitle("PID MR");

  TH2D * PID_1RMRXdWall = new TH2D("PID_1RMRXdWall","",dwallNBins,dwallMin,dwallMax,PID_1RMRNBins,PID_1RMRMin,PID_1RMRMax);
  PID_1RMRXdWall->GetXaxis()->SetTitle("PIDXdWall");
  PID_1RMRXdWall->SetLineWidth(2);

  TH2D * PID_1RMRXtoWall = new TH2D("PID_1RMRXtoWall","",towallNBins,towallMin,towallMax,PID_1RMRNBins,PID_1RMRMin,PID_1RMRMax);
  PID_1RMRXtoWall->GetXaxis()->SetTitle("PIDXtoWall");
  PID_1RMRXtoWall->SetLineWidth(2);

  TH1D * VertexResolutionPrimary_1RMR = new TH1D("VertexResolutionPrimary_1RMR","",50,0,250);
  VertexResolutionPrimary_1RMR->GetXaxis()->SetTitle("primary vtx_{rec} - vtx_{true} (cm)");
  VertexResolutionPrimary_1RMR->SetLineWidth(2);

  TH1D * DirectionResolutionPrimary_1RMR = new TH1D("DirectionResolutionPrimary_1RMR","",180,0,180);
  DirectionResolutionPrimary_1RMR->GetXaxis()->SetTitle("primary #theta_{rec} - #theta_{true} (cm)");
  DirectionResolutionPrimary_1RMR->SetLineWidth(2);
  
  TH1D * MomentumResolutionPrimary_1RMR = new TH1D("MomentumResolutionPrimary_1RMR","",100,-250,250);
  //TH1D * MomentumResolutionPrimary_1RMR = new TH1D("MomentumResolutionPrimary_1RMR","",100,0,1000);
  MomentumResolutionPrimary_1RMR->GetXaxis()->SetTitle("primary p_{rec} - p_{true} (cm)");
  MomentumResolutionPrimary_1RMR->SetLineWidth(2);

  //TH2D * MomentumPrimaryRecXTrue_1RMR = new TH2D("MomentumPrimaryRecXTrue_1RMR","",100,0,1000,100,0,1000);
  TH2D * MomentumPrimaryRecXTrue_1RMR = new TH2D("MomentumPrimaryRecXTrue_1RMR","",100,0,1000,100,-250,250);
  MomentumPrimaryRecXTrue_1RMR->GetXaxis()->SetTitle("primary p_{rec} - p_{true} (MeV/c)");
  MomentumPrimaryRecXTrue_1RMR->GetYaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");

  //TH1D * MomentumSum_1RMR = new TH1D("MomentumSum_1RMR","",100,-250,250);
  TH1D * MomentumSum_1RMR = new TH1D("MomentumSum_1RMR","",100,0,1000);
  MomentumSum_1RMR->GetXaxis()->SetTitle("primary p_{rec} - p_{true} (cm)");
  MomentumSum_1RMR->SetLineWidth(2);
  //MomentumSum_1RMR->Write();

  TH2D * MomentumSumXMomentumPrimary_1RMR = new TH2D("MomentumSumXMomentumPrimary_1RMR","",100,0,1000,100,0,1000);
  MomentumSumXMomentumPrimary_1RMR->GetXaxis()->SetTitle("primary p_{rec} - p_{true} (cm)");
  MomentumSumXMomentumPrimary_1RMR->GetYaxis()->SetTitle("primary p_{rec} - p_{true} (cm)");
  MomentumSumXMomentumPrimary_1RMR->SetLineWidth(2);
  //MomentumSum_1RMR->Write();

  TH2D * Likelihood_1RMRXMomentumSum = new TH2D("Likelihood_1RMRXMomentumSum","",100,0,1000,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXMomentumSum->GetXaxis()->SetTitle("Sum p_{rec} (MeV/c)");
  Likelihood_1RMRXMomentumSum->GetYaxis()->SetTitle("likelihood");
  Likelihood_1RMRXMomentumSum->SetLineWidth(2);
  //MomentumSum_1RMR->Write();
  
  TH2D * Likelihood_1RMRXMomentumPrimary = new TH2D("Likelihood_1RMRXMomentumPrimary","",100,0,1000,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXMomentumPrimary->GetXaxis()->SetTitle("primary p_{rec} (MeV/c)");
  Likelihood_1RMRXMomentumPrimary->GetYaxis()->SetTitle("likelihood");
  Likelihood_1RMRXMomentumPrimary->SetLineWidth(2);
  //MomentumPrimary_1RMR->Write();
  
  TH2D * Likelihood_1RMRXMomentumSecondary = new TH2D("Likelihood_1RMRXMomentumSecondary","",100,0,1000,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXMomentumSecondary->GetXaxis()->SetTitle("secondary p_{rec} (MeV/c)");
  Likelihood_1RMRXMomentumSecondary->GetYaxis()->SetTitle("likelihood");
  Likelihood_1RMRXMomentumSecondary->SetLineWidth(2);
  //MomentumPrimary_1RMR->Write();

  TH2D * Likelihood_1RMRXConvLengthPrimary = new TH2D("Likelihood_1RMRXConvLengthPrimary","",ConvL_1RMRNBins,ConvL_1RMRMin,ConvL_1RMRMax,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXConvLengthPrimary->GetXaxis()->SetTitle("primary conv L (cm)");
  Likelihood_1RMRXConvLengthPrimary->GetYaxis()->SetTitle("likelihood");
  Likelihood_1RMRXConvLengthPrimary->SetLineWidth(2);
  //MomentumPrimary_1RMR->Write();
  
  TH2D * Likelihood_1RMRXConvLengthSecondary = new TH2D("Likelihood_1RMRXConvLengthSecondary","",ConvL_1RMRNBins,ConvL_1RMRMin,ConvL_1RMRMax,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXConvLengthSecondary->GetXaxis()->SetTitle("secondary conv L (cm)");
  Likelihood_1RMRXConvLengthSecondary->GetYaxis()->SetTitle("likelihood");
  Likelihood_1RMRXConvLengthSecondary->SetLineWidth(2);
  //MomentumPrimary_1RMR->Write();
  
  TH2D * ConvLPrimaryRecXTrue_1RMR = new TH2D("ConvLPrimaryRecXTrue_1RMR","",ConvL_1RMRNBins,ConvL_1RMRMin,ConvL_1RMRMax,ConvL_1RMRNBins,ConvL_1RMRMin,ConvL_1RMRMax);
  ConvLPrimaryRecXTrue_1RMR->GetXaxis()->SetTitle("true primary conv L (cm)");
  ConvLPrimaryRecXTrue_1RMR->GetYaxis()->SetTitle("rec primary conv L (cm)");
  ConvLPrimaryRecXTrue_1RMR->SetLineWidth(2);
  //MomentumPrimary_1RMR->Write();

  TH2D * ConvLSecondaryRecXTrue_1RMR = new TH2D("ConvLSecondaryRecXTrue_1RMR","",ConvL_1RMRNBins,ConvL_1RMRMin,ConvL_1RMRMax,ConvL_1RMRNBins,ConvL_1RMRMin,ConvL_1RMRMax);
  ConvLSecondaryRecXTrue_1RMR->GetXaxis()->SetTitle("true secondary conv L (cm)");
  ConvLSecondaryRecXTrue_1RMR->GetYaxis()->SetTitle("rec secondary conv L (cm)");
  ConvLSecondaryRecXTrue_1RMR->SetLineWidth(2);
  //MomentumSecondary_1RMR->Write();

  TH1D * VertexResolutionSecondary_1RMR = new TH1D("VertexResolutionSecondary_1RMR","",50,0,250);
  VertexResolutionSecondary_1RMR->GetXaxis()->SetTitle("secondary vtx_{rec} - vtx_{true} (cm)");
  VertexResolutionSecondary_1RMR->SetLineWidth(2);

  TH1D * DirectionResolutionSecondary_1RMR = new TH1D("DirectionResolutionSecondary_1RMR","",180,0,180);
  DirectionResolutionSecondary_1RMR->GetXaxis()->SetTitle("secondary #theta_{rec} - #theta_{true} (#theta)");
  DirectionResolutionSecondary_1RMR->SetLineWidth(2);
  
  TH1D * MomentumResolutionSecondary_1RMR = new TH1D("MomentumResolutionSecondary_1RMR","",100,-250,250);
  //TH1D * MomentumResolutionSecondary_1RMR = new TH1D("MomentumResolutionSecondary_1RMR","",100,0,1000);
  MomentumResolutionSecondary_1RMR->GetXaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV)");
  MomentumResolutionSecondary_1RMR->SetLineWidth(2);

  //TH2D * MomentumSecondaryRecXTrue_1RMR = new TH2D("MomentumSecondaryRecXTrue_1RMR","",100,0,1000,100,0,1000);
  TH2D * MomentumSecondaryRecXTrue_1RMR = new TH2D("MomentumSecondaryRecXTrue_1RMR","",100,0,1000,100,-250,250);
  MomentumSecondaryRecXTrue_1RMR->GetXaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");
  MomentumSecondaryRecXTrue_1RMR->GetYaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");

  TH2D * VertexResolutionPrimaryXSecondary_1RMR = new TH2D("VertexResolutionPrimaryXSecondary_1RMR","",50,0,250,50,0,250);
  VertexResolutionPrimaryXSecondary_1RMR->GetXaxis()->SetTitle("primary vtx_{rec} - vtx_{true} (cm)");
  VertexResolutionPrimaryXSecondary_1RMR->GetYaxis()->SetTitle("secondary vtx_{rec} - vtx_{true} (cm)");

  TH2D * DirectionResolutionPrimaryXSecondary_1RMR = new TH2D("DirectionResolutionPrimaryXSecondary_1RMR","",180,0,180,180,0,180);
  DirectionResolutionPrimaryXSecondary_1RMR->GetXaxis()->SetTitle("primary #theta_{rec} - #theta_{true} (#circ)");
  DirectionResolutionPrimaryXSecondary_1RMR->GetYaxis()->SetTitle("secondary #theta_{rec} - #theta_{true} (#circ)");

  TH2D * MomentumResolutionPrimaryXSecondary_1RMR = new TH2D("MomentumResolutionPrimaryXSecondary_1RMR","",100,-250,250,100,-250,250);
  //TH2D * MomentumResolutionPrimaryXSecondary_1RMR = new TH2D("MomentumResolutionPrimaryXSecondary_1RMR","",100,0,1000,100,0,1000);
  MomentumResolutionPrimaryXSecondary_1RMR->GetXaxis()->SetTitle("primary p_{rec} - p_{true} (MeV/c)");
  MomentumResolutionPrimaryXSecondary_1RMR->GetYaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");

  TH2D * VertexXDirectionResolutionPrimary_1RMR = new TH2D("VertexXDirectionResolutionPrimary_1RMR","",50,0,250,180,0,180);
  VertexXDirectionResolutionPrimary_1RMR->GetXaxis()->SetTitle("primary vtx_{rec} - vtx_{true} (cm)");
  VertexXDirectionResolutionPrimary_1RMR->GetYaxis()->SetTitle("primary #theta_{rec} - #theta_{true} (#circ)");

  TH2D * VertexXMomentumResolutionPrimary_1RMR = new TH2D("VertexXMomentumResolutionPrimary_1RMR","",50,0,250,100,-250,250);
  //TH2D * VertexXMomentumResolutionPrimary_1RMR = new TH2D("VertexXMomentumResolutionPrimary_1RMR","",50,0,250,100,0,1000);
  VertexXMomentumResolutionPrimary_1RMR->GetXaxis()->SetTitle("primary vtx_{rec} - vtx_{true} (cm)");
  VertexXMomentumResolutionPrimary_1RMR->GetYaxis()->SetTitle("primary p_{rec} - p_{true} (MeV/c)");

  TH2D * DirectionXMomentumResolutionPrimary_1RMR = new TH2D("DirectionXMomentumResolutionPrimary_1RMR","",180,0,180,100,-250,250);
  //TH2D * DirectionXMomentumResolutionPrimary_1RMR = new TH2D("DirectionXMomentumResolutionPrimary_1RMR","",180,0,180,100,0,1000);
  DirectionXMomentumResolutionPrimary_1RMR->GetXaxis()->SetTitle("primary #theta_{rec} - #theta_{true} (cm)");
  DirectionXMomentumResolutionPrimary_1RMR->GetYaxis()->SetTitle("primary p_{rec} - p_{true} (MeV/c)");

  TH2D * VertexXDirectionResolutionSecondary_1RMR = new TH2D("VertexXDirectionResolutionSecondary_1RMR","",50,0,250,180,0,180);
  VertexXDirectionResolutionSecondary_1RMR->GetXaxis()->SetTitle("secondary vtx_{rec} - vtx_{true} (cm)");
  VertexXDirectionResolutionSecondary_1RMR->GetYaxis()->SetTitle("secondary #theta_{rec} - #theta_{true} (#circ)");

  TH2D * VertexXMomentumResolutionSecondary_1RMR = new TH2D("VertexXMomentumResolutionSecondary_1RMR","",50,0,250,100,-250,250);
  //TH2D * VertexXMomentumResolutionSecondary_1RMR = new TH2D("VertexXMomentumResolutionSecondary_1RMR","",50,0,250,100,0,1000);
  VertexXMomentumResolutionSecondary_1RMR->GetXaxis()->SetTitle("secondary vtx_{rec} - vtx_{true} (cm)");
  VertexXMomentumResolutionSecondary_1RMR->GetYaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");

  TH2D * DirectionXMomentumResolutionSecondary_1RMR = new TH2D("DirectionXMomentumResolutionSecondary_1RMR","",180,0,180,100,-250,250);
  //TH2D * DirectionXMomentumResolutionSecondary_1RMR = new TH2D("DirectionXMomentumResolutionSecondary_1RMR","",180,0,180,100,0,1000);
  DirectionXMomentumResolutionSecondary_1RMR->GetXaxis()->SetTitle("secondary #theta_{rec} - #theta_{true} (cm)");
  DirectionXMomentumResolutionSecondary_1RMR->GetYaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");

  TH2D * trueEnergyDifferenceXtruePhotonAngle = new TH2D("trueEnergyDifferenceXtruePhotonAngle","",180,0,180,100,0.,1.);
  trueEnergyDifferenceXtruePhotonAngle->GetXaxis()->SetTitle("LikelihoodXtruePhotonAngle");
  //trueEnergyDifferenceXphotonAngle->SetLineWidth(2);

  TH2D * Likelihood_1RMRXtruePhotonAngle = new TH2D("Likelihood_1RMRXtruePhotonAngle","",180,0,180,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXtruePhotonAngle->GetXaxis()->SetTitle("LikelihoodXtruePhotonAngle");
  //Likelihood_1RMRXphotonAngle->SetLineWidth(2);

  TH2D * Likelihood_1RMRXrecPhotonAngle = new TH2D("Likelihood_1RMRXrecPhotonAngle","",180,0,180,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  Likelihood_1RMRXrecPhotonAngle->GetXaxis()->SetTitle("LikelihoodXrecPhotonAngle");
  //Likelihood_1RMRXphotonAngle->SetLineWidth(2);

  TH2D * VertexResolutionPrimaryXLikelihood_1RMR = new TH2D("VertexResolutionPrimaryXLikelihood_1RMR","",50,0,250,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * VertexResolutionPrimaryXLikelihood_1RMR = new TH2D("VertexResolutionPrimaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  VertexResolutionPrimaryXLikelihood_1RMR->GetXaxis()->SetTitle("primary vtx_{rec} - vtx_{true} (cm)");
  VertexResolutionPrimaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  TH2D * DirectionResolutionPrimaryXLikelihood_1RMR = new TH2D("DirectionResolutionPrimaryXLikelihood_1RMR","",180,-0,290,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * DirectionResolutionPrimaryXLikelihood_1RMR = new TH2D("DirectionResolutionPrimaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  DirectionResolutionPrimaryXLikelihood_1RMR->GetXaxis()->SetTitle("primary #theta_{rec} - #theta_{true} (#circ)");
  DirectionResolutionPrimaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  TH2D * MomentumResolutionPrimaryXLikelihood_1RMR = new TH2D("MomentumResolutionPrimaryXLikelihood_1RMR","",100,-250,250,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * MomentumResolutionPrimaryXLikelihood_1RMR = new TH2D("MomentumResolutionPrimaryXLikelihood_1RMR","",100,0,1000,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * MomentumResolutionPrimaryXLikelihood_1RMR = new TH2D("MomentumResolutionPrimaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  MomentumResolutionPrimaryXLikelihood_1RMR->GetXaxis()->SetTitle("primary p_{rec} - p_{true} (MeV/c)");
  MomentumResolutionPrimaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  TH2D * VertexResolutionSecondaryXLikelihood_1RMR = new TH2D("VertexResolutionSecondaryXLikelihood_1RMR","",50,0,250,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * VertexResolutionSecondaryXLikelihood_1RMR = new TH2D("VertexResolutionSecondaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  VertexResolutionSecondaryXLikelihood_1RMR->GetXaxis()->SetTitle("secondary vtx_{rec} - vtx_{true} (cm)");
  VertexResolutionSecondaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  TH2D * DirectionResolutionSecondaryXLikelihood_1RMR = new TH2D("DirectionResolutionSecondaryXLikelihood_1RMR","",180,-0,290,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * DirectionResolutionSecondaryXLikelihood_1RMR = new TH2D("DirectionResolutionSecondaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  DirectionResolutionSecondaryXLikelihood_1RMR->GetXaxis()->SetTitle("secondary #theta_{rec} - #theta_{true} (#circ)");
  DirectionResolutionSecondaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  TH2D * MomentumResolutionSecondaryXLikelihood_1RMR = new TH2D("MomentumResolutionSecondaryXLikelihood_1RMR","",100,-250,250,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * MomentumResolutionSecondaryXLikelihood_1RMR = new TH2D("MomentumResolutionSecondaryXLikelihood_1RMR","",100,0,1000,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * MomentumResolutionSecondaryXLikelihood_1RMR = new TH2D("MomentumResolutionSecondaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  MomentumResolutionSecondaryXLikelihood_1RMR->GetXaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");
  MomentumResolutionSecondaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  //TH2D * MomentumResolutionSecondaryXLikelihood_1RMR = new TH2D("MomentumResolutionSecondaryXLikelihood_1RMR","",100,-250,250,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  TH2D * DirectionResolutionSecondaryVSPrimaryXLikelihood_1RMR = new TH2D("DirectionResolutionSecondaryVSPrimaryXLikelihood_1RMR","",180,0,180,Likelihood_1RMRNBins,Likelihood_1RMRMin,Likelihood_1RMRMax);
  //TH2D * MomentumResolutionSecondaryXLikelihood_1RMR = new TH2D("MomentumResolutionSecondaryXLikelihood_1RMR","",100,0,1000,100,0,1000);
  DirectionResolutionSecondaryVSPrimaryXLikelihood_1RMR->GetXaxis()->SetTitle("secondary p_{rec} - p_{true} (MeV/c)");
  DirectionResolutionSecondaryVSPrimaryXLikelihood_1RMR->GetYaxis()->SetTitle("Likelihood 2R vs 1R");

  TH2D * Likelihood_1RMR_eeXep_2d = new TH2D("Likelihood_1RMR_eeXep_2d","",Likelihood_1RMR_simpleNBins,Likelihood_1RMR_simpleMin,Likelihood_1RMR_simpleMax,Likelihood_1RMR_simpleNBins,Likelihood_1RMR_simpleMin,Likelihood_1RMR_simpleMax);
  Likelihood_1RMR_eeXep_2d->GetXaxis()->SetTitle("Likelihood");
  Likelihood_1RMR_eeXep_2d->SetLineWidth(2);

  TH1D * Likelihood_1RMR_eeXep = new TH1D("Likelihood_1RMR_eeXep","",Likelihood_1RMR_simpleNBins,Likelihood_1RMR_simpleMin,Likelihood_1RMR_simpleMax);
  Likelihood_1RMR_eeXep->GetXaxis()->SetTitle("Likelihood");
  Likelihood_1RMR_eeXep->SetLineWidth(2);
      
  TH2D * Likelihood_1RMR_eeXepXdWall = new TH2D("Likelihood_1RMR_eeXepXdWall","",dwallNBins,dwallMin,dwallMax,Likelihood_1RMR_simpleNBins,Likelihood_1RMR_simpleMin,Likelihood_1RMR_simpleMax);
  Likelihood_1RMR_eeXepXdWall->GetXaxis()->SetTitle("LikelihoodXdWall");
  Likelihood_1RMR_eeXepXdWall->SetLineWidth(2);
  
  TH2D * Likelihood_1RMR_eeXepXtoWall = new TH2D("Likelihood_1RMR_eeXepXtoWall","",towallNBins,towallMin,towallMax,Likelihood_1RMR_simpleNBins,Likelihood_1RMR_simpleMin,Likelihood_1RMR_simpleMax);
  Likelihood_1RMR_eeXepXtoWall->GetXaxis()->SetTitle("LikelihoodXtoWall");
  Likelihood_1RMR_eeXepXtoWall->SetLineWidth(2);
  
  TH1D * Vertex[3];
  for(int i=0;i<3;i++){
    Vertex[i] = new TH1D(Form("Vertex%d",i),"",1000,-500,500);
    Vertex[i]->GetXaxis()->SetTitle("r (cm)");
    Vertex[i]->SetLineWidth(2);
  }
  
  TH1D * VertexResolution = new TH1D("VertexResolution","",500,0,500);
  VertexResolution->GetXaxis()->SetTitle("r (cm)");
  VertexResolution->SetLineWidth(2);

  TH2D * VertexResolutionXdWall = new TH2D("VertexResolutionXdWall","",dwallNBins,dwallMin,dwallMax,500,0,500);
  VertexResolutionXdWall->GetXaxis()->SetTitle("dWall (cm)");
  VertexResolutionXdWall->GetYaxis()->SetTitle("r (cm)");

  TH2D * VertexResolutionXtoWall = new TH2D("VertexResolutionXtoWall","",towallNBins,towallMin,towallMax,500,0,500);
  VertexResolutionXtoWall->GetXaxis()->SetTitle("toWall (cm)");
  VertexResolutionXtoWall->GetYaxis()->SetTitle("r (cm)");

  TH1D * DirectionResolution = new TH1D("DirectionResolution","",80,0,20);
  DirectionResolution->GetXaxis()->SetTitle("Angle (#circ)");
  DirectionResolution->SetLineWidth(2);

  TH2D * DirectionResolutionXdWall = new TH2D("DirectionResolutionXdWall","",dwallNBins,dwallMin,dwallMax,80,0,20);
  DirectionResolutionXdWall->GetXaxis()->SetTitle("dWall (cm)");
  DirectionResolutionXdWall->GetYaxis()->SetTitle("Angle (#circ)");

  double IDMomentumEvent0 = 0;//Just to identify the first event of the chain.
  cout<<"Starting the event loop"<<endl;
  int evtFake=0;
  double vtxTrue[3];
  double dirTrue[3];
  for(int ievt=0;ievt<nevt;ievt++){
    if(ievt%100 == 0) cout<<ievt<<"/"<<nevt<<endl;
    if(verbose){
      cout<<"##############################################################################################################"<<endl;
      cout<<"##############################################################################################################"<<endl;
      cout<<"##############################################EVENT #"<<ievt<<" ######################################################"<<endl;
      cout<<"##############################################################################################################"<<endl;
      cout<<"##############################################################################################################"<<endl;

    }
    FQTree->GetEntry(ievt);
    cout<<"Nevt / fitqun file = "<<fqnse<<endl;
    int evtWCSim=evtId-1;
    //std::cout<<"Evt = "<<ievt<<", id evt = "<<evtId<<std::endl;
    //if(evtId!=ievt) continue;
    WCSimTree->GetEntry(evtWCSim);      
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    for (int j=0; j<3; j++) vtxTrue[j] = wcsimrootevent->GetVtx(j);
    if(verbose){
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
      std::cout<<"Mode = "<<wcsimrootevent->GetMode()<<", Jmu = "<<wcsimrootevent->GetJmu()<<", Jp = "<<wcsimrootevent->GetJp()<<", Npar = "<<wcsimrootevent->GetNpar()<<", ? = "<<wcsimrootevent->GetVecRecNumber()<<std::endl;
      std::cout<<"Mode = "<<wcsimrootevent2->GetMode()<<", Jmu = "<<wcsimrootevent2->GetJmu()<<", Jp = "<<wcsimrootevent2->GetJp()<<", Npar = "<<wcsimrootevent2->GetNpar()<<", ? = "<<wcsimrootevent2->GetVecRecNumber()<<std::endl;
    }
// Get the number of tracks
  int ntrack = wcsimrootevent->GetNtrack();
  if(verbose) printf("ntracks=%d\n",ntrack);

  //if(ntrack > 10) continue;//I do not understand, but very strange thing happen not using this condition :(. To investigate.
  double particleStart[3];
  double particleStop[3];
  double particleP;
  double particleE;
  double particleDir[3];
  bool isSelected=false;
  const int max_gammas = 20;
  double gammaStart[max_gammas][3]={{0.}};
  double gammaStop[max_gammas][3]={{0.}};
  double gammaDir[max_gammas][3]={{0.}};
  double gammaNorm[max_gammas]={0.};
  double gammaP[max_gammas]={0.};
  int gammaID_secondary = 0;
  int gammaID_primary = 0;
  int ngammas=0;
  double pi0_convl_dist = 0;
  double gamma_convl_primary=0., gamma_convl_secondary=0.;
  double costh_gamma=-2;
  
  // Loop through elements in the TClonesArray of WCSimTracks
  for (int i=0; i<ntrack; i++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(i);
      
      //WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
      WCSimRootTrack *wcsimroottrack = (WCSimRootTrack*) (element);

      if(wcsimroottrack->GetFlag() == 0){
	for (int j=0; j<3; j++) dirTrue[j] = wcsimroottrack->GetDir(j);
      }
      if(verbose){
	printf("Track #%d\n",i);
	printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	cout<<"Track vtx: "<<wcsimroottrack->GetStart(0)<<", "<<wcsimroottrack->GetStart(1)<<", "<<wcsimroottrack->GetStart(2)<<" and stop: "<<wcsimroottrack->GetStop(0)<<", "<<wcsimroottrack->GetStop(1)<<", "<<wcsimroottrack->GetStop(2)<<endl;
	cout<<"Track dir: "<<wcsimroottrack->GetDir(0)<<", "<<wcsimroottrack->GetDir(1)<<", "<<wcsimroottrack->GetDir(2)<<endl;
	printf("Track energy: %f\n", wcsimroottrack->GetE());
	printf("Track momentum: %f\n", wcsimroottrack->GetP());
	printf("Track mass: %f\n", wcsimroottrack->GetM());
      std::cout<<"Ipnu = "<<wcsimroottrack->GetIpnu()<<", Flag = "<<wcsimroottrack->GetFlag()<<std::endl;
      }
      //if(i==(ntrack-1) && ){//Mother particle
      if(wcsimroottrack->GetIpnu() == pdgID[ParticleType] && TMath::Abs(wcsimroottrack->GetStart(0)-vtxTrue[0])<1.){
	for(int j=0; j<3; j++){
	  particleStart[j] = wcsimroottrack->GetStart(j);
	  particleStop[j] = wcsimroottrack->GetStop(j);
	  particleDir[j] = wcsimroottrack->GetDir(j);
	}
	particleP = wcsimroottrack->GetP();
	particleE = wcsimroottrack->GetE();
	double rad=TMath::Sqrt(pow(particleStop[0],2)+pow(particleStop[1],2));
	double height=TMath::Abs(particleStop[2]);
	bool trueStopTank=false;
	if(rad<TankRadius && height<TankHalfHeight) trueStopTank=true;
	if(verbose){
	  cout<<"Selected particle, but is stopping in the tank truly?"<<trueStopTank<<endl;
	  cout<<"Track Length = "<<TMath::Sqrt(pow(particleStop[0]-particleStart[0],2)+pow(particleStop[1]-particleStart[1],2)+pow(particleStop[2]-particleStart[2],2))<<endl;
	}
	isSelected=true;
      }
      else{
	cout<<"Not selected because ID ="<<wcsimroottrack->GetIpnu()<<", ptype = "<<ParticleType<<", ID searched for = "<<pdgID[ParticleType]<<", test = "<< pdgID[0] <<", particle type = "<<ParticleType<<", dist = "<<TMath::Abs(wcsimroottrack->GetStart(0)-vtxTrue[0])<<endl;
      }
      if(wcsimroottrack->GetParenttype()==111 && wcsimroottrack->GetIpnu()==22){//If parent is a pi0 and the particle is one of its 2 daughter gamma
	for(int j=0; j<3; j++){
	  gammaStart[ngammas][j]=wcsimroottrack->GetStart(j);
	  gammaStop[ngammas][j]=wcsimroottrack->GetStop(j);
	  gammaDir[ngammas][j]=wcsimroottrack->GetDir(j);
	  gammaNorm[ngammas]+=pow(gammaDir[ngammas][j],2);
	}
	gammaP[ngammas]=wcsimroottrack->GetP();
	ngammas++;
      }
    }  // End of loop over tracks
  
  if(ParticleType == 4){
    if(gammaP[0]>gammaP[1]){
      gammaID_primary = 0;
      gammaID_secondary = 1;
    }
    else{
      gammaID_primary = 1;
      gammaID_secondary = 0;
    }
    if(ngammas==2){
      gamma_convl_primary = TMath::Sqrt( pow(gammaStop[0][0]-gammaStart[0][0],2) + pow(gammaStop[0][1]-gammaStart[0][1],2) + pow(gammaStop[0][2]-gammaStart[0][2],2) );
      gamma_convl_secondary = TMath::Sqrt( pow(gammaStop[1][0]-gammaStart[1][0],2) + pow(gammaStop[1][1]-gammaStart[1][1],2) + pow(gammaStop[1][2]-gammaStart[1][2],2) );
      pi0_convl_dist = TMath::Sqrt( pow(gammaStop[0][0]-gammaStop[1][0],2) + pow(gammaStop[0][1]-gammaStop[1][1],2) + pow(gammaStop[0][2]-gammaStop[1][2],2) );
      double ps = gammaDir[0][0]*gammaDir[1][0]+gammaDir[0][1]*gammaDir[1][1]+gammaDir[0][2]*gammaDir[1][2];
      costh_gamma = ps/(TMath::Sqrt(gammaNorm[0])*TMath::Sqrt(gammaNorm[1]));
    }
    cout<<"Ngammas = "<<ngammas<<", costh = "<<costh_gamma<<", norm 1st gamma = "<<TMath::Sqrt(gammaNorm[0])<<", norm 2nd gamma = "<<TMath::Sqrt(gammaNorm[1])<<endl;
  }
  //if(pi0_convl_dist>25.) continue;
  //if(ParticleType_true==3) ParticleType=1; 
  //FQTree->GetEntry(ievt+evtFake);
  //FQTree->GetEntry(ievt);
    if(verbose){
      cout<<"Number of subevents = "<<fqnse<<endl;
      for(int s=0;s<fqnse;s++){
	cout<<"total charge = "<<fqtotq[s]<<endl;
      }
    }
    //if(ievt == 0) IDMomentumEvent0=fq1rmom[0][ParticleType];
    //else{
    //while(fq1rmom[0][ParticleType] == IDMomentumEvent0){
    //evtFake++;
    //FQTree->GetEntry(ievt+evtFake);
    //}
    //}
    if(!isSelected){
      cout<<endl<<endl<<endl<<endl;
      cout<<"ISSUE WITH THE TRUE VERTEX"<<endl;
      for(int i=0;i<3;i++){
     	cout<<"particle start = "<<particleStart[i]<<", vtx start="<<wcsimrootevent->GetVtx(i)<<endl;
      }	      //if(particleDir[i]!=
      cout<<endl<<endl<<endl<<endl;
      //break;
    }
    if(verbose){
      cout << "Rec vtx = (" << fq1rpos[0][ParticleType][0] << ", " << fq1rpos[0][ParticleType][1] << ", " << fq1rpos[0][ParticleType][2] << ")" << endl; 
      cout << "Rec dir = (" << fq1rdir[0][ParticleType][0] << ", " << fq1rdir[0][ParticleType][1] << ", " << fq1rdir[0][ParticleType][2] << ")" << endl;
      cout << "Rec mom = " << fq1rmom[0][ParticleType] << endl;
    }
    //if(TMath::Sqrt( pow(wcsimrootevent->GetVtx(0),2)+pow(wcsimrootevent->GetVtx(1),2)+pow(wcsimrootevent->GetVtx(2),2)) > 200) continue;
    double dwall, dwallBarrel, dwallCap;
    double StartRadius,StartHeight;
    //StartRadius = TMath::Sqrt(wcsimrootevent->GetVtx(0)*wcsimrootevent->GetVtx(0) +wcsimrootevent->GetVtx(1)*wcsimrootevent->GetVtx(1));
    //StartHeight = TMath::Abs(wcsimrootevent->GetVtx(2));
    StartRadius = TMath::Sqrt(particleStart[0]*particleStart[0]+particleStart[1]*particleStart[1]);
    StartHeight = TMath::Abs(particleStart[2]);
    dwallBarrel = TankRadius - StartRadius;
    dwallCap = TankHalfHeight - StartHeight;
    dwall = TMath::Min(dwallBarrel,dwallCap);

    //double stopRadius = TMath::Sqrt(pow(particleStop[0],2)+pow(particleStop[1],2));
    //double stopHalfHeight = TMath::Abs(particleStop[2]);
    //if(stopRadius > TankRadius || stopHalfHeight > TankHalfHeight) stopOffTank = true;
    //if(fq1rpcflg[0][ParticleType]!=stopOffTank){
    //cout<<"Particle stopped in tank or not??"<<endl<<"Stopping radius = "<<stopRadius<<", height = "<<stopHalfHeight<<", dwall = "<<dwall<<endl;
    //}
    if(containedOnly){
      if(stopReconstructed){
	if(fq1rpcflg[0][ParticleType]) continue;
      }
      else{
	bool stopOffTank = false;
	if(ParticleType == 4){
	  stopOffTank = (outOfTank(gammaStop[0]) || outOfTank(gammaStop[1]));
	}
	else{
	  stopOffTank = outOfTank(particleStop);
	}
	if(stopOffTank /*fq1rpcflg[0][ParticleType]*/){
	  if(verbose) cout<<"Not FC"<<endl;
	  continue;
	}
      }
    }
    if(dwallCut && dwall<dwallCutValue){//To remove interactions in PMTs
      continue;
    }

    
  //dwallBarrel = TankRadius - StartRadius;
  //dwallCap = TankHalfHeight - StartHeight;
  //dwall = TMath::Min(dwallBarrel,dwallCap);
    double mom_1r, totq_1r, vtx_1r[3], dir_1r[3];
    int ptype = ParticleType==4?1:ParticleType;//for pi0, use the electron hypothesis.
    mom_1r = fq1rmom[0][ptype];
    //if( (ParticleType==4 && mom_1r<500) || mom_1r<400) continue;
    
    totq_1r = fq1rtotmu[0][ptype];
    for(int j=0;j<3;j++){
      vtx_1r[j] = fq1rpos[0][ptype][j];
      dir_1r[j] = fq1rdir[0][ptype][j];
    }

    double pid_emu=fq1rnll[0][2]-fq1rnll[0][1]-0.2*fq1rmom[0][1];//fq1rnll[0][1]
    double pid=fq1rnll[0][2]-fq1rnll[0][1]-0.2*fq1rmom[0][1];//fq1rnll[0][1]
    double pid_pi0=-fqpi0nll[0]+fq1rnll[0][1];//Ln(pi0/e)
    double dist_vtx_pi0=0;
    for(int i=0;i<3;i++){
      dist_vtx_pi0+=pow(fqpi0pos[0][i]-fqpi0pos[1][i],2);
      cout<<fqpi0pos[0][i]<<" vs "<<fqpi0pos[1][i]<<endl;
    }
    dist_vtx_pi0 = TMath::Sqrt(dist_vtx_pi0);
    cout<<"pid pi0 = "<<pid_pi0<<", distance between pi0 vtx = "<<dist_vtx_pi0<<endl;
    //if(fqpi0mass[0]>300){
    //cout<<fq1rnll[0][2]-fq1rnll[0][1]<<endl;
    double vRes=0;double vDir=0;
    double normTrue=0;double normRec=0;

    for(int i=0;i<3;i++){
      vRes += pow(vtx_1r[i]-particleStart[i],2);
      vDir += dir_1r[i]*dirTrue[i];
      normTrue += pow(dirTrue[i],2);
      normRec += pow(dir_1r[i],2);
    }
    vDir /= TMath::Sqrt(normTrue*normRec);

    if(TMath::Sqrt(vRes) > 100 && verbose){
      cout << "Reso = " << TMath::Sqrt(vRes) << endl;
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
      cout << "Rec vtx = (" << vtx_1r[0] << ", " << vtx_1r[1] << ", " << vtx_1r[2] << ")" << endl;
    }
    double towall = FindToWall(particleStart,dirTrue);
#ifdef ONLYGOODEVENTS
    if(TMath::Abs(mom_1r - particleP)>100. || vRes>50.) continue;  
#endif
    if(towall<100){
    //if(fq1rnll[0][2]-fq1rnll[0][1]>0){
    //if(fqpi0mass[0]>300){
      cout<<"Evt fq #"<<ievt+evtFake<<"towall = "<<towall<<", dWall = "<<dwall<<endl;
      cout<<"e fitter results: vertex resolution = "<<TMath::Sqrt(vRes)<<", rec mom e = "<<fq1rmom[0][1]<<", rec mom mu = "<<fq1rmom[0][2]<<", likelihood ratio = "<<fq1rnll[0][2]-fq1rnll[0][1]<<", pid = "<<pid<<", Ln mu = "<<fq1rnll[0][2]<<", flag e = "<<fq1rpcflg[0][1]<<", flag mu = "<<fq1rpcflg[0][2]<<", number of subevents = "<<fqnse<<", total charge = "<<fqtotq[0]<<", total number of hits = "<<fqnhitpmt[0]<<", true radius = "<<StartRadius<<", height = "<<StartHeight<<endl;
      cout<<"Pi0 fitter results: Mass = "<<fqpi0mass[0]<<", rec mom of 1st photon = "<<fqpi0mom1[0]<<", rec mom 2nd photon = "<<fqpi0mom2[0]<<", rec mom of pi0 = "<<fqpi0momtot[0]<<", charge = "<<fqpi0totmu[0]<<endl<<endl;
    }
    //if(towall<100){
      //cout<<"Evt fq #"<<ievt+evtFake<<"towall = "<<towall<<", dWall = "<<dwall<<", vertex resolution = "<<TMath::Sqrt(vRes)<<", rec mom = "<<fq1rmom[0][1]<<", likelihood ratio = "<<fq1rnll[0][2]-fq1rnll[0][1]<<", pid = "<<pid<<", flag = "<<fq1rpcflg[0][1]<<", number of subevents = "<<fqnse<<", total charge = "<<fqtotq[0]<<", total number of hits = "<<fqnhitpmt[0]<<", radius = "<<StartRadius<<", height = "<<StartHeight<<endl;
    //}
    Momentum->Fill(mom_1r);
    Likelihood->Fill(pid_emu);
    Likelihood_pi0->Fill(fqpi0mass[0],pid_pi0);
    Likelihood_e->Fill(fq1rnll[0][1]);

    VertexResolutionXtoWall->Fill(towall,TMath::Sqrt(vRes));
    VertexResolutionXdWall->Fill(dwall,TMath::Sqrt(vRes));
    VertexResolution->Fill(TMath::Sqrt(vRes));
    for(int i=0;i<3;i++) Vertex[i]->Fill(vtx_1r[i]-particleStart[i]);
    //for(int i=0;i<3;i++) Vertex[i]->Fill(vtx_1r[i]-wcsimrootevent->GetVtx(i));
    //VertexResolution->Fill(TMath::Sqrt(vtx_1r[0]*vtx_1r[0]+vtx_1r[1]*vtx_1r[1]+vtx_1r[2]*vtx_1r[2]));
    DirectionResolutionXdWall->Fill(dwall,TMath::ACos(vDir)*180/TMath::Pi());
    DirectionResolution->Fill(TMath::ACos(vDir)*180/TMath::Pi());//Since beam is now along x direction
    MomentumXdWall->Fill(dwall,mom_1r);
    MomentumXtoWall->Fill(towall,mom_1r);
    LikelihoodXMomentum->Fill(mom_1r,pid);
    LikelihoodXdWall->Fill(dwall,pid);
    LikelihoodXtoWall->Fill(towall,pid);
    Likelihood_pi0XdWall->Fill(dwall,fqpi0mass[0],pid_pi0);
    Likelihood_pi0XtoWall->Fill(towall,fqpi0mass[0],pid_pi0);
    TotalCharge->Fill(totq_1r);
    TotalChargeXdWall->Fill(dwall,totq_1r);
    TotalChargeXtoWall->Fill(towall,totq_1r);
      //LikelihoodXMomentum->Fill(mom_1r,fq1rnll[0][2]-fq1rnll[0][1]-0.2*fq1rmom[0][1]);
    //LikelihoodXdWall->Fill(dwall,fq1rnll[0][2]-fq1rnll[0][1]-0.2*fq1rmom[0][1]);
    //LikelihoodXtoWall->Fill(towall,fq1rnll[0][2]-fq1rnll[0][1]-0.2*fq1rmom[0][1]);
    if(ievt%100 == 0 && verbose){
      cout << "Vx true = " << wcsimrootevent->GetVtx(0) << "Vx = " << vtx_1r[0] << ", dWall = " << dwall << ", toWall = " << towall << endl;
      cout<<"Mom e = "<<fq1rmom[0][1]<<", mu = "<<fq1rmom[0][2]<<", nll e = "<<fq1rnll[0][1]<<", mu = "<<fq1rnll[0][2]<<endl;
      cout<<"Number of tracks = " << ntrack << endl;
    }
    double LnLE = -fq1rnll[0][1];
    double LnLMu = -fq1rnll[0][2];
    //cout << LnLE - LnLMu << " vs " << 0.2*fq1rmom[0][1] << endl;
    //cout<<test<<endl;
    //cout<<"Nll e = "<<fq1rnll[0][1]<<", mu = "<<fq1rnll[0][2]<<endl;

    //EMom->Fill(fq1rmom[0][1]);
    //EPos->Fill(fq1rpos[0][1]);
    //EDir->Fill(fq1rdir[0][1]);
    NLL->Fill(fq1rnll[0][2],fq1rnll[0][1]);

    //Multi-ring part:
    int nring=0;
    double nll_diff=0;
    double nll_min=1e9;
    double nll_min_2rings_mainRingElectron=1e9;
    double dist_vtx=0;
    int pid_nll_min = 1e9;
    int pid_1rmr = 0;//Lower than 4 rings
    //for(int iFit=0;iFit<TMath::Min(fqnmrfit,max_nfit);iFit++){
    if(MRTuneCut){
      bool check_2rings=true;//false;
      //For 1R hypothesis
      double e_pid=0, e_mom=0, e_pos[3] = {0.}, e_dir[3] = {0.};
      //For 2R hypothesis. First entry is the most energetic ring of the 2, second is the less energetic.
      double ee_pid=0, ee_mom[2]={0.}, ee_pos[2][3] = {{0.}}, ee_dir[2][3] = {{0.}};
      //
      double pi_pid=0;
      double epi_pid=0;
      double pipi_pid=0;
      double pos_diff_primary =0., pos_diff_secondary = 0., dir_diff_primary = 0., dir_diff_secondary = 0., mom_diff_primary = 0., mom_diff_secondary = 0.;
      double mom_primary = 0., mom_secondary = 0.;
      double mom_sum = 0;
      double pos_diff_secondaryVSprimary = 0;
      double convl_primary = 0., convl_secondary = 0.;
      cout<<"Number of ring fits : "<<fqnmrfit<<endl;
      for(int iFit=0;iFit<fqnmrfit;iFit++){
	//cout<<"Fit "<<iFit<<", nrings = "<<fqmrnring[iFit]<<", MINUIT converged = "<<fqmrpcflg[iFit]<<", check fake rings = "<<fqmrifit[iFit]<<", ll = "<<fqmrnll[iFit]<<", conv L 1st ring = "<<fqmrdconv[iFit][0]<<", 2nd = "<<fqmrdconv[iFit][1]<<endl;
	cout<<endl<<"Fit "<<iFit<<", nrings = "<<fqmrnring[iFit]<<", MINUIT converged = "<<fqmrpcflg[iFit]<<", check fake rings = "<<fqmrifit[iFit]<<", ll = "<<fqmrnll[iFit]<<endl;
	cout<<"1st ring: PID = "<<fqmrpid[iFit][0]<<", vtx pos = "<<fqmrpos[iFit][0][0]<<", "<<fqmrpos[iFit][0][1]<<", "<<fqmrpos[iFit][0][2]<<", "<<endl;
	cout<<"2nd ring: PID = "<<fqmrpid[iFit][1]<<", vtx pos = "<<fqmrpos[iFit][1][0]<<", "<<fqmrpos[iFit][1][1]<<", "<<fqmrpos[iFit][1][2]<<", "<<endl;
	if(fqmrpcflg[iFit]>=0){//Means that minuit converged
	  //We wish to tune here the 2 ring cut for e-like event.
	  //We will therefore compare:
	  //Here, we wish to tune the case where we reconstruct 2 e-like ring vs 1 e-like. So we can compare the result :
	  //a. 8 : 1***1
	  //b. 8 : 2***11
	  //if(fqmrifit[iFit]>=0){//To avoid having rings counted as fake rings by the first part of MR algorithm
	  //if(((int) (fqmrifit[iFit]/1e8))==0){//Case where there is no final sequential fitter which merges rings.
	  //if( (((int) (fqmrifit[iFit]/1e8))==2) && (((int) (fqmrifit[iFit]/1e9))==0) ){//Case where there is no final sequential fitter which merges rings.
	  if( ((((int) (fqmrifit[iFit]/1e8))==2) && (((int) (fqmrifit[iFit]/1e9))==0)) || ((((int) (fqmrifit[iFit]/1e8))==0) && (fqmrnring[iFit]==1)) ){//Case where there is no final sequential fitter which merges rings.

	      //Access the two last digits of the index:
	      //int lastdigits=TMath::Abs(fqmrifit[iFit])%1000000;
	      int lastdigits = 0;
	      
	      if(((int) (fqmrifit[iFit]/1e8))==0){
		lastdigits=TMath::Abs(fqmrifit[iFit])%1000000;
	      }
	      else{
		if(fqmrnring[iFit]==1) lastdigits = fqmrpid[iFit][0];
		else if(fqmrnring[iFit]==2) lastdigits = fqmrpid[iFit][1]*10 + fqmrpid[iFit][0];
		else lastdigits = 100;
	      }

	    std::cout<<"Considered to fill plots, last digits = "<<lastdigits<<std::endl;
	      if(lastdigits==1){
		e_pid=fqmrnll[iFit];
		e_mom=fqmrmom[iFit][0];
		for(int k=0;k<3;k++){
		  e_pos[k]=fqmrpos[iFit][0][k];
		  e_dir[k]=fqmrdir[iFit][0][k];
		}
	      }
	      else if(lastdigits==3) pi_pid=fqmrnll[iFit];
	      else if(lastdigits==11){
		ee_pid=fqmrnll[iFit];
		if(fqmrnll[iFit]<nll_min_2rings_mainRingElectron){
		  nll_min_2rings_mainRingElectron = fqmrnll[iFit];
		  ee_mom[0]=fqmrmom[iFit][0];
		  ee_mom[1]=fqmrmom[iFit][1];
		  convl_primary=fqmrdconv[iFit][0];
		  convl_secondary=fqmrdconv[iFit][1];
		  
		  for(int k=0;k<3;k++){
		    ee_pos[0][k]=fqmrpos[iFit][0][k];
		    ee_dir[0][k]=fqmrdir[iFit][0][k];
		    ee_pos[1][k]=fqmrpos[iFit][1][k];
		    ee_dir[1][k]=fqmrdir[iFit][1][k];
		  }
		}
	      }
	      else if(lastdigits==31){
		epi_pid=fqmrnll[iFit];
		if(fqmrnll[iFit]<nll_min_2rings_mainRingElectron){
		  nll_min_2rings_mainRingElectron = fqmrnll[iFit];
		  ee_mom[0]=fqmrmom[iFit][0];
		  ee_mom[1]=fqmrmom[iFit][1];
		  convl_primary=fqmrdconv[iFit][0];
		  convl_secondary=fqmrdconv[iFit][1];
		  for(int k=0;k<3;k++){
		    ee_pos[0][k]=fqmrpos[iFit][0][k];
		    ee_dir[0][k]=fqmrdir[iFit][0][k];
		    ee_pos[1][k]=fqmrpos[iFit][1][k];
		    ee_dir[1][k]=fqmrdir[iFit][1][k];
		  }
		}
	      }
	      //else if(lastdigits==13 || lastdigits==31) epi_pid=fqmrnll[iFit];
	      else if(lastdigits==33) pipi_pid=fqmrnll[iFit];

	      if(fqmrnll[iFit]<nll_min && lastdigits<100){
		pid_nll_min = lastdigits;
		nll_min = fqmrnll[iFit];
		std::cout<<"New minimum"<<std::endl;
	      }
	  }
	    //}
	  //else cout<<"
	}
      }
      //nll_diff=e_pid-epi_pid;
      if(epi_pid==0) nll_diff=e_pid-ee_pid;
      else  nll_diff=e_pid-std::min(ee_pid,epi_pid);
      //if(nll_diff>100) continue;
      if(nll_diff>100 && nll_diff<200) cout<<"ISSUE 666"<<endl;
      //nll_diff=fq1rnll[0][1]-ee_pid;
      //if(nll_diff>1e3) cout<<"This is an error 666"<<endl;
      pid_1rmr = DeterminePID(pid_nll_min);//Lower than 4 rings
      cout<<"1 ring PID ="<<fq1rnll[0][1]<<", e pid ="<<e_pid<<", pi pid ="<<pi_pid<<", epi pid ="<<epi_pid<<", ee_pid="<<ee_pid<<", pipi_pid="<<pipi_pid<<", PID = "<<pid_1rmr<<endl;
      Likelihood_1RMR_e->Fill(e_pid);
      Likelihood_1RMR_ee->Fill(ee_pid);
      double nll_eeXe = e_pid - ee_pid;
      double nll_epXe = e_pid - epi_pid;
      Likelihood_1RMR_eeXep_2d->Fill(nll_eeXe,nll_epXe);
      Likelihood_1RMR_eeXep->Fill(nll_eeXe-nll_epXe);
      Likelihood_1RMR_eeXepXtoWall->Fill(towall,nll_eeXe-nll_epXe);
      Likelihood_1RMR_eeXepXdWall->Fill(dwall,nll_eeXe-nll_epXe);
      //double angle_gamma = TMath::ACos(scalar_product(gammaDir[gammaID_primary],particleDir,true))*180./TMath::Pi();
      double angle_gamma = TMath::ACos(scalar_product(gammaDir[gammaID_primary],gammaDir[gammaID_secondary],true))*180./TMath::Pi();
      Likelihood_1RMRXtruePhotonAngle->Fill(angle_gamma,nll_diff);
      double primary_energy = gammaP[gammaID_primary]/particleE;
      trueEnergyDifferenceXtruePhotonAngle->Fill(angle_gamma,primary_energy);
      double angle_gamma_rec = TMath::ACos(scalar_product(ee_dir[1],ee_dir[0],true))*180./TMath::Pi();
      Likelihood_1RMRXrecPhotonAngle->Fill(angle_gamma_rec,nll_diff);
      if(check_2rings){
	if(ParticleType == 4){
	  pos_diff_primary = vector_difference(ee_pos[0],gammaStop[gammaID_primary]);
	  pos_diff_secondary = vector_difference(ee_pos[1],gammaStop[gammaID_secondary]);
	  dir_diff_primary = TMath::ACos( scalar_product(ee_dir[0],gammaDir[gammaID_primary],true) )*180./TMath::Pi();
	  dir_diff_secondary = TMath::ACos( scalar_product(ee_dir[1],gammaDir[gammaID_secondary],true) )*180./TMath::Pi();
	  mom_diff_primary = ee_mom[0] - gammaP[gammaID_primary];
	  mom_diff_secondary = ee_mom[1] - gammaP[gammaID_secondary];
	  mom_primary = ee_mom[0];// - gammaP[gammaID_primary];
	  mom_secondary = ee_mom[1];// - gammaP[gammaID_secondary];
	  pos_diff_secondaryVSprimary = vector_difference(ee_pos[0],gammaStop[gammaID_secondary]);
	}
	else{
	  pos_diff_primary = vector_difference(ee_pos[0],particleStart);
	  pos_diff_secondary = vector_difference(ee_pos[1],particleStart);
	  dir_diff_primary = TMath::ACos( scalar_product(ee_dir[0],particleDir,true) )*180./TMath::Pi();
	  dir_diff_secondary = TMath::ACos( scalar_product(ee_dir[1],particleDir,true) )*180./TMath::Pi();
	  mom_diff_primary = ee_mom[0] - particleP;
	  mom_diff_secondary = ee_mom[1] - particleP;
	  mom_primary = ee_mom[0];// - particleP;
	  mom_secondary = ee_mom[1];// - particleP;
	  pos_diff_secondaryVSprimary = vector_difference(ee_pos[0],particleDir);
	}
	mom_sum = TMath::Sqrt( pow(ee_mom[0],2) + pow(ee_mom[1],2) + 2*ee_mom[0]*ee_mom[1]*scalar_product(ee_dir[0],ee_dir[1],true) );
      }
      else{
	  pos_diff_primary = vector_difference(e_pos,particleStart);
	  pos_diff_secondary = 0;//vector_difference(e_pos[1],particleStart);
	  dir_diff_primary = TMath::ACos( scalar_product(e_dir,particleDir,true) )*180./TMath::Pi();
	  dir_diff_secondary = 0;//TMath::ACos( scalar_product(ee_dir[1],particleDir,true) )*180./TMath::Pi();
	  mom_diff_primary = e_mom - particleP;
	  mom_diff_secondary = 0;//ee_mom[1] - particleP;
	  mom_primary = e_mom;// - particleP;
	  mom_secondary = 0;//ee_mom[1] - particleP;
	  mom_sum = e_mom;//TMath::Sqrt( pow(ee_mom[0],2) + pow(ee_mom[1],2) + 2*ee_mom[0]*ee_mom[1]*scalar_product(ee_dir[0],ee_dir[1],true) );
      }

      //if(TMath::Abs(convl_secondary)<200.) continue;
      cout<<"Difference in primary position = "<<pos_diff_primary<<"cm, secondary = "<<pos_diff_secondary<<"cm, difference in primary dirition = "<<dir_diff_primary<<"deg, secondary = "<<dir_diff_secondary<<"deg, momentum diff primary = "<<mom_diff_primary<<", momentum diff secondary = "<<mom_diff_secondary<<", convl primary = "<<convl_primary<<", secondary = "<<convl_secondary<<endl;

      VertexResolutionPrimary_1RMR->Fill(pos_diff_primary);
      DirectionResolutionPrimary_1RMR->Fill(dir_diff_primary);
      MomentumResolutionPrimary_1RMR->Fill(mom_diff_primary);
      MomentumPrimaryRecXTrue_1RMR->Fill(gammaP[gammaID_primary],mom_diff_primary);
      MomentumSum_1RMR->Fill(mom_sum);
      MomentumSumXMomentumPrimary_1RMR->Fill(mom_primary,mom_sum);
      Likelihood_1RMRXMomentumSum->Fill(mom_sum,nll_diff);
      Likelihood_1RMRXMomentumPrimary->Fill(mom_primary,nll_diff);
      Likelihood_1RMRXMomentumSecondary->Fill(mom_secondary,nll_diff);
      Likelihood_1RMRXConvLengthPrimary->Fill(convl_primary,nll_diff);
      Likelihood_1RMRXConvLengthSecondary->Fill(convl_secondary,nll_diff);
      ConvLPrimaryRecXTrue_1RMR->Fill(gamma_convl_primary,convl_primary);
      ConvLSecondaryRecXTrue_1RMR->Fill(gamma_convl_secondary,convl_secondary);
      //MomentumSum_1RMR->Fill(mom_sum);

      VertexResolutionSecondary_1RMR->Fill(pos_diff_secondary);
      DirectionResolutionSecondary_1RMR->Fill(dir_diff_secondary);
      MomentumResolutionSecondary_1RMR->Fill(mom_diff_secondary);
      MomentumSecondaryRecXTrue_1RMR->Fill(gammaP[gammaID_secondary],mom_diff_secondary);

      VertexResolutionPrimaryXSecondary_1RMR->Fill(pos_diff_primary,pos_diff_secondary);
      DirectionResolutionPrimaryXSecondary_1RMR->Fill(dir_diff_primary,dir_diff_secondary);
      MomentumResolutionPrimaryXSecondary_1RMR->Fill(mom_diff_primary,mom_diff_secondary); 

      VertexXDirectionResolutionPrimary_1RMR->Fill(pos_diff_primary,dir_diff_primary);
      VertexXMomentumResolutionPrimary_1RMR->Fill(pos_diff_primary,mom_diff_primary);
      DirectionXMomentumResolutionPrimary_1RMR->Fill(dir_diff_primary,mom_diff_primary);

      VertexXDirectionResolutionSecondary_1RMR->Fill(pos_diff_secondary,dir_diff_secondary);
      VertexXMomentumResolutionSecondary_1RMR->Fill(pos_diff_secondary,mom_diff_secondary);
      DirectionXMomentumResolutionSecondary_1RMR->Fill(dir_diff_secondary,mom_diff_secondary);

      VertexResolutionPrimaryXLikelihood_1RMR->Fill(pos_diff_primary,nll_diff);
      DirectionResolutionPrimaryXLikelihood_1RMR->Fill(dir_diff_primary,nll_diff);
      MomentumResolutionPrimaryXLikelihood_1RMR->Fill(mom_diff_primary,nll_diff);

      VertexResolutionSecondaryXLikelihood_1RMR->Fill(pos_diff_secondary,nll_diff);
      DirectionResolutionSecondaryXLikelihood_1RMR->Fill(dir_diff_secondary,nll_diff);
      MomentumResolutionSecondaryXLikelihood_1RMR->Fill(mom_diff_secondary,nll_diff);
      
      DirectionResolutionSecondaryVSPrimaryXLikelihood_1RMR->Fill(pos_diff_secondaryVSprimary,nll_diff);
      
    }
    else{
      //for(int iFit=0;iFit<TMath::Min(fqnmrfit,1);iFit++){
      //if(fqmrpcflg[iFit]>=0){//Means that minuit converged
	  //if(fqmrnll[iFit]<nll_min){
	  //nll_min=fqmrnll[iFit];
	  //nring=fqmrnring[iFit];
	  //}
	  //Search for min NLL w/o Dequential Fitter
      for(int iFit=0;iFit<fqnmrfit;iFit++){
	cout<<endl<<"Fit "<<iFit<<", nrings = "<<fqmrnring[iFit]<<", MINUIT converged = "<<fqmrpcflg[iFit]<<", check fake rings = "<<fqmrifit[iFit]<<", ll = "<<fqmrnll[iFit]<<endl;
	if(fqmrpcflg[iFit]>=0){//Means that minuit converged
	  if(fqmrifit[iFit]>=0){//To avoid having rings counted as fake rings by the first part of MR algorithm
	    //if(((int) (fqmrifit[iFit]/1e8))==0){//Case where there is no final sequential fitter which merges rings.
	    if( (((int) (fqmrifit[iFit]/1e8))==3) && (((int) (fqmrifit[iFit]/1e9))==0) ){//Case where there is no final sequential fitter which merges rings.
	    //if( (((int) (fqmrifit[iFit]/1e9))==0) && ){//Case where there is sequential fitter which merges rings.
	      //Access the two last digits of the index:
	      //int lastdigits=fqmrifit[iFit]%1000000;
	      //if(lastdigits==1 || lastdigits==11){
		if(fqmrnll[iFit]<nll_min){
		  cout<<"Nes best hypothesis"<<endl;
		  cout<<"1st ring: PID = "<<fqmrpid[iFit][0]<<", vtx pos = "<<fqmrpos[iFit][0][0]<<", "<<fqmrpos[iFit][0][1]<<", "<<fqmrpos[iFit][0][2]<<", "<<endl;
		  cout<<"2nd ring: PID = "<<fqmrpid[iFit][1]<<", vtx pos = "<<fqmrpos[iFit][1][0]<<", "<<fqmrpos[iFit][1][1]<<", "<<fqmrpos[iFit][1][2]<<", "<<endl;
		  nll_min=fqmrnll[iFit];
		  nring=fqmrnring[iFit];
		  if(nring==2){
		    dist_vtx = TMath::Sqrt(pow(fqmrpos[iFit][0][0]-fqmrpos[iFit][1][0],2)+pow(fqmrpos[iFit][0][1]-fqmrpos[iFit][1][1],2)+pow(fqmrpos[iFit][0][2]-fqmrpos[iFit][1][2],2));
		    cout<<"dist vtx = "<<dist_vtx<<endl;
		  }
		}
		//}
	    }
	  }
	  if(ievt%1 == 0 && verbose) cout<<"Fit = "<<iFit<<", flag = "<<fqmrpcflg[iFit]<<", number of rings ="<<fqmrnring[iFit]<<", nll = "<<fqmrnll[iFit]<<", min NLL = "<<nll_min<<endl;
	}
      }
      nll_diff=fq1rnll[0][1]-nll_min;
    }
    std::cout<<"Selected hypothesis has min nll = "<<nll_min<<", and nring = "<<nring<<std::endl;    
    NRing->Fill(nring);
    NRingXdWall->Fill(dwall,nring);
    NRingXtoWall->Fill(towall,nring);
    cout<<nll_diff<<endl;
    Likelihood_1RMR->Fill(nll_diff);
    cout<<pi0_convl_dist<<", "<<dist_vtx<<", nll min = "<<nll_diff<<endl;
    Likelihood_1RMRXdistVtx->Fill(pi0_convl_dist,nll_diff);
    Likelihood_1RMRXcosth->Fill(costh_gamma,nll_diff);
    DistVtxRecXTrue->Fill(pi0_convl_dist,dist_vtx);
    NRingXLikelihood_1RMR->Fill(nll_diff,nring);
    Likelihood_1RMRXdWall->Fill(dwall,nll_diff);
    Likelihood_1RMRXtoWall->Fill(towall,nll_diff);

    PID_1RMR->Fill(pid_1rmr);
    PID_1RMRXdWall->Fill(dwall,pid_1rmr);
    PID_1RMRXtoWall->Fill(towall,pid_1rmr);

    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    wcsimrootsuperevent2->ReInitialize();

  }

  if(Input) fOutput = new TFile(outfilename,"recreate");
  //if(Input) fOutput = new TFile(Form("plots/plots_%s",filename),"recreate");
  else  fOutput = new TFile(Form("plots/plots_fiTQun_%s.root",Suffix),"recreate");
  //else  fOutput = new TFile(Form("/disk01/usr5/bquilain/plots_fiTQun_%s.root",Suffix),"recreate");
  Momentum->Write();
  MomentumXdWall->Write();
  MomentumXtoWall->Write();
  Likelihood->Write();
  Likelihood_e->Write();
  Likelihood_pi0->Write();
  LikelihoodXtoWall->Write();
  LikelihoodXdWall->Write();
  Likelihood_pi0XdWall->Write();
  Likelihood_pi0XtoWall->Write();
  LikelihoodXMomentum->Write();
  TotalCharge->Write();
  TotalChargeXdWall->Write();
  TotalChargeXtoWall->Write();
  
  for(int i=0;i<3;i++) Vertex[i]->Write();
  VertexResolutionXtoWall->Write();
  VertexResolutionXdWall->Write();
  VertexResolution->Write();
  DirectionResolution->Write();
  DirectionResolutionXdWall->Write();
  NLL->Write();
  NRing->Write();
  NRingXdWall->Write();
  NRingXtoWall->Write();

  Likelihood_1RMR->Write();
  Likelihood_1RMR_e->Write();
  Likelihood_1RMR_ee->Write();
  NRingXLikelihood_1RMR->Write();
  Likelihood_1RMRXdistVtx->Write();
  Likelihood_1RMRXcosth->Write();
  DistVtxRecXTrue->Write();
  Likelihood_1RMRXtoWall->Write();
  Likelihood_1RMRXdWall->Write();
  Likelihood_1RMRXrecPhotonAngle->Write();
  Likelihood_1RMRXtruePhotonAngle->Write();
  trueEnergyDifferenceXtruePhotonAngle->Write();
  Likelihood_1RMRXMomentumSum->Write();
  Likelihood_1RMRXMomentumPrimary->Write();
  Likelihood_1RMRXMomentumSecondary->Write();
  Likelihood_1RMRXConvLengthPrimary->Write();
  Likelihood_1RMRXConvLengthSecondary->Write();
  
  ConvLPrimaryRecXTrue_1RMR->Write();
  ConvLSecondaryRecXTrue_1RMR->Write();
  PID_1RMR->Write();//Fill(pid_1rmr);
  PID_1RMRXdWall->Write();//Fill(dwall,pid_1rmr);
  PID_1RMRXtoWall->Write();//Fill(towall,pid_1rmr);

  VertexResolutionPrimary_1RMR->Write();
  DirectionResolutionPrimary_1RMR->Write();
  MomentumResolutionPrimary_1RMR->Write();
  MomentumPrimaryRecXTrue_1RMR->Write();
  MomentumSum_1RMR->Write();
  MomentumSumXMomentumPrimary_1RMR->Write();
  
  VertexResolutionSecondary_1RMR->Write();
  DirectionResolutionSecondary_1RMR->Write();
  MomentumResolutionSecondary_1RMR->Write();
  MomentumSecondaryRecXTrue_1RMR->Write();
  
  VertexResolutionPrimaryXSecondary_1RMR->Write();
  DirectionResolutionPrimaryXSecondary_1RMR->Write();//Fill(dir_diff_primary,dir_diff_secondary);
  MomentumResolutionPrimaryXSecondary_1RMR->Write();//Fill(mom_diff_primary,mom_diff_secondary); 
  
  VertexXDirectionResolutionPrimary_1RMR->Write();//Fill(pos_diff_primary,dir_diff_primary);
  VertexXMomentumResolutionPrimary_1RMR->Write();//Fill(pos_diff_primary,mom_diff_primary);
  DirectionXMomentumResolutionPrimary_1RMR->Write();//Fill(dir_diff_primary,mom_diff_primary);
  
  VertexXDirectionResolutionSecondary_1RMR->Write();//Fill(pos_diff_secondary,dir_diff_secondary);
  VertexXMomentumResolutionSecondary_1RMR->Write();//Fill(pos_diff_secondary,mom_diff_secondary);
  DirectionXMomentumResolutionSecondary_1RMR->Write();//Fill(dir_diff_secondary,mom_diff_secondary);
  
  VertexResolutionPrimaryXLikelihood_1RMR->Write();
  DirectionResolutionPrimaryXLikelihood_1RMR->Write();//Fill(dir_diff_primary,dir_diff_secondary);
  MomentumResolutionPrimaryXLikelihood_1RMR->Write();//Fill(mom_diff_primary,mom_diff_secondary); 

  VertexResolutionSecondaryXLikelihood_1RMR->Write();
  DirectionResolutionSecondaryXLikelihood_1RMR->Write();//Fill(dir_diff_secondary,dir_diff_secondary);
  MomentumResolutionSecondaryXLikelihood_1RMR->Write();//Fill(mom_diff_secondary,mom_diff_secondary); 
  DirectionResolutionSecondaryVSPrimaryXLikelihood_1RMR->Write();//Fill(mom_diff_secondary,mom_diff_secondary); 
  
  Likelihood_1RMR_eeXep_2d->Write();//(nll_eeXe,nll_epXe);
  Likelihood_1RMR_eeXep->Write();//(nll_eeXe-nll_epXe);
  Likelihood_1RMR_eeXepXtoWall->Write();//(towall,nll_eeXe-nll_epXe);
  Likelihood_1RMR_eeXepXdWall->Write();//(dwall,nll_eeXe-nll_epXe);
      
  fOutput->Close();
  //gApplication->Write();Terminate();
  return 0;
}
