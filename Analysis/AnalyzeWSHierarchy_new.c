#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
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
#include <TMath.h>
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"
//#define REMOVEMOST
//#define NEUTRONCAPTURE
//#define LEAFPDFONLY
bool separatedTriggers=false;//Assume two independent triggers, one for mPMT, one for B&L
bool produceDirectionPDF=true;//false;//true;
bool thetaOnly=false;//true;
WCSimRootGeom *geo = 0; 

bool eventDisplay_closeWall = false;
int stopCounterPlots=0;
int nPlotsStopping = 50;
int topPMT=19;
int nPMTpermPMT=19;
using namespace std;
// Simple example of reading a generated Root file
const int nPMTtypes = 1;
double PMTradius[nPMTtypes];
const int nGroupsPMTs=3;//Different PMT config in an mPMT. I assumed here a rotational symetry of the mPMT, so PMT 1 to 12 are the same, 13 to 18 are the same and 19 is separated

//2D plots
TH2D * ChargeProfile2D[nPMTtypes];
TH2D * ChargeProfile2D_onlyFront[nPMTtypes];
TH2D * ChargeProfile2DTop[nPMTtypes];
TH2D * ChargeProfile2DCylindricCoordinate[nPMTtypes];
TH2D * ChargeProfile2DRelativeAngle[nPMTtypes];
TH2D * TimeTOFProfileXTOF[nPMTtypes];
TH2D * ChargeProfile2D_event[nPMTtypes];

//1D plots of charge & hit
TH1D * ChargeProfile[nPMTtypes];
TH1D * HitProfile[nPMTtypes];

//1D time profiles
TH1D * TimeProfile[nPMTtypes];
TH1D * TimeHitProfile[nPMTtypes];
TH1D * TimeTOFProfile[nPMTtypes];
TH1D * TriggerTimeProfile[nPMTtypes];

//Charge per PMT
TH1D * ChargePerPMT[nPMTtypes];
TH1D * ChargePerPMT_decay[nPMTtypes];
TH1D * HitPerPMT[nPMTtypes];
TH1D * HitPerPMT_decay[nPMTtypes];

TH1D * TimeDifferenceHitPerPMT[nPMTtypes];
TH1D * TimeDifferenceHitPerPMT_decay[nPMTtypes];

//Total charge
TH1D * TotalCharge[nPMTtypes];
TH1D * TotalHit[nPMTtypes];
TH1D * TotalHit_decay[nPMTtypes];
TH1D * TotalHit_decay_nonshared[nPMTtypes];

TH2D * TotalChargeXdWall[nPMTtypes];

//new
TH1D * HitTimeTOFProfile[nPMTtypes];
TH3D * angularResponseXPMT_hits[nPMTtypes][nGroupsPMTs];
TH3D * angularResponseXPMT_allPMTs[nPMTtypes][nGroupsPMTs];
TH2D * distanceXweight_hits[nPMTtypes];
//Tests new
//TH3D * ChargeProfile2DXPMT[nPMTtypes];
TH2D * HitTimeTOFProfileXTOF[nPMTtypes];
//TH2D * ChargeProfile2D_onemPMT[nPMTtypes];

TH1D * neutrinoEnergy;

//2D plots
TH2D * ChargeProfileSpread[nPMTtypes];

double * crossProduct(double* v0, double* v1, double * v2){
  v2[0] = v0[1]*v1[2]-v0[2]*v1[1];
  v2[1] = v0[2]*v1[0]-v0[0]*v1[2];
  v2[2] = v0[0]*v1[1]-v0[1]*v1[0];
  return v2;
}

double scalarProduct(double* v0, double* v1){
  double total=0;
  for(int j=0;j<3;j++){
    total+=v0[j]*v1[j];
  }
  return total;
}

double groupPMTs(int pmtType, int pmt_number){
  if(pmtType==0) return 0;
  else{
    if(pmt_number<=12) return 0;
    else if(pmt_number<=18) return 1;
    else return nGroupsPMTs-1;
  }
}

double calculateWeight(double dist, double PMTrad, double partAng, int verbose=0){
  double L0=TMath::Sqrt(dist*dist+PMTrad*PMTrad-dist*2*PMTrad*TMath::Sin(partAng));
  double L1=TMath::Sqrt(dist*dist+PMTrad*PMTrad+dist*2*PMTrad*TMath::Sin(partAng));
  double num=pow(L0,2)+pow(L1,2)-pow(PMTrad*2,2);
  double denom=2*L0*L1;
  //double C=dist*dist+PMTrad*PMTrad;
  //double num=2*C-2*PMTrad*2*PMTrad;
  //double denom=TMath::Sqrt( (C+dist*2*PMTrad*TMath::Sin(partAng)) * (C-dist*2*PMTrad*TMath::Sin(partAng)) );
  double cosAlpha=num/denom;
  double cosHalfAlpha=TMath::Sqrt((1+cosAlpha)/2);
  //double cosAlpha=(2*C-2*PMTrad*2*PMTrad)/( C*TMath::Sqrt( 1 - pow((dist*PMTrad*2*TMath::Sin(partAng)/C),2) ) );
  double surfaceEff=dist*dist*2*TMath::Pi()*(1-cosHalfAlpha);//solid angle
  double fracSurface=surfaceEff/(4*TMath::Pi()*dist*dist);
  if(verbose==3){
    cout<<"Particle angle = "<<partAng*180/TMath::Pi()<<", Alpha = "<<TMath::ACos(cosAlpha)*180/TMath::Pi()<<", rad = "<<PMTrad<<", num = "<<num<<", denom = "<<denom<<endl;
    cout<<"Distance = "<<dist<<", area seen ="<<surfaceEff<<", fraction obs = "<<fracSurface<<endl;
  }
  return 1./fracSurface;
}

void vectorVertexPMT(double * vertex, int tubeNumber, int pmtType, double angles[4],int pmt_number_in_mpmt,int verbose=0){
  clock_t timeStart=clock();
  double ex[3]={0.};
  double ey[3]={0.};
  double ez[3]={0.};

  
  //1. Get the PMT position
  WCSimRootPMT pmt;
  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
  else pmt  = geo->GetPMT(tubeNumber-1,true);
  int mPMT_number = pmt.GetmPMTNo();
  if(verbose==3) cout<<"Original PMT number = "<<pmt_number_in_mpmt<<endl;
  double PMTpos[3];
  //double PMTradius=pmt.GetWCPMTRadius(pmtType);
  for(int j=0;j<3;j++){
    PMTpos[j] = pmt.GetPosition(j);
    ez[j] = pmt.GetOrientation(j);
  }
  if(verbose==3){
    cout<<"Position of the PMT = "<<PMTpos[0]<<", "<<PMTpos[1]<<", "<<PMTpos[2]<<endl<<"Orientation of the PMT = "<<ez[0]<<", "<<ez[1]<<", "<<ez[2]<<endl;
  }
  //2. Get the mPMT position and its direction -> Provides orgin and first axis of the referential (ez)
  //We do not have the mPMT position, so we will take the position of the PMT on top of the mPMT: number:
  double mPMTpos[3];
  double mPMTdir[3];

  int tubeNumber2 = tubeNumber + ( topPMT - pmt_number_in_mpmt );
  WCSimRootPMT pmt2;
  if(pmtType == 0) pmt2 = geo->GetPMT(tubeNumber2-1,false);
  else pmt2  = geo->GetPMT(tubeNumber2-1,true);
  int mPMT_number2 = pmt2.GetmPMTNo();
  int pmt_number_in_mpmt2 = pmt2.GetmPMT_PMTNo();
  for(int j=0;j<3;j++){
    mPMTpos[j]=pmt2.GetPosition(j);
    mPMTdir[j]=pmt2.GetOrientation(j);
  }
  if(verbose == 3){
    cout<<"Top PMT position = "<<mPMTpos[0]<<", "<<mPMTpos[1]<<", "<<mPMTpos[2]<<endl;
  }

  //3. Define the second vector of the referential in the orthogonal plane to ez and colinear to the center PMT -> Active PMT plane
  //a. PMT position in the mPMT referencec plane
  double PMTpos_mPMTref[3];
  for(int j=0;j<3;j++) PMTpos_mPMTref[j] = PMTpos[j] - mPMTpos[j];
  double totlen_PMTpos = TMath::Sqrt(pow(PMTpos_mPMTref[0],2)+pow(PMTpos_mPMTref[1],2)+pow(PMTpos_mPMTref[2],2));
  for(int j=0;j<3;j++) PMTpos_mPMTref[j]/=totlen_PMTpos;
  if(verbose==3) cout<<"Vertex position in the mPMT referential = "<<PMTpos_mPMTref[0]<<", "<<PMTpos_mPMTref[1]<<", "<<PMTpos_mPMTref[2]<<endl;
  //b. Define ey, perpendicular to ez and the center PMT -> Active PMT vector  
  crossProduct(ez,PMTpos_mPMTref,ey);
  double totlen_ey=TMath::Sqrt(pow(ey[0],2)+pow(ey[1],2)+pow(ey[2],2));
    for(int j=0;j<3;j++){
    ey[j]/=totlen_ey;
  }
  if(verbose==3) cout<<"Scal product test = "<<TMath::ACos(scalarProduct(ey,ez))*180/TMath::Pi()<<", "<<TMath::Sqrt(pow(ey[0],2)+pow(ey[1],2)+pow(ey[2],2))<<", "<<TMath::Sqrt(pow(ez[0],2)+pow(ez[1],2)+pow(ez[2],2))<<endl;

  //4. And then create ex which should be orthogonal to the 2 remaining vectors
  crossProduct(ey,ez,ex);
  //As ez and ey should be unitary, so does ex. So, let's check as a debug
  if(TMath::Sqrt(pow(ex[0],2)+pow(ex[1],2)+pow(ex[2],2))!=1){
    if(verbose==3) cout<<"There is a referential vector unitarity issue, length is ="<<TMath::Sqrt(pow(ex[0],2)+pow(ex[1],2)+pow(ex[2],2))<<endl;
  }

  //Conclusion: we now have our referential.
  
  //5. Now we have our referential, we should just calculate the angles of the PMT to vertex position vector in this referential.
  //a. calculate the PMT to vertex position vector.
  double particlepos_PMTref[3];
  for(int j=0;j<3;j++) particlepos_PMTref[j] = vertex[j] - PMTpos[j];
  double totlen_particlepos = TMath::Sqrt(pow(particlepos_PMTref[0],2)+pow(particlepos_PMTref[1],2)+pow(particlepos_PMTref[2],2));
  for(int j=0;j<3;j++) particlepos_PMTref[j]/=totlen_particlepos;
  if(verbose==3){
    cout<<"Vertex position = "<<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<endl;
    cout<<"Vertex position from PMT = "<<particlepos_PMTref[0]<<", "<<particlepos_PMTref[1]<<", "<<particlepos_PMTref[2]<<endl;
  }
  //b. Then extract Theta and Phi:
  double cosTheta=scalarProduct(particlepos_PMTref,ez);
  double Theta=TMath::ACos(cosTheta);
  //We know x=cosPhi x sinTheta and y=sinPhi x sinTheta
  double x=scalarProduct(particlepos_PMTref,ex);
  double y=scalarProduct(particlepos_PMTref,ey);
  double tanPhi=y/x;
  //tan is symetric from -pi/2 to +pi/2
  double Phi=TMath::ATan(tanPhi);
  if(x==0){
    if(y<0) Phi=-TMath::Pi()/2.;
    else Phi=TMath::Pi()/2.;
  }
  if(x<0) Phi+=TMath::Pi();//With this, angle become defines between -pi/2 to -pi/2 +2pi.
  //We wish to bring this from 0 to 2pi:
  if(Phi<0) Phi+=2*TMath::Pi();
  //Actually, we have a symmetry in Phi between [0,pi] and [pi,2pi]. So, we will just define Phi in [0,pi] modulo pi
  if(Phi > TMath::Pi()) Phi = TMath::Pi() - (Phi - TMath::Pi());
  //Phi=TMath::Abs(Phi);
  
  if(pmt_number_in_mpmt==topPMT) Phi=0.;//Phi is not defined in that case...
  angles[0]=Phi*180./TMath::Pi();
  angles[1]=Theta*180./TMath::Pi();
  angles[2]=totlen_particlepos;
  angles[3]=1;//calculateWeight(totlen_particlepos,PMTradius[pmtType],Theta,verbose);
  
  if(verbose==3) cout<<"Angles phi = "<<angles[0]<<", Theta = "<<angles[1]<<endl;
  clock_t timeEnd=clock();
  if(verbose == 3) cout<<"Time for the loop = "<<timeEnd-timeStart<<endl;
  
}


double findDirectionTheta(double * vertex,int tubeNumber,int pmtType,int pmt_number_in_mpmt,int verbose=0){
  clock_t timeStart=clock();
  double PMTpos[3];
  double PMTdir[3];
  WCSimRootPMT pmt;
  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
  else pmt  = geo->GetPMT(tubeNumber-1,true);
  for(int j=0;j<3;j++){
    PMTpos[j] = pmt.GetPosition(j);
    PMTdir[j] = pmt.GetOrientation(j);
  }

  int mPMT_number = pmt.GetmPMTNo();
  pmt_number_in_mpmt = pmt.GetmPMT_PMTNo();
  double particleRelativePMTpos[3];
  for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vertex[j];
  if(verbose == 3){
    cout<<"Original PMT number = "<<pmt_number_in_mpmt<<endl;  
    cout<<"Direction = "<<PMTdir[0]<<", "<<PMTdir[1]<<", "<<PMTdir[2]<<endl;
  }

  double norm=0;
  double ps=0;
  for(int j =0;j<3;j++){
    ps+=(-particleRelativePMTpos[j])*PMTdir[j];
    norm+=pow(particleRelativePMTpos[j],2);
  }
  double cosTheta=ps/TMath::Sqrt(norm);
  double Theta=TMath::ACos(cosTheta)*180/TMath::Pi();
  
  if(verbose == 3){
    cout<<"PMT #"<<pmt_number_in_mpmt<< ", Theta = "<< Theta << endl;
  }
  
  clock_t timeEnd=clock();
  if(verbose == 3) cout<<"Time for the loop = "<<timeEnd-timeStart<<endl;
  return Theta;
  //delete vPMTvertex;
  
}

int main(int argc, char **argv){

  //A test:
  TString test = "WCSim";  
  char * filename=NULL;
  char * outfilename=NULL;
  bool verbose = false;
  bool hybrid = true;//false;//true;
  bool Gamma = false;//Is the mother particle a gamma or another particle?
  float cvacuum = 3e8 / 1e9;//speed of light, in meter per ns.
  float nindex = 1.373;//refraction index of water
  bool plotDigitized = true;//false;

  bool HK=false;
  int startEvent=0;
  int endEvent=0;
  int c=-1;
  while( (c = getopt(argc,argv,"f:o:s:e:hv")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
      switch(c){
      case 'f':
	filename = optarg;
	break;
      case 'h':
	HK = true;
	break;
      case 'o':
	outfilename = optarg;
	break;
      case 's':
	startEvent = atoi(optarg);
	break;
      case 'e':
	endEvent = atoi(optarg);
	break;
      case 'v':
	verbose = true;
	break;
      }
  }
  
  // Clear global scope
  //gROOT->Reset();
  /*
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetStatColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleSize(0.04);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetPadBorderMode(0);
  gStylern->SetPalette(1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(.5);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetHatchesLineWidth(2);
  gStyle->SetLineWidth(1.5);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  */
  // Load the library with class dictionary info
  // (create with "gmake shared")
  //char* wcsimdirenv;
  //wcsimdirenv = getenv ("WCSIMDIR");
  //if(wcsimdirenv !=  NULL) gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  //else gSystem->Load("../libWCSimRoot.so");

  double TotalNHits = 1e4;
  double TankSize = 1100/2;//in cm, the maximal size of the plots for vertices etc...
  double TankRadius = 742/2;//in cm, the maximal size of the plots for vertices etc...
  double TankHalfHeight = 1042/2;//in cm, the maximal size of the plots for vertices etc...

  if(HK){
    TotalNHits = 3e4;
    TankSize = 7100;//in cm, the maximal size of the plots for vertices etc...
    TankRadius = 7080/2;//in cm, the maximal size of the plots for vertices etc...
    TankHalfHeight = 5480/2;//in cm, the maximal size of the plots for vertices etc...
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );
  }
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );

  
  TFile *file;

  // Open the file
  if (filename==NULL){
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hk_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hkmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_trash.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_2muniformsphere.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_2muniformsphere.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_1000events.root","read");   
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_10hitstrigger.root","read");
    //wcsim_hkhybridmpmt_e10_center_nominal_fulltank_0hitstrigger_1000events.root
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_0hitstrigger_1000events.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_0hitstrigger_nodn.root","read");


    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374_e10_center_nominal_fulltank_0hitstrigger_nodnno950.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374100Hz_e10_center_nominal_fulltank_0hitstrigger_10000.root","read");
    //TChain * 
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_4.2kHzAnd50Hz_14374100Hz_e10_center_nominal_fulltank_0hitstrigger_10000.root","read");

    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374100Hz_2.2MeVgamma_center_nominal_fulltank_0hitstrigger_100.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_2.2MeVgamma_center_nominal_fulltank_0hitstrigger_100.root","read");
    
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374100Hz_e3_center_nominal_fulltank_0hitstrigger.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374_e10_center_nominal_fulltank_0hitstrigger_noqe.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374_e10_trash.root","read");
    
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_2muniformsphere_100events.root","read");
    //file = new TFile("/disk01/usr5/bquilain/trash.root","read");

    
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_onlybal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_2.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hkmpmt_mu500_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_100_0_0_0/11_100_0_0_0.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_900_0_6_0/11_900_0_6_0.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkmpmt100HzDN_e4_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkmpmt_e4_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hk_e3_center_nominal_fulltank.root","read");
    //file = new TFile("../test.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_900_0_3_0/11_900_0_3_0.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_900_0_3_0/11_900_0_3_0.root","read");


    //tree = (TTree*)file->Get("wcsimT");
    //geotree = (TTree*)file->Get("wcsimGeoT");
    //opttree = (TTree*)file->Get("wcsimRootOptionsT");
    
  }else{
    file = new TFile(filename,"read");

  }
  if (!file->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }
  

  // Get the a pointer to the tree from the file
  TTree * tree = (TTree*)file->Get("wcsimT");
  TTree * geotree = (TTree*)file->Get("wcsimGeoT");
  TTree * opttree = (TTree*)file->Get("wcsimRootOptionsT");

  // Get the number of events
  int nevent = ((int)tree->GetEntries());//std::min(((int)tree->GetEntries()),100000);
  if(endEvent!=0 && endEvent<=nevent) nevent = endEvent;
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the tree in

  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  TBranch *branch2;
  if(hybrid){
    branch2 = tree->GetBranch("wcsimrootevent2");
    branch2->SetAddress(&wcsimrootsuperevent2);
  // Force deletion to prevent memory leak 
    tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
  }

  // Geometry tree - only need 1 "event"
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
      exit(9);
  }
  geotree->GetEntry(0);
  PMTradius[0]=geo->GetWCPMTRadius();
  PMTradius[1]=geo->GetWCPMTRadius(true);
  cout << "Number of PMTs of 1st type = " << geo->GetWCNumPMT() << ", radius = " << PMTradius[0] << endl;
  cout << "Number of PMTs of 2nd type = " << geo->GetWCNumPMT(true) << ", radius = " << PMTradius[1] <<endl;
  
  // Options tree - only need 1 "event"
  WCSimRootOptions *opt = 0; 
  opttree->SetBranchAddress("wcsimrootoptions", &opt);
  if(verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
  if (opttree->GetEntries() == 0) {
    exit(9);
  }
  opttree->GetEntry(0);
  opt->Print();

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;

  TH1F *h1 = new TH1F("PMT Hits", "PMT Hits", 8000, 0, 8000);
  TH1F *hvtx0 = new TH1F("Event VTX0", "Event VTX0", 200, -3500, 3500);
  TH1F *hvtx1 = new TH1F("Event VTX1", "Event VTX1", 200, -3500, 3500);
  TH1F *hvtx2 = new TH1F("Event VTX2", "Event VTX2", 200, -3500, 3500);

  // Alper's modification for position plots
  TH1F *hdir0 = new TH1F("Event DIR0", "Event DIR0", 200, -1.2, 1.2);
  TH1F *hdir1 = new TH1F("Event DIR1", "Event DIR1", 200, -1.2, 1.2);
  TH1F *hdir2 = new TH1F("Event DIR2", "Event DIR2", 200, -1.2, 1.2);
  //
  
  int num_trig=0;

    //TFile * fOutput = new TFile(OutputFile,"recreate");
#ifdef LEAFPDFONLY
  int nbins_position = 1;
  int nbins_angle = 720;
  int nbins_time = 1e4;
  int nbins_TOF = 1e3;
  int nbins_pe = 1;
  int nbins_totalcharge = 1;
  int nbins_distancePMT = 700;
  int nbins_phi = 18;
  int nbins_theta = 18;
  int nbins_weight = 1;
#else
  int nbins_position = TankSize/2;
  int nbins_angle = 720;
  int nbins_time = 1e4;
  int nbins_TOF = 1e3;
  int nbins_pe = 500;
  int nbins_totalcharge = 3e4;
  int nbins_distancePMT = 700;
  int nbins_phi = 18;
  int nbins_theta = 18;
  int nbins_weight = 200;  
#endif
  for(int i=0;i<nPMTtypes;i++){
    ChargeProfile2D_onlyFront[i] = new TH2D(Form("ChargeProfile2D_onlyFront_pmtType%d",i),"",nbins_position,-TankSize/2.,TankSize/2.,nbins_position,-TankSize/2.,TankSize/2.);
    ChargeProfile2D[i] = new TH2D(Form("ChargeProfile2D_pmtType%d",i),"",nbins_position,-TankSize/2.,TankSize/2.,nbins_position,-TankSize/2.,TankSize/2.);
    ChargeProfile2DTop[i] = new TH2D(Form("ChargeProfile2DTop_pmtType%d",i),"",nbins_position,-TankSize/2.,TankSize/2.,nbins_position,-TankSize/2.,TankSize/2.);
    ChargeProfile2DCylindricCoordinate[i] = new TH2D(Form("ChargeProfile2DCylindricCoordinate_pmtType%d",i),"",nbins_angle,0,360,nbins_position,-TankSize/2.,TankSize/2.);
    ChargeProfile2DRelativeAngle[i] = new TH2D(Form("ChargeProfile2DRelativeAngle_pmtType%d",i),"",nbins_angle,0,360,nbins_position,-TankSize/2.,TankSize/2.);
    ChargeProfile2D_event[i] = new TH2D(Form("ChargeProfile2D_event_pmtType%d",i),"",nbins_position/10,-TankSize/2.,TankSize/2.,nbins_position/10,-TankSize/2.,TankSize/2.);
    
    ChargeProfileSpread[i] = new TH2D(Form("ChargeProfileSpread_pmtType%d",i),"",1000,0,1000.,nbins_angle,0,180);
    ChargeProfile[i] = new TH1D(Form("ChargeProfile_pmtType%d",i),"",nbins_angle,0,180);
    HitProfile[i] = new TH1D(Form("HitProfile_pmtType%d",i),"",nbins_angle,0,180);

#ifdef NEUTRONCAPTURE
    TimeTOFProfileXTOF[i] = new TH2D(Form("TimeTOFProfileXTOF_pmtType%d",i),"",nbins_time,-1.5e3,5e5,nbins_TOF,0,5e2);TimeTOFProfileXTOF[i]->Sumw2();
    TimeProfile[i] = new TH1D(Form("TimeProfile_pmtType%d",i),"",nbins_time,0,5e5);
    TimeHitProfile[i] = new TH1D(Form("TimeHitProfile_pmtType%d",i),"",nbins_time,0,5e5);
    TimeTOFProfile[i] = new TH1D(Form("TimeTOFProfile_pmtType%d",i),"",nbins_time,-1.5e3,5e5);TimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfile[i] = new TH1D(Form("HitTimeTOFProfile_pmtType%d",i),"",nbins_time,-1.5e3,5e5);HitTimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfileXTOF[i] = new TH2D(Form("HitTimeTOFProfileXTOF_pmtType%d",i),"",nbins_time,-1.5e3,5e5,nbins_TOF,0,1e3);HitTimeTOFProfileXTOF[i]->Sumw2();
    TriggerTimeProfile[i] = new TH1D(Form("TriggerTimeProfile_pmtType%d",i),"",nbins_time,-100,5e5);
#else
    TimeTOFProfileXTOF[i] = new TH2D(Form("TimeTOFProfileXTOF_pmtType%d",i),"",nbins_time,-1.5e3,1.5e3,nbins_TOF,0,5e2);TimeTOFProfileXTOF[i]->Sumw2();
    TimeProfile[i] = new TH1D(Form("TimeProfile_pmtType%d",i),"",nbins_time,0,5e3);
    TimeHitProfile[i] = new TH1D(Form("TimeHitProfile_pmtType%d",i),"",nbins_time,0,5e3);
    TimeTOFProfile[i] = new TH1D(Form("TimeTOFProfile_pmtType%d",i),"",nbins_time,-1.5e3,1.5e3);TimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfile[i] = new TH1D(Form("HitTimeTOFProfile_pmtType%d",i),"",nbins_time,-1.5e3,1.5e3);HitTimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfileXTOF[i] = new TH2D(Form("HitTimeTOFProfileXTOF_pmtType%d",i),"",nbins_time,-1.5e3,1.5e3,nbins_TOF,0,1e3);HitTimeTOFProfileXTOF[i]->Sumw2();
    TriggerTimeProfile[i] = new TH1D(Form("TriggerTimeProfile_pmtType%d",i),"",nbins_time,-100,5e3);
#endif
    
    ChargePerPMT[i] = new TH1D(Form("ChargePerPMT_pmtType%d",i),"",nbins_pe,0,500);
    ChargePerPMT_decay[i] = new TH1D(Form("ChargePerPMT_decay_pmtType%d",i),"",nbins_pe,0,500);
    HitPerPMT[i] = new TH1D(Form("HitPerPMT_pmtType%d",i),"",10,0,10);
    HitPerPMT_decay[i] = new TH1D(Form("HitPerPMT_decay_pmtType%d",i),"",10,0,10);
    TimeDifferenceHitPerPMT[i] = new TH1D(Form("TimeDifferenceHitPerPMT_pmtType%d",i),"",nbins_time,-1.5e-3,1.5e3);
    TimeDifferenceHitPerPMT_decay[i] = new TH1D(Form("TimeDifferenceHitPerPMT_decay_pmtType%d",i),"",nbins_time,-1.5e-3,1.5e3);

    TotalCharge[i] = new TH1D(Form("TotalCharge_pmtType%d",i),"",nbins_totalcharge,0,3e4);
    TotalHit[i] = new TH1D(Form("TotalHit_pmtType%d",i),"",nbins_totalcharge,0,TotalNHits);
    TotalHit_decay[i] = new TH1D(Form("TotalHit_decay_pmtType%d",i),"",nbins_totalcharge,0,TotalNHits);
    TotalHit_decay_nonshared[i] = new TH1D(Form("TotalHit_decay_nonshared_pmtType%d",i),"",nbins_totalcharge,0,TotalNHits);
 
    TotalChargeXdWall[i] = new TH2D(Form("TotalChargeXdWall_pmtType%d",i),"",nbins_position,0,TankSize/2.,nbins_time,0,3e4);

    ChargeProfile2D[i]->Sumw2();
    ChargeProfile2D_onlyFront[i]->Sumw2();
    ChargeProfile2DTop[i]->Sumw2();
    ChargeProfile2DCylindricCoordinate[i]->Sumw2();
    ChargeProfile2DRelativeAngle[i]->Sumw2();
    ChargeProfile[i]->Sumw2(); HitProfile[i]->Sumw2();TimeProfile[i]->Sumw2();TimeHitProfile[i]->Sumw2();ChargePerPMT[i]->Sumw2();
    TotalCharge[i]->Sumw2();TotalHit[i]->Sumw2();
    TotalChargeXdWall[i]->Sumw2();
    ChargeProfile2D_event[i]->Sumw2();

    //new
    for(int p=0;p<nGroupsPMTs;p++){
      angularResponseXPMT_hits[i][p] = new TH3D(Form("angularResponseXPMT_hits_pmtType%d_pmtgroup%d",i,p),"",nbins_distancePMT,0,TankSize,nbins_phi,0,180,nbins_theta,0,180);
      angularResponseXPMT_allPMTs[i][p] = new TH3D(Form("angularResponseXPMT_allPMTs_pmtType%d_pmtgroup%d",i,p),"",nbins_distancePMT,0,TankSize,nbins_phi,0,180,nbins_theta,0,180);
    }
    distanceXweight_hits[i] = new TH2D(Form("distanceXweight_hits_pmtType%d",i),"",nbins_distancePMT,0,TankSize,nbins_weight,0,20);
    //ChargeProfile2DXPMT[i] = new TH3D(Form("ChargeProfile2DXPMT_pmtType%d",i),"",19,1,20,TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);

  }
  neutrinoEnergy = new TH1D("neutrinoEnergy","Neutrino Energy (MeV)",200,0.,60.);
  /*
  
  TH2D * HitProfile2D = new TH2D("HitProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * TimeProfile2D = new TH2D("TimeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * ChargeXPositionXTime = new TH2D("ChargeXPositionXTime","",36,0,180,TankSize,0,TankSize);
  TH2D * ChargeProfileXdWall = new TH2D("ChargeProfileXdWall","",nbins_angle,0,180,100,0,TankSize);
  TH2D * TotalChargeXdWall = new TH2D("TotalChargeXdWall","",100,0,TankSize,1e3,0,1e5);

  TH2D * TimeAngleProfile = new TH2D("TimeAngleHitProfile","",nbins_angle,0,180,1e4,0,5e3);TimeAngleProfile->Sumw2();
  TH2D * HitTimeAngleProfile = new TH2D("HitTimeAngleHitProfile","",nbins_angle,0,180,1e4,0,5e3);HitTimeAngleProfile->Sumw2();
  TH2D * TimeTOFAngleProfile = new TH2D("TimeTOFAngleProfile","",nbins_angle,0,180,1e4,-1e2,5e3);TimeTOFAngleProfile->Sumw2();
  TH2D * HitTimeTOFAngleProfile = new TH2D("HitTimeTOFAngleProfile","",nbins_angle,0,180,1e4,-1e2,5e3);HitTimeTOFAngleProfile->Sumw2();
  TH2D * TimeTOFTOF = new TH2D("TimeTOFTOF","",1e4,-1e2,5e3,1e2,0,5e2);TimeTOFTOF->Sumw2();
  TH2D * HitTimeTOFTOF = new TH2D("HitTimeTOFTOF","",1e4,-1e2,5e3,1e2,0,5e2);HitTimeTOFTOF->Sumw2();


  TH3D * VertexPosition = new TH3D("VertexPosition","",100,-TankSize,TankSize,100,-TankSize,TankSize,100,-TankSize,TankSize);
  TH3D * VertexDirection = new TH3D("VertexDirection","",100,-1.1,1.1,100,-1.1,1.1,100,-1.1,1.1);
  TH1D * VertexXdWall = new TH1D("VertexXdWall","",100,0,TankSize);
  TH1D * ParentFlyingDistance = new TH1D("ParentFlyingDistance","",10000,0,TankSize);

  TH2D * EGammaSeparationXdWall = new TH2D("EGammaSeparationxdWall","",100,0,TankSize,1000,0.,5.);
  TH1D * EGammaSeparation = new TH1D("EGammaSeparation","",1000,0.,5.);

  //Set the 2D histogram to 0 to increase lisibility of plots
  for(int ibinx=1;ibinx<=ChargeProfile2D->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=ChargeProfile2D->GetNbinsY();ibiny++){
      HitProfile2D->SetBinContent(ibinx,ibiny,0);
      ChargeProfile2D->SetBinContent(ibinx,ibiny,0);     
      ChargeProfile2DXZ->SetBinContent(ibinx,ibiny,0);     
    }
  }
  //
  double ConversionLength = 36;//cm
  TH2D * ChargeProfile2D_PerEvent[nevt];
  */
  
  //TH2D * ChargeProfile2D = new TH2D("ChargeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);

  // Now loop over events
  for (int ev=startEvent; ev<nevent; ev++)
  {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);

    //To store the decay-e ID.
    std::vector <int> decayID; decayID.clear();
    
    //New: Step 0 -> Create a map of the position of one PMT inside the mPMT
    /*
    //a. Determine the number of PMT per mPMT
    const int limit_maxPMT_per_mPMT = 100;
    int maxPMT_per_mPMT[nPMTtypes];
    
    //b. At the same time, for mPMT=0, fill the PMT position container.
    double position_PMT_local[nPMTtypes][limit_maxPMT_per_mPMT][3];
    double center_mPMT[nPMTtypes][3];//Position of the center of the mPMT
  
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      maxPMT_per_mPMT[pmtType]=0;
      for(int j=0;j<3;j++){
	center_mPMT[pmtType][j]=0;
      }
      //Ok this loop seems complicated but the idea is the following:
      //Using the fact that limit_max >> max_mPMT, we should loop over the PMTs of several mPMT.
      //For the first mPMT, we do not know at which PMT we will start. It may be at #3, and not #1.
      //To avoid that, we will only fill the table for PMT belonging to the second mPMT we will probe.
      //This is the purpose of this mPMT_counter, same_PMT etc.
      bool same_mPMT=true;
      int mPMT_counter=0;
      int previous_PMT_mPMT=-999;
      
      for(int p=0;p<limit_maxPMT_per_mPMT;p++){
	cout<<"debug, pmt #"<<p<<endl;
	WCSimRootPMT pmt_temp;
	if(pmtType == 0) pmt_temp = geo->GetPMT(p,false);
	else pmt_temp  = geo->GetPMT(p,true);
	int mPMT_number2 = pmt_temp.GetmPMTNo();
	int pmt_number_in_mPMT2 = pmt_temp.GetmPMT_PMTNo();
	cout << "Debug:, mPMT number = "<<mPMT_number2<<", number in the mPMT = "<<pmt_number_in_mPMT2<<endl;

	if(p==0) previous_PMT_mPMT=mPMT_number2;
	if(mPMT_number2 != previous_PMT_mPMT){
	  same_mPMT=false;
	  mPMT_counter++;
	  previous_PMT_mPMT=mPMT_number2;
	}

	if(mPMT_counter==0) continue;
	else if(mPMT_counter==1){
	  if(pmt_number_in_mPMT2>maxPMT_per_mPMT[pmtType]) maxPMT_per_mPMT[pmtType] = pmt_number_in_mPMT2;
	  for(int j=0;j<3;j++){
	    position_PMT_local[pmtType][pmt_number_in_mPMT2][j]=pmt_temp.GetPosition(j);
	    center_mPMT[pmtType][j]+=pmt_temp.GetPosition(j);
	  }
	}
      }
      for(int j=0;j<3;j++){
	if(maxPMT_per_mPMT[pmtType]!=0) center_mPMT[pmtType][j]/=maxPMT_per_mPMT[pmtType];
      }
    //c. What shall we keep in mind? That most of the position_PMT_local is non filled, but it is filled normally everywhere that we need.
    // Also, let's be careful. We do know that PMTs are numered from 1 to 19 in mPMT, not from 0 to 18. Let's keep this numbering to avoid any danger when someone may change the numbering in WCSim (e.g., moving to 0 -> 8

    //d. Now, we should find the central PMT of the mPMT. What can it be?
    // the best would be that the code find it by itself. 
    // There may be even some cases where there is no central PMT. So, we should find the module center. We can only average over the different component I guess.
      same_mPMT=true;
      mPMT_counter=0;
      previous_PMT_mPMT=-999;

      for(int p=0;p<limit_maxPMT_per_mPMT;p++){
	WCSimRootPMT pmt_temp;
	if(pmtType == 0) pmt_temp = geo->GetPMT(p,false);
	else pmt_temp  = geo->GetPMT(p,true);
	int mPMT_number2 = pmt_temp.GetmPMTNo();
	int pmt_number_in_mPMT2 = pmt_temp.GetmPMT_PMTNo();

	if(p==0) previous_PMT_mPMT=mPMT_number2;
	if(mPMT_number2 != previous_PMT_mPMT){
	  same_mPMT=false;
	  mPMT_counter++;
	  previous_PMT_mPMT=mPMT_number2;
	}
	if(mPMT_counter==0) continue;
	else if(mPMT_counter==1){
	  for(int j=0;j<3;j++){
	    position_PMT_local[pmtType][pmt_number_in_mPMT2][j]-=center_mPMT[pmtType][j];
	    cout<<j<<"="<<position_PMT_local[pmtType][pmt_number_in_mPMT2][j]<<", ";
	  }
	  cout<<endl;
	}
      }
      
    }
    */    
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    //wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());
      
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    }

    /*
    // Two ways to get the total number of tracks from the WCSimRootEvent
    cout << "Number of tracks = " << wcsimrootevent->GetTracks()->GetEntries() << endl;
    cout << "Number of tracks = " << wcsimrootevent->GetNtrack() << endl;

    // Just a print out of the direction of all tracks in the event
    cout << "***************************" << endl;
    for (int i=0 ; i<wcsimrootevent->GetNtrack() ; i++) {
      WCSimRootTrack *myTrack = (WCSimRootTrack*)(wcsimrootevent->GetTracks()->At(i));
      cout << "direction of Track #" << i << ": (" 
	   << myTrack->GetDir(0) << ", " 
	   << myTrack->GetDir(1) << ", " 
	   << myTrack->GetDir(2) << ")" << endl;
    }
    cout << "***************************" << endl;
    
    // Correct way to chose the real primary track
    cout << "***************************" << endl;
    for (int i=0 ; i<wcsimrootevent->GetNtrack() ; i++) {
      WCSimRootTrack *myTrack = (WCSimRootTrack*)(wcsimrootevent->GetTracks()->At(i));
      cout << myTrack->GetIpnu() << "   " << myTrack->GetParenttype() << " " 
	   << myTrack->GetP() << " " << myTrack->GetParentId() << " " << myTrack->GetFlag() << endl; 
      cout << "direction of Track #" << i << ": (" 
	   << myTrack->GetDir(0) << ", " 
	   << myTrack->GetDir(1) << ", " 
	   << myTrack->GetDir(2) << ")" << endl;
    }
    cout << "***************************" << endl;
    */
    
    // As one can see from previous lines of output, the correct way to identify the track
    // if to look for myTrack->GetFlag() = -1
    // for what I have seen, it is always track number 0
    for (int i=0 ; i<wcsimrootevent->GetNtrack() ; i++) {
      WCSimRootTrack *myTrack = (WCSimRootTrack*)(wcsimrootevent->GetTracks()->At(i));
      if ( myTrack->GetFlag() == -1 ) {
	if (verbose){
	  cout << "True direction : (" 
	     << myTrack->GetDir(0) << ", " 
	     << myTrack->GetDir(1) << ", " 
	     << myTrack->GetDir(2) << ")" << endl;
	}
      // Alper's modification
	hdir0->Fill(myTrack->GetDir(0));
	hdir1->Fill(myTrack->GetDir(1));
	hdir2->Fill(myTrack->GetDir(2));
      //
      }
    }


    hvtx0->Fill(wcsimrootevent->GetVtx(0));
    hvtx1->Fill(wcsimrootevent->GetVtx(1));
    hvtx2->Fill(wcsimrootevent->GetVtx(2));

    if(verbose){
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
      printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
      
    }

    std::vector<double> triggerInfo;
    triggerInfo.clear();
    triggerInfo = wcsimrootevent->GetTriggerInfo();

    std::vector<double> triggerInfo2;
    triggerInfo2.clear();
    if(hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();

#ifdef REMOVEMOST
    if(triggerInfo.size() == 0 && triggerInfo2.size() == 0) continue;
#endif

    if(verbose){
      for(int v=0;v<triggerInfo.size();v++){
	cout << "Trigger entry #" << v << ", info = " << triggerInfo[v] << endl;
      }
      if(hybrid){
	for(int v=0;v<triggerInfo2.size();v++){
	  cout << "Trigger2 entry #" << v << ", info = " << triggerInfo2[v] << endl;
	}
      }
    }

    double triggerShift[nPMTtypes];
    double triggerTime[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      triggerShift[pmtType]=0;
      triggerTime[pmtType]=0;
      if(triggerInfo.size()>=3){
	if(pmtType==0){
	  triggerShift[pmtType] = triggerInfo[1];
	  triggerTime[pmtType] = triggerInfo[2];
	}
      }
      if(triggerInfo2.size()>=3){
	if(pmtType==1 && hybrid){
	  triggerShift[pmtType] = triggerInfo2[1];
	  triggerTime[pmtType] = triggerInfo2[2];
	}
      }
      TriggerTimeProfile[pmtType]->Fill(triggerTime[pmtType]);
    }    
    // Now read the tracks in the event
    
    // Get the number of tracks
    int ntrack = wcsimrootevent->GetNtrack();
    if(verbose) printf("ntracks=%d\n",ntrack);

    double particleStart[3];
    double particleStop[3];
    double particleDir[3];
    //For decay-e
    double particleStart_decay[3];
    double particleStop_decay[3];
    double particleDir_decay[3];

    int i;
    // Loop through elements in the TClonesArray of WCSimTracks
    for (i=0; i<ntrack; i++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(i);
      
      //WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
      WCSimRootTrack *wcsimroottrack = (WCSimRootTrack*) (element);


      if(wcsimroottrack->GetParenttype() == 0 && wcsimroottrack->GetIpnu() != 0){//Mother particle
	//if(i==(ntrack-1)){//Mother particle
	for(int j=0; j<3; j++){
	  particleStart[j] = wcsimroottrack->GetStart(j);
	  particleStop[j] = wcsimroottrack->GetStop(j);
	  particleDir[j] = wcsimroottrack->GetDir(j);
	}
	//if(particleDir[0] != 1 || particleDir[1] != 0 || particleDir[2] != 0) cout << endl << endl << endl << "###############################################" << endl << "Particle direction has an issue" << endl << "###############################################" << endl;
	//if(particleStart[0] != 0 || particleStart[1] != 0 || particleStart[2] != 0) cout << endl << endl << endl << "###############################################" << endl << "Particle position has an issue" << endl << "###############################################" << endl;
      }

      //For decay-e
      if(wcsimroottrack->GetIpnu() == 11 && wcsimroottrack->GetParenttype() == 13){
	decayID.push_back(wcsimroottrack->GetId());
	for(int j=0; j<3; j++){
	  particleStart_decay[j] = wcsimroottrack->GetStart(j);
	  particleStop_decay[j] = wcsimroottrack->GetStop(j);
	}
	//decayID.push_back(wcsimroottrack->GetId());
      }

      if(verbose){
	printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	printf("Track ID: %d\n",wcsimroottrack->GetId());
	printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	printf("Track energy: %f\n", wcsimroottrack->GetE());
	printf("Track momentum: %f\n", wcsimroottrack->GetP());
	printf("Track mass: %f\n", wcsimroottrack->GetM());
      
	for (int j=0; j<3; j++){
	  printf("Track start: %d %f\n",j, wcsimroottrack->GetStart(j));
	  printf("Track dir: %d %f\n",j, wcsimroottrack->GetDir(j));
	}
      }

      if(wcsimroottrack->GetFlag()==-1){//This is the neutrino                                                                                                                            
	//ipnu=wcsimroottrack->GetIpnu();
	//for(int j=0; j<3; j++){
	//vnu[j] = wcsimroottrack->GetStop(j);
	//dirnu[j] = wcsimroottrack->GetDir(j);
	//}
	double pnu = wcsimroottrack->GetE();
	neutrinoEnergy->Fill(pnu);
	//cout<<"Energy = "<<pnu<<endl;
      }
    }  // End of loop over tracks
  
    if(verbose) cout<<"Particle start position = "<<particleStart[0]<<" "<<particleStart[1]<<" "<<particleStart[2]<<endl;

    double dwall, dwallBarrel, dwallCap;
    double StartRadius,StartHeight;
    double dwall_decay;
    if(Gamma){
      StartRadius = TMath::Sqrt(particleStop[0]*particleStop[0]+particleStop[1]*particleStop[1]);
      StartHeight = TMath::Abs(particleStop[2]);
    }
    else{
      StartRadius = TMath::Sqrt(particleStart[0]*particleStart[0]+particleStart[1]*particleStart[1]);
      StartHeight = TMath::Abs(particleStart[2]);
    }
 
   dwallBarrel = TankRadius - StartRadius;
   dwallCap = TankHalfHeight - StartHeight;
   dwall = min(dwallBarrel,dwallCap);

   //For decay-e
    double StartRadius_decay = TMath::Sqrt(particleStart_decay[0]*particleStart_decay[0]+particleStart_decay[1]*particleStart_decay[1]);
    double StartHeight_decay = TMath::Abs(particleStart_decay[2]);
    double dwallBarrel_decay = TankRadius - StartRadius_decay;
    double dwallCap_decay = TankHalfHeight - StartHeight_decay;
    dwall_decay = min(dwallBarrel_decay,dwallCap_decay);
    //

    //if(dwallCap>200) continue;
    if(ev%1 == 0){
      cout<<"Evt #"<<ev<<endl;
      //cout<<"**************************"<<endl<<"Vertex:"<<Start_x[0]<<", "<<Start_y[0]<<", "<<Start_z[0]<<endl;
      cout<<"Dwall = "<<dwall<<", radius = "<<StartRadius<<", Height = "<<StartHeight<<endl;
      //VertexXdWall->Fill(dwall);
    }
    // Now look at the Cherenkov hits
    
    // Get the number of Cherenkov hits.
    // Note... this is *NOT* the number of photons that hit tubes.
    // It is the number of tubes hit with Cherenkov photons.
    // The number of digitized tubes will be smaller because of the threshold.
    // Each hit "raw" tube has several photon hits.  The times are recorded.
    // See chapter 5 of ../doc/DetectorDocumentation.pdf
    // for more information on the structure of the root file.
    //  
    // The following code prints out the hit times for the first 10 tubes and also
    // adds up the total pe.
    // 
    // For digitized info (one time/charge tube after a trigger) use
    // the digitized information.
    //

    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits(); 
    int ncherenkovhits2 = 0; if(hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();
    int ncherenkovdigihits2 = 0;if(hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits(); 
    
    h1->Fill(ncherenkovdigihits);
    //if(verbose){
      printf("node id: %i\n", ev);
      printf("Ncherenkovhits %d\n",     ncherenkovhits);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      printf("Ncherenkovhits2 %d\n",     ncherenkovhits2);
      printf("Ncherenkovdigihits2 %d\n", ncherenkovdigihits2);
      cout << "RAW HITS:" << endl;
      //}

    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(separatedTriggers){
	if(triggerInfo2.size()!=0 && pmtType==0) continue;
	if(triggerInfo.size()!=0 && pmtType==1) continue;
      }
      if(verbose) cout << "PMT Type = " << pmtType << endl;
      int nPMTs = 0;
      if(produceDirectionPDF){
	nPMTs = geo->GetWCNumPMT(pmtType);
	for (int p=0; p< nPMTs ; p++)
	  {
	    //if(p%100 == 0) cout<<"PMT tested = "<<p<<endl;
	    int tubeNumber=p+1;
	    WCSimRootPMT pmt_temp;
	    if(pmtType == 0) pmt_temp = geo->GetPMT(p,false);
	    else pmt_temp  = geo->GetPMT(p,true);
	    int mPMT_number = pmt_temp.GetmPMTNo();
	    int pmt_number_in_mPMT = pmt_temp.GetmPMT_PMTNo();
	  //TMP
	    double dist_tmp =0;
	    for(int i=0;i<3;i++) dist_tmp+=pow(particleStart[i]-pmt_temp.GetPosition(i),2);
	    dist_tmp=TMath::Sqrt(dist_tmp);
	    //if(dist_tmp>500) continue;
	    //
	    // void vectorVertexPMT(double * vertex, int tubeNumber, int pmtType, double * vDirection,int pmt_number_in_mpmt,int verbose=0){
	    //double * vPMTvertex_temp = new double[3];vectorVertexPMT(particleStart,p+1,pmtType,vPMTvertex_temp,pmt_number_in_mPMT_temp);
	    double Phi=0;//=vPMTvertex_temp[0];
	    double Theta=0;//=vPMTvertex_temp[1];
	    double weight=1;
	    double distance_pmt_vertex=0;
	    if(thetaOnly){
	      Theta=findDirectionTheta(particleStart,p+1,pmtType,pmt_number_in_mPMT);
	    }
	    else{
	      if(pmtType==1){
		double vPMTvertex[4];
		vectorVertexPMT(particleStart,tubeNumber,pmtType,vPMTvertex,pmt_number_in_mPMT);
		Phi=vPMTvertex[0];
		Theta=vPMTvertex[1];
		distance_pmt_vertex=vPMTvertex[2];
		weight=vPMTvertex[3];
	      }
	    }
	      
	    //We just want phi:
	    //double cosTheta_temp = vPMTvertex_temp[2];
	    //double tanPhi_temp = vPMTvertex_temp[1]/vPMTvertex_temp[0];
	    //atan is defined between -pi/2 and pi/2
	    //double Theta_temp = TMath::ACos(cosTheta_temp);
	    //double Phi_temp;
	    //if(vPMTvertex_temp[0]>=0) Phi_temp=TMath::ATan(tanPhi_temp);
	    //else Phi_temp=-TMath::ATan(tanPhi_temp);
	    //cout<<"Phi angle = "<<Phi*180/TMath::Pi()<<", theta = "<<Theta*180/TMath::Pi()<<endl;
	    int pmtGroup=groupPMTs(pmtType,pmt_number_in_mPMT);
	    angularResponseXPMT_allPMTs[pmtType][pmtGroup]->Fill(distance_pmt_vertex,Phi,Theta,weight);
	  }
      }
 
      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray;//An array of pointers on CherenkovHitsTimes.
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();

      double particleRelativePMTpos[3];
      double totalPe = 0;
      int totalHit = 0;
      int totalHit_decay = 0;
      int totalHitPerPMT[nPMTs] = {0};
      int totalHitPerPMT_decay[nPMTs] = {0};
      vector <double> timePerPMT[nPMTs];//contains: time of each hit.
      vector <double> timePerPMT_decay[nPMTs];//contains: time of each hit only for decay-e hits.
      //vector <bool> isdecayPerPMT[nPMTs];//contains: proof it is a hit from a decay-e
      for(int ihit = 0; ihit<nPMTs;ihit++){
	timePerPMT[i].clear();
	timePerPMT_decay[i].clear();
      }

      int nhits;
      if(pmtType == 0) nhits = ncherenkovhits;
      else nhits = ncherenkovhits2;

      //cout<<"Number of hits non digitized = "<<nhits<<endl;
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (i=0; i< nhits ; i++)
	{
	  //if(verbose) cout << "Hit #" << i << endl;

	  WCSimRootCherenkovHit *wcsimrootcherenkovhit;
	  if(pmtType==0) wcsimrootcherenkovhit = (wcsimrootevent->GetCherenkovHits())->At(i);
	  else wcsimrootcherenkovhit = (wcsimrootevent2->GetCherenkovHits())->At(i);
	  
	  int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	  int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);

	  WCSimRootPMT pmt;
	  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
	  else pmt  = geo->GetPMT(tubeNumber-1,true);

	  double PMTpos[3];
	  double PMTdir[3];
	  for(int j=0;j<3;j++){
	    PMTpos[j] = pmt.GetPosition(j);
	    PMTdir[j] = pmt.GetOrientation(j);
	  }
	  
	  if ( i < 10 ) // Only print first XX=10 tubes
	    {
	      if(verbose){
		cout << "Number of entries in the array = " << timeArray->GetEntries() << endl;
		cout << "time array index = " << timeArrayIndex << endl;
		printf("Total pe: %d times( ",peForTube);
	      }
	      for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
		{
		  //cout << j << endl;
		  //WCSimRootCherenkovHitTime HitTime = (WCSimRootCherenkovHitTime) timeArray->At(j);		  
		  WCSimRootCherenkovHitTime* HitTime = (WCSimRootCherenkovHitTime*) timeArray->At(j);//Take the element j of an array of pointer on CherenkovHitTimes, so should be a pointer on a CherenkovHitsTimes		  
		  //  dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(j));
		  if(verbose) printf("%6.2f ", HitTime->GetTruetime() );
		}
	      if(verbose) cout << ")" << endl;
	      if(verbose) std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << endl;
	    }
	  
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
	  double phiAngle;
	  if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
	  else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));
	  
	  if(Gamma){//ProduceDirectionPDF the profile at the gamma conversion point, since it takes ~36 cm to convert/be visible
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStop[j];
	  }
	  else{
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStart[j];
	  }
		  
	  double vDir[3];double vOrientation[3];
	  double vDir2D[3];double vOrientation2D[3];
	  for(int j=0;j<3;j++){
	    vDir[j] = particleRelativePMTpos[j];
	    vOrientation[j] = PMTdir[j];
	    if(j<2){
	      vDir2D[j] = particleRelativePMTpos[j];
	      vOrientation2D[j] = PMTdir[j];
	    }	  
	  }
	  double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
	  double tof = Norm*1e-2/(cvacuum/nindex);
	  //if(verbose) cout << "Time of flight = " << tof << endl;
	  double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
	  for(int j=0;j<3;j++){
	    vDir[j] /= Norm;
	    vOrientation[j] /= NormOrientation;
	  }
	  double Norm2D = TMath::Sqrt(vDir2D[0]*vDir2D[0]+vDir2D[1]*vDir2D[1]);
	  double NormOrientation2D = TMath::Sqrt(vOrientation2D[0]*vOrientation2D[0]+vOrientation2D[1]*vOrientation2D[1]);
	  for(int j=0;j<2;j++){
	    vDir2D[j] /= Norm2D;
	    vOrientation2D[j] /= NormOrientation2D;
	  }
	  double HorizontalScalarProduct = TMath::ACos(vDir2D[0]*vOrientation2D[0] + vDir2D[1]*vOrientation2D[1])*180./TMath::Pi();
	  double relativeAngle = TMath::ACos(vDir[0]*particleDir[0]+vDir[1]*particleDir[1]+vDir[2]*particleDir[2])*180./TMath::Pi();
	  //WCSimRootCherenkovHitTime * HitTime = dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(timeArrayIndex));//Assumes that the earliest time the PMT is hit is in the first element -> Everything is ordered.

	  double Theta = 0;
	  double Phi = 0;
	  double weight = 1;
	  double distance_pmt_vertex=0;
	  int mPMT_number = pmt.GetmPMTNo();
	  int pmt_number_in_mPMT = pmt.GetmPMT_PMTNo();

	  //TMP
	    double dist_tmp =0;
	    for(int i=0;i<3;i++) dist_tmp+=pow(particleStart[i]-pmt.GetPosition(i),2);
	    dist_tmp=TMath::Sqrt(dist_tmp);
	    ////if(dist_tmp>1000) continue;
	    //if(dist_tmp>500) continue;
	    //

	  if(produceDirectionPDF){
	    //double * vPMTvertex = new double[3];vectorVertexPMT(geo,particleStart,pmt,tubeNumber,pmtType,vPMTvertex);
	    if(thetaOnly){
	      Theta=findDirectionTheta(particleStart,tubeNumber,pmtType,pmt_number_in_mPMT);
	    }
	    else{
	      if(pmtType==1){
		double vPMTvertex[4];
		vectorVertexPMT(particleStart,tubeNumber,pmtType,vPMTvertex,pmt_number_in_mPMT);
		//cout<<endl;
		Phi=vPMTvertex[0];
		Theta=vPMTvertex[1];
		distance_pmt_vertex=vPMTvertex[2];
		weight=vPMTvertex[3];
	      }
	    }

	    //We just want phi:
	    //double cosTheta = vPMTvertex[2];
	    //double tanPhi = vPMTvertex[1]/vPMTvertex[0];
	    //atan is defined between -pi/2 and pi/2
	    //Theta = TMath::ACos(cosTheta);
	    //if(vPMTvertex[0]>=0) Phi=TMath::ATan(tanPhi);
	    //else Phi=-TMath::ATan(tanPhi);
	    //cout<<"Phi angle = "<<Phi*180/TMath::Pi()<<", theta = "<<Theta*180/TMath::Pi()<<endl;	  
	  }
	  
	  WCSimRootCherenkovHitTime * HitTime = (WCSimRootCherenkovHitTime*) timeArray->At(0);//Takes the first hit of the array as the timing, It should be the earliest hit
	  //WCSimRootCherenkovHitTime HitTime = (WCSimRootCherenkovHitTime) timeArray->At(j);		  
	  //double time = HitTime->GetTruetime();
	
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  
	  //std::cout << "Relative angle = " << relativeAngle << ", time = " << time << endl;
	
	  if(!plotDigitized){
	    //plotBQ(PMTpos[3], PMTdir[3],phiAngle,relativeAngle,peForTube,pmtType);
	    //if(TMath::Abs(PMTpos[2]) < 4500){//Temp
	      if(phiAngle < TMath::Pi()/2 || phiAngle > 3*TMath::Pi()/2) ChargeProfile2D_onlyFront[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2D[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2DTop[pmtType]->Fill(PMTpos[0],PMTpos[1],peForTube);
	      ChargeProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],peForTube);
	      ChargeProfile2DRelativeAngle[pmtType]->Fill(relativeAngle,PMTpos[2],peForTube);
	      ChargeProfile2D_event[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	 
	      ChargeProfileSpread[pmtType]->Fill(ev,relativeAngle,peForTube);
	      ChargeProfile[pmtType]->Fill(relativeAngle,peForTube);
	      HitProfile[pmtType]->Fill(relativeAngle,1);

	      if(i%1000 == 0) cout<<"tube ID = "<<tubeNumber<<", time="<<time<<", pe = "<<peForTube<<endl;
    
	      for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
		{
		  WCSimRootCherenkovHitTime* HitTime = (WCSimRootCherenkovHitTime*) timeArray->At(j);//Take the element j of an array of pointer on CherenkovHitTimes, so should be a pointer on a CherenkovHitsTimes
		  double time = HitTime->GetTruetime();
		  int parentID = HitTime->GetParentID();
		  //For decay-e
		  //bool isdecay = false;
		  for(int id = 0; id<decayID.size(); id++){
		    if(parentID == decayID[id]){
		      totalHitPerPMT_decay[tubeNumber-1]++; 
		      timePerPMT_decay[tubeNumber-1].push_back(time);
		      //cout<<"This is a decay-e hit, with parent ID = "<<parentID<<", and PMT ID = "<<tubeNumber - 1<<", and charge = "<<peForTube<<", test hit = "<<totalHitPerPMT_decay[tubeNumber-1]<<endl;
		      //isdecay = true;
		      //isdecayPerPMT[tubeNumber-1].push_back(time);
		    }
		  }
		  //if(parentID > 50) continue;
		  if((i%100 == 0) && (j==timeArrayIndex)) cout<<"Time = "<<time<<", parent ID = "<<parentID<<endl;
		  TimeProfile[pmtType]->Fill(time,peForTube);
		  TimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],peForTube);
		  TimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,peForTube);
		  
		  HitTimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],1);
		  HitTimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,1);
		  timePerPMT[tubeNumber-1].push_back(time);
		  //isdecayPerPMT[tubeNumber-1].push_back(isdecay);
		  //cout<<timePerPMT[tubeNumber-1].size()<<endl;
		  totalHitPerPMT[tubeNumber-1]++;
		}
	      
	      int pmtGroup=groupPMTs(pmtType,pmt_number_in_mPMT);
	      angularResponseXPMT_hits[pmtType][pmtGroup]->Fill(distance_pmt_vertex,Phi,Theta,weight);
	      //angularResponseXPMT_hits[pmtType]->Fill(pmtGroup,Phi,Theta,weight);
	      distanceXweight_hits[pmtType]->Fill(distance_pmt_vertex,TMath::Log(weight));
	    //angularResponseXPMT_allPMTs[pmtType]->Fill(pmt_number_in_mPMT,Phi*180/TMath::Pi(),Theta*180/TMath::Pi(),1);
	    //ChargeProfile2DXPMT[pmtType]->Fill(pmt_number_in_mPMT,PMTpos[1],PMTpos[2],peForTube);
	      
	      ChargePerPMT[pmtType]->Fill(peForTube);
	      if(decayID.size()!=0 && dwall_decay > 0){
		ChargePerPMT_decay[pmtType]->Fill(totalHitPerPMT_decay[tubeNumber-1]);
		if(totalHitPerPMT_decay[tubeNumber-1] != 0){
		  totalHit_decay ++;
		  //cout<<"Decay e- received in Tube number "<<tubeNumber-1<<endl;
		}
	      }
	      totalPe += peForTube;
	      totalHit ++;
	      //}
	    //if(i%50 == 0) cout << "Hit angle wrt particle direction = " << relativeAngle << ", Phi = " << phiAngle*180/TMath::Pi() << " Z position = " << PMTpos[2] << endl;
	  }
	} // End of loop over Cherenkov hits
      if(verbose) cout << "Total Pe : " << totalPe << endl;
      if(!plotDigitized){
	TotalCharge[pmtType]->Fill(totalPe);
	TotalHit[pmtType]->Fill(totalHit);
	if(decayID.size()!=0 && dwall_decay > 0){
	  TotalHit_decay[pmtType]->Fill(totalHit_decay);
	}
	TotalChargeXdWall[pmtType]->Fill(dwall,totalPe);
	int decay_hitpmt = 0;
	int decay_hitpmt_independent = 0;
	for (int p=0; p< nPMTs ; p++)
	  {
	    int nh = totalHitPerPMT[p];
	    HitPerPMT[pmtType]->Fill(nh);
	    int nh_decay = totalHitPerPMT_decay[p];
	    if(decayID.size()!=0 && dwall_decay > 0){
	      HitPerPMT_decay[pmtType]->Fill(nh_decay);
	      if(nh_decay > 0){ 
		//cout<<"Check: Decay e- received in Tube number "<<p<<endl;
		//cout<<"Number of hits in this PMT = "<<nh<<" vs only decay-e hits = "<<nh_decay<<endl;
		if(nh == nh_decay) decay_hitpmt_independent++;
		decay_hitpmt ++;
	      }
	    }
	    if(nh>=2){
	      //sort vector by time
	      sort(timePerPMT[p].begin(),timePerPMT[p].end());
	      sort(timePerPMT_decay[p].begin(),timePerPMT_decay[p].end());
	      //Get earlier hit.
	      double start_hit = timePerPMT[p][0];
	      //Only for decay-e:
	      if(nh_decay!=0 && (nh_decay != nh)){//Ensure that there is one non-decay e hit also.
		double tdif_decay = timePerPMT_decay[p][0] - start_hit;
		if(decayID.size()!=0 && dwall_decay > 0) TimeDifferenceHitPerPMT_decay[pmtType]->Fill(tdif_decay);
	      }
	      for(int ihit = 1; ihit < nh; ihit++){
		double tdif = timePerPMT[p][ihit] - start_hit;
		//cout<<
		TimeDifferenceHitPerPMT[pmtType]->Fill(tdif);
    	      }
	    }
	  }
	if(decayID.size()!=0 && dwall_decay > 0){
	  //int decay_hitpmt_nonshared = decay_hitpmt - decay_hitpmt_shared;
	  TotalHit_decay_nonshared[pmtType]->Fill(decay_hitpmt_independent);
	}
      }
    }
    // Get the number of digitized hits
    // Loop over sub events


    
    if(verbose) cout << "DIGITIZED HITS:" << endl;
    //for (int index = 0 ; index < wcsimrootsuperevent->GetNumberOfEvents(); index++) 
    //{
    //wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
    //if(verbose) cout << "Sub event number = " << index << "\n";
    
    //int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
    //if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);

    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(separatedTriggers){
	if(triggerInfo2.size()!=0 && pmtType==0) continue;
	if(triggerInfo.size()!=0 && pmtType==1) continue;
      }
      if(verbose) cout << "PMT Type = " << pmtType << endl;
      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray;
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();
      
      double particleRelativePMTpos[3];
      double totalPe = 0;
      int totalHit = 0;

      int nhits;
      if(pmtType == 0) nhits = ncherenkovdigihits;
      else nhits = ncherenkovdigihits2;
      
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (i=0; i< nhits ; i++)
	{
	  TObject *Hit;
	  if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	  else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

	  WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
	    dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
	  
	  int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId();
	  double peForTube      = wcsimrootcherenkovdigihit->GetQ();

	  //int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  //int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	  //int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	  WCSimRootPMT pmt;
	  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
	  else pmt  = geo->GetPMT(tubeNumber-1,true); 
	  
	  double PMTpos[3];
	  double PMTdir[3];
	  for(int j=0;j<3;j++){
	    PMTpos[j] = pmt.GetPosition(j);
	    PMTdir[j] = pmt.GetOrientation(j);
	  }
	  
	  if ( i < 10 ) // Only print first XX=10 tubes
	    {
	      if(verbose) printf("Total pe: %1.1f times( ",peForTube);
	      /*	      for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
		{
		  WCSimRootCherenkovHitTime HitTime = 
		    dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(j));
		  
		  if(verbose) printf("%6.2f ", HitTime.GetTruetime() );	     
		  }*/
	      if(verbose) cout << ")" << endl;
	      if(verbose) std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << endl;
	    }
	  
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
	  double phiAngle;
	  if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
	  else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));
	  
	  if(Gamma){//ProduceDirectionPDF the profile at the gamma conversion point, since it takes ~36 cm to convert/be visible
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStop[j];
	  }
	  else{
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStart[j];
	  }
	  
	  double vDir[3];double vOrientation[3];
	  double vDir2D[3];double vOrientation2D[3];
	  for(int j=0;j<3;j++){
	    vDir[j] = particleRelativePMTpos[j];
	    vOrientation[j] = PMTdir[j];
	    if(j<2){
	      vDir2D[j] = particleRelativePMTpos[j];
	      vOrientation2D[j] = PMTdir[j];
	    }	  
	  }
	  double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
	  double tof = Norm*1e-2/(cvacuum/nindex);
	  double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
	  for(int j=0;j<3;j++){
	    vDir[j] /= Norm;
	    vOrientation[j] /= NormOrientation;
	  }
	  double Norm2D = TMath::Sqrt(vDir2D[0]*vDir2D[0]+vDir2D[1]*vDir2D[1]);
	  double NormOrientation2D = TMath::Sqrt(vOrientation2D[0]*vOrientation2D[0]+vOrientation2D[1]*vOrientation2D[1]);
	  for(int j=0;j<2;j++){
	    vDir2D[j] /= Norm2D;
	    vOrientation2D[j] /= NormOrientation2D;
	  }
	  double HorizontalScalarProduct = TMath::ACos(vDir2D[0]*vOrientation2D[0] + vDir2D[1]*vOrientation2D[1])*180./TMath::Pi();
	  double relativeAngle = TMath::ACos(vDir[0]*particleDir[0]+vDir[1]*particleDir[1]+vDir[2]*particleDir[2])*180./TMath::Pi();
	  //WCSimRootCherenkovHitTime HitTime = dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(timeArrayIndex));//Assumes that the earliest time the PMT is hit is in the first element -> Everything is ordered.
	  //double time = HitTime.GetTruetime();
	  double time = wcsimrootcherenkovdigihit->GetT();
	 
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////

	  double Theta = 0;
	  double Phi = 0;
	  double weight = 1;
	  double distance_pmt_vertex=0;
	  int mPMT_number = pmt.GetmPMTNo();
	  int pmt_number_in_mPMT = pmt.GetmPMT_PMTNo();

	  //TMP
	    double dist_tmp =0;
	    for(int i=0;i<3;i++) dist_tmp+=pow(particleStart[i]-pmt.GetPosition(i),2);
	    dist_tmp=TMath::Sqrt(dist_tmp);
	    //if(dist_tmp>500) continue;
	    //if(dist_tmp>1000) continue;
	    //

	  if(produceDirectionPDF){
	    //double * vPMTvertex = new double[3];vectorVertexPMT(geo,particleStart,pmt,tubeNumber,pmtType,vPMTvertex);
	    if(thetaOnly){
	      Theta=findDirectionTheta(particleStart,tubeNumber,pmtType,pmt_number_in_mPMT);
	    }
	    else{
	      if(pmtType==1){
		double vPMTvertex[4];
		vectorVertexPMT(particleStart,tubeNumber,pmtType,vPMTvertex,pmt_number_in_mPMT);
		Phi=vPMTvertex[0];
		Theta=vPMTvertex[1];
		distance_pmt_vertex=vPMTvertex[2];
		weight=vPMTvertex[3];
	      }
	    }
	  }
	  //std::cout << "Relative angle = " << relativeAngle << ", time = " << time << endl;
	  
	  if(plotDigitized){
	    //plotBQ(PMTpos[3], PMTdir[3],phiAngle,relativeAngle,peForTube,pmtType);
	    if(TMath::Abs(PMTpos[2]) < 4500){//Temp
	      if(phiAngle < TMath::Pi()/2 || phiAngle > 3*TMath::Pi()/2) ChargeProfile2D_onlyFront[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2D[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2DTop[pmtType]->Fill(PMTpos[0],PMTpos[1],peForTube);
	      ChargeProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],peForTube);
	      ChargeProfile2DRelativeAngle[pmtType]->Fill(relativeAngle,PMTpos[2],peForTube);
	      ChargeProfile2D_event[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);

	      ChargeProfileSpread[pmtType]->Fill(ev,relativeAngle,peForTube);
	      ChargeProfile[pmtType]->Fill(relativeAngle,peForTube);
	      HitProfile[pmtType]->Fill(relativeAngle,1);
#ifdef NEUTRONCAPTURE
	      //In that case, we do not care about the trigger. We just want to check in the whole time window.
	      time+=triggerTime[pmtType];
#endif
	    
	      if(ev%100==0 && i<10){
		cout<<"PMT type = "<<pmtType<<", pmt number = "<<tubeNumber<<", charge = "<<peForTube<<", time = "<<time<<", after shift by trigger = "<<time+triggerTime[pmtType]<<", charge in the PMT = "<<peForTube<<endl;
	      }
		
	      
	      TimeProfile[pmtType]->Fill(time,peForTube);
	      TimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],peForTube);
	      TimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,peForTube);

	      int pmtGroup=groupPMTs(pmtType,pmt_number_in_mPMT);
	      angularResponseXPMT_hits[pmtType][pmtGroup]->Fill(distance_pmt_vertex,Phi,Theta,weight);
	      //angularResponseXPMT_hits[pmtType]->Fill(pmtGroup,Phi,Theta,weight);
	      distanceXweight_hits[pmtType]->Fill(distance_pmt_vertex,TMath::Log(weight));
	      //angularResponseXPMT_allPMTs[pmtType]->Fill(pmt_number_in_mPMT,Phi*180/TMath::Pi(),Theta*180/TMath::Pi(),1);
	      //ChargeProfile2DXPMT[pmtType]->Fill(pmt_number_in_mPMT,PMTpos[1],PMTpos[2],peForTube);
	      //cout<<"DIST = "<<distance_pmt_vertex<<", "<<weight<<endl;
	      HitTimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],1);
	      HitTimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,1);
	      
	      //if(verbose && i<10){
	      //cout<<"Hit time="<<time<<", tof="<<tof<<", vertex time=0"<<endl;
	      //cout<<"Vertex pos = "<<particleStart[0]<<", "<<particleStart[1]<<", "<<particleStart[2]<<", distance PMT to vertex="<<Norm<<"cm"<<", pmt pos X="<<PMTpos[0]<<","<<PMTpos[1]<<","<<PMTpos[2]<<endl;
	      //}
	      ChargePerPMT[pmtType]->Fill(peForTube);

	      totalPe += peForTube;
	      totalHit ++;
	    }
	    //if(i%50 == 0) cout << "Hit angle wrt particle direction = " << relativeAngle << ", Phi = " << phiAngle*180/TMath::Pi() << " Z position = " << PMTpos[2] << endl;
	  }
	} // End of loop over Cherenkov hits
      if(plotDigitized){
	if(verbose) cout << "Total Pe : " << totalPe << ", total hit : " << totalHit << endl;
	TotalCharge[pmtType]->Fill(totalPe);
	TotalHit[pmtType]->Fill(totalHit);
	TotalChargeXdWall[pmtType]->Fill(dwall,totalPe);
      }
    }
    // Get the number of digitized hits
    // Loop over sub events
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(eventDisplay_closeWall && dwall<200 && dwall>100){
	ChargeProfile2D_event[pmtType]->SaveAs(Form("trashPlots/ChargeProfile2D_closewall_evt%d_pmtType%d.root",stopCounterPlots,pmtType));
	stopCounterPlots++;
	if(stopCounterPlots >= nPlotsStopping) return 0;
      }
      ChargeProfile2D_event[pmtType]->Reset();
    }
    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    if(hybrid) wcsimrootsuperevent2->ReInitialize();
    
    } // End of loop over events
  //  TCanvas c1("c1"); 
  /*
  float win_scale = 0.75;
  int n_wide(2);
  int n_high(2);
  TCanvas* c1 = new TCanvas("c1", "First canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
  c1->Draw();
  c1->Divide(2,2);
  c1->cd(1); hvtx0->Draw();
  c1->cd(2); hvtx1->Draw();
  c1->cd(3); hvtx2->Draw();
  c1->cd(4); h1->Draw();
  */
  std::cout<<"num_trig "<<num_trig<<"\n";
  if(outfilename==NULL) sprintf(outfilename,"out.root");
  
  TFile * outfile = new TFile(outfilename,"RECREATE");
  cout<<"File "<<outfilename<<" is open for writing"<<endl;
  hvtx0->Write();
  hvtx1->Write();
  hvtx2->Write();
  hdir0->Write();
  hdir1->Write();
  hdir2->Write();
  for(int i=0;i<nPMTtypes;i++){
    
    ChargeProfile2D_onlyFront[i]->Write();
    
    ChargeProfile2D[i]->Write();
    
    ChargeProfile2DTop[i]->Write();
    
    ChargeProfile2DCylindricCoordinate[i]->Write();
    
    ChargeProfile2DRelativeAngle[i]->Write();
    
    ChargeProfile[i]->Write();
    
    HitProfile[i]->Write();
    
    TimeProfile[i]->Write();
    
    TimeTOFProfile[i]->Write();
    
    TimeTOFProfileXTOF[i]->Write();
    
    HitTimeTOFProfile[i]->Write();
    
    HitTimeTOFProfileXTOF[i]->Write();
    
    TriggerTimeProfile[i]->Write();

    ChargePerPMT[i]->Write();
    ChargePerPMT_decay[i]->Write();
    
    TotalCharge[i]->Write();
    
    TotalHit[i]->Write();
    TotalHit_decay[i]->Write();
    
    TotalChargeXdWall[i]->Write();//Fill(dWall,totalPe);
    
    HitPerPMT[i]->Write();
    HitPerPMT_decay[i]->Write();
    TotalHit_decay_nonshared[i]->Write();

    TimeDifferenceHitPerPMT[i]->Write();
    TimeDifferenceHitPerPMT_decay[i]->Write();
    //ChargeProfile2DXPMT[i]->Write();

    for(int p=0;p<nGroupsPMTs;p++){
      angularResponseXPMT_hits[i][p]->Write();    
      angularResponseXPMT_allPMTs[i][p]->Write();
    }
    distanceXweight_hits[i]->Write();
    
  ChargeProfileSpread[i]->Write();
  }
  neutrinoEnergy->Write();
  outfile->Close();
  
  return 0;
 }
