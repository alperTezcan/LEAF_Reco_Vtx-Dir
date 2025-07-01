#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMinuit.h>
#include <TFitter.h>

using namespace std;

//use array of c++ 11 as we cannot create a vector of traditional arrays e.g. double[5].
const int ntuple=5;
vector < array < double,ntuple > > dataID;

double CredibleInterval=0.68;

double ComputeR(double x, double y){
	return TMath::Sqrt(x*x + y*y);
}

double ComputeDwall(double r, double z, double R, double Z){
	double absz = TMath::Abs(z);
	double signflg = 1.;
	if(absz > Z || r > R) signflg = -1.;
	double distz = TMath::Abs(Z - absz);
	double distr = TMath::Abs(R - r);
	return signflg * fmin(distr, distz);
}

double ComputeTowall(double x, double y, double z, double dx, double dy, double dz, double R, double Z){
	double tbar = 1000000000.; // Placeholder, impossibly big value. This way, even dx = dy = dz returns something.
	double tcap = 1000000000.; // Placeholder, impossibly big value. This way, even dx = dy = dz returns something.
	if(dx!=0 || dy!=0){ 
		// Distance for hitting the barrel, assuming there is some radial direction. Found by solving the quadratic X + t*DX = R.
		tbar = 1/(dx*dx + dy*dy) * (-1 * (x*dx + y*dy) + TMath::Sqrt((x*dx + y*dy)*(x*dx + y*dy) + (R*R - x*x - y*y))); 
	}
	if(dz == 0){
		// Directly returns the radial distance if particle is on the radial plane
		return tbar;
	}else{
		// Distance for hitting one of the caps, analytical result of (signOf(dz) * Z - z)/dz
	  tcap = Z/TMath::Abs(dz) - z/dz;
	}
	return (tcap > tbar ? tbar:tcap);
}

bool migrateFromOVtoFV(double dwall, double dwall_reconstructed, double dwall_FV_definition, bool verbose){
  if(verbose) cout<<"FV cut is set at a distance = "<<dwall_FV_definition<<", Event dwall = "<<dwall<<", reconstructed dwall = "<<dwall_reconstructed<<endl;
  // Event is truly out of FV AND reconstructed inside the FV
  return ((dwall<dwall_FV_definition) && (dwall_reconstructed >= dwall_FV_definition));
}

double vRes(TH1D * hvertex, double CI = CredibleInterval){
  double integralSoFar = 0;
  double integralTotal = hvertex->Integral();
  double integral = 0;
  int ibin1D=1;
  while(integral<CI && ibin1D<1e4){
    integralSoFar += hvertex->GetBinContent(ibin1D);
    integral = integralSoFar/(integralTotal==0?1:integralTotal);
    ibin1D++;
  }
  return hvertex->GetBinCenter(ibin1D);
}

double vResAndError(TH1D * hvertex, double CI = CredibleInterval, bool verbose = false){
  double integralSoFar = 0;
  double integralTotal = hvertex->Integral();
  if(integralTotal==0) return 0;
  double integral = 0;
  int ibin1D=1;
  double error=0;
  while(integral<CI && ibin1D<1e4){
    integralSoFar += hvertex->GetBinContent(ibin1D);
    integral = integralSoFar/(integralTotal==0?1:integralTotal);
    ibin1D++;
  }
  double res=hvertex->GetBinCenter(ibin1D);
  if(verbose) cout<<"res="<<res<<", "<<integralSoFar<<", "<<integralTotal<<endl;
  error=integral*TMath::Sqrt(pow(1/integralSoFar,2)+pow(1/integralTotal,2));
  if(verbose) cout<<"Error on integral = "<<error<<endl;
  double error_vertex_inf=vRes(hvertex,CI-error);
  double error_vertex_sup=vRes(hvertex,CI+error);
  double error_full=TMath::Max(TMath::Abs(res-error_vertex_inf),TMath::Abs(res-error_vertex_sup));
  if (verbose) cout<<"Error on position inf = "<<error_vertex_inf<<", "<<error_vertex_sup<<", full = "<<error_full<<endl;
  error_full+=hvertex->GetBinWidth(ibin1D);//Add bin width
  return error_full;
}

int vertexResolution(double energy, double vertex_resolution, double efficiency, double misIDFV, bool verbose, bool wrt_to_ROOT, char * fname, char * ofname, bool monoEnergy=true, double photocov=20.2150576375662, double darkrate=4.2){
	//Tank constants
  double tankRadius = 3240.0;//In cm
  double tankHeight = 3287.55;//In cm
  
	//Cut valuesvertexRe
  bool FVCut = true;
  double FVCutValue = 0;//in cm
  bool ExtCut = false;
  double ExtCutvalue = 1200;//in cm  

	// The input file
	TFile * file = new TFile(fname, "read");
	if(verbose) cout<<"We will try to open file = "<<file->GetName()<<endl;
  file->cd();

	// The TTree linking to variables
	TTree *bstree = (TTree*)file->Get("Reduced");
  int nentries = bstree->GetEntries();
  if(verbose) cout<<"Number of entries = "<<nentries<<endl;

  //Preparation of host variables
  std::vector<float> *true_vertex_X = 0;
	std::vector<float> *true_vertex_Y = 0;
	std::vector<float> *true_vertex_Z = 0;
	std::vector<float> *true_vertex_T = 0;

  std::vector<float> *true_dir_X = 0;
	std::vector<float> *true_dir_Y = 0;
	std::vector<float> *true_dir_Z = 0;

  // double lf_time = 0.;
  double lf_vertex[4] = {0.};
  double lf_dir[3] = {0.};
  double lf_wall = 0;
  int lf_n50[3] = {0};
  double lf_good = 0.;
  double lf_dirKS[3] = {0.};
  int lf_intime = 0;
  int eventId = 0, triggerId = 0, ID_hits_200 = 0, ID_hits = 0;

	//True variables
  bstree->SetBranchAddress("true_vertex_X", &true_vertex_X);
  bstree->SetBranchAddress("true_vertex_Y", &true_vertex_Y);
  bstree->SetBranchAddress("true_vertex_Z", &true_vertex_Z);
  bstree->SetBranchAddress("true_vertex_T", &true_vertex_T);
  bstree->SetBranchAddress("true_dir_X", &true_dir_X);
  bstree->SetBranchAddress("true_dir_Y", &true_dir_Y);
  bstree->SetBranchAddress("true_dir_Z", &true_dir_Z);  

  //LEAF reconstructed variables			     
  //bstree->SetBranchAddress("lf_time", lf_time);
  bstree->SetBranchAddress("lf_vertex", lf_vertex);
  bstree->SetBranchAddress("lf_dir", lf_dir);
  bstree->SetBranchAddress("lf_good", &lf_good);
  bstree->SetBranchAddress("lf_dirKS", lf_dirKS);
  bstree->SetBranchAddress("lf_wall",&lf_wall);
  bstree->SetBranchAddress("lf_n50", lf_n50);

  //Pure detector raw variables
  bstree->SetBranchAddress("eventId",&eventId);
  bstree->SetBranchAddress("triggerId",&triggerId);
  bstree->SetBranchAddress("ID_hits", &ID_hits);
  bstree->SetBranchAddress("ID_hits_200",&ID_hits_200);

	// Declaration of output histograms
	//Vertex related
  TH2D * hvertexLongitudinalXTime;
  TH1D * hvertexLongitudinal;
	TH1D * hvertexOrtho;
  TH1D * hvertex;
  TH1D * hradius;
  double nVertexFound = 0., nEvents = 0., nEventsTotal = 0.;
  TH2D * hvertex2D;
  TH2D * htruevertex2D;
  TH2D * hvertexXdwall;
  TH2D * hvertexXtowall;
  TH2D * hdwallXdwallrec; 

  //Nhits related
  TH1D * hnhits;
  TH2D * hdwallXnhits;
  TH2D * htowallXnhits;

  //Direction related
  TH1D * hdirection;
  TH1D * hOpeningAngle;
  TH1D * hdirX;
  TH1D * hdirY;
  TH1D * hdirZ;
  
  //Goodness of fit & dirKS related
  TH1D * hgood;
  TH1D * hdirKS;
  
  //Mix of previous variables
  TH2D * hvertexXnhits_intime;
  TH2D * hvertexXgood;
  TH2D * hvertexXnhits;
 
  //Cut related
  TH1D * hvertexMisIDFV; TH1D * hvertexMisIDFV_total;
  TH1D * hvertexMisIDExtCut; TH1D * hvertexMisIDExtCut_total;
  double nMisIDExtCut = 0., nMisIDExtCut_total = 0.;
  double nMisIDFV = 0., nMisIDFV_total = 0.;

	//Alper's Extra Plots
	TH1D * hvtx0;
	TH1D * hvtx1;
	TH1D * hvtx2;
	TH1D * hvtx3;

	TH2D * hvtx0Xdwall;
	TH2D * hvtx1Xdwall;
	TH2D * hvtx2Xdwall;
	TH2D * hvtx3Xdwall;

	TH2D * hvtx0Xtowall;
	TH2D * hvtx1Xtowall;
	TH2D * hvtx2Xtowall;
	TH2D * hvtx3Xtowall;
	
	TH2D * hvertexLongitudinalXdwall;
	TH2D * hvertexOrthoXdwall;
	TH2D * hradiusXdwall;

	TH2D * hvertexLongitudinalXtowall;
	TH2D * hvertexOrthoXtowall;
	TH2D * hradiusXtowall;

	// Binning choice
	int nbin_vertex = 6000;//Number of bins used for vertex resolution
  double max_vertex = 6000;
  int nbin_direction = 180;//Number of bins used for angle resolution
  double max_direction = 180;
  int nbin_nhits = 800;
  double max_nhits = 800;
  int nbin_time = 400;
  double max_time = 50;
  int nbin_distance = 200;
  double min_distance = 0.;
  double max_distance = 10000;
  double nbin_tanksize = 760;
  double min_tanksize = -3800;
  double max_tanksize = 3800;
  int nbin_good = 100;//Number of bins used for vertex resolution
  double max_good = 1000;
  int nbin_dirKS = 50;//Number of bins used for vertex resolution
  double max_dirKS = 1;

	// Histogram initialization
	//Vertex related
  hvertexLongitudinalXTime = new TH2D(Form("hvertexLongitudinalXTime"),"",nbin_vertex,0,max_vertex,nbin_time,-max_time,max_time);
  hvertexLongitudinalXTime->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexLongitudinalXTime->GetYaxis()->SetTitle("dwall (cm)");			     

  hvertexLongitudinal = new TH1D(Form("hvertexLongitudinal"),"",nbin_vertex,0,max_vertex);
  hvertexLongitudinal->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");

  hvertexOrtho = new TH1D(Form("hvertexOrtho"),"",nbin_vertex,0,max_vertex);
  hvertexOrtho->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");

  hvertex = new TH1D(Form("hvertex"),"",nbin_vertex,0,max_vertex);
  hvertex->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");

  hradius = new TH1D(Form("hradius"),"",nbin_vertex,0,max_vertex*max_vertex);
  hradius->GetXaxis()->SetTitle("R2 (cm2)");
 
  hvertex2D = new TH2D(Form("hvertex2D"),"",nbin_tanksize,min_tanksize,max_tanksize,nbin_tanksize,min_tanksize,max_tanksize);
  hvertex2D->GetXaxis()->SetTitle("x position (cm)");
  hvertex2D->GetYaxis()->SetTitle("y position (cm)");

  htruevertex2D = new TH2D(Form("htruevertex2D"),"",nbin_tanksize,min_tanksize,max_tanksize,nbin_tanksize,min_tanksize,max_tanksize);
  htruevertex2D->GetXaxis()->SetTitle("x position (cm)");
  htruevertex2D->GetYaxis()->SetTitle("y position (cm)");

  hdwallXdwallrec = new TH2D(Form("hdwallXdwallrec"),"",nbin_distance,0,max_distance,20,-1,1);
  hdwallXdwallrec->GetXaxis()->SetTitle("dwall true (cm)");
  hdwallXdwallrec->GetYaxis()->SetTitle("dwall rec  - dwall true / dwall true");

  hvertexXdwall = new TH2D(Form("hvertexXdwall"),"",nbin_vertex,0,max_vertex,nbin_distance,0,max_distance);
  hvertexXdwall->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXdwall->GetYaxis()->SetTitle("dwall (cm)");

  hvertexXtowall = new TH2D(Form("hvertexXtowall"),"",nbin_vertex,0,max_vertex,nbin_distance,0,max_distance);
  hvertexXtowall->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXtowall->GetYaxis()->SetTitle("towall (cm)");
  
  //Nhits related
  hnhits = new TH1D(Form("hnhits"),"",nbin_nhits,0,max_nhits);
  hnhits->GetXaxis()->SetTitle("number of hits");

  hdwallXnhits = new TH2D(Form("hdwallXnhits"),"",nbin_distance,0,max_distance,nbin_nhits,0,max_nhits);
  hdwallXnhits->GetXaxis()->SetTitle("dwall (cm)");
  hdwallXnhits->GetYaxis()->SetTitle("number of hits");

  htowallXnhits = new TH2D(Form("htowallXnhits"),"",nbin_distance,0,max_distance,nbin_nhits,0,max_nhits);
  htowallXnhits->GetXaxis()->SetTitle("towall (cm)");
  htowallXnhits->GetYaxis()->SetTitle("number of hits");

  //Direction related
  hdirection = new TH1D(Form("hdirection"),"",nbin_direction,0,max_direction);
  hdirection->GetXaxis()->SetTitle("angle from rec. to mc. direction (#circ)");

  hOpeningAngle = (TH1D*) file->Get("hOpeningAngle");

  hdirX = new TH1D("hdirX","",100,-1,1);
  hdirX->GetXaxis()->SetTitle("(rec dir - mc dir)/ mc dir for X");
  hdirY = new TH1D("hdirY","",100,-1,1);
  hdirY->GetXaxis()->SetTitle("(rec dir - mc dir)/ mc dir for Y");
  hdirZ = new TH1D("hdirZ","",100,-1,1);
  hdirZ->GetXaxis()->SetTitle("(rec dir - mc dir)/ mc dir for Z");


  //Goodness of fit & dirKS related
  hgood = new TH1D(Form("hgood"),"",nbin_good,0,max_good);
  hgood->GetXaxis()->SetTitle("good of fit");
  hdirKS = new TH1D(Form("hdirKS"),"",nbin_dirKS,0,max_dirKS);
  hdirKS->GetXaxis()->SetTitle("dirKS");

  //Mix of previous variables
  hvertexXnhits_intime = new TH2D(Form("hvertexXnhits_intime"),"",nbin_vertex,0,max_vertex,nbin_nhits,0.,max_nhits);
  hvertexXnhits_intime->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXnhits_intime->GetYaxis()->SetTitle("number of hits");

  hvertexXgood = new TH2D(Form("hvertexXgoddness"),"",nbin_vertex,0,max_vertex,nbin_good,0.,max_good);
  hvertexXgood->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXgood->GetYaxis()->SetTitle("good of fit");

  hvertexXnhits = new TH2D(Form("hvertexXnhits"),"",nbin_vertex,0,max_vertex,nbin_nhits,0,max_nhits);
  hvertexXnhits->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXnhits->GetYaxis()->SetTitle("number of hits");

  //Cut related
  hvertexMisIDFV = new TH1D(Form("hvertexMisIDFV"),"",20,0,200);
  hvertexMisIDFV->GetXaxis()->SetTitle("FV distance to the wall");
  hvertexMisIDFV->GetYaxis()->SetTitle("Probability of migration from out to in FV");
  hvertexMisIDFV_total = new TH1D(Form("hvertexMisIDFV_total"),"",20,0,200);

  hvertexMisIDExtCut = new TH1D(Form("hvertexMisIDExtCut"),"",20,0,2000);
  hvertexMisIDExtCut->GetXaxis()->SetTitle("ExtCut distance to the wall");
  hvertexMisIDExtCut->GetYaxis()->SetTitle("Probability of migration from out to in ExtCut");
  hvertexMisIDExtCut_total = new TH1D(Form("hvertexMisIDExtCut_total"),"",20,0,2000);

	// Alper's plots
	hvtx0 = new TH1D(Form("hvtx0"), "", nbin_vertex, 0, max_vertex);
  hvtx0->GetXaxis()->SetTitle("X resolution (cm)");
	hvtx1 = new TH1D(Form("hvtx1"), "", nbin_vertex, 0, max_vertex);
  hvtx1->GetXaxis()->SetTitle("Y resolution (cm)");
	hvtx2 = new TH1D(Form("hvtx2"), "", nbin_vertex, 0, max_vertex);
  hvtx2->GetXaxis()->SetTitle("Z resolution (cm)");
	hvtx3 = new TH1D(Form("hvtx3"), "", nbin_vertex, 0, max_vertex);
  hvtx3->GetXaxis()->SetTitle("T resolution (cm)");

	hvtx0Xdwall = new TH2D(Form("hvtx0Xdwall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx0Xdwall->GetXaxis()->SetTitle("X resolution (cm)");
  hvtx0Xdwall->GetYaxis()->SetTitle("dwall (cm)");
	hvtx1Xdwall = new TH2D(Form("hvtx1Xdwall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx1Xdwall->GetXaxis()->SetTitle("Y resolution (cm)");
  hvtx1Xdwall->GetYaxis()->SetTitle("dwall (cm)");
	hvtx2Xdwall = new TH2D(Form("hvtx2Xdwall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx2Xdwall->GetXaxis()->SetTitle("Z resolution (cm)");
  hvtx2Xdwall->GetYaxis()->SetTitle("dwall (cm)");
	hvtx3Xdwall = new TH2D(Form("hvtx3Xdwall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx3Xdwall->GetXaxis()->SetTitle("T resolution (cm)");
  hvtx3Xdwall->GetYaxis()->SetTitle("dwall (cm)");

	hvtx0Xtowall = new TH2D(Form("hvtx0Xtowall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx0Xtowall->GetXaxis()->SetTitle("X resolution (cm)");
  hvtx0Xtowall->GetYaxis()->SetTitle("towall (cm)");
	hvtx1Xtowall = new TH2D(Form("hvtx1Xtowall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx1Xtowall->GetXaxis()->SetTitle("Y resolution (cm)");
  hvtx1Xtowall->GetYaxis()->SetTitle("towall (cm)");
	hvtx2Xtowall = new TH2D(Form("hvtx2Xtowall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx2Xtowall->GetXaxis()->SetTitle("Z resolution (cm)");
  hvtx2Xtowall->GetYaxis()->SetTitle("towall (cm)");
	hvtx3Xtowall = new TH2D(Form("hvtx3Xtowall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvtx3Xtowall->GetXaxis()->SetTitle("T resolution (cm)");
  hvtx3Xtowall->GetYaxis()->SetTitle("towall (cm)");

  hvertexLongitudinalXdwall = new TH2D(Form("hvertexLongitudinalXdwall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvertexLongitudinalXdwall->GetXaxis()->SetTitle("VTX Long");
  hvertexLongitudinalXdwall->GetYaxis()->SetTitle("dwall (cm)");  
	hvertexLongitudinalXtowall = new TH2D(Form("hvertexLongitudinalXtowall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvertexLongitudinalXtowall->GetXaxis()->SetTitle("VTX Long");
  hvertexLongitudinalXtowall->GetYaxis()->SetTitle("towall (cm)");

  hvertexOrthoXdwall = new TH2D(Form("hvertexOrthoXdwall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvertexOrthoXdwall->GetXaxis()->SetTitle("VTX Long");
  hvertexOrthoXdwall->GetYaxis()->SetTitle("dwall (cm)");  
	hvertexOrthoXtowall = new TH2D(Form("hvertexOrthoXtowall"), "", nbin_vertex, 0, max_vertex, nbin_distance, 0, max_distance);
  hvertexOrthoXtowall->GetXaxis()->SetTitle("VTX Long");
  hvertexOrthoXtowall->GetYaxis()->SetTitle("towall (cm)");

  hradiusXdwall = new TH2D(Form("hradiusXdwall"), "", nbin_vertex, 0, max_vertex*max_vertex, nbin_distance, 0, max_distance);
  hradiusXdwall->GetXaxis()->SetTitle("R2 (cm2)");
  hradiusXdwall->GetYaxis()->SetTitle("dwall (cm)");
  hradiusXtowall = new TH2D(Form("hradiusXtowall"), "", nbin_vertex, 0, max_vertex*max_vertex, nbin_distance, 0, max_distance);
  hradiusXtowall->GetXaxis()->SetTitle("R2 (cm2)");
  hradiusXtowall->GetYaxis()->SetTitle("towall (cm)");

	// Loop over each DAQ/true vertex event
	for (int entry=0; entry<nentries; entry++){
    bstree->GetEntry(entry);
    nEventsTotal++;
		// One DAQ/Sim event can contain several vertices, and several trigger per vertices.
    // Only "one vertex per event and one trigger for each" is implemented.
		float mc_vertex[3] = {true_vertex_X->at(0),true_vertex_Y->at(0),true_vertex_Z->at(0)};
    float mc_time = true_vertex_T->at(0);
    double lf_time = lf_vertex[3];
    float mc_dir[3] = {true_dir_X->at(0),true_dir_Y->at(0),true_dir_Z->at(0)};
    if(verbose){
			std::cout<<"############ True event ############"<<std::endl;
			std::cout<<"Vertex pos = "<<mc_vertex[0]<<", "<<mc_vertex[1]<<", "<<mc_vertex[2]<<endl;
      std::cout<<"Vertex time = "<<mc_time<<std::endl;
			std::cout<<"Particle dir = "<<mc_dir[0]<<", "<<mc_dir[1]<<", "<<mc_dir[2]<<endl;
			std::cout<<"############ Reco event ############"<<std::endl;
			std::cout<<"Vertex pos = "<<lf_vertex[0]<<", "<<lf_vertex[1]<<", "<<lf_vertex[2]<<endl;
      std::cout<<"Vertex time = "<<lf_time<<std::endl;
			std::cout<<"Particle dir PMT 1 = "<<lf_dir[0]<<", "<<lf_dir[1]<<", "<<lf_dir[2]<<endl;
			//std::cout<<"Particle dir PMT 2 = "<<lf_dir[1][0]<<", "<<lf_dir[1][1]<<", "<<lf_dir[1][2]<<endl;
			//std::cout<<"Particle dir PMT 3 = "<<lf_dir[2][0]<<", "<<lf_dir[2][1]<<", "<<lf_dir[2][2]<<endl;
			std::cout<<"####################################"<<std::endl<<std::endl;
    }
		// Calculation of dwall and towall
		double r_true_vertex = ComputeR(mc_vertex[0], mc_vertex[1]);
		double r_reco_vertex = ComputeR(lf_vertex[0], lf_vertex[1]);
		double dwall = ComputeDwall(r_true_vertex, mc_vertex[2], tankRadius, tankHeight);
		double dwall_reco = ComputeDwall(r_reco_vertex, lf_vertex[2], tankRadius, tankHeight);
		double towall = ComputeTowall(mc_vertex[0], mc_vertex[1], mc_vertex[2], mc_dir[0], mc_dir[1], mc_dir[2], tankRadius, tankHeight);
		double pmt_to_vertex = ComputeTowall(mc_vertex[0], mc_vertex[1], mc_vertex[2], (-1)*mc_dir[0], (-1)*mc_dir[1], (-1)*mc_dir[2], tankRadius, tankHeight);
		double pmt_to_vertex_reco = ComputeTowall(lf_vertex[0], lf_vertex[1], lf_vertex[2], (-1)*mc_dir[0], (-1)*mc_dir[1], (-1)*mc_dir[2], tankRadius, tankHeight);
	
		// Analysis cuts
		if(FVCut && lf_wall<FVCutValue) continue;
    if(ExtCut && pmt_to_vertex<ExtCutvalue) continue;
    nEvents++;

		// Vertex resolution
		//Simple vertex resolution
		double vertex_distance = sqrt(pow(lf_vertex[0]-mc_vertex[0],2) + pow(lf_vertex[1]-mc_vertex[1],2) + pow(lf_vertex[2]-mc_vertex[2],2));
    hvertex->Fill(vertex_distance);
    hradius->Fill(pow(r_reco_vertex,2));
    hnhits->Fill(ID_hits);
    hvertexXnhits->Fill(vertex_distance,ID_hits);
    hdwallXnhits->Fill(dwall,ID_hits);
    htowallXnhits->Fill(towall,ID_hits);
    hvertexXnhits_intime->Fill(vertex_distance,lf_intime);
    nVertexFound++;

    //Vertex resolution along particle direction and transverse
    double vector_director[3];double vector_director_norm;
    for(int i=0;i<3;i++){
			vector_director[i]=lf_vertex[i]-mc_vertex[i];
			vector_director_norm+=pow(vector_director[i],2);
    }
    vector_director_norm = TMath::Sqrt(vector_director_norm);
    double vector_director_longitudinal = 0;//Component parralel to the direction of the particle
    double vector_director_ortho = 0;//Component ortho to the direction of the particle
    for(int i=0;i<3;i++){
			vector_director_longitudinal += vector_director[i]*mc_dir[i];
    }
    double angle=TMath::Abs(TMath::ACos(vector_director_longitudinal/vector_director_norm));
    vector_director_ortho = vector_director_norm*TMath::Sin(angle);
    hvertexLongitudinal->Fill(TMath::Abs(vector_director_longitudinal));
    hvertexOrtho->Fill(vector_director_ortho);
    hvertexLongitudinalXTime->Fill(TMath::Abs(vector_director_longitudinal),lf_time-mc_time);


    //Vertex resolution vs dwall and towall
		hvertexXdwall->Fill(vertex_distance,dwall);
    hvertexXtowall->Fill(vertex_distance,towall);
    hdwallXdwallrec->Fill(dwall,(dwall_reco-dwall)/dwall);

		//Alper's plots
		hvtx0->Fill(TMath::Abs(lf_vertex[0]-mc_vertex[0]));
		hvtx1->Fill(TMath::Abs(lf_vertex[1]-mc_vertex[1]));
		hvtx2->Fill(TMath::Abs(lf_vertex[2]-mc_vertex[2]));
		hvtx3->Fill(TMath::Abs(lf_time-mc_time));
		
		hvtx0Xdwall->Fill(TMath::Abs(lf_vertex[0]-mc_vertex[0]), dwall);
		hvtx1Xdwall->Fill(TMath::Abs(lf_vertex[1]-mc_vertex[1]), dwall);
		hvtx2Xdwall->Fill(TMath::Abs(lf_vertex[2]-mc_vertex[2]), dwall);
		hvtx3Xdwall->Fill(TMath::Abs(lf_time-mc_time), dwall);
		
		hvtx0Xtowall->Fill(TMath::Abs(lf_vertex[0]-mc_vertex[0]), towall);
		hvtx1Xtowall->Fill(TMath::Abs(lf_vertex[1]-mc_vertex[1]), towall);
		hvtx2Xtowall->Fill(TMath::Abs(lf_vertex[2]-mc_vertex[2]), towall);
		hvtx3Xtowall->Fill(TMath::Abs(lf_time-mc_time), towall);

		hradiusXdwall->Fill(pow(r_reco_vertex,2), dwall);
		hvertexLongitudinalXdwall->Fill(TMath::Abs(vector_director_longitudinal), dwall);
    hvertexOrthoXdwall->Fill(vector_director_ortho, dwall);
		hradiusXtowall->Fill(pow(r_reco_vertex,2), towall);
		hvertexLongitudinalXtowall->Fill(TMath::Abs(vector_director_longitudinal), towall);
    hvertexOrthoXtowall->Fill(vector_director_ortho, towall);

		for(int ibinx=1;ibinx<=hvertexMisIDFV->GetNbinsX();ibinx++){
			double dwall_FV_definition = hvertexMisIDFV->GetBinLowEdge(ibinx) + hvertexMisIDFV->GetBinWidth(ibinx);
			if(dwall<dwall_FV_definition){//check only if the event is out of FV
			  bool misID = migrateFromOVtoFV(dwall, dwall_reco, FVCutValue, verbose);
	    	if(misID) hvertexMisIDFV->Fill(dwall_FV_definition-hvertexMisIDFV->GetBinWidth(ibinx));
  	  	hvertexMisIDFV_total->Fill(dwall_FV_definition-hvertexMisIDFV->GetBinWidth(ibinx));
			}
    }
    
    if(dwall<FVCutValue){
			bool migration = migrateFromOVtoFV(dwall,dwall_reco,FVCutValue, verbose);
			if(migration) nMisIDFV++;
			nMisIDFV_total++;
    }

    for(int ibinx=1; ibinx<=hvertexMisIDExtCut->GetNbinsX(); ibinx++){
			double ExtCut_definition=hvertexMisIDExtCut->GetBinLowEdge(ibinx)+hvertexMisIDExtCut->GetBinWidth(ibinx);
			bool misID=migrateFromOVtoFV(pmt_to_vertex, pmt_to_vertex_reco, ExtCut_definition, verbose);
			if(misID) hvertexMisIDExtCut->Fill(ExtCut_definition-hvertexMisIDExtCut->GetBinWidth(ibinx));
			hvertexMisIDExtCut_total->Fill(ExtCut_definition-hvertexMisIDExtCut->GetBinWidth(ibinx));
		}
		bool migration = migrateFromOVtoFV(pmt_to_vertex, pmt_to_vertex_reco, ExtCutvalue, verbose);
    if(migration) nMisIDExtCut++;
    nMisIDExtCut_total++;

    // Direction resolution
		//Simple direction resolution
    double direction_angle = TMath::ACos(lf_dir[0]*mc_dir[0] + lf_dir[1]*mc_dir[1] + lf_dir[2]*mc_dir[2])*180/TMath::Pi();
    hdirection->Fill(direction_angle);
    hdirX->Fill((lf_dir[0] - mc_dir[0])/mc_dir[0]);
    hdirY->Fill((lf_dir[1] - mc_dir[1])/mc_dir[1]);
    hdirZ->Fill((lf_dir[2] - mc_dir[2])/mc_dir[2]);
    
    if(verbose){
	    std::cout<<"####################################"<<std::endl;
      std::cout<<"Direction difference rec vs true = "<<direction_angle<<"degrees"<<std::endl;
      std::cout<<"############ True event ############"<<std::endl;
	    std::cout<<"Vertex pos = "<<mc_vertex[0]<<", "<<mc_vertex[1]<<", "<<mc_vertex[2]<<endl;
	    std::cout<<"Particle dir = "<<mc_dir[0]<<", "<<mc_dir[1]<<", "<<mc_dir[2]<<endl;
	    std::cout<<"############ Reco event ############"<<std::endl;
	    std::cout<<"Vertex pos = "<<lf_vertex[0]<<", "<<lf_vertex[1]<<", "<<lf_vertex[2]<<endl;
	    std::cout<<"Particle dir PMT 1 = "<<lf_dir[0]<<", "<<lf_dir[1]<<", "<<lf_dir[2]<<endl;
	    //std::cout<<"Particle dir PMT 2 = "<<lf_dir[1][0]<<", "<<lf_dir[1][1]<<", "<<lf_dir[1][2]<<endl;
	    //std::cout<<"Particle dir PMT 3 = "<<lf_dir[2][0]<<", "<<lf_dir[2][1]<<", "<<lf_dir[2][2]<<endl;
	    std::cout<<"dwall = "<<lf_wall<<" vs rec dwall = "<<dwall_reco<<", n50 hits = "<<lf_n50[2]<<std::endl;
    	std::cout<<"####################################"<<std::endl<<std::endl;
    }

    // Goodness of fit & dirKS
    //Simple goodness of fit and dirKS
    hgood->Fill(lf_good);
    hdirKS->Fill(lf_dirKS[2]);
	}
  
  efficiency = (nEvents==0 ? 0:nVertexFound/nEvents);
  misIDFV = (nEvents==0 ? 0:nMisIDFV/nMisIDFV_total);
  cout << "Efficiency in loop : "<< efficiency << " and misID = " << misIDFV << endl;

  double misIDExtCut;
  misIDExtCut = (nEvents==0 ? 0:nMisIDExtCut/nMisIDExtCut_total);
  cout<<"Efficiency in loop : "<< efficiency << " and misID = "<<misIDExtCut<<endl;
  cout<<"Number of events passing the FV and other cuts = "<<nEvents/nEventsTotal<<endl;

  vertex_resolution = vRes(hvertex);
  double error = vResAndError(hvertex,CredibleInterval,verbose);
  cout<<"Vertex resolution = "<<vertex_resolution<<endl;

  if(!wrt_to_ROOT){
    TCanvas * cvertex2D = new TCanvas("cvertex2D","");
    hvertex2D->Draw("colz");
    cvertex2D->SaveAs("plots/vertex2D.eps");

    TCanvas * ctruevertex2D = new TCanvas("ctruevertex2D","");
    htruevertex2D->Draw("colz");
    ctruevertex2D->SaveAs("plots/truevertex2D.eps");
    
    TCanvas * cvertex = new TCanvas("cvertex","");
    hvertex->Rebin(10);
    hvertex->GetXaxis()->SetRangeUser(0,300);
    hvertex->Scale(1./hvertex->Integral());
    hvertex->SetLineColor(kBlue);
    hvertex->Draw();
    hvertex->SaveAs("root_output/hvertex.root");
    cvertex->SaveAs("plots/vertex.eps");

    TCanvas * cvertexLongitudinal = new TCanvas("cvertexLongitudinal","");
    hvertexLongitudinal->Rebin(10);
    hvertexLongitudinal->GetXaxis()->SetRangeUser(0,300);
    hvertexLongitudinal->SetLineColor(kBlue);
    hvertexLongitudinal->Draw();
    cvertexLongitudinal->SaveAs("plots/vertexLongitudinal.eps");
    
    TCanvas * cvertexOrtho = new TCanvas("cvertexOrtho","");
    hvertexOrtho->Rebin(10);
    hvertexOrtho->GetXaxis()->SetRangeUser(0,300);
    hvertexOrtho->SetLineColor(kBlue);
    hvertexOrtho->Draw();
    cvertexOrtho->SaveAs("plots/vertexOrtho.eps");

    TCanvas * cvertexLongitudinalXTime = new TCanvas("cvertexLongitudinalXTime","");
    hvertexLongitudinalXTime->Draw("colz");

    TCanvas * cvertexXdwall = new TCanvas("cvertexXdwall","");
    hvertexXdwall->Draw("colz");

    TCanvas * cvertexXtowall = new TCanvas("cvertexXtowall","");
    hvertexXtowall->Draw("colz");   

    TCanvas * cvertexXdwall1DAll = new TCanvas("cvertexXdwall1DAll","");
    const int nSlicesdwall = nbin_distance;//16;//hvertexXdwall[0]->GetNbinsX();
    cvertexXdwall1DAll->Divide(TMath::Sqrt(nSlicesdwall),TMath::Sqrt(nSlicesdwall));
    TH1D * hvertexXdwall1D[nSlicesdwall];
    double dwallPosition[nSlicesdwall];
    double resolution68[nSlicesdwall];
    double error_resolution68[nSlicesdwall];
    double error_x[nSlicesdwall];
    TGraphErrors * gdwall;
    for(int ibinx=1;ibinx<=nSlicesdwall;ibinx++){
      hvertexXdwall1D[ibinx-1] = (TH1D*) hvertexXdwall->ProjectionX(Form("hvertexXdwall1D_%d",ibinx-1),ibinx,ibinx);
      dwallPosition[ibinx-1]=hvertexXdwall->GetYaxis()->GetBinCenter(ibinx);
      resolution68[ibinx-1]=vRes(hvertexXdwall1D[ibinx-1]);
      error_resolution68[ibinx-1]=vResAndError(hvertexXdwall1D[ibinx-1],CredibleInterval,verbose);
      //cout<<"Error = "<<error_resolution68[ibinx-1]<<endl;
      error_x[ibinx-1]=0;
      //hvertexXdwall1D[ibinx-1]->Scale(1./hvertexXdwall1D[ibinx-1]->Integral());
      hvertexXdwall1D[ibinx-1]->GetXaxis()->SetTitle(hvertexXdwall->GetYaxis()->GetTitle());
      cvertexXdwall1DAll->cd(ibinx);
      hvertexXdwall1D[ibinx-1]->SetLineColor(kBlue);
      hvertexXdwall1D[ibinx-1]->Draw();
    }
    gdwall = new TGraphErrors(nSlicesdwall,dwallPosition,resolution68,error_x,error_resolution68);
    gdwall->GetXaxis()->SetTitle(hvertexXdwall->GetYaxis()->GetTitle());
    gdwall->GetYaxis()->SetTitle(hvertexXdwall->GetXaxis()->GetTitle());
    gdwall->SaveAs(Form("plots/gdwall_%s.root",file->GetName()));

    TCanvas * cvertexXdwall1D = new TCanvas("cvertexXdwall1D","");
    gdwall->SetMarkerStyle(20);
    gdwall->GetXaxis()->SetRangeUser(min_distance,max_distance);
    gdwall->GetYaxis()->SetRangeUser(50,100);
    gdwall->SetLineColor(kBlue);
    gdwall->SetMarkerColor(kBlue);
    gdwall->Draw("ALP");
    gdwall->SaveAs(Form("plots/gdwall_%s.root",file->GetName()));
    cvertexXdwall1D->SaveAs("plots/vertexXdwall.eps");

    TCanvas * cvertexXtowall1DAll = new TCanvas("cvertexXtowall1DAll","");
    const int nSlicestowall = nbin_distance;//16;//hvertexXtowall[0]->GetNbinsX();
    cvertexXtowall1DAll->Divide(TMath::Sqrt(nSlicestowall),TMath::Sqrt(nSlicestowall));
    TH1D * hvertexXtowall1D[nSlicestowall];
    double towallPosition[nSlicestowall];
    double resolution68_towall[nSlicestowall];
    double error_resolution68_towall[nSlicestowall];
    double error_x_towall[nSlicestowall];
    TGraphErrors * gtowall;
    for(int ibinx=1;ibinx<=nSlicestowall;ibinx++){
      hvertexXtowall1D[ibinx-1] = (TH1D*) hvertexXtowall->ProjectionX(Form("hvertexXtowall1D_%d",ibinx-1),ibinx,ibinx);
      towallPosition[ibinx-1]=hvertexXtowall->GetYaxis()->GetBinCenter(ibinx);
      resolution68_towall[ibinx-1]=vRes(hvertexXtowall1D[ibinx-1]);
      error_resolution68_towall[ibinx-1]=vResAndError(hvertexXtowall1D[ibinx-1],CredibleInterval,verbose);
      //cout<<"Error = "<<error_resolution68_towall[ibinx-1]<<endl;
      error_x_towall[ibinx-1]=0;
      //hvertexXtowall1D[ibinx-1]->Scale(1./hvertexXtowall1D[ibinx-1]->Integral());
      hvertexXtowall1D[ibinx-1]->GetXaxis()->SetTitle(hvertexXtowall->GetYaxis()->GetTitle());
      cvertexXtowall1DAll->cd(ibinx);
      hvertexXtowall1D[ibinx-1]->SetLineColor(kBlue);
      hvertexXtowall1D[ibinx-1]->Draw();
      gtowall = new TGraphErrors(nSlicestowall,towallPosition,resolution68_towall,error_x_towall,error_resolution68_towall);
      gtowall->GetXaxis()->SetTitle(hvertexXtowall->GetYaxis()->GetTitle());
      gtowall->GetYaxis()->SetTitle(hvertexXtowall->GetXaxis()->GetTitle());
      gtowall->SaveAs(Form("plots/gtowall_%s.root",file->GetName()));
    }

    TCanvas * cvertexXtowall1D = new TCanvas("cvertexXtowall1D","");
    gtowall->SetMarkerStyle(20);
    gtowall->GetXaxis()->SetRangeUser(min_distance,max_distance);
    gtowall->GetYaxis()->SetRangeUser(50,100);
    gtowall->SetLineColor(kBlue);
    gtowall->SetMarkerColor(kBlue);
    gtowall->Draw("ALP");
    gtowall->SaveAs(Form("plots/gtowall_%s.root",file->GetName()));
    cvertexXtowall1D->SaveAs("plots/vertexXtowall.eps");

    TCanvas * cradius = new TCanvas("cradius","");
    hradius->Rebin(10);
    hradius->Scale(1./hradius->Integral());
    hradius->SetLineColor(kBlue);
    hradius->Draw();
    hradius->SaveAs(Form("root_output/hradius_%s.root",file->GetName()));
    cradius->SaveAs("plots/radius.eps");

    TCanvas * cdwallXdwallrec = new TCanvas("cdwallXdwallrec","");
    hdwallXdwallrec->Draw("colz");

    TCanvas * cnhits = new TCanvas("cnhits","");
    hnhits->SetLineColor(kBlue);
    hnhits->Draw();

    TCanvas * cvertexXnhits = new TCanvas("cvertexXnhits","");
    hvertexXnhits->Rebin2D(1,2);
    hvertexXnhits->Draw("colz");    

    TCanvas * cdwallXnhits = new TCanvas("cdwallXnhits","");
    hdwallXnhits->Rebin2D(1,2);
    hdwallXnhits->Draw("colz");

    TCanvas * ctowallXnhits = new TCanvas("ctowallXnhits","");
    htowallXnhits->Rebin2D(1,2);
    htowallXnhits->Draw("colz");

    TCanvas * cvertexXnhits_intime = new TCanvas("cvertexXnhits_intime","");
    hvertexXnhits_intime->Rebin2D(10,1);
    hvertexXnhits_intime->GetYaxis()->SetRangeUser(0,10);
    hvertexXnhits_intime->Draw("colz");  

    TCanvas * cvertexXgood = new TCanvas("cvertexXgood","");
    hvertexXgood->Rebin2D(100,2);
    hvertexXgood->GetYaxis()->SetRangeUser(0.,1.);
    hvertexXgood->Draw("colz");

    TCanvas * cdirection = new TCanvas("cdirection","");
    cout<<"Direction resolution = "<<vRes(hdirection)<<endl;
    //hdirection->Rebin(2);
    hdirection->GetXaxis()->SetRangeUser(0,100);
    hdirection->Scale(1./hdirection->Integral());
    hdirection->SetLineColor(kBlue);
    hdirection->Draw();
    hdirection->SaveAs(Form("plots/hdirection_%s.root",file->GetName()));
    cdirection->SaveAs("plots/direction.eps");

    TCanvas * cgood = new TCanvas("cgood","");
    //hgood->Rebin(2);
    hgood->Scale(1./hgood->Integral());
    hgood->SetLineColor(kBlue);
    hgood->Draw();
    hgood->SaveAs(Form("plots/hgood_%s.root",file->GetName()));
    cgood->SaveAs("plots/good.eps");
  
    TCanvas * cdirKS = new TCanvas("cdirKS","");
    //hdirKS->Rebin(2);
    hdirKS->Scale(1./hdirKS->Integral());
    hdirKS->SetLineColor(kBlue);
    hdirKS->Draw();
    hdirKS->SaveAs(Form("plots/hdirKS_%s.root",file->GetName()));
    cdirKS->SaveAs("plots/dirKS.eps");

    TCanvas * cvertexMisIDFV = new TCanvas("cvertexMisIDFV","");
    gStyle->SetOptStat(0.);
    cvertexMisIDFV->SetGridx(true);
    cvertexMisIDFV->SetGridy(true);
    hvertexMisIDFV->Divide(hvertexMisIDFV_total);
    hvertexMisIDFV->SetLineWidth(2);
    hvertexMisIDFV->SetLineColor(kBlue);
    hvertexMisIDFV->GetXaxis()->SetRangeUser(0,200); 
    hvertexMisIDFV->GetYaxis()->SetRangeUser(0,1.); 
    hvertexMisIDFV->Draw();
    hvertexMisIDFV->SaveAs(Form("root_output/hvertexMisIDFV_%s.root",file->GetName()));
    cvertexMisIDFV->SaveAs("plots/vertexMisIDFV.eps");
  
    TCanvas * cvertexMisIDExtCut = new TCanvas("cvertexMisIDExtCut","");
    gStyle->SetOptStat(0.);
    cvertexMisIDExtCut->SetGridx(true);
    cvertexMisIDExtCut->SetGridy(true);
    hvertexMisIDExtCut->Divide(hvertexMisIDExtCut_total);
    hvertexMisIDExtCut->SetLineWidth(2);
    hvertexMisIDExtCut->SetLineColor(kBlue);
    hvertexMisIDExtCut->GetXaxis()->SetRangeUser(0,2000); 
    hvertexMisIDExtCut->GetYaxis()->SetRangeUser(0,1.); 
    hvertexMisIDExtCut->Draw();
    hvertexMisIDExtCut->SaveAs(Form("root_output/hvertexMisIDExtCut_%s.root",file->GetName()));
    cvertexMisIDExtCut->SaveAs("plots/vertexMisIDExtCut.eps");
  }else{
    // Open the outputfile
    TFile * outfile = new TFile(ofname, "RECREATE");
    cout<<"File "<<ofname<<" is open for writing"<<endl;
    
    //Vertex related (1D)
    hvertexLongitudinal->Write();
  	hvertexOrtho->Write();
    hvertex->Write();
    hradius->Write();
    //Vertex related (2D)
    hvertexLongitudinalXTime->Draw("colz");
    hvertexLongitudinalXTime->Write();
    hvertex2D->Draw("colz");
    hvertex2D->Write();
    htruevertex2D->Draw("colz");
    htruevertex2D->Write();
    hvertexXdwall->Draw("colz");
    hvertexXdwall->Write();
    hvertexXtowall->Draw("colz");
    hvertexXtowall->Write();
    hdwallXdwallrec->Draw("colz"); 
    hdwallXdwallrec->Write();
    
    //Nhits related (1D)
    hnhits->Write();
    //Nhits related (2D)
    hdwallXnhits->Draw("colz");
    hdwallXnhits->Write();
    htowallXnhits->Draw("colz");
    htowallXnhits->Write();

    //Direction related
    hdirection->Write();
    hOpeningAngle->Write();
    hdirX->Write();
    hdirY->Write();
    hdirZ->Write();

    //Goodness of fit & dirKS related
    hgood->Write();
    hdirKS->Write();

    //Mix of previous variables
    hvertexXnhits_intime->Draw("colz");
    hvertexXnhits_intime->Write();
    hvertexXgood->Draw("colz");
    hvertexXgood->Write();
    hvertexXnhits->Draw("colz");
    hvertexXnhits->Write();

    //Cut related (1D)
    hvertexMisIDFV->Write(); 
    hvertexMisIDFV_total->Write();
    hvertexMisIDExtCut->Write(); 
    hvertexMisIDExtCut_total->Write();

    //Alper's plots (1D)
    hvtx0->Write();
    hvtx1->Write();
    hvtx2->Write();
    hvtx3->Write();
    //Alper's plots (2D)
    hvtx0Xdwall->Draw("colz");
    hvtx0Xdwall->Write();
    hvtx1Xdwall->Draw("colz");
    hvtx1Xdwall->Write();
    hvtx2Xdwall->Draw("colz");
    hvtx2Xdwall->Write();
    hvtx3Xdwall->Draw("colz");
    hvtx3Xdwall->Write();

    hvtx0Xtowall->Draw("colz");
    hvtx0Xtowall->Write();
    hvtx1Xtowall->Draw("colz");
    hvtx1Xtowall->Write();
    hvtx2Xtowall->Draw("colz");
    hvtx2Xtowall->Write();
    hvtx3Xtowall->Draw("colz");
    hvtx3Xtowall->Write();

    hvertexLongitudinalXdwall->Draw("colz");
    hvertexLongitudinalXdwall->Write();
    hvertexOrthoXdwall->Draw("colz");
    hvertexOrthoXdwall->Write();
    hradiusXdwall->Draw("colz");
    hradiusXdwall->Write();

    hvertexLongitudinalXtowall->Draw("colz");
    hvertexLongitudinalXtowall->Write();
    hvertexOrthoXtowall->Draw("colz");
    hvertexOrthoXtowall->Write();
    hradiusXtowall->Draw("colz");
    hradiusXtowall->Write();

    outfile->Close();
  }

  cout<<"End of this file"<<endl;
  return 0;
}

int main(int argc, char **argv){
  // Optional arguments for infile, outfile: adapted from Analysis code for WCSim
  char * filename=NULL;
  char * outfilename=NULL;
  double NRG = 0.;

  bool verb = false;
	bool wrt_ROOT = false;
  int c = -1;

  while( (c = getopt(argc,argv,"f:o:e:v:w")) != -1 ){
    // Input in c the argument (-f etc...) and in optarg the next argument.
    // When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
    case 'f':
      filename = optarg;
      break;
    case 'o':
      outfilename = optarg;
      break;
    case 'e':
      NRG = atoi(optarg);
      break;
    case 'v':
      verb = true;
      break;
		case 'w':
			wrt_ROOT = true;
			break;
    }
  }

  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  
  TApplication * theApp;
  if(!wrt_ROOT) theApp = new TApplication("theApp", &argc, argv);

  double vertex_resolution, efficiency, migrationProbability;

  // Implementation done only for mono energetic case
  vertexResolution(NRG, vertex_resolution, efficiency, migrationProbability, verb, wrt_ROOT, filename, outfilename);

	if(!wrt_ROOT) theApp->Run();
  return 0;
}
