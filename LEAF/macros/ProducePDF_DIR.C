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
#include <TSpline.h>

using namespace std;

/*
  This piece of code is a adapted from ProducePDF_new.C code.
  It creates the DIR probability density function needed for LEAF to work.
  It uses a root input file that contains necessary information extracted 
  from an initial WCSim MC simulation.
  Alper wrote it as an input-taking executable form instead of pure root version
*/

int producePDF(char *file, char *ofile, int factor, bool verbose){
  TFile* infile = new TFile(file, "read");
  
  if(verbose) cout << "Reading file: " << infile->GetName() << endl;
  
  const int nPMTtypes = 1;
  int colorScale[nPMTtypes];
  for(int i=0;i<nPMTtypes;i++){
    if(i==0) colorScale[i] = 4;
    else if(i==1) colorScale[i] = 2;
    else colorScale[i] = 3;
  }

  TH1D* ChargeProfile  [nPMTtypes];
  TH1D* HitProfile     [nPMTtypes];

  TSpline3 * splineCharge[nPMTtypes];
  TSpline3 * splineHit[nPMTtypes];
  TGraph * graphCharge[nPMTtypes];
  TGraph * graphHit[nPMTtypes];

  for(int i=0; i<nPMTtypes; i++){
    if(verbose) cout << "PMT type = " << i << endl;

    ChargeProfile[i] = (TH1D*) infile->Get(Form("ChargeProfile_pmtType%d",i));
    ChargeProfile[i]->SetLineWidth(2);
    ChargeProfile[i]->SetLineColor(colorScale[i]);
    ChargeProfile[i]->GetXaxis()->SetTitle("Angle (degrees)");
    ChargeProfile[i]->GetYaxis()->SetTitle("Charge (p.e) (normalized)");
    
    HitProfile[i] = (TH1D*) infile->Get(Form("HitProfile_pmtType%d",i));
    HitProfile[i]->SetLineWidth(2);
    HitProfile[i]->SetLineColor(colorScale[i]);
    HitProfile[i]->GetXaxis()->SetTitle("Angle (degrees)");
    HitProfile[i]->GetYaxis()->SetTitle("Hits (normalized)");
  }

  for(int i=0; i<nPMTtypes; i++){
    ChargeProfile[i]->Rebin(factor);
    HitProfile[i]->Rebin(factor);
    
    double scaleValueC = ChargeProfile[i]->Integral();
    double scaleValueH = HitProfile[i]->Integral();
    ChargeProfile[i]->Scale(1/scaleValueC);
    HitProfile[i]->Scale(1/scaleValueH);

    vector <double> xPosC, yPosC, xPosH, yPosH;
    xPosC.clear();
    yPosC.clear();
    xPosH.clear();
    yPosH.clear();

    int nBinsC_graph = ChargeProfile[i]->GetNbinsX();
    int nBinsH_graph = HitProfile[i]->GetNbinsX();

    //Loop over all bins of the histogram (0-180 degrees with bins of 0.25 * factor degrees)
    for(int ibinx = 0; ibinx<=nBinsC_graph; ibinx++){
	    xPosC.push_back(ChargeProfile[i]->GetBinCenter(ibinx));
	    yPosC.push_back(ChargeProfile[i]->GetBinContent(ibinx));
    }

    for(int ibinx = 0; ibinx<=nBinsH_graph; ibinx++){
      xPosH.push_back(HitProfile[i]->GetBinCenter(ibinx));
      yPosH.push_back(HitProfile[i]->GetBinContent(ibinx));
    }

    double * xPosC_graph = new double[nBinsC_graph];
    double * yPosC_graph = new double[nBinsC_graph];

    double * xPosH_graph = new double[nBinsH_graph];
    double * yPosH_graph = new double[nBinsH_graph];
    
    for(int ibin = 0; ibin<nBinsC_graph;ibin++){
      xPosC_graph[ibin] = xPosC.at(ibin);
      yPosC_graph[ibin] = yPosC.at(ibin);
    }

    for(int ibin = 0; ibin<nBinsH_graph;ibin++){
      xPosH_graph[ibin] = xPosH.at(ibin);
      yPosH_graph[ibin] = yPosH.at(ibin);
    }
    
    graphCharge[i] = new TGraph(nBinsC_graph,xPosC_graph,yPosC_graph);
    graphCharge[i]->SetMarkerColor(kRed);
    graphCharge[i]->SetLineColor(kRed);
    graphCharge[i]->SetMarkerStyle(20);
    
    splineCharge[i] = new TSpline3(Form("splineCharge%d",i), graphCharge[i]);
    splineCharge[i]->SetLineColor(kGreen);
    if(verbose) splineCharge[i]->Draw("lcsame");
    
    graphHit[i] = new TGraph(nBinsH_graph,xPosH_graph,yPosH_graph);
    graphHit[i]->SetMarkerColor(kRed);
    graphHit[i]->SetLineColor(kRed);
    graphHit[i]->SetMarkerStyle(20);
    
    splineHit[i] = new TSpline3(Form("splineHit%d",i), graphHit[i]);
    splineHit[i]->SetLineColor(kGreen);
    if(verbose) splineHit[i]->Draw("lcsame");
  }

  if(verbose){
    HitProfile[0]->Draw();  
    graphHit[0]->Draw("LPsame");
    splineHit[0]->Draw("Csame");
  }

  TFile * outfile = new TFile(ofile, "recreate");
  if(verbose) cout<<"Opened "<<outfile->GetName()<<" for writing."<<endl;

  for ( int i=0 ; i<nPMTtypes ; i++ ){
    splineCharge[i]->SetName(Form("SplineCharge_%d",i));
    splineHit[i]->SetName(Form("SplineHit_%d",i));
    graphCharge[i]->SetName(Form("GraphCharge_%d",i));
    graphHit[i]->SetName(Form("GraphHit_%d",i));
    
    ChargeProfile[i]->Write();
    HitProfile   [i]->Write();
    
    splineCharge[i]->Write();
    graphCharge [i]->Write();
    splineHit   [i]->Write();
    graphHit    [i]->Write();
  }
  
  outfile->Close();
  return 0;
}

int main(int argc, char **argv){
  // Optional arguments for infile, outfile: adapted from Analysis code for WCSim
  char * filename=NULL;
  char * outfilename=NULL;
  bool verb = false;
  int factor = 1;
  int c = -1;

  while( (c = getopt(argc,argv,"f:o:rb:v")) != -1 ){
    // Input in c the argument (-f etc...) and in optarg the next argument.
    // When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
    case 'f':
      filename = optarg;
      break;
    case 'o':
      outfilename = optarg;
      break;
    case 'r':
      if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
   			optarg = argv[optind++];
			}
			if (optarg == NULL) {
				std::cerr << "Error: Argument for -r option is missing." << std::endl;
				exit(EXIT_FAILURE);
			}
			factor = atoi(optarg);
      break;
    case 'v':
      verb = true;
      break;
    }
  }
  
  producePDF(filename, outfilename, factor, verb);
  
  return 0;
}
