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
  This piece of code is a simplified and commented version of Benjamin's original
  ProducePDF.cpp code.
  It is meant to produce the t-tof probability density function needed for LEAF to work.
  It uses a root input file that contains necessary information extracted from an initial
  WCSim MC simulation
  Alper wrote it as an input-taking executable form instead of pure root version
*/

int producePDF(char *file, char *ofile, bool verbose){
  TFile* infile = new TFile(file, "read");
  
  if(verbose) cout << "Reading file: " << infile->GetName() << endl;
  
  const int nPMTtypes = 1;
  double DRTotalPerNS[nPMTtypes];
  int colorScale[nPMTtypes];
  for(int i=0;i<nPMTtypes;i++){
    if(i==0) colorScale[i] = 4;
    else if(i==1) colorScale[i] = 2;
    else colorScale[i] = 3;
  }

  TH1D* TimeProfile      [nPMTtypes];
  TH1D* TimeTOFProfile   [nPMTtypes];
  TH1D* HitTimeTOFProfile[nPMTtypes];
  TH1D* HitTimeTOFDR     [nPMTtypes];

  double pdfLowerLimit[nPMTtypes];
  for (int i=0 ; i<nPMTtypes ; i++ ) {
    (i==0) ? pdfLowerLimit[i] = -800 : pdfLowerLimit[i] = -800;
  }

  TSpline3* splineExpoQueue[nPMTtypes];
  TSpline3* splineDR       [nPMTtypes];

  TGraph * graphExpoQueue[nPMTtypes];
  TGraph * graphDR       [nPMTtypes];

  for(int i=0; i<nPMTtypes; i++){
    if(verbose) cout << "PMT type = " << i << endl;

    TimeProfile[i] = (TH1D*) infile->Get(Form("TimeProfile_pmtType%d",i));
    TimeProfile[i]->SetLineWidth(2);
    TimeProfile[i]->SetLineColor(colorScale[i]);
    TimeProfile[i]->GetXaxis()->SetTitle("Time (ns)");
    TimeProfile[i]->GetYaxis()->SetTitle("Charge (p.e) (normalized)");
    
    TimeTOFProfile[i] = (TH1D*) infile->Get(Form("TimeTOFProfile_pmtType%d",i));
    TimeTOFProfile[i]->SetLineWidth(2);
    TimeTOFProfile[i]->SetLineColor(colorScale[i]);
    TimeTOFProfile[i]->GetXaxis()->SetTitle("Time - TOF (ns)");
    TimeTOFProfile[i]->GetYaxis()->SetTitle("Charge (p.e) (normalized)");
    
    HitTimeTOFProfile[i] = (TH1D*) infile->Get(Form("HitTimeTOFProfile_pmtType%d",i));
    HitTimeTOFProfile[i]->SetLineWidth(2);
    HitTimeTOFProfile[i]->SetLineColor(colorScale[i]);
    HitTimeTOFProfile[i]->GetXaxis()->SetTitle("HitTime - TOF (ns)");
    HitTimeTOFProfile[i]->GetYaxis()->SetTitle("nhits (normalized)");
    
    HitTimeTOFDR[i] = (TH1D*) HitTimeTOFProfile[i]->Clone(Form("HitTimeTOFDR_pmtType%d",i));
    HitTimeTOFDR[i]->Reset();
  }

  for(int i=0; i<nPMTtypes; i++){
    double scaleValue = TimeTOFProfile[i]->Integral();
    TimeTOFProfile[i]->Scale(1/scaleValue);
    
    double startDR = -100;  // ns
    double endDR   =  -20;  // ns
    double timeWindowDR = endDR - startDR ; // ns
    
    DRTotalPerNS[i] = HitTimeTOFProfile[i]->Integral( HitTimeTOFProfile[i]->FindBin(startDR), HitTimeTOFProfile[i]->FindBin(endDR));
    DRTotalPerNS[i] /= timeWindowDR;
    
    for(int ibinx=1; ibinx <= HitTimeTOFProfile[i]->GetNbinsX(); ibinx++){
      double timeWindow = HitTimeTOFProfile[i]->GetBinWidth(ibinx);
      HitTimeTOFDR[i]->SetBinContent(ibinx,DRTotalPerNS[i]*timeWindow);
    }
    
    HitTimeTOFProfile[i]->Scale(1/HitTimeTOFProfile[i]->Integral());
    HitTimeTOFDR     [i]->Scale(1/scaleValue);
    
    if(verbose){
      cout<<"Integrated hit from -100 to 500ns = "<< HitTimeTOFProfile[i]->Integral(HitTimeTOFProfile[i]->FindBin(-100),HitTimeTOFProfile[i]->FindBin(500))<< ", Dark Rate = "<< HitTimeTOFDR[i]->Integral(HitTimeTOFDR[i]->FindBin(-100),HitTimeTOFDR[i]->FindBin(500))<<endl;
      cout << "Start to merge bins in higher size bins to optimize PDF" << endl;
    }

    vector <double> xPos, yPos, drPos;
    xPos.clear();
    yPos.clear();
    drPos.clear();
    double xAverage=0;
    double yAverage=0;
    double drAverage=0;
    const int nLimits = 8 ;
    double minValue = pdfLowerLimit[i];
    double maxValue = HitTimeTOFProfile[i]->GetXaxis()->GetXmax();
    double lowLimits[nLimits] = { minValue, -100., -30., -10., 8., 30., 100., maxValue };
    double binLowLimits[nLimits];
    for(int il = 0;il<nLimits;il++){
      binLowLimits[il] = HitTimeTOFProfile[i]->FindBin(lowLimits[il]);
    }
    int rebinLowLimits[nLimits] = {100, 10, 5, 2, 10, 30, 100, 999};
    int currentLimit = 0;
    int nBinsAverage = 0;
    int nBins = 0;//150;//Number of big bins

    //Loop over all bins of the histogram starting at pdfLowerLimit
    for(int ibinx = HitTimeTOFProfile[i]->FindBin(pdfLowerLimit[i]); ibinx<=HitTimeTOFProfile[i]->GetNbinsX(); ibinx++){
      if(verbose){
        cout<<"Bin = "<<ibinx<<" i.e. value of low edge = "<< HitTimeTOFProfile[i]->GetBinLowEdge(ibinx)<<", current low edge limit bin = "<<binLowLimits[currentLimit]<<endl;
        cout<<"nBinsAverage = "<< nBinsAverage<< ", rebinning factor = "<<rebinLowLimits[currentLimit]<<endl;
      }
      
      if( (ibinx >= binLowLimits[currentLimit] && ibinx < binLowLimits[currentLimit+1]) && nBinsAverage < rebinLowLimits[currentLimit]){
	      //Sum over bin content to average them until we reach the last bin 
	      xAverage+=HitTimeTOFProfile[i]->GetBinCenter(ibinx);
	      yAverage+=HitTimeTOFProfile[i]->GetBinContent(ibinx);
	      drAverage+=HitTimeTOFDR[i]->GetBinContent(ibinx);
	      nBinsAverage++;
      }

      if( (ibinx+1 >= binLowLimits[currentLimit+1]) || (nBinsAverage >= rebinLowLimits[currentLimit])){
        // If we reached the number of binning to rebin: store information.
	      // Same if we reached the limit of the region to rebin of a given factor.
	      // What if both happens at the same time, or worse: we reach limit of bins on bin n, and at n+1,
	      // we overcome the limit. In that case, we do not have anyting to fill our average? 
	      // So, we understand that if we are in the last bin below the limit, we should store and then pass
	      // to the next step of limit. So, the check of the limit should always be on the next bin. 
	      xPos.push_back(xAverage/nBinsAverage);
	      yPos.push_back((yAverage/nBinsAverage));
	      drPos.push_back(drAverage/nBinsAverage);
	      if(verbose) cout<<"Position = "<<xPos[nBins]<<", value="<<yPos[nBins]<<", DR = "<<drPos[nBins]<<endl;
      	nBins++;
	      nBinsAverage = 0;
	      xAverage = 0;
	      yAverage = 0;
	      drAverage = 0;
	      if(ibinx+1 >= binLowLimits[currentLimit+1]) currentLimit++;
      }
    }

    double nBins_graph = nBins;
    double * xPos_graph = new double[nBins];
    double * yPos_graph = new double[nBins];
    double * drPos_graph = new double[nBins];
    for(int ibin = 0; ibin<nBins;ibin++){
      xPos_graph[ibin] = xPos.at(ibin);
      yPos_graph[ibin] = yPos.at(ibin);
      drPos_graph[ibin] = drPos.at(ibin);
    }
    graphExpoQueue[i] = new TGraph(nBins_graph,xPos_graph,yPos_graph);
    graphExpoQueue[i]->SetMarkerColor(kRed);
    graphExpoQueue[i]->SetLineColor(kRed);
    graphExpoQueue[i]->SetMarkerStyle(20);
    
    splineExpoQueue[i] = new TSpline3(Form("splineExpoQueue%d",i),graphExpoQueue[i]);
    splineExpoQueue[i]->SetLineColor(kGreen);
    if(verbose) splineExpoQueue[i]->Draw("lcsame");
    
    graphDR[i]  = new TGraph(nBins_graph,xPos_graph,drPos_graph);
    splineDR[i] = new TSpline3(Form("splineDR%d",i),graphDR[i]);
    splineDR[i]->SetLineColor(kCyan);
  }

  if(verbose){
    HitTimeTOFProfile[0]->Draw();  
    graphExpoQueue[0]->Draw("LPsame");
    splineExpoQueue[0]->Draw("Csame");
  }

  TFile *outfile = new TFile(ofile, "recreate");
  if(verbose) cout<<"Opened "<<outfile->GetName()<<" for writing."<<endl;
  for ( int i=0 ; i<nPMTtypes ; i++ ){
    splineExpoQueue[i]->SetName(Form("SplineExpoQueue_%d",i));
    splineDR[i]->SetName(Form("SplineDR_%d",i));
    graphExpoQueue[i]->SetName(Form("GraphExpoQueue_%d",i));
    graphDR[i]->SetName(Form("GraphDR_%d",i));
    
    TimeProfile       [i]->Write();
    TimeTOFProfile    [i]->Write();
    HitTimeTOFProfile [i]->Write();
    HitTimeTOFDR      [i]->Write();
    
    // To write graph and splines
    splineExpoQueue   [i]->Write();
    splineDR          [i]->Write();
    graphExpoQueue    [i]->Write();
    graphDR           [i]->Write();
  }
  
  outfile->Close();
  return 0;
}

int main(int argc, char **argv){
  // Optional arguments for infile, outfile: adapted from Analysis code for WCSim
  char * filename=NULL;
  char * outfilename=NULL;
  bool verb = false;
  int c = -1;

  while( (c = getopt(argc,argv,"f:o:v")) != -1 ){
    // Input in c the argument (-f etc...) and in optarg the next argument.
    // When the above test becomes -1, it means it fails to find a new argument.
    switch(c){
    case 'f':
      filename = optarg;
      break;
    case 'o':
      outfilename = optarg;
      break;
    case 'v':
      verb = true;
      break;
    }
  }
  
  producePDF(filename, outfilename, verb);
  
  return 0;
}
