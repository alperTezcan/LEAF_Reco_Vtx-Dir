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

void analyze2DHistograms(const char* outputFileName) {
  // Open the output CSV file
  ofstream outputFile(outputFileName);

  // Energy list
  int E[6] = {2, 3, 5, 10, 20, 50};

  // Loop over each ROOT file
  for (int i = 0; i<6; i++) {
  	// Open the ROOT file and get the 2D histogram
    TFile* file = new TFile(Form("../Outputs/RND/%dMeV/wcsim_output_e%d_rnd_reco_plots.root", E[i], E[i]), "READ");
    TH2* histD = (TH2*) file->Get("hvertexXdwall");
		TH2* histT = (TH2*) file->Get("hvertexXtowall");

    // Rebin the histogram along the Y-axis
    histD->RebinY(25);
		histT->RebinY(100);

    // Get the number of Y bins
    int numYBinsD = 7;
    int numYBinsT = histT->GetNbinsY();
    
    // Write the header to the CSV file
    if(i == 0){
      outputFile << "E,";
      for(int binD = 0; binD < numYBinsD; binD++){
        outputFile << "binD" << binD << "entries," << "binD" << binD << "mean," << "binD" << binD << "std,";
        /*
        outputFile << "binD" << binD << "meanY," << "binD" << binD << "stdY,";
        outputFile << "binD" << binD << "meanZ," << "binD" << binD << "stdZ,";
        outputFile << "binD" << binD << "meanT," << "binD" << binD << "stdT,";
        */
      }
      for(int binT = 0; binT < numYBinsT; binT++){
        outputFile << "binT" << binT << "entries," << "binT" << binT << "mean," << "binT" << binT << "std,";
        /*
        outputFile << "binT" << binT << "meanY," << "binT" << binT << "stdY,";
        outputFile << "binT" << binT << "meanZ," << "binT" << binT << "stdZ,";
        outputFile << "binT" << binT << "meanT," << "binT" << binT << "stdT,";
        */
      }
      outputFile << endl;
    }
    
    outputFile << E[i] << ",";

    // Loop over each Y bin for Dwall
    for (int yBinD = 1; yBinD < 8; yBinD++){
      // Get the statistics for each X-axis bin within the current Y-axis bin
      double mean = histD->ProjectionX("proj", yBinD, yBinD)->GetMean();
      double rms = histD->ProjectionX("proj", yBinD, yBinD)->GetRMS();
      double entries = histD->ProjectionX("proj", yBinD, yBinD)->GetEntries();

      // Write the statistics to the CSV file
      outputFile << entries << "," << mean << "," << rms << ",";
    }

    // Loop over each Y bin for Towall
    for (int yBinT = 1; yBinT <= numYBinsT; yBinT++){
      // Get the statistics for each X-axis bin within the current Y-axis bin
      double mean = histT->ProjectionX("proj", yBinT, yBinT)->GetMean();
      double rms = histT->ProjectionX("proj", yBinT, yBinT)->GetRMS();
      double entries = histT->ProjectionX("proj", yBinT, yBinT)->GetEntries();

      // Write the statistics to the CSV file
      outputFile << entries << "," << mean << "," << rms << ",";
    }
    
    // Write a new line to the CSV file to move to the next file
    outputFile << endl;

    // Close the ROOT file
    file->Close();
    delete file;
  }
  // Close the output CSV file
  outputFile.close();
}

int main() {
  const char* outputFileName = "../Outputs/2dBin_output.csv"; // Name of the output CSV file

  // Call the function to analyze the 2D histograms and write to CSV
  analyze2DHistograms(outputFileName);

  return 0;
}