#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

void plotHistograms(const char* file1, const char* file2, const char* file3, const char* outputFileName) {
  // Open the Root files
  TFile* f1 = new TFile(file1, "READ");
  TFile* f2 = new TFile(file2, "READ");
	TFile* f3 = new TFile(file3, "READ");

  // Get the histograms
  TH1* hC1 = (TH1*)f1->Get("ChargeProfile_pmtType0");
  TH1* hC2 = (TH1*)f2->Get("ChargeProfile_pmtType0");
  TH1* hC3 = (TH1*)f3->Get("ChargeProfile_pmtType0");
  TH1* hH1 = (TH1*)f1->Get("HitProfile_pmtType0");
  TH1* hH2 = (TH1*)f2->Get("HitProfile_pmtType0");
  TH1* hH3 = (TH1*)f3->Get("HitProfile_pmtType0");

  // Rebin them to 1deg bins
  hC1->Rebin(4);
  hC2->Rebin(4);
  hC3->Rebin(4);
  hH1->Rebin(4);
  hH2->Rebin(4);
  hH3->Rebin(4);

  // Normalize the histograms
  hC1->Scale(1.0 / hC1->Integral());
  hC2->Scale(1.0 / hC2->Integral());
  hC3->Scale(1.0 / hC3->Integral());
  hH1->Scale(1.0 / hH1->Integral());
  hH2->Scale(1.0 / hH2->Integral());
  hH3->Scale(1.0 / hH3->Integral());

  // Plot the histograms
  TCanvas* Ccanvas = new TCanvas("Ccanvas", "Normalized Charge Profile Histograms", 800, 600);
  hC1->SetLineColor(kRed);
  hC2->SetLineColor(kBlue);
  hC3->SetLineColor(kGreen);
  hC1->Draw();
  hC2->Draw("same");
  hC3->Draw("same");

  TCanvas* Hcanvas = new TCanvas("Hcanvas", "Normalized Hit Profile Histograms", 800, 600);
  hH1->SetLineColor(kRed);
  hH2->SetLineColor(kBlue);
  hH3->SetLineColor(kGreen);
  hH1->Draw();
  hH2->Draw("same");
  hH3->Draw("same");

  // Create histograms for relative changes
  TH1* Crel_change1 = new TH1F("Crel_change1", "(hC1 - hC2) / hC2", hC2->GetNbinsX(), hC2->GetXaxis()->GetXmin(), hC2->GetXaxis()->GetXmax());
  TH1* Crel_change2 = new TH1F("Crel_change2", "(hC3 - hC2) / hC2", hC2->GetNbinsX(), hC2->GetXaxis()->GetXmin(), hC2->GetXaxis()->GetXmax());

	TH1* Hrel_change1 = new TH1F("Hrel_change1", "(hH1 - hH2) / hH2", hH2->GetNbinsX(), hH2->GetXaxis()->GetXmin(), hH2->GetXaxis()->GetXmax());
  TH1* Hrel_change2 = new TH1F("Hrel_change2", "(hH3 - hH2) / hH2", hH2->GetNbinsX(), hH2->GetXaxis()->GetXmin(), hH2->GetXaxis()->GetXmax());

  // Calculate relative changes
  for (int i = 1; i <= hC1->GetNbinsX(); ++i) {
    double content1 = hC1->GetBinContent(i);
    double content2 = hC2->GetBinContent(i);
    double content3 = hC3->GetBinContent(i);

    Crel_change1->SetBinContent(i, (content1 - content2) / content2);
    Crel_change2->SetBinContent(i, (content3 - content2) / content2);
  }

	for (int i = 1; i <= hH1->GetNbinsX(); ++i) {
    double content1 = hH1->GetBinContent(i);
    double content2 = hH2->GetBinContent(i);
    double content3 = hH3->GetBinContent(i);

    Hrel_change1->SetBinContent(i, (content1 - content2) / content2);
    Hrel_change2->SetBinContent(i, (content3 - content2) / content2);
  }

  // Plot the relative changes
  TCanvas* Ccanvas2 = new TCanvas("Ccanvas_rel", "Relative Changes for Charge Profile", 800, 600);
  Crel_change1->SetLineColor(kRed);
  Crel_change2->SetLineColor(kBlue);
  Crel_change1->Draw();
  Crel_change2->Draw("same");

	TCanvas* Hcanvas2 = new TCanvas("Hcanvas_rel", "Relative Changes for Hit Profile", 800, 600);
  Hrel_change1->SetLineColor(kRed);
  Hrel_change2->SetLineColor(kBlue);
  Hrel_change1->Draw();
  Hrel_change2->Draw("same");

  // Create a new ROOT file to store the histograms and canvases
  TFile* outputFile = new TFile(outputFileName, "RECREATE");

  // Write histograms and canvases to the output ROOT file
  hC1->Write();
  hC2->Write();
  hC3->Write();
  Crel_change1->Write();
  Crel_change2->Write();
  Ccanvas->Write();
  Ccanvas2->Write();

	hH1->Write();
  hH2->Write();
  hH3->Write();
  Hrel_change1->Write();
  Hrel_change2->Write();
  Hcanvas->Write();
  Hcanvas2->Write();

  // Close the output ROOT file
  outputFile->Close();

  // Cleanup: delete all allocated objects
  delete f1;
  delete f2;
  delete f3;
  delete Ccanvas;
  delete Ccanvas2;
	delete Hcanvas;
  delete Hcanvas2;
  delete outputFile;
}

int main() {
  const char* file1 = "../Outputs/RND/2MeV/wcsim_output_e2_rnd_plots.root";
  const char* file2 = "../Outputs/RND/10MeV/wcsim_output_e10_rnd_plots.root";
  const char* file3 = "../Outputs/RND/50MeV/wcsim_output_e50_rnd_plots.root";
  const char* outputFileName = "../Outputs/RND/compare_profiles_high.root";

  plotHistograms(file1, file2, file3, outputFileName);

  return 0;
}