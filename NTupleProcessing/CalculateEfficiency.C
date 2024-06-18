#include "FilePropertiesManager.hh"
#include "EventCategory.hh"

#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

#include <vector>
#include <iostream>
#include <iomanip>

// For a given analysis TTree (containing the events to use), compute the
// efficiency and purity given definitions for what events count as signal
// (based on MC truth information) and what events are selected (based on
// reconstructed information) in terms of our analysis TTree branch variables.
// Store the results in the output variables eff and pur.
void compute_eff_pur( TTree* stv_tree, const std::string signal_cuts,
		      const std::string selection_cuts, double& eff, double& pur )
{
  // For computing efficiencies and purities, we need to only use MC events.
  // Unconditionally add this requirement to the cuts defined above.
  std::string signal = signal_cuts + " && is_mc ";
  std::string selection = selection_cuts + " && is_mc ";

  // These are actually integer counts, but since we will use them below to
  // compute ratios, intrinsically cast them to double-precision values for
  // convenience.
  double num_signal = stv_tree->GetEntries(signal.c_str());
  double num_selected = stv_tree->GetEntries(selection.c_str());
  double num_selected_signal = stv_tree->GetEntries((signal + " && " + selection).c_str());

  eff = num_selected_signal / num_signal;
  pur = num_selected_signal / num_selected;

  std::cout << "selection = " << selection << '\n';
  std::cout << "signal = " << num_signal << '\n';
  std::cout << "selected = " << num_selected << '\n';
  std::cout << "selected_signal = " << num_selected_signal << '\n';

}

int Calculate(std::vector<int> runs) {

  // Get access to the singleton utility class that manages the processed
  // ntuple files
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();
  
  std::vector<NtupleFileType> FileTypes = {NtupleFileType::kExtBNB,NtupleFileType::kNumuMC,NtupleFileType::kIntrinsicNueMC,NtupleFileType::kDirtMC};
  //std::vector<NtupleFileType> FileTypes = {NtupleFileType::kNumuMC};
  std::string SignalString = "CC1mu2p0pi_MC_Signal";
  std::vector<std::string> CutStrings = {"1","CC1mu2p0pi_sel_reco_vertex_in_FV","CC1mu2p0pi_sel_npfps_eq_3","CC1mu2p0pi_sel_ntracks_eq_3","CC1mu2p0pi_sel_correctparticles","CC1mu2p0pi_sel_containedparticles","CC1mu2p0pi_sel_momentum_threshold_passed"};

  TChain* Chain = new TChain("stv_tree");

  const auto& ntuple_map = fpm.ntuple_file_map();
  for ( const auto& filetype : FileTypes ) {

    for ( int run : runs ) {
      const auto& ntuple_files = ntuple_map.at( run ).at( filetype );
      
      for ( const auto& file_name : ntuple_files ) {
	TString NTupleType_Str = fpm.ntuple_type_to_string(fpm.get_ntuple_file_type(file_name));
	std::cout << "Adding: " << file_name << "(" << NTupleType_Str << ")" << std::endl;
	Chain->Add(file_name.c_str());
      }

    }
  }

  size_t num_points = CutStrings.size();
  TGraph* eff_graph = new TGraph( num_points );
  TGraph* pur_graph = new TGraph( num_points );

  double eff, pur;
  std::string CutString = "";

  for (int iCut=0;iCut<num_points;iCut++) {
    CutString += CutStrings[iCut];
    compute_eff_pur(Chain, SignalString, CutString, eff, pur );
    CutString += " && ";

    eff_graph->SetPoint( iCut, iCut + 1, eff );
    pur_graph->SetPoint( iCut, iCut + 1, pur );


    std::cout << "selection = " << CutString << '\n';
    std::cout << "eff = " << eff << '\n';
    std::cout << "pur = " << pur << '\n';
    std::cout << "\n\n";

  }

  const std::vector< std::string > bin_labels = { "no cuts", "in FV", "3 PFPs", "3 Tracks", "Correct Particles", "Containment", "Momentum thresh"};
  for ( int b = 1; b <= bin_labels.size(); ++b ) {
    eff_graph->GetHistogram()->GetXaxis()->SetBinLabel( eff_graph->GetHistogram()->FindBin(b),bin_labels.at(b - 1).c_str() );
  }

  TCanvas* c1 = new TCanvas;
  c1->SetBottomMargin(0.21);
  eff_graph->SetTitle( ";; efficiency or purity" );
  eff_graph->SetLineColor(kBlue);
  eff_graph->SetMarkerColor(kBlue);
  eff_graph->SetLineWidth(3);
  eff_graph->SetMarkerStyle(20);
  eff_graph->GetYaxis()->SetRangeUser( 0., 1. );
  eff_graph->Draw( "alp" );
  
  pur_graph->SetLineColor(kRed);
  pur_graph->SetMarkerColor(kRed);
  pur_graph->SetLineWidth(3);
  pur_graph->SetMarkerStyle(20);
  pur_graph->Draw("same lp");

  TLegend* lg = new TLegend(0.7, 0.7, 0.9, 0.9);
  lg->AddEntry( eff_graph, "efficiency", "lp" );
  lg->AddEntry( pur_graph, "purity", "lp" );
  lg->Draw("same");

  TString OutputName = "Performance";
  for (int run : runs) {
    OutputName += TString("_run")+TString(Form("%i",run));
  }
  c1->Print(OutputName+".png");
}

int main() {
  std::vector<int> runs;

  runs = {1};
  Calculate(runs);

  runs = {2};
  Calculate(runs);

  runs = {3};
  Calculate(runs);

  runs = {1,2,3};
  Calculate(runs);
}
