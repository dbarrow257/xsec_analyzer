#include "FilePropertiesManager.hh"
#include "EventCategory.hh"

#include "TChain.h"

#include <vector>
#include <iostream>
#include <iomanip>

void PrintCuts(TChain* Chain, TString Run, TString NTupleType_Str, std::vector<TString> CutStrings) {
  float TotalEvents = Chain->GetEntries();

  std::cout<< "\n=============================================" << std::endl;
  std::cout << "Run:" << Run << std::endl;
  TString TotalCut = "";
  for ( const auto& cut : CutStrings) {
    if (TotalCut != "") {
      TotalCut += " && ";
    }

    TotalCut += cut + "==1";
    float Entries = Chain->GetEntries(TotalCut);
    std::cout << std::setw(15) << NTupleType_Str << " | " << std::setw(40) << cut << " : " << Entries << Form(" (%3.2f%)",100.0*Entries/TotalEvents) << std::endl;
  }
  std::cout << "\n\n" << std::endl;
}

int main() {
  std::cout << std::setprecision(10);
  std::vector<int> runs = {1,2,3};

  // Get access to the singleton utility class that manages the processed
  // ntuple files
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();
  
  // Get the number of defined EventCategory values by checking the number
  // of elements in the "label map" managed by the EventCategoryInterpreter
  // singleton class
  const auto& eci = EventCategoryInterpreter::Instance();
  std::map< EventCategory, std::string > category_map = eci.label_map();

  std::map< int, std::string > interaction_map = {
    {0, "QE"},
    {3, "Coh"},
    {10, "MEC"},
    {1, "RES"},
    {2, "DIS"}
  };
  
  std::vector<NtupleFileType> FileTypes = {NtupleFileType::kOnBNB,NtupleFileType::kExtBNB,NtupleFileType::kNumuMC,NtupleFileType::kIntrinsicNueMC,NtupleFileType::kDirtMC};
  std::vector<TString> CutStrings = {"1","CC1mu2p0pi_sel_nslice_eq_1","CC1mu2p0pi_sel_reco_vertex_in_FV","CC1mu2p0pi_sel_npfps_eq_3","CC1mu2p0pi_sel_ntracks_eq_3","CC1mu2p0pi_sel_correctparticles","CC1mu2p0pi_sel_containedparticles","CC1mu2p0pi_sel_momentum_threshold_passed","CC1mu2p0pi_Selected"};

  TChain* TotalNuMuOverlaySample = new TChain("stv_tree");

  const auto& ntuple_map = fpm.ntuple_file_map();
  for ( const auto& filetype : FileTypes ) {
    TChain* AllRunsChain = new TChain( "stv_tree" );
    TString NTupleType_Str = "";

    for ( int run : runs ) {
      const auto& ntuple_files = ntuple_map.at( run ).at( filetype );
      TChain* SingleRunChain = new TChain( "stv_tree" );
      
      for ( const auto& file_name : ntuple_files ) {
	NTupleType_Str = fpm.ntuple_type_to_string(fpm.get_ntuple_file_type(file_name));
	
	AllRunsChain->Add(file_name.c_str());
	SingleRunChain->Add(file_name.c_str());
      }
      
      PrintCuts(SingleRunChain,Form("%i",run),NTupleType_Str,CutStrings);
    }

    if (NTupleType_Str == "numuMC") {
      TotalNuMuOverlaySample = (TChain*)AllRunsChain->Clone();
    }
    
    PrintCuts(AllRunsChain,"All",NTupleType_Str,CutStrings);
  }
  double TotalNuMuOverlayCC2PEntries = TotalNuMuOverlaySample->GetEntries(CutStrings[CutStrings.size()-1]+"==1");

  std::cout << "\n=============================================" << std::endl;
  std::cout << "True Event Cateogry for selected events" << std::endl;

  std::cout << std::setw(3) << "-1" << " " << std::setw(30) << "Total" << " = " << std::setw(10) << TotalNuMuOverlayCC2PEntries << Form(" (%2.2f)",100.*TotalNuMuOverlayCC2PEntries/TotalNuMuOverlayCC2PEntries) << std::endl;

  for (auto const& EC : category_map) {
    TString CutStr = CutStrings[CutStrings.size()-1]+Form("==1 && CC1mu2p0pi_EventCategory==%i",EC.first);
    double Entries = TotalNuMuOverlaySample->GetEntries(CutStr);
    std::cout << std::setw(3) << EC.first << " " << std::setw(30) << EC.second << " = " << std::setw(10) << Entries << Form(" (%2.2f)",100.*Entries/TotalNuMuOverlayCC2PEntries) << std::endl;;
  }

  std::cout << "\n=============================================" << std::endl;
  std::cout << "Mode breakdown for selected events" << std::endl;

  std::cout << std::setw(3) << "-1" << " " << std::setw(30) << "Total" << " = " << std::setw(10) << TotalNuMuOverlayCC2PEntries << Form(" (%2.2f)",100.*TotalNuMuOverlayCC2PEntries/TotalNuMuOverlayCC2PEntries) << std::endl;

  //NC
  TString NC_CutStr = CutStrings[CutStrings.size()-1]+"==1 && mc_ccnc==1";
  double NC_Entries = TotalNuMuOverlaySample->GetEntries(NC_CutStr);
  std::cout << std::setw(3) << -1 << " " << std::setw(30) << "NC" << " = " << std::setw(10) << NC_Entries << Form(" (%2.2f)",100.*NC_Entries/TotalNuMuOverlayCC2PEntries) << std::endl;

  //Nue
  TString Nue_CutStr = CutStrings[CutStrings.size()-1]+"==1 && mc_ccnc==0 && mc_nu_pdg==12";
  double Nue_Entries = TotalNuMuOverlaySample->GetEntries(Nue_CutStr);
  std::cout << std::setw(3) << -1 << " " << std::setw(30) << "Nue" << " = " << std::setw(10) << Nue_Entries << Form(" (%2.2f)",100.*Nue_Entries/TotalNuMuOverlayCC2PEntries) << std::endl;
  
  for ( auto& mode : interaction_map ) {
    TString CutStr = CutStrings[CutStrings.size()-1]+Form("==1 && mc_ccnc==0 && mc_nu_pdg!=12 && mc_interaction==%i",mode.first);    
    double Entries = TotalNuMuOverlaySample->GetEntries(CutStr);
    std::cout << std::setw(3) << mode.first << " " << std::setw(30) << mode.second << " = " << std::setw(10) << Entries << Form(" (%2.2f)",100.*Entries/TotalNuMuOverlayCC2PEntries) << std::endl;
  }

  std::cout << "\n=============================================" << std::endl;
  std::cout << "Mode breakdown for signal events" << std::endl;

  double TotalNuMuOverlaySignalEntries = TotalNuMuOverlaySample->GetEntries( CutStrings[CutStrings.size()-1]+"==1 & CC1mu2p0pi_MC_Signal==1");

  std::cout << std::setw(3) << "-1" << " " << std::setw(30) << "Total" << " = " << std::setw(10) << TotalNuMuOverlaySignalEntries << Form(" (%2.2f)",100.*TotalNuMuOverlaySignalEntries/TotalNuMuOverlaySignalEntries) << std::endl;
  for ( auto& mode : interaction_map ) {
    TString CutStr = CutStrings[CutStrings.size()-1]+Form("==1 && CC1mu2p0pi_MC_Signal==1 && mc_interaction==%i",mode.first);    
    double Entries = TotalNuMuOverlaySample->GetEntries(CutStr);
    std::cout << std::setw(3) << mode.first << " " << std::setw(30) << mode.second << " = " << std::setw(10) << Entries << Form(" (%2.2f)",100.*Entries/TotalNuMuOverlaySignalEntries) << std::endl;
  }

  std::cout << "\n=============================================" << std::endl;
  std::cout << "Total Number of Signal Events:" << std::endl;

  double TotalSignalEntries = TotalNuMuOverlaySample->GetEntries("CC1mu2p0pi_MC_Signal==1");
  std::cout << std::setw(3) << "-1" << " " << std::setw(30) << "Total" << " = " << std::setw(10) << TotalSignalEntries << std::endl;

}
