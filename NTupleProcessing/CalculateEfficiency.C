#include "FilePropertiesManager.hh"
#include "HistUtils.hh"

#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TLine.h"
#include "TStyle.h"

#include <vector>
#include <iostream>
#include <iomanip>

struct Variable {
  std::string BranchName_;
  std::string Name_;
  std::string DrawCut_;
  std::string ShortDrawCut_;
  TH1F Binning_;
  double LowBoundSelCut_ = -999999.;
  double UppBoundSelCut_ = 999999.;
};

int CalculateEfficiency() {
  gStyle->SetOptStat(false);
  gStyle->SetEndErrorSize(5);

  //================================================================================
  //Define event weights, cuts and ntuple file types considered

  std::vector<int> runs = {1,2,3};

  const std::string mc_event_weight = DEFAULT_MC_EVENT_WEIGHT;

  std::vector<NtupleFileType> FileTypes = {NtupleFileType::kExtBNB,NtupleFileType::kNumuMC,NtupleFileType::kIntrinsicNueMC,NtupleFileType::kDirtMC};

  //================================================================================
  //Variables we want to plot:

  std::vector<Variable> Variables;

  std::vector<double> MuonMomentumBinEdges = {0.1,0.3,0.5,0.7,0.9,1.2};
  Variables.emplace_back(Variable{"CC1mu2p0pi_MuonMomentumVector_Reco.Mag()","MuonMomentum","1==1","All",TH1F("CC1mu2p0piMuonMomentumVectorReco",";CC1mu2p0piMuonMomentumVectorReco.Mag();events",MuonMomentumBinEdges.size()-1,MuonMomentumBinEdges.data())});

  std::vector<double> ProtonMomentumBinEdges = {0.3,0.4,0.5,0.6,0.7,0.8,1.0};
  Variables.emplace_back(Variable{"CC1mu2p0pi_LeadingProtonMomentumVector_Reco.Mag()","LeadingProtonMomentum","1==1","All",TH1F("CC1mu2p0piLeadingProtonMomentumVectorReco",";CC1mu2p0piLeadingProtonMomentumVectorReco.Mag();events",ProtonMomentumBinEdges.size()-1,ProtonMomentumBinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_RecoilProtonMomentumVector_Reco.Mag()","RecoilProtonMomentum","1==1","All",TH1F("CC1mu2p0piRecoilProtonMomentumVectorReco",";CC1mu2p0piRecoilProtonMomentumVectorReco.Mag();events",ProtonMomentumBinEdges.size()-1,ProtonMomentumBinEdges.data())});

  std::vector<double> CosThetaBinEdges = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0};
  Variables.emplace_back(Variable{"TMath::Cos(CC1mu2p0pi_MuonMomentumVector_Reco.Theta())","MuonCosTheta","1==1","All",TH1F("CC1mu2p0piMuonMomentumVectorReco_CosTheta",";CosTheta;events",CosThetaBinEdges.size()-1,CosThetaBinEdges.data())});
  Variables.emplace_back(Variable{"TMath::Cos(CC1mu2p0pi_LeadingProtonMomentumVector_Reco.Theta())","LeadingProtonCosTheta","1==1","All",TH1F("CC1mu2p0piLeadingProtonMomentumVectorReco_CosTheta",";CosTheta;events",CosThetaBinEdges.size()-1,CosThetaBinEdges.data())});
  Variables.emplace_back(Variable{"TMath::Cos(CC1mu2p0pi_RecoilProtonMomentumVector_Reco.Theta())","RecoilProtonCosTheta","1==1","All",TH1F("CC1mu2p0piRecoilProtonMomentumVectorReco_CosTheta",";CosTheta;events",CosThetaBinEdges.size()-1,CosThetaBinEdges.data())});

  std::vector<double> PhiBinEdges = {-3.14,-2.51,-1.88,-1.26,-0.63,0.00,0.63,1.26,1.88,2.51,3.14};
  Variables.emplace_back(Variable{"CC1mu2p0pi_MuonMomentumVector_Reco.Phi()","MuonPhi","1==1","All",TH1F("CC1mu2p0piMuonMomentumVectorReco_Phi",";Phi;events",PhiBinEdges.size()-1,PhiBinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_LeadingProtonMomentumVector_Reco.Phi()","LeadingProtonPhi","1==1","All",TH1F("CC1mu2p0piLeadingProtonMomentumVectorReco_Phi",";Phi;events",PhiBinEdges.size()-1,PhiBinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_RecoilProtonMomentumVector_Reco.Phi()","RecoilProtonPhi","1==1","All",TH1F("CC1mu2p0piRecoilProtonMomentumVectorReco_Phi",";Phi;events",PhiBinEdges.size()-1,PhiBinEdges.data())});

  std::vector<double> OpeningAngleBinEdges = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0};
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_CosPlPr","CosPlPr","1==1","All",TH1F("CosPlPr",";CosPlPr;events",OpeningAngleBinEdges.size()-1,OpeningAngleBinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_CosMuPsum","CosMuPsum","1==1","All",TH1F("CosMuPsum",";CosMuPsum;events",OpeningAngleBinEdges.size()-1,OpeningAngleBinEdges.data())});

  std::vector<double> DeltaPTBinEdges = {0.,0.1,0.2,0.3,0.4,0.5,0.6,1.0};
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_Pt","DeltaPT","1==1","All",TH1F("DeltaPT",";DeltaPT;events",DeltaPTBinEdges.size()-1,DeltaPTBinEdges.data())});

  std::vector<double> OtherTKIBinEdges = {0.,30.,60.,90.,120.,150.,180.};
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_DeltaAlphaT","DeltaAlphaT","1==1","All",TH1F("DeltaAlphaT",";DeltaAlphaT;events",OtherTKIBinEdges.size()-1,OtherTKIBinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_DeltaPhiT","DeltaPhiT","1==1","All",TH1F("DeltaPhiT",";DeltaPhiT;events",OtherTKIBinEdges.size()-1,OtherTKIBinEdges.data())});

  //================================================================================

  int nCats = 3;
  int Selected = 0;
  int Signal = 1;
  int SignalSelected = 2;

  std::vector<std::string> Category(nCats);
  std::vector<std::string> CategoryCutString(nCats);

  Category[Selected] = "Selected";
  CategoryCutString[Selected] = "CC1mu2p0pi_Selected==1";

  Category[Signal] = "Signal";
  CategoryCutString[Signal] = "CC1mu2p0pi_MC_Signal==1";

  Category[SignalSelected] = "Selected Signal";
  CategoryCutString[SignalSelected] = "(CC1mu2p0pi_Selected==1 && CC1mu2p0pi_MC_Signal==1)";

  //HistogramsToFill[NTupleType][PlotVariable]
  std::vector< std::vector< std::vector<TH1F*> > > HistogramsToFill;
  HistogramsToFill.resize(FileTypes.size());
  for (size_t iFT=0;iFT<FileTypes.size();iFT++) {
    HistogramsToFill[iFT].resize(Variables.size());
    for (size_t iVar=0;iVar<Variables.size();iVar++) {
      HistogramsToFill[iFT][iVar].resize(Category.size());
      for (int iCat=0;iCat<Category.size();iCat++) {
	HistogramsToFill[iFT][iVar][iCat] = (TH1F*)(Variables[iVar].Binning_).Clone();
	HistogramsToFill[iFT][iVar][iCat]->Sumw2();
      }
    }
  }

  // Get access to the singleton utility class that manages the processed
  // ntuple files
  const FilePropertiesManager& fpm = FilePropertiesManager::Instance();
  const auto& ntuple_map = fpm.ntuple_file_map();
  std::map< std::string, FilePropertiesManager::TriggersAndPOT > data_norm_map = fpm.data_norm_map();

  int histcounter = 0;

  for ( int run : runs ) {
    
    std::string data_file_name = "";
    for ( const auto& file_name : ntuple_map.at( run ).at( NtupleFileType::kOnBNB ) ) {
      if (data_file_name != "") {
	std::cerr << "Two data files found for run :" << run << std::endl;
	throw;
      }
      data_file_name = file_name;
    }
    double DataNTrigs = (data_norm_map.at(data_file_name)).trigger_count_;
    double DataPOT = (data_norm_map.at(data_file_name)).pot_;
			
    for (int iFT=0;iFT<FileTypes.size();iFT++) {
      const auto& ntuple_files = ntuple_map.at( run ).at( FileTypes[iFT] );
				       
      for ( const auto& file_name : ntuple_files ) {
        NtupleFileType NTupleType = fpm.get_ntuple_file_type(file_name);
	double ScaleFactor = -1;

	if (NTupleType == NtupleFileType::kOnBNB) {
	  ScaleFactor = 1.;
	}
	if (NTupleType == NtupleFileType::kExtBNB) {
	  FilePropertiesManager::TriggersAndPOT TrigPot = data_norm_map.at(file_name);
	  ScaleFactor = (double)DataNTrigs/(double)TrigPot.trigger_count_;
	}
	if (NTupleType == NtupleFileType::kNumuMC || NTupleType == NtupleFileType::kIntrinsicNueMC || NTupleType == NtupleFileType::kDirtMC) {
	  TFile temp_file( file_name.c_str(), "read" );
	  TParameter<float>* temp_pot = nullptr;
	  temp_file.GetObject( "summed_pot", temp_pot );
	  double pot = temp_pot->GetVal();

	  ScaleFactor = DataPOT/pot;
	}

	std::cout << "Adding: " << std::setw(140) << file_name << " (" << std::setw(15) << fpm.ntuple_type_to_string(NTupleType) << ")" << " | Factor: " << std::setw(10) << ScaleFactor << std::endl;
	TChain Chain("stv_tree");
	Chain.Add(file_name.c_str());

	for (int iVar=0;iVar<Variables.size();iVar++) {

	  TString CutString;
	  if (NTupleType == NtupleFileType::kOnBNB) {
	    CutString = "(" + Variables[iVar].DrawCut_ + ")";
	  } else if (NTupleType == NtupleFileType::kExtBNB) {
	    CutString = TString(Form("%2.6f * ",ScaleFactor)) + "(" + Variables[iVar].DrawCut_ + ")";
	  } else {
	    CutString = Form("%2.6f * ",ScaleFactor) + mc_event_weight + " * (" + Variables[iVar].DrawCut_ + ")";
	  }

	  for (int iCat=0;iCat<Category.size();iCat++) {
	    TString HistName = Form("%i",histcounter);
	    TH1F* Hist = (TH1F*)(Variables[iVar].Binning_).Clone(HistName);
	    Hist->Sumw2();
	    histcounter += 1;

	    TString FullCutString  = CutString + " * (" + CategoryCutString[iCat] + ")";
	    
	    Chain.Draw((Variables[iVar].BranchName_+">>"+HistName),FullCutString,"goff");
	    HistogramsToFill[iFT][iVar][iCat]->Add(Hist);	    

	  }
	}

      }

    }
  }

  std::cout << "============================================================" << std::endl;

  for (size_t iFT=0;iFT<FileTypes.size();iFT++) {
    for (size_t iVar=0;iVar<Variables.size();iVar++) {
      for (size_t iCat=0;iCat<Category.size();iCat++) {
	std::cout << std::setw(15) << fpm.ntuple_type_to_string(FileTypes[iFT]) << " " << std::setw(65) << Variables[iVar].BranchName_ << " " << std::setw(35) << Category[iCat] << " | " << std::setw(15) << HistogramsToFill[iFT][iVar][iCat]->Integral() << " " << std::setw(15) << HistogramsToFill[iFT][iVar][iCat]->GetEntries() << std::endl;
      }
    }
  }

  //========================================================================================================
  //Integrate over the different file types

  std::vector< std::vector<TH1F*> > IntegratedHistograms;
  IntegratedHistograms.resize(Variables.size());
  for (size_t iVar=0;iVar<Variables.size();iVar++) {
    IntegratedHistograms[iVar].resize(Category.size());
    
    for (int iCat=0;iCat<Category.size();iCat++) {
      IntegratedHistograms[iVar][iCat] = (TH1F*)(Variables[iVar].Binning_).Clone();
      IntegratedHistograms[iVar][iCat]->Sumw2();

      for (size_t iFT=0;iFT<FileTypes.size();iFT++) {
	IntegratedHistograms[iVar][iCat]->Add(HistogramsToFill[iFT][iVar][iCat]);
      }
    }
  }

  //========================================================================================================
  //Build the stats histograms

  int nStats = 2;
  int PurityIndex = 0;
  int EfficiencyIndex = 1;

  std::vector< std::vector<TH1F*> > StatsHistograms;
  StatsHistograms.resize(Variables.size());
  for (size_t iVar=0;iVar<Variables.size();iVar++) {
    StatsHistograms[iVar].resize(nStats);
    for (int iStat=0;iStat<nStats;iStat++) {
      StatsHistograms[iVar][iStat] = new TH1F(Variables[iVar].Binning_);
    }
  }

  //========================================================================================================
  //Fill the stats histograms

  for (size_t iVar=0;iVar<Variables.size();iVar++) {
    for (int xBin=1;xBin<=Variables[iVar].Binning_.GetNbinsX();xBin++) {
      double nSelected = IntegratedHistograms[iVar][Selected]->GetBinContent(xBin);
      double nSignal = IntegratedHistograms[iVar][Signal]->GetBinContent(xBin);
      double nSignalSelected = IntegratedHistograms[iVar][SignalSelected]->GetBinContent(xBin);

      double Purity = 0.;
      double Purity_Err = 0.;
      if (nSignalSelected > 0. && nSignal > 0.) {
	Purity = nSignalSelected/nSignal;
	Purity_Err = Purity*std::sqrt((1./nSignalSelected)+(1./nSignal));
      }

      double Efficiency = 0.;
      double Efficiency_Err = 0.;
      if (nSignalSelected > 0. && nSelected > 0.) {
	Efficiency = nSignalSelected/nSelected;
	Efficiency_Err = Efficiency*std::sqrt((1./nSignalSelected)+(1./nSelected));
      }

      StatsHistograms[iVar][PurityIndex]->SetBinContent(xBin,Purity);
      StatsHistograms[iVar][PurityIndex]->SetBinError(xBin,Purity_Err);

      StatsHistograms[iVar][EfficiencyIndex]->SetBinContent(xBin,Efficiency);
      StatsHistograms[iVar][EfficiencyIndex]->SetBinError(xBin,Efficiency_Err);
    }

    double Min = 1e8;
    double Max = -1e8;

    for (int xBin=1;xBin<=Variables[iVar].Binning_.GetNbinsX();xBin++) {

      double Purity = StatsHistograms[iVar][PurityIndex]->GetBinContent(xBin);
      double Purity_Err = StatsHistograms[iVar][PurityIndex]->GetBinError(xBin);
      double Efficiency = StatsHistograms[iVar][EfficiencyIndex]->GetBinContent(xBin);
      double Efficiency_Err = StatsHistograms[iVar][EfficiencyIndex]->GetBinError(xBin);

      std::cout << std::setw(30) << Variables[iVar].Name_ << " " << std::setw(3) << xBin << " || Eff. = " << std::setw(10) << Efficiency << " +/- " << std::setw(10) << Efficiency_Err << " | Pur. = " << std::setw(10) << Purity << " +/- " << Purity_Err << std::endl;

      double MinBinContent = ((Purity-Purity_Err) < (Efficiency-Efficiency_Err)) ? Purity-Purity_Err : Efficiency-Efficiency_Err;
      if (MinBinContent < Min) Min = MinBinContent;

      double MaxBinContent = ((Purity+Purity_Err) > (Efficiency+Efficiency_Err)) ? Purity+Purity_Err : Efficiency+Efficiency_Err;
      if (MaxBinContent > Max) Max = MaxBinContent;
    }

    TCanvas* Canv = new TCanvas(Variables[iVar].Name_.c_str());

    StatsHistograms[iVar][PurityIndex]->GetYaxis()->SetTitle("%");
    StatsHistograms[iVar][PurityIndex]->GetXaxis()->SetTitle(Variables[iVar].Name_.c_str());
    StatsHistograms[iVar][PurityIndex]->SetLineColor(kBlack);
    StatsHistograms[iVar][EfficiencyIndex]->SetLineColor(kBlue);
    StatsHistograms[iVar][PurityIndex]->SetLineWidth(2);
    StatsHistograms[iVar][EfficiencyIndex]->SetLineWidth(2);

    Min = Min - 0.3*(Max-Min);
    Min = (Min < 0) ? 0 : Min;
    Max = Max + 0.3*(Max-Min);

    StatsHistograms[iVar][PurityIndex]->GetYaxis()->SetRangeUser(Min,Max);

    StatsHistograms[iVar][PurityIndex]->Draw("E1");
    StatsHistograms[iVar][EfficiencyIndex]->Draw("SAME E1");

    TLegend Legend = TLegend(0.7,0.8,0.9,0.9);
    Legend.AddEntry(StatsHistograms[iVar][PurityIndex],"Purity","l");
    Legend.AddEntry(StatsHistograms[iVar][EfficiencyIndex],"Efficiency","l");
    Legend.Draw();

    Canv->Print((Variables[iVar].Name_+".pdf").c_str());
  }

  //========================================================================================================

  TFile* File = new TFile("EfficiencyAndPurity.root","RECREATE");
  for (size_t iVar=0;iVar<Variables.size();iVar++) {
    File->cd();
    File->mkdir(Variables[iVar].Name_.c_str());
    File->cd(Variables[iVar].Name_.c_str());

    StatsHistograms[iVar][PurityIndex]->Write("Purity");
    StatsHistograms[iVar][EfficiencyIndex]->Write("Efficiency");
  }
  File->Close();
  
}

int main() {
  CalculateEfficiency();
}
