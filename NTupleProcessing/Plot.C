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

#include <vector>
#include <iostream>
#include <iomanip>

struct Variable {
  std::string BranchName_;
  std::string DrawCut_;
  std::string ShortDrawCut_;
  TH1F Binning_;
  double LowBoundSelCut_ = -999999.;
  double UppBoundSelCut_ = 999999.;
};

int Plot() {

  //================================================================================
  //Define event weights, cuts and ntuple file types considered

  std::vector<int> runs = {1,2,3};
  //std::vector<int> runs = {1};
  
  const std::string mc_event_weight = DEFAULT_MC_EVENT_WEIGHT;

  std::vector<NtupleFileType> FileTypes = {NtupleFileType::kOnBNB,NtupleFileType::kExtBNB,NtupleFileType::kNumuMC,NtupleFileType::kIntrinsicNueMC,NtupleFileType::kDirtMC};
  //std::vector<NtupleFileType> FileTypes = {NtupleFileType::kOnBNB,NtupleFileType::kNumuMC};

  //================================================================================
  //Variables we want to plot:

  std::vector<Variable> Variables;

  /*
  Variables.emplace_back(Variable{"nslice","1==1","All",TH1F("nslice","1==1;nslice;events",4,-1,3),1.,2.});
  */

  Variables.emplace_back(Variable{"mc_nu_vtx_sce_x","nslice==1","nslice==1",TH1F("mc_nu_vtx_sce_x",";mc_nu_vtx_sce_x;events",50,-50.,300.),10.,246.35});
  Variables.emplace_back(Variable{"mc_nu_vtx_sce_y","nslice==1","nslice==1",TH1F("mc_nu_vtx_sce_y",";mc_nu_vtx_sce_y;events",50,-150.,150.),-106.5,106.5});
  Variables.emplace_back(Variable{"mc_nu_vtx_sce_z","nslice==1","nslice==1",TH1F("mc_nu_vtx_sce_z",";mc_nu_vtx_sce_z;events",50,-50.,1100.),10.,1026.8});

  /*
  Variables.emplace_back(Variable{"trk_score_v","nslice==1 && CC1mu2p0pi_sel_reco_vertex_in_FV==1","1 slice, in FV",TH1F("trk_score_v",";trk_score_v;events",30,0.,1.),0.8});
  Variables.emplace_back(Variable{"trk_distance_v","nslice==1 && CC1mu2p0pi_sel_reco_vertex_in_FV==1","1 slice, in FV",TH1F("trk_distance_v",";trk_distance_v;events",30,0.,10.),-999999.,4.0});

  Variables.emplace_back(Variable{"trk_llr_pid_score_v","nslice==1 && CC1mu2p0pi_sel_npfps_eq_3==1 && CC1mu2p0pi_sel_ntracks_eq_3==1","1 slice, in FV, 3 Tracks",TH1F("trk_llr_pid_score_v_3PFPs",";trk_llr_pid_score_v;events",30,-1.,1.),0.2,999999.});

  Variables.emplace_back(Variable{"CC1mu2p0pi_MuonMomentumVector_Reco.Mag()","nslice==1 && CC1mu2p0pi_sel_reco_vertex_in_FV==1 && CC1mu2p0pi_sel_npfps_eq_3==1 && CC1mu2p0pi_sel_ntracks_eq_3==1 && CC1mu2p0pi_sel_correctparticles==1","1 slice, in FV, 2 protons, 1 muon",TH1F("CC1mu2p0piMuonMomentumVectorReco",";CC1mu2p0piMuonMomentumVectorReco.Mag();events",60,0.05,2.0),0.1,1.2});
  Variables.emplace_back(Variable{"CC1mu2p0pi_LeadingProtonMomentumVector_Reco.Mag()","nslice==1 && CC1mu2p0pi_sel_reco_vertex_in_FV==1 && CC1mu2p0pi_sel_npfps_eq_3==1 && CC1mu2p0pi_sel_ntracks_eq_3==1 && CC1mu2p0pi_sel_correctparticles==1","1 slice, in FV, 2 protons, 1 muon",TH1F("CC1mu2p0piLeadingProtonMomentumVectorReco",";CC1mu2p0piLeadingProtonMomentumVectorReco.Mag();events",60,0.05,2.0),0.3,1.0});
  Variables.emplace_back(Variable{"CC1mu2p0pi_RecoilProtonMomentumVector_Reco.Mag()","nslice==1 && CC1mu2p0pi_sel_reco_vertex_in_FV==1 && CC1mu2p0pi_sel_npfps_eq_3==1 && CC1mu2p0pi_sel_ntracks_eq_3==1 && CC1mu2p0pi_sel_correctparticles==1","1 slice, in FV, 2 protons, 1 muon",TH1F("CC1mu2p0piRecoilProtonMomentumVectorReco",";CC1mu2p0piRecoilProtonMomentumVectorReco.Mag();events",60,0.05,2.0),0.3,1.0});

  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_Q2","CC1mu2p0pi_Selected==1","CC1mu2p0pi_Selected==1",TH1F("CC1mu2p0piRecoQ2",";RecoQ2;events;",20,0.,2.)});
  Variables.emplace_back(Variable{"CC1mu2p0pi_True_Q2","CC1mu2p0pi_Selected==1","CC1mu2p0pi_Selected==1",TH1F("CC1mu2p0piTrueQ2",";TrueQ2;events;",20,0.,2.)});
  */
  
  /*
  std::vector<double> BinEdges = {0.,0.1,0.2,0.3,0.4,0.5,0.6,1.0};
  Variables.emplace_back(Variable{"CC1mu2p0pi_True_Pt","CC1mu2p0pi_MC_Signal==1","CC1mu2p0pi_MC_Signal==1",TH1F("CC1mu2p0pi_True_Pt",";CC1mu2p0pi_True_Pt;events",BinEdges.size()-1,BinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_Pt","CC1mu2p0pi_Selected==1","CC1mu2p0pi_Selected==1",TH1F("CC1mu2p0pi_Reco_Pt",";CC1mu2p0pi_Reco_Pt;events",BinEdges.size()-1,BinEdges.data())});

  Variables.emplace_back(Variable{"CC1mu2p0pi_True_Pt","CC1mu2p0pi_MC_Signal==1 && CC1mu2p0pi_Selected==1","CC1mu2p0pi_MC_Signal==1 && CC1mu2p0pi_Selected==1",TH1F("CC1mu2p0pi_True_Pt_SigSel",";CC1mu2p0pi_True_Pt;events",BinEdges.size()-1,BinEdges.data())});
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_Pt","CC1mu2p0pi_MC_Signal==1 && CC1mu2p0pi_Selected==1","CC1mu2p0pi_MC_Signal==1 && CC1mu2p0pi_Selected==1",TH1F("CC1mu2p0pi_Reco_Pt_SigSel",";CC1mu2p0pi_Reco_Pt;events",BinEdges.size()-1,BinEdges.data())});
  */

  /*
  std::vector<double> BinEdges = {-1.0,-0.75,-0.50,-0.25,0.,0.25,0.50,0.75,1.0};
  Variables.emplace_back(Variable{"CC1mu2p0pi_Reco_CosPlPr","CC1mu2p0pi_Selected==1","CC1mu2p0pi_Selected==1",TH1F("CC1mu2p0pi_Reco_CosPlPr",";cos(#theta_{#vec{P_{L}} #bullet #vec{P_{R}}});events",BinEdges.size()-1,BinEdges.data())});
  */

  //================================================================================
  //Cateogries we want to plot by

  std::vector<std::string> CategoryName;
  std::vector<int> CategoryColor;
  std::vector<std::string> CategoryCutString;

  /*
  CategoryName.emplace_back("Proton");
  CategoryColor.emplace_back(46);
  CategoryCutString.emplace_back("backtracked_pdg==2212");

  CategoryName.emplace_back("Muon");
  CategoryColor.emplace_back(kCyan+2);
  CategoryCutString.emplace_back("TMath::Abs(backtracked_pdg)==13");

  CategoryName.emplace_back("Pion");
  CategoryColor.emplace_back(38);
  CategoryCutString.emplace_back("TMath::Abs(backtracked_pdg)==211");

  CategoryName.emplace_back("#gamma");
  CategoryColor.emplace_back(24);
  CategoryCutString.emplace_back("backtracked_pdg==22");

  CategoryName.emplace_back("Electron");
  CategoryColor.emplace_back(28);
  CategoryCutString.emplace_back("TMath::Abs(backtracked_pdg)==11");

  CategoryName.emplace_back("Kaon");
  CategoryColor.emplace_back(42);
  CategoryCutString.emplace_back("TMath::Abs(backtracked_pdg)==321");

  CategoryName.emplace_back("Other");
  CategoryColor.emplace_back(kGray);
  CategoryCutString.emplace_back("(backtracked_pdg==2112 || backtracked_pdg==0)");
  */

  CategoryName.emplace_back("Other");
  CategoryColor.emplace_back(kBlack);
  CategoryCutString.emplace_back("(CC1mu2p0pi_EventCategory==0 || CC1mu2p0pi_EventCategory==19 || CC1mu2p0pi_EventCategory==20 || CC1mu2p0pi_EventCategory==22)");

  CategoryName.emplace_back("OOFV");
  CategoryColor.emplace_back(kGray);
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==21");

  CategoryName.emplace_back("CC1#mu CCOth");
  CategoryColor.emplace_back(46);
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==18");

  CategoryName.emplace_back("CC1#muN#pi");
  CategoryColor.emplace_back(42);
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==17");

  CategoryName.emplace_back("CC1#mu0p");
  CategoryColor.emplace_back(21);
  CategoryCutString.emplace_back("(CC1mu2p0pi_EventCategory==1 || CC1mu2p0pi_EventCategory==2 || CC1mu2p0pi_EventCategory==3 || CC1mu2p0pi_EventCategory==4)");

  CategoryName.emplace_back("CC1#mu1p");
  CategoryColor.emplace_back(24);
  CategoryCutString.emplace_back("(CC1mu2p0pi_EventCategory==5 || CC1mu2p0pi_EventCategory==6 || CC1mu2p0pi_EventCategory==7 || CC1mu2p0pi_EventCategory==8)");

  CategoryName.emplace_back("CC1#mu(M>2)p");
  CategoryColor.emplace_back(28);
  CategoryCutString.emplace_back("(CC1mu2p0pi_EventCategory==13 || CC1mu2p0pi_EventCategory==14 || CC1mu2p0pi_EventCategory==15 || CC1mu2p0pi_EventCategory==16)");

  CategoryName.emplace_back("CC1#mu2p Other");
  CategoryColor.emplace_back(9);
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==12");

  CategoryName.emplace_back("CC1#mu2p QE");
  CategoryColor.emplace_back(38);
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==9");

  CategoryName.emplace_back("CC1#mu2p RES");
  CategoryColor.emplace_back(kCyan-3);
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==11");

  CategoryName.emplace_back("CC1#mu2p MEC");
  CategoryColor.emplace_back(kCyan+2);  
  CategoryCutString.emplace_back("CC1mu2p0pi_EventCategory==10");

  /*
  CategoryName.emplace_back("CCQE");
  CategoryColor.emplace_back(46);
  CategoryCutString.emplace_back("abs(mc_nu_pdg) == 14 && mc_ccnc == 0 && mc_interaction == 0");

  CategoryName.emplace_back("CCMEC");
  CategoryColor.emplace_back(kCyan+2);
  CategoryCutString.emplace_back("abs(mc_nu_pdg) == 14 && mc_ccnc == 0 && mc_interaction == 10");

  CategoryName.emplace_back("CCRES");
  CategoryColor.emplace_back(42);
  CategoryCutString.emplace_back("abs(mc_nu_pdg) == 14 && mc_ccnc == 0 && mc_interaction == 1");

  CategoryName.emplace_back("CCDIS");
  CategoryColor.emplace_back(21);
  CategoryCutString.emplace_back("abs(mc_nu_pdg) == 14 && mc_ccnc == 0 && mc_interaction == 2");

  CategoryName.emplace_back("#nu_{e}");
  CategoryColor.emplace_back(kOrange);
  CategoryCutString.emplace_back("abs(mc_nu_pdg) == 12 && mc_ccnc == 0");

  CategoryName.emplace_back("NC");
  CategoryColor.emplace_back(kGray);
  CategoryCutString.emplace_back("mc_ccnc == 1");
  */

  /*
  CategoryName.emplace_back("All");
  CategoryColor.emplace_back(46);
  CategoryCutString.emplace_back("1");
  */

  //================================================================================

  if (CategoryName.size() != CategoryCutString.size()) throw;

  //HistogramsToFill[NTupleType][PlotVariable][Category]
  std::vector< std::vector< std::vector<TH1F*> > > HistogramsToFill;
  HistogramsToFill.resize(FileTypes.size());
  for (size_t iFT=0;iFT<FileTypes.size();iFT++) {
    HistogramsToFill[iFT].resize(Variables.size());
    for (size_t iVar=0;iVar<Variables.size();iVar++) {
      HistogramsToFill[iFT][iVar].resize(CategoryName.size());
      for (int iCat=0;iCat<CategoryName.size();iCat++) {
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
	    //CutString = mc_event_weight + " * (" + Variables[iVar].DrawCut_ + ") * FDSWeight";
	  } else if (NTupleType == NtupleFileType::kExtBNB) {
	    CutString = TString(Form("%2.6f * ",ScaleFactor)) + "(" + Variables[iVar].DrawCut_ + ")";
	  } else {
	    CutString = Form("%2.6f * ",ScaleFactor) + mc_event_weight + " * (" + Variables[iVar].DrawCut_ + ")";
	  }

	  for (int iCat=0;iCat<CategoryName.size();iCat++) {
	    
	    //TString HistName = Variables[iVar]+"_"+Form("%i",histcounter);
	    TString HistName = Form("%i",histcounter);
	    TH1F* Hist = (TH1F*)(Variables[iVar].Binning_).Clone(HistName);
	    Hist->Sumw2();
	    histcounter += 1;

	    TString FullCutString;
	    
	    if (NTupleType == NtupleFileType::kOnBNB || NTupleType == NtupleFileType::kExtBNB || NTupleType == NtupleFileType::kDirtMC) {
	      FullCutString = CutString;
	    } else {
	      FullCutString = CutString + " * (" + CategoryCutString[iCat] + ")";
	    }

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
      for (size_t iCat=0;iCat<CategoryName.size();iCat++) {
	std::cout << std::setw(15) << fpm.ntuple_type_to_string(FileTypes[iFT]) << " " << std::setw(65) << Variables[iVar].BranchName_ << " " << std::setw(15) << CategoryName[iCat] << " | " << std::setw(15) << HistogramsToFill[iFT][iVar][iCat]->Integral() << " " << std::setw(15) << HistogramsToFill[iFT][iVar][iCat]->GetEntries() << std::endl;
      }
    }
  }

  //========================================================================================================

  //HistogramsToFill[NTupleType][PlotVariable][Category]
  for (size_t iVar=0;iVar<Variables.size();iVar++) {
    std::cout << std::setw(15) << Variables[iVar].BranchName_ << std::endl;

    //First let's get the BNBon and EXT hists
    //Don't integrate over categories because the cutstring is not differernt for each category 
    TH1F* BNBOnHist = (TH1F*)(Variables[iVar].Binning_).Clone();
    TH1F* ExtHist = (TH1F*)(Variables[iVar].Binning_).Clone();

    for (size_t iFT=0;iFT<FileTypes.size();iFT++) {
      if (FileTypes[iFT] == NtupleFileType::kOnBNB) {
	BNBOnHist->Add(HistogramsToFill[iFT][iVar][0]);
      }
      if (FileTypes[iFT] == NtupleFileType::kExtBNB || FileTypes[iFT] == NtupleFileType::kDirtMC) {
	ExtHist->Add(HistogramsToFill[iFT][iVar][0]);
      }
    }

    std::cout << std::setw(15) << "Data" << ": ";
    for (int xBin=1;xBin<=BNBOnHist->GetNbinsX();xBin++) {
      std::cout << std::setw(15) << BNBOnHist->GetBinContent(xBin) << ", ";
    }
    std::cout << ", " << BNBOnHist->Integral() << std::endl;

    std::cout << std::setw(15) << "Ext" << ": ";
    for (int xBin=1;xBin<=ExtHist->GetNbinsX();xBin++) {
      std::cout << std::setw(15) << ExtHist->GetBinContent(xBin) << ", ";
    }
    std::cout << ", " << ExtHist->Integral() << std::endl;

    TCanvas* Canv = new TCanvas(TString(Variables[iVar].Binning_.GetName())+"_Canv","");
    THStack* Stack = new THStack((Variables[iVar].BranchName_+"_Stack").c_str(),(";"+std::string(Variables[iVar].Binning_.GetXaxis()->GetTitle())+";Events").c_str());
    TLegend* Legend = new TLegend(0.1,0.9,0.9,1.0);
    Legend->SetNColumns(std::ceil(float(CategoryName.size()+2)/2.0));

    ExtHist->SetFillColor( 28 );
    ExtHist->SetLineColor( 28 );
    ExtHist->SetLineWidth( 2 );
    ExtHist->SetFillStyle( 3005 );
    Stack->Add(ExtHist);
    Legend->AddEntry(ExtHist,"EXT","f");

    std::vector<TH1F*> CatHists(CategoryName.size());

    for (size_t iCat=0;iCat<CategoryName.size();iCat++) {
      CatHists[iCat] = new TH1F(Variables[iVar].Binning_);
      CatHists[iCat]->Sumw2();

      for (size_t iFT=0;iFT<FileTypes.size();iFT++) {
	if (FileTypes[iFT] == NtupleFileType::kOnBNB || FileTypes[iFT] == NtupleFileType::kExtBNB || FileTypes[iFT] == NtupleFileType::kDirtMC) continue;
	CatHists[iCat]->Add(HistogramsToFill[iFT][iVar][iCat]);
      }
    }

    for (size_t iCat=0;iCat<CategoryName.size();iCat++) {
      CatHists[iCat]->SetFillColor(CategoryColor[iCat]);
      Legend->AddEntry(CatHists[iCat],CategoryName[iCat].c_str(),"f");
      Stack->Add(CatHists[iCat]);

      std::cout << std::setw(15) << CategoryName[iCat] << ": ";
      for (int xBin=1;xBin<=CatHists[iCat]->GetNbinsX();xBin++) {
	std::cout << std::setw(15) << CatHists[iCat]->GetBinContent(xBin) << ", ";
      }
      std::cout << ", " << CatHists[iCat]->Integral() << std::endl;
    }

    std::cout << std::setw(15) << "Total" << ": ";
    TH1 *StackTotal = (TH1*)Stack->GetStack()->Last();
    for (int xBin=1;xBin<=StackTotal->GetNbinsX();xBin++) {
      std::cout << std::setw(15) << StackTotal->GetBinContent(xBin) << ", ";
    }
    std::cout << ", " << StackTotal->Integral() << std::endl;
    
    Stack->Draw("hist");
    double Max = (Stack->GetMaximum() > BNBOnHist->GetMaximum()) ? Stack->GetMaximum() : BNBOnHist->GetMaximum();
    double Min = Stack->GetMinimum();
    Stack->SetMaximum(1.3*Max);

    BNBOnHist->SetMarkerColor(kBlack);
    BNBOnHist->SetLineColor(kBlack);
    BNBOnHist->SetLineWidth( 3 );

    BNBOnHist->Draw("SAME E");

    Legend->AddEntry(BNBOnHist,"Data","l");
    Legend->Draw("SAME");

    TLatex Text(.13,.85,Variables[iVar].ShortDrawCut_.c_str());
    Text.SetNDC(kTRUE);
    //Text.Draw("SAME");

    if (Variables[iVar].LowBoundSelCut_ != -999999. || Variables[iVar].UppBoundSelCut_ != 999999.) {

      double LowVal = (Variables[iVar].LowBoundSelCut_ < BNBOnHist->GetXaxis()->GetBinLowEdge(1)) ? BNBOnHist->GetXaxis()->GetBinLowEdge(1) : Variables[iVar].LowBoundSelCut_;
      double UppVal = (Variables[iVar].UppBoundSelCut_ > BNBOnHist->GetXaxis()->GetBinLowEdge(BNBOnHist->GetNbinsX()+1)) ? BNBOnHist->GetXaxis()->GetBinLowEdge(BNBOnHist->GetNbinsX()+1) : Variables[iVar].UppBoundSelCut_;

      if (Variables[iVar].LowBoundSelCut_ > BNBOnHist->GetXaxis()->GetBinLowEdge(1)) {
	TLine* LowLine = new TLine(LowVal,Min,LowVal,1.1*Max);
	LowLine->Draw();
	
	TArrow* LowArr = new TArrow(LowVal,1.1*Max,LowVal+0.3*(UppVal-LowVal),1.1*Max,0.03,"|>");
	LowArr->Draw();
      }

      if (Variables[iVar].UppBoundSelCut_ < BNBOnHist->GetXaxis()->GetBinLowEdge(BNBOnHist->GetNbinsX()+1)) {
	TLine* UppLine = new TLine(UppVal,Min,UppVal,1.1*Max);
	UppLine->Draw();
	
	TArrow* UppArr = new TArrow(UppVal,1.1*Max,UppVal-0.3*(UppVal-LowVal),1.1*Max,0.03,"|>");
	UppArr->Draw();
      }
    }

    Canv->Update();
    Canv->Print(TString(Variables[iVar].Binning_.GetName())+".pdf");
  }

}

int main() {
  Plot();
}
