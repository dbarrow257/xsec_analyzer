// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// STV analysis includes
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "PlotUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

using NFT = NtupleFileType;

void tutorial_slice_plots(std::string FPM_Config, std::string SYST_Config, std::string SLICE_Config, std::string Univ_Output, std::string Plot_OutputDir) {

  std::string Plot_Prefix = "SlicePlots";
  std::string Plot_Suffix = ".pdf";

  std::cout << "\nRunning Slice_Plots with options:" << std::endl;
  std::cout << "\tFPM_Config: " << FPM_Config << std::endl;
  std::cout << "\tSYST_Config: " << SYST_Config << std::endl;
  std::cout << "\tSLICE_Config: " <<  SLICE_Config << std::endl;
  std::cout << "\tUniv_Output: " << Univ_Output << std::endl;
  std::cout << "\tPlot_OutputDir: " << Plot_OutputDir << std::endl;
  std::cout << "\n" << std::endl;

  std::vector< std::pair<TString, std::vector<int> > > CategoriesToPlot;
  std::vector<int> FillColor;

  CategoriesToPlot.emplace_back(TString("Other"),std::vector<int>{0,19,20,22});
  FillColor.emplace_back(kBlack);

  CategoriesToPlot.emplace_back(TString("OOFV"),std::vector<int>{21});
  FillColor.emplace_back(kGray);

  CategoriesToPlot.emplace_back(TString("CC1#mu CCOth"),std::vector<int>{18});
  FillColor.emplace_back(46);

  CategoriesToPlot.emplace_back(TString("CC1#muN#pi"),std::vector<int>{17});
  FillColor.emplace_back(42);

  CategoriesToPlot.emplace_back(TString("CC1#mu0p"),std::vector<int>{1,2,3,4});
  FillColor.emplace_back(21);

  CategoriesToPlot.emplace_back(TString("CC1#mu1p"),std::vector<int>{5,6,7,8});
  FillColor.emplace_back(24);

  CategoriesToPlot.emplace_back(TString("CC1#mu(M>2)p"),std::vector<int>{13,14,15,16});
  FillColor.emplace_back(28);

  CategoriesToPlot.emplace_back(TString("CC1#mu2p Other"),std::vector<int>{12});
  FillColor.emplace_back(9);

  CategoriesToPlot.emplace_back(TString("CC1#mu2p QE"),std::vector<int>{9});
  FillColor.emplace_back(38);

  CategoriesToPlot.emplace_back(TString("CC1#mu2p RES"),std::vector<int>{11});
  FillColor.emplace_back(kCyan-3);

  CategoriesToPlot.emplace_back(TString("CC1#mu2p MEC"),std::vector<int>{10});
  FillColor.emplace_back(kCyan+2);
  
  // Show fractional uncertainties computed using these covariance matrices
  // in the ROOT plot. All configured fractional uncertainties will be
  // included in the output pgfplots file regardless of whether they appear
  // in this vector.
  std::vector< std::pair<std::string, std::string> > cov_mat_keys;

  /*  
  cov_mat_keys.emplace_back("total","Total");
  cov_mat_keys.emplace_back("detVar_total","Det. Var.");
  cov_mat_keys.emplace_back("flux","Flux");
  cov_mat_keys.emplace_back("reint","Reint.");
  cov_mat_keys.emplace_back("xsec_total","XSec");
  cov_mat_keys.emplace_back("POT","POT");
  cov_mat_keys.emplace_back("numTargets","Num. Targets");
  cov_mat_keys.emplace_back("MCstats","MC Stats.");
  cov_mat_keys.emplace_back("EXTstats","EXT Stats.");
  cov_mat_keys.emplace_back("BNBstats","BNB Stats.");

  std::string FinalPlot = "Total";
  */

  /*
  cov_mat_keys.emplace_back("xsec_AxFFCCQEshape","AxFFCCQEshape");
  cov_mat_keys.emplace_back("xsec_DecayAngMEC","DecayAngMEC");  
  cov_mat_keys.emplace_back("xsec_NormCCCOH","NormCCCOH");
  cov_mat_keys.emplace_back("xsec_NormNCCOH","NormNCCOH");
  cov_mat_keys.emplace_back("xsec_RPA_CCQE","RPA_CCQE");
  cov_mat_keys.emplace_back("xsec_ThetaDelta2NRad","ThetaDelta2NRad");
  cov_mat_keys.emplace_back("xsec_Theta_Delta2Npi","Theta_Delta2Npi");
  cov_mat_keys.emplace_back("xsec_VecFFCCQEshape","VecFFCCQEshape");
  cov_mat_keys.emplace_back("xsec_XSecShape_CCMEC","XSecShape_CCMEC");
  cov_mat_keys.emplace_back("xsec_xsr_scc_Fa3_SCC","xsr_scc_Fa3_SCC");
  cov_mat_keys.emplace_back("xsec_xsr_scc_Fv3_SCC","xsr_scc_Fv3_SCC");
  cov_mat_keys.emplace_back("xsec_multi","Multi-Sim");
  cov_mat_keys.emplace_back("xsec_total","XSec Total");

  std::string FinalPlot = "XSec Total";
  */

  cov_mat_keys.emplace_back("detVarLYatten","LYatten");
  cov_mat_keys.emplace_back("detVarLYdown","LYdown");
  cov_mat_keys.emplace_back("detVarLYrayl","LYrayl");
  cov_mat_keys.emplace_back("detVarRecomb2","Recomb2");
  cov_mat_keys.emplace_back("detVarSCE","SCE");
  cov_mat_keys.emplace_back("detVarWMAngleXZ","WMAngleXZ");
  cov_mat_keys.emplace_back("detVarWMAngleYZ","WMAngleYZ");
  cov_mat_keys.emplace_back("detVarWMX","WMX");
  cov_mat_keys.emplace_back("detVarWMYZ","WMYZ");
  cov_mat_keys.emplace_back("detVar_total","Det. Var. Total");

  std::string FinalPlot = "Det. Var. Total";

//===========================================================================================================================================================================

  std::cout << "Cateogries grouped together ----" << std::endl;
  for (const auto& Cat: CategoriesToPlot) {
    std::cout << Cat.first << ":" << "\n\t";
    for (const auto& iCat: Cat.second) {
      std::cout << iCat << ",";
    }
    std::cout << std::endl;
  }
  std::cout << "\n\n" << std::endl;

#ifdef USE_FAKE_DATA
  // Initialize the FilePropertiesManager and tell it to treat the NuWro
  // MC ntuples as if they were data
  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( FPM_Config );
#endif
  
  // Check that we can read the universe output file
  TFile* temp_file = new TFile(Univ_Output.c_str(), "read");
  if (!temp_file || temp_file->IsZombie()) {
    std::cerr << "Could not read file: " << Univ_Output << std::endl;
    throw;
  }
  delete temp_file;
  
  auto* syst_ptr = new MCC9SystematicsCalculator(Univ_Output, SYST_Config);
  auto& syst = *syst_ptr;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

  #ifdef USE_FAKE_DATA
    // Add the EXT to the "data" when working with fake data
    reco_bnb_hist->Add( reco_ext_hist );
  #endif

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  std::string PlotFileName;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    
    const auto& slice = sb.slices_.at( sl_idx );
    auto& SliceVar = sb.slice_vars_.at( slice.active_var_indices_.front() );
    std::string SliceLatexName = SliceVar.latex_name_;
    std::string SliceLatexUnits = SliceVar.latex_units_;

    if (SliceLatexUnits != "") {
      SliceLatexName += " ("+SliceLatexUnits+")";
    }

    std::string SliceVariableName = SliceVar.name_;
    SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());

    //===========================================================================================================================================================================
    //Plot the slices
    
    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    const auto& eci = EventCategoryInterpreter::Instance();
    eci.set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& cat_map = eci.label_map();

    TLegend* Legend = new TLegend(0.1,0.9,0.9,1.0);
    Legend->SetNColumns(std::ceil(float(CategoriesToPlot.size()+1)/2.0));
    Legend->AddEntry((slice_ext->hist_).get(),"EXT","f");

    int counter = 0;
    for (const auto& CatGroup : CategoriesToPlot) {
      TString CatName = CatGroup.first;
      std::vector<int> Cats = CatGroup.second;
      int nCats = Cats.size();

      if (nCats == 0) {
	std::cerr << "Given a category group with no categories" << std::endl;
	throw;
      }

      TH1D* temp_mc_hist = category_hist->ProjectionY("temp_mc_hist", Cats[0]+1, Cats[0]+1);
      temp_mc_hist->SetDirectory( nullptr );

      if (nCats > 1) {
	for (int iCat=1;iCat<nCats;iCat++) {
	  TH1D* temp_mc_hist_iCat = category_hist->ProjectionY(Form("temp_mc_hist_%i",iCat), Cats[iCat]+1,Cats[iCat]+1);
	  temp_mc_hist->Add(temp_mc_hist_iCat);
	}
      }

      SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(*temp_mc_hist, slice);
      eci.set_mc_histogram_style((EventCategory)Cats[0],temp_slice_mc->hist_.get());
      temp_slice_mc->hist_->SetFillColor(FillColor[counter]);
      temp_slice_mc->hist_->SetLineColor(FillColor[counter]);
      counter ++ ;

      slice_pred_stack->Add(temp_slice_mc->hist_.get());

      Legend->AddEntry((temp_slice_mc->hist_).get(),CatName,"f");
    }
    //=========================================================

    TCanvas* c1 = new TCanvas;
    slice_bnb->hist_->SetLineColor( kBlack );
    slice_bnb->hist_->SetLineWidth( 3 );
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 1.2 );
    slice_bnb->hist_->SetStats( false );
    double ymax = std::max(slice_bnb->hist_->GetMaximum(),slice_mc_plus_ext->hist_->GetMaximum() ) * 1.5;
    slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );
    
    slice_bnb->hist_->SetTitle("");
    slice_bnb->hist_->GetYaxis()->SetTitle("Events/Bin");
    slice_bnb->hist_->GetXaxis()->SetTitle(SliceLatexName.c_str());

    slice_bnb->hist_->Draw( "e" );

    slice_pred_stack->Draw( "hist same" );

    slice_mc_plus_ext->hist_->SetLineWidth( 3 );
    slice_mc_plus_ext->hist_->SetLineColor(kBlack);
    slice_mc_plus_ext->hist_->SetFillColor(kBlack);
    slice_mc_plus_ext->hist_->SetFillStyle(3353);
    slice_mc_plus_ext->hist_->Draw( "same e2" );

    slice_bnb->hist_->Draw( "same e" );
    slice_bnb->hist_->Draw( "SAME AXIS" );

    Legend->Draw("SAME");

    c1->Update();
    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + "_" + SliceVariableName + Plot_Suffix;
    c1->SaveAs(PlotFileName.c_str());

    //===========================================================================================================================================================================
    //Plot the fractional uncertainty

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(slice.hist_->Clone("slice_hist") );
    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    double Max = -1e8;

    // Loop over the various systematic uncertainties
    int color = 1;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(*reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // Check whether the current covariance matrix name is present in
      // the vector defined above this loop. If it isn't, don't bother to
      // plot it, and just move on to the next one.
      bool AddToMap = false;
      std::string FancyName = "";
      for (size_t iCov=0;iCov<cov_mat_keys.size();iCov++) {
	if (cov_mat_keys[iCov].first == key) {
	  AddToMap = true;
	  FancyName = cov_mat_keys[iCov].second;
	  break;
	}
      }
      if (!AddToMap) continue;

      frac_uncertainty_hists[ FancyName ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 10;

      slice_for_syst->hist_->SetLineColor( color );
      slice_for_syst->hist_->SetLineWidth( 3 );

      for (int iBin=1;iBin<=slice_for_syst->hist_->GetNbinsX();iBin++) {
	double Content = slice_for_syst->hist_->GetBinContent(iBin);
	if (Content > Max) Max = Content;
      }
    }

    TCanvas* c2 = new TCanvas;
    TLegend* lg2 = new TLegend( 0.1, 0.9, 0.9, 1.0 );
    lg2->SetNColumns(std::ceil(float(cov_mat_keys.size()+1)/2.0));

    bool isFirst = true;
    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;

      hist->SetTitle("");
      hist->GetYaxis()->SetTitle("Fractional Uncertainty (%)");
      hist->GetXaxis()->SetTitle(SliceLatexName.c_str());
      hist->GetYaxis()->SetRangeUser(0.,Max*1.3);
      hist->SetStats( false );
      hist->SetLineWidth( 3 );

      if (name == FinalPlot) {
	hist->SetLineColor( kBlack );
      } else {
	lg2->AddEntry( hist, name.c_str(), "l" );
      }

      if (isFirst) {
	isFirst = false;
	hist->Draw( "hist" );
      } else {
	hist->Draw( "same hist" );
      }
    }

    lg2->AddEntry( frac_uncertainty_hists.at(FinalPlot), FinalPlot.c_str(), "l" );
    lg2->Draw( "same" );

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + "_" + SliceVariableName + "_FracUncer" + Plot_Suffix;
    c2->SaveAs(PlotFileName.c_str());

  } // slices

}

int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    std::cout << "Usage: Slice_Plots FPM_CONFIG"
	      << " SYST_Config SLICE_Config Univ_Output Plot_OutputDir\n";
    return 1;
  }

  std::string list_file_name( argv[1] );
  std::string univmake_config_file_name( argv[2] );
  std::string output_file_name( argv[3] );

  std::string FPM_Config( argv[1] );
  std::string SYST_Config( argv[2] );
  std::string SLICE_Config( argv[3] );
  std::string Univ_Output( argv[4] );
  std::string Plot_OutputDir( argv[5] );

  tutorial_slice_plots(FPM_Config, SYST_Config, SLICE_Config, Univ_Output, Plot_OutputDir);
  return 0;
}
