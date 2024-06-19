// Standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

#include <algorithm>

// ROOT includes
#include "TCanvas.h"
#include "TLegend.h"

// STV analysis includes
#include "CrossSectionExtractor.hh"
#include "PGFPlotsDumpUtils.hh"
#include "PlotUtils.hh"
#include "SliceBinning.hh"
#include "SliceHistogram.hh"

#include "Functions.h"

//Useful DEBUG options which can be turned on/off
std::string PlotExtension = ".pdf";
std::string TextExtension = ".txt";
bool DumpToText = false;
bool DumpToPlot = false;

double MigrationMatrixThreshold = 0.50;

void Unfolder(std::string XSEC_Config, std::string SLICE_Config, std::string OutputDirectory, std::string OutputFileName) {

  // Make results in both Event Count units and then the XSec units
  std::vector<std::string> ResultTypes(2);
  ResultTypes[0] = "EventCountUnits";
  ResultTypes[1] = "XsecUnits";
  size_t nResultTypes = ResultTypes.size();

  // Make predictions in both AC smeared and truth distributions
  std::vector<std::string> SmearTypes(2);
  SmearTypes[0] = "Unsmeared";
  SmearTypes[1] = "Smeared";
  size_t nSmearTypes = SmearTypes.size();

  std::cout << "\nRunning Unfolder.C with options:" << std::endl;
  std::cout << "\tXSEC_Config: " << XSEC_Config << std::endl;
  std::cout << "\tSLICE_Config: " << SLICE_Config << std::endl;
  std::cout << "\tOutputDirectory: " << OutputDirectory << std::endl;
  std::cout << "\tOutputFileName: " << OutputFileName << std::endl;
  std::cout << "\n" << std::endl;

  OutputFileName = OutputDirectory+"/"+OutputFileName;

  // Use a CrossSectionExtractor object to handle the systematics and unfolding
  auto extr = std::make_unique< CrossSectionExtractor >( XSEC_Config );

  //Grab the systematics calculator from the extractor
  auto syst_calc = extr->get_syst_calc();

  // Plot slices of the unfolded result                                                                                                                                                                             
  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr; 

  // Check that the output file can be written to
  TFile* File = new TFile(OutputFileName.c_str(), "RECREATE");
  if (!File || File->IsZombie()) {
    std::cerr << "Could not write to output file:" << OutputFileName << std::endl;
    throw;
  }

  // Get the tuned GENIE CV prediction in each true bin (including the
  // background true bins)
  TH1D* genie_cv_truth = syst_calc->cv_universe().hist_true_.get();
  
  int num_true_signal_bins = 0;
  for ( int t = 0; t < syst_calc->true_bins_.size(); ++t ) {
    const auto& tbin = syst_calc->true_bins_.at( t );
    if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
  }

  // Save the GENIE CV model (before A_C multiplication) using a column vector
  // of event counts
  TMatrixD genie_cv_truth_vec( num_true_signal_bins, 1 );
  for ( int b = 0; b < num_true_signal_bins; ++b ) {
    double true_evts = genie_cv_truth->GetBinContent( b + 1 );
    genie_cv_truth_vec( b, 0 ) = true_evts;
  }

  // Figure out if we are using fake data
  // If present, then get the fake data event counts in each true bin
  // (including the background true bins). We hope to approximately reproduce
  // these event counts in the signal true bins via unfolding the fake data.
  const auto& fake_data_univ = syst_calc->fake_data_universe();
  TH1D* fake_data_truth_hist = nullptr;

  bool using_fake_data = false;
  if ( fake_data_univ ) {
    using_fake_data = true;
    fake_data_truth_hist = fake_data_univ->hist_true_.get();
  }

  // Save the fake data truth (before A_C multiplication) using a column vector
  // of event counts
  TMatrixD fake_data_truth( num_true_signal_bins, 1 );
  if ( using_fake_data ) {
    for ( int b = 0; b < num_true_signal_bins; ++b ) {
      double true_evts = fake_data_truth_hist->GetBinContent( b + 1 );
      fake_data_truth( b, 0 ) = true_evts;
    }
  }

  std::cout << "\n\nSaving inputs -----------------" << std::endl;
  File->cd();
  File->mkdir("Inputs");
  File->cd("Inputs");

  auto smearcept = syst_calc->get_cv_smearceptance_matrix();
  smearcept->Write("smearcept");

  auto true_signal = syst_calc->get_cv_true_signal();
  true_signal->Write("true_signal");

  auto meas = syst_calc->get_measured_events();
  const auto& data_signal = meas.reco_signal_;
  const auto& data_covmat = meas.cov_matrix_;

  const auto& MC_selected = meas.reco_mc_plus_ext_;
  MC_selected->Write("true_selected");

  data_signal->Write("data_signal");
  data_covmat->Write("data_covmat");

  File->mkdir("Inputs/CovarianceMatrix");
  File->cd("Inputs/CovarianceMatrix");
  const auto covariances = syst_calc->get_covariances();
  for ( const auto& matrix_pair : *covariances ) {
    const std::string& matrix_key = matrix_pair.first;
    auto temp_cov_mat = matrix_pair.second.get_matrix();
    
    temp_cov_mat->Write(matrix_key.c_str());
  }

  //========================================================================================================================================
  //Grab the pre-unfolding MC prediction slices, also grab the total covariance matrix slices

  File->mkdir("Inputs/MCSlices");
  File->cd("Inputs/MCSlices");

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    auto& Slice = sb.slices_.at( sl_idx );

    //The following line will fall over if several active variables are used per Slice                                                                                                                              
    auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );

    std::string SliceVariableName = SliceVar.name_;
    SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());

    File->mkdir(("Inputs/MCSlices/"+SliceVariableName).c_str());
    File->cd(("Inputs/MCSlices/"+SliceVariableName).c_str());

    SliceHistogram* Slice_MC = SliceHistogram::make_slice_histogram( *MC_selected, Slice, nullptr );
    TH1* SliceHist = Slice_MC->hist_.get();

    SliceHist->Write(SliceVariableName.c_str());
  }

  File->mkdir("Inputs/TotalCovarianceMatrix");
  File->cd("Inputs/TotalCovarianceMatrix");

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    auto& Slice = sb.slices_.at( sl_idx );
    //The following line will fall over if several active variables are used per Slice
    auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );

    std::string SliceVariableName = SliceVar.name_;
    SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());

    for ( const auto& matrix_pair : *covariances ) {
      const std::string& matrix_key = matrix_pair.first;
      auto temp_cov_mat = matrix_pair.second.get_matrix();
      
      if (matrix_key != "total") continue;

      SliceHistogram* Slice_Cov = SliceHistogram::make_slice_histogram( *true_signal, Slice, temp_cov_mat.get() );
      auto SliceCovMat_Ptr = Slice_Cov->cmat_.get_matrix();
      TMatrixD* SliceCovMat_Mat = SliceCovMat_Ptr.get();

      SliceCovMat_Mat->Write(SliceVariableName.c_str());
    }
  }
  
  //========================================================================================================================================

  /*
    ToDo: Add this to output
    syst_->total_bnb_data_pot_;    
  */

  std::cout << "Total POT:" << syst_calc->total_bnb_data_pot_ << std::endl;

  //========================================================================================================================================
  std::cout << "\n\nCalculating Migration Matrix -----------------" << std::endl;
  std::cout << "Migration matrix threshold set to:" << MigrationMatrixThreshold << std::endl;

  File->cd();
  File->mkdir("Inputs/MigrationMatrix");
  File->cd("Inputs/MigrationMatrix");

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    auto& Slice = sb.slices_.at( sl_idx );
    //The following line will fall over if several active variables are used per Slice                                                                                                                              
    auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );

    std::string SliceVariableName = SliceVar.name_;
    SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());
    if (SliceVariableName == "recobinnumber") continue;

    TMatrixD TrueSignal_Mat = *(syst_calc->get_cv_true_signal().get());
    const TMatrixD* smearcept_Mat = smearcept.get();

    SliceHistogram* SliceHist = SliceHistogram::make_slice_histogram( TrueSignal_Mat, Slice, smearcept_Mat );
    TH1* TrueSignal_Hist = SliceHist->hist_.get();
    auto SmearceptSlice = SliceHist->cmat_.get_matrix();
    TMatrixD* SmearceptSlice_Mat = SmearceptSlice.get();

    TMatrixD* MigrationMatrix = new TMatrixD(*SmearceptSlice_Mat);
    for (int yBin=0;yBin<TrueSignal_Hist->GetNbinsX();yBin++) {
      double Sum = 0.;
      for (int xBin=0;xBin<TrueSignal_Hist->GetNbinsX();xBin++) {
	MigrationMatrix->operator()(xBin,yBin) = MigrationMatrix->operator()(xBin,yBin)/TrueSignal_Hist->GetBinContent(yBin+1);
	Sum += MigrationMatrix->operator()(xBin,yBin);
      }
      for (int xBin=0;xBin<TrueSignal_Hist->GetNbinsX();xBin++) {
	MigrationMatrix->operator()(xBin,yBin) = MigrationMatrix->operator()(xBin,yBin)/Sum;
      } 
    }    

    for (int xBin=0;xBin<TrueSignal_Hist->GetNbinsX();xBin++) {
      if (MigrationMatrix->operator()(xBin,xBin) < MigrationMatrixThreshold) {
	std::cout << "Migration Matrix Diagonal Bin Below Threshold! : Slice:" << sl_idx << "(" << std::setw(20) << SliceVariableName << ")" << " , Bin:" << xBin  << " , Value:" << MigrationMatrix->operator()(xBin,xBin) << std::endl;
      }
    }

    MigrationMatrix->Write(SliceVariableName.c_str());
  }

  //========================================================================================================================================
  std::cout << "\n\nUnfolding -----------------" << std::endl;

  auto xsec = extr->get_unfolded_events();
  double conv_factor = extr->conversion_factor();
  std::cout << "Conversion factor:" << 1.0/conv_factor << std::endl;
  const auto& pred_map = extr->get_prediction_map();

  //Grab the additional smearing matrix to apply to results:
  TMatrixD AC_matrix = *xsec.result_.add_smear_matrix_;
  TH2D* AC_matrix_TH2 = new TH2D("ACMatrix","",AC_matrix.GetNrows(),0,AC_matrix.GetNrows(),AC_matrix.GetNrows(),0,AC_matrix.GetNrows());
  TM2TH2(AC_matrix, AC_matrix_TH2);

  //Smear GENIE CV truth
  TMatrixD genie_cv_truth_smeared( AC_matrix, TMatrixD::kMult, genie_cv_truth_vec );

  //Smear FDS
  TMatrixD fake_data_truth_smeared( AC_matrix, TMatrixD::kMult, fake_data_truth );

  std::cout << "\n\nSaving results -----------------" << std::endl;

  std::cout << "Output file - " << OutputFileName << std::endl;
  if (DumpToText) std::cout << "\tDumping plots to " << TextExtension << " files" << std::endl;
  if (DumpToPlot) std::cout << "\tDumping plots to " << PlotExtension << " files" << std::endl;
  std::cout << "\n" << std::endl;
  
  //Loop over the two ResultTypes (Event Counts and Xsec Units)
  for (size_t iRT=0;iRT<nResultTypes;iRT++) {

    std::string RT = ResultTypes[iRT];

    //======================================================================================
    //Loop over all covariance matrices stored in the unfolded result and save them

    File->cd();
    File->mkdir((std::string("Outputs/")+RT+"/CovarianceMatrix").c_str());
    File->cd((std::string("Outputs/")+RT+"/CovarianceMatrix").c_str());

    // Convert units on the covariance matrices one-by-one and dump them
    for ( const auto& cov_pair : xsec.unfolded_cov_matrix_map_ ) {
      const auto& name = cov_pair.first;
      TMatrixD temp_cov_matrix = *cov_pair.second;
      // Note that we need to square the unit conversion factor for the
      // covariance matrix elements

      if (RT == "XsecUnits") {
	temp_cov_matrix *= std::pow( 1.0 / conv_factor, 2 );
      }

      temp_cov_matrix.Write(name.c_str());
      if (DumpToText) dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_" + name + TextExtension, temp_cov_matrix );
      if (DumpToPlot) draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_" + name + PlotExtension, temp_cov_matrix, (name+" matrix").c_str(), "Bin Number", "Bin Number", "COLZ");
    }

    // No unit conversions are necessary for the unfolding, error propagation,
    // and additional smearing matrices since they are dimensionless
    TMatrixD temp_unfolding_matrix = *xsec.result_.unfolding_matrix_;
    TMatrixD temp_err_prop_matrix = *xsec.result_.err_prop_matrix_;
    TMatrixD temp_add_smear_matrix = *xsec.result_.add_smear_matrix_;
    TMatrixD temp_response_matrix = *xsec.result_.response_matrix_;
    temp_unfolding_matrix.Write("UnfoldingMatrix");
    temp_err_prop_matrix.Write("ErrorPropagationMatrix");
    temp_add_smear_matrix.Write("AdditionalSmearingMatrix");
    temp_response_matrix.Write("ResponseMatrix");

    if (DumpToText) {
      dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_unfolding"+TextExtension, temp_unfolding_matrix );
      dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_err_prop"+TextExtension, temp_err_prop_matrix );
      dump_text_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_add_smear"+TextExtension, temp_add_smear_matrix );
    }
    if (DumpToPlot) {
      draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_unfolding"+PlotExtension, temp_unfolding_matrix, "Unfolding matrix", "Bin Number", "Bin Number", "COLZ");
      draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_err_prop"+PlotExtension, temp_err_prop_matrix, "Error Propagation matrix", "Bin Number", "Bin Number", "COLZ");
      draw_matrix( OutputDirectory+"/"+RT+"_mat_table_cov_add_smear"+PlotExtension, temp_add_smear_matrix, "Addition Smearing matrix", "Bin Number", "Bin Number", "COLZ");
    }

    //====================================================================================== 
    //Loop over all the slices taken from the config and save the unfolded distribution
    
    for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
      auto& Slice = sb.slices_.at( sl_idx );
      //The following line will fall over if several active variables are used per Slice
      auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );
      
      std::string SliceVariableName = SliceVar.name_;
      SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());
      
      File->cd();
      File->mkdir((std::string("Outputs/")+RT+"/"+SliceVariableName).c_str());
      File->cd((std::string("Outputs/")+RT+"/"+SliceVariableName).c_str());
      
      //====================================================================================== 
      //Save unfolded distribution for each covariance matrix
      
      // Make a histogram showing the unfolded counts in the current slice
      // for a particular covariance matrix being used to define the uncertainties
      for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
	const auto& uc_name = uc_pair.first;
	const auto& uc_matrix = uc_pair.second;
	
	SliceHistogram* Slice_unf = SliceHistogram::make_slice_histogram( *xsec.result_.unfolded_signal_, Slice, uc_matrix.get() );
	TH1* SliceHist = Slice_unf->hist_.get();
	
	SliceHist->Scale(1.0,"width");
	SliceHist->GetYaxis()->SetTitle("Events/Bin Width");
	
	if (RT == "XsecUnits") {
	  SliceHist->Scale(1.0 / conv_factor);
	  SliceHist->GetYaxis()->SetTitle("Differential Cross Section [10^{-38} cm^{2}/Ar/Bin Width]");
	}
	SliceHist->Write((SliceVariableName+"_"+uc_name).c_str());
	
	if (DumpToText) dump_text_column_vector( OutputDirectory+"/"+RT+"_vec_table_unfolded_signal_"+uc_name+TextExtension, *xsec.result_.unfolded_signal_ );
	if (DumpToPlot) draw_column_vector( OutputDirectory+"/"+RT+"_vec_table_unfolded_signal_"+uc_name+PlotExtension, *xsec.result_.unfolded_signal_, "Unfolded Signal", "Bin Number", "Cross Section [#times 10^{-38} cm^{2}]");
      }
      
      //====================================================================================== 
      //Loop over all the prediction slices
      
      for (size_t iST=0;iST<nSmearTypes;iST++) {
	
	std::string ST = SmearTypes[iST];

	TMatrixD* GENIE_CV;
	TMatrixD* FDS;

	if (ST == "Smeared") {
	  GENIE_CV = &genie_cv_truth_smeared;
	  if (using_fake_data) FDS = &fake_data_truth_smeared;
	} else {
	  GENIE_CV = &genie_cv_truth_vec;
	  if (using_fake_data) FDS = &fake_data_truth;
	}

	for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
	  const auto& uc_name = uc_pair.first;
	  const auto& uc_matrix = uc_pair.second;
	  
	  if (!(uc_name == "total" || uc_name == "FDS_Uncer")) continue;
	  
	  File->cd();
	  File->mkdir(("Predictions/"+RT+"/"+ST+"/"+SliceVariableName+"/"+uc_name).c_str());
	  File->cd(("Predictions/"+RT+"/"+ST+"/"+SliceVariableName+"/"+uc_name).c_str());
	  
	  //====================================================================================== 
	  //Get the unfolded Slice Histogram	  
	  SliceHistogram* Slice_unf = SliceHistogram::make_slice_histogram( *xsec.result_.unfolded_signal_, Slice, uc_matrix.get() );
	  TH1* Slice_unf_hist = Slice_unf->hist_.get();
	  
	  auto SliceCovMat_Ptr = Slice_unf->cmat_.get_matrix();
	  TMatrixD* SliceCovMat_Mat = SliceCovMat_Ptr.get();
	  
	  //====================================================================================== 
	  //Get the MC Predictions 
	  
	  //CV Slice
	  SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram( *GENIE_CV, Slice, nullptr );
	  TH1* genie_cv_truth_smeared_slice = slice_cv->hist_.get();
	  
	  const auto& cv_chi2_result = Slice_unf->get_chi2( *slice_cv );
	  if (iRT == 0) {
	    std::cout << std::setw(15) << uc_name << " | " << std::setw(15) << ST << " | " << std::setw(30) << SliceVariableName << " | " << std::setw(15) << "MicroBooNE Tune" << " | " << std::setw(15) << cv_chi2_result.chi2_ << "/" << std::setw(3) << cv_chi2_result.num_bins_ << " | " << std::setw(15) <<  cv_chi2_result.p_value_ << std::endl;
	  }
	  genie_cv_truth_smeared_slice->SetTitle(Form("%4.4f",cv_chi2_result.chi2_));
	  
	  //====================================================================================== 
	  //Get the FDS Predictions 
	  
	  TH1* FDS_slice;
	  
	  if ( using_fake_data ) {
	    SliceHistogram* slice_FDS = SliceHistogram::make_slice_histogram( *FDS, Slice, nullptr );
	    FDS_slice = slice_FDS->hist_.get();
	    const auto& FDS_chi2_result = Slice_unf->get_chi2( *slice_FDS );
	    
	    //Only print to console once as Chi2 for Events units and Xsec units are the same
	    if (iRT == 0) {
	      std::cout << std::setw(15) << uc_name << " | " << std::setw(15) << ST << " | " << std::setw(30) << SliceVariableName << " | " << std::setw(15) << "FDS" << " | " << std::setw(15) << FDS_chi2_result.chi2_ << "/" << std::setw(3) << FDS_chi2_result.num_bins_ << " | " << std::setw(15) <<  FDS_chi2_result.p_value_ << std::endl;
	    }
	    
	    FDS_slice->SetTitle(Form("%4.4f",FDS_chi2_result.chi2_));
	  }

	  //====================================================================================== 
	  //Scale by bin width and conversion factor if required
	  
	  Slice_unf_hist->Scale(1.0,"width");
	  genie_cv_truth_smeared_slice->Scale(1.0,"width");
	  
	  int nBins = Slice_unf_hist->GetNbinsX();
	  for (int xBin=0;xBin<nBins;xBin++) {
	    for (int yBin=0;yBin<nBins;yBin++) {
	      SliceCovMat_Mat->operator()(xBin,yBin) = SliceCovMat_Mat->operator()(xBin,yBin)/(Slice_unf_hist->GetXaxis()->GetBinWidth(xBin+1)*Slice_unf_hist->GetXaxis()->GetBinWidth(yBin+1));
	    }
	  }
	  
	  if (RT == "XsecUnits") {
	    Slice_unf_hist->Scale(1.0/conv_factor);
	    genie_cv_truth_smeared_slice->Scale(1.0/conv_factor);
	    SliceCovMat_Mat->operator*=(std::pow( 1.0 / conv_factor, 2 ));

	    Slice_unf_hist->GetYaxis()->SetTitle("Differential Cross Section [10^{-38} cm^{2}/Ar/Bin Width]");
	    genie_cv_truth_smeared_slice->GetYaxis()->SetTitle("Differential Cross Section [10^{-38} cm^{2}/Ar/Bin Width]");
	  }

	  if ( using_fake_data ) {
	    FDS_slice->Scale(1.0,"width");
	    
	    if (RT == "XsecUnits") {
	      FDS_slice->Scale(1.0/conv_factor);
	      FDS_slice->GetYaxis()->SetTitle("Differential Cross Section [10^{-38} cm^{2}/Ar/Bin Width]");
	    }
	  }

	  //====================================================================================== 
	  //Write to File
	  
	  SliceCovMat_Mat->Write("UncerCovMat");
	  Slice_unf_hist->Write("UnfoldedData");
	  genie_cv_truth_smeared_slice->Write("CVPred_Smeared");
	  
	  if ( using_fake_data ) {
	    FDS_slice->Write("FakeDataPred_Smeared");
	  }
	}
      }
    }
  }

  File->Close();
}

int main( int argc, char* argv[] ) {

  if ( argc != 4 ) {
    std::cout << "Usage: Unfolder.C XSEC_Config"
	      << " SLICE_Config OUTPUT_FILE";
    return 1;
  }

  std::string XSEC_Config( argv[1] );
  std::string SLICE_Config( argv[2] );
  std::string OutputFile( argv[3] );

  //Take the output directory from the file handed as the expected output
  //Only used for dumping to text or plot, if that option is requested in the hardcoded options at start of file
  std::string OutputDirectory = OutputFile.substr(0, OutputFile.find_last_of("/") + 1);
  std::string OutputFileName = OutputFile.substr(OutputFile.find_last_of("/") + 1);

  Unfolder(XSEC_Config, SLICE_Config, OutputDirectory, OutputFileName);
  return 0;
}
