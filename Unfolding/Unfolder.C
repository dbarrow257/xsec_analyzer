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

void Unfolder(std::string XSEC_Config, std::string SLICE_Config, std::string OutputDirectory, std::string OutputFileName) {

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

  /*
    ToDo: Add this to output
    syst_->total_bnb_data_pot_;    
  */

  std::cout << "Total POT:" << syst_calc->total_bnb_data_pot_ << std::endl;

  std::cout << "\n\nUnfolding -----------------" << std::endl;

  auto xsec = extr->get_unfolded_events();
  double conv_factor = extr->conversion_factor();
  std::cout << "Conversion factor:" << 1.0/conv_factor << std::endl;
  const auto& pred_map = extr->get_prediction_map();

  //Grab the additional smearing matrix to apply to results:
  TMatrixD AC_matrix = *xsec.result_.add_smear_matrix_;
  TH2D* AC_matrix_TH2 = new TH2D("ACMatrix","",AC_matrix.GetNrows(),0,AC_matrix.GetNrows(),AC_matrix.GetNrows(),0,AC_matrix.GetNrows());
  TM2TH2(AC_matrix, AC_matrix_TH2);

  std::cout << "\n\nSaving results -----------------" << std::endl;

  std::cout << "Output file - " << OutputFileName << std::endl;
  if (DumpToText) std::cout << "\tDumping plots to " << TextExtension << " files" << std::endl;
  if (DumpToPlot) std::cout << "\tDumping plots to " << PlotExtension << " files" << std::endl;
  std::cout << "\n" << std::endl;

  // Make results in both Event Count units and then the XSec units
  std::vector<std::string> ResultTypes(2);
  ResultTypes[0] = "EventCountUnits";
  ResultTypes[1] = "XsecUnits";
  size_t nResultTypes = ResultTypes.size();

  // Make predictions in both AC smeared and truth distributions
  std::vector<std::string> SmearTypes(2);
  SmearTypes[0] = "Truth";
  SmearTypes[1] = "ACSmeared";
  size_t nSmearTypes = SmearTypes.size();

  //Loop over the two ResultTypes (Event Counts and Xsec Units)
  for (size_t iRT=0;iRT<nResultTypes;iRT++) {

    std::string RT = ResultTypes[iRT];
    File->cd();
    File->mkdir((std::string("Outputs/")+RT).c_str());
    File->cd((std::string("Outputs/")+RT).c_str());

    //======================================================================================
    //Loop over all covariance matrices stored in the unfolded result and save them

    File->cd((std::string("Outputs/")+RT).c_str());
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
    temp_unfolding_matrix.Write("UnfoldingMatrix");
    temp_err_prop_matrix.Write("ErrorPropagationMatrix");
    temp_add_smear_matrix.Write("AdditionalSmearingMatrix");

    //DB - Need to figure out if these need to be converted - It's a thing which compares truth and reco
    TMatrixD temp_response_matrix = *xsec.result_.response_matrix_;
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
    //Loop over all the slices taken from the config and save the unfolded distribution/generator prediction

    for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
      auto& Slice = sb.slices_.at( sl_idx );
      //The following line will fall over if several active variables are used per Slice
      auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );
      
      std::string SliceVariableName = SliceVar.name_;
      SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());

      File->cd((std::string("Outputs/")+RT).c_str());      
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

	divide_TH1_by_bin_width(SliceHist,true);
	SliceHist->GetYaxis()->SetTitle("Events/GeV");

	if (RT == "XsecUnits") {
	  SliceHist->Scale(1.0 / conv_factor);
	  SliceHist->GetYaxis()->SetTitle("Differential Cross Section [10^{-38} cm^{2}/Ar/GeV]");
	}
	SliceHist->Write((SliceVariableName+"_"+uc_name).c_str());

	if (DumpToText) dump_text_column_vector( OutputDirectory+"/"+RT+"_vec_table_unfolded_signal_"+uc_name+TextExtension, *xsec.result_.unfolded_signal_ );
	if (DumpToPlot) draw_column_vector( OutputDirectory+"/"+RT+"_vec_table_unfolded_signal_"+uc_name+PlotExtension, *xsec.result_.unfolded_signal_, "Unfolded Signal", "Bin Number", "Cross Section [#times 10^{-38} cm^{2}]");
      }

    }

    File->cd();
    File->mkdir((std::string("Predictions/")+RT).c_str());
    File->cd((std::string("Predictions/")+RT).c_str());
    
    //======================================================================================
    //Loop over all generator predictions and save them to the same output
    
    for (int iST=0;iST<nSmearTypes;iST++) {
      std::string ST = SmearTypes[iST];

      File->cd();
      File->mkdir((std::string("Predictions/")+RT+"/"+SmearTypes[iST]).c_str());
      File->cd((std::string("Predictions/")+RT+"/"+SmearTypes[iST]).c_str());

      for ( const auto& gen_pair : extr->get_prediction_map()) {
	std::string gen_short_name = gen_pair.second->name();
	TMatrixD temp_gen = gen_pair.second->get_prediction();
	if (RT == "XsecUnits") {
	  temp_gen *= (1.0 / conv_factor);
	}
	
	TH1D* temp_gen_hist = Matrix_To_TH1(temp_gen,gen_short_name,"","Events");

	if (ST == "ACSmeared") {
	  multiply_1d_hist_by_matrix(&AC_matrix,temp_gen_hist);
	}

	temp_gen_hist->Write(("Pred_"+gen_short_name+"_"+ST).c_str());
	
	if (DumpToText) dump_text_column_vector( OutputDirectory+"/"+RT+"_vec_table_" + gen_short_name + TextExtension, temp_gen );
	if (DumpToPlot) draw_column_vector( OutputDirectory+"/"+RT+"_vec_table_" + gen_short_name + PlotExtension, temp_gen, (gen_short_name + " Prediction").c_str(), "Bin Number", "Cross Section [#times 10^{-38} cm^{2}]");
      }
    }

  }

  std::cout << "\n\nCalculating Chi2 -----------------" << std::endl;
  TString MatrixNameToUse = "FDS_Uncer";
  std::cout << "\n\tUsing " << MatrixNameToUse << " unfolded covariance matrix to calculate Chi2 metrics" << std::endl;

  bool CalculatedChi2 = false;

  for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
    const auto& uc_name = uc_pair.first;
    const auto& uc_matrix = uc_pair.second;

    if (uc_name == MatrixNameToUse) {
      
      //====================================================================================== 
      //Loop over the different slices
      //for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
      for ( size_t sl_idx = 0u; sl_idx < 1; ++sl_idx ) {
	auto& Slice = sb.slices_.at( sl_idx );
	//The following line will fall over if several active variables are used per Slice
	auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );
	
	std::string SliceVariableName = SliceVar.name_;
	SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());
	
	//====================================================================================== 
	//Loop over the two ResultTypes (Event Counts and Xsec Units)
	//for (size_t iRT=0;iRT<nResultTypes;iRT++) {
	for (size_t iRT=0;iRT<1;iRT++) {
	  std::string RT = ResultTypes[iRT];
	  
	  //====================================================================================== 
	  //Get the unfolded Slice Histogram	  
	  SliceHistogram* Slice_unf = SliceHistogram::make_slice_histogram( *xsec.result_.unfolded_signal_, Slice, uc_matrix.get() );

	  //====================================================================================== 
	  //Get the Predictions 

	  //CV Slice
	  SliceHistogram* slice_cv = SliceHistogram::make_slice_histogram(genie_cv_truth_vec, Slice, nullptr );

	  //Fake data slice
	  SliceHistogram* slice_truth = nullptr;
	  if ( using_fake_data ) {
	    slice_truth = SliceHistogram::make_slice_histogram( fake_data_truth, Slice, nullptr );
	  }

	  // Keys are legend labels, values are SliceHistogram objects containing
	  // true-space predictions from the corresponding generator models
	  auto* slice_gen_map_ptr = new std::map< std::string, SliceHistogram* >();
	  auto& slice_gen_map = *slice_gen_map_ptr;

	  slice_gen_map[ "unfolded data" ] = Slice_unf;
	  if ( using_fake_data ) {
	    slice_gen_map[ "truth" ] = slice_truth;
	  }
	  slice_gen_map[ "MicroBooNE Tune" ] = slice_cv;

	  /*
	  for ( const auto& pair : generator_truth_map ) {
	    const auto& model_name = pair.first;
	    TMatrixD* truth_mat = pair.second;
	    
	    SliceHistogram* temp_slice = SliceHistogram::make_slice_histogram(*truth_mat, Slice, nullptr );
	    slice_gen_map[ model_name ] = temp_slice;
	  }
	  */

	  // Keys are generator legend labels, values are the results of a chi^2
	  // test compared to the unfolded data (or, in the case of the unfolded
	  // data, to the fake data truth)
	  std::map< std::string, SliceHistogram::Chi2Result > chi2_map;
	  for ( const auto& pair : slice_gen_map ) {
	    const auto& name = pair.first;
	    const auto* slice_h = pair.second;


	    // Decide what other slice histogram should be compared to this one,
	    // then calculate chi^2
	    SliceHistogram* other = nullptr;
	    // We don't need to compare the unfolded data to itself, so just skip to
	    // the next SliceHistogram and leave a dummy Chi2Result object in the map
	    if ( name == "unfolded data" ) {
	      chi2_map[ name ] = SliceHistogram::Chi2Result();
	      continue;
	    }
	    // Compare all other distributions to the unfolded data
	    else {
	      other = slice_gen_map.at( "unfolded data" );
	    }

	    // Store the chi^2 results in the map
	    const auto& chi2_result = chi2_map[ name ] = slice_h->get_chi2( *other );


	    std::cout << "Slice " << sl_idx << ", " << name << ": \u03C7\u00b2 = "
		      << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bin";
	    if ( chi2_result.num_bins_ > 1 ) std::cout << 's';
	    std::cout << ", p-value = " << chi2_result.p_value_ << '\n';
	  }

	}
      }

      CalculatedChi2 = true;
    }
  }

  if (!CalculatedChi2) {
    std::cerr << "Could not find matrix:" << MatrixNameToUse << " in unfolded covariance matrix map. Did not calculate Chi2 metrics" << std::endl;
    throw;
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
