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

  std::cout << "\nRunning Validation.C with options:" << std::endl;
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

  std::cout << "\n\nStoring inputs -----------------" << std::endl;

  auto smearcept = syst_calc->get_cv_smearceptance_matrix();
  auto true_signal = syst_calc->get_cv_true_signal();

  TMatrixD* true_signal_tmat = true_signal.get();

  auto meas = syst_calc->get_measured_events();
  const auto& data_signal = meas.reco_signal_;
  const auto& data_covmat = meas.cov_matrix_;

  TMatrixD* data_signal_tmat = data_signal.get();
  TMatrixD* data_covmat_tmat = data_covmat.get();

  std::cout << "\n\nUnfolding -----------------" << std::endl;

  auto xsec = extr->get_unfolded_events();
  double conv_factor = extr->conversion_factor();
  std::cout << "Conversion factor:" << 1.0/conv_factor << std::endl;
  const auto& pred_map = extr->get_prediction_map();

  TMatrixD AC_matrix = *xsec.result_.add_smear_matrix_;
  TMatrixD MCPred_Smeared(AC_matrix, TMatrixD::kMult, *true_signal_tmat );

  TMatrixT<double>* TotalCovarianceMatrix;
  for ( const auto& uc_pair : xsec.unfolded_cov_matrix_map_ ) {
    const auto& uc_name = uc_pair.first;
    const auto& uc_matrix = uc_pair.second;
    if (uc_name == "total") {
      TotalCovarianceMatrix = uc_matrix.get();
    }
  }

  TMatrixD response_matrix = *xsec.result_.response_matrix_;
  TMatrixD true_signal_respmat(response_matrix, TMatrixD::kMult, *true_signal_tmat);

  std::cout << "\n\nChecking results -----------------" << std::endl;

  //====================================================================================== 
  //Loop over all the slices taken from the config and save the unfolded distribution/generator prediction
  
  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {
    auto& Slice = sb.slices_.at( sl_idx );
    //The following line will fall over if several active variables are used per Slice
    auto& SliceVar = sb.slice_vars_.at( Slice.active_var_indices_.front() );
    
    std::string SliceVariableName = SliceVar.name_;
    SliceVariableName.erase(std::remove(SliceVariableName.begin(), SliceVariableName.end(), ' '), SliceVariableName.end());
    
    //====================================================================================== 
    //Save unfolded distribution for each covariance matrix
    
    SliceHistogram* Slice_Data_BeforeUnfolding = SliceHistogram::make_slice_histogram( *data_signal_tmat, Slice, data_covmat_tmat );
    SliceHistogram* Slice_MCPred_BeforeSmearing = SliceHistogram::make_slice_histogram( true_signal_respmat, Slice, nullptr );
    
    const auto& chi2_result_beforesmearing = Slice_Data_BeforeUnfolding->get_chi2( *Slice_MCPred_BeforeSmearing );
    std::cout << "Slice " << sl_idx << ", Before Smearing/Unfolding:- " << SliceVariableName << ": \u03C7\u00b2 = " << chi2_result_beforesmearing.chi2_ << '/' << chi2_result_beforesmearing.num_bins_ << " bin" << std::endl;

    SliceHistogram* Slice_Data_Unfolded = SliceHistogram::make_slice_histogram( *xsec.result_.unfolded_signal_, Slice, TotalCovarianceMatrix );
    SliceHistogram* Slice_MCPred_Smeared = SliceHistogram::make_slice_histogram( MCPred_Smeared, Slice, nullptr );

    const auto& chi2_result_afterunfolding = Slice_Data_Unfolded->get_chi2( *Slice_MCPred_Smeared );
    std::cout << "Slice " << sl_idx << ", After Smearing/Unfolding:- " << SliceVariableName << ": \u03C7\u00b2 = " << chi2_result_afterunfolding.chi2_ << '/' << chi2_result_afterunfolding.num_bins_ << " bin" << std::endl;

    TH1* Hist_Data_BeforeUnfolding = Slice_Data_BeforeUnfolding->hist_.get();
    TH1* Hist_MCPred_BeforeSmearing = Slice_MCPred_BeforeSmearing->hist_.get();
    TH1* Hist_Data_Unfolded = Slice_Data_Unfolded->hist_.get();
    TH1* Hist_MCPred_Smeared = Slice_MCPred_Smeared->hist_.get();

    Hist_MCPred_BeforeSmearing->SetTitle(Form("%4.2f",chi2_result_beforesmearing.chi2_));
    Hist_MCPred_Smeared->SetTitle(Form("%4.2f",chi2_result_afterunfolding.chi2_));

    std::cout << std::setw(10) << "Raw Data" << " | " << std::setw(10) << "Raw MC" << " | " << std::setw(10) << "Unf. Data" << " | " << std::setw(10) << "Smr. MC" << std::endl;
    for (int iBin=1;iBin<=Hist_Data_BeforeUnfolding->GetNbinsX();iBin++) {
      std::cout << std::setw(10) << Hist_Data_BeforeUnfolding->GetBinContent(iBin) << " | " << std::setw(10) << Hist_MCPred_BeforeSmearing->GetBinContent(iBin) << " | " << std::setw(10) << Hist_Data_Unfolded->GetBinContent(iBin) << " | " << Hist_MCPred_Smeared->GetBinContent(iBin) << std::endl;
    }

    File->cd();
    File->mkdir(SliceVariableName.c_str());
    File->cd(SliceVariableName.c_str());

    Hist_Data_BeforeUnfolding->Write("Data_BeforeUnfolding");
    Hist_Data_Unfolded->Write("Data_AfterUnfolding");
    Hist_MCPred_BeforeSmearing->Write("MC_NoSmearing");
    Hist_MCPred_Smeared->Write("MC_Smeared");
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
