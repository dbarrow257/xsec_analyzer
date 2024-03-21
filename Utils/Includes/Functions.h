#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include "Constants.h"
#include "EventCategory.hh"
#include "AnalysisEvent.h"
#include "FiducialVolume.hh"

#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TH2D.h"
#include "TMatrixD.h"

// Helper function that avoids NaNs when taking square roots of negative
// numbers
inline double real_sqrt( double x ) {
  if ( x < 0. ) return 0.;
  else return std::sqrt( x );
}

// Helper function that returns true if a given PDG code represents a meson or
// antimeson. Otherwise returns false. Based on points 10, 12, and 13 of the
// Particle Data Group's "Monte Carlo Particle Numbering Scheme"
// (2019 revision).

inline bool is_meson_or_antimeson( int pdg_code ) {
  // Ignore differences between mesons and antimesons for this test. Mesons
  // will have positive PDG codes, while antimesons will have negative ones.
  int abs_pdg = std::abs( pdg_code );

  // Meson PDG codes have no more than seven digits. Seven-digit
  // codes beginning with "99" are reserved for generator-specific
  // particles
  if ( abs_pdg >= 9900000 ) return false;

  // Mesons have a value of zero for $n_{q1}$, the thousands digit
  int thousands_digit = ( abs_pdg / 1000 ) % 10;
  if ( thousands_digit != 0 ) return false;

  // They also have a nonzero value for $n_{q2}$, the hundreds digit
  int hundreds_digit = ( abs_pdg / 100 ) % 10;
  if ( hundreds_digit == 0 ) return false;

  // Reserved codes for Standard Model parton distribution functions
  if ( abs_pdg >= 901 && abs_pdg <= 930 ) return false;

  // Reggeon and pomeron
  if ( abs_pdg == 110 || abs_pdg == 990 ) return false;

  // Reserved codes for GEANT tracking purposes
  if ( abs_pdg == 998 || abs_pdg == 999 ) return false;

  // Reserved code for generator-specific pseudoparticles
  if ( abs_pdg == 100 ) return false;

  // If we've passed all of the tests above, then the particle is a meson
  return true;
}

// Function that defines the track-length-dependent proton PID cut
inline double proton_pid_cut( double track_length ) {

  double cut = DEFAULT_PROTON_PID_CUT;

  // Piecewise cut removed 27 June 2021
  //// All track length values are in cm
  //if ( track_length >= 0. && track_length <= 10.5 ) {
  //  cut = -0.0034219305*std::pow( track_length, 2 );
  //  cut += 0.018436866*track_length + 0.062718401;
  //}
  //else if ( track_length > 10.5 && track_length <= 33.1776508 ) {
  //  cut = 0.014153245*( track_length - 10.5 ) - 0.12096235;
  //}
  
  return cut;
}

// Helper function for computing STVs (either reco or true)
inline void compute_stvs( const TVector3& p3mu, const TVector3& p3p, double& delta_pT,
  double& delta_phiT, double& delta_alphaT, double& delta_pL, double& pn,
  double& delta_pTx, double& delta_pTy )
{
  delta_pT = (p3mu + p3p).Perp();

  delta_phiT = std::acos( (-p3mu.X()*p3p.X() - p3mu.Y()*p3p.Y())
    / (p3mu.XYvector().Mod() * p3p.XYvector().Mod()) );

  TVector2 delta_pT_vec = (p3mu + p3p).XYvector();
  delta_alphaT = std::acos( (-p3mu.X()*delta_pT_vec.X()
    - p3mu.Y()*delta_pT_vec.Y())
    / (p3mu.XYvector().Mod() * delta_pT_vec.Mod()) );

  double Emu = std::sqrt(std::pow(MUON_MASS, 2) + p3mu.Mag2());
  double Ep = std::sqrt(std::pow(PROTON_MASS, 2) + p3p.Mag2());
  double R = TARGET_MASS + p3mu.Z() + p3p.Z() - Emu - Ep;

  // Estimated mass of the final remnant nucleus (CCQE assumption)
  double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;
  delta_pL = 0.5*R - (std::pow(mf, 2) + std::pow(delta_pT, 2)) / (2.*R);

  pn = std::sqrt( std::pow(delta_pL, 2) + std::pow(delta_pT, 2) );

  // Components of the 2D delta_pT vector (see arXiv:1910.08658)
  
  // We assume that the neutrino travels along the +z direction (also done
  // in the other expressions above)
  TVector3 zUnit( 0., 0., 1. );

  // Defines the x direction for the components of the delta_pT vector
  TVector2 xTUnit = zUnit.Cross( p3mu ).XYvector().Unit();

  delta_pTx = xTUnit.X()*delta_pT_vec.X() + xTUnit.Y()*delta_pT_vec.Y();

  // Defines the y direction for the components of the delta_T vector
  TVector2 yTUnit = ( -p3mu ).XYvector().Unit();

  delta_pTy = yTUnit.X()*delta_pT_vec.X() + yTUnit.Y()*delta_pT_vec.Y();
}

inline void TM2TH2(const TMatrixD Matrix, TH2D* Hist) {
  int nBins = Matrix.GetNrows();
  if (Hist->GetNbinsX() != nBins) {
    std::cerr << "Matrix and given histogram have different binning" << std::endl;
    throw;
  }  

  for (int xBin=0;xBin<nBins;xBin++) {
    for (int yBin=0;yBin<nBins;yBin++) {
      Hist->SetBinContent(xBin+1,yBin+1,(Matrix)(yBin,xBin));
    }
  }
}

inline void TV2TH(const TVectorD vec, TH1D* histo) {
  // Fill vector to histogram,
  for(Int_t i=0; i<vec.GetNrows(); i++)
    {
      histo->SetBinContent(i+1, vec(i));
    }
}

inline void TH2TM(const TH2D* histo, TMatrixD& mat, bool rowcolumn) {
  // Fill 2D histogram into matrix
  // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE

  for (Int_t i=0; i<histo->GetNbinsX(); i++) {

    for (Int_t j=0; j<histo->GetNbinsY(); j++) {

      if (rowcolumn) { mat(i, j) = histo->GetBinContent(i+1, j+1); }
      else { mat(j, i) = histo->GetBinContent(i+1, j+1); }

    }

  }

}

inline void TH2TV(const TH1D* histo, TVectorD& vec)
{
  // Fill 1D histogram into matrix
  for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
      vec(i) = histo->GetBinContent(i+1);
    }
}

inline void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval, double &sigma) {

  // Clone them so we can scale them 
  TH1D* h_model_clone = (TH1D*)h_model->Clone();
  TH1D* h_data_clone  = (TH1D*)h_data->Clone();
  TH2D* h_cov_clone   = (TH2D*)cov->Clone();
  int NBins = h_cov_clone->GetNbinsX();

  // Getting covariance matrix in TMatrix form
  TMatrixD cov_m;
  cov_m.Clear();
  cov_m.ResizeTo(NBins,NBins);

  // loop over rows

  for (int i = 0; i < NBins; i++) {

    // loop over columns

    for (int j = 0; j < NBins; j++) {

      cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
 
    }
    
  }

  TMatrixD copy_cov_m = cov_m;

  // Inverting the covariance matrix
  TMatrixD inverse_cov_m = cov_m.Invert();

  // Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
  // x = data, mu = model, E^(-1) = inverted covariance matrix 
  chi = 0.;
  
  for (int i = 0; i < NBins; i++) {
    //double XWidth = h_data_clone->GetBinWidth(i+1);
    for (int j = 0; j < NBins; j++) {
      //double YWidth = h_data_clone->GetBinWidth(i+1);
      double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
      double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
      double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
      chi += LocalChi;
    }

  }

  ndof = h_data_clone->GetNbinsX();
  pval = TMath::Prob(chi, ndof);
  sigma = TMath::Sqrt( TMath::ChisquareQuantile( 1-pval, 1 ) ); 

  delete h_model_clone;
  delete h_data_clone;
  delete h_cov_clone;
}

inline TH1D* Multiply(TH1D* True, TH2D* SmearMatrix) {
  TH1D* TrueClone = (TH1D*)(True->Clone());

  int XBins = SmearMatrix->GetXaxis()->GetNbins();
  int YBins = SmearMatrix->GetYaxis()->GetNbins();

  if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

  TVectorD signal(XBins);
  TMatrixD response(XBins,XBins);

  TH2TV(TrueClone, signal);
  TH2TM(SmearMatrix, response, kTRUE);

  TVectorD RecoSpace = response * signal;
  TV2TH(RecoSpace, TrueClone);

  return TrueClone;
}

#endif
