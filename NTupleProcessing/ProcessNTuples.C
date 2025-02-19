// Analysis macro for use in the CCNp0pi single transverse variable analysis
// Designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 22 April 2023
// Steven Gardiner <gardiner@fnal.gov>
//
// Update September 2024
// Daniel Barrow <daniel.barrow@physics.ox.ac.uk>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

//Dan's includes
#include "AnalysisEvent.h"
#include "Constants.h"
#include "Functions.h"
#include "Branches.h"

#include "SelectionBase.h"
#include "SelectionFactory.h"

void analyze(const std::vector<std::string>& in_file_names,
  const std::string& output_filename)
{
  std::cout << "\nRunning ProcessNTuples with options:" << std::endl;
  std::cout << "\toutput_filename: " << output_filename << std::endl;
  std::cout << "\tinput_file_names: " << std::endl;
  for (size_t i=0;i<in_file_names.size();i++) {
    std::cout << "\t\t- " << in_file_names[i] << std::endl;
  }
  std::cout << "\n" << std::endl;
  
  // Get the TTrees containing the event ntuples and subrun POT information
  // Use TChain objects for simplicity in manipulating multiple files
  TChain events_ch( "nuselection/NeutrinoSelectionFilter" );
  TChain subruns_ch( "nuselection/SubRun" );

  for ( const auto& f_name : in_file_names ) {
    events_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
  }
  
  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( output_filename.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
  }

  TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot",
    summed_pot );

  summed_pot_param->Write();

  std::vector<SelectionBase*> Selections;

  SelectionFactory* SelFactory = new SelectionFactory();
  Selections.push_back(SelFactory->CreateSelection("CC1mu1p0pi"));
  Selections.push_back(SelFactory->CreateSelection("CC1mu2p0pi"));
  Selections.push_back(SelFactory->CreateSelection("CC1muNp0pi"));

  out_file->cd();
  for (size_t i=0;i<Selections.size();i++) {
    Selections[i]->Setup(out_tree);
  }

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;

  while ( true ) {
    
    //if ( events_entry > 1000) break;
    
    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }
    
    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }
    set_event_output_branch_addresses(*out_tree, cur_event, create_them );

    for (size_t i=0;i<Selections.size();i++) {
      Selections[i]->ApplySelection(&(cur_event));
    }
    
    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }
  
  for (size_t i=0;i<Selections.size();i++) {
    Selections[i]->Summary();
  }
  std::cout << "Wrote output to:" << output_filename << std::endl;

  for (size_t i=0;i<Selections.size();i++) {
    Selections[i]->FinalTasks();
  }
  
  out_tree->Write();
  out_file->Close();
  delete out_file;
}

void analyzer(const std::string& in_file_name,
 const std::string& output_filename)
{
  std::vector<std::string> in_files = { in_file_name };
  analyze( in_files, output_filename );
}

int main( int argc, char* argv[] ) {

  if ( argc != 3 ) {
    std::cout << "Usage: analyzer INPUT_PELEE_NTUPLE_FILE OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[1] );
  std::string output_file_name( argv[2] );

  analyzer( input_file_name, output_file_name );

  return 0;
}
