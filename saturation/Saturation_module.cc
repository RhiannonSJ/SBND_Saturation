////////////////////////////////////////////////////////////////////////
/// \file  Saturation_module.cc
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
#ifndef SATURATION_H
#define SATURATION_H 

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "cetlib_except/exception.h"

// LArSoft Includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TTree.h>

// C++ Includes
#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cstring>
#include <sys/stat.h>

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

///Geant4 interface 
namespace larg4 {  

  class Saturation : public art::EDAnalyzer{
    public:

      /// Standard constructor and destructor for an FMWK module.
      explicit Saturation(fhicl::ParameterSet const& pset);
      virtual ~Saturation();

      void analyze (const art::Event& evt); 
      void beginJob();
      void endJob();
      void reconfigure(fhicl::ParameterSet const& pset);

    private:

      std::string fG4ModuleLabel;     ///< module label for the Geant
      std::string fCorsikaModuleLabel;  ///< module label for corkisa
      std::string fGenieModuleLabel;  ///< module label for Genie

      TTree *corsika_tree;  ///< Cosmic information 
      TTree *neutrino_tree; ///< Neutrino products information
      TTree *event_tree;    ///< Event information

      // Multi-TTree variables (particle &/or coriska &/or event)
      unsigned int bEventId;    ///< Event ID
      unsigned int bParticleId; ///< ID of the particle to count wires crossed
      unsigned int bWireNumber; ///< wire number of the deposited charge
      unsigned int bWirePlane;  ///< wire plane of the deposited charge
      float bTDCStart;          ///< start tdc of the track/shower

      // Cosmics only
      float bEnergy;            ///< Energy of the particles we are looking at 
      float bAngleZ;            ///< Angle to the neutrino direction of the particles we are looking at
      float bAngleY;            ///< Angle to the vertical of the particles we are looking at
      float bCharge;            ///< Max charge branch

      // Event only (BNB-only maximum charge on a single wire over the entire event)
      unsigned int bInteraction;      ///< ScatterCode of the event
      unsigned int bNeutrinoStartTDC; ///< StartTDC of the neutrino

      // Vectors to carry above variables
      std::vector<geo::PlaneID::PlaneID_t> wire_planes; ///< Vector of plane definitions

      // Algorithm counters and module variables
      unsigned int event, valid_neutrino, prior_cosmic, prior_cosmic_muon;
      std::ofstream file;
  };

} // namespace larg4

namespace larg4 {

  //-----------------------------------------------------------------------
  // Constructor
  Saturation::Saturation(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  Saturation::~Saturation() 
  {
  }

  //-----------------------------------------------------------------------
  void Saturation::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;

    // 1 TDC = 0.5 us
    wire_planes = std::vector<geo::PlaneID::PlaneID_t>({0, 1, 2});

    // Write particle tree too for comparison purposes (overlays)
    event_tree    = tfs->make<TTree>("event_tree","Tree to hold saturation information of events");
    corsika_tree  = tfs->make<TTree>("corsika_tree","Tree to hold saturation information of cosmics");
    neutrino_tree = tfs->make<TTree>("neutrino_tree","Tree to hold saturation information of products of the neutrino");

    event_tree->Branch("bEventId",          &bEventId,          "bEventId/i");
    event_tree->Branch("bNeutrinoStartTDC", &bNeutrinoStartTDC, "bNeutrinoStartTDC/i");
    event_tree->Branch("bInteraction",      &bInteraction,      "bInteraction/i");

    corsika_tree->Branch("bEventId",    &bEventId,    "bEventId/i");
    corsika_tree->Branch("bParticleId", &bParticleId, "bParticleId/i");
    corsika_tree->Branch("bTDCStart",   &bTDCStart,   "bTDCStart/F");
    corsika_tree->Branch("bCharge",     &bCharge,     "bCharge/F");
    corsika_tree->Branch("bWirePlane",  &bWirePlane,  "bWirePlane/i");
    corsika_tree->Branch("bWireNumber", &bWireNumber, "bWireNumber/i");
    corsika_tree->Branch("bEnergy",     &bEnergy,     "bEnergy/F");
    corsika_tree->Branch("bAngleZ",     &bAngleZ,     "bAngleZ/F");
    corsika_tree->Branch("bAngleY",     &bAngleY,     "bAngleY/F");

    neutrino_tree->Branch("bEventId",    &bEventId,    "bEventId/i");
    neutrino_tree->Branch("bParticleId", &bParticleId, "bParticleId/i");
    neutrino_tree->Branch("bTDCStart",   &bTDCStart,   "bTDCStart/F");
    neutrino_tree->Branch("bCharge",     &bCharge,     "bCharge/F");
    neutrino_tree->Branch("bWirePlane",  &bWirePlane,  "bWirePlane/i");
    neutrino_tree->Branch("bWireNumber", &bWireNumber, "bWireNumber/i");

    event             = 0;
    valid_neutrino    = 0;
    prior_cosmic      = 0;
    prior_cosmic_muon = 0;

    file.open("studies.txt");
  }
  //-----------------------------------------------------------------------
  void Saturation::endJob()
  {
    std::cout << " Number of events : " << event << std::endl;
    std::cout << " Events with 1 valid neutrino, and possibility for prior cosmics  : " << valid_neutrino    << std::endl;
    std::cout << " Events with 1 valid neutrino, and at least one prior cosmic      : " << prior_cosmic      << std::endl;
    std::cout << " Events with 1 valid neutrino, and at least one prior cosmic muon : " << prior_cosmic_muon << std::endl;
    file.close();
  }

  //-----------------------------------------------------------------------
  void Saturation::reconfigure(fhicl::ParameterSet const& p)
  {
    fG4ModuleLabel      = p.get< std::string >("GeantModuleLabel");
    fCorsikaModuleLabel = p.get< std::string >("CorsikaModuleLabel"); 
    fGenieModuleLabel   = p.get< std::string >("GenieModuleLabel"); 

    return;
  }

  //-----------------------------------------------------------------------
  void Saturation::analyze(const art::Event& evt) 
  {

    // Increment
    event++;

    //get the list of particles from this event
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all sim::SimChannels in the event and make sure there are no
    // sim::IDEs with trackID values that are not in the sim::ParticleList
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);

    // Get the Genie MCTruth information to access the corresponding MCParticle list
    art::Handle< std::vector< simb::MCTruth > > genHandle;
    evt.getByLabel(fGenieModuleLabel, genHandle);

    // Get the Corsika MCTruth information to access the corresponding MCParticle list
    art::Handle< std::vector< simb::MCTruth > > corHandle;
    evt.getByLabel(fCorsikaModuleLabel, corHandle);

    // Get the GTruth information for the interaction codes
    art::Handle< std::vector< simb::GTruth > > gtHandle;
    evt.getByLabel(fGenieModuleLabel, gtHandle);

    // Check that the size of GTruth is 1 or 0, and if it is 0 return since 
    // we don't have the required information
    // Currently not bothered about pileup
    if(gtHandle->size() == 0 || genHandle->size() == 0 || corHandle->size() == 0) return;
    if(gtHandle->size() != 1){
      mf::LogWarning("Saturation") << "GTruth size is " << gtHandle->size() << ", only looking at events with a single neutrino";
      return;      
    }
    if(genHandle->size() != 1){
      mf::LogWarning("Saturation") << "Genie MCTruth size is " << genHandle->size() << ", only looing at events with a single neutrino";
      return;      
    }
    if(corHandle->size() < 1){
      mf::LogWarning("Saturation") << "Corsika MCTruth size is " << corHandle->size() << ", need at least 1 cosmic for the study";
      return;      
    }

    // Set the scattering code of the event
    art::Ptr< simb::GTruth > gt( gtHandle, 0 );
    unsigned int scatterCode = gt->fGscatter;

    // Get the vectors of Corsika and Genie MCParticles 
    art::FindManyP<simb::MCParticle> fgenpart(genHandle, evt, fG4ModuleLabel);
    art::FindManyP<simb::MCParticle> fcorpart(corHandle, evt, fG4ModuleLabel);

    std::vector< art::Ptr<simb::MCParticle> > genParticles = fgenpart.at(0);
    std::vector< art::Ptr<simb::MCParticle> > corParticles = fcorpart.at(0);

    // Need at least 1 neutrino product in the event
    if(genParticles.size() == 0) return;

    // Vector to hold the smallest TDC's given by products of the neutrino interaction
    // in the event
    // This way we can find the first recorded TDC of the neutrino event and set
    // the cosmic search 400us before that (800TDC)
    std::vector<unsigned short> smallest_neutrino_tdcs;

    file << "-----------------------------------" << std::endl;

    if(!sccol.size()) {
      mf::LogWarning("Saturation") << "No SimChannels";
      return;
    }

    // Loop over all the sim channels to find various charge-based things
    for(size_t sc = 0; sc < sccol.size(); ++sc){

      // Get the plane ID of the current channel
      //raw::ChannelID_t channel        = sccol[sc]->Channel();

      // Get the map between TDCs and IDEs, ordered by increasing TDC
      const auto & tdcidemap    = sccol[sc]->TDCIDEMap();

      // Loop over the map of TDC to IDE
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){

        // Get the TDC and the vector of IDEs
        //const unsigned int tdc = static_cast<unsigned int>(mapitr->first);

        // Looking at the neutrinos
        // Loop over IDEs and check there is a track associated with each one
        const std::vector<sim::IDE> idevec = mapitr->second;
        for(size_t iv = 0; iv < idevec.size(); ++iv){
          if(plist.find(idevec[iv].trackID) == plist.end() && idevec[iv].trackID != sim::NoParticleId){
            mf::LogWarning("Saturation") << idevec[iv].trackID << " is not in particle list";
            continue;
          }

          // Get the current particle from the back tracker
          const simb::MCParticle *part = plist.at(idevec[iv].trackID);

          // Get the MCTruth information from the particle
          art::Ptr<simb::MCTruth> currentTruth = pi_serv->ParticleToMCTruth_P(part);

          // Make sure we are looking at the neutrino
          //if(currentTruth->Origin() != simb::kBeamNeutrino) continue;
          if(currentTruth->Origin() == simb::kBeamNeutrino){

            // Get the start TDC for the channel to use when finding the relative end tdcs 
            // for the different time intervals
            // .front() returns the smallest TDC value for the simchannel 
            // .back() returns the largest TDC value for the simchannel
            smallest_neutrino_tdcs.push_back(mapitr->first);
          }
        }// IDE vector
      }// TDCIDEMap
    }// SimChannels
    // Find the smallest TDC given by the neutrino part of the event
    file << " Event                 : " << event << std::endl;
    if(!smallest_neutrino_tdcs.size()) {
      mf::LogWarning("Saturation") << "There are no neutrino TDCs in the event";
      return;
    }
    unsigned short smallest_neutrino_tdc = *std::min_element(begin(smallest_neutrino_tdcs), end(smallest_neutrino_tdcs));
    if(smallest_neutrino_tdc == 0){
      mf::LogWarning("Saturation") << "Smallest neutrino TDC is 0, so the cosmics cannot have been recorded before this";
      return;
    }
    valid_neutrino++;

    unsigned short valid_cosmic      = 0;
    unsigned short valid_cosmic_muon = 0;

    // If we have found a neutrino with a start TDC greater than 0 then look for the cosmics which have come before it
    for(size_t sc = 0; sc < sccol.size(); ++sc){

      // Get the plane ID of the current channel
      raw::ChannelID_t channel        = sccol[sc]->Channel();
      std::vector<geo::WireID> wires  = geom->ChannelToWire(channel);             
      geo::PlaneID::PlaneID_t planeID = wires.front().planeID().deepestIndex();  
      
      // Only looking at the collection plane
      if(planeID != 2) continue;

      // Get the map between TDCs and IDEs, ordered by increasing TDC
      const auto & tdcidemap    = sccol[sc]->TDCIDEMap();

      // Loop over the map of TDC to IDE
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){

        // Get the TDC and the vector of IDEs
        //const unsigned int tdc = static_cast<unsigned int>(mapitr->first);

        // Looking at the cosmics
        // Loop over IDEs and check there is a track associated with each one
        const std::vector<sim::IDE> idevec = mapitr->second;
        for(size_t iv = 0; iv < idevec.size(); ++iv){
          if(plist.find(idevec[iv].trackID) == plist.end() && idevec[iv].trackID != sim::NoParticleId){
            mf::LogWarning("Saturation") << idevec[iv].trackID << " is not in particle list";
            continue;
          }

          // Get the current particle from the back tracker
          const simb::MCParticle *part = plist.at(idevec[iv].trackID);

          // Get the MCTruth information from the particle
          art::Ptr<simb::MCTruth> currentTruth = pi_serv->ParticleToMCTruth_P(part);

          // Having determined that we are looking at an event with a valid neutrino and some prior cosmics, start writing information to the root output

          if(currentTruth->Origin() == simb::kCosmicRay){
            // If the cosmic TDC is between 
            // (The smallest neutrino tdc - 800) & the smallest neutrino TDC
            // Call it a valid cosmic
            if(mapitr->first < smallest_neutrino_tdc && mapitr->first > (smallest_neutrino_tdc - 2000)){
              valid_cosmic++;
            
              if(part->PdgCode() != 13) continue;
              valid_cosmic_muon++;

              TVector3 direction(part->Px()/part->P(), part->Py()/part->P(), part->Pz()/part->P());
              TVector3 z(0,0,1);
              TVector3 y(0,1,0);
              double costhetaZ = direction.Dot(z);
              double costhetaY = direction.Dot(y);

              bEventId       = event;
              bParticleId    = part->TrackId();
              bCharge        = idevec[iv].numElectrons;
              bEnergy        = idevec[iv].energy;
              bAngleY        = costhetaY;
              bAngleZ        = costhetaZ;
              bWireNumber    = channel;
              bWirePlane     = planeID;
              bTDCStart      = tdcidemap.front().first;
              corsika_tree->Fill();
            } // Valid cosmic
          } // Cosmic rays
        } // IDE vector
      } // TDC - IDE map
    } // SimChannels
    if(valid_cosmic == 0){
      mf::LogWarning("Saturation") << "There are no cosmic TDCs before the neutrino interaction in the event";
      return;
    }
    prior_cosmic++;
    if(valid_cosmic_muon == 0){
      mf::LogWarning("Saturation") << "There are no cosmic muon TDCs before the neutrino interaction in the event";
      return;
    }
    prior_cosmic_muon++;
    
    // If we have found a neutrino with a start TDC greater than 0 then look for the cosmics which have come before it
    for(size_t sc = 0; sc < sccol.size(); ++sc){

      // Get the plane ID of the current channel
      raw::ChannelID_t channel        = sccol[sc]->Channel();
      std::vector<geo::WireID> wires  = geom->ChannelToWire(channel);             
      geo::PlaneID::PlaneID_t planeID = wires.front().planeID().deepestIndex();  

      // Get the map between TDCs and IDEs, ordered by increasing TDC
      const auto & tdcidemap    = sccol[sc]->TDCIDEMap();

      // Loop over the map of TDC to IDE
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){

        // Get the TDC and the vector of IDEs
        //const unsigned int tdc = static_cast<unsigned int>(mapitr->first);

        // Looking at the cosmics
        // Loop over IDEs and check there is a track associated with each one
        const std::vector<sim::IDE> idevec = mapitr->second;
        for(size_t iv = 0; iv < idevec.size(); ++iv){
          if(plist.find(idevec[iv].trackID) == plist.end() && idevec[iv].trackID != sim::NoParticleId){
            mf::LogWarning("Saturation") << idevec[iv].trackID << " is not in particle list";
            continue;
          }

          // Get the current particle from the back tracker
          const simb::MCParticle *part = plist.at(idevec[iv].trackID);

          // Get the MCTruth information from the particle
          art::Ptr<simb::MCTruth> currentTruth = pi_serv->ParticleToMCTruth_P(part);

          // Make sure we are looking at the neutrino now
          if(currentTruth->Origin() != simb::kBeamNeutrino) continue;
          
          // Only look at primary neutrino final state particles 
          // Don't want to look at elements (high pdgcode)
          if(part->Process() != "primary" || part->PdgCode() >= 1000018039) continue;

          // Write to the neutrino tree
          bEventId       = event;
          bParticleId    = part->TrackId();
          bCharge        = idevec[iv].numElectrons;
          bWireNumber    = channel;
          bWirePlane     = planeID;
          bTDCStart      = tdcidemap.front().first;
          neutrino_tree->Fill();
        } // IDE vector
      } // TDC - IDE map
    } // SimChannels
    bEventId          = event;
    bNeutrinoStartTDC = smallest_neutrino_tdc;
    bInteraction      = scatterCode;
    event_tree->Fill();
  } // analyze
} // namespace larg4

namespace larg4 {

  DEFINE_ART_MODULE(Saturation)

} // namespace larg4

#endif // SATURATION_H

