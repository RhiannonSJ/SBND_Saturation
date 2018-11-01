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
#include "cetlib_except/exception.h"

// LArSoft Includes
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include <TTree.h>

// C++ Includes
#include <iostream>
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
    void reconfigure(fhicl::ParameterSet const& pset);

  private:

    std::string fG4ModuleLabel;     ///< module label for the Geant
    std::string fTruthModuleLabel;  ///< module label for the Geant

    TTree *tree; ///< Tree to hold 45 entries per event for the time, planes and particle types
    unsigned int bTDC; ///< TDC branch
    unsigned int bWire; ///< Wire plane branch
    unsigned int bParticle; ///< Particle branch
    float bMaxCharge; ///< Max charge branch

    std::vector<unsigned int> tdc_diffs; ///< Vector to hold the standard differences in tdcs for the different time intervals
    std::vector<geo::PlaneID::PlaneID_t> wire_planes; ///< Vector of plane definitions
    std::vector<std::string> particle_types; ///< Vector of particle types, MIP, non-MIP, shower
    
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

    tdc_diffs = std::vector<unsigned int>({4, 10, 20, 40, 100});
    wire_planes = std::vector<geo::PlaneID::PlaneID_t>({0, 1, 2});
    particle_types = std::vector<std::string>({"MIP", "non-MIP", "shower"});

    tree = tfs->make<TTree>("tree","Tree to hold saturation information");
    tree->Branch("bTDC", &bTDC, "bTDC/i");
    tree->Branch("bWire", &bWire, "bWire/i");
    tree->Branch("bParticle", &bParticle, "bParticle/i");
    tree->Branch("bMaxCharge", &bMaxCharge, "bMaxCharge/F");

  }

  //-----------------------------------------------------------------------
  void Saturation::reconfigure(fhicl::ParameterSet const& p)
  {
    fG4ModuleLabel    = p.get< std::string >("GeantModuleLabel");
    fTruthModuleLabel = p.get< std::string >("TruthModuleLabel"); 
    
    return;
  }

  //-----------------------------------------------------------------------
  void Saturation::analyze(const art::Event& evt) 
  {
    
    //get the list of particles from this event
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all sim::SimChannels in the event and make sure there are no
    // sim::IDEs with trackID values that are not in the sim::ParticleList
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);

    // Vector to find the maximum charge on any wire in the time intervals
    std::vector< std::vector< std::vector<float> > > max_charges;

    // For the time ranges
    //    0 = 2us, 1 = 5us, 2 = 10us, 3 = 20us, 4 = 50us
    for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
      // Temp vector of vectors of float (wire planes)
      std::vector< std::vector<float> > temp_wire;
      for(unsigned int j = 0; j < wire_planes.size(); ++j){
        // Temporary vector of float (particle type)
        std::vector<float> temp_particle;
        for(unsigned int k = 0; k < particle_types.size(); ++k){
          temp_particle.push_back(-std::numeric_limits<float>::max());
        }
        temp_wire.push_back(temp_particle);
      }
      max_charges.push_back(temp_wire);
    }

    // Loop over all the sim channels to find various charge-based things
    for(size_t sc = 0; sc < sccol.size(); ++sc){

      // Define a temporary charge vector for the current channel, to be compared
      // to the maximum at the end of the loop, initialize with 0
      std::vector< std::vector< std::vector<float> > > charges;

      // For the time ranges
      //    0 = 2us, 1 = 5us, 2 = 10us, 3 = 20us, 4 = 50us
      for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
        // Temp vector of vectors of float (wire planes)
        std::vector< std::vector<float> > zero_wire;
        for(unsigned int j = 0; j < wire_planes.size(); ++j){
          // Temporary vector of float (particle type)
          std::vector<float> zero_particle;
          for(unsigned int k = 0; k < particle_types.size(); ++k){
            zero_particle.push_back(0);
          }
          zero_wire.push_back(zero_particle);
        }
        charges.push_back(zero_wire);
      }

      // Get the plane ID of the current channel
      raw::ChannelID_t channel        = sccol[sc]->Channel();
      std::vector<geo::WireID> wires  = geom->ChannelToWire(channel);
      geo::PlaneID::PlaneID_t planeID = wires.front().planeID().deepestIndex();

      // Get the map between TDCs and IDEs, ordered by increasing TDC
      const auto & tdcidemap    = sccol[sc]->TDCIDEMap();

      // Get the start TDC for the channel to use when finding the relative end tdcs 
      // for the different time intervals
      unsigned int start_tdc    = tdcidemap.front().first;

      // Loop over the map of TDC to IDE
      for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
        // Get the TDC and the vector of IDEs
        const unsigned int tdc = static_cast<unsigned int>(mapitr->first);
        
        // Loop over each of the defined time ranges
        for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
          // If we are outside of the current range, continue
          if(tdc >= start_tdc + tdc_diffs[i]) continue;

          for(unsigned int j = 0; j < wire_planes.size(); ++j){
            // If the channel is not on the current wire, continue
            if(planeID != wire_planes[j]) continue;

            // Loop over IDEs and check there is a track associated with each one
            const std::vector<sim::IDE> idevec = mapitr->second;
            for(size_t iv = 0; iv < idevec.size(); ++iv){
              if(plist.find( idevec[iv].trackID ) == plist.end() && idevec[iv].trackID != sim::NoParticleId){
                mf::LogWarning("Saturation") << idevec[iv].trackID << " is not in particle list";
                continue;
              }
             
              // Get the MCParticle
              const simb::MCParticle *part = plist[idevec[iv].trackID];

              // If we have the particle, check if it is one of the desired types 
              // and then calculate the charge accumulated on the wire
              for(unsigned int k = 0; k < particle_types.size(); ++k){
                // For MIPs, look for pdgcodes = 211, -211 or 13
                if(k == 0){ // MIP, muons or pions
                  if(part->PdgCode() == 211 || part->PdgCode() == -211 || part->PdgCode() == 13)
                    charges[i][j][k] += idevec[iv].numElectrons; 
                }
                if(k == 1){ // non-MIP, stopping protons
                  if(part->PdgCode() == 2212)
                    charges[i][j][k] += idevec[iv].numElectrons; 
                }
                if(k == 2){ // Showers, electrons or photons
                  if(part->PdgCode() == 11 || part->PdgCode() == 22)
                    charges[i][j][k] += idevec[iv].numElectrons; 
                }
              }
            }
          }
        }
      }
      for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
        for(unsigned int j = 0; j < wire_planes.size(); ++j){
          for(unsigned int k = 0; k < particle_types.size(); ++k){
            if(charges[i][j][k] > max_charges[i][j][k])
              max_charges[i][j][k] = charges[i][j][k];
          }
        }
      }
    } // SimChannels
    for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
      for(unsigned int j = 0; j < wire_planes.size(); ++j){
        for(unsigned int k = 0; k < particle_types.size(); ++k){
          if(max_charges[i][j][k] > 0){
            bTDC = tdc_diffs[i];
            bWire = wire_planes[j];
            bParticle = k;
            bMaxCharge = max_charges[i][j][k];
            tree->Fill();
          }
        }
      }
    }
  } // analyze
} // namespace larg4

namespace larg4 {

  DEFINE_ART_MODULE(Saturation)

} // namespace larg4

#endif // SATURATION_H

