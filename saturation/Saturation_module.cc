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
    std::string fTruthModuleLabel;  ///< module label for the Geant

    TTree *particle_tree; ///< Tree to hold 45 entries per event for the time, planes and particle types
    TTree *event_tree; ///< Tree to hold 15 entries per event for the time and planes
    unsigned int bTDC; ///< TDC branch
    unsigned int bWire; ///< Wire plane branch
    unsigned int bParticle; ///< Particle branch
    unsigned int bInteraction; ///< ScatterCode of the event 
    float bMaxCharge; ///< Max charge branch
    float bEnergy; ///< Energy of the particles we are looking at 
    float bAngle; ///< Angle to the neutrino direction of the particles we are looking at
    float bEventCharge; ///< Charge accumulated on a wire over the whole event

    std::vector<unsigned int> tdc_diffs; ///< Vector to hold the standard differences in tdcs for the different time intervals
    std::vector<geo::PlaneID::PlaneID_t> wire_planes; ///< Vector of plane definitions
    std::vector<std::string> particle_types; ///< Vector of particle types, MIP, non-MIP, shower
   
    // Counters for finding interesting events
    unsigned int event;
    unsigned int event_high_mip;
    unsigned int event_low_shower;
    
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

    tdc_diffs = std::vector<unsigned int>({4, 10, 20, 40, 100});
    wire_planes = std::vector<geo::PlaneID::PlaneID_t>({0, 1, 2});
    particle_types = std::vector<std::string>({"MIP", "non-MIP", "shower"});

    particle_tree = tfs->make<TTree>("particle_tree","Tree to hold saturation information of particles");
    event_tree    = tfs->make<TTree>("event_tree","Tree to hold saturation information of events");
    particle_tree->Branch("bTDC", &bTDC, "bTDC/i");
    particle_tree->Branch("bWire", &bWire, "bWire/i");
    particle_tree->Branch("bParticle", &bParticle, "bParticle/i");
    particle_tree->Branch("bMaxCharge", &bMaxCharge, "bMaxCharge/F");
    particle_tree->Branch("bEnergy", &bEnergy, "bEnergy/F");
    particle_tree->Branch("bAngle", &bAngle, "bAngle/F");
    particle_tree->Branch("bInteraction", &bInteraction, "bInteraction/i");
    event_tree->Branch("bEventCharge", &bEventCharge, "bEventCharge/F");
    event_tree->Branch("bEventTDC", &bTDC, "bEventTDC/i");
    event_tree->Branch("bEventWire", &bWire, "bEventWire/i");
    event_tree->Branch("bEventInteraction", &bInteraction, "bEventInteraction/i");

    event            = 0;

    // File to hold interesting events
    file.open("evd.txt");
    file << std::setw(20) << " High charge MIP" << std::setw(20) << " High charge shower " << std::setw(20) << " Low charge shower " << std::endl;
  }
  //-----------------------------------------------------------------------
  void Saturation::endJob()
  {
    file.close();
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
  
    event++;

    //get the list of particles from this event
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    const sim::ParticleList& plist = pi_serv->ParticleList();
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all sim::SimChannels in the event and make sure there are no
    // sim::IDEs with trackID values that are not in the sim::ParticleList
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);

    // Get the GTruth information for the interaction codes
    art::Handle< std::vector< simb::GTruth > > gtHandle;
    evt.getByLabel(fTruthModuleLabel, gtHandle);

    // Check that the size of GTruth is 1 or 0, and if it is 0 return since 
    // we don't have the required information
    if(gtHandle->size() == 0) return;
    if(gtHandle->size() != 1){
      mf::LogWarning("Saturation") << "GTruth size is " << gtHandle->size() << ", should be 0 or 1";
      return;      
    }

    // Set the scattering code of the event
    art::Ptr< simb::GTruth > gt( gtHandle, 0 );
    unsigned int scatterCode = gt->fGscatter;

    // Define the maximum charge on a wire in an event 
    std::vector< std::vector<float> > max_event_charge;

    // Vector to find the maximum charge on any wire in the time intervals
    std::vector< std::vector< std::vector<float> > > max_charges;
    std::vector< std::vector< std::vector<float> > > max_energy;
    std::vector< std::vector< std::vector<float> > > max_angle;

    // For the time ranges
    //    0 = 2us, 1 = 5us, 2 = 10us, 3 = 20us, 4 = 50us
    for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
      // Temp vector of vectors of float (wire planes)
      std::vector< std::vector<float> > temp_wire;
      std::vector<float> temp_event_wire;
      for(unsigned int j = 0; j < wire_planes.size(); ++j){
        // Temporary vector of float (particle type)
        std::vector<float> temp_particle;
        for(unsigned int k = 0; k < particle_types.size(); ++k){
          temp_particle.push_back(-std::numeric_limits<float>::max());
        }
        temp_wire.push_back(temp_particle);
        temp_event_wire.push_back(-std::numeric_limits<float>::max());
      }
      max_charges.push_back(temp_wire);
      max_energy.push_back(temp_wire);
      max_angle.push_back(temp_wire);
      max_event_charge.push_back(temp_event_wire);
    }

    // Loop over all the sim channels to find various charge-based things
    for(size_t sc = 0; sc < sccol.size(); ++sc){

      // Define a temporary charge vector for the current channel, to be compared
      // to the maximum at the end of the loop, initialize with 0
      std::vector< std::vector< std::vector<float> > > charges;
      std::vector< std::vector< std::vector<float> > > energy;
      std::vector< std::vector< std::vector<float> > > angle;

      // The total charge on the current wire is 0 so far
      std::vector< std::vector<float> > event_charge;

      // For the time ranges
      //    0 = 2us, 1 = 5us, 2 = 10us, 3 = 20us, 4 = 50us
      for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
        // Temp vector of vectors of float (wire planes)
        std::vector< std::vector<float> > zero_wire;
        std::vector<float> zero_event_wire;
        for(unsigned int j = 0; j < wire_planes.size(); ++j){
          // Temporary vector of float (particle type)
          std::vector<float> zero_particle;
          for(unsigned int k = 0; k < particle_types.size(); ++k){
            zero_particle.push_back(0);
          }
          zero_wire.push_back(zero_particle);
          zero_event_wire.push_back(0);
        }
        charges.push_back(zero_wire);
        energy.push_back(zero_wire);
        angle.push_back(zero_wire);
        event_charge.push_back(zero_event_wire);
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
              TVector3 direction(part->Px()/part->P(), part->Py()/part->P(), part->Pz()/part->P());
              TVector3 z(0,0,1);
              double costheta = direction.Dot(z);

              // If we have the particle, check if it is one of the desired types 
              // and then calculate the charge accumulated on the wire
              for(unsigned int k = 0; k < particle_types.size(); ++k){
                // For MIPs, look for pdgcodes = 211, -211 or 13
                if(k == 0){ // MIP, muons or pions
                  if(part->PdgCode() == 211 || part->PdgCode() == -211 || part->PdgCode() == 13){
                    charges[i][j][k] += idevec[iv].numElectrons;
                    energy[i][j][k]   = idevec[iv].energy;
                    angle[i][j][k]    = costheta;
                  }
                }
                if(k == 1){ // non-MIP, stopping protons
                  if(part->PdgCode() == 2212){
                    charges[i][j][k] += idevec[iv].numElectrons; 
                    energy[i][j][k]   = idevec[iv].energy;
                    angle[i][j][k]    = costheta;
                  }
                }
                if(k == 2){ // Showers, electrons or photons
                  if(part->PdgCode() == 11 || part->PdgCode() == 22){
                    charges[i][j][k] += idevec[iv].numElectrons; 
                    energy[i][j][k]   = idevec[iv].energy;
                    angle[i][j][k]    = costheta;
                  }
                }
              }
              // Fill the event_charge vector for the current wire plane and tdc tick 
              event_charge[i][j] += idevec[iv].numElectrons;
            }
          }
        }
      }

      // Set the maximum charges and respective energies and angles
      for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
        for(unsigned int j = 0; j < wire_planes.size(); ++j){
          for(unsigned int k = 0; k < particle_types.size(); ++k){
            if(charges[i][j][k] > max_charges[i][j][k]){
              max_charges[i][j][k] = charges[i][j][k];
              max_energy[i][j][k]  = energy[i][j][k];
              max_angle[i][j][k]   = angle[i][j][k];
            }
          }
          if(event_charge[i][j] > max_event_charge[i][j])
            max_event_charge[i][j] = event_charge[i][j];
        }
      }
    } // SimChannels
    for(unsigned int i = 0; i < tdc_diffs.size(); ++i){
      for(unsigned int j = 0; j < wire_planes.size(); ++j){
        if(max_event_charge[i][j] > 0){
          bEventCharge = max_event_charge[i][j];
          bTDC = tdc_diffs[i];
          bWire = wire_planes[j];
          bInteraction = scatterCode;
          event_tree->Fill();
        }
        for(unsigned int k = 0; k < particle_types.size(); ++k){
          if(max_charges[i][j][k] > 0){
            bParticle = k;
            bMaxCharge = max_charges[i][j][k];
            bEnergy = max_energy[i][j][k];
            bAngle  = max_angle[i][j][k];
            particle_tree->Fill();

            // High charge MIP
            if(max_charges[i][j][k] * 1.6e-19 * 1e15 > 250 && k == 0)
              file << std::setw(8) << event << std::setw(8) << max_charges[i][j][k] * 1.6e-19 *1e15 << std::setw(4) << j << std::setw(20) << " " << std::setw(20) << " " << std::endl;
            if(max_charges[i][j][k] * 1.6e-19 * 1e15 >= 100 && k == 2)
              file << std::setw(20) << " " << std::setw(8) << event << std::setw(8) << max_charges[i][j][k] * 1.6e-19 *1e15 << std::setw(4) << j << std::setw(20) << " " << std::endl;
            if(max_charges[i][j][k] * 1.6e-19 * 1e15 > 50 && max_charges[i][j][k] * 1e-19 * 1e15 < 100 && k == 2)
              file << std::setw(20) << " " << std::setw(20) << " " << std::setw(8) << event << std::setw(8) << max_charges[i][j][k] * 1.6e-19 *1e15 << std::setw(4) << j << std::endl;
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

