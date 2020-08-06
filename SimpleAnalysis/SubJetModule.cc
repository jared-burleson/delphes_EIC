//Module focused on studying SubJet elements from EIC ROOT events
//Currently very simple and under-developed, this is merely the beginnings of the work to study SubJet variables

#include "SubJetModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"

#include "AnalysisFunctions.cc" //Contains the reference to functions to produce SubJet variables
#include "TreeHandler.h"

SubJetModule::SubJetModule(ExRootTreeReader* data) : Module(data)
{

}

SubJetModule::~SubJetModule()
{

}

void SubJetModule::initialize()
{
  TreeHandler *tree_handler = tree_handler->getInstance();

  //std::cout << "Inside SubJetModule::intialize()" << std::endl;

  if(tree_handler->getTree() != nullptr)
  {
    _jet_P4 = std::vector<std::vector<float>>();
    _jet_flavor = std::vector<int>();
    _jet_sIP3D = std::vector<std::vector<float>>();
    _jet_TaggedsIP3D = std::vector<int>();
    
    _subjet1_P4 = std::vector<std::vector<float>>();
    _subjet1_sIP3D = std::vector<std::vector<float>>();
    _subjet1_TaggedsIP3D = std::vector<std::vector<int>>();

    _subjet2_P4 = std::vector<std::vector<float>>();
    _subjet2_sIP3D = std::vector<std::vector<float>>();
    _subjet2_TaggedsIP3D = std::vector<std::vector<int>>();

    _subjet3_P4 = std::vector<std::vector<float>>();
    _subjet3_sIP3D = std::vector<std::vector<float>>();
    _subjet3_TaggedsIP3D = std::vector<std::vector<int>>();

    _subjet4_P4 = std::vector<std::vector<float>>();
    _subjet4_sIP3D = std::vector<std::vector<float>>();
    _subjet4_TaggedsIP3D = std::vector<std::vector<int>>();

    _individual_jet_P4 = std::vector<float>();
    _individual_jet_sIP3D = std::vector<float>();
    _whole_jet_subjet1_P4 = std::vector<float>();
    _whole_jet_subjet2_P4 = std::vector<float>();
    _whole_jet_subjet3_P4 = std::vector<float>();
    _whole_jet_subjet4_P4 = std::vector<float>();
    _whole_jet_subjet1_sIP3D = std::vector<float>();
    _whole_jet_subjet2_sIP3D = std::vector<float>();
    _whole_jet_subjet3_sIP3D = std::vector<float>();
    _whole_jet_subjet4_sIP3D = std::vector<float>();
    _whole_jet_subjet1_TaggedsIP3D = std::vector<int>();
    _whole_jet_subjet2_TaggedsIP3D = std::vector<int>();
    _whole_jet_subjet3_TaggedsIP3D = std::vector<int>();
    _whole_jet_subjet4_TaggedsIP3D = std::vector<int>();
  
    tree_handler->getTree()->Branch("jet.P4", "std::vector<std::vector<float>>", &_jet_P4);
    tree_handler->getTree()->Branch("jet.flavor", "std::vector<int>", &_jet_flavor);
    tree_handler->getTree()->Branch("jet.sIP3D", "std::vector<std::vector<float>>", &_jet_sIP3D);
    tree_handler->getTree()->Branch("jet.TaggedsIP3D", "std::vector<int>", &_jet_TaggedsIP3D);
    tree_handler->getTree()->Branch("subjet1.P4", "std::vector<std::vector<float>>", &_subjet1_P4);
    tree_handler->getTree()->Branch("subjet1.sIP3D","std::vector<std::vector<float>>", &_subjet1_sIP3D);
    tree_handler->getTree()->Branch("subjet1.TaggedsIP3D", "std::vector<std::vector<int>>", &_subjet1_TaggedsIP3D);
    tree_handler->getTree()->Branch("subjet2.P4", "std::vector<std::vector<float>>", &_subjet2_P4);
    tree_handler->getTree()->Branch("subjet2.sIP3D","std::vector<std::vector<float>>", &_subjet2_sIP3D);
    tree_handler->getTree()->Branch("subjet2.TaggedsIP3D", "std::vector<std::vector<int>>", &_subjet2_TaggedsIP3D);
    tree_handler->getTree()->Branch("subjet3.P4", "std::vector<std::vector<float>>", &_subjet3_P4);
    tree_handler->getTree()->Branch("subjet3.sIP3D","std::vector<std::vector<float>>", &_subjet3_sIP3D);
    tree_handler->getTree()->Branch("subjet3.TaggedsIP3D", "std::vector<std::vector<int>>", &_subjet3_TaggedsIP3D);
    tree_handler->getTree()->Branch("subjet4.P4", "std::vector<std::vector<float>>", &_subjet4_P4);
    tree_handler->getTree()->Branch("subjet4.sIP3D","std::vector<std::vector<float>>", &_subjet4_sIP3D);
    tree_handler->getTree()->Branch("subjet4.TaggedsIP3D", "std::vector<std::vector<int>>", &_subjet4_TaggedsIP3D);
  }
  else
  {
    std::cout << "ERROR in SubJetModule::intialize()! tree_handler->getTree() returns nullptr" << std:: endl;
  }
}

void SubJetModule::finalize()
{
  //Empty for the moment. Not sure what should go here yet
}

bool SubJetModule::execute(std::map<std::string, std::any>* DataScore)
{
  //Clearing all memory from previous execute calls

  _jet_P4.resize(0);
  _jet_P4.shrink_to_fit();
  _jet_flavor.resize(0);
  _jet_flavor.shrink_to_fit();
  _jet_sIP3D.resize(0);
  _jet_sIP3D.shrink_to_fit();
  _jet_TaggedsIP3D.resize(0);
  _jet_TaggedsIP3D.shrink_to_fit();

  _subjet1_P4.resize(0);
  _subjet1_P4.shrink_to_fit();
  _subjet1_sIP3D.resize(0);
  _subjet1_sIP3D.shrink_to_fit();
  _subjet1_TaggedsIP3D.resize(0);
  _subjet1_TaggedsIP3D.shrink_to_fit();
  
  _subjet2_P4.resize(0);
  _subjet2_P4.shrink_to_fit();
  _subjet2_sIP3D.resize(0);
  _subjet2_sIP3D.shrink_to_fit();
  _subjet2_TaggedsIP3D.resize(0);
  _subjet2_TaggedsIP3D.shrink_to_fit();
  
  _subjet3_P4.resize(0);
  _subjet3_P4.shrink_to_fit();
  _subjet3_sIP3D.resize(0);
  _subjet3_sIP3D.shrink_to_fit();
  _subjet3_TaggedsIP3D.resize(0);
  _subjet3_TaggedsIP3D.shrink_to_fit();
  
  _subjet4_P4.resize(0);
  _subjet4_P4.shrink_to_fit();
  _subjet4_sIP3D.resize(0);
  _subjet4_sIP3D.shrink_to_fit();
  _subjet4_TaggedsIP3D.resize(0);
  _subjet4_TaggedsIP3D.shrink_to_fit();

  //Populating all vectors and variables with event-information
  
  //Getting all Tracks
  auto tracks = getTracks();
  std::vector<Track*> all_tracks;
  for (auto obj_track : *tracks) 
  {
    auto track = static_cast<Track*>(obj_track);
    all_tracks.push_back( track );
  }

  //Getting all Jets
  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++)
  {
    //Take each individual jet in Event
    Jet *jet = (Jet*) getJets()->At(ijet);
    all_jets.push_back(jet);

    //Clearing all memory from previous jet calls in loop
    _individual_jet_P4.resize(0);
    _individual_jet_P4.shrink_to_fit();
    _individual_jet_sIP3D.resize(0);
    _individual_jet_sIP3D.shrink_to_fit();
    
    _whole_jet_subjet1_P4.resize(0);
    _whole_jet_subjet1_P4.shrink_to_fit();
    _whole_jet_subjet1_sIP3D.resize(0);
    _whole_jet_subjet1_sIP3D.shrink_to_fit();
    _whole_jet_subjet1_TaggedsIP3D.resize(0);
    _whole_jet_subjet1_TaggedsIP3D.shrink_to_fit();
  
    _whole_jet_subjet2_P4.resize(0);
    _whole_jet_subjet2_P4.shrink_to_fit();
    _whole_jet_subjet2_sIP3D.resize(0);
    _whole_jet_subjet2_sIP3D.shrink_to_fit();
    _whole_jet_subjet2_TaggedsIP3D.resize(0);
    _whole_jet_subjet2_TaggedsIP3D.shrink_to_fit();
  
    _whole_jet_subjet3_P4.resize(0);
    _whole_jet_subjet3_P4.shrink_to_fit();
    _whole_jet_subjet3_sIP3D.resize(0);
    _whole_jet_subjet3_sIP3D.shrink_to_fit();
    _whole_jet_subjet3_TaggedsIP3D.resize(0);
    _whole_jet_subjet3_TaggedsIP3D.shrink_to_fit();
  
    _whole_jet_subjet4_P4.resize(0);
    _whole_jet_subjet4_P4.shrink_to_fit();
    _whole_jet_subjet4_sIP3D.resize(0);
    _whole_jet_subjet4_sIP3D.shrink_to_fit();
    _whole_jet_subjet4_TaggedsIP3D.resize(0);
    _whole_jet_subjet4_TaggedsIP3D.shrink_to_fit();

    //Filling all variables/vectors with Jet Data (and SubJet Data)
 
    //Jet Variables
    _individual_jet_P4.push_back( jet->PT );
    _individual_jet_P4.push_back( jet->Eta );
    _individual_jet_P4.push_back( jet->Phi );
    _individual_jet_P4.push_back( jet->Mass );

    _jet_P4.push_back( _individual_jet_P4 );
    _jet_flavor.push_back( jet->Flavor );

    //Tagging Parameters
    float minSignif_jet    = 3.00;
    float minSignif_subjet = 3.00;
    float minPT_jet        = 0.50;
    float minPT_subjet     = 0.50;
    int minTracks_jet      = 2;
    int minTracks_subjet   = 2;
    float delRval_jet      = 0.5;
    float delRval_subjet   = 0.1;
    
    //Jet sIP3D Tagging
    _jet_TaggedsIP3D.push_back(Tagged_sIP3D(jet, *tracks, minSignif_jet, minPT_jet, minTracks_jet, delRval_jet));

    for(int i = 0; i < tracks->GetEntries(); i++)
    {
      auto constituent = tracks->At(i);
      if (constituent->IsA() == Track::Class()) 
      {
        auto track = static_cast<Track*>(constituent);
        const TLorentzVector &trkMomentum = track->P4();
          
        if(trkMomentum.Pt() > minPT_jet && trkMomentum.DeltaR(jet->P4()) < delRval_jet && TMath::Abs(track->D0) < 3.0)
        {
          _individual_jet_sIP3D.push_back(sIP3D(jet, track));
        }
      }
    }
    
    _jet_sIP3D.push_back( _individual_jet_sIP3D );

    //SubJet Section
    auto subjets_trimmed = jet->SubJetsTrimmed;

    //NOTE: This code will ONLY work if subjets_trimmed has a maximum 4 subjets of the jet
    for(int i = 0; i < subjets_trimmed.size(); i++)
    {
      Jet individual_subjet = subjets_trimmed[i];

      if(i == 0)
      {
        _whole_jet_subjet1_P4.push_back( individual_subjet.PT );
        _whole_jet_subjet1_P4.push_back( individual_subjet.Eta );
        _whole_jet_subjet1_P4.push_back( individual_subjet.Phi );
        _whole_jet_subjet1_P4.push_back( individual_subjet.Mass );
        _whole_jet_subjet1_TaggedsIP3D.push_back( Tagged_sIP3D(&individual_subjet, *tracks, minSignif_subjet, minPT_subjet, minTracks_subjet, delRval_subjet) );
        
        for(int i = 0; i < tracks->GetEntries(); i++)                                                                                                                                                                      {                                                                                                                                                                                                                    auto constituent = tracks->At(i);                                                                                                                                                                                  if (constituent->IsA() == Track::Class())
          {
            auto track = static_cast<Track*>(constituent);
            const TLorentzVector &trkMomentum = track->P4();

            if(trkMomentum.Pt() > minPT_subjet && trkMomentum.DeltaR(individual_subjet.P4()) < delRval_subjet && TMath::Abs(track->D0) < 3.0)
            {
              _whole_jet_subjet1_sIP3D.push_back( sIP3D(&individual_subjet, track) );
            }
          }
        }
      }
      else if(i == 1)
      {
        _whole_jet_subjet2_P4.push_back( individual_subjet.PT );
        _whole_jet_subjet2_P4.push_back( individual_subjet.Eta );
        _whole_jet_subjet2_P4.push_back( individual_subjet.Phi );
        _whole_jet_subjet2_P4.push_back( individual_subjet.Mass );
        _whole_jet_subjet2_TaggedsIP3D.push_back( Tagged_sIP3D(&individual_subjet, *tracks, minSignif_subjet, minPT_subjet, minTracks_subjet, delRval_subjet) );  
        
        for(int i = 0; i < tracks->GetEntries(); i++)                                                                                                                                                                      {                                                                                                                                                                                                                    auto constituent = tracks->At(i);                                                                                                                                                                                  if (constituent->IsA() == Track::Class())
          {
            auto track = static_cast<Track*>(constituent);
            const TLorentzVector &trkMomentum = track->P4();

            if(trkMomentum.Pt() > minPT_subjet && trkMomentum.DeltaR(individual_subjet.P4()) < delRval_subjet && TMath::Abs(track->D0) < 3.0)
            {
              _whole_jet_subjet2_sIP3D.push_back( sIP3D(&individual_subjet, track) );
            }
          }
        }
      }
      else if(i == 2)
      {
        _whole_jet_subjet3_P4.push_back( individual_subjet.PT );
        _whole_jet_subjet3_P4.push_back( individual_subjet.Eta );
        _whole_jet_subjet3_P4.push_back( individual_subjet.Phi );
        _whole_jet_subjet3_P4.push_back( individual_subjet.Mass );
        _whole_jet_subjet3_TaggedsIP3D.push_back( Tagged_sIP3D(&individual_subjet, *tracks, minSignif_subjet, minPT_subjet, minTracks_subjet, delRval_subjet) );
        
        for(int i = 0; i < tracks->GetEntries(); i++)                                                                                                                                                                      {                                                                                                                                                                                                                    auto constituent = tracks->At(i);                                                                                                                                                                                  if (constituent->IsA() == Track::Class())
          {
            auto track = static_cast<Track*>(constituent);
            const TLorentzVector &trkMomentum = track->P4();

            if(trkMomentum.Pt() > minPT_subjet && trkMomentum.DeltaR(individual_subjet.P4()) < delRval_subjet && TMath::Abs(track->D0) < 3.0)
            {
              _whole_jet_subjet3_sIP3D.push_back( sIP3D(&individual_subjet, track) );
            }
          }
        }
      }
      else if(i == 3)
      {
        _whole_jet_subjet4_P4.push_back( individual_subjet.PT );
        _whole_jet_subjet4_P4.push_back( individual_subjet.Eta );
        _whole_jet_subjet4_P4.push_back( individual_subjet.Phi );
        _whole_jet_subjet4_P4.push_back( individual_subjet.Mass );
        _whole_jet_subjet4_TaggedsIP3D.push_back( Tagged_sIP3D(&individual_subjet, *tracks, minSignif_subjet, minPT_subjet, minTracks_subjet, delRval_subjet) );    
        
        for(int i = 0; i < tracks->GetEntries(); i++)                                                                                                                                                                      {                                                                                                                                                                                                                    auto constituent = tracks->At(i);                                                                                                                                                                                  if (constituent->IsA() == Track::Class())
          {
            auto track = static_cast<Track*>(constituent);
            const TLorentzVector &trkMomentum = track->P4();

            if(trkMomentum.Pt() > minPT_subjet && trkMomentum.DeltaR(individual_subjet.P4()) < delRval_subjet && TMath::Abs(track->D0) < 3.0)
            {
              _whole_jet_subjet4_sIP3D.push_back( sIP3D(&individual_subjet, track) );
            }
          }
        }
      }
    }

    _subjet1_P4.push_back( _whole_jet_subjet1_P4 );
    _subjet2_P4.push_back( _whole_jet_subjet2_P4 );
    _subjet3_P4.push_back( _whole_jet_subjet3_P4 );
    _subjet4_P4.push_back( _whole_jet_subjet4_P4 );
    _subjet1_sIP3D.push_back( _whole_jet_subjet1_sIP3D );
    _subjet2_sIP3D.push_back( _whole_jet_subjet2_sIP3D );
    _subjet3_sIP3D.push_back( _whole_jet_subjet3_sIP3D );
    _subjet4_sIP3D.push_back( _whole_jet_subjet4_sIP3D );
    _subjet1_TaggedsIP3D.push_back( _whole_jet_subjet1_TaggedsIP3D );
    _subjet2_TaggedsIP3D.push_back( _whole_jet_subjet2_TaggedsIP3D );
    _subjet3_TaggedsIP3D.push_back( _whole_jet_subjet3_TaggedsIP3D );
    _subjet4_TaggedsIP3D.push_back( _whole_jet_subjet4_TaggedsIP3D );

  }
 
  return true;
}
