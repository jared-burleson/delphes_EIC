//Header File for first attempt at a SubJet Module for SimpleAnalysis for Jet Studies at EIC

#ifndef SUBJETMODULE_HH
#define SUBJETMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"
#include <vector>
#include <map>
#include <utility>

class SubJetModule : public Module
{
  public:

    SubJetModule(ExRootTreeReader* data); //Constructor

    ~SubJetModule(); //Destructor

    void initialize() override;
    bool execute(std::map<std::string, std::any>* DataScore) override;
    void finalize() override;

  private:

    //Jet-Level Variables, will be duplicated for each subjet that comes from the same jet
    std::vector<std::vector<float>> _jet_P4;
    std::vector<int>                _jet_flavor;
    std::vector<std::vector<float>> _jet_sIP3D;
    std::vector<int>                _jet_TaggedsIP3D;
    
    //SubJet-Level Variables, only using 4 hardest (ranked by pT), so will contain 
    std::vector<std::vector<float>> _subjet1_P4;
    std::vector<std::vector<float>> _subjet1_sIP3D;
    std::vector<std::vector<int>>   _subjet1_TaggedsIP3D;

    std::vector<std::vector<float>> _subjet2_P4;
    std::vector<std::vector<float>> _subjet2_sIP3D;
    std::vector<std::vector<int>>   _subjet2_TaggedsIP3D;

    std::vector<std::vector<float>> _subjet3_P4;
    std::vector<std::vector<float>> _subjet3_sIP3D;
    std::vector<std::vector<int>>   _subjet3_TaggedsIP3D;

    std::vector<std::vector<float>> _subjet4_P4;
    std::vector<std::vector<float>> _subjet4_sIP3D;
    std::vector<std::vector<int>>   _subjet4_TaggedsIP3D;

    //Temporary Vectors that are not outputted but rather used in code to store data for later usage
    std::vector<float> _individual_jet_P4;
    std::vector<float> _individual_jet_sIP3D;
    std::vector<float> _whole_jet_subjet1_P4;
    std::vector<float> _whole_jet_subjet2_P4;
    std::vector<float> _whole_jet_subjet3_P4;
    std::vector<float> _whole_jet_subjet4_P4;                                                                                                                                                                          std::vector<float> _whole_jet_subjet1_sIP3D;                                                                                                                                                                       std::vector<float> _whole_jet_subjet2_sIP3D;                                                                                                                                                                       std::vector<float> _whole_jet_subjet3_sIP3D;
    std::vector<float> _whole_jet_subjet4_sIP3D;
    std::vector<int>   _whole_jet_subjet1_TaggedsIP3D;
    std::vector<int>   _whole_jet_subjet2_TaggedsIP3D;
    std::vector<int>   _whole_jet_subjet3_TaggedsIP3D;
    std::vector<int>   _whole_jet_subjet4_TaggedsIP3D;
};

#endif
