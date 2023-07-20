#ifndef FILE_LOOPER_HH_
#define FILE_LOOPER_HH_

// C++
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <stdexcept>
#include <regex>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiE4D.h>
#include <Math/PxPyPzE4D.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

// Plugins
#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"
#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"

const double E_MASS  = 0.0005109989; //GeV
const double MU_MASS = 0.1056583715; //GeV

class FileLooper {
	/* Class for processing data in a ROOT file event by event */

private:
	// Names
	using LorentzVectorPEP = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>;
	using LorentzVector    = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all, _use_deep_csv, _inc_other_regions, _inc_all_jets, _inc_data, _only_kl1, _only_sm_vbf;
    std::set<std::string> _requested;
    unsigned int _n_feats;
    std::vector<std::string> _feat_names;
    EvtProc* _evt_proc;

	// Methods
    inline int _get_split(const unsigned long int&);
    void _prep_file(TTree* tree, const std::vector<std::unique_ptr<float>>& feat_vals, double* weight, int* sample, int* region, int* jet_cat,
                    int* class_id, unsigned long long int* strat_key);
    Channel _get_channel(std::string);
    Year _get_year(std::string);
    unsigned long long int _get_strat_key(const int& sample, const int& jet_cat, const Channel& channel, const Year& year, const int& region);
    std::vector<std::string> _get_evt_names(const std::map<unsigned long, std::string>&, const std::vector<unsigned long>&);
    int _jet_cat_lookup(const bool has_b_pair, const bool has_vbf_pair, const bool is_boosted, const int num_btag_loose, const int num_btag_medium);
    int _region_lookup(const std::string&);
    void _sample_lookup(std::string& sample, int& sample_id, Spin& spin, float& res_mass);
    int _sample2class_lookup(const int&);
    bool _accept_evt(const int& region, const int& jet_cat, const int& class_id, const float& klambda,
                     const float& cv, const float& c2v, const float& c3);
    bool _apply_baseline(std::string channel, int c_event, int pairType, int nleps, int nbjetscand, int isLeptrigger);
    int _get_region(std::string channel, int isOS, float dau1_deepTauVSJet, float dau2_deepTauVSJet, float dau1_iso, float dau1_eleMVAiso ); 
public:
    // Methods
	FileLooper(bool return_all=true, std::vector<std::string> requested={}, bool use_deep_bjet_wps=true,
               bool inc_all_jets=true, bool inc_other_regions=false, bool inc_data=false,
               bool only_kl1=true, bool only_sm_vbf=true);
	~FileLooper();
    bool loop_file(const std::string& fname, const std::string& out_dir, const std::string& channel,
                   const std::string& year, std::string &sample_str, double sum_w, const long int& n_events,
                   const long int& start_evt, const long int& end_evt);
    std::map<unsigned, std::string> build_dataset_id_map(TFile* in_file);
    std::map<unsigned, std::string> build_region_id_map(TFile* in_file);
};

#endif /* FILE_LOOPER_HH_ */
