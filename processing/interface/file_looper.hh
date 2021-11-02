#ifndef FILE_LOOPER_HH_
#define FILE_LOOPER_HH_

// C++
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <stdexcept>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/PxPyPzM4D.h>
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
	using LorentzVectorPEP = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
	using LorentzVector    = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all, _use_deep_csv, _apply_cut, _inc_other_regions, _inc_all_jets, _inc_unc, _inc_data, _only_kl1, _only_sm_vbf;
    std::set<std::string> _requested;
    unsigned int _n_feats;
    std::vector<std::string> _feat_names;
    EvtProc* _evt_proc;

	// Methods
    inline int _get_split(const unsigned long int&);
    void _prep_file(TTree* tree, const std::vector<std::unique_ptr<float>>& feat_vals, double* weight, int* sample, int* region, int* jet_cat,
                    bool* scale, bool* central_unc, int* class_id, unsigned long long int* strat_key, float* mva_score);
    Channel _get_channel(std::string);
    Year _get_year(std::string);
    unsigned long long int _get_strat_key(const int& sample, const int& jet_cat, const Channel& channel, const Year& year, const int& region,
                                          const int& central_unc, const int& cut_pass);
    std::vector<std::string> _get_evt_names(const std::map<unsigned long, std::string>&, const std::vector<unsigned long>&);
    void _extract_flags(const std::vector<std::string>& names, int& sample, int& region, bool& central_unc, bool& scale,
                        int& jet_cat, bool& cut_pass, int& class_id, Spin& spin, float& klambda, float& res_mass,
                        bool& is_boosted, bool& accept, std::vector<unsigned int>& idxs , float& cv, float& c2v, float& c3);
    int _jet_cat_lookup(const std::string&);
    int _region_lookup(const std::string&);
    void _sample_lookup(const std::string& sample, int& sample_id, Spin& spin, float& klambda, float& res_mass, float& cv, float& c2v, float& c3);
    int _sample2class_lookup(const int&);
    bool _accept_evt(const int& region, const bool& central_unc, const int& jet_cat, const bool& cut_pass, const int& class_id, const float& klambda,
                     const float& cv, const float& c2v, const float& c3);
    double _get_weight(TTreeReaderValue<std::vector<double>>& rv_weight, const std::vector<unsigned int>& idxs,
                               const std::vector<std::string>& names);
    float _get_mva_score(TTreeReaderValue<std::vector<float>>& rv_mva_score, const std::vector<unsigned int>& idxs,
                                 const std::vector<std::string>& names);

public:
    // Methods
	FileLooper(bool return_all=true, std::vector<std::string> requested={}, bool use_deep_bjet_wps=true,
               bool apply_cut=false, bool inc_all_jets=true, bool inc_other_regions=false, bool inc_data=false, bool inc_unc=false,
               bool only_kl1=true, bool only_sm_vbf=true);
	~FileLooper();
	bool loop_file(const std::string&, const std::string&, const std::string&, const std::string&, const long int&);
    std::map<unsigned long, std::string> build_id_map(TFile*);
};

#endif /* FILE_LOOPER_HH_ */
