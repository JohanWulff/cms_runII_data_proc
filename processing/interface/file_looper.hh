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
    bool _all, _use_deep_csv, _apply_cut, _inc_other_regions, _inc_all_jets, _inc_unc, _inc_data;
    std::set<std::string> _requested;
    unsigned int _n_feats;
    std::vector<std::string> _feat_names;
    EvtProc* _evt_proc;

	// Methods
    inline int _get_split(const unsigned long int&);
    void _prep_file(TTree*, const std::vector<float>&, const double&, const int&, const int&, const int&, const int&, const bool&, const bool&, const int&,
                    const unsigned long long int&);
    Channel _get_channel(std::string);
    Year _get_year(std::string);
    unsigned long long int _get_strat_key(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&);
    std::vector<std::string> _get_evt_names(const std::map<unsigned long, std::string>&, const std::vector<unsigned long>&);
    void _extract_flags(const std::vector<std::string>&, int&, int&, bool&, bool&, int&, int&, int&, Spin&, float&, float&, bool&);
    int _jet_cat_lookup(const std::string&);
    int _region_lookup(const std::string&);
    void _sample_lookup(const std::string&, int&, Spin&, float&, float&);
    int _sample2class_lookup(const int&);
    bool _accept_evt(const int&, const bool&, const int&, const int&, const int&);
    int _cut_lookup(const std::string&);

public:
    // Methods
	FileLooper(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true,
               bool apply_cut=true, bool inc_all_jets=true, bool inc_other_regions=false, bool inc_data=false, bool inc_unc=false);
	~FileLooper();
	bool loop_file(const std::string&, const std::string&, const std::string&, const std::string&, const long int&);
    std::map<unsigned long, std::string> build_id_map(TFile*);
};

#endif /* FILE_LOOPER_HH_ */
