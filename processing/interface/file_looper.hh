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
	using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

	// Variables
    bool _all, _use_deep_csv, _apply_cut, _inc_other_regions, _inc_all_jets;
    std::set<std::string> _requested;
    unsigned int _n_feats;
    std::vector<std::string> _feat_names;
    EvtProc* _evt_proc;

	// Methods
    inline int _get_split(const unsigned long int&);
    void FileLooper::_prep_file(TTree*, std::vector<float*>& feats);
    Channel _get_channel(std::string);
    Year _get_year(std::string);
    // int _get_strat_key()
    std::string _get_evt_name_get_evt_name(TTreeReader&, TTreeReaderValue&, TTreeReaderValue&, const unsigned long int&);
    void _extract_flags(const std::string&, int&, int&, bool&, bool&, int&, bool&);
    int _jet_cat_lookup(const std::string&);
    int _region_lookup(const std::string&);
    int _sample_lookup(const std::string&);
    int _sample2class_lookup(const std::string&);
    bool _accept_evt(const int&, const bool&, const bool&, const int&, const bool&, const int&);

public:
    // Methods
	FileLooper(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true,
               bool apply_cut=true, inc_all_jets=true, inc_other_regions=false, bool inc_data=false);
	~FileLooper();
	bool loop_file(const& std::string&, const std::string&, const std::string&, const std::string&, const long int&, const bool&);
};

#endif /* FILE_LOOPER_HH_ */
