#ifndef FILE_LOOPER_HH_
#define FILE_LOOPER_HH_

// C++
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>
#include <TFile.h>

// Delphes
#include <ExRootAnalysis/ExRootTreeReader.h>

// Plugins
#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"
#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"

const double eMass = 0.0005109989; //GeV
const double muMass = 0.1056583715; //GeV

class FileLooper {
	/* Class for processing data in a ROOT file event by event */

private:
	// Names
	using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all, _use_deep_csv;
    std::set<std::string> _requested;
    TFile* _output_file;


	// Methods
	void _load_file(const& std::string);
    inline int _get_split(unsigned long int);
    // int _get_strat_key()

public:
    // Methods
	FileLooper(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true);
	~FileLooper();
	bool loop_file(const& std::string, const& std::string, unsigned long int, bool);
};

#endif /* FILE_LOOPER_HH_ */
