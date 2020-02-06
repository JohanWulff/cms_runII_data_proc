#include "cms_runII_data_proc/processing/interface/file_looper.hh"

FileLooper::FileLooper(bool return_all, std::set<std::string> requested, bool use_deep_csv,
                       bool apply_cut, bool inc_all_jets, bool inc_other_regions, bool inc_data, bool inc_unc) {
    _evt_proc = new EvtProc(return_all, requested, use_deep_csv);
    _feat_names = _evt_proc->get_feats();
    _n_feats = _feat_names.size();
    _apply_cut = apply_cut;
    _inc_all_jets = inc_all_jets;
    _inc_other_regions = inc_other_regions;
    _inc_unc = inc_unc;
    _inc_data = inc_data;
}

FileLooper::~FileLooper() {
    delete _evt_proc;
}

bool FileLooper::loop_file(const std::string& in_dir, const std::string& out_dir, const std::string& channel, const std::string& year,
                           const long int& n_events) {
    /*
    Loop though file {in_dir}/{year}_{channel}.root processing {n_events} (all events if n_events < 0).
    Processed events will be saved to one of two trees inside {out_dir}/{year}_{channel}.root:
    Even event IDs will be saved to data_0 and odd to data_1.
    */

    std::string fname = in_dir+"/"+year+"_"+channel+".root";
    std::cout << "Reading from file: " << fname << "\n";
    TFile* in_file = TFile::Open(fname.c_str());
    TTreeReader reader(channel.c_str(), in_file);

    // Enums
    Channel e_channel = FileLooper::_get_channel(channel);
    Year e_year = FileLooper::_get_year(year);
    Spin spin(nonres);
    float klambda;

    // Meta info
    std::cout << "Extracting auxiliary data...";
    TTreeReaderValue<std::vector<double>> rv_weight(reader, "all_weights");
    TTreeReaderValue<std::vector<float>> rv_weight(reader, "evt");
    TTreeReaderValue<std::vector<unsigned long>> rv_id(reader, "dataIds");
    std::map<unsigned long, std::string> id2name = FileLooper::build_id_map(in_file);
    std::cout << " Extracted\n";
    double weight;
    float res_mass, evt;
    std::vector<std::string> names;
    int sample, region, jet_cat, cut;
    unsigned long long int strat_key;
    bool scale, syst_unc;
    std::vector<unsigned long> ids;
    int class_id;

    // HL feats
    TTreeReaderValue<float> rv_kinfit_mass(reader, "m_ttbb_kinfit");
    TTreeReaderValue<float> rv_kinfit_chi2(reader, "chi2_kinFit");
    TTreeReaderValue<float> rv_mt2(reader, "MT2");
    TTreeReaderValue<float> rv_mt_tot(reader, "mt_tot");
    TTreeReaderValue<float> rv_top_1_mass(reader, "mass_top1");
    TTreeReaderValue<float> rv_top_2_mass(reader, "mass_top2");
    TTreeReaderValue<float> rv_p_zetavisible(reader, "p_zetavisible");
    TTreeReaderValue<float> rv_p_zeta(reader, "p_zeta");
    float kinfit_mass, kinfit_chi2, mt2, mt_tot, top_1_mass, top_2_mass, p_zetavisible, p_zeta;

    // Tagging
    TTreeReaderValue<float> rv_b_1_csv(reader, "csv_b1");
    TTreeReaderValue<float> rv_b_2_csv(reader, "csv_b2");
    TTreeReaderValue<float> rv_b_1_deepcsv(reader, "deepcsv_b1");
    TTreeReaderValue<float> rv_b_2_deepcsv(reader, "deepcsv_b2");
    float b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv;
    bool is_boosted;

    // SVFit feats
    //TTreeReaderValue<float> rv_svfit_pT(reader, "pt_sv");
    //TTreeReaderValue<float> rv_svfit_eta(reader, "eta_sv");
    //TTreeReaderValue<float> rv_svfit_phi(reader, "phi_sv");
    TTreeReaderValue<float> rv_svfit_mass(reader, "m_sv");
    LorentzVectorPEP pep_svfit;
    LorentzVector svfit;

    // l1 feats
    TTreeReaderValue<float> rv_l_1_pT(reader, "pt_1");
    TTreeReaderValue<float> rv_l_1_eta(reader, "eta_1");
    TTreeReaderValue<float> rv_l_1_phi(reader, "phi_1");
    TTreeReaderValue<float> rv_l_1_mass(reader, "m_1");
    TTreeReaderValue<float> rv_l_1_mt(reader, "mt_1");
    float l_1_mt;
    LorentzVectorPEP pep_l_1;
    LorentzVector l_1;


    // l2 feats
    TTreeReaderValue<float> rv_l_2_pT(reader, "pt_2");
    TTreeReaderValue<float> rv_l_2_eta(reader, "eta_2");
    TTreeReaderValue<float> rv_l_2_phi(reader, "phi_2");
    TTreeReaderValue<float> rv_l_2_mass(reader, "m_2");
    TTreeReaderValue<float> rv_l_2_mt(reader, "mt_2");
    float l_2_mass, l_2_mt;
    LorentzVectorPEP pep_l_2;\
    LorentzVector l_2;

    // MET feats
    TTreeReaderValue<float> rv_met_pT(reader, "pt_MET");
    TTreeReaderValue<float> rv_met_phi(reader, "phiMET");
    LorentzVectorPEP pep_met;
    LorentzVector met;

    // b1 feats
    TTreeReaderValue<float> rv_b_1_pT(reader, "pt_b1");
    TTreeReaderValue<float> rv_b_1_eta(reader, "eta_b1");
    TTreeReaderValue<float> rv_b_1_phi(reader, "phi_b1");
    TTreeReaderValue<float> rv_b_1_mass(reader, "m_b1");
    LorentzVectorPEP pep_b_1;
    LorentzVector b_1;

    // b2 feats
    TTreeReaderValue<float> rv_b_2_pT(reader, "pt_b2");
    TTreeReaderValue<float> rv_b_2_eta(reader, "eta_b2");
    TTreeReaderValue<float> rv_b_2_phi(reader, "phi_b2");
    TTreeReaderValue<float> rv_b_2_mass(reader, "m_b2");
    LorentzVectorPEP pep_b_2;
    LorentzVector b_2;

    // vbf1 feats
    TTreeReaderValue<float> rv_vbf_1_pT(reader, "pt_VBF_1");
    TTreeReaderValue<float> rv_vbf_1_eta(reader, "eta_VBF_1");
    TTreeReaderValue<float> rv_vbf_1_phi(reader, "phi_VBF_1");
    TTreeReaderValue<float> rv_vbf_1_mass(reader, "m_VBF_1");
    LorentzVectorPEP pep_vbf_1;
    LorentzVector vbf_1;

    // vbf2 feats
    TTreeReaderValue<float> rv_vbf_2_pT(reader, "pt_VBF_2");
    TTreeReaderValue<float> rv_vbf_2_eta(reader, "eta_VBF_2");
    TTreeReaderValue<float> rv_vbf_2_phi(reader, "phi_VBF_2");
    TTreeReaderValue<float> rv_vbf_2_mass(reader, "m_VBF_2");
    LorentzVectorPEP pep_vbf_2;
    LorentzVector vbf_2;

    std::vector<std::unique_ptr<float>> feat_vals;
    feat_vals.reserve(_n_feats);
    for (unsigned int i = 0; i < _n_feats; i++) feat_vals.emplace_back(new float(0));
    
    // Outfiles
    std::string oname = out_dir+"/"+year+"_"+channel+".root";
    std::cout << "Preparing output file: " << oname << " ...";
    TFile* out_file  = new TFile(oname.c_str(), "recreate");
    TTree* data_even = new TTree("data_0", "Even id data");
    TTree* data_odd  = new TTree("data_1", "Odd id data");
    FileLooper::_prep_file(data_even, feat_vals, &weight, &sample, &region, &jet_cat, &cut, &scale, &syst_unc, &class_id, &strat_key);
    FileLooper::_prep_file(data_odd,  feat_vals, &weight, &sample, &region, &jet_cat, &cut, &scale, &syst_unc, &class_id, &strat_key);
    std::cout << "\tprepared.\nBeginning loop.\n";

    long int c_event(0), n_tot_events(reader.GetEntries(true));
    while (reader.Next()) {
        c_event++;
        if (c_event%1000 == 0) std::cout << c_event << " / " << n_tot_events << "\n";
        ids = *rv_id;

        names = FileLooper::_get_evt_names(id2name, ids);
        FileLooper::_extract_flags(names, sample, region, syst_unc, scale, jet_cat, cut, class_id, spin, klambda, res_mass, is_boosted);
        if (!FileLooper::_accept_evt(region, syst_unc, jet_cat, cut, class_id)) continue;
        strat_key = FileLooper::_get_strat_key(sample, static_cast<int>(jet_cat), region, static_cast<int>(spin), static_cast<int>(syst_unc), cut);

        // Load meta
        weight = (*rv_weight)[0];
        evt    =  *evt;

        // Load HL feats
        kinfit_mass   = *rv_kinfit_mass;
        kinfit_chi2   = *rv_kinfit_chi2;
        mt2           = *rv_mt2;
        mt_tot        = *rv_mt_tot;
        top_1_mass    = *rv_top_1_mass;
        top_2_mass    = *rv_top_2_mass;
        p_zetavisible = *rv_p_zetavisible;
        p_zeta        = *rv_p_zeta;
        l_1_mt        = *rv_l_1_mt;
        l_2_mt        = *rv_l_2_mt;

        // Load tagging
        b_1_csv     = *rv_b_1_csv;
        b_2_csv     = *rv_b_2_csv;
        b_1_deepcsv = *rv_b_1_deepcsv;
        b_2_deepcsv = *rv_b_2_deepcsv;

        // Load vectors
        // pep_svfit.SetCoordinates(*rv_svfit_pT, *rv_svfit_eta, *rv_svfit_phi, *rv_svfit_mass);
        pep_svfit.SetCoordinates(0, 0, 0, *rv_svfit_mass); // TODO remove once svfit present
        pep_l_1.SetCoordinates(*rv_l_1_pT, *rv_l_1_eta, *rv_l_1_phi, *rv_l_1_mass);
        if (channel == "muTau") {  // Fix mass for light leptons
            l_2_mass = MU_MASS;
        } else if (channel == "eTau") {
            l_2_mass = E_MASS;
        } else {
            l_2_mass = *rv_l_2_mass;
        }
        pep_l_2.SetCoordinates(*rv_l_2_pT, *rv_l_2_eta, *rv_l_2_phi, l_2_mass);
        pep_met.SetCoordinates(*rv_met_pT, 0,           *rv_met_phi, 0);
        pep_b_1.SetCoordinates(*rv_b_1_pT, *rv_b_1_eta, *rv_b_1_phi, *rv_b_1_mass);
        pep_b_2.SetCoordinates(*rv_b_2_pT, *rv_b_2_eta, *rv_b_2_phi, *rv_b_2_mass);
        pep_vbf_1.SetCoordinates(*rv_vbf_1_pT, *rv_vbf_1_eta, *rv_vbf_1_phi, *rv_vbf_1_mass);
        pep_vbf_2.SetCoordinates(*rv_vbf_2_pT, *rv_vbf_2_eta, *rv_vbf_2_phi, *rv_vbf_2_mass);

        svfit.SetCoordinates(pep_svfit.Px(), pep_svfit.Py(), pep_svfit.Pz(), pep_svfit.M());
        l_1.SetCoordinates(pep_l_1.Px(),     pep_l_1.Py(),   pep_l_1.Pz(),   pep_l_1.M());
        l_2.SetCoordinates(pep_l_2.Px(),     pep_l_2.Py(),   pep_l_2.Pz(),   pep_l_2.M());
        met.SetCoordinates(pep_met.Px(),     pep_met.Py(),   0,              0);
        b_1.SetCoordinates(pep_b_1.Px(),     pep_b_1.Py(),   pep_b_1.Pz(),   pep_b_1.M());
        b_2.SetCoordinates(pep_b_2.Px(),     pep_b_2.Py(),   pep_b_2.Pz(),   pep_b_2.M());
        vbf_1.SetCoordinates(pep_vbf_1.Px(), pep_vbf_1.Py(), pep_vbf_1.Pz(), pep_vbf_1.M());
        vbf_2.SetCoordinates(pep_vbf_2.Px(), pep_vbf_2.Py(), pep_vbf_2.Pz(), pep_vbf_2.M());

        _evt_proc->process_to_vec(feat_vals, b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, kinfit_mass, kinfit_chi2, mt2, mt_tot, p_zetavisible, p_zeta, top_1_mass,
                                  top_2_mass, l_1_mt, l_2_mt, is_boosted, b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv, e_channel, e_year, res_mass, spin,
                                  klambda);

        std::cout << evt << "\n";
        if (c_event%2 == 0) {  // TODO: Replace with evt once implemented
            data_even->Fill();
        } else {
            data_odd->Fill();
        }        
        if (n_events > 0 && c_event >= n_events) {
            std::cout << "Exiting after " << c_event << " events.\n";
            break;
        }
    }

    std::cout << "Loop complete, saving results.\n";
    data_even->Write();
    data_odd->Write();
    delete data_even;
    delete data_odd;
    in_file->Close();
    out_file->Close();
    return true;
}

void FileLooper::_prep_file(TTree* tree, const std::vector<std::unique_ptr<float>>& feat_vals, double* weight, int* sample, int* region, int* jet_cat,
                            int* cut, bool* scale, bool* syst_unc, int* class_id, unsigned long long int* strat_key) {
    /* Add branches to tree and set addresses for values */

    for (unsigned int i = 0; i < _n_feats; i++) tree->Branch(_feat_names[i].c_str(), feat_vals[i].get());
    tree->Branch("weight",    weight);
    tree->Branch("sample",    sample);
    tree->Branch("region",    region);
    tree->Branch("jet_cat",   jet_cat);
    tree->Branch("cut",       cut);
    tree->Branch("scale",     scale);
    tree->Branch("syst_unc",  syst_unc);
    tree->Branch("class_id",  class_id);
    tree->Branch("strat_key", strat_key);
}

Channel FileLooper::_get_channel(std::string channel) {
    /* COnvert channel to enum */

    if (channel == "tauTau") {
        return Channel(tauTau);
    } else if (channel == "muTau") {
        return Channel(muTau);
    } else if (channel == "eTau") {
        return Channel(eTau);
    }
    throw std::invalid_argument("Invalid channel: options are tauTau, muTau, eTau");
    return Channel(tauTau);
}

Year FileLooper::_get_year(std::string year) {
    /* Convert year to enum */

    if (year == "2016") {
        return Year(y16);
    } else if (year == "2017") {
        return Year(y17);
    } else if (year == "2018") {
        return Year(y18);
    }
    throw std::invalid_argument("Invalid year: options are 2016, 2017, 2018");
    return Year(y16);
}

std::vector<std::string> FileLooper::_get_evt_names(const std::map<unsigned long, std::string>& id2name, const std::vector<unsigned long>& ids) {
    /* Match data IDs to aux names */
    
    std::vector<std::string> names(ids.size());
    for (unsigned int i = 0; i < ids.size(); i++) names[i] = id2name.at(ids[i]);
    return names;
}

void FileLooper::_extract_flags(const std::vector<std::string>& name, int& sample, int& region, bool& syst_unc, bool& scale, int& jet_cat, int& cut, int& class_id,
                                Spin& spin, float& klambda, float& res_mass, bool& is_boosted) {
    /*
    Extract event flags from name 
    Example: "2j/NoCuts/SS_AntiIsolated/None/Central/DY_MC_M-10-50"
    */

    std::string val;
    int tmp;
    jet_cat = -1;
    cut = -1;
    for (unsigned int n = 0; n < name.size(); n++) {
        std::istringstream iss(name[n]);
        int i = 0;
        while (std::getline(iss, val, '/')) {   
            if (i == 0) {
                tmp = FileLooper::_jet_cat_lookup(val);
                if (tmp > jet_cat) jet_cat = tmp;
                is_boosted = (jet_cat == 5);  // TODO: update this
            } else if (i == 1) {
                tmp = FileLooper::_cut_lookup(val);
                if (tmp > cut) cut = tmp;
            } else if (i == 2 && n == 0) {
                region = FileLooper::_region_lookup(val);
            } else if (i == 3 && n == 0) {
                syst_unc = (val == "None");
            } else if (i == 4 && n == 0) {
                scale = (val == "Central");
            } else if (i == 5 && n == 0) {
                FileLooper::_sample_lookup(val, sample, spin, klambda, res_mass);
            }
            i++;
        }
        class_id = FileLooper::_sample2class_lookup(sample);
    }
}

int FileLooper::_cut_lookup(const std::string& cut) {
    if (cut == "NoCuts") return 0;
    if (cut == "mhVis")  return 1;
    throw std::invalid_argument("Unrecognised cut category: " + cut);
    return -1;
}

int FileLooper::_jet_cat_lookup(const std::string& jet_cat) {
    if (jet_cat == "2j")            return 0;
    if (jet_cat == "2j0bR_noVBF")   return 1;
    if (jet_cat == "2j1bR_noVBF")   return 2;
    if (jet_cat == "2j2b+R_noVBF")  return 3;
    if (jet_cat == "4j1b+_VBF")     return 4;
    if (jet_cat == "2j2Lb+B_noVBF") return 5;
    throw std::invalid_argument("Unrecognised jet category: " + jet_cat);
    return -1;
}

int FileLooper::_region_lookup(const std::string& region) {
    if (region == "OS_Isolated")      return 0;
    if (region == "OS_AntiIsolated")  return 1;
    if (region == "SS_Isolated")      return 2;
    if (region == "SS_AntiIsolated")  return 3;
    if (region == "SS_LooseIsolated") return 4;
    throw std::invalid_argument("Unrecognised region: " + region);
    return -1;
}

void FileLooper::_sample_lookup(const std::string& sample, int& sample_id, Spin& spin, float& klambda, float& res_mass) {
    spin = nonres;
    res_mass = 125;
    klambda = 1;
    
    if (sample.find("Signal_NonRes") != std::string::npos) {
        if (sample.find("_kl") != std::string::npos) {
            sample_id = -12;
            try {
                klambda = std::stof(sample.substr(sample.find("_kl")+3));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_kl")+3) << "\n";
                assert(false);
            }
        } else {  // VBF
            sample_id = -13;
        }
    } else if (sample.find("Signal_Radion") != std::string::npos) {
        spin = radion;
        try {
            res_mass = std::stof(sample.substr(sample.find("_M")+2));
        } catch (...) {
            std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
            assert(false);
        }
        if (res_mass <= 400) {
            sample_id = -14;
        } else if (res_mass <= 600) {
            sample_id = -15;
        } else {
            sample_id = -16;
        }
    } else if (sample.find("Signal_Graviton") != std::string::npos) {
        spin = graviton;
        try {
            res_mass = std::stof(sample.substr(sample.find("_M")+2));
        } catch (...) {
            std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
            assert(false);
        }
        if (res_mass <= 400) {
            sample_id = -17;
        } else if (res_mass <= 600) {
            sample_id = -18;
        } else {
            sample_id = -19;
        }
    } else if (sample.find("Data") != std::string::npos) {
        sample_id = 0;
    } else if (sample.find("TT") != std::string::npos) {
        sample_id = 1;
    } else if (sample.find("ttH") != std::string::npos) {
        sample_id = 2;
    } else if (sample.find("DY") != std::string::npos) {
        sample_id = 3;
    } else if (sample.find("Wjets") != std::string::npos) {
        sample_id = 4;
    } else if (sample.find("SM_Higgs") != std::string::npos) {
        sample_id = 5;
    } else if (sample.find("VH") != std::string::npos) {
        sample_id = 6;
    } else if (sample.find("VVV") != std::string::npos) {
        sample_id = 7;
    } else if (sample.find("EWK") != std::string::npos) {
        sample_id = 8;
    } else if (sample.find("VV") != std::string::npos) {
        sample_id = 9;
    } else if (sample.find("ST") != std::string::npos) {
        sample_id = 10;
    } else if (sample.find("ttV") != std::string::npos) {
        sample_id = 11;
    } else{
        throw std::invalid_argument("Unrecognised sample: " + sample);
    }
}

int FileLooper::_sample2class_lookup(const int& sample) {
    if (sample < 0)  return 1;   // Signal
    if (sample == 0) return -1;  // Collider data
    return 0;                    // Background
}

bool FileLooper::_accept_evt(const int& region, const bool& syst_unc, const int& jet_cat, const int& cut, const int& class_id) {
    if (_apply_cut && cut == 0) return false;  // Require cut and cut failed
    if (!_inc_data && class_id == -1) return false; // Don't inlclude data and event is data
    if (!_inc_other_regions && region != 0) return false;  //Don't include other regions and event is not SS Iso
    if (!_inc_unc && !syst_unc) return false;  //Don't systematicss and event is a systematic
    if (!_inc_all_jets && jet_cat == 0) return false;  // Only use inference category jets and event is non-inference category
    return true;
}

unsigned long long int FileLooper::_get_strat_key(const int& sample, const int& jet_cat, const int& region, const int& spin, const int& syst_unc,
                                                  const int& cut) {
    unsigned long long int strat_key = std::pow(2,  std::abs(sample))*
                                       std::pow(3,  jet_cat)*
                                       std::pow(5, region)*
                                       std::pow(7, spin)*
                                       std::pow(11, cut)*
                                       std::pow(13, syst_unc);
    if (strat_key == 0) throw std::overflow_error("Strat key overflow\n");    
    return strat_key;
}

std::map<unsigned long, std::string> FileLooper::build_id_map(TFile* in_file) {
    TTreeReader aux_reader("aux", in_file);
    TTreeReaderValue<std::vector<unsigned long>> rv_aux_id(aux_reader, "dataIds");
    TTreeReaderValue<std::vector<std::string>> rv_aux_name(aux_reader, "dataId_names");
    std::vector<unsigned long> ids;
    std::vector<std::string> names;
    while (aux_reader.Next()) {
        ids   = *rv_aux_id;
        names = *rv_aux_name;
    }
    std::map<unsigned long, std::string> id2name;
    for (unsigned int i = 0; i < ids.size(); i++) id2name[ids[i]] = names[i];
    return id2name;
}
