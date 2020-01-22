#include "cms_runII_data_proc/processing/interface/file_looper.hh"

FileLooper::FileLooper(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true,
                       bool apply_cut=true, inc_all_jets=true, inc_other_regions=false, inc_data=false) {
    _evt_proc = new EvtProc(return_all, requested, use_deep_csv);
    _feat_names = EvtProc.get_feats();
    _n_feats = _feats.size();
    _apply_cut = apply_cut;
    _inc_all_jets = inc_all_jets;
    _inc_other_regions = inc_other_regions;
}

~FileLooper() {}

bool FileLooper::loop_file(const std::string& in_dir, const std::string& out_dir, const std::string& channel, const std::year& year, const long int& n_events) {
    /*
    Loop though file {in_dir}/{year}_{channel}.root processing {n_events} (all events if n_events < 0).
    Processed events will be saved to one of two trees inside {out_dir}/{year}_{channel}.root:
    Even event IDs will be saved to data_0 and odd to data_1.
    */

    TFile* in_file = TFile::Open((in_dir+"/"+year+"_"+channel+".root").c_str());
    TTreeReader reader(channel, in_file);
    TTreeReader aux_reader('"aux', in_file);

    // Enums
    Channel e_channel = FileLooper::_get_channel(channel);
    Year e_year = FileLooper::_get_year(year);

    // Meta info
    TTreeReaderValue<float> rv_weight(reader, "all_weights");  // TODO: Check this
    TTreeReaderValue<unsigned long int> rv_id(reader, "dataIds");
    TTreeReaderValue<unsigned long int> rv_aux_id(aux_reader, "dataIds");
    TTreeReaderValue<std::string> rv_aux_name(aux_reader, "dataId_names");
    float weight;
    std::string name;
    int sample, region, jet_cat;
    bool cut, scale, syst_unc;
    unsigned long int id;
    unsigned int class;

    // HL feats
    TTreeReaderValue<float> rv_kinfit_mass(reader, "m_ttbb_kinfit");
    TTreeReaderValue<float> rv_kinfit_chi2(reader, "chi2_kinFit");
    TTreeReaderValue<float> rv_mt2(reader, "MT2");
    TTreeReaderValue<float> rv_mt_tot(reader, "mt_tot");
    TTreeReaderValue<float> rv_top_1_mass(reader, "mass_top1");
    TTreeReaderValue<float> rv_top_2_mass(reader, "mass_top2");
    TTreeReaderValue<float> rv_p_zetavisible(reader, "p_zetavisible");
    TTreeReaderValue<float> rv_p_zeta(reader, "p_zeta");
    float kinfit_mass, kinfit_chi2, mt2, top_1_mass, top_2_mass, p_zetavisible, p_zeta;

    // Tagging
    TTreeReaderValue<float> rv_b_1_csv(reader, "csv_b1");
    TTreeReaderValue<float> rv_b_2_csv(reader, "csv_b2");
    TTreeReaderValue<float> rv_b_1_deepcsv(reader, "deepcsv_b1");
    TTreeReaderValue<float> rv_b_2_deepcsv(reader, "deepcsv_b2");
    //TODO: is boosted
    float b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv;

    // SVFit feats
    TTreeReaderValue<float> rv_svfit_pT(reader, "pt_sv");
    TTreeReaderValue<float> rv_svfit_eta(reader, "eta_sv");
    TTreeReaderValue<float> rv_svfit_phi(reader, "phi_sv");
    TTreeReaderValue<float> rv_svfit_mass(reader, "m_sv");
    LorentzVector svfit;

    // l1 feats
    TTreeReaderValue<float> rv_l_1_pT(reader, "pt_1");
    TTreeReaderValue<float> rv_l_1_eta(reader, "eta_1");
    TTreeReaderValue<float> rv_l_1_phi(reader, "phi_1");
    TTreeReaderValue<float> rv_l_1_mass(reader, "m_1");
    TTreeReaderValue<float> rv_l_1_mt(reader, "mt_1");
    float l_1_mt;
    LorentzVector l_1;

    // l2 feats
    TTreeReaderValue<float> rv_l_2_pT(reader, "pt_2");
    TTreeReaderValue<float> rv_l_2_eta(reader, "eta_2");
    TTreeReaderValue<float> rv_l_2_phi(reader, "phi_2");
    TTreeReaderValue<float> rv_l_2_mass(reader, "m_2");
    TTreeReaderValue<float> rv_l_2_mt(reader, "mt_2");
    float l_2_mass, l_2_mt;
    LorentzVector l_2;

    // MET feats
    TTreeReaderValue<float> rv_met_pT(reader, "pt_MET");
    TTreeReaderValue<float> rv_met_phi(reader, "phiMET");
    LorentzVector met;

    // b1 feats
    TTreeReaderValue<float> rv_b_1_pT(reader, "pt_b1");
    TTreeReaderValue<float> rv_b_1_eta(reader, "eta_b1");
    TTreeReaderValue<float> rv_b_1_phi(reader, "phi_b1");
    TTreeReaderValue<float> rv_b_1_mass(reader, "m_b1");
    LorentzVector b_1;

    // b2 feats
    TTreeReaderValue<float> rv_b_2_pT(reader, "pt_b2");
    TTreeReaderValue<float> rv_b_2_eta(reader, "eta_b2");
    TTreeReaderValue<float> rv_b_2_phi(reader, "phi_b2");
    TTreeReaderValue<float> rv_b_2_mass(reader, "m_b2");
    LorentzVector b_2;

    std::vector<float*> feat_vals(_n_feats);
    
    // Outfiles
    TFile* out_file  = new TFile((out_dir+"/"+year+"_"+channel+".root").c_str(), "recreate");
    TTree* data_even = new TTree("data_0", "Even id data");
    TTree* data_odd  = new TTree("data_1", "Odd id data");
    FileLooper::_prep_file(data_even, feat_vals);
    FileLooper::_prep_file(data_odd,  feat_vals);

    unsigned long c_event = 0;
    while (reader.Next()) {
        id = *rv_data_id;
        name = FileLooper::_get_evt_name(aux_reader, rv_aux_id, rv_aux_name, id);
        FileLooper::_extract_flags(std::string& name, sample, region, syst_unc, scale, jet_cat, cut, class);
        if (!FileLooper::_accept_evt(region, syst_unc, scale, jet_cat, cut, class)) continue;

        // Load meta
        weight = *rv_weight;

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
        svfit.SetCoordinates(*rv_svfit_pT, *rv_svfit_eta, *rv_svfit_phi, *rv_svfit_mass);
        l_1.SetCoordinates(*rv_l_1_pT, *rv_l_1_eta, *rv_l_1_phi, *rv_l_1_mass);
        if (channel == "muTau") {  // Fix mass for light leptons
            l_2_mass == MU_MASS;
        } else if (channel == "eTau") {
            l_2_mass == E_MASS;
        } else {
            l_2_mass = *rv_l_2_mass;
        }
        l_2.SetCoordinates(*rv_l_2_pT, *rv_l_2_eta, *rv_l_2_phi, l_2_mass);
        met.SetCoordinates(*rv_met_pT, 0,           *rv_met_phi, 0);
        b_1.SetCoordinates(*rv_b_1_pT, *rv_b_1_eta, *rv_b_1_phi, *rv_b_1_mass);
        b_2.SetCoordinates(*rv_b_2_pT, *rv_b_2_eta, *rv_b_2_phi, *rv_b_2_mass);

        _evt_proc.process_to_vec(feat_vals, b_1, b_2, l_1, l_2, met, svfit, kinfit_mass, kinfit_chi2, mt2, mt_tot, p_zetavisible, p_zeta, 
                                 top_1_mass, top_2_mass, l_1_mt, l_2_mt, is_boosted, csv_1, csv_2, deepcsv_1, deepcsv_2, e_channel, e_year, res_mass);

        c_event++;
        if (c_event >= n_events) break;
    }

    return true;
}

void FileLooper::_prep_file(TFile* tree, std::vector<float*>& feat_vals) {
    /* Add branches to tree and set addresses for values */

    for (unsigned int=0; i < _n_feats; i++) tree->Branch(_feat_names[i], feat_val[i]);
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
    return Channel;
}

Year FileLooper::_get_year(year) {
    /* Convert year to enum */

    if (channel == "y16") {
        return Year(y16);
    } else if (channel == "y17") {
        return Year(y17);
    } else if (channel == "y18") {
        return Year(y18);
    }
    throw std::invalid_argument("Invalid year: options are y16, y17, y18");
    return Year;
}

bool FileLooper::_get_evt_name(TTreeReader& aux_reader, TTreeReaderValue& rv_aux_id, TTreeReaderValue& rv_aux_name, const unsigned long int& id) {
    /* Match data ID to aux name */
    
    aux_reader.reset();
    std::string name;
    while (aux_reader.next()) {
        if (*rv_aux_id == id) {
            name = *rv_aux_name;
            break;
        }
    }
    return name;
}

void FileLooper::_extract_flags(const std::string& name, int& sample, int& region, bool& syst_unc, bool& scale, int& jet_cat, bool& cut, int& class) {
    /*
    Extract event flags from name 
    Example: "2j/NoCuts/SS_AntiIsolated/None/Central/DY_MC_M-10-50"
    */

    std::string val;
    std::istringstream iss(name);
    int i = 0;
    while (std::getline(iss, val, '/')) {   
        if (i == 0) {
            jet_cat = FileLooper::_jet_cat_lookup(val);
        } else if (i == 1) {
            cut = (val == "NoCuts");
        } else if (i == 2) {
            region = FileLooper::_region_lookup(val);
        } else if (i == 3) {
            syst_unc = (val == "None");
        } else if (i == 4) {
            scale = (val == "Central");
        } else if (i == 5) {
            sample = FileLooper::_sample_lookup(val);
        }
        i++;
    }
    class = FileLooper::_sample2class_lookup(sample);
}

int FileLooper::_jet_cat_lookup(const std::string& jet_cat) {
    if (val == "2j")           return 0;
    if (val == "2j0bR_noVBF")  return 1;
    if (val == "2j1bR_noVBF")  return 2;
    if (val == "2j2b+R_noVBF") return 3;
    if (val == "4j1b+_VBF")    return 4;
    // TODO: Boosted?
}

int FileLooper::_region_lookup(const std::string& region) {
    if (val == "OS_Isolated")      return 0;
    if (val == "OS_AntiIsolated")  return 1;
    if (val == "SS_Isolated")      return 2;
    if (val == "SS_AntiIsolated")  return 3;
}

int FileLooper::_sample_lookup(const std::string& sample) {
    if //HERE
    
    Signal_NonRes = -125
    Signal_Radion = [-260, -270, -280, -300, -320, -340, -350, -400, -450, -500, -550, -600, -650, -750, -800, -900]
    Data  = 0
    TT    = 1
    DY    = 2
    Wjets = 3
    ZH    = 4
    WW    = 5
    WZ    = 6
    ZZ    = 7
    EWK   = 8
    tW    = 9
}

int FileLooper::_sample2class_lookup(const int& sample) {
    if (sample < 0)  return 1;   // Signal
    if (sample == 0) return -1;  // Collider data
    return 0;                    // Background
}

bool FileLooper::_accept_evt(const int& region, const bool& syst_unc, const bool& scale, const int& jet_cat, const bool& cut, const int& class) {
    if (_apply_cut && cut) return false;  // Require cut and cut failed
    if (!_inc_data && class == -1) return false; // Don't inlclude data and 
    if (!_inc_other_regions && region != 0)
    _inc_other_regions = inc_other_regions
}