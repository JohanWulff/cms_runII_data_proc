#include "cms_runII_data_proc/processing/interface/file_looper.hh"
#include "cms_runII_data_proc/processing/interface/kinfitter.hh"

int use_kl = 1;


FileLooper::FileLooper(bool return_all, std::vector<std::string> requested, bool use_deep_bjet_wps,
                       bool inc_all_jets, bool inc_other_regions, bool inc_data, bool only_kl1, bool only_sm_vbf) {
    _evt_proc = new EvtProc(return_all, requested, use_deep_bjet_wps);
    _feat_names = _evt_proc->get_feats();
    _n_feats = _feat_names.size();
    _inc_all_jets = inc_all_jets;
    _inc_other_regions = inc_other_regions;
    _inc_data = inc_data;
    _only_kl1 = only_kl1;
    _only_sm_vbf = only_sm_vbf;
}

FileLooper::~FileLooper() {
    delete _evt_proc;
}

bool FileLooper::loop_file(const std::string& in_dir, const std::string& out_dir, const std::string& channel, const bool add_zz_zh_feats,
                           const std::string& year, const long int& n_events, const long int& start_evt, const long int& end_evt) {
    /*
    Loop though file {in_dir}/{year}_{channel}.root processing {n_events} (all events if n_events < 0).
    Processed events will be saved to one of two trees inside {out_dir}/{year}_{channel}_{evt_range}.root:
    Even event IDs will be saved to data_0 and odd to data_1.
    */

    std::string fname = in_dir+"/"+year+"_"+channel+"_Central.root";
    std::cout << "Reading from file: " << fname << "\n";
    TFile* in_file = TFile::Open(fname.c_str());
    TTreeReader reader(channel.c_str(), in_file);
    TTreeReader aux_reader("aux", in_file);

    // Enums
    Channel e_channel = FileLooper::_get_channel(channel);
    Year e_year = FileLooper::_get_year(year);
    Spin spin(nonres);
    float klambda, cv, c2v, c3;

    // Meta info
    std::cout << "Extracting auxiliary data...";

    TTreeReaderValue<unsigned long long> rv_evt(reader, "evt");
    TTreeReaderValue<float> rv_weight(reader, "weight");
    TTreeReaderValue<UInt_t> rv_dataset_id(reader, "dataset");
    TTreeReaderValue<UInt_t> rv_region_id(reader, "event_region");
    std::map<unsigned, std::string> id2dataset = FileLooper::build_dataset_id_map(in_file);
    std::map<unsigned, std::string> id2region = FileLooper::build_region_id_map(in_file);
    std::cout << " Extracted\n";
    
    float weight, res_mass;
    int sample, region, jet_cat, n_vbf, class_id;
    unsigned long long int strat_key, evt;
    bool svfit_conv, hh_kinfit_conv;

    // Gen Info
    TTreeReaderValue<int> rv_tau1_gen_match(reader, "tau1_gen_match");
    TTreeReaderValue<int> rv_tau2_gen_match(reader, "tau2_gen_match");
    TTreeReaderValue<int> rv_b1_hadronFlavour(reader, "b1_hadronFlavour");
    TTreeReaderValue<int> rv_b2_hadronFlavour(reader, "b2_hadronFlavour");
    int tau1_gen_match, tau2_gen_match, b1_hadronFlavour, b2_hadronFlavour;

    // HL feats
    TTreeReaderValue<float> rv_kinfit_mass(reader, "kinFit_m");
    TTreeReaderValue<float> rv_kinfit_chi2(reader, "kinFit_chi2");
    TTreeReaderValue<float> rv_mt2(reader, "MT2");
    float kinfit_mass, kinfit_chi2, mt2;

    // ZZ & ZH KinFit;
    std::pair<float,float> kinfit_ZZ, kinfit_ZH;
    std::vector<float> kinINinfo;

    // Tagging
    TTreeReaderValue<float> rv_b_1_csv(reader, "b1_DeepFlavour");
    TTreeReaderValue<float> rv_b_2_csv(reader, "b2_DeepFlavour");
    TTreeReaderValue<bool> rv_is_boosted(reader, "is_boosted");
    TTreeReaderValue<bool> rv_has_b_pair(reader, "has_b_pair");
    TTreeReaderValue<bool> rv_has_vbf_pair(reader, "has_VBF_pair");
    TTreeReaderValue<int> rv_num_btag_loose(reader, "num_btag_Loose");
    TTreeReaderValue<int> rv_num_btag_medium(reader, "num_btag_Medium");
    float b_1_csv, b_2_csv;
    bool is_boosted, has_vbf_pair, has_b_pair;
    int num_btag_loose, num_btag_medium;

    // SVFit feats
    TTreeReaderValue<float> rv_svfit_pT(reader, "SVfit_pt");
    TTreeReaderValue<float> rv_svfit_eta(reader, "SVfit_eta");
    TTreeReaderValue<float> rv_svfit_phi(reader, "SVfit_phi");
    TTreeReaderValue<float> rv_svfit_mass(reader, "SVfit_m");
    LorentzVectorPEP pep_svfit;
    LorentzVector svfit;

    // l1 feats
    TTreeReaderValue<float> rv_l_1_pT(reader, "tau1_pt");
    TTreeReaderValue<float> rv_l_1_eta(reader, "tau1_eta");
    TTreeReaderValue<float> rv_l_1_phi(reader, "tau1_phi");
    TTreeReaderValue<float> rv_l_1_mass(reader, "tau1_m");
    float l_1_mass;
    LorentzVectorPEP pep_l_1;
    LorentzVector l_1;

    // l2 feats
    TTreeReaderValue<float> rv_l_2_pT(reader, "tau2_pt");
    TTreeReaderValue<float> rv_l_2_eta(reader, "tau2_eta");
    TTreeReaderValue<float> rv_l_2_phi(reader, "tau2_phi");
    TTreeReaderValue<float> rv_l_2_mass(reader, "tau2_m");
    LorentzVectorPEP pep_l_2;
    LorentzVector l_2;

    // MET feats
    TTreeReaderValue<float> rv_met_pT(reader, "MET_pt");
    TTreeReaderValue<float> rv_met_phi(reader, "MET_phi");
    TTreeReaderValue<float> rv_met_cov_00(reader, "MET_cov_00");
    TTreeReaderValue<float> rv_met_cov_01(reader, "MET_cov_01");
    TTreeReaderValue<float> rv_met_cov_11(reader, "MET_cov_11");
    LorentzVectorPEP pep_met;
    LorentzVector met;

    // b1 feats
    TTreeReaderValue<float> rv_b_1_pT(reader, "b1_pt");
    TTreeReaderValue<float> rv_b_1_eta(reader, "b1_eta");
    TTreeReaderValue<float> rv_b_1_phi(reader, "b1_phi");
    TTreeReaderValue<float> rv_b_1_mass(reader, "b1_m");
    TTreeReaderValue<float> rv_b_1_hhbtag(reader, "b1_HHbtag");
    TTreeReaderValue<float> rv_b_1_cvsl(reader, "b1_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_b_1_cvsb(reader, "b1_DeepFlavour_CvsB");
    float b_1_hhbtag, b_1_cvsl, b_1_cvsb;
    LorentzVectorPEP pep_b_1;
    LorentzVector b_1;

    // b2 feats
    TTreeReaderValue<float> rv_b_2_pT(reader, "b2_pt");
    TTreeReaderValue<float> rv_b_2_eta(reader, "b2_eta");
    TTreeReaderValue<float> rv_b_2_phi(reader, "b2_phi");
    TTreeReaderValue<float> rv_b_2_mass(reader, "b2_m");
    TTreeReaderValue<float> rv_b_2_hhbtag(reader, "b2_HHbtag");
    TTreeReaderValue<float> rv_b_2_cvsl(reader, "b2_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_b_2_cvsb(reader, "b2_DeepFlavour_CvsB");
    float b_2_hhbtag, b_2_cvsl, b_2_cvsb;
    LorentzVectorPEP pep_b_2;
    LorentzVector b_2;

    // vbf1 feats
    TTreeReaderValue<float> rv_vbf_1_pT(reader, "VBF1_pt");
    TTreeReaderValue<float> rv_vbf_1_eta(reader, "VBF1_eta");
    TTreeReaderValue<float> rv_vbf_1_phi(reader, "VBF1_phi");
    TTreeReaderValue<float> rv_vbf_1_mass(reader, "VBF1_m");
    TTreeReaderValue<float> rv_vbf_1_hhbtag(reader, "VBF1_HHbtag");
    TTreeReaderValue<float> rv_vbf_1_cvsl(reader, "VBF1_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_vbf_1_cvsb(reader, "VBF1_DeepFlavour_CvsB");
    float vbf_1_hhbtag, vbf_1_cvsl, vbf_1_cvsb;
    LorentzVectorPEP pep_vbf_1;
    LorentzVector vbf_1;

    // vbf2 feats
    TTreeReaderValue<float> rv_vbf_2_pT(reader, "VBF2_pt");
    TTreeReaderValue<float> rv_vbf_2_eta(reader, "VBF2_eta");
    TTreeReaderValue<float> rv_vbf_2_phi(reader, "VBF2_phi");
    TTreeReaderValue<float> rv_vbf_2_mass(reader, "VBF2_m");
    TTreeReaderValue<float> rv_vbf_2_hhbtag(reader, "VBF2_HHbtag");
    TTreeReaderValue<float> rv_vbf_2_cvsl(reader, "VBF2_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_vbf_2_cvsb(reader, "VBF2_DeepFlavour_CvsB");
    float vbf_2_hhbtag, vbf_2_cvsl, vbf_2_cvsb;
    LorentzVectorPEP pep_vbf_2;
    LorentzVector vbf_2;

    std::vector<std::unique_ptr<float>> feat_vals;
    feat_vals.reserve(_n_feats);
    for (unsigned int i = 0; i < _n_feats; i++) feat_vals.emplace_back(new float(0));
    
    // Outfiles
    std::string oname = out_dir+"/"+year+"_"+channel+"_"+std::to_string(start_evt)+"-"+std::to_string(end_evt)+".root";
    std::cout << "Preparing output file: " << oname << " ...";
    TFile* out_file  = new TFile(oname.c_str(), "recreate");
    TTree* data_even = new TTree("data_0", "Even id data");
    TTree* data_odd  = new TTree("data_1", "Odd id data");
    FileLooper::_prep_file(data_even, feat_vals, &weight, &sample, &region, &jet_cat, &class_id, &strat_key,
                           &tau1_gen_match, &tau2_gen_match, &b1_hadronFlavour, &b2_hadronFlavour);
    FileLooper::_prep_file(data_odd,  feat_vals, &weight, &sample, &region, &jet_cat, &class_id, &strat_key,
                           &tau1_gen_match, &tau2_gen_match, &b1_hadronFlavour, &b2_hadronFlavour);

    
    std::vector<std::unique_ptr<float>> zz_feat_vals;
    TTree* data_zz_even;
    TTree* data_zz_odd;
    std::vector<std::unique_ptr<float>> zh_feat_vals;
    TTree* data_zh_even;
    TTree* data_zh_odd;
    if (add_zz_zh_feats) {
        std::cout << "\tIncluding extra features for ZZ and ZH searches...";
        zz_feat_vals.reserve(_n_feats);
        for (unsigned int i = 0; i < _n_feats; i++) zz_feat_vals.emplace_back(new float(0));
        data_zz_even = new TTree("data_ZZ_0", "Even id ZZ data");
        data_zz_odd  = new TTree("data_ZZ_1", "Odd id ZZ data");
        FileLooper::_prep_file(data_zz_even, zz_feat_vals, &weight, &sample, &region, &jet_cat, &class_id, &strat_key,
                               &tau1_gen_match, &tau2_gen_match, &b1_hadronFlavour, &b2_hadronFlavour);
        FileLooper::_prep_file(data_zz_odd,  zz_feat_vals, &weight, &sample, &region, &jet_cat, &class_id, &strat_key,
                               &tau1_gen_match, &tau2_gen_match, &b1_hadronFlavour, &b2_hadronFlavour);

        zh_feat_vals.reserve(_n_feats);
        for (unsigned int i = 0; i < _n_feats; i++) zh_feat_vals.emplace_back(new float(0));
        data_zh_even = new TTree("data_ZH_0", "Even id ZH data");
        data_zh_odd  = new TTree("data_ZH_1", "Odd id ZH data");
        FileLooper::_prep_file(data_zh_even, zh_feat_vals, &weight, &sample, &region, &jet_cat, &class_id, &strat_key,
                               &tau1_gen_match, &tau2_gen_match, &b1_hadronFlavour, &b2_hadronFlavour);
        FileLooper::_prep_file(data_zh_odd,  zh_feat_vals, &weight, &sample, &region, &jet_cat, &class_id, &strat_key,
                               &tau1_gen_match, &tau2_gen_match, &b1_hadronFlavour, &b2_hadronFlavour);
    }
    std::cout << "\tprepared.\nBeginning loop.\n";

    long int c_event(0), n_saved_events(0), n_tot_events(reader.GetEntries(true));
    while (reader.Next()) {
        c_event++;
        if (c_event < start_evt) continue;
        if (end_evt > 0 && c_event >= end_evt) break;
        if (c_event%1000 == 0) std::cout << c_event << " / " << n_tot_events << "\n";

        // Load meta
        weight = *rv_weight;
        evt    = *rv_evt;
        
        FileLooper::_sample_lookup(id2dataset[*rv_dataset_id], sample, spin, klambda, res_mass, cv, c2v, c3);
        class_id = FileLooper::_sample2class_lookup(sample);
        
        region = FileLooper::_region_lookup(id2region[*rv_region_id]);
        
        is_boosted = *rv_is_boosted;
        has_vbf_pair = *rv_has_vbf_pair;
        has_b_pair = *rv_has_b_pair;
        num_btag_loose = *rv_num_btag_loose;
        num_btag_medium = *rv_num_btag_medium;

        jet_cat = FileLooper::_jet_cat_lookup(has_b_pair, has_vbf_pair, is_boosted, num_btag_loose, num_btag_medium);
        
        if (!FileLooper::_accept_evt(region, jet_cat, class_id, klambda, cv, c2v, c3)) continue;
        n_saved_events++;

        strat_key = FileLooper::_get_strat_key(sample, jet_cat, e_channel, e_year, region);

        // Gen info
        tau1_gen_match = *rv_tau1_gen_match;
        tau2_gen_match = *rv_tau2_gen_match;
        b1_hadronFlavour = *rv_b1_hadronFlavour;
        b2_hadronFlavour = *rv_b2_hadronFlavour;

        // Load HL feats
        kinfit_mass   = *rv_kinfit_mass;
        kinfit_chi2   = *rv_kinfit_chi2;
        mt2           = *rv_mt2;
        b_1_hhbtag    = *rv_b_1_hhbtag;
        b_2_hhbtag    = *rv_b_2_hhbtag;
        vbf_1_hhbtag  = *rv_vbf_1_hhbtag;
        vbf_2_hhbtag  = *rv_vbf_2_hhbtag;
        b_1_cvsl      = *rv_b_1_cvsl;
        b_2_cvsl      = *rv_b_2_cvsl;
        vbf_1_cvsl    = *rv_vbf_1_cvsl;
        vbf_2_cvsl    = *rv_vbf_2_cvsl;
        b_1_cvsb      = *rv_b_1_cvsb;
        b_2_cvsb      = *rv_b_2_cvsb;
        vbf_1_cvsb    = *rv_vbf_1_cvsb;
        vbf_2_cvsb    = *rv_vbf_2_cvsb;

        // Load tagging
        b_1_csv     = *rv_b_1_csv;
        b_2_csv     = *rv_b_2_csv;

        // Load vectors
        pep_svfit.SetCoordinates(*rv_svfit_pT, *rv_svfit_eta, *rv_svfit_phi, *rv_svfit_mass);
        if (channel == "muTau") {  // Fix mass for light leptons
            l_1_mass = MU_MASS;
        } else if (channel == "eTau") {
            l_1_mass = E_MASS;
        } else {
            l_1_mass = *rv_l_1_mass;
        }
        pep_l_1.SetCoordinates(*rv_l_1_pT, *rv_l_1_eta, *rv_l_1_phi, l_1_mass);
        pep_l_2.SetCoordinates(*rv_l_2_pT, *rv_l_2_eta, *rv_l_2_phi, *rv_l_2_mass);
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

        // VBF
        n_vbf = has_vbf_pair ? 2 : 0;

        // Convergence
        svfit_conv     = *rv_svfit_mass > 0;
        hh_kinfit_conv = kinfit_chi2    > 0;

        _evt_proc->process_to_vec(feat_vals, b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, kinfit_mass, kinfit_chi2, mt2, is_boosted, b_1_csv, b_2_csv,
                                  e_channel, e_year, res_mass, spin, klambda, n_vbf, svfit_conv, hh_kinfit_conv, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag,
                                  vbf_2_hhbtag, b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb, cv, c2v, c3, true);
        
        if (evt%2 == 0) {
            data_even->Fill();
        } else {
            data_odd->Fill();
        }

        if (add_zz_zh_feats) {
            // KinFit for ZZ/ZH
            // create a single object with all the needed info to give kinfit { 4 lep1 coords, 4 lep2 coords, 4 bjet1 coords, 4 bjet2 coords, 2 MET coors, 3 MET cov entries }
            kinINinfo = { *rv_l_1_pT, *rv_l_1_eta, *rv_l_1_phi, l_1_mass, *rv_l_2_pT, *rv_l_2_eta, *rv_l_2_phi, *rv_l_2_mass ,*rv_b_1_pT, *rv_b_1_eta, *rv_b_1_phi,
                         *rv_b_1_mass ,*rv_b_2_pT, *rv_b_2_eta, *rv_b_2_phi, *rv_b_2_mass ,*rv_met_pT, *rv_met_phi, *rv_met_cov_00, *rv_met_cov_01, *rv_met_cov_11 };
            // compute KinFit 
            KinFitter fitter(kinINinfo);
            kinfit_ZZ = fitter.fit("ZZ");
            kinfit_ZH = fitter.fit("ZH");
            kinINinfo.clear();

            _evt_proc->process_to_vec(zz_feat_vals, b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, kinfit_ZZ.first, kinfit_ZZ.second, mt2, is_boosted, b_1_csv, b_2_csv,
                                      e_channel, e_year, res_mass, spin, klambda, n_vbf, svfit_conv, kinfit_ZZ.second > 0, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag,
                                      vbf_2_hhbtag, b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb, cv, c2v, c3, true);
            
            _evt_proc->process_to_vec(zh_feat_vals, b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, kinfit_ZH.first, kinfit_ZH.second, mt2, is_boosted, b_1_csv, b_2_csv,
                                      e_channel, e_year, res_mass, spin, klambda, n_vbf, svfit_conv, kinfit_ZH.second > 0, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag,
                                      vbf_2_hhbtag, b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb, cv, c2v, c3, true);

            if (evt%2 == 0) {
                data_zz_even->Fill();
                data_zh_even->Fill();
            } else {
                data_zz_odd->Fill();
                data_zh_odd->Fill();
            }
        }

        if (n_events > 0 && n_saved_events >= n_events) {
            std::cout << "Exiting after " << n_saved_events << " events.\n";
            break;
        }
    }

    std::cout << "Loop complete, saving results.\n";
    data_even->Write();
    data_odd->Write();
    if (add_zz_zh_feats) {
        data_zz_even->Write();
        data_zz_odd->Write();
        data_zh_even->Write();
        data_zh_odd->Write();
        delete data_zz_even;
        delete data_zz_odd;
        delete data_zh_even;
        delete data_zh_odd;
    }
    delete data_even;
    delete data_odd;
    in_file->Close();
    out_file->Close();
    return true;
}

std::map<unsigned, std::string> FileLooper::build_dataset_id_map(TFile* in_file) {
    TTreeReader aux_reader("aux", in_file);
    TTreeReaderValue<std::vector<std::string>> rv_dataset_names(aux_reader, "dataset_names");
    TTreeReaderValue<std::vector<unsigned>> rv_dataset_hashes(aux_reader, "dataset_hashes");
    
    std::vector<unsigned> ids;
    std::vector<std::string> names;
    while (aux_reader.Next()) {
        ids   = *rv_dataset_hashes;
        names = *rv_dataset_names;
    }
    std::map<unsigned, std::string> id2name;
    for (unsigned int i = 0; i < ids.size(); i++) id2name[ids[i]] = names[i];
    return id2name;
}

std::map<unsigned, std::string> FileLooper::build_region_id_map(TFile* in_file) {
    TTreeReader aux_reader("aux", in_file);
    TTreeReaderValue<std::vector<std::string>> rv_region_names(aux_reader, "region_names");
    TTreeReaderValue<std::vector<unsigned>> rv_region_hashes(aux_reader, "region_hashes");
    
    std::vector<unsigned> ids;
    std::vector<std::string> names;
    while (aux_reader.Next()) {
        ids   = *rv_region_hashes;
        names = *rv_region_names;
    }
    std::map<unsigned, std::string> id2name;
    for (unsigned int i = 0; i < ids.size(); i++) id2name[ids[i]] = names[i];
    return id2name;
}

void FileLooper::_prep_file(TTree* tree, const std::vector<std::unique_ptr<float>>& feat_vals, float* weight, int* sample, int* region, int* jet_cat,
                            int* class_id, unsigned long long int* strat_key,
                            int* tau1_gen_match, int* tau2_gen_match, int* b1_hadronFlavour, int* b2_hadronFlavour) {
    /* Add branches to tree and set addresses for values */

    for (unsigned int i = 0; i < _n_feats; i++) tree->Branch(_feat_names[i].c_str(), feat_vals[i].get());
    tree->Branch("weight",      weight);
    tree->Branch("sample",      sample);
    tree->Branch("region",      region);
    tree->Branch("jet_cat",     jet_cat);
    tree->Branch("tau1_gen_match", tau1_gen_match);
    tree->Branch("tau2_gen_match", tau2_gen_match);
    tree->Branch("b1_hadronFlavour", b1_hadronFlavour);
    tree->Branch("b2_hadronFlavour", b2_hadronFlavour);
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

int FileLooper::_jet_cat_lookup(const bool has_b_pair, const bool has_vbf_pair, const bool is_boosted, const int num_btag_Loose, const int num_btag_Medium) {
    if (!has_b_pair)  return -1;
    if (has_vbf_pair && num_btag_Loose >= 1) return 5;  // 2j1b+_VBFL, 2j1b+_VBF, 2j1b+_VBFT
    if (!has_vbf_pair && is_boosted && num_btag_Loose >= 2)   return 4;  // 2j2Lb+B_noVBF
    if (!has_vbf_pair && num_btag_Medium >= 2) return 3;  // 2j2b+R_noVBF
    if (!has_vbf_pair && num_btag_Loose >= 1)  return 2;  // 2j1bR_noVBF
    if (!has_vbf_pair && num_btag_Loose == 0)  return 1;  // 2j0bR_noVBF
    return 0;  // 2j 
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

void FileLooper::_sample_lookup(std::string& sample, int& sample_id, Spin& spin, float& klambda, float& res_mass, float& cv, float& c2v, float& c3) {
    // TODO Update this
    
    spin = nonres;
    res_mass = 125;
    klambda = use_kl;
    cv = 1;
    c2v = 1;
    c3 = 1;
    
    if (sample.find("GluGluSignal") != std::string::npos) { 
        if (sample.find("NonRes_klScan") != std::string::npos) {
            sample_id = (sample.find("nlo") != std::string::npos) ? -26 : -12;
            try {
                klambda = std::stof(sample.substr(sample.find("_klScan_kl")+10));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_klScan_kl")+10) << "\n";
                assert(false);
            }
        } else if (sample.find("Radion") != std::string::npos) {
            spin = radion;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -13;
            } else if (res_mass <= 600) {
                sample_id = -14;
            } else {
                sample_id = -15;
            }
        } else if (sample.find("Graviton") != std::string::npos) {
            spin = graviton;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -16;
            } else if (res_mass <= 600) {
                sample_id = -17;
            } else {
                sample_id = -18;
            }
        } else {
            sample_id = -999; 
        }
    } else if (sample.find("VBFSignal") != std::string::npos) {
        if (sample.find("NonRes") != std::string::npos) {
            sample_id = (sample.find("nlo") != std::string::npos) ? -27 : -19;
            try {
                cv = std::stof(sample.substr(sample.find("CV_")+3, sample.find("_C2V")));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("CV")+2, sample.find("_C2V")) << "\n";
                assert(false);
            }
            try {
                c2v = std::stof(sample.substr(sample.find("C2V_")+4, sample.find("_C3")));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("C2V_")+4, sample.find("_C3")) << "\n";
                assert(false);
            }
            try {
                c3 = std::stof(sample.substr(sample.find("C3_")+3));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("C3_")+3) << "\n";
                assert(false);
            }
        } else if (sample.find("Radion") != std::string::npos) {
            spin = radion;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -20;
            } else if (res_mass <= 600) {
                sample_id = -21;
            } else {
                sample_id = -22;
            }
        } else if (sample.find("Graviton") != std::string::npos) {
            spin = graviton;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -23;
            } else if (res_mass <= 600) {
                sample_id = -24;
            } else {
                sample_id = -25;
            }
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
    } else if (sample.find("GluGluH") != std::string::npos || sample.find("VBFH") != std::string::npos) {
        sample_id = 5;
    } else if (sample.find("ZHToTauTau_M125") != std::string::npos) {
        sample_id = 6;
    } else if (sample.find("ZH_HToBB_ZToLL_M125") != std::string::npos) {
        sample_id = 7;
    } else if (sample.find("ZH_HToBB_ZToQQ_M125") != std::string::npos) {
        sample_id = 8;
    } else if (sample.find("WminusH") != std::string::npos || sample.find("WplusH") != std::string::npos) {
        sample_id = 9;
    } else if (sample.find("WWW") != std::string::npos) {
        sample_id = 10;
    } else if (sample.find("WWZ") != std::string::npos) {
        sample_id = 11;
    } else if (sample.find("WZZ") != std::string::npos) {
        sample_id = 12;
    } else if (sample.find("ZZZ") != std::string::npos) {
        sample_id = 13;
    } else if (sample.find("EWK") != std::string::npos) {
        sample_id = 14;
    } else if (sample.find("WW") != std::string::npos) {
        sample_id = 15;
    } else if (sample.find("WZ") != std::string::npos) {
        sample_id = 16;
    } else if (sample.find("ZZTo2L2Nu") != std::string::npos) {
        sample_id = 17;
    } else if (sample.find("ZZTo2L2Q") != std::string::npos) {
        sample_id = 18;
    } else if (sample.find("ZZTo2Q2Nu") != std::string::npos) {
        sample_id = 19;
    } else if (sample.find("ZZTo4L") != std::string::npos) {
        sample_id = 20;
    } else if (sample.find("ZZTo4Q") != std::string::npos) {
        sample_id = 21;
    } else if (sample.find("ST") != std::string::npos) {
        sample_id = 22;
    } else{
        throw std::invalid_argument("Unrecognised sample: " + sample);
    }
}

int FileLooper::_sample2class_lookup(const int& sample) {
    if (sample == -999)  return -999;  // Sample to reject
    if (sample < 0)      return 1;     // Signal
    if (sample == 0)     return -1;    // Collider data
    return 0;                          // Background
}

bool FileLooper::_accept_evt(const int& region, const int& jet_cat, const int& class_id, const float& klambda,
                             const float& cv, const float& c2v, const float& c3) {
    if (_only_kl1 && klambda != use_kl) {
        // std::cout << "Rejecting due to klambda = " << klambda << "\n";
        return false; // Only consider klambda at SM point
    }
    if (!_inc_data && class_id == -1) {
        // std::cout << "Rejecting due to class ID = " << class_id << "\n";
        return false; // Don't include data and event is data
    }
    if (!_inc_other_regions && region != 0) {
        // std::cout << "Rejecting due to region = " << region << "\n";
        return false;  //Don't include other regions and event is not SS Iso
    }
    if (!_inc_all_jets && jet_cat == 0) {
        // std::cout << "Rejecting due to jet_cat = " << jet_cat << "\n";
        return false;  // Only use inference category jets and event is non-inference category
    }
    if (jet_cat < 0) {
        // std::cout << "Rejecting due to jet_cat = " << jet_cat << "\n";
        return false;  // Don't include experimental jet categories
    }
    if (_only_sm_vbf && (cv != 1 || c2v != 1 || c3 != 1)) {
        // std::cout << "Rejecting due to cv = " << cv << " c2v = " << c2v << " c3 = " << c3 << "\n";
        return false; // Only consider SM VBF
    } if (class_id == -999) {
        // std::cout << "Rejecting due to class ID = -999\n";
        return false; // Reject sample
    }
    // std::cout << "Accepting\n";
    return true;
}

unsigned long long int FileLooper::_get_strat_key(const int& sample, const int& jet_cat, const Channel& channel, const Year& year, const int& region) {
    unsigned long long int strat_key = std::pow(2,  std::abs(sample))*
                                       std::pow(3,  jet_cat)*
                                       std::pow(5,  (float)channel)*
                                       std::pow(7,  (float)year)*
                                       std::pow(11, region);
    if (strat_key == 0) {
        std::cout << "sample " << sample << " jet_cat " << jet_cat << " channel " << channel << " year " << year << " region " << region << "\n";
        throw std::overflow_error("Strat key overflow\n");  
    }  
    return strat_key;
}
