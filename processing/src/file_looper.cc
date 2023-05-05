#include "cms_runII_data_proc/processing/interface/file_looper.hh"
#include "cms_runII_data_proc/processing/interface/kinfitter.hh"
#include "TH1D.h"

int use_kl = 1;

FileLooper::FileLooper(bool return_all, std::vector<std::string> requested, bool use_deep_bjet_wps,
                       bool inc_all_jets, bool inc_other_regions, bool inc_data, bool only_kl1, bool only_sm_vbf)
{
    _evt_proc = new EvtProc(return_all, requested, use_deep_bjet_wps);
    _feat_names = _evt_proc->get_feats();
    _n_feats = _feat_names.size();
    _inc_all_jets = inc_all_jets;
    _inc_other_regions = inc_other_regions;
    _inc_data = inc_data;
    _only_kl1 = only_kl1;
    _only_sm_vbf = only_sm_vbf;
}

FileLooper::~FileLooper()
{
    delete _evt_proc;
}

bool FileLooper::loop_file(const std::string &fname, const std::string &oname, const std::string &channel,
                           const std::string &year, std::string &sample_str, double sum_w,
                           const long int &n_events, const long int &start_evt, const long int &end_evt)
{
    /*
    Loop though file {in_dir}/{year}_{channel}.root processing {n_events} (all events if n_events < 0).
    Processed events will be saved to one of two trees inside {out_dir}/{year}_{channel}_{evt_range}.root:
    Even event IDs will be saved to data_0 and odd to data_1.
    */

    std::cout << "Reading from file: " << fname << "\n";
    TFile *in_file = TFile::Open(fname.c_str());
    TTreeReader reader("HTauTauTree", in_file);

    // Enums
    Channel e_channel = FileLooper::_get_channel(channel);
    Year e_year = FileLooper::_get_year(year);
    Spin spin(nonres);
    float klambda, cv, c2v, c3;

    // Meta info
    std::cout << "Extracting auxiliary data...";

    TTreeReaderValue<unsigned long long> rv_evt(reader, "EventNumber");

    TTreeReaderValue<float> rv_MC_weight(reader, "MC_weight");
    TTreeReaderValue<float> rv_prescaleWeight(reader, "prescaleWeight");
    TTreeReaderValue<float> rv_L1pref_weight(reader, "L1pref_weight");
    TTreeReaderValue<float> rv_PUjetID_SF(reader, "PUjetID_SF");
    TTreeReaderValue<float> rv_PUReweight(reader, "PUReweight");
    TTreeReaderValue<float> rv_bTagweightReshape(reader, "bTagweightReshape");
    TTreeReaderValue<float> rv_trigSF(reader, "trigSF");
    TTreeReaderValue<float> rv_IdAndIsoAndFakeSF_deep_pt(reader, "IdAndIsoAndFakeSF_deep_pt");
    TTreeReaderValue<float> rv_DYscale_MTT(reader, "DYscale_MTT");

    double weight;
    float bTagweightReshape, PUReweight, PUjetID_SF, L1pref_weight, prescaleWeight, MC_weight;
    float trigSF, DYscale_MTT, IdAndIsoAndFakeSF_deep_pt;

    std::cout << " Extracted\n";

    int sample_id, jet_cat, n_vbf, class_id;
    unsigned long long int strat_key, evt;
    bool svfit_conv, hh_kinfit_conv;

    // HL feats
    TTreeReaderValue<float> rv_tauH_mass(reader, "tauH_mass");
    TTreeReaderValue<float> rv_tauH_pt(reader, "tauH_pt");
    TTreeReaderValue<float> rv_tauH_eta(reader, "tauH_eta");
    TTreeReaderValue<float> rv_tauH_phi(reader, "tauH_phi");
    TTreeReaderValue<float> rv_tauH_e(reader, "tauH_e");
    float tauH_mass, tauH_pt, tauH_eta, tauH_phi, tauH_e;

    TTreeReaderValue<float> rv_bH_mass(reader, "bH_mass");
    TTreeReaderValue<float> rv_bH_pt(reader, "bH_pt");
    TTreeReaderValue<float> rv_bH_eta(reader, "bH_eta");
    TTreeReaderValue<float> rv_bH_phi(reader, "bH_phi");
    TTreeReaderValue<float> rv_bH_e(reader, "bH_e");
    float bH_mass, bH_pt, bH_eta, bH_phi, bH_e;

    TTreeReaderValue<float> rv_HH_mass(reader, "HH_mass");
    TTreeReaderValue<float> rv_HH_mass_raw(reader, "HH_mass_raw");
    TTreeReaderValue<float> rv_HH_pt(reader, "HH_pt");
    TTreeReaderValue<float> rv_HH_eta(reader, "HH_eta");
    TTreeReaderValue<float> rv_HH_phi(reader, "HH_phi");
    TTreeReaderValue<float> rv_HH_e(reader, "HH_e");
    float HH_mass, H_mass_raw, HH_pt, HH_eta, HH_phi, HH_e;
    
    TTreeReaderValue<float> rv_kinfit_mass(reader, "HHKin_mass_raw");
    TTreeReaderValue<float> rv_kinfit_chi2(reader, "HHKin_mass_raw_chi2");
    TTreeReaderValue<float> rv_mt2(reader, "MT2");
    float kinfit_mass, kinfit_chi2, mt2;

    // Selection Stuff
    TTreeReaderValue<int> rv_pairType(reader, "pairType");
    TTreeReaderValue<int> rv_nleps(reader, "nleps");
    TTreeReaderValue<int> rv_nbjetscand(reader, "nbjetscand");
    TTreeReaderValue<int> rv_isLeptrigger(reader, "isLeptrigger");
    int pairType, nleps, nbjetscand, isLeptrigger;

    // Tagging
    TTreeReaderValue<float> rv_b_1_csv(reader, "bjet1_bID_deepFlavor");
    TTreeReaderValue<float> rv_b_2_csv(reader, "bjet2_bID_deepFlavor");
    TTreeReaderValue<int> rv_is_boosted(reader, "isBoosted");
    float b_1_csv, b_2_csv;
    int is_boosted;

    // SVFit feats
    TTreeReaderValue<float> rv_svfit_pT(reader, "tauH_SVFIT_pt");
    TTreeReaderValue<float> rv_svfit_eta(reader, "tauH_SVFIT_eta");
    TTreeReaderValue<float> rv_svfit_phi(reader, "tauH_SVFIT_phi");
    TTreeReaderValue<float> rv_svfit_mass(reader, "tauH_SVFIT_mass");
    LorentzVectorPEP pep_svfit;
    LorentzVector svfit;

    // l1 feats
    TTreeReaderValue<float> rv_l_1_pT(reader, "dau1_pt");
    TTreeReaderValue<float> rv_l_1_eta(reader, "dau1_eta");
    TTreeReaderValue<float> rv_l_1_phi(reader, "dau1_phi");
    TTreeReaderValue<float> rv_l_1_e(reader, "dau1_e");
    TTreeReaderValue<float> rv_l_1_iso(reader, "dau1_iso");
    TTreeReaderValue<int> rv_l_1_eleMVAiso(reader, "dau1_eleMVAiso");
    LorentzVectorPEP pep_l_1;
    LorentzVector l_1;
    float dau1_iso;
    int dau1_eleMVAiso;

    // l2 feats
    TTreeReaderValue<float> rv_l_2_pT(reader, "dau2_pt");
    TTreeReaderValue<float> rv_l_2_eta(reader, "dau2_eta");
    TTreeReaderValue<float> rv_l_2_phi(reader, "dau2_phi");
    TTreeReaderValue<float> rv_l_2_e(reader, "dau2_e");
    LorentzVectorPEP pep_l_2;
    LorentzVector l_2;

    // MET feats
    TTreeReaderValue<float> rv_met_pT(reader, "met_et"); // seems strange.. I know
    TTreeReaderValue<float> rv_met_phi(reader, "met_phi");
    TTreeReaderValue<float> rv_met_cov_00(reader, "met_cov00");
    TTreeReaderValue<float> rv_met_cov_01(reader, "met_cov01");
    TTreeReaderValue<float> rv_met_cov_11(reader, "met_cov11");
    LorentzVectorPEP pep_met;
    LorentzVector met;

    // b1 feats
    TTreeReaderValue<float> rv_b_1_pT(reader, "bjet1_pt");
    TTreeReaderValue<float> rv_b_1_eta(reader, "bjet1_eta");
    TTreeReaderValue<float> rv_b_1_phi(reader, "bjet1_phi");
    TTreeReaderValue<float> rv_b_1_e(reader, "bjet1_e");
    TTreeReaderValue<float> rv_b_1_hhbtag(reader, "bjet1_HHbtag");
    TTreeReaderValue<float> rv_b_1_cvsl(reader, "bjet1_CvsL");
    TTreeReaderValue<float> rv_b_1_cvsb(reader, "bjet1_CvsB");
    float b_1_hhbtag, b_1_cvsl, b_1_cvsb;
    LorentzVectorPEP pep_b_1;
    LorentzVector b_1;

    // b2 feats
    TTreeReaderValue<float> rv_b_2_pT(reader, "bjet2_pt");
    TTreeReaderValue<float> rv_b_2_eta(reader, "bjet2_eta");
    TTreeReaderValue<float> rv_b_2_phi(reader, "bjet2_phi");
    TTreeReaderValue<float> rv_b_2_e(reader, "bjet2_e");
    TTreeReaderValue<float> rv_b_2_hhbtag(reader, "bjet2_HHbtag");
    TTreeReaderValue<float> rv_b_2_cvsl(reader, "bjet2_CvsL");
    TTreeReaderValue<float> rv_b_2_cvsb(reader, "bjet2_CvsB");
    float b_2_hhbtag, b_2_cvsl, b_2_cvsb;
    LorentzVectorPEP pep_b_2;
    LorentzVector b_2;

    // VBF trigger
    TTreeReaderValue<int> rv_isVBFtrigger(reader, "isVBFtrigger");
    TTreeReaderValue<int> rv_isVBF(reader, "isVBF");

    // VBF jj
    TTreeReaderValue<float> rv_VBFjj_mass(reader, "VBFjj_mass");
    TTreeReaderValue<float> rv_VBFjj_deltaEta(reader, "VBFjj_deltaEta");

    // vbf1 feats
    TTreeReaderValue<float> rv_vbf_1_pT(reader, "VBFjet1_pt");
    TTreeReaderValue<float> rv_vbf_1_eta(reader, "VBFjet1_eta");
    TTreeReaderValue<float> rv_vbf_1_phi(reader, "VBFjet1_phi");
    TTreeReaderValue<float> rv_vbf_1_e(reader, "VBFjet1_e");
    TTreeReaderValue<float> rv_vbf_1_hhbtag(reader, "VBFjet1_HHbtag");
    TTreeReaderValue<float> rv_vbf_1_cvsl(reader, "VBFjet1_CvsL");
    TTreeReaderValue<float> rv_vbf_1_cvsb(reader, "VBFjet1_CvsB");
    float vbf_1_hhbtag, vbf_1_cvsl, vbf_1_cvsb;
    LorentzVectorPEP pep_vbf_1;
    LorentzVector vbf_1;

    // vbf2 feats
    TTreeReaderValue<float> rv_vbf_2_pT(reader, "VBFjet2_pt");
    TTreeReaderValue<float> rv_vbf_2_eta(reader, "VBFjet2_eta");
    TTreeReaderValue<float> rv_vbf_2_phi(reader, "VBFjet2_phi");
    TTreeReaderValue<float> rv_vbf_2_e(reader, "VBFjet2_e");
    TTreeReaderValue<float> rv_vbf_2_hhbtag(reader, "VBFjet2_HHbtag");
    TTreeReaderValue<float> rv_vbf_2_cvsl(reader, "VBFjet2_CvsL");
    TTreeReaderValue<float> rv_vbf_2_cvsb(reader, "VBFjet2_CvsB");
    float vbf_2_hhbtag, vbf_2_cvsl, vbf_2_cvsb;
    LorentzVectorPEP pep_vbf_2;
    LorentzVector vbf_2;

    // region vars
    TTreeReaderValue<int> rv_isOS(reader, "isOS");
    TTreeReaderValue<int> rv_dau1_deepTauVsJet(reader, "dau1_deepTauVsJet");
    TTreeReaderValue<int> rv_dau2_deepTauVsJet(reader, "dau2_deepTauVsJet");
    int dau1_deepTauVsJet, dau2_deepTauVsJet;
    int region_id, sample, isOS;

    std::vector<std::unique_ptr<float>> feat_vals;
    feat_vals.reserve(_n_feats);
    for (unsigned int i = 0; i < _n_feats; i++)
        feat_vals.emplace_back(new float(0));

    // Outfiles
    std::cout << "Preparing output file: " << oname << " ...";
    TFile *out_file = new TFile(oname.c_str(), "recreate");
    TTree *data_even = new TTree("data_0", "Even id data");
    TTree *data_odd = new TTree("data_1", "Odd id data");
    FileLooper::_prep_file(data_even, feat_vals, &weight, &sample_id, &region_id, &jet_cat, &class_id, &strat_key);

    data_even->Branch("tauH_mass", &tauH_mass);
    data_even->Branch("tauH_e", &tauH_e);
    data_even->Branch("tauH_phi", &tauH_phi);
    data_even->Branch("tauH_pt", &tauH_pt);
    data_even->Branch("tauH_eta", &tauH_eta);

    data_even->Branch("bH_mass", &bH_mass);
    data_even->Branch("bH_e", &bH_e);
    data_even->Branch("bH_phi", &bH_phi);
    data_even->Branch("bH_pt", &bH_pt);
    data_even->Branch("bH_eta", &bH_eta);

    data_even->Branch("HH_mass", &HH_mass);
    data_even->Branch("HH_e", &HH_e);
    data_even->Branch("HH_phi", &HH_phi);
    data_even->Branch("HH_pt", &HH_pt);
    data_even->Branch("HH_eta", &HH_eta);

    data_even->Branch("kinfit_mass", &kinfit_mass);
    FileLooper::_prep_file(data_odd, feat_vals, &weight, &sample_id, &region_id, &jet_cat, &class_id, &strat_key);

    data_odd->Branch("tauH_mass", &tauH_mass);
    data_odd->Branch("tauH_e", &tauH_e);
    data_odd->Branch("tauH_phi", &tauH_phi);
    data_odd->Branch("tauH_pt", &tauH_pt);
    data_odd->Branch("tauH_eta", &tauH_eta);

    data_odd->Branch("bH_mass", &bH_mass);
    data_odd->Branch("bH_e", &bH_e);
    data_odd->Branch("bH_phi", &bH_phi);
    data_odd->Branch("bH_pt", &bH_pt);
    data_odd->Branch("bH_eta", &bH_eta);

    data_odd->Branch("HH_mass", &HH_mass);
    data_odd->Branch("HH_e", &HH_e);
    data_odd->Branch("HH_phi", &HH_phi);
    data_odd->Branch("HH_pt", &HH_pt);
    data_odd->Branch("HH_eta", &HH_eta);

    data_odd->Branch("kinfit_mass", &kinfit_mass);

    std::cout << "\tprepared.\nBeginning loop.\n";

    long int c_event(0), n_saved_events(0), n_tot_events(reader.GetEntries(true));
    while (reader.Next())
    {
        float res_mass;
        FileLooper::_sample_lookup(sample_str, sample_id, spin, res_mass);

        bTagweightReshape = *rv_bTagweightReshape;
        PUReweight = *rv_PUReweight;
        PUjetID_SF = *rv_PUjetID_SF;
        L1pref_weight = *rv_L1pref_weight;
        prescaleWeight = *rv_prescaleWeight;
        MC_weight = *rv_MC_weight;
        trigSF = *rv_trigSF;
        DYscale_MTT = *rv_DYscale_MTT;
        IdAndIsoAndFakeSF_deep_pt = *rv_IdAndIsoAndFakeSF_deep_pt;
        // calc weight
        weight = 1.0;
        if (sample_id != 0)
        {
            weight *= bTagweightReshape * PUReweight * PUjetID_SF * L1pref_weight * prescaleWeight;
            weight *= MC_weight * trigSF * IdAndIsoAndFakeSF_deep_pt * DYscale_MTT;
            weight /= sum_w;
        }

        // baseline selection
        pairType = *rv_pairType;
        nleps = *rv_nleps;
        nbjetscand = *rv_nbjetscand;
        isLeptrigger = *rv_isLeptrigger;
        bool pass_baseline;
        bool has_b_pair = true;

        pass_baseline = FileLooper::_apply_baseline(channel, c_event, pairType, nleps, nbjetscand, isLeptrigger);
        if (pass_baseline == 0)
            continue;
        c_event++;
        if (c_event < start_evt)
            continue;
        if (end_evt > 0 && c_event >= end_evt)
            break;
        if (c_event % 1000 == 0)
            std::cout << c_event << " / " << n_tot_events << "\n";
        // Load meta
        evt = *rv_evt;

        class_id = FileLooper::_sample2class_lookup(sample_id);

        // determine region depending on channel
        dau1_deepTauVsJet = *rv_dau1_deepTauVsJet;
        dau2_deepTauVsJet = *rv_dau2_deepTauVsJet;
        isOS = *rv_isOS;
        dau1_iso = *rv_l_1_iso;
        dau1_eleMVAiso = *rv_l_1_eleMVAiso;
        region_id = FileLooper::_get_region(channel, isOS, dau1_deepTauVsJet, dau2_deepTauVsJet, dau1_iso, dau1_eleMVAiso);

        if (region_id == -1){
            continue;
        }

        is_boosted = *rv_is_boosted;
        bool boosted = is_boosted != 0;

        // Load tagging
        b_1_csv = *rv_b_1_csv;
        b_2_csv = *rv_b_2_csv;
        int num_btag_loose = 0;
        int num_btag_medium = 0;
        if ((b_1_csv > 0.049) ^ (b_2_csv > 0.049))
            num_btag_loose = 1;
        if ((b_1_csv > 0.049) && (b_2_csv > 0.049))
            num_btag_loose = 2;
        if ((b_1_csv > 0.2783) ^ (b_2_csv > 0.2783))
            num_btag_medium = 1;
        if ((b_1_csv > 0.2783) && (b_2_csv > 0.2783))
            num_btag_medium = 2;

        // VBF

        bool has_vbf_pair = false;
        if (*rv_isVBF == 1 && *rv_VBFjj_mass > 500 && *rv_VBFjj_deltaEta > 3 &&
            (((*rv_l_1_pT > 25 && *rv_l_2_pT > 25 && (*rv_l_1_pT <= 40 || *rv_l_2_pT <= 40)) &&
              *rv_VBFjj_mass > 800 && *rv_vbf_1_pT > 140 && *rv_vbf_2_pT > 60) ||
             *rv_isVBFtrigger == 0))
        {
            has_vbf_pair = true;
        }
        n_vbf = has_vbf_pair ? 2 : 0;
        jet_cat = FileLooper::_jet_cat_lookup(has_b_pair, has_vbf_pair, boosted, num_btag_loose, num_btag_medium);

        // if (!FileLooper::_accept_evt(region, jet_cat, class_id, klambda, cv, c2v, c3)) continue;
        n_saved_events++;

        strat_key = FileLooper::_get_strat_key(sample_id, jet_cat, e_channel, e_year, region_id);

        // Load HL feats
        kinfit_mass = *rv_kinfit_mass;
        kinfit_chi2 = *rv_kinfit_chi2;
        mt2 = *rv_mt2;
        b_1_hhbtag = *rv_b_1_hhbtag;
        b_2_hhbtag = *rv_b_2_hhbtag;
        vbf_1_hhbtag = *rv_vbf_1_hhbtag;
        vbf_2_hhbtag = *rv_vbf_2_hhbtag;
        b_1_cvsl = *rv_b_1_cvsl;
        b_2_cvsl = *rv_b_2_cvsl;
        vbf_1_cvsl = *rv_vbf_1_cvsl;
        vbf_2_cvsl = *rv_vbf_2_cvsl;
        b_1_cvsb = *rv_b_1_cvsb;
        b_2_cvsb = *rv_b_2_cvsb;
        vbf_1_cvsb = *rv_vbf_1_cvsb;
        vbf_2_cvsb = *rv_vbf_2_cvsb;

        // Load vectors
        pep_svfit.SetCoordinates(*rv_svfit_pT, *rv_svfit_eta, *rv_svfit_phi, *rv_svfit_mass);
        pep_l_1.SetCoordinates(*rv_l_1_pT, *rv_l_1_eta, *rv_l_1_phi, *rv_l_1_e);
        pep_l_2.SetCoordinates(*rv_l_2_pT, *rv_l_2_eta, *rv_l_2_phi, *rv_l_2_e);
        pep_met.SetCoordinates(*rv_met_pT, 0, *rv_met_phi, 0);
        pep_b_1.SetCoordinates(*rv_b_1_pT, *rv_b_1_eta, *rv_b_1_phi, *rv_b_1_e);
        pep_b_2.SetCoordinates(*rv_b_2_pT, *rv_b_2_eta, *rv_b_2_phi, *rv_b_2_e);
        pep_vbf_1.SetCoordinates(*rv_vbf_1_pT, *rv_vbf_1_eta, *rv_vbf_1_phi, *rv_vbf_1_e);
        pep_vbf_2.SetCoordinates(*rv_vbf_2_pT, *rv_vbf_2_eta, *rv_vbf_2_phi, *rv_vbf_2_e);

        //svfit.SetCoordinates(pep_svfit.Px(), pep_svfit.Py(), pep_svfit.Pz(), pep_svfit.M());
        //l_1.SetCoordinates(pep_l_1.Px(), pep_l_1.Py(), pep_l_1.Pz(), pep_l_1.M());
        //l_2.SetCoordinates(pep_l_2.Px(), pep_l_2.Py(), pep_l_2.Pz(), pep_l_2.M());
        //met.SetCoordinates(pep_met.Px(), pep_met.Py(), 0, 0);
        //b_1.SetCoordinates(pep_b_1.Px(), pep_b_1.Py(), pep_b_1.Pz(), pep_b_1.M());
        //b_2.SetCoordinates(pep_b_2.Px(), pep_b_2.Py(), pep_b_2.Pz(), pep_b_2.M());
        //vbf_1.SetCoordinates(pep_vbf_1.Px(), pep_vbf_1.Py(), pep_vbf_1.Pz(), pep_vbf_1.M());
        //vbf_2.SetCoordinates(pep_vbf_2.Px(), pep_vbf_2.Py(), pep_vbf_2.Pz(), pep_vbf_2.M());

        // Convergence
        svfit_conv = *rv_svfit_mass > 0;
        hh_kinfit_conv = kinfit_chi2 > 0;

        _evt_proc->process_to_vec(feat_vals, pep_b_1, pep_b_2, pep_l_1, pep_l_2, pep_met, pep_svfit, pep_vbf_1, pep_vbf_2, kinfit_mass, kinfit_chi2, mt2, boosted, b_1_csv, b_2_csv,
                                  e_channel, e_year, res_mass, spin, klambda, n_vbf, svfit_conv, hh_kinfit_conv, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag,
                                  vbf_2_hhbtag, b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb, cv, c2v, c3, true);
        tauH_mass = *rv_tauH_mass;
        tauH_pt = *rv_tauH_pt;
        tauH_eta = *rv_tauH_eta;
        tauH_phi = *rv_tauH_phi;
        tauH_e = *rv_tauH_e;

        bH_mass = *rv_bH_mass;
        bH_pt = *rv_bH_pt;
        bH_eta = *rv_bH_eta;
        bH_phi = *rv_bH_phi;
        bH_e = *rv_bH_e;

        HH_mass = *rv_HH_mass;
        HH_pt = *rv_HH_pt;
        HH_eta = *rv_HH_eta;
        HH_phi = *rv_HH_phi;
        HH_e = *rv_HH_e;
        if (evt % 2 == 0)
        {
            data_even->Fill();
        }
        else
        {
            data_odd->Fill();
        }

        if (n_events > 0 && n_saved_events >= n_events)
        {
            std::cout << "Exiting after " << n_saved_events << " events.\n";
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

std::map<unsigned, std::string> FileLooper::build_dataset_id_map(TFile *in_file)
{
    TTreeReader aux_reader("aux", in_file);
    TTreeReaderValue<std::vector<std::string>> rv_dataset_names(aux_reader, "dataset_names");
    TTreeReaderValue<std::vector<unsigned>> rv_dataset_hashes(aux_reader, "dataset_hashes");

    std::vector<unsigned> ids;
    std::vector<std::string> names;
    while (aux_reader.Next())
    {
        ids = *rv_dataset_hashes;
        names = *rv_dataset_names;
    }
    std::map<unsigned, std::string> id2name;
    for (unsigned int i = 0; i < ids.size(); i++)
        id2name[ids[i]] = names[i];
    return id2name;
}

std::map<unsigned, std::string> FileLooper::build_region_id_map(TFile *in_file)
{
    TTreeReader aux_reader("aux", in_file);
    TTreeReaderValue<std::vector<std::string>> rv_region_names(aux_reader, "region_names");
    TTreeReaderValue<std::vector<unsigned>> rv_region_hashes(aux_reader, "region_hashes");

    std::vector<unsigned> ids;
    std::vector<std::string> names;
    while (aux_reader.Next())
    {
        ids = *rv_region_hashes;
        names = *rv_region_names;
    }
    std::map<unsigned, std::string> id2name;
    for (unsigned int i = 0; i < ids.size(); i++)
        id2name[ids[i]] = names[i];
    return id2name;
}

void FileLooper::_prep_file(TTree *tree, const std::vector<std::unique_ptr<float>> &feat_vals, double *weight, int *sample, int *region_id, int *jet_cat,
                            int *class_id, unsigned long long int *strat_key)
{
    /* Add branches to tree and set addresses for values */

    for (unsigned int i = 0; i < _n_feats; i++)
        tree->Branch(_feat_names[i].c_str(), feat_vals[i].get());
    tree->Branch("weight", weight);
    tree->Branch("sample", sample);
    tree->Branch("region", region_id);
    tree->Branch("jet_cat", jet_cat);
    tree->Branch("strat_key", strat_key);
}

Channel FileLooper::_get_channel(std::string channel)
{
    /* COnvert channel to enum */

    if (channel == "tauTau")
    {
        return Channel(tauTau);
    }
    else if (channel == "muTau")
    {
        return Channel(muTau);
    }
    else if (channel == "eTau")
    {
        return Channel(eTau);
    }
    throw std::invalid_argument("Invalid channel: options are tauTau, muTau, eTau");
    return Channel(tauTau);
}

Year FileLooper::_get_year(std::string year)
{
    /* Convert year to enum */

    if (year == "2016")
    {
        return Year(y16);
    }
    else if (year == "2017")
    {
        return Year(y17);
    }
    else if (year == "2018")
    {
        return Year(y18);
    }
    throw std::invalid_argument("Invalid year: options are 2016, 2017, 2018");
    return Year(y16);
}

int FileLooper::_jet_cat_lookup(const bool has_b_pair, const bool has_vbf_pair, const bool is_boosted, const int num_btag_Loose, const int num_btag_Medium)
{
    if (!has_b_pair)
        return -1;
    if (has_vbf_pair && num_btag_Loose >= 1)
        return 5; // 2j1b+_VBFL, 2j1b+_VBF, 2j1b+_VBFT
    if (!has_vbf_pair && is_boosted && num_btag_Loose >= 2)
        return 4; // 2j2Lb+B_noVBF
    if (!has_vbf_pair && num_btag_Medium >= 2)
        return 3; // 2j2b+R_noVBF
    if (!has_vbf_pair && num_btag_Loose >= 1)
        return 2; // 2j1bR_noVBF
    if (!has_vbf_pair && num_btag_Loose == 0)
        return 1; // 2j0bR_noVBF
    return 0;     // 2j
}

bool FileLooper::_apply_baseline(std::string channel, int c_event, int pairType, int nleps, int nbjetscand, int isLeptrigger)
{
    // channel-dependant baseline selection
    if (channel == "tauTau")
    {
        if (nleps != 0 || nbjetscand < 2 || pairType != 2 || isLeptrigger != 1)
        {
            c_event++;
            return 0;
        }
        else
            return 1;
    }
    else if (channel == "muTau")
    {
        if (nleps != 0 || nbjetscand < 2 || pairType != 0 || isLeptrigger != 1)
        {
            c_event++;
            return 0;
        }
        else
            return 1;
    }
    else if (channel == "eTau")
    {
        if (nleps != 0 || nbjetscand < 2 || pairType != 1 || isLeptrigger != 1)
        {
            c_event++;
            return 0;
        }
        else
            return 1;
    }
    else
    {
        std::cout << "Specified channel: " << channel << std::endl;
        throw std::invalid_argument("Channel should either be tauTau, muTau or eTau!");
    }
}

int FileLooper::_get_region(std::string channel, int isOS, float dau1_deepTauVsJet, float dau2_deepTauVsJet, float dau1_iso, float dau1_eleMVAiso)
{
    if (channel == "tauTau")
    {
        if (isOS != 0 && dau1_deepTauVsJet >= 5 && dau2_deepTauVsJet >= 5){
            //region = "SR"; // signal region: opposite sign, isolated taus
            return 0;
        }
        else if (isOS == 0 && dau1_deepTauVsJet >= 5 && dau2_deepTauVsJet >= 5){
            //region = "SStight"; // B region
            return 2;
        }
        else if (isOS != 0 && dau1_deepTauVsJet >= 5 && dau2_deepTauVsJet >= 1 && dau2_deepTauVsJet < 5){
            //region = "OSinviso"; //  # C region
            return 1;
        }
        else if (isOS == 0 && dau1_deepTauVsJet >= 5 && dau2_deepTauVsJet >= 1 && dau2_deepTauVsJet < 5){
            //region = "SSinviso"; //  # D region
            return 3;
        }
        else{
            //region = "";
            return -1;
        }
    }
    else if (channel == "muTau")
    {
        if (isOS != 0 && dau1_iso < 0.15 && dau2_deepTauVsJet >= 5){
            //region = "SR"; // signal region: opposite sign, isolated taus
            return 0;
        }
        else if (isOS == 0 && dau1_iso < 0.15 && dau2_deepTauVsJet >= 5){
            //region = "SStight"; // B region
            return 2;
        }
        else if (isOS != 0 && dau1_iso < 0.15 && dau2_deepTauVsJet >= 1 && dau2_deepTauVsJet < 5){
            //region = "OSinviso"; //  # C region
            return 1;
        }
        else if (isOS == 0 && dau1_iso < 0.15 && dau2_deepTauVsJet >= 1 && dau2_deepTauVsJet < 5){
            //region = "SSinviso"; //  # D region
            return 3;
        }
        else{
            //region = "";
            return -1;
        }
    }
    else if (channel == "eTau")
    {
        if (isOS != 0 && dau1_eleMVAiso == 1 && dau2_deepTauVsJet >= 5){
            //region = "SR"; // signal region: opposite sign, isolated taus
            return 0;
        }
        else if (isOS == 0 && dau1_eleMVAiso == 1 && dau2_deepTauVsJet >= 5){
            //region = "SStight"; // B region
            return 2;
        }
        else if (isOS != 0 && dau1_eleMVAiso == 1 && dau2_deepTauVsJet >= 1 && dau2_deepTauVsJet < 5){
            //region = "OSinviso"; //  # C region
            return 1;
        }
        else if (isOS == 0 && dau1_eleMVAiso == 1 && dau2_deepTauVsJet >= 1 && dau2_deepTauVsJet < 5){
            //region = "SSinviso"; //  # D region
            return 3;
        }
        else{
            //region = "";
            return -1;
        }
    }
    else
    {
        std::cout << "Specified channel: " << channel << std::endl;
        throw std::invalid_argument("Channel should either be tauTau, muTau or eTau!");
    }
}

int FileLooper::_region_lookup(const std::string &region)
{
    /* PISA             <----> LLR
    ------------------------------------
    OS_Isolated         <----> SR
    OS_AntiIsolated     <----> OSinviso
    SS_Isolated         <----> SStight
    SS_Antiisolated     <----> SSinviso
    SS_LooseIsolated    <----> SSrlx */
    if (region == "SR")
        return 0;
    if (region == "OSinviso")
        return 1;
    if (region == "SStight")
        return 2;
    if (region == "SSinviso")
        return 3;
    if (region == "SSrlx")
        return 4;
    if (region == "OSrlx")
        return 5;
    //throw std::invalid_argument("Unrecognised region: " + region);
    else{
        return -1;
    }
}

void FileLooper::_sample_lookup(std::string &sample, int &sample_id, Spin &spin, float &res_mass)
{
    // TODO Update this

    spin = nonres;
    res_mass = 125;
    // matches anything up to '_M-' or '_m-' or '_M' or '_m', then captures one or more digits in a group
    std::regex mass_regex(".*_[Mm]-?([0-9]+)");
    std::smatch match;

    if (sample.find("_ggF_") != std::string::npos || sample.find("GluGluTo") != std::string::npos )
    {
        if (sample.find("Radion") != std::string::npos)
        {
            spin = radion;
            if (std::regex_search(sample, match, mass_regex))
            {
                std::string number_str = match[1];
                if (!number_str.empty())
                {
                    res_mass = std::stoi(number_str);
                }
                else
                {
                    throw std::runtime_error("No mass captured in the sample!");
                }
            }
            else
            {
                throw std::runtime_error("No '_M-' or '_m' found in sample!");
            }
            if (res_mass <= 400)
            {
                sample_id = -15;
            }
            else if (res_mass <= 600)
            {
                sample_id = -16;
            }
            else
            {
                sample_id = -17;
            }
        }
        else if (sample.find("Graviton") != std::string::npos)
        {
            spin = graviton;
            if (std::regex_search(sample, match, mass_regex))
            {
                std::string number_str = match[1];
                if (!number_str.empty())
                {
                    res_mass = std::stoi(number_str);
                }
                else
                {
                    throw std::runtime_error("No mass captured in the sample!");
                }
            }
            else
            {
                throw std::runtime_error("No '_M-' or '_m' found in sample!");
            }
            if (res_mass <= 400)
            {
                sample_id = -18;
            }
            else if (res_mass <= 600)
            {
                sample_id = -19;
            }
            else
            {
                sample_id = -20;
            }
        }
        else
        {
            sample_id = -999;
        }
    }
    else if (sample.find("_VBF_") != std::string::npos || sample.find("VBFTo") != std::string::npos)
    {
        if (sample.find("Radion") != std::string::npos)
        {
            spin = radion;
            if (std::regex_search(sample, match, mass_regex))
            {
                std::string number_str = match[1];
                if (!number_str.empty())
                {
                    res_mass = std::stoi(number_str);
                }
                else
                {
                    throw std::runtime_error("No mass captured in the sample!");
                }
            }
            else
            {
                throw std::runtime_error("No '_M-' or '_m' found in sample!");
            }
            if (res_mass <= 400)
            {
                sample_id = -21;
            }
            else if (res_mass <= 600)
            {
                sample_id = -22;
            }
            else
            {
                sample_id = -23;
            }
        }
        else if (sample.find("Graviton") != std::string::npos)
        {
            spin = graviton;
            if (std::regex_search(sample, match, mass_regex))
            {
                std::string number_str = match[1];
                if (!number_str.empty())
                {
                    res_mass = std::stoi(number_str);
                }
                else
                {
                    throw std::runtime_error("No mass captured in the sample!");
                }
            }
            else
            {
                throw std::runtime_error("No '_M-' or '_m' found in sample!");
            }
            if (res_mass <= 400)
            {
                sample_id = -24;
            }
            else if (res_mass <= 600)
            {
                sample_id = -25;
            }
            else
            {
                sample_id = -26;
            }
        }
    }
    else if (sample.find("_Run") != std::string::npos)
    {
        sample_id = 0;
    }
    else if (sample.find("_TT_fully") != std::string::npos || sample.find("_TT_semi") != std::string::npos || sample.find( "TTTo2L2Nu" ) != std::string::npos || sample.find( "TTToSemiLeptonic" ) != std::string::npos || sample.find( "TTToHadronic" ) != std::string::npos)
    {
        sample_id = 1;
    }
    else if (sample.find("TTWJets") != std::string::npos || sample.find("TTZTo") != std::string::npos)
    {
        sample_id = 2;
    }
    else if (sample.find("TTWW") != std::string::npos || sample.find("TTWZ") != std::string::npos || sample.find("TTZZ") != std::string::npos)
    {
        sample_id = 3;
    }
    else if (sample.find("ttH") != std::string::npos)
    {
        sample_id = 4;
    }
    else if (sample.find("DY") != std::string::npos)
    {
        sample_id = 5;
    }
    else if (sample.find("WJets") != std::string::npos)
    {
        sample_id = 6;
    }
    else if (sample.find("_ggH") != std::string::npos || sample.find("VBFH") != std::string::npos || sample.find("GGHH") != std::string::npos || sample.find("GluGluHToTauTau") != std::string::npos)
    {
        sample_id = 7;
    }
    else if (sample.find("ZH") != std::string::npos)
    {
        sample_id = 8;
    }
    else if (sample.find("WminusH") != std::string::npos)
    {
        sample_id = 9;
    }
    else if (sample.find("WplusHT") != std::string::npos)
    {
        sample_id = 10;
    }
    else if (sample.find("EWK") != std::string::npos)
    {
        sample_id = 11;
    }
    else if (sample.find("_WW") != std::string::npos)
    {
        sample_id = 12;
    }
    else if (sample.find("_WZ") != std::string::npos)
    {
        sample_id = 13;
    }
    else if (sample.find("_ZH") != std::string::npos)
    {
        sample_id = 14;
    }
    else if (sample.find("_ZZ") != std::string::npos)
    {
        sample_id = 15;
    }
    else if (sample.find("_ST_") != std::string::npos || sample.find("ST_") != std::string::npos)
    {
        sample_id = 16;
    }
    else if (sample.find("WWW") != std::string::npos || sample.find("WWZ") != std::string::npos || sample.find("WZZ") != std::string::npos || sample.find("ZZZ") != std::string::npos )
    {
        sample_id = 17;
    }
    else
    {
        throw std::invalid_argument("Unrecognised sample: " + sample);
    }
}

int FileLooper::_sample2class_lookup(const int &sample)
{
    if (sample == -999)
        return -999; // Sample to reject
    if (sample < 0)
        return 1; // Signal
    if (sample == 0)
        return -1; // Collider data
    return 0;      // Background
}

bool FileLooper::_accept_evt(const int &region_id, const int &jet_cat, const int &class_id, const float &klambda,
                             const float &cv, const float &c2v, const float &c3)
{
    if (_only_kl1 && klambda != use_kl)
    {
        // std::cout << "Rejecting due to klambda = " << klambda << "\n";
        return false; // Only consider klambda at SM point
    }
    if (!_inc_data && class_id == -1)
    {
        // std::cout << "Rejecting due to class ID = " << class_id << "\n";
        return false; // Don't include data and event is data
    }
    if (!_inc_other_regions && region_id != 0)
    {
        // std::cout << "Rejecting due to region = " << region << "\n";
        return false; // Don't include other regions and event is not SS Iso
    }
    if (!_inc_all_jets && jet_cat == 0)
    {
        // std::cout << "Rejecting due to jet_cat = " << jet_cat << "\n";
        return false; // Only use inference category jets and event is non-inference category
    }
    if (jet_cat < 0)
    {
        // std::cout << "Rejecting due to jet_cat = " << jet_cat << "\n";
        return false; // Don't include experimental jet categories
    }
    if (_only_sm_vbf && (cv != 1 || c2v != 1 || c3 != 1))
    {
        // std::cout << "Rejecting due to cv = " << cv << " c2v = " << c2v << " c3 = " << c3 << "\n";
        return false; // Only consider SM VBF
    }
    if (class_id == -999)
    {
        // std::cout << "Rejecting due to class ID = -999\n";
        return false; // Reject sample
    }
    // std::cout << "Accepting\n";
    return true;
}

unsigned long long int FileLooper::_get_strat_key(const int &sample, const int &jet_cat, const Channel &channel, const Year &year, const int &region)
{
    unsigned long long int strat_key = std::pow(2, std::abs(sample)) *
                                       std::pow(3, jet_cat) *
                                       std::pow(5, (float)channel) *
                                       std::pow(7, (float)year) *
                                       std::pow(11, region);
    if (strat_key == 0)
    {
        std::cout << "sample " << sample << " jet_cat " << jet_cat << " channel " << channel << " year " << year << " region " << region << "\n";
        throw std::overflow_error("Strat key overflow\n");
    }
    return strat_key;
}
