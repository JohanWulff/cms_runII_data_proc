#include "cms_runII_data_proc/processing/interface/file_looper.hh"
#include <iostream>
#include <string>


void show_help() {
    /* Show help for input arguments */

    std::cout << "-y : Year\n";
    std::cout << "-c : Channel\n";
    std::cout << "--sample : sample name\n";
    std::cout << "-n : # events, default = -1 (all)\n";
    std::cout << "-s : start event number, default = 0\n";
    std::cout << "-e : end event number, default = -1 (last event)\n";
    std::cout << "-i : input file";
    std::cout << "-o : output file";
}

std::map<std::string, std::string> get_options(int argc, char* argv[]) {
    /*Interpret input arguments*/

    std::map<std::string, std::string> options;
    options.insert(std::make_pair("-y", "")); // Year
    options.insert(std::make_pair("-c", "")); // Channel
    options.insert(std::make_pair("--sample", "")); // sample str. 
    options.insert(std::make_pair("--sum_w", "1.0")); // sample sum_w. 
    options.insert(std::make_pair("-n", "-1")); // # events
    options.insert(std::make_pair("-s", "0")); // start event number
    options.insert(std::make_pair("-e", "-1")); // end event number
    options.insert(std::make_pair("-i", "")); // input file name
    options.insert(std::make_pair("-o", "")); // output file name

    if (argc >= 2) { //Check if help was requested
        std::string option(argv[1]);
        if (option == "-h" || option == "--help") {
            show_help();
            options.clear();
            return options;
        }
    }

    for (int i = 1; i < argc; i = i+2) {
        std::string option(argv[i]);
        std::string argument(argv[i+1]);
        if (option == "-h" || option == "--help" || argument == "-h" || argument == "--help") { // Check if help was requested
            show_help();
            options.clear();
            return options;
        }
        options[option] = argument;
    }
    return options;
}

int main(int argc, char *argv[]) {
    std::map<std::string, std::string> options = get_options(argc, argv); // Parse arguments
    if (options.size() == 0) return 1;
    std::vector<std::string> requested={"sv_mt",
                                        "met_pT",
                                        "dR_b1_b2",
                                        "sv_E",
                                        "bjet1_HHbtag",
                                        "bjet2_HHbtag",
                                        "ll_mt",
                                        "sv_mass",
                                        "dR_l1_l2_x_sv_pT",
                                        "l_2_pT",
                                        "l_1_mt",
                                        "hh_kinfit_chi2",
                                        "hh_kinfit_m",
                                        "diH_mass_met",
                                        "h_bb_mass",
                                        "top_1_mass",
                                        "hh_pT",
                                        "dphi_hbb_met",
                                        "res_mass",
                                        "boosted",
                                        "channel",
                                        "is_vbf",
                                        "spin",
                                        "year",
                                        "dR_l1_l2",
                                        "deta_hbb_httvis",
                                        "met_px",
                                        "met_py",
                                        "dmet_resp_px",
                                        "dmet_resp_py",
                                        "dmet_reso_px",
                                        "dmet_reso_py",
                                        "dphi_l1_l2",
                                        "dphi_b1_b2",
                                        "deta_l1_l2",
                                        "deta_b1_b2",
                                        "dau1_px",
                                        "dau1_py",
                                        "dau1_pz",
                                        "l_1_E",
                                        "dau1_dxy",
                                        "dau1_dz",
                                        "dau1_iso",
                                        "dau2_px",
                                        "dau2_py",
                                        "dau2_pz",
                                        "l_2_E",
                                        "dau2_dxy",
                                        "dau2_dz",
                                        "dau2_iso",
                                        "met_cov00",
                                        "met_cov01",
                                        "met_cov11",
                                        "bjet1_px",
                                        "bjet1_py",
                                        "bjet1_pz",
                                        "b_1_E",
                                        "bjet1_bID_deepFlavor",
                                        "bjet1_cID_deepFlavor",
                                        "bjet1_pnet_bb",
                                        "bjet1_pnet_cc",
                                        "bjet1_pnet_b",
                                        "bjet1_pnet_c",
                                        "bjet1_pnet_g",
                                        "bjet1_pnet_uds",
                                        "bjet1_pnet_pu",
                                        "bjet1_pnet_undef",
                                        "bjet2_px",
                                        "bjet2_py",
                                        "bjet2_pz",
                                        "b_2_E",
                                        "bjet2_bID_deepFlavor",
                                        "bjet2_cID_deepFlavor",
                                        "bjet2_pnet_bb",
                                        "bjet2_pnet_cc",
                                        "bjet2_pnet_b",
                                        "bjet2_pnet_c",
                                        "bjet2_pnet_g",
                                        "bjet2_pnet_uds",
                                        "bjet2_pnet_pu",
                                        "bjet2_pnet_undef",
                                        "pairType",
                                        "dau1_decayMode",
                                        "dau2_decayMode",
                                        "dau1_charge",
                                        "dau2_charge"};
    FileLooper file_looper(false, requested, true, true, false, false, true, true);
    std::cout << "Using sumw of " << std::stof(options["--sum_w"]) << std::endl;
    bool ok = file_looper.loop_file(options["-i"], options["-o"], options["-y"],
                                    options["--sample"], std::stof(options["--sum_w"]), 
                                    std::stoi(options["-n"]), std::stoi(options["-s"]), std::stoi(options["-e"]));
    if (ok) std::cout << "File loop ran ok!\n";
    return 0;
}