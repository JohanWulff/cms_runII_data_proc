#include "cms_runII_data_proc/processing/interface/file_looper.hh"
#include <iostream>
#include <string>

std::string root_dir = "/eos/home-k/kandroso/cms-it-hh-bbtautau/anaTuples/2020-01-10/";
std::string out_dir = "/eos/user/g/gstrong/cms_runII_dnn/data/";


void show_help() {
    /* Show help for input arguments */

    std::cout << "-y : Year\n";
    std::cout << "-c : Channel\n";
    std::cout << "-n : # events, default = -1 (all)\n";
    std::cout << "-i : input dir, default " << root_dir << "data/set\n";
    std::cout << "-o : out dir, default = " << out_dir << "\n";
}

std::map<std::string, std::string> get_options(int argc, char* argv[]) {
    /*Interpret input arguments*/

    std::map<std::string, std::string> options;
    options.insert(std::make_pair("-y", "")); // Year
    options.insert(std::make_pair("-c", "")); // Channel
    options.insert(std::make_pair("-n", "-1")); // # events
    options.insert(std::make_pair("-i", root_dir)); // input dir name
    options.insert(std::make_pair("-o", out_dir)); // output name

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

    FileLooper file_looper;
    bool ok = file_looper.loop_file(options["-i"], options["-o"], options["-c"], options["-y"], options["-n"]);
    if (ok) std::cout << "File loop ran ok!";
    return 0;
}