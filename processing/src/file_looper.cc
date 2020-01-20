#include "cms_runII_data_proc/processing/interface/file_looper.hh"

FileLooper::FileLooper(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true) {
    _all = return_all;
    _requested = requested;
    _use_deep_csv = use_deep_csv;
}

~FileLooper() {}

bool FileLooper::loop_file(const& std::string, const& std::string, unsigned long int, bool) {
    return true;
}