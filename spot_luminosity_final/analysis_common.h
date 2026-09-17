#ifndef ANALYSIS_COMMON_H
#define ANALYSIS_COMMON_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <string>
#include <cctype>
#include <limits>

#include "TSystem.h"
#include "TString.h"

using std::string;
using std::vector;
using std::map;

struct AnalysisRow {
    int spot = -1;
    double x = 0.0, y = 0.0;
    double luminosity = 0.0, error = 0.0;
    double T = 0.0, v = 0.0, v_fin = 0.0;
    string phase;
    bool detected = true;
    // Symmetric luminosity systematic uncertainty:
    // deltaL = max(|L_R24-L_R20|, |L_R20-L_R16|).
    double deltaL = 0.0;
    vector<string> raw_fields;
};

struct CsvTable {
    vector<string> header;
    map<string,int> col;
    vector<AnalysisRow> rows;
};

inline string trim_copy(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == string::npos) return "";
    return s.substr(a, b-a+1);
}

inline vector<string> split_csv_simple(const string& line) {
    vector<string> out;
    string field;
    std::stringstream ss(line);
    while (std::getline(ss, field, ',')) out.push_back(trim_copy(field));
    return out;
}

inline bool parse_bool_safe(string s) {
    s = trim_copy(s);
    for (char& c : s) c = std::tolower(static_cast<unsigned char>(c));
    return (s=="true" || s=="1" || s=="yes" || s=="y");
}

inline bool finite_number(double x) { return std::isfinite(x); }

inline bool has_col(const CsvTable& t, const string& name) {
    return t.col.find(name) != t.col.end();
}

inline CsvTable read_analysis_csv(const string& filename, bool require_vfin=false) {
    CsvTable t;
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return t;
    }

    string line;
    if (!std::getline(fin, line)) {
        std::cerr << "Error: empty CSV " << filename << std::endl;
        return t;
    }
    t.header = split_csv_simple(line);
    for (int i=0; i<(int)t.header.size(); ++i) t.col[t.header[i]] = i;

    vector<string> required = {"spot","x","y","luminosity","error","T","v","phase"};
    for (const auto& c : required) {
        if (!has_col(t,c)) {
            std::cerr << "Error: missing required column '" << c << "' in " << filename << std::endl;
            t.rows.clear();
            return t;
        }
    }
    if (require_vfin && !has_col(t,"v_fin")) {
        std::cerr << "Warning: v_fin is missing in " << filename
                  << ". The macro will use the absolute voltage column v as x variable." << std::endl;
    }

    while (std::getline(fin,line)) {
        if (trim_copy(line).empty()) continue;
        vector<string> f = split_csv_simple(line);
        if (f.size() < t.header.size()) f.resize(t.header.size(), "");
        try {
            AnalysisRow r;
            r.raw_fields = f;
            r.spot = (int)std::llround(std::stod(f[t.col["spot"]]));
            r.x = std::stod(f[t.col["x"]]);
            r.y = std::stod(f[t.col["y"]]);
            r.luminosity = std::stod(f[t.col["luminosity"]]);
            r.error = std::stod(f[t.col["error"]]);
            r.T = std::stod(f[t.col["T"]]);
            r.v = std::stod(f[t.col["v"]]);
            r.phase = f[t.col["phase"]];
            r.v_fin = has_col(t,"v_fin") && !f[t.col["v_fin"]].empty() ? std::stod(f[t.col["v_fin"]]) : r.v;
            if (has_col(t,"detection") && !f[t.col["detection"]].empty()) r.detected = parse_bool_safe(f[t.col["detection"]]);
            else if (has_col(t,"detected") && !f[t.col["detected"]].empty()) r.detected = parse_bool_safe(f[t.col["detected"]]);
            // New pipeline convention: a single symmetric luminosity systematic.
            if (has_col(t,"deltaL") && !f[t.col["deltaL"]].empty()) {
                r.deltaL = std::fabs(std::stod(f[t.col["deltaL"]]));
            }
            // Backward-compatible fallback for old CSV files. New outputs must use deltaL only.
            else {
                double legacy_plus = 0.0, legacy_minus = 0.0;
                if (has_col(t,"deltaL_plus") && !f[t.col["deltaL_plus"]].empty())
                    legacy_plus = std::fabs(std::stod(f[t.col["deltaL_plus"]]));
                if (has_col(t,"deltaL_minus") && !f[t.col["deltaL_minus"]].empty())
                    legacy_minus = std::fabs(std::stod(f[t.col["deltaL_minus"]]));
                r.deltaL = std::max(legacy_plus, legacy_minus);
            }
            t.rows.push_back(r);
        } catch (const std::exception& e) {
            std::cerr << "Warning: skipped malformed CSV row: " << e.what() << "\n" << line << std::endl;
        }
    }
    return t;
}

inline vector<AnalysisRow> filter_phase(const vector<AnalysisRow>& in, const string& phase) {
    vector<AnalysisRow> out;
    for (const auto& r : in) if (r.phase == phase) out.push_back(r);
    return out;
}

inline string safe_token(string s) {
    for (char& c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c=='_' || c=='=')) c = '_';
    }
    return s;
}

inline void ensure_dir(const string& d) {
    if (!d.empty()) gSystem->mkdir(d.c_str(), kTRUE);
}

inline bool file_exists(const string& p) {
    return gSystem->AccessPathName(p.c_str()) == kFALSE;
}

inline double positive_or_zero(double x) { return x > 0.0 ? x : 0.0; }

inline vector<string> standard_phases() {
    return {
        "before_annealing",
        "annealing_T=75_h=5",
        "annealing_T=75_h=25",
        "annealing_T=100_h=5",
        "annealing_T=100_h=25",
        "annealing_T=125_h=5",
        "annealing_T=125_h=25",
        "annealing_T=150_h=5",
        "annealing_T=150_h=25"
    };
}

#endif
