#ifndef GOLD_ANALYSIS_COMMON_H
#define GOLD_ANALYSIS_COMMON_H

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"

using std::map;
using std::string;
using std::vector;

struct GoldRow {
    int spot = -1;
    double x = 0.0;
    double y = 0.0;
    double luminosity = 0.0;
    double error = 0.0;
    double T = 0.0;
    double v = 0.0;
    double v_fin = std::numeric_limits<double>::quiet_NaN();
    string dataset_key;
    string run_name;
    string run_number;
    bool detected = true;
    double deltaL = 0.0;
};

struct GoldTable {
    vector<string> header;
    map<string,int> col;
    vector<GoldRow> rows;
};

inline string gold_trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

inline vector<string> gold_split_csv(const string& line) {
    vector<string> out;
    string field;
    std::stringstream ss(line);
    while (std::getline(ss, field, ',')) out.push_back(gold_trim(field));
    return out;
}

inline bool gold_parse_bool(string s) {
    s = gold_trim(s);
    for (char& c : s) c = std::tolower(static_cast<unsigned char>(c));
    return s=="true" || s=="1" || s=="yes" || s=="y" ||
           s=="present" || s=="detected";
}

inline bool gold_finite(double x) { return std::isfinite(x); }

inline bool gold_has_col(const GoldTable& t, const string& name) {
    return t.col.find(name) != t.col.end();
}

inline GoldTable read_gold_csv(const string& filename, bool require_deltaL=false) {
    GoldTable t;
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

    t.header = gold_split_csv(line);
    for (int i=0; i<(int)t.header.size(); ++i) t.col[t.header[i]] = i;

    const vector<string> required = {
        "spot","x","y","luminosity","error","T","v","dataset_key","detected"
    };
    for (const auto& c : required) {
        if (!gold_has_col(t,c)) {
            std::cerr << "Error: missing required column '" << c
                      << "' in " << filename << std::endl;
            t.rows.clear();
            return t;
        }
    }
    if (require_deltaL && !gold_has_col(t,"deltaL")) {
        std::cerr << "Error: missing required column 'deltaL' in "
                  << filename << std::endl;
        t.rows.clear();
        return t;
    }

    while (std::getline(fin,line)) {
        if (gold_trim(line).empty()) continue;
        auto f = gold_split_csv(line);
        if (f.size() < t.header.size()) f.resize(t.header.size(), "");

        try {
            GoldRow r;
            r.spot = (int)std::llround(std::stod(f[t.col["spot"]]));
            r.x = std::stod(f[t.col["x"]]);
            r.y = std::stod(f[t.col["y"]]);
            r.luminosity = std::stod(f[t.col["luminosity"]]);
            r.error = std::fabs(std::stod(f[t.col["error"]]));
            r.T = std::stod(f[t.col["T"]]);
            r.v = std::stod(f[t.col["v"]]);
            r.dataset_key = f[t.col["dataset_key"]];
            r.detected = gold_parse_bool(f[t.col["detected"]]);

            if (gold_has_col(t,"v_fin") && !f[t.col["v_fin"]].empty())
                r.v_fin = std::stod(f[t.col["v_fin"]]);
            if (gold_has_col(t,"run_name")) r.run_name = f[t.col["run_name"]];
            if (gold_has_col(t,"run_number")) r.run_number = f[t.col["run_number"]];
            if (gold_has_col(t,"deltaL") && !f[t.col["deltaL"]].empty())
                r.deltaL = std::fabs(std::stod(f[t.col["deltaL"]]));

            t.rows.push_back(r);
        } catch (const std::exception& e) {
            std::cerr << "Warning: skipped malformed CSV row: " << e.what()
                      << "\n" << line << std::endl;
        }
    }
    return t;
}

inline vector<GoldRow> gold_select_dataset(const vector<GoldRow>& rows,
                                           const string& key,
                                           bool detected_only=true) {
    vector<GoldRow> out;
    for (const auto& r : rows) {
        if (r.dataset_key != key) continue;
        if (detected_only && !r.detected) continue;
        out.push_back(r);
    }
    return out;
}

inline map<int,vector<GoldRow>> gold_group_spots(const vector<GoldRow>& rows) {
    map<int,vector<GoldRow>> out;
    for (const auto& r : rows) out[r.spot].push_back(r);
    return out;
}

inline void gold_ensure_dir(const string& d) {
    if (!d.empty()) gSystem->mkdir(d.c_str(), kTRUE);
}

inline bool gold_file_exists(const string& p) {
    return !p.empty() && gSystem->AccessPathName(p.c_str()) == kFALSE;
}

inline string gold_safe_token(string s) {
    for (char& c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c=='_' || c=='=' || c=='-')) c='_';
    }
    return s;
}

inline void gold_save_canvas(TCanvas* c, const string& base) {
    c->SaveAs((base+".png").c_str());
    c->SaveAs((base+".pdf").c_str());
}

// Draw systematic uncertainty as short horizontal brackets, preserving the
// visual convention used in the previous analysis macros.
inline void gold_draw_syst_brackets(const vector<double>& x,
                                    const vector<double>& y,
                                    const vector<double>& delta,
                                    double xmin,
                                    double xmax,
                                    int color=kBlack,
                                    int line_width=4,
                                    double width_fraction=0.005) {
    if (x.empty() || y.size()!=x.size() || delta.size()!=x.size()) return;
    double span = xmax-xmin;
    if (!(span>0.0)) span=1.0;
    const double half_width = width_fraction*span;

    for (size_t i=0; i<x.size(); ++i) {
        if (!gold_finite(delta[i]) || delta[i]<=0.0) continue;
        const double top_y = y[i]+delta[i];
        const double bot_y = y[i]-delta[i];
        const double hook = 0.12*delta[i];

        TLine* top = new TLine(x[i]-half_width,top_y,x[i]+half_width,top_y);
        top->SetLineColor(color); top->SetLineWidth(line_width); top->Draw("SAME");
        TLine* tl = new TLine(x[i]-half_width,top_y,x[i]-half_width,top_y-hook);
        tl->SetLineColor(color); tl->SetLineWidth(line_width); tl->Draw("SAME");
        TLine* tr = new TLine(x[i]+half_width,top_y,x[i]+half_width,top_y-hook);
        tr->SetLineColor(color); tr->SetLineWidth(line_width); tr->Draw("SAME");

        TLine* bot = new TLine(x[i]-half_width,bot_y,x[i]+half_width,bot_y);
        bot->SetLineColor(color); bot->SetLineWidth(line_width); bot->Draw("SAME");
        TLine* bl = new TLine(x[i]-half_width,bot_y,x[i]-half_width,bot_y+hook);
        bl->SetLineColor(color); bl->SetLineWidth(line_width); bl->Draw("SAME");
        TLine* br = new TLine(x[i]+half_width,bot_y,x[i]+half_width,bot_y+hook);
        br->SetLineColor(color); br->SetLineWidth(line_width); br->Draw("SAME");
    }
}

inline const GoldRow* gold_find_vfin(const vector<GoldRow>& rows,double x,double tol=1e-9) {
    for (const auto& r:rows)
        if (gold_finite(r.v_fin) && std::fabs(r.v_fin-x)<=tol) return &r;
    return nullptr;
}

inline const GoldRow* gold_find_T(const vector<GoldRow>& rows,double x,double tol=1e-9) {
    for (const auto& r:rows) if (std::fabs(r.T-x)<=tol) return &r;
    return nullptr;
}

inline void gold_estimate_exp(const vector<GoldRow>& rows,double& A0,double& lambda0) {
    vector<GoldRow> pos;
    for (const auto& r:rows) if (r.luminosity>0.0) pos.push_back(r);
    std::sort(pos.begin(),pos.end(),[](const GoldRow&a,const GoldRow&b){return a.T<b.T;});
    if (pos.size()>=2 && std::fabs(pos.back().T-pos.front().T)>1e-12) {
        lambda0=(std::log(pos.back().luminosity)-std::log(pos.front().luminosity)) /
                (pos.back().T-pos.front().T);
        A0=std::exp(std::log(pos.front().luminosity)-lambda0*pos.front().T);
    } else if (pos.size()==1) {
        A0=pos.front().luminosity;
        lambda0=0.0;
    } else {
        A0=1.0;
        lambda0=0.0;
    }
}

#endif
