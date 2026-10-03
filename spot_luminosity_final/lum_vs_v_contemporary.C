#include "analysis_common.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::set;
using std::string;
using std::vector;

namespace {

struct StoredVFit {
    double A = std::numeric_limits<double>::quiet_NaN();
    double B = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
};

struct PhaseKey {
    int group = 2;
    double T = 1e99;
    double h = 1e99;
    string name;
};

PhaseKey phase_key_Vcontemporary(const string& phase) {
    PhaseKey k;
    k.name = phase;
    if (phase == "before_annealing") {
        k.group = 0;
        k.T = -1e99;
        k.h = -1e99;
        return k;
    }

    double T = 0.0, h = 0.0;
    if (std::sscanf(phase.c_str(), "annealing_T=%lf_h=%lf", &T, &h) == 2) {
        k.group = 1;
        k.T = T;
        k.h = h;
    }
    return k;
}

bool phase_less_Vcontemporary(const string& a, const string& b) {
    const PhaseKey ka = phase_key_Vcontemporary(a);
    const PhaseKey kb = phase_key_Vcontemporary(b);
    if (ka.group != kb.group) return ka.group < kb.group;
    if (ka.T != kb.T) return ka.T < kb.T;
    if (ka.h != kb.h) return ka.h < kb.h;
    return ka.name < kb.name;
}

string pretty_phase_Vcontemporary(const string& phase) {
    if (phase == "before_annealing") return "Before annealing";
    double T = 0.0, h = 0.0;
    if (std::sscanf(phase.c_str(), "annealing_T=%lf_h=%lf", &T, &h) == 2)
        return Form("%g °C, %g h", T, h);
    return phase;
}

int phase_color_Vcontemporary(size_t i) {
    static const int colors[] = {
        kBlack, kRed+1, kOrange+7, kGreen+2, kCyan+2,
        kAzure+1, kBlue+1, kViolet+1, kMagenta+1, kGray+2,
        kPink+7, kTeal+3
    };
    return colors[i % (sizeof(colors)/sizeof(colors[0]))];
}

int phase_marker_Vcontemporary(size_t i) {
    static const int markers[] = {20,21,22,23,24,25,26,27,28,30,33,34};
    return markers[i % (sizeof(markers)/sizeof(markers[0]))];
}

int find_col_Vcontemporary(const vector<string>& header, const string& name) {
    for (int i=0; i<(int)header.size(); ++i)
        if (trim_copy(header[i]) == name) return i;
    return -1;
}

map<int,StoredVFit> read_stored_V_fit(const string& filename) {
    map<int,StoredVFit> out;
    std::ifstream fin(filename);
    if (!fin.is_open()) return out;

    string line;
    if (!std::getline(fin,line)) return out;
    const vector<string> header = split_csv_simple(line);

    const int cspot = find_col_Vcontemporary(header,"spot");
    const int cA = find_col_Vcontemporary(header,"A");
    const int cB = find_col_Vcontemporary(header,"B");
    const int cconv = find_col_Vcontemporary(header,"converged");
    const int cstatus = find_col_Vcontemporary(header,"fit_status");

    if (cspot<0 || cA<0 || cB<0) return out;

    while (std::getline(fin,line)) {
        if (trim_copy(line).empty()) continue;
        vector<string> f = split_csv_simple(line);
        if ((int)f.size() < (int)header.size()) f.resize(header.size(),"");

        try {
            const int spot = (int)std::llround(std::stod(f[cspot]));
            StoredVFit s;
            s.A = std::stod(f[cA]);
            s.B = std::stod(f[cB]);

            if (cconv>=0 && !trim_copy(f[cconv]).empty())
                s.converged = parse_bool_safe(f[cconv]);
            else if (cstatus>=0 && !trim_copy(f[cstatus]).empty())
                s.converged = (std::stoi(f[cstatus]) == 0);
            else
                s.converged = finite_number(s.A) && finite_number(s.B);

            out[spot] = s;
        } catch (...) {
            // Ignore malformed fit-result rows without stopping the full plot production.
        }
    }
    return out;
}

void draw_syst_brackets_Vcontemporary(const vector<double>& x,
                                      const vector<double>& y,
                                      const vector<double>& d,
                                      double xmin, double xmax,
                                      int color) {
    if (x.empty() || y.size()!=x.size() || d.size()!=x.size()) return;
    double span = xmax-xmin;
    if (!(span>0.0)) span=1.0;
    const double half_width = 0.009*span;

    for (size_t i=0; i<x.size(); ++i) {
        if (!finite_number(d[i]) || d[i]<=0.0) continue;
        const double ytop = y[i]+d[i];
        const double ybottom = y[i]-d[i];
        const double hook = 0.15*d[i];

        TLine* top = new TLine(x[i]-half_width,ytop,x[i]+half_width,ytop);
        top->SetLineColor(color); top->SetLineWidth(2); top->Draw("SAME");
        TLine* tl = new TLine(x[i]-half_width,ytop,x[i]-half_width,ytop-hook);
        tl->SetLineColor(color); tl->SetLineWidth(2); tl->Draw("SAME");
        TLine* tr = new TLine(x[i]+half_width,ytop,x[i]+half_width,ytop-hook);
        tr->SetLineColor(color); tr->SetLineWidth(2); tr->Draw("SAME");

        TLine* bot = new TLine(x[i]-half_width,ybottom,x[i]+half_width,ybottom);
        bot->SetLineColor(color); bot->SetLineWidth(2); bot->Draw("SAME");
        TLine* bl = new TLine(x[i]-half_width,ybottom,x[i]-half_width,ybottom+hook);
        bl->SetLineColor(color); bl->SetLineWidth(2); bl->Draw("SAME");
        TLine* br = new TLine(x[i]+half_width,ybottom,x[i]+half_width,ybottom+hook);
        br->SetLineColor(color); br->SetLineWidth(2); br->Draw("SAME");
    }
}

} // namespace

// ============================================================================
// Draw, for every global hotspot, all annealing phases on the same Lum-vs-Vover
// canvas. NO FIT IS PERFORMED HERE.
//
// Data points are read from the nominal R=20 all-phases master CSV.
// For each phase, A and B are read from the *_B_fit_results.csv produced
// previously by lum_vs_v_fit.C and are used only to reconstruct the already-
// obtained curve Lum(Vover) = A * Vover^B.
// ============================================================================
void lum_vs_v_contemporary(const char* all_phases_csv,
                           const char* analysis_base_dir,
                           const char* output_dir,
                           const char* prefix) {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const string master = all_phases_csv ? all_phases_csv : "";
    const string analysis_base = analysis_base_dir ? analysis_base_dir : "";
    const string outdir = output_dir ? output_dir : "";
    const string pref = prefix ? prefix : "";

    if (master.empty() || analysis_base.empty() || outdir.empty() || pref.empty()) {
        std::cerr << "ERROR: lum_vs_v_contemporary received an empty required argument." << std::endl;
        return;
    }

    gSystem->mkdir(outdir.c_str(), true);

    CsvTable table = read_analysis_csv(master,true);
    if (table.rows.empty()) {
        std::cerr << "WARNING: no rows found in " << master << std::endl;
        return;
    }

    set<string> phase_set;
    for (const auto& r : table.rows)
        if (!r.phase.empty()) phase_set.insert(r.phase);

    vector<string> phases(phase_set.begin(),phase_set.end());
    std::sort(phases.begin(),phases.end(),phase_less_Vcontemporary);

    // Keep exactly the genuinely detected points, as done by lum_vs_v_fit.C.
    map<int,map<string,vector<AnalysisRow>>> data;
    set<int> spot_ids;
    for (const string& ph : phases) {
        vector<AnalysisRow> selected = filter_phase_detected(table.rows,ph);
        for (const auto& r : selected) {
            data[r.spot][ph].push_back(r);
            spot_ids.insert(r.spot);
        }
    }

    // Read once the stored single-hotspot fit result for every phase.
    map<string,map<int,StoredVFit>> stored_fits;
    for (const string& ph : phases) {
        const string fitfile = analysis_base + "/" + ph + "/" +
                               pref + "_" + ph + "_B_fit_results.csv";
        stored_fits[ph] = read_stored_V_fit(fitfile);
        if (stored_fits[ph].empty())
            std::cerr << "WARNING: no stored V-fit results read from " << fitfile << std::endl;
    }

    int ncreated = 0;

    for (int spot : spot_ids) {
        double xmin =  1e99, xmax = -1e99;
        double ymin =  1e99, ymax = -1e99;
        bool have_points = false;
        double spot_x = 0.0, spot_y = 0.0;
        bool have_coordinates = false;

        for (const string& ph : phases) {
            auto itd = data.find(spot);
            if (itd==data.end()) continue;
            auto itp = itd->second.find(ph);
            if (itp==itd->second.end()) continue;

            for (const auto& r : itp->second) {
                have_points = true;
                if (!have_coordinates) {
                    spot_x = r.x; spot_y = r.y; have_coordinates = true;
                }
                xmin = std::min(xmin,r.v_fin);
                xmax = std::max(xmax,r.v_fin);
                const double e = std::max(std::fabs(r.error),std::fabs(r.deltaL));
                ymin = std::min(ymin,r.luminosity-e);
                ymax = std::max(ymax,r.luminosity+e);
            }
        }

        if (!have_points) continue;

        double xspan = xmax-xmin;
        if (!(xspan>0.0)) xspan=1.0;
        double frame_xmin = xmin-0.08*xspan;
        const double frame_xmax = xmax+0.08*xspan;
        if (frame_xmin<0.0) frame_xmin=0.0;

        double yspan = ymax-ymin;
        if (!(yspan>0.0)) yspan=std::max(1.0,std::fabs(ymax)*0.20);
        double frame_ymin = ymin-0.12*yspan;
        double frame_ymax = ymax+0.25*yspan;
        if (ymin>=0.0 && frame_ymin<0.0) frame_ymin=0.0;

        TCanvas* c = new TCanvas(Form("c_v_contemporary_spot_%d",spot),
                                 Form("Spot %d - contemporary Luminosity vs Vover",spot),
                                 1800,1100);
        c->SetLeftMargin(0.10);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.11);
        c->SetTopMargin(0.08);

        TH1* frame = c->DrawFrame(
            frame_xmin,frame_ymin,frame_xmax,frame_ymax,
            Form("%s - Spot %d (x = %.2f, y = %.2f): luminosity vs overvoltage for all annealing phases;Overvoltage (V);Luminosity",
                 pref.c_str(),spot,spot_x,spot_y)
        );
        frame->GetXaxis()->SetTitleSize(0.045);
        frame->GetYaxis()->SetTitleSize(0.045);
        frame->GetYaxis()->SetTitleOffset(1.05);

        TLegend* leg = new TLegend(0.12,0.73,0.88,0.91);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.022);
        leg->SetNColumns(3);

        vector<TGraphErrors*> graphs;
        vector<TF1*> funcs;

        for (size_t ip=0; ip<phases.size(); ++ip) {
            const string& ph = phases[ip];
            auto itd = data.find(spot);
            if (itd==data.end()) continue;
            auto itp = itd->second.find(ph);
            if (itp==itd->second.end() || itp->second.empty()) continue;

            vector<AnalysisRow> rows = itp->second;
            std::sort(rows.begin(),rows.end(),[](const AnalysisRow& a,const AnalysisRow& b){return a.v_fin<b.v_fin;});

            const int color = phase_color_Vcontemporary(ip);
            const int marker = phase_marker_Vcontemporary(ip);
            const int n = (int)rows.size();

            vector<double> xv(n),yv(n),exv(n,0.0),eyv(n),esyst(n);
            for (int i=0;i<n;++i) {
                xv[i]=rows[i].v_fin;
                yv[i]=rows[i].luminosity;
                eyv[i]=rows[i].error;
                esyst[i]=rows[i].deltaL;
            }

            TGraphErrors* gr = new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
            gr->SetName(Form("gr_v_contemporary_spot_%d_phase_%zu",spot,ip));
            gr->SetMarkerColor(color);
            gr->SetLineColor(color);
            gr->SetMarkerStyle(marker);
            gr->SetMarkerSize(1.25);
            gr->SetLineWidth(2);
            gr->Draw("PE SAME");

            draw_syst_brackets_Vcontemporary(xv,yv,esyst,frame_xmin,frame_xmax,color);
            gr->Draw("P SAME");
            graphs.push_back(gr);

            bool have_fit = false;
            StoredVFit sf;
            auto itfphase = stored_fits.find(ph);
            if (itfphase!=stored_fits.end()) {
                auto itf = itfphase->second.find(spot);
                if (itf!=itfphase->second.end()) {
                    sf=itf->second;
                    have_fit = sf.converged && finite_number(sf.A) && finite_number(sf.B);
                }
            }

            if (have_fit) {
                // No call to Fit(): use only the parameters already obtained by lum_vs_v_fit.C.
                const double fitmin = std::max(1e-6,frame_xmin);
                const double fitmax = frame_xmax;
                TF1* f = new TF1(Form("stored_power_spot_%d_phase_%zu",spot,ip),
                                 "[0]*TMath::Power(x,[1])",fitmin,fitmax);
                f->SetParameters(sf.A,sf.B);
                f->SetLineColor(color);
                f->SetLineWidth(3);
                f->SetLineStyle(1);
                f->Draw("SAME");
                funcs.push_back(f);
            }

            string label = pretty_phase_Vcontemporary(ph);
            if (have_fit) label += Form("  (B = %.3f)",sf.B);
            else          label += "  (fit unavailable)";
            leg->AddEntry(gr,label.c_str(),"pe");
        }

        leg->AddEntry((TObject*)0,"Solid curves: stored fits (no refit)","");
        leg->AddEntry((TObject*)0,"Vertical bars: stat.; brackets: syst.","");
        leg->Draw();

        c->Modified();
        c->Update();
        c->SaveAs(Form("%s/%s_lum_vs_v_contemporary_spot%d.png",outdir.c_str(),pref.c_str(),spot));
        c->SaveAs(Form("%s/%s_lum_vs_v_contemporary_spot%d.pdf",outdir.c_str(),pref.c_str(),spot));
        delete c;
        ++ncreated;
    }

    std::cout << "lum_vs_v_contemporary: created " << ncreated
              << " hotspot canvases in " << outdir << std::endl;
}
