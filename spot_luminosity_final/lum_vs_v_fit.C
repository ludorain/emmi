#include "analysis_common.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TH1D.h"
#include "TPaveText.h"
#include <limits>
#include <cmath>
#include <algorithm>

struct BSystematicPair {
    bool has16=false, has24=false;
    double B16=0.0, B24=0.0;
    bool conv16=false, conv24=false;
    bool has_delta=false;
    double deltaB=0.0;
};

struct VNominalFit {
    int spot=-1;
    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TF1* func=nullptr;

    double A=0.0, Aerr=0.0;
    double B=0.0, Berr=0.0;
    double chi2=0.0, chi2ndf=0.0, prob=0.0, edm=-1.0;
    int ndf=0, status=-999, covstatus=-999;

    bool fit_done=false, converged=false;
    bool has_param_syst=false;
    double deltaB=0.0;
};

struct VRadiusOverlayFit {
    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TF1* func=nullptr;
    double B=0.0, Berr=0.0;
    int status=-999;
    bool fit_done=false, converged=false;
};

// -----------------------------------------------------------------------------
// Append the statistics of the accepted individual B values for one annealing
// phase.  The file is shared by all phases and is reset by 5analysis.sh before
// a complete analysis run.
// -----------------------------------------------------------------------------
static void append_B_phase_summary(const string& filename,
                                   const string& phase,
                                   double B_mean,
                                   double B_mean_stat_error,
                                   double B_rms,
                                   size_t n_spots) {
    if (filename.empty()) return;

    bool write_header=true;
    {
        std::ifstream check(filename,std::ios::binary | std::ios::ate);
        if (check.is_open() && check.tellg()>0) write_header=false;
    }

    std::ofstream fout(filename,std::ios::app);
    if (!fout.is_open()) {
        std::cerr << "Warning: cannot open B-vs-phase summary CSV: "
                  << filename << std::endl;
        return;
    }

    if (write_header) {
        fout << "phase,B_mean,B_mean_stat_error,B_rms,n_spots\n";
    }

    fout << phase << ','
         << B_mean << ','
         << B_mean_stat_error << ','
         << B_rms << ','
         << n_spots << '\n';
}

// -----------------------------------------------------------------------------
// Read the R=16 and R=24 fit values produced by lum_vs_v_systematics.C.
// ROOT fit status == 0 is the convergence criterion;
// -----------------------------------------------------------------------------
static map<int,BSystematicPair> read_B_systematics(const string& filename) {
    map<int,BSystematicPair> out;
    std::ifstream fin(filename);
    if (!fin.is_open()) return out;

    string line;
    if (!std::getline(fin,line)) return out;
    auto header=split_csv_simple(line);
    map<string,int> col;
    for (int i=0;i<(int)header.size();++i) col[header[i]]=i;

    while (std::getline(fin,line)) {
        if (trim_copy(line).empty()) continue;
        auto f=split_csv_simple(line);
        if (f.size()<header.size()) f.resize(header.size(),"");

        try {
            int id=(int)std::llround(std::stod(f[col["spot"]]));
            BSystematicPair s;

            if (col.count("B_R16") && !f[col["B_R16"]].empty()) {
                s.B16=std::stod(f[col["B_R16"]]);
                s.has16=finite_number(s.B16);
            }
            if (col.count("B_R24") && !f[col["B_R24"]].empty()) {
                s.B24=std::stod(f[col["B_R24"]]);
                s.has24=finite_number(s.B24);
            }

            if (col.count("fit_status_R16") && !f[col["fit_status_R16"]].empty())
                s.conv16=(std::stoi(f[col["fit_status_R16"]])==0 && s.has16);
            else if (col.count("converged_R16"))
                s.conv16=parse_bool_safe(f[col["converged_R16"]]);

            if (col.count("fit_status_R24") && !f[col["fit_status_R24"]].empty())
                s.conv24=(std::stoi(f[col["fit_status_R24"]])==0 && s.has24);
            else if (col.count("converged_R24"))
                s.conv24=parse_bool_safe(f[col["converged_R24"]]);

            if (col.count("deltaB") && !f[col["deltaB"]].empty()) {
                s.deltaB=std::fabs(std::stod(f[col["deltaB"]]));
                s.has_delta=finite_number(s.deltaB);
            }
            out[id]=s;
        } catch (...) {
            // Skip malformed systematic rows without aborting the full analysis.
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
// Draw a symmetric systematic uncertainty as two black horizontal brackets.
// The bracket width is defined as a small fraction of the visible x range, so it
// remains only slightly wider than the experimental marker for any axis scale.
// -----------------------------------------------------------------------------
static void draw_horizontal_syst_brackets(const vector<double>& x,
                                           const vector<double>& y,
                                           const vector<double>& delta,
                                           double xmin,
                                           double xmax,
                                           int color=kBlack,
                                           int line_width=4) {
    if (x.empty() || y.size()!=x.size() || delta.size()!=x.size()) return;

    double span=xmax-xmin;
    if (!(span>0.0)) span=1.0;

    // Total cap width = 0.9% of the visible x range.
    // On the 1800-pixel canvas this is slightly wider than MarkerSize(1.5).
    const double half_width=0.009*span;

    for (size_t i=0;i<x.size();++i) {
        if (!finite_number(delta[i]) || delta[i]<=0.0) continue;

        const double ytop=y[i]+delta[i];
        const double ybottom=y[i]-delta[i];
        const double hook=0.20*delta[i];

        // Upper horizontal bracket with short downward hooks.
        TLine* top=new TLine(x[i]-half_width,ytop,x[i]+half_width,ytop);
        top->SetLineColor(color); top->SetLineWidth(line_width); top->Draw("SAME");
        TLine* top_l=new TLine(x[i]-half_width,ytop,x[i]-half_width,ytop-hook);
        top_l->SetLineColor(color); top_l->SetLineWidth(line_width); top_l->Draw("SAME");
        TLine* top_r=new TLine(x[i]+half_width,ytop,x[i]+half_width,ytop-hook);
        top_r->SetLineColor(color); top_r->SetLineWidth(line_width); top_r->Draw("SAME");

        // Lower horizontal bracket with short upward hooks.
        TLine* bottom=new TLine(x[i]-half_width,ybottom,x[i]+half_width,ybottom);
        bottom->SetLineColor(color); bottom->SetLineWidth(line_width); bottom->Draw("SAME");
        TLine* bottom_l=new TLine(x[i]-half_width,ybottom,x[i]-half_width,ybottom+hook);
        bottom_l->SetLineColor(color); bottom_l->SetLineWidth(line_width); bottom_l->Draw("SAME");
        TLine* bottom_r=new TLine(x[i]+half_width,ybottom,x[i]+half_width,ybottom+hook);
        bottom_r->SetLineColor(color); bottom_r->SetLineWidth(line_width); bottom_r->Draw("SAME");
    }
}

// -----------------------------------------------------------------------------
// Write an analysis CSV using the symmetric-systematic convention.
// -----------------------------------------------------------------------------
static void write_augmented_v_csv(const CsvTable& original,
                                  const string& phase,
                                  const map<int,VNominalFit>& fits,
                                  const string& output_csv) {
    std::ofstream fout(output_csv);
    if (!fout.is_open()) return;

    vector<int> kept_columns;
    bool input_has_deltaL=false;
    for (int i=0;i<(int)original.header.size();++i) {
        const string& name=original.header[i];
        if (name=="deltaL_plus" || name=="deltaL_minus") continue;
        if (name=="deltaL") input_has_deltaL=true;
        kept_columns.push_back(i);
    }

    bool first=true;
    for (int idx:kept_columns) {
        if (!first) fout << ',';
        fout << original.header[idx];
        first=false;
    }
    if (!input_has_deltaL) fout << ",deltaL";
    fout << ",B,B_stat_error,deltaB,A,A_stat_error,fit_chi2,fit_ndf,fit_chi2ndf,fit_status,covmatrix_status,edm,fit_converged\n";

    for (const auto& r:original.rows) {
        if (r.phase!=phase) continue;

        first=true;
        for (int idx:kept_columns) {
            if (!first) fout << ',';
            if (idx<(int)r.raw_fields.size()) fout << r.raw_fields[idx];
            first=false;
        }
        if (!input_has_deltaL) fout << ',' << r.deltaL;

        auto it=fits.find(r.spot);
        if (it==fits.end()) {
            fout << ",,,,,,,,,,,,0\n";
            continue;
        }

        const auto& f=it->second;
        fout << ',' << f.B
             << ',' << f.Berr
             << ',' << (f.has_param_syst ? f.deltaB : 0.0)
             << ',' << f.A
             << ',' << f.Aerr
             << ',' << f.chi2
             << ',' << f.ndf
             << ',' << f.chi2ndf
             << ',' << f.status
             << ',' << f.covstatus
             << ',' << f.edm
             << ',' << (f.converged ? 1 : 0)
             << '\n';
    }
}


// -----------------------------------------------------------------------------
// Build the optional R=16/R=24 overlays. The same power-law model and starting
// parameters used for the nominal R=20 fit are used here. "N" prevents ROOT
// from attaching extra fit-statistics boxes to the comparison graphs; the fit
// functions are drawn manually as dashed lines.
// -----------------------------------------------------------------------------
static map<int,VRadiusOverlayFit> build_v_radius_overlays(
    const string& filename,
    const string& phase,
    int color,
    int marker,
    const string& tag,
    double fit_xmin,
    double fit_xmax,
    double A_init,
    double B_init) {

    map<int,VRadiusOverlayFit> out;
    if (filename.empty()) return out;

    CsvTable table=read_analysis_csv(filename,true);
    auto selected=filter_phase_detected(table.rows,phase);
    if (selected.empty()) {
        std::cerr << "Warning: no comparison rows for phase " << phase
                  << " in " << filename << std::endl;
        return out;
    }

    map<int,vector<AnalysisRow>> spots;
    for (const auto& r:selected) spots[r.spot].push_back(r);

    for (auto& entry:spots) {
        int spot=entry.first;
        auto rows=entry.second;
        std::sort(rows.begin(),rows.end(),
                  [](const AnalysisRow& a,const AnalysisRow& b){return a.v_fin<b.v_fin;});

        VRadiusOverlayFit info;
        info.rows=rows;

        int n=(int)rows.size();
        if (n<2) {
            out[spot]=info;
            continue;
        }

        vector<double> xv(n),yv(n),exv(n,0.0),eyv(n);
        for (int i=0;i<n;++i) {
            xv[i]=rows[i].v_fin;
            yv[i]=rows[i].luminosity;
            eyv[i]=rows[i].error;
        }

        // Logarithmic prefit 
        double A_start = A_init;
        double B_start = B_init;

        vector<double> log_v;
        vector<double> log_l;
        vector<double> log_ev;
        vector<double> log_el;

        for (int i=0; i<n; ++i) {

            if (xv[i] <= 0.0 || yv[i] <= 0.0)
                continue;

            if (!finite_number(xv[i]) ||
                !finite_number(yv[i]) ||
                !finite_number(eyv[i]))
                continue;

            log_v.push_back(std::log(xv[i]));
            log_l.push_back(std::log(yv[i]));
            log_ev.push_back(0.0);
            log_el.push_back(std::fabs(eyv[i] / yv[i]));
        }

        if (log_v.size() >= 2) {

            TGraphErrors gr_log(
                (int)log_v.size(),
                log_v.data(),
                log_l.data(),
                log_ev.data(),
                log_el.data()
            );

            TF1 f_log(
                Form("f_log_%s_spot_%d",tag.c_str(),spot),
                "[0] + [1]*x",
                log_v.front(),
                log_v.back()
            );

            TFitResultPtr prefit = gr_log.Fit(&f_log,"QRSN0");

            if ((int)prefit == 0) {

                const double A_prefit = std::exp(f_log.GetParameter(0));

                const double B_prefit = f_log.GetParameter(1);

                if (finite_number(A_prefit) &&
                    finite_number(B_prefit) &&
                    A_prefit > 0.0) {

                    A_start = A_prefit;
                    B_start = B_prefit;
                }
            }
        }

        // Main non-linear fit

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
        info.graph->SetName(Form("gr_%s_spot_%d",tag.c_str(),spot));
        info.graph->SetMarkerColor(color);
        info.graph->SetLineColor(color);
        info.graph->SetMarkerStyle(marker);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetLineWidth(2);

        info.func=new TF1(Form("fit_%s_spot_%d",tag.c_str(),spot),
                          "[0]*TMath::Power(x,[1])",
                          fit_xmin,fit_xmax);
        info.func->SetParameters(A_start,B_start);
        info.func->SetParNames("A","B");
        info.func->SetLineColor(color);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(2);

        TFitResultPtr fit=info.graph->Fit(info.func,"QRSN");
        info.fit_done=true;
        info.status=(int)fit;
        info.B=info.func->GetParameter(1);
        info.Berr=info.func->GetParError(1);
        info.converged=(info.status==0 && finite_number(info.B));

        out[spot]=info;
    }

    return out;
}

void lum_vs_v_fit(const char* all_phases_csv,
                  const char* phase,
                  const char* systematic_values_csv,
                  const char* output_dir,
                  const char* prefix,
                  const char* r16_all_phases_csv = "",
                  const char* r24_all_phases_csv = "",
                  const char* B_summary_csv = "") {
    //gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    const double FIT_XMIN=0.0;
    const double FIT_XMAX=8.0;
    const double A_INIT=1.0;
    const double B_INIT=2.0;

    string ph=phase, outdir=output_dir, pref=prefix;
    const string summary_file=(B_summary_csv ? B_summary_csv : "");
    const string r16_file=(r16_all_phases_csv ? r16_all_phases_csv : "");
    const string r24_file=(r24_all_phases_csv ? r24_all_phases_csv : "");
    const bool comparison_mode=(!r16_file.empty() || !r24_file.empty());
    ensure_dir(outdir);

    CsvTable table=read_analysis_csv(all_phases_csv,true);
    auto selected=filter_phase_detected(table.rows,ph);
    if (selected.empty()) {
        std::cerr << "No rows for phase " << ph << std::endl;
        return;
    }

    auto sys=read_B_systematics(systematic_values_csv);

    map<int,vector<AnalysisRow>> spots;
    for (const auto& r:selected) spots[r.spot].push_back(r);

    map<int,VNominalFit> fits;

    // =====================================================================
    // FIT STAGE: each hotspot is fitted once and the result is reused below.
    // =====================================================================
    for (auto& entry:spots) {
        int spot=entry.first;
        auto rows=entry.second;
        std::sort(rows.begin(),rows.end(),
                  [](const AnalysisRow& a,const AnalysisRow& b){return a.v_fin<b.v_fin;});

        VNominalFit info;
        info.spot=spot;
        info.rows=rows;

        int n=(int)rows.size();
        if (n<2) {
            fits[spot]=info;
            continue;
        }

        vector<double> xv(n), yv(n), exv(n,0.0), eyv(n);
        for (int i=0;i<n;++i) {
            xv[i]=rows[i].v_fin;
            yv[i]=rows[i].luminosity;
            eyv[i]=rows[i].error;
        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
        info.graph->SetName(Form("gr_spot_%d",spot));
        info.graph->SetTitle(Form("%s - Spot %d: x = %.2f, y = %.2f, T = %.1f #circC, %s;Overvoltage (V);Luminosity",
                                pref.c_str(),spot,rows[0].x,rows[0].y,rows[0].T,ph.c_str()));
        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetLineWidth(2);


        // =====================================================================
        // LOGARITHMIC PREFIT
        //
        // Lum = A * V^B
        //
        // ln(Lum) = ln(A) + B * ln(V)
        //
        // The linear fit therefore provides:
        //   intercept = ln(A)
        //   slope     = B
        // =====================================================================

        double A_start = A_INIT;
        double B_start = B_INIT;

        vector<double> log_v;
        vector<double> log_l;
        vector<double> log_ev;
        vector<double> log_el;

        for (int i=0; i<n; ++i) {

            // Logarithms require positive V and luminosity.
            if (xv[i] <= 0.0 || yv[i] <= 0.0)
                continue;

            // Also reject non-finite quantities.
            if (!finite_number(xv[i]) ||
                !finite_number(yv[i]) ||
                !finite_number(eyv[i]))
                continue;

            log_v.push_back(std::log(xv[i]));
            log_l.push_back(std::log(yv[i]));

            // No uncertainty on V.
            log_ev.push_back(0.0);

            // d(ln L) = dL / L
            log_el.push_back(std::fabs(eyv[i] / yv[i]));
        }


        // At least two valid points are required for the linear prefit.
        if (log_v.size() >= 2) {

            TGraphErrors gr_log(
                (int)log_v.size(),
                log_v.data(),
                log_l.data(),
                log_ev.data(),
                log_el.data()
            );

            TF1 f_log(
                Form("f_log_spot_%d",spot),
                "[0] + [1]*x",
                log_v.front(),
                log_v.back()
            );

            TFitResultPtr prefit = gr_log.Fit(&f_log,"QRSN0");

            const int prefit_status = (int)prefit;

            if (prefit_status == 0) {

                const double lnA_prefit = f_log.GetParameter(0);
                const double B_prefit   = f_log.GetParameter(1);

                const double A_prefit = std::exp(lnA_prefit);

                if (finite_number(A_prefit) &&
                    finite_number(B_prefit) &&
                    A_prefit > 0.0) {

                    A_start = A_prefit;
                    B_start = B_prefit;
                }
            }
        }


        // =====================================================================
        // MAIN NON-LINEAR FIT
        // =====================================================================

        info.func=new TF1(Form("f_quad_spot_%d",spot),"[0]*TMath::Power(x,[1])",FIT_XMIN,FIT_XMAX);

        // Starting parameters obtained from the logarithmic prefit.
        // If the prefit failed, A_INIT and B_INIT are retained.
        info.func->SetParameters(A_start,B_start);

        info.func->SetParNames("A","B");

        TFitResultPtr fit=info.graph->Fit(info.func,"QRS0");

        info.fit_done=true;
        info.status=(int)fit;
        info.covstatus=fit->CovMatrixStatus();
        info.edm=fit->Edm();
        info.A=info.func->GetParameter(0);
        info.Aerr=info.func->GetParError(0);
        info.B=info.func->GetParameter(1);
        info.Berr=info.func->GetParError(1);
        info.chi2=info.func->GetChisquare();
        info.ndf=info.func->GetNDF();
        info.chi2ndf=(info.ndf>0) ? info.chi2/info.ndf : 0.0;
        info.prob=info.func->GetProb();

        // A point enters Canvas 9 when the ROOT minimisation converged and B is
        // finite. No chi2/ndf or covariance-status cut is applied here.
        info.converged=(info.status==0 && finite_number(info.B));

        // Symmetric systematic on B:
        // deltaB = max(|B24-B20|, |B20-B16|).
        if (sys.count(spot) && sys[spot].has_delta) {
            // deltaB is computed upstream from R16/R20/R24 fits performed on
            // the exact same common set of detected overvoltage points.
            info.deltaB=sys[spot].deltaB;
            info.has_param_syst=true;
        }

        fits[spot]=info;
    }

    map<int,VRadiusOverlayFit> fits16;
    map<int,VRadiusOverlayFit> fits24;
    if (!r16_file.empty()) {
        fits16=build_v_radius_overlays(r16_file,ph,kP6Blue,24,"R16",
                                       FIT_XMIN,FIT_XMAX,A_INIT,B_INIT);
    }
    if (!r24_file.empty()) {
        fits24=build_v_radius_overlays(r24_file,ph,kP6Red,25,"R24",
                                       FIT_XMIN,FIT_XMAX,A_INIT,B_INIT);
    }

    // =====================================================================
    // SINGLE-HOTSPOT CANVASES
    // =====================================================================
    string spotdir=outdir+"/lum_vs_v_all_spots";
    ensure_dir(spotdir);

    for (auto& entry:fits) {
        auto& info=entry.second;
        int id=entry.first;
        if (!info.fit_done || !info.graph || !info.func) continue;

        TCanvas* c=new TCanvas(Form("c_v_spot_%d",id),
                               Form("Spot %d - Luminosity vs overvoltage",id),
                               1600,1200);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.13);
        c->SetTopMargin(0.08);

        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetMarkerColor(kBlack);
        info.graph->SetLineColor(kBlack);
        info.graph->SetLineWidth(2);
        info.func->SetLineColor(kBlack);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(1);
        info.graph->GetXaxis()->SetTitleSize(0.045);
        info.graph->GetYaxis()->SetTitleSize(0.045);
        info.graph->GetYaxis()->SetTitleOffset(1.4);

        // Make sure R=20 (stat+syst) and the optional R=16/R=24
        // statistical error bars all fit inside the visible y range.
        double ymin=1e99, ymax=-1e99;
        for (const auto& r:info.rows) {
            const double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));
            ymin=std::min(ymin,r.luminosity-emax);
            ymax=std::max(ymax,r.luminosity+emax);
        }

        auto expand_overlay_range=[&](const map<int,VRadiusOverlayFit>& overlays) {
            auto it=overlays.find(id);
            if (it==overlays.end()) return;
            for (const auto& r:it->second.rows) {
                ymin=std::min(ymin,r.luminosity-std::fabs(r.error));
                ymax=std::max(ymax,r.luminosity+std::fabs(r.error));
            }
        };
        expand_overlay_range(fits16);
        expand_overlay_range(fits24);

        if (finite_number(ymin) && finite_number(ymax) && ymax>ymin) {
            double span=ymax-ymin;
            info.graph->SetMinimum(ymin-0.08*span);
            info.graph->SetMaximum(ymax+0.15*span);
        }

        info.graph->Draw("APZ");
        c->Update();

        vector<double> sx, sy, sd;
        for (const auto& r:info.rows) {
            sx.push_back(r.v_fin);
            sy.push_back(r.luminosity);
            sd.push_back(r.deltaL);
        }
        draw_horizontal_syst_brackets(sx,sy,sd,
                                      info.graph->GetXaxis()->GetXmin(),
                                      info.graph->GetXaxis()->GetXmax());

        // Experimental points and fit remain on top of the systematic brackets.
        info.graph->Draw("P SAME");
        info.func->Draw("SAME");

        VRadiusOverlayFit* ov16=nullptr;
        VRadiusOverlayFit* ov24=nullptr;

        auto it16=fits16.find(id);
        if (it16!=fits16.end() && it16->second.graph) {
            ov16=&it16->second;
            ov16->graph->Draw("PZ SAME");
            if (ov16->fit_done && ov16->func) ov16->func->Draw("SAME");
        }

        auto it24=fits24.find(id);
        if (it24!=fits24.end() && it24->second.graph) {
            ov24=&it24->second;
            ov24->graph->Draw("PZ SAME");
            if (ov24->fit_done && ov24->func) ov24->func->Draw("SAME");
        }

        TLine* syst_proxy=new TLine(0,0,1,0);
        syst_proxy->SetLineColor(kBlack);
        syst_proxy->SetLineWidth(4);

        TLegend* leg=nullptr;
        if (!comparison_mode) {
            // Original 5analysis legend.
            leg=new TLegend(0.12,0.73,0.48,0.91);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(info.graph,"Data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"Systematic uncertainty","l");
            leg->AddEntry(info.func,"Fit: Lum = A#upointV_{over}^{B}","l");
        } else {
            // Radius-comparison legend: R=20 remains the nominal dataset.
            leg=new TLegend(0.12,0.67,0.66,0.91);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.030);

            leg->AddEntry(info.graph,"R = 20 px data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"R = 20 px luminosity systematic uncertainty","l");

            string nominal_fit_label;
            if (info.has_param_syst) {
                nominal_fit_label=Form("R = 20 px fit: B = %.3f #pm %.3f (stat) #pm %.3f (syst)",
                                       info.B,info.Berr,info.deltaB);
            } else {
                nominal_fit_label=Form("R = 20 px fit: B = %.3f #pm %.3f",
                                       info.B,info.Berr);
            }
            leg->AddEntry(info.func,nominal_fit_label.c_str(),"l");

            if (ov16) {
                leg->AddEntry(ov16->graph,"R = 16 px data (statistical uncertainty)","lep");
                if (ov16->fit_done && ov16->func) {
                    string label16=Form("R = 16 px fit: B = %.3f #pm %.3f",
                                        ov16->B,ov16->Berr);
                    leg->AddEntry(ov16->func,label16.c_str(),"l");
                }
            }

            if (ov24) {
                leg->AddEntry(ov24->graph,"R = 24 px data (statistical uncertainty)","lep");
                if (ov24->fit_done && ov24->func) {
                    string label24=Form("R = 24 px fit: B = %.3f #pm %.3f",
                                        ov24->B,ov24->Berr);
                    leg->AddEntry(ov24->func,label24.c_str(),"l");
                }
            }
        }
        leg->Draw();

        c->Update();
        
        if (!comparison_mode) {
            TPaveStats* stats=(TPaveStats*)info.graph->FindObject("stats");
            if (!stats) {
                stats=new TPaveStats(0.12, 0.51, 0.48, 0.71, "brNDC");
                stats->SetName(Form("fit_stats_spot_%d",id));
                stats->AddText(Form("A = %.4g #pm %.3g",info.A,info.Aerr));
                stats->AddText(Form("B = %.4g #pm %.3g",info.B,info.Berr));
                stats->AddText(Form("#chi^{2}/ndf = %.3g / %d",info.chi2,info.ndf));
                stats->SetTextSize(0.020);
                stats->SetBorderSize(1);
                stats->SetFillStyle(0);
                stats->Draw(); }
        }

        c->Modified();
        c->Update();
        c->SaveAs(Form("%s/%s_lum_vs_v_spot%d.png",spotdir.c_str(),pref.c_str(),id));
        c->SaveAs(Form("%s/%s_lum_vs_v_spot%d.pdf",spotdir.c_str(),pref.c_str(),id));
        delete c;
    }

    if (comparison_mode) {
        // The comparison pipeline is intended only to redraw the individual
        // luminosity curves with R=16/20/24 overlays. The nominal B-vs-spot
        // and CSV products remain the responsibility of 5analysis.sh.
        return;
    }

    // =====================================================================
    // CANVAS 9: B vs global hotspot ID.
    //
    // There is NO simultaneous/common-B fit here.  We keep the B value from
    // each accepted individual hotspot fit and calculate simple distribution
    // statistics over exactly the same hotspot population shown in the plot.
    // =====================================================================
    vector<double> x, ex, B, Bstat, Bsyst;

    for (auto& kv:fits) {
        const auto& f=kv.second;

        // Keep the same quality selection for the B-vs-ID plot, the B
        // distribution histogram and the phase-summary statistics.
        if (!f.fit_done ||
            !f.converged ||
            !finite_number(f.B) ||
            !finite_number(f.Berr) ||
            f.ndf <= 0 ||
            f.chi2ndf >= 3.5) {
            continue;
        }

        x.push_back((double)kv.first);
        ex.push_back(0.0);
        B.push_back(f.B);
        Bstat.push_back(f.Berr);
        Bsyst.push_back(f.has_param_syst ? f.deltaB : 0.0);
    }

    if (!B.empty()) {
        // -----------------------------------------------------------------
        // Distribution statistics of the accepted individual B_i values.
        //
        // Mean:
        //     B_mean = (1/N) sum_i B_i
        //
        // Statistical error on the arithmetic mean, propagated from the
        // individual fit errors:
        //     sigma_mean = sqrt(sum_i sigma_i^2) / N
        //
        // RMS of the hotspot-to-hotspot distribution:
        //     RMS = sqrt[(1/N) sum_i (B_i-B_mean)^2]
        // -----------------------------------------------------------------
        const size_t nB=B.size();

        double B_mean=0.0;
        for (double b:B) B_mean+=b;
        B_mean/=(double)nB;

        double sum_stat2=0.0;
        for (double eb:Bstat) {
            if (finite_number(eb)) sum_stat2+=eb*eb;
        }
        const double B_mean_stat_error=std::sqrt(sum_stat2)/(double)nB;

        double sum_dev2=0.0;
        for (double b:B) {
            const double d=b-B_mean;
            sum_dev2+=d*d;
        }
        const double B_rms=std::sqrt(sum_dev2/(double)nB);

        std::cout << "B distribution for phase " << ph
                  << ": mean = " << B_mean
                  << " +/- " << B_mean_stat_error << " (stat on mean)"
                  << ", RMS = " << B_rms
                  << ", N = " << nB
                  << std::endl;

        // Append one row for this annealing phase.  5analysis.sh removes this
        // file before starting the phase loop, so rows accumulate only within
        // the current complete analysis run.
        append_B_phase_summary(summary_file,
                               ph,
                               B_mean,
                               B_mean_stat_error,
                               B_rms,
                               nB);

        const string sensor=sensor_from_prefix(pref);
        const double xmin=*std::min_element(x.begin(),x.end())-1.0;
        const double xmax=*std::max_element(x.begin(),x.end())+1.0;

        double ymin=1e99,ymax=-1e99;
        for (size_t i=0;i<B.size();++i) {
            const double emax=std::max(std::fabs(Bstat[i]),std::fabs(Bsyst[i]));
            ymin=std::min(ymin,B[i]-emax);
            ymax=std::max(ymax,B[i]+emax);
        }
        double yspan=ymax-ymin;
        if (!(yspan>0.0)) yspan=std::max(0.2,std::fabs(ymax)*0.2);

        // -------------------------------------------------------------
        // B vs global hotspot ID.  Two versions are retained, as in the
        // previous macro: full y range and focused 1<B<3 range.
        // -------------------------------------------------------------
        auto draw_B_canvas=[&](const string& tag,bool focus) {
            TCanvas* c9=new TCanvas(Form("c9_%s",tag.c_str()),
                                    "B vs global hotspot ID",1800,1000);

            TGraphErrors* grB=new TGraphErrors((int)B.size(),
                                               x.data(),B.data(),
                                               ex.data(),Bstat.data());
            grB->SetTitle(Form("Sensor %s - Power-law exponent B - %s;Global spot ID;B",
                               sensor.c_str(),ph.c_str()));
            grB->SetMarkerStyle(21);
            grB->SetMarkerSize(1.2);
            grB->SetMarkerColor(kBlack);
            grB->SetLineColor(kBlack);
            grB->SetLineWidth(2);
            grB->SetMinimum(focus ? -10.0 : ymin-0.15*yspan);
            grB->SetMaximum(focus ? 20.0 : ymax+0.35*yspan);
            grB->Draw("APZ");
            grB->GetXaxis()->SetLimits(xmin,xmax);

            // Systematic uncertainty on each individual B_i.
            draw_horizontal_syst_brackets(x,B,Bsyst,xmin,xmax,kBlack,4);
            grB->Draw("P SAME");

            // Reference B=2.
            TLine* line_B2=new TLine(xmin,2.0,xmax,2.0);
            line_B2->SetLineColor(kP6Red);
            line_B2->SetLineStyle(2);
            line_B2->SetLineWidth(2);
            line_B2->Draw("SAME");

            // NEW: dashed horizontal line at the arithmetic mean of the
            // accepted individual B values.
            TLine* mean_line=new TLine(xmin,B_mean,xmax,B_mean);
            mean_line->SetLineColor(kP10Gray);
            mean_line->SetLineStyle(2);
            mean_line->SetLineWidth(2);
            mean_line->Draw("SAME");

            TLine* syst_proxy=new TLine(0,0,1,0);
            syst_proxy->SetLineColor(kBlack);
            syst_proxy->SetLineWidth(4);

            TLegend* leg=new TLegend(.12,.67,.48,.91);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(grB,"B statistical uncertainty","lep");
            leg->AddEntry(syst_proxy,"B systematic uncertainty","l");

            string mean_label=Form("Mean B = %.4f #pm %.4f (stat)",
                                   B_mean,B_mean_stat_error);
            leg->AddEntry(mean_line,mean_label.c_str(),"l");
            leg->AddEntry(line_B2,"Reference: B = 2","l");
            leg->Draw();

            // Display the RMS separately from the statistical error on the mean.
            TPaveText* stats=new TPaveText(.64,.73,.91,.90,"NDC");
            stats->SetBorderSize(1);
            stats->SetFillStyle(0);
            stats->SetTextAlign(12);
            stats->SetTextSize(.028);
            stats->AddText(Form("N_{spots} = %d",(int)nB));
            stats->AddText(Form("Mean B = %.5g",B_mean));
            stats->AddText(Form("Stat. error on mean = %.4g",B_mean_stat_error));
            stats->AddText(Form("RMS = %.4g",B_rms));
            stats->Draw();

            c9->Modified();
            c9->Update();

            string base=outdir+"/"+pref+"_B_vs_spot"+(focus?"_focus_1_3":"");
            c9->SaveAs((base+".png").c_str());
            c9->SaveAs((base+".pdf").c_str());
            delete c9;
        };

        draw_B_canvas("full",false);
        draw_B_canvas("focus",true);

        /*
        // -------------------------------------------------------------
        // SECOND PLOT: binned distribution of the individual B_i values.
        // This is a diagnostic plot intended to show the shape of the B
        // distribution.  No Gaussian fit is imposed here.
        // -------------------------------------------------------------
        double Bmin=*std::min_element(B.begin(),B.end());
        double Bmax=*std::max_element(B.begin(),B.end());
        double Bspan=Bmax-Bmin;
        if (!(Bspan>0.0)) Bspan=std::max(0.2,std::fabs(Bmax)*0.2);

        const double hist_min=Bmin-0.10*Bspan;
        const double hist_max=Bmax+0.10*Bspan;
        const int nbins=std::max(8,
                         std::min(30,
                         (int)std::ceil(std::sqrt((double)nB))));

        TCanvas* ch=new TCanvas("c_B_distribution",
                                "Distribution of B values",1500,1000);
        ch->SetLeftMargin(0.12);
        ch->SetBottomMargin(0.12);
        ch->SetRightMargin(0.05);
        ch->SetTopMargin(0.08);

        TH1D* hB=new TH1D("h_B_distribution",
                          Form("Sensor %s - Distribution of power-law exponent B - %s;B;Number of hotspots",
                               sensor.c_str(),ph.c_str()),
                          nbins,hist_min,hist_max);
        hB->SetLineWidth(2);
        for (double b:B) hB->Fill(b);
        hB->Draw("HIST");

        const double hist_ymax=std::max(1.0,hB->GetMaximum()*1.12);
        hB->SetMaximum(hist_ymax);

        // Mean of the B_i distribution.
        TLine* hist_mean_line=new TLine(B_mean,0.0,B_mean,hist_ymax);
        hist_mean_line->SetLineColor(kP10Gray);
        hist_mean_line->SetLineStyle(2);
        hist_mean_line->SetLineWidth(3);
        hist_mean_line->Draw("SAME");

        // B=2 reference, retained also in the distribution plot.
        TLine* hist_B2_line=new TLine(2.0,0.0,2.0,hist_ymax);
        hist_B2_line->SetLineColor(kP6Red);
        hist_B2_line->SetLineStyle(2);
        hist_B2_line->SetLineWidth(2);
        hist_B2_line->Draw("SAME");

        TPaveText* hstats=new TPaveText(.62,.68,.91,.90,"NDC");
        hstats->SetBorderSize(1);
        hstats->SetFillStyle(0);
        hstats->SetTextAlign(12);
        hstats->SetTextSize(.030);
        hstats->AddText(Form("N_{spots} = %d",(int)nB));
        hstats->AddText(Form("Mean B = %.5g",B_mean));
        hstats->AddText(Form("Stat. error on mean = %.4g",B_mean_stat_error));
        hstats->AddText(Form("RMS = %.4g",B_rms));
        hstats->Draw();

        TLegend* hleg=new TLegend(.14,.73,.45,.89);
        hleg->SetBorderSize(0);
        hleg->SetFillStyle(0);
        hleg->AddEntry(hB,"Individual hotspot B values","l");
        hleg->AddEntry(hist_mean_line,"Mean B","l");
        hleg->AddEntry(hist_B2_line,"Reference: B = 2","l");
        hleg->Draw();

        ch->Modified();
        ch->Update();
        ch->SaveAs(Form("%s/%s_B_distribution.png",outdir.c_str(),pref.c_str()));
        ch->SaveAs(Form("%s/%s_B_distribution.pdf",outdir.c_str(),pref.c_str()));
        delete ch;

    } else {
        std::cerr << "Warning: no accepted B fit is available for the B-vs-ID statistics in phase "
                  << ph << std::endl;
    } */
        // -------------------------------------------------------------
        // SECOND PLOT: binned distribution of the individual B_i values.
        // This plot is meant to inspect the shape of the hotspot-to-hotspot
        // distribution; no Gaussian fit is imposed here.
        // -------------------------------------------------------------
        if (!B.empty()) {

            double hmin = *std::min_element(B.begin(),B.end());
            double hmax = *std::max_element(B.begin(),B.end());

            if (!(hmax > hmin)) {
                const double half = std::max(0.1,0.10*std::fabs(hmax));
                hmin -= half;
                hmax += half;
            } else {
                const double margin = 0.10*(hmax-hmin);
                hmin -= margin;
                hmax += margin;
            }

            // Use a finer binning than the previous sqrt(N) choice.
            // This keeps enough granularity to inspect the distribution shape.
            int nbins = std::max(
                15,
                (int)std::ceil(2.0*std::sqrt((double)B.size()))
            );
            nbins = std::min(nbins,60);

            TCanvas* ch = new TCanvas(
                "c_B_distribution",
                "B distribution",
                1900,900
            );

            // Reserve a right-side area for the summary box so that it does
            // not overlap the histogram.
            ch->SetRightMargin(0.28);

            TH1D* hB = new TH1D(
                "h_B_distribution",
                Form("Sensor %s - Distribution of fitted B - %s;B;Hotspots",
                     sensor.c_str(),ph.c_str()),
                nbins,hmin,hmax
            );

            hB->SetLineColor(kP6Blue);
            hB->SetFillColor(kP6Blue);
            hB->SetFillStyle(1001);
            hB->SetLineWidth(2);

            for (double B_i : B)
                hB->Fill(B_i);

            hB->Draw("HIST");

            const double ymax_hist = hB->GetMaximum();

            // ---------------------------------------------------------
            // Mean B.
            // ---------------------------------------------------------
            TLine* mean_hist_line = new TLine(
                B_mean,0.0,
                B_mean,1.05*ymax_hist
            );
            mean_hist_line->SetLineColor(kBlack);
            mean_hist_line->SetLineWidth(2);
            mean_hist_line->SetLineStyle(2);
            mean_hist_line->Draw("SAME");

            // ---------------------------------------------------------
            // Reference B = 2.
            // ---------------------------------------------------------
            TLine* B2_hist_line = new TLine(
                2.0,0.0,
                2.0,1.05*ymax_hist
            );
            B2_hist_line->SetLineColor(kP6Red);
            B2_hist_line->SetLineWidth(2);
            B2_hist_line->SetLineStyle(2);
            B2_hist_line->Draw("SAME");

            // ---------------------------------------------------------
            // Statistics box in the reserved right-side area.
            // ---------------------------------------------------------
            TPaveText* stats_box = new TPaveText(
                0.74,0.66,
                0.97,0.88,
                "NDC"
            );

            stats_box->SetBorderSize(0);
            stats_box->SetFillStyle(0);
            stats_box->SetTextAlign(12);

            stats_box->AddText(Form("Mean B = %.5f",B_mean));
            stats_box->AddText(
                Form("Stat. error on mean = %.5f",
                     B_mean_stat_error)
            );
            stats_box->AddText(Form("RMS = %.5f",B_rms));
            stats_box->AddText(Form("N = %zu",B.size()));

            stats_box->Draw();

            ch->Modified();
            ch->Update();

            ch->SaveAs(
                Form("%s/%s_B_distribution.png",
                     outdir.c_str(),pref.c_str())
            );
            ch->SaveAs(
                Form("%s/%s_B_distribution.pdf",
                     outdir.c_str(),pref.c_str())
            );

            delete ch;
        }

    } else {
        std::cerr
            << "Warning: no accepted B fit is available for the B-vs-ID "
               "statistics in phase "
            << ph << std::endl;
    }

    // =====================================================================
    // CSV OUTPUTS
    // =====================================================================
    string result_csv=outdir+"/"+pref+"_B_fit_results.csv";
    std::ofstream fout(result_csv);
    fout << "spot,x,y,phase,A,A_stat_error,B,B_stat_error,deltaB,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";

    for (auto& kv:fits) {
        const auto& f=kv.second;
        double xx=f.rows.empty()?0.0:f.rows[0].x;
        double yy=f.rows.empty()?0.0:f.rows[0].y;
        fout << kv.first << ',' << xx << ',' << yy << ',' << ph
             << ',' << f.A << ',' << f.Aerr
             << ',' << f.B << ',' << f.Berr
             << ',' << (f.has_param_syst ? f.deltaB : 0.0)
             << ',' << f.chi2 << ',' << f.ndf << ',' << f.chi2ndf
             << ',' << f.prob << ',' << f.status << ',' << f.covstatus
             << ',' << f.edm << ',' << (f.converged?1:0) << '\n';
    }
    fout.close();

    write_augmented_v_csv(table,ph,fits,outdir+"/"+pref+"_analysis.csv");
}
