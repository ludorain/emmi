//root lum_vs_T_fit.C()

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
#include "TLine.h"
#include "TH1D.h"
#include "TPaveText.h"

struct LambdaSystematicPair {

    bool has16=false,has24=false;
    double l16=0.0,l24=0.0;
    bool conv16=false,conv24=false;
    bool has_delta=false;
    double deltaLambda=0.0;

};

struct TNominalFit {

    int spot=-1;

    vector<AnalysisRow> rows;

    TGraphErrors* graph=nullptr;

    TGraphErrors* syst_graph=nullptr;

    TF1* func=nullptr;

    double A=0.0,Aerr=0.0,lambda=0.0,lambdaerr=0.0;

    double chi2=0.0,chi2ndf=0.0,prob=0.0,edm=-1.0;

    int ndf=0,status=-999,covstatus=-999;

    bool fit_done=false,converged=false;

    bool has_param_syst=false;

    double deltaLambda=0.0;

};

struct TRadiusOverlayFit {

    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TF1* func=nullptr;
    double lambda=0.0,lambdaerr=0.0;
    int status=-999;
    bool fit_done=false,converged=false;

};

static void estimate_exp_parameters_nom(const vector<AnalysisRow>& rows, double& A0, double& lambda0) {

    vector<AnalysisRow> pos;

    // Mantieni solo i punti con luminosità positiva,

    // perché dobbiamo calcolare log(L).

    for (const auto& r : rows) {

        if (r.luminosity > 0.0 && finite_number(r.luminosity)) {

            pos.push_back(r);

        }

    }

    std::sort(pos.begin(), pos.end(),

              [](const AnalysisRow& a, const AnalysisRow& b) {

                  return a.T < b.T;

              });

    if (pos.size() >= 2 &&

        std::fabs(pos.back().T - pos.front().T) > 1e-12) {

        // -------------------------------------------------------------

        // Pre-parametrizzazione logaritmica:

        //

        //     L(T) = A exp(lambda T)

        //

        // diventa

        //

        //     ln(L) = ln(A) + lambda T

        //

        // quindi:

        //     p0 = ln(A)

        //     p1 = lambda

        // -------------------------------------------------------------

        TGraph gr_log;

        for (size_t i = 0; i < pos.size(); ++i) {

            gr_log.SetPoint(i,

                            pos[i].T,

                            std::log(pos[i].luminosity));

        }

        TF1 f_lin("f_lin_init",

                  "pol1",

                  pos.front().T,

                  pos.back().T);

        TFitResultPtr fit_lin = gr_log.Fit(&f_lin, "Q0SN");

        const int status = (int)fit_lin;

        if (status == 0 &&

            finite_number(f_lin.GetParameter(0)) &&

            finite_number(f_lin.GetParameter(1))) {

            lambda0 = f_lin.GetParameter(1);

            A0 = std::exp(f_lin.GetParameter(0));

        } else {

            // Fallback alla vecchia inizializzazione

            // usando il primo e l'ultimo punto.

            lambda0 =

                (std::log(pos.back().luminosity)

                 - std::log(pos.front().luminosity))

                /

                (pos.back().T - pos.front().T);

            A0 =

                std::exp(std::log(pos.front().luminosity)

                         - lambda0 * pos.front().T);

        }

    } else if (pos.size() == 1) {

        A0 = pos[0].luminosity;

        lambda0 = 0.0;

    } else {

        A0 = 2.0;

        lambda0 = 0.06;

    }

}

// -----------------------------------------------------------------------------

// Build optional R=16/R=24 temperature overlays. They use the same exponential

// model and initial-parameter estimate as the nominal R=20 analysis. "N" keeps

// ROOT from creating extra fit-statistics boxes for the comparison datasets.

// -----------------------------------------------------------------------------

static map<int,TRadiusOverlayFit> build_T_radius_overlays(const string& filename, const string& phase, int color, int marker, const string& tag) {

    map<int,TRadiusOverlayFit> out;

    if (filename.empty()) return out;

    CsvTable table=read_analysis_csv(filename,false);

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

                  [](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});

        TRadiusOverlayFit info;

        info.rows=rows;

        int n=(int)rows.size();

        if (n<3) {

            out[spot]=info;

            continue;

        }

        double Tmin=rows.front().T,Tmax=rows.back().T,dT=Tmax-Tmin;

        if (dT<=0) dT=1.0;

        double fitmin=Tmin-0.10*dT,fitmax=Tmax+0.10*dT;

        vector<double> xv(n),yv(n),exv(n,0.0),eyv(n);

        for (int i=0;i<n;++i) {

            xv[i]=rows[i].T;

            yv[i]=rows[i].luminosity;

            eyv[i]=rows[i].error;

        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());

        info.graph->SetName(Form("gr_%s_spot_%d",tag.c_str(),spot));

        info.graph->SetMarkerColor(color);

        info.graph->SetLineColor(color);

        info.graph->SetMarkerStyle(marker);

        info.graph->SetMarkerSize(1.25);

        info.graph->SetLineWidth(2);

        double A0,l0;

        estimate_exp_parameters_nom(rows,A0,l0);

        info.func=new TF1(Form("fit_%s_spot_%d",tag.c_str(),spot),

                          "[0]*exp([1]*x)",fitmin,fitmax);

        info.func->SetParNames("A","#lambda");

        info.func->SetParameters(A0,l0);

        info.func->SetLineColor(color);

        info.func->SetLineWidth(2);

        info.func->SetLineStyle(2);

        TFitResultPtr fit=info.graph->Fit(info.func,"QRSN");

        info.fit_done=true;

        info.status=(int)fit;

        info.lambda=info.func->GetParameter(1);

        info.lambdaerr=info.func->GetParError(1);

        info.converged=(info.status==0 && finite_number(info.lambda));

        out[spot]=info;

    }

    return out;

}

static map<int,LambdaSystematicPair> read_lambda_systematics(const string& filename) {

    map<int,LambdaSystematicPair> out;

    std::ifstream fin(filename);

    if (!fin.is_open()) return out;

    string line;

    if (!std::getline(fin,line)) return out;

    auto h=split_csv_simple(line);

    map<string,int> c;

    for (int i=0;i<(int)h.size();++i) c[h[i]]=i;

    while (std::getline(fin,line)) {

        if (trim_copy(line).empty()) continue;

        auto f=split_csv_simple(line);

        if (f.size()<h.size()) f.resize(h.size(),"");

        try {

            int id=(int)std::llround(std::stod(f[c["spot"]]));

            LambdaSystematicPair s;

            if (c.count("lambda_R16") && !f[c["lambda_R16"]].empty()) {

                s.l16=std::stod(f[c["lambda_R16"]]);

                s.has16=finite_number(s.l16);

            }

            if (c.count("lambda_R24") && !f[c["lambda_R24"]].empty()) {

                s.l24=std::stod(f[c["lambda_R24"]]);

                s.has24=finite_number(s.l24);

            }

            if (c.count("fit_status_R16") && !f[c["fit_status_R16"]].empty())

                s.conv16=(std::stoi(f[c["fit_status_R16"]])==0 && s.has16);

            else if (c.count("converged_R16"))

                s.conv16=parse_bool_safe(f[c["converged_R16"]]);

            if (c.count("fit_status_R24") && !f[c["fit_status_R24"]].empty())

                s.conv24=(std::stoi(f[c["fit_status_R24"]])==0 && s.has24);

            else if (c.count("converged_R24"))

                s.conv24=parse_bool_safe(f[c["converged_R24"]]);

            if (c.count("deltaLambda") && !f[c["deltaLambda"]].empty()) {

                s.deltaLambda=std::fabs(std::stod(f[c["deltaLambda"]]));

                s.has_delta=finite_number(s.deltaLambda);

            }

            out[id]=s;

        } catch (...) {

            // Skip malformed rows and continue with the remaining hotspots.

        }

    }

    return out;

}

// Write a phase CSV using only the new symmetric-systematic columns.

static void write_augmented_T_csv(const CsvTable& original, const string& phase, const map<int,TNominalFit>& fits, const string& output_csv) {

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

    fout << ",lambda,lambda_stat_error,deltaLambda,A,A_stat_error,fit_chi2,fit_ndf,fit_chi2ndf,fit_status,covmatrix_status,edm,fit_converged\n";

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

        fout << ',' << f.lambda

             << ',' << f.lambdaerr

             << ',' << (f.has_param_syst?f.deltaLambda:0.0)

             << ',' << f.A

             << ',' << f.Aerr

             << ',' << f.chi2

             << ',' << f.ndf

             << ',' << f.chi2ndf

             << ',' << f.status

             << ',' << f.covstatus

             << ',' << f.edm

             << ',' << (f.converged?1:0)

             << '\n';

    }

}

void lum_vs_T_fit(const char* all_phases_csv,
                  const char* phase,
                  const char* systematic_values_csv,
                  const char* output_dir,
                  const char* prefix,
                  const char* r16_all_phases_csv = "",
                  const char* r24_all_phases_csv = "",
                  const char* lambda_summary_csv = "") {

    gStyle->SetOptStat(0);

    gStyle->SetOptFit(0);

    // Statistical errors are drawn with option Z (no end caps).
    // Systematic errors are drawn with ROOT bracket option [].


    string ph=phase,outdir=output_dir,pref=prefix;

    const string r16_file=(r16_all_phases_csv ? r16_all_phases_csv : "");

    const string r24_file=(r24_all_phases_csv ? r24_all_phases_csv : "");

    const string lambda_summary_file=(lambda_summary_csv ? lambda_summary_csv : "");

    const bool comparison_mode=(!r16_file.empty() || !r24_file.empty());

    ensure_dir(outdir);

    CsvTable table=read_analysis_csv(all_phases_csv,false);

    auto selected=filter_phase_detected(table.rows,ph);

    if (selected.empty()) {

        std::cerr << "No rows for phase " << ph << std::endl;

        return;

    }

    auto sys=read_lambda_systematics(systematic_values_csv);

    map<int,vector<AnalysisRow>> spots;

    for (const auto& r:selected) spots[r.spot].push_back(r);

    map<int,TNominalFit> fits;

    for (auto& entry:spots) {

        int spot=entry.first;

        auto rows=entry.second;

        std::sort(rows.begin(),rows.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});

        TNominalFit info;

        info.spot=spot;

        info.rows=rows;

        int n=(int)rows.size();

        if (n<3) {

            fits[spot]=info;

            continue;

        }

        double Tmin=rows.front().T,Tmax=rows.back().T,dT=Tmax-Tmin;

        if (dT<=0) dT=1.0;

        double fitmin=Tmin-0.10*dT,fitmax=Tmax+0.10*dT;

        vector<double> xv(n),yv(n),exv(n,0.0),eyv(n),esyst(n);

        for (int i=0;i<n;++i) {

            xv[i]=rows[i].T;

            yv[i]=rows[i].luminosity;

            eyv[i]=rows[i].error;

            esyst[i]=rows[i].deltaL;

        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());

        info.graph->SetTitle(Form("%s - Spot %d: x = %.2f, y = %.2f, v = %.1f V, %s;T (#circC);Luminosity",

                                  pref.c_str(),spot,rows[0].x,rows[0].y,rows[0].v,ph.c_str()));

        info.graph->SetMarkerStyle(20);

        info.graph->SetMarkerSize(1.4);

        info.graph->SetLineWidth(2);

        info.syst_graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),esyst.data());

        info.syst_graph->SetLineColor(kBlack);

        info.syst_graph->SetLineWidth(2);

        info.syst_graph->SetMarkerSize(0);

        double A0,l0;

        estimate_exp_parameters_nom(rows,A0,l0);
        // Main exponential fit 

        info.func=new TF1(Form("fit_exp_spot_%d",spot),"[0]*exp([1]*x)",fitmin,fitmax);
        info.func->SetParNames("A","#lambda");
        info.func->SetParameters(A0,l0);
        info.func->SetLineColor(kBlack);
        info.func->SetLineWidth(2);
        TFitResultPtr fit=info.graph->Fit(info.func,"QRS0");
        info.fit_done=true;
        info.status=(int)fit;

        info.covstatus=fit->CovMatrixStatus();

        info.edm=fit->Edm();

        info.A=info.func->GetParameter(0);

        info.Aerr=info.func->GetParError(0);

        info.lambda=info.func->GetParameter(1);

        info.lambdaerr=info.func->GetParError(1);

        info.chi2=info.func->GetChisquare();

        info.ndf=info.func->GetNDF();

        info.chi2ndf=(info.ndf>0)?info.chi2/info.ndf:0.0;

        info.prob=info.func->GetProb();

        info.converged=(info.status==0 && finite_number(info.lambda));

        // Symmetric systematic on lambda:

        // deltaLambda=max(|lambda24-lambda20|, |lambda20-lambda16|).

        if (sys.count(spot) && sys[spot].has_delta) {

            // Computed upstream from R16/R20/R24 fits using exactly the same

            // common set of genuinely-detected temperature points.

            info.deltaLambda=sys[spot].deltaLambda;

            info.has_param_syst=true;

        }

        fits[spot]=info;

    }

    map<int,TRadiusOverlayFit> fits16;

    map<int,TRadiusOverlayFit> fits24;

    if (!r16_file.empty()) {

        fits16=build_T_radius_overlays(r16_file,ph,kP6Blue,24,"R16");

    }

    if (!r24_file.empty()) {

        fits24=build_T_radius_overlays(r24_file,ph,kP6Red,25,"R24");

    }

    // Individual hotspot canvases. The original T-fit aesthetics are preserved.

    string spotdir=outdir+"/lum_vs_T_all_spots_expfit";

    ensure_dir(spotdir);

    for (auto& kv:fits) {

        auto& info=kv.second;

        int spot=kv.first;

        if (!info.fit_done || !info.graph) continue;

        TCanvas* c=new TCanvas(Form("c_T_spot_%d",spot),

                               Form("Spot %d - Luminosity vs T",spot),

                               1200,800);

        info.graph->SetMarkerColor(kBlack);
        info.graph->SetLineColor(kBlack);
        info.func->SetLineColor(kBlack);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(1);

        double ymin=1e99,ymax=-1e99;

        for (const auto& r:info.rows) {

            const double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));
            ymin=std::min(ymin,r.luminosity-emax);
            ymax=std::max(ymax,r.luminosity+emax);

        }

        auto expand_overlay_range=[&](const map<int,TRadiusOverlayFit>& overlays) {

            auto it=overlays.find(spot);

            if (it==overlays.end()) return;

            for (const auto& r:it->second.rows) {

                ymin=std::min(ymin,r.luminosity-std::fabs(r.error));

                ymax=std::max(ymax,r.luminosity+std::fabs(r.error));

            }

        };

        expand_overlay_range(fits16);
        expand_overlay_range(fits24);

        if (ymin>0) ymin*=0.80; else ymin=0.0;

        ymax*=1.25;

        info.graph->SetMinimum(ymin);

        info.graph->SetMaximum(ymax);

        // Statistical uncertainty: vertical bars only (ROOT option Z removes end caps).
        info.graph->Draw("APZ");

        c->Update();

        // Systematic uncertainty: ROOT bracket representation.
        // Option [] draws only the bracket ends of the error bars.
        info.syst_graph->Draw("[] SAME");

        // Redraw markers and fit above the systematic brackets.
    

        info.func->Draw("SAME");

        TRadiusOverlayFit* ov16=nullptr;

        TRadiusOverlayFit* ov24=nullptr;

        auto it16=fits16.find(spot);

        if (it16!=fits16.end() && it16->second.graph) {

            ov16=&it16->second;

            ov16->graph->Draw("PZ SAME");

            if (ov16->fit_done && ov16->func) ov16->func->Draw("SAME");

        }

        auto it24=fits24.find(spot);

        if (it24!=fits24.end() && it24->second.graph) {

            ov24=&it24->second;

            ov24->graph->Draw("PZ SAME");

            if (ov24->fit_done && ov24->func) ov24->func->Draw("SAME");

        }

        TLegend* leg=nullptr;

        if (!comparison_mode) {

            // Original 5analysis legend.

            leg=new TLegend(0.14,0.68,0.60,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(info.graph,"Data (statistical uncertainty)","pez");
            leg->AddEntry(info.syst_graph,"Systematic uncertainty","e[]");
            leg->AddEntry(info.func,"Fit: Lum = A e^{#lambda T}","l");

        } else {

            leg=new TLegend(0.14,0.48,0.66,0.88);

            leg->SetBorderSize(0);

            leg->SetFillStyle(0);

            leg->SetTextSize(0.030);

            leg->AddEntry(info.graph,"R = 20 px data (statistical uncertainty)","pez");

            leg->AddEntry(info.syst_graph,"R = 20 px luminosity systematic uncertainty","e[]");

            string nominal_fit_label;

            if (info.has_param_syst) {

                nominal_fit_label=Form("R = 20 px fit: #lambda = %.5f #pm %.5f (stat) #pm %.5f (syst)",

                                       info.lambda,info.lambdaerr,info.deltaLambda);

            } else {

                nominal_fit_label=Form("R = 20 px fit: #lambda = %.5f #pm %.5f",

                                       info.lambda,info.lambdaerr);

            }

            leg->AddEntry(info.func,nominal_fit_label.c_str(),"l");

            if (ov16) {

                leg->AddEntry(ov16->graph,"R = 16 px data (statistical uncertainty)","pez");

                if (ov16->fit_done && ov16->func) {

                    string label16=Form("R = 16 px fit: #lambda = %.5f #pm %.5f",

                                        ov16->lambda,ov16->lambdaerr);

                    leg->AddEntry(ov16->func,label16.c_str(),"l");

                }

            }

            if (ov24) {

                leg->AddEntry(ov24->graph,"R = 24 px data (statistical uncertainty)","pez");

                if (ov24->fit_done && ov24->func) {

                    string label24=Form("R = 24 px fit: #lambda = %.5f #pm %.5f",

                                        ov24->lambda,ov24->lambdaerr);

                    leg->AddEntry(ov24->func,label24.c_str(),"l");

                }

            }

        }

        leg->Draw();

        c->SaveAs(Form("%s/%s_lum_vs_T_spot%d.png",spotdir.c_str(),pref.c_str(),spot));

        c->SaveAs(Form("%s/%s_lum_vs_T_spot%d.pdf",spotdir.c_str(),pref.c_str(),spot));

        delete c;

    }

    if (comparison_mode) {

        // Only the individual R=16/20/24 comparison canvases are generated in

        // comparison mode. Nominal lambda-vs-spot and CSV products remain in

        // the standard 5analysis.sh workflow.

        return;

    }

    // -------------------------------------------------------------------------
    // Lambda summary over hotspots.
    //
    // IMPORTANT: there is no global/simultaneous fit here.  Each lambda value
    // comes exclusively from the independent exponential fit of one hotspot.
    // The summary quantities are then computed from the selected lambda_i:
    //
    //   mean(lambda)      = (1/N) sum_i lambda_i
    //   stat. error mean  = sqrt(sum_i sigma_i^2) / N
    //   RMS               = sqrt[(1/N) sum_i (lambda_i-mean)^2]
    //
    // The same hotspot-quality selection used for the old Canvas 5 is kept.
    // -------------------------------------------------------------------------
    vector<double> x, ex, L, Lstat, Lsyst;

    for (auto& kv : fits) {
        const auto& f = kv.second;

        if (!f.fit_done ||
            !finite_number(f.lambda) ||
            !finite_number(f.lambdaerr) ||
            !finite_number(f.A) ||
            f.ndf <= 0 ||
            f.chi2ndf >= 3.5) {
            continue;
        }

        x.push_back((double)kv.first);
        ex.push_back(0.0);
        L.push_back(f.lambda);
        Lstat.push_back(f.lambdaerr);
        Lsyst.push_back(f.has_param_syst ? f.deltaLambda : 0.0);
    }

    double lambda_mean = std::numeric_limits<double>::quiet_NaN();
    double lambda_mean_stat_error = std::numeric_limits<double>::quiet_NaN();
    double lambda_rms = std::numeric_limits<double>::quiet_NaN();

    if (!L.empty()) {
        const double N = (double)L.size();

        double sum_lambda = 0.0;
        double sum_stat2 = 0.0;

        for (size_t i = 0; i < L.size(); ++i) {
            sum_lambda += L[i];
            sum_stat2 += Lstat[i] * Lstat[i];
        }

        lambda_mean = sum_lambda / N;
        lambda_mean_stat_error = std::sqrt(sum_stat2) / N;

        double sum_sq_dev = 0.0;
        for (double lambda_i : L) {
            const double d = lambda_i - lambda_mean;
            sum_sq_dev += d * d;
        }
        lambda_rms = std::sqrt(sum_sq_dev / N);
    }

    // -------------------------------------------------------------------------
    // Save one summary row for this annealing phase.
    // -------------------------------------------------------------------------
    if (!lambda_summary_file.empty()) {
        bool write_header = true;

        {
            std::ifstream check(lambda_summary_file,
                                std::ios::binary | std::ios::ate);
            if (check.is_open() && check.tellg() > 0)
                write_header = false;
        }

        std::ofstream gout(lambda_summary_file, std::ios::app);

        if (!gout.is_open()) {
            std::cerr << "Warning: cannot open lambda summary file: "
                      << lambda_summary_file << std::endl;
        } else {
            if (write_header) {
                gout << "phase,lambda_mean,lambda_mean_stat_error,lambda_rms,n_spots\n";
            }

            gout << ph << ','
                 << lambda_mean << ','
                 << lambda_mean_stat_error << ','
                 << lambda_rms << ','
                 << L.size() << '\n';
        }
    }

    // -------------------------------------------------------------------------
    // CANVAS 5: lambda_i versus global hotspot ID.
    // The dashed horizontal line is the arithmetic mean of the selected lambda_i.
    // -------------------------------------------------------------------------
    if (!L.empty()) {
        TCanvas* c5 = new TCanvas("c5_lambda_vs_spot_mean",
                                  "lambda vs spot ID",
                                  1800,1000);

        TGraphErrors* gr = new TGraphErrors((int)L.size(),
                                             x.data(), L.data(),
                                             ex.data(), Lstat.data());

        TGraphErrors* gr_syst = new TGraphErrors((int)L.size(),
                                                  x.data(), L.data(),
                                                  ex.data(), Lsyst.data());
        gr_syst->SetLineColor(kBlack);
        gr_syst->SetLineWidth(2);
        gr_syst->SetMarkerSize(0);
        gr->SetTitle(Form("%s - Exponential fit parameter #lambda - %s;Global spot ID;#lambda",
                          pref.c_str(), ph.c_str()));
        gr->SetMarkerStyle(21);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);
        // Statistical uncertainty: vertical bars only, without horizontal end caps.
        gr->Draw("APZ");

        double xmin = *std::min_element(x.begin(),x.end()) - 1.0;
        double xmax = *std::max_element(x.begin(),x.end()) + 1.0;
        gr->GetXaxis()->SetLimits(xmin,xmax);

        double lmin = *std::min_element(L.begin(),L.end());
        double lmax = *std::max_element(L.begin(),L.end());
        double span = lmax-lmin;
        if (span<=0.0) span=std::max(0.1,std::fabs(lmax)*0.2);

        // Include both the individual errors and the mean line in the visible range.
        double display_min = lmin;
        double display_max = lmax;
        for (size_t i=0;i<L.size();++i) {
            const double emax = std::max(std::fabs(Lstat[i]),std::fabs(Lsyst[i]));
            display_min = std::min(display_min,L[i]-emax);
            display_max = std::max(display_max,L[i]+emax);
        }
        display_min = std::min(display_min,lambda_mean);
        display_max = std::max(display_max,lambda_mean);

        double display_span = display_max-display_min;
        if (display_span<=0.0)
            display_span=std::max(0.1,std::fabs(display_max)*0.2);

        gr->GetYaxis()->SetRangeUser(display_min-0.20*display_span,
                                     display_max+0.45*display_span);

        c5->Update();

        // Systematic uncertainty: ROOT bracket representation.
        gr_syst->Draw("[] SAME");

        // Keep data markers clearly visible on top.
        //gr->Draw("P SAME");

        TLine* mean_line = new TLine(xmin,lambda_mean,xmax,lambda_mean);
        mean_line->SetLineColor(kBlack);
        mean_line->SetLineWidth(2);
        mean_line->SetLineStyle(2);
        mean_line->Draw("SAME");

        TLegend* leg = new TLegend(0.10,0.68,0.62,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr,"Individual #lambda (statistical uncertainty)","pez");
        leg->AddEntry(gr_syst,"Systematic uncertainty","e[]");

        string mean_label = Form("Mean #lambda = %.5f #pm %.5f (stat), RMS = %.5f",
                                 lambda_mean,
                                 lambda_mean_stat_error,
                                 lambda_rms);
        leg->AddEntry(mean_line,mean_label.c_str(),"l");
        leg->Draw();

        c5->SaveAs(Form("%s/%s_lambda_vs_spot.png",outdir.c_str(),pref.c_str()));
        c5->SaveAs(Form("%s/%s_lambda_vs_spot.pdf",outdir.c_str(),pref.c_str()));
        delete c5;
    }

    // -------------------------------------------------------------------------
    // CANVAS 6: binned distribution of the lambda_i values.
    // This plot is meant to inspect the shape of the hotspot-to-hotspot
    // distribution; no Gaussian fit is imposed here.
    // -------------------------------------------------------------------------
    if (!L.empty()) {
        double hmin = *std::min_element(L.begin(),L.end());
        double hmax = *std::max_element(L.begin(),L.end());

        if (!(hmax > hmin)) {
            const double half = std::max(0.01,0.10*std::fabs(hmax));
            hmin -= half;
            hmax += half;
        } else {
            const double margin = 0.10*(hmax-hmin);
            hmin -= margin;
            hmax += margin;
        }

        // Use a finer binning than the previous sqrt(N) choice.
        // This keeps enough granularity to inspect the distribution shape.
        int nbins = std::max(15,(int)std::ceil(2.0*std::sqrt((double)L.size())));
        nbins = std::min(nbins,60);

        TCanvas* c6 = new TCanvas("c6_lambda_distribution",
                                  "lambda distribution",
                                  1900,900);

        // Reserve a right-side area for the summary box so it does not overlap
        // the histogram.
        c6->SetRightMargin(0.28);

        TH1D* hLambda = new TH1D("h_lambda_distribution",
                                  Form("%s - Distribution of fitted #lambda - %s;#lambda;Hotspots",
                                       pref.c_str(),ph.c_str()),
                                  nbins,hmin,hmax);
        hLambda->SetLineColor(kP6Blue);
        hLambda->SetFillColor(kP6Blue);
        hLambda->SetFillStyle(1001);
        hLambda->SetLineWidth(2);

        for (double lambda_i : L)
            hLambda->Fill(lambda_i);

        hLambda->Draw("HIST");

        const double ymax_hist = hLambda->GetMaximum();
        TLine* mean_hist_line = new TLine(lambda_mean,0.0,
                                           lambda_mean,1.05*ymax_hist);
        mean_hist_line->SetLineColor(kBlack);
        mean_hist_line->SetLineWidth(2);
        mean_hist_line->SetLineStyle(2);
        mean_hist_line->Draw("SAME");

        TPaveText* stats_box = new TPaveText(0.74,0.66,0.97,0.88,"NDC");
        stats_box->SetBorderSize(0);
        stats_box->SetFillStyle(0);
        stats_box->SetTextAlign(12);
        stats_box->AddText(Form("Mean #lambda = %.5f",lambda_mean));
        stats_box->AddText(Form("Stat. error on mean = %.5f",lambda_mean_stat_error));
        stats_box->AddText(Form("RMS = %.5f",lambda_rms));
        stats_box->AddText(Form("N = %zu",L.size()));
        stats_box->Draw();

        c6->SaveAs(Form("%s/%s_lambda_distribution.png",outdir.c_str(),pref.c_str()));
        c6->SaveAs(Form("%s/%s_lambda_distribution.pdf",outdir.c_str(),pref.c_str()));
        delete c6;
    }

    string result_csv=outdir+"/"+pref+"_lambda_fit_results.csv";

    std::ofstream fout(result_csv);

    fout << "spot,x,y,phase,A,A_stat_error,lambda,lambda_stat_error,deltaLambda,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";

    for (auto& kv:fits) {

        const auto& f=kv.second;

        double xx=f.rows.empty()?0.0:f.rows[0].x;

        double yy=f.rows.empty()?0.0:f.rows[0].y;

        fout << kv.first << ',' << xx << ',' << yy << ',' << ph

             << ',' << f.A << ',' << f.Aerr

             << ',' << f.lambda << ',' << f.lambdaerr

             << ',' << (f.has_param_syst?f.deltaLambda:0.0)

             << ',' << f.chi2 << ',' << f.ndf << ',' << f.chi2ndf

             << ',' << f.prob << ',' << f.status << ',' << f.covstatus

             << ',' << f.edm << ',' << (f.converged?1:0) << '\n';

    }

    fout.close();

    write_augmented_T_csv(table,ph,fits,outdir+"/"+pref+"_analysis.csv");

}
