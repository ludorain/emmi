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

struct BSystematicPair {
    bool has16=false, has24=false;
    double B16=0.0, B24=0.0;
    bool conv16=false, conv24=false;
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
// Read the R=16 and R=24 fit values produced by lum_vs_v_systematics.C.
// ROOT fit status == 0 is the convergence criterion; ndf is intentionally not
// required because a two-point/two-parameter fit can converge with ndf == 0.
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

            out[id]=s;
        } catch (...) {
            // Skip malformed systematic rows without aborting the full analysis.
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
// Draw a symmetric systematic uncertainty as two green horizontal brackets.
// The bracket width is defined as a small fraction of the visible x range, so it
// remains only slightly wider than the experimental marker for any axis scale.
// This is deliberately different from a conventional statistical error bar.
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
// Write an analysis CSV using the new symmetric-systematic convention.
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
    auto selected=filter_phase(table.rows,phase);
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
        info.func->SetParameters(A_init,B_init);
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

// -----------------------------------------------------------------------------
// Nominal R=20 power-law fit: Lum = A * V_over^B.
// Luminosity systematic uncertainties are displayed but are NOT fit weights.
// -----------------------------------------------------------------------------
void lum_vs_v_fit(const char* all_phases_csv,
                  const char* phase,
                  const char* systematic_values_csv,
                  const char* output_dir,
                  const char* prefix,
                  const char* r16_all_phases_csv = "",
                  const char* r24_all_phases_csv = "") {
    gStyle->SetOptFit(000);

    const double FIT_XMIN=0.0;
    const double FIT_XMAX=8.0;
    const double A_INIT=1.0;
    const double B_INIT=2.0;

    string ph=phase, outdir=output_dir, pref=prefix;
    const string r16_file=(r16_all_phases_csv ? r16_all_phases_csv : "");
    const string r24_file=(r24_all_phases_csv ? r24_all_phases_csv : "");
    const bool comparison_mode=(!r16_file.empty() || !r24_file.empty());
    ensure_dir(outdir);

    CsvTable table=read_analysis_csv(all_phases_csv,true);
    auto selected=filter_phase(table.rows,ph);
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
        info.graph->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f #circC %s;Overvoltage (V);Luminosity",
                                  spot,rows[0].x,rows[0].y,rows[0].T,ph.c_str()));
        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetLineWidth(2);

        info.func=new TF1(Form("f_quad_spot_%d",spot),
                          "[0]*TMath::Power(x,[1])",
                          FIT_XMIN,FIT_XMAX);
        info.func->SetParameters(A_INIT,B_INIT);
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
        if (sys.count(spot) && sys[spot].has16 && sys[spot].has24 &&
            sys[spot].conv16 && sys[spot].conv24) {
            const double deltaB1=std::fabs(sys[spot].B24-info.B);
            const double deltaB2=std::fabs(info.B-sys[spot].B16);
            info.deltaB=std::max(deltaB1,deltaB2);
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

        info.graph->Draw("AP");
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
            ov16->graph->Draw("PE SAME");
            if (ov16->fit_done && ov16->func) ov16->func->Draw("SAME");
        }

        auto it24=fits24.find(id);
        if (it24!=fits24.end() && it24->second.graph) {
            ov24=&it24->second;
            ov24->graph->Draw("PE SAME");
            if (ov24->fit_done && ov24->func) ov24->func->Draw("SAME");
        }

        TLine* syst_proxy=new TLine(0,0,1,0);
        syst_proxy->SetLineColor(kBlack);
        syst_proxy->SetLineWidth(4);

        TLegend* leg=nullptr;
        if (!comparison_mode) {
            // Original 5analysis legend.
            leg=new TLegend(0.14,0.63,0.43,0.80);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(info.graph,"Data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"Systematic uncertainty","l");
            leg->AddEntry(info.func,"Fit: Lum = A#upointV_{over}^{B}","l");
        } else {
            // Radius-comparison legend: R=20 remains the nominal dataset.
            leg=new TLegend(0.14,0.48,0.64,0.80);
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
        
        TPaveStats* stats=(TPaveStats*)info.graph->FindObject("stats");
        if (stats) {
            // Compact nominal R=20 fit-statistics box, aligned at the top with
            // the legend. Comparison-fit values are reported in the legend.
            if (comparison_mode) {
                stats->SetX1NDC(0.50);
                stats->SetX2NDC(0.80);
            } else {
                stats->SetX1NDC(0.36);
                stats->SetX2NDC(0.60);
            }
            stats->SetY1NDC(0.79);
            stats->SetY2NDC(0.99);
            stats->SetTextSize(0.020);
            stats->Draw();
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
    // IMPORTANT: no chi2/ndf cut. Every converged ROOT fit is shown.
    // =====================================================================
    vector<double> x, ex, B, Bstat, Bsyst;
    for (auto& kv:fits) {
        const auto& f=kv.second;
        if (!f.converged) continue;

        x.push_back((double)kv.first);
        ex.push_back(0.0);
        B.push_back(f.B);
        Bstat.push_back(f.Berr);
        Bsyst.push_back(f.has_param_syst ? f.deltaB : 0.0);
    }

    if (!B.empty()) {
        TCanvas* c9=new TCanvas("c9","B vs spot",1800,1000);
        TGraphErrors* grB=new TGraphErrors((int)B.size(),x.data(),B.data(),ex.data(),Bstat.data());
        grB->SetTitle(Form("Power law exponent: B - %s;Spot;B",ph.c_str()));
        grB->SetMarkerStyle(21);
        grB->SetMarkerSize(1.2);
        grB->SetLineWidth(2);

        double xmin=*std::min_element(x.begin(),x.end())-1.0;
        double xmax=*std::max_element(x.begin(),x.end())+1.0;

        double ymin=1e99, ymax=-1e99;
        for (size_t i=0;i<B.size();++i) {
            const double emax=std::max(std::fabs(Bstat[i]),std::fabs(Bsyst[i]));
            ymin=std::min(ymin,B[i]-emax);
            ymax=std::max(ymax,B[i]+emax);
        }
        double yspan=ymax-ymin;
        if (!(yspan>0.0)) yspan=std::max(0.2,std::fabs(ymax)*0.2);

        grB->SetMinimum(ymin-0.15*yspan);
        grB->SetMaximum(ymax+0.35*yspan);
        grB->Draw("AP");
        grB->GetXaxis()->SetLimits(xmin,xmax);
        c9->Update();

        draw_horizontal_syst_brackets(x,B,Bsyst,xmin,xmax,kGreen+2,4);
        grB->Draw("P SAME");

        TF1* line_B=new TF1("line_B","2",xmin,xmax);
        line_B->SetLineColor(kRed+1);
        line_B->SetLineStyle(2);
        line_B->SetLineWidth(2);
        line_B->Draw("SAME");

        TF1* fit_B=new TF1("fit_B","[0]",xmin,xmax);
        fit_B->SetLineColor(kP6Red);
        fit_B->SetLineWidth(2);
        if (B.size()>=2) {
            grB->Fit(fit_B,"RQ");
            fit_B->Draw("SAME");
        }

        TLine* syst_proxy=new TLine(0,0,1,0);
        syst_proxy->SetLineColor(kGreen+2);
        syst_proxy->SetLineWidth(4);

        TLegend* leg=new TLegend(0.11,0.68,0.43,0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(grB,"B statistical uncertainty","lep");
        leg->AddEntry(syst_proxy,"B systematic uncertainty","l");
        if (B.size()>=2) leg->AddEntry(fit_B,"Constant fit: B = const","l");
        leg->AddEntry(line_B,"Line: B = 2","l");
        leg->Draw();

        c9->SaveAs(Form("%s/%s_B_vs_spot.png",outdir.c_str(),pref.c_str()));
        delete c9;
    } else {
        std::cerr << "Warning: no ROOT-converged B fit is available for Canvas 9 in phase "
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
