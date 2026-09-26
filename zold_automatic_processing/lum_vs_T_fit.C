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

struct LambdaSystematicPair {
    bool has16=false,has24=false;
    double l16=0.0,l24=0.0;
    bool conv16=false,conv24=false;
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

// -----------------------------------------------------------------------------
// Draw a symmetric systematic uncertainty as two green horizontal brackets.
// The horizontal cap is intentionally only slightly wider than the data marker,
// so statistical and systematic uncertainties remain visually distinct.
// -----------------------------------------------------------------------------
static void draw_horizontal_syst_brackets_T(const vector<double>& x,
                                             const vector<double>& y,
                                             const vector<double>& delta,
                                             double xmin,
                                             double xmax,
                                             int color=kBlack,
                                             int line_width=4) {
    if (x.empty() || y.size()!=x.size() || delta.size()!=x.size()) return;

    double span=xmax-xmin;
    if (!(span>0.0)) span=1.0;

    const double half_width=0.0045*span;

    for (size_t i=0;i<x.size();++i) {
        if (!finite_number(delta[i]) || delta[i]<=0.0) continue;

        const double ytop=y[i]+delta[i];
        const double ybottom=y[i]-delta[i];
        const double hook=0.10*delta[i];

        TLine* top=new TLine(x[i]-half_width,ytop,x[i]+half_width,ytop);
        top->SetLineColor(color); top->SetLineWidth(line_width); top->Draw("SAME");
        TLine* top_l=new TLine(x[i]-half_width,ytop,x[i]-half_width,ytop-hook);
        top_l->SetLineColor(color); top_l->SetLineWidth(line_width); top_l->Draw("SAME");
        TLine* top_r=new TLine(x[i]+half_width,ytop,x[i]+half_width,ytop-hook);
        top_r->SetLineColor(color); top_r->SetLineWidth(line_width); top_r->Draw("SAME");

        TLine* bottom=new TLine(x[i]-half_width,ybottom,x[i]+half_width,ybottom);
        bottom->SetLineColor(color); bottom->SetLineWidth(line_width); bottom->Draw("SAME");
        TLine* bottom_l=new TLine(x[i]-half_width,ybottom,x[i]-half_width,ybottom+hook);
        bottom_l->SetLineColor(color); bottom_l->SetLineWidth(line_width); bottom_l->Draw("SAME");
        TLine* bottom_r=new TLine(x[i]+half_width,ybottom,x[i]+half_width,ybottom+hook);
        bottom_r->SetLineColor(color); bottom_r->SetLineWidth(line_width); bottom_r->Draw("SAME");
    }
}

static void estimate_exp_parameters_nom(const vector<AnalysisRow>& rows,
                                        double& A0,
                                        double& lambda0) {
    vector<AnalysisRow> pos;
    for (const auto& r:rows) if (r.luminosity>0) pos.push_back(r);
    std::sort(pos.begin(),pos.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});

    if (pos.size()>=2 && std::fabs(pos.back().T-pos.front().T)>1e-12) {
        lambda0=(std::log(pos.back().luminosity)-std::log(pos.front().luminosity))/(pos.back().T-pos.front().T);
        A0=std::exp(std::log(pos.front().luminosity)-lambda0*pos.front().T);
    } else if (pos.size()==1) {
        A0=pos[0].luminosity;
        lambda0=0.0;
    } else {
        A0=1.0;
        lambda0=0.0;
    }
}


// -----------------------------------------------------------------------------
// Build optional R=16/R=24 temperature overlays. They use the same exponential
// model and initial-parameter estimate as the nominal R=20 analysis. "N" keeps
// ROOT from creating extra fit-statistics boxes for the comparison datasets.
// -----------------------------------------------------------------------------
static map<int,TRadiusOverlayFit> build_T_radius_overlays(
    const string& filename,
    const string& phase,
    int color,
    int marker,
    const string& tag) {

    map<int,TRadiusOverlayFit> out;
    if (filename.empty()) return out;

    CsvTable table=read_analysis_csv(filename,false);
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

            out[id]=s;
        } catch (...) {
            // Skip malformed rows and continue with the remaining hotspots.
        }
    }
    return out;
}

// Write a phase CSV using only the new symmetric-systematic columns.
static void write_augmented_T_csv(const CsvTable& original,
                                  const string& phase,
                                  const map<int,TNominalFit>& fits,
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
                  const char* r24_all_phases_csv = "") {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);

    string ph=phase,outdir=output_dir,pref=prefix;
    const string r16_file=(r16_all_phases_csv ? r16_all_phases_csv : "");
    const string r24_file=(r24_all_phases_csv ? r24_all_phases_csv : "");
    const bool comparison_mode=(!r16_file.empty() || !r24_file.empty());
    ensure_dir(outdir);

    CsvTable table=read_analysis_csv(all_phases_csv,false);
    auto selected=filter_phase(table.rows,ph);
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
        info.graph->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, v = %.1f, %s;T (#circC);Luminosity",
                                  spot,rows[0].x,rows[0].y,rows[0].v,ph.c_str()));
        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.4);
        info.graph->SetLineWidth(2);

        info.syst_graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),esyst.data());
        info.syst_graph->SetLineColor(kBlack);
        info.syst_graph->SetLineWidth(4);
        info.syst_graph->SetMarkerSize(0);

        double A0,l0;
        estimate_exp_parameters_nom(rows,A0,l0);
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
        if (sys.count(spot) && sys[spot].has16 && sys[spot].has24 &&
            sys[spot].conv16 && sys[spot].conv24) {
            double d1=std::fabs(sys[spot].l24-info.lambda);
            double d2=std::fabs(info.lambda-sys[spot].l16);
            info.deltaLambda=std::max(d1,d2);
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
        info.graph->Draw("AP");
        c->Update();

        vector<double> sx,sy,sd;
        for (const auto& r:info.rows) {
            sx.push_back(r.T);
            sy.push_back(r.luminosity);
            sd.push_back(r.deltaL);
        }
        draw_horizontal_syst_brackets_T(sx,sy,sd,
                                        info.graph->GetXaxis()->GetXmin(),
                                        info.graph->GetXaxis()->GetXmax());

        // Redraw the experimental markers and the fit above the systematic brackets.
        info.graph->Draw("P SAME");
        info.func->Draw("SAME");

        TRadiusOverlayFit* ov16=nullptr;
        TRadiusOverlayFit* ov24=nullptr;

        auto it16=fits16.find(spot);
        if (it16!=fits16.end() && it16->second.graph) {
            ov16=&it16->second;
            ov16->graph->Draw("PE SAME");
            if (ov16->fit_done && ov16->func) ov16->func->Draw("SAME");
        }

        auto it24=fits24.find(spot);
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
            leg=new TLegend(0.14,0.68,0.60,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(info.graph,"Data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"Systematic uncertainty","l");
            leg->AddEntry(info.func,"Fit: Lum = A e^{#lambda T}","l");
        } else {
            leg=new TLegend(0.14,0.48,0.66,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.030);

            leg->AddEntry(info.graph,"R = 20 px data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"R = 20 px luminosity systematic uncertainty","l");

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
                leg->AddEntry(ov16->graph,"R = 16 px data (statistical uncertainty)","lep");
                if (ov16->fit_done && ov16->func) {
                    string label16=Form("R = 16 px fit: #lambda = %.5f #pm %.5f",
                                        ov16->lambda,ov16->lambdaerr);
                    leg->AddEntry(ov16->func,label16.c_str(),"l");
                }
            }

            if (ov24) {
                leg->AddEntry(ov24->graph,"R = 24 px data (statistical uncertainty)","lep");
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

    // lambda vs spot: show every fit that produced a finite lambda value and
    // statistical uncertainty. Fit status, ndf and chi2/ndf are kept only as
    // diagnostics in the CSV and do not remove points from this canvas.
    vector<double> x,ex,L,Lstat,Lsyst;
    for (auto& kv:fits) {
        const auto& f=kv.second;
        if (!f.fit_done || !finite_number(f.lambda) || !finite_number(f.lambdaerr)) continue;
        x.push_back(kv.first);
        ex.push_back(0.0);
        L.push_back(f.lambda);
        Lstat.push_back(f.lambdaerr);
        Lsyst.push_back(f.has_param_syst?f.deltaLambda:0.0);
    }

    if (!L.empty()) {
        TCanvas* c5=new TCanvas("c5_lambda_vs_spot_constfit","lambda vs spot ID",1800,1000);
        TGraphErrors* gr=new TGraphErrors((int)L.size(),x.data(),L.data(),ex.data(),Lstat.data());
        gr->SetTitle(Form("Exponential fit parameter lambda vs spot ID - %s;Spot;lambda",ph.c_str()));
        gr->SetMarkerStyle(21);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);
        gr->Draw("AP");

        double xmin=*std::min_element(x.begin(),x.end())-1.0;
        double xmax=*std::max_element(x.begin(),x.end())+1.0;
        gr->GetXaxis()->SetLimits(xmin,xmax);

        double lmin=*std::min_element(L.begin(),L.end());
        double lmax=*std::max_element(L.begin(),L.end());
        double span=lmax-lmin;
        if (span<=0) span=std::max(0.1,std::fabs(lmax)*0.2);
        gr->GetYaxis()->SetRangeUser(lmin-0.25*span,lmax+0.90*span);

        c5->Update();
        draw_horizontal_syst_brackets_T(x,L,Lsyst,xmin,xmax,kBlack,4);
        gr->Draw("P SAME");

        TF1* fc=new TF1("fit_lambda_const","[0]",xmin,xmax);
        fc->SetParNames("lambda_{const}");
        fc->SetLineColor(kBlack);
        fc->SetLineWidth(2);
        if (L.size()>=2) {
            gr->Fit(fc,"RQ");
            fc->Draw("SAME");
        }

        TLine* syst_proxy_lambda=new TLine(0,0,1,0);
        syst_proxy_lambda->SetLineColor(kBlack);
        syst_proxy_lambda->SetLineWidth(4);

        TLegend* leg=new TLegend(0.10,0.70,0.54,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr,"Statistical uncertainty","lep");
        leg->AddEntry(syst_proxy_lambda,"Systematic uncertainty (R=16/24)","l");
        if (L.size()>=2) leg->AddEntry(fc,"Constant fit: lambda = const","l");
        leg->Draw();

        c5->SaveAs(Form("%s/%s_lambda_vs_spot.png",outdir.c_str(),pref.c_str()));
        c5->SaveAs(Form("%s/%s_lambda_vs_spot.pdf",outdir.c_str(),pref.c_str()));
        delete c5;
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
