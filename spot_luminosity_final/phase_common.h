#ifndef PHASE_COMMON_H
#define PHASE_COMMON_H

#include "analysis_common.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TPad.h"

struct PhaseFitPoint {
    bool ok=false;
    double val=0.0;
    double err=0.0;
};

// Select the representative point used in luminosity-vs-phase plots:
// maximum voltage for T=const data, maximum temperature for v=const data.
inline AnalysisRow choose_phase_representative(const vector<AnalysisRow>& rows,
                                               bool T_const) {
    AnalysisRow best=rows.front();
    for (const auto& r:rows) {
        if (T_const) {
            if (r.v>best.v) best=r;
        } else {
            if (r.T>best.T) best=r;
        }
    }
    return best;
}

// Fit B or lambda inside one annealing phase using the same nominal models used
// in the analysis macros. Systematic luminosity errors are displayed elsewhere
// and are not used as fit weights.
inline PhaseFitPoint fit_phase_parameter(vector<AnalysisRow> rows,
                                         bool T_const,
                                         bool has_vfin) {
    PhaseFitPoint p;

    if (T_const) {
        if (!has_vfin || rows.size()<2) return p;

        std::sort(rows.begin(),rows.end(),
                  [](const AnalysisRow&a,const AnalysisRow&b){return a.v_fin<b.v_fin;});

        vector<double> x,y,ex,ey;
        for (const auto& r:rows) {
            x.push_back(r.v_fin);
            y.push_back(r.luminosity);
            ex.push_back(0.0);
            ey.push_back(r.error);
        }

        TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
        TF1 f("phase_B_tmp","[0]*TMath::Power(x,[1])",0,8);
        f.SetParameters(1,2);
        auto res=gr.Fit(&f,"QRS0");

        if ((int)res==0 && finite_number(f.GetParameter(1))) {
            p.ok=true;
            p.val=f.GetParameter(1);
            p.err=f.GetParError(1);
        }
        return p;
    }

    if (rows.size()<3) return p;

    std::sort(rows.begin(),rows.end(),
              [](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});

    double Tmin=rows.front().T;
    double Tmax=rows.back().T;
    double d=Tmax-Tmin;
    if (d<=0) d=1.0;

    vector<double> x,y,ex,ey;
    for (const auto& r:rows) {
        x.push_back(r.T);
        y.push_back(r.luminosity);
        ex.push_back(0.0);
        ey.push_back(r.error);
    }

    double A0=1.0,l0=0.0;
    vector<AnalysisRow> pos;
    for (const auto& r:rows) if (r.luminosity>0) pos.push_back(r);
    if (pos.size()>=2 && std::fabs(pos.back().T-pos.front().T)>1e-12) {
        l0=(std::log(pos.back().luminosity)-std::log(pos.front().luminosity)) /
           (pos.back().T-pos.front().T);
        A0=std::exp(std::log(pos.front().luminosity)-l0*pos.front().T);
    }

    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
    TF1 f("phase_lambda_tmp","[0]*exp([1]*x)",Tmin-.1*d,Tmax+.1*d);
    f.SetParameters(A0,l0);
    auto res=gr.Fit(&f,"QRS0");

    if ((int)res==0 && finite_number(f.GetParameter(1))) {
        p.ok=true;
        p.val=f.GetParameter(1);
        p.err=f.GetParError(1);
    }
    return p;
}

inline void run_phase_analysis(const char* csvfile,
                               const char* output_dir,
                               const char* prefix,
                               bool T_const) {
    gStyle->SetOptStat(0);
    const int MAX_SPOTS_PER_CANVAS=15;

    string out=output_dir,pref=prefix;
    ensure_dir(out);
    ensure_dir(out+"/single_spots");
    ensure_dir(out+"/grouped");
    ensure_dir(out+"/ratios");
    ensure_dir(out+"/fit_parameter_vs_phase");

    CsvTable table=read_analysis_csv(csvfile,T_const);
    if (table.rows.empty()) return;
    bool has_vfin=has_col(table,"v_fin");

    vector<string> phases=standard_phases();
    map<string,int> pidx;
    for (int i=0;i<(int)phases.size();++i) pidx[phases[i]]=i;

    // Store all operating points for per-phase fits.
    map<int,map<string,vector<AnalysisRow>>> all_by_spot;
    for (const auto& r:table.rows) {
        if (pidx.count(r.phase)) all_by_spot[r.spot][r.phase].push_back(r);
    }

    // Store one representative point per phase for luminosity-vs-phase plots.
    map<int,map<string,AnalysisRow>> rep;
    for (auto& skv:all_by_spot) {
        for (auto& pkv:skv.second) {
            if (!pkv.second.empty())
                rep[skv.first][pkv.first]=choose_phase_representative(pkv.second,T_const);
        }
    }

    vector<int> colors={
        kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+1,
        kOrange+7,kCyan+2,kViolet+1,kAzure+1,kPink+7,
        kSpring+5,kTeal+3,kGray+2,kRed-4,kBlue-4
    };
    vector<int> markers={20,21,22,23,24,25,26,27,28,29,30,33,34,43,47};

    // =====================================================================
    // SINGLE-HOTSPOT LUMINOSITY VS PHASE
    // =====================================================================
    for (auto& skv:rep) {
        int spot=skv.first;
        vector<double> x,y,ex,estat,esyst;
        double ymin=1e99,ymax=-1e99;

        for (int i=0;i<(int)phases.size();++i) {
            if (!skv.second.count(phases[i])) continue;
            const auto& r=skv.second[phases[i]];

            x.push_back(i);
            y.push_back(r.luminosity);
            ex.push_back(0.0);
            estat.push_back(r.error);
            esyst.push_back(r.deltaL);

            double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));
            ymin=std::min(ymin,r.luminosity-emax);
            ymax=std::max(ymax,r.luminosity+emax);
        }

        if (x.empty()) continue;
        if (ymin>0) ymin*=.95;
        ymax*=1.05;
        if (ymax<=ymin) ymax=ymin+1.0;

        TCanvas* c=new TCanvas(Form("c_phase_single_%d",spot),
                               Form("Spot %d - Luminosity vs phase",spot),
                               1200,800);
        c->SetGrid();

        TH1D* frame=new TH1D(Form("frame_single_%d",spot),
                             Form("Spot %d: luminosity vs phase",spot),
                             (int)phases.size(),-.5,(double)phases.size()-.5);
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle("phase");
        frame->GetYaxis()->SetTitle("luminosity");
        for (int i=0;i<(int)phases.size();++i)
            frame->GetXaxis()->SetBinLabel(i+1,phases[i].c_str());
        frame->GetXaxis()->LabelsOption("h");
        frame->GetXaxis()->SetLabelSize(.035);
        frame->GetYaxis()->SetTitleOffset(1.25);
        frame->Draw();

        TGraphErrors* gs=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),estat.data());
        gs->SetLineColor(kBlue+1);
        gs->SetMarkerColor(kBlue+1);
        gs->SetMarkerStyle(20);
        gs->SetMarkerSize(1.2);
        gs->SetLineWidth(2);
        gs->Draw("PL SAME");

        TGraphErrors* gy=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),esyst.data());
        gy->SetLineColor(kGray+2);
        gy->SetLineWidth(2);
        gy->SetMarkerSize(0);
        gy->Draw("[] SAME");
        gs->Draw("PL SAME");

        TLegend* leg=new TLegend(.62,.68,.90,.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gs,Form("spot %d - statistical",spot),"pl");
        leg->AddEntry(gy,"systematic","l");
        leg->Draw();

        c->SaveAs(Form("%s/single_spots/%s_luminosity_vs_phase_spot%d.png",
                       out.c_str(),pref.c_str(),spot));
        delete c;
    }

    // =====================================================================
    // GROUPED LUMINOSITY VS PHASE
    // =====================================================================
    map<int,vector<int>> groups;
    for (auto& skv:rep) {
        double m=-1.0;
        for (auto& pkv:skv.second) m=std::max(m,pkv.second.luminosity);
        groups[m>0 ? (int)std::floor(std::log10(m)) : -999].push_back(skv.first);
    }

    int counter=0;
    for (auto& gkv:groups) {
        auto ids=gkv.second;
        std::sort(ids.begin(),ids.end());

        for (int start=0;start<(int)ids.size();start+=MAX_SPOTS_PER_CANVAS) {
            int end=std::min(start+MAX_SPOTS_PER_CANVAS,(int)ids.size());
            vector<int> sub(ids.begin()+start,ids.begin()+end);

            double ymin=1e99,ymax=-1e99;
            for (int id:sub) {
                for (auto& pkv:rep[id]) {
                    const auto& r=pkv.second;
                    double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));
                    ymin=std::min(ymin,r.luminosity-emax);
                    ymax=std::max(ymax,r.luminosity+emax);
                }
            }

            if (ymin>0) ymin*=.8; else ymin*=1.2;
            ymax*=1.25;
            if (ymax<=ymin) ymax=ymin+1.0;

            string title=(gkv.first==-999)
                ? "Luminosity vs phase - zero or negative luminosity"
                : Form("Luminosity vs phase - order 10^{%d}",gkv.first);

            TCanvas* c=new TCanvas(Form("c_phase_group_%d",counter),title.c_str(),1100,750);
            c->SetGrid();

            TH1D* frame=new TH1D(Form("frame_group_%d",counter),title.c_str(),
                                 (int)phases.size(),-.5,(double)phases.size()-.5);
            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);
            frame->GetXaxis()->SetTitle("phase");
            frame->GetYaxis()->SetTitle("luminosity");
            for (int i=0;i<(int)phases.size();++i)
                frame->GetXaxis()->SetBinLabel(i+1,phases[i].c_str());
            frame->GetXaxis()->LabelsOption("v");
            frame->GetXaxis()->SetLabelSize(.035);
            frame->GetYaxis()->SetTitleOffset(1.25);
            frame->Draw();

            TLegend* leg=new TLegend(.72,.58,.90,.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            int ig=0;
            for (int id:sub) {
                vector<double> x,y,ex,estat,esyst;
                for (int i=0;i<(int)phases.size();++i) {
                    if (!rep[id].count(phases[i])) continue;
                    const auto& r=rep[id][phases[i]];
                    x.push_back(i);
                    y.push_back(r.luminosity);
                    ex.push_back(0.0);
                    estat.push_back(r.error);
                    esyst.push_back(r.deltaL);
                }
                if (x.empty()) continue;

                int col=colors[ig%colors.size()];
                int mar=markers[ig%markers.size()];

                TGraphErrors* gs=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),estat.data());
                gs->SetLineColor(col);
                gs->SetMarkerColor(col);
                gs->SetMarkerStyle(mar);
                gs->SetMarkerSize(1.1);
                gs->SetLineWidth(2);
                gs->Draw("PL SAME");

                TGraphErrors* gy=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),esyst.data());
                gy->SetLineColor(col);
                gy->SetLineStyle(2);
                gy->SetLineWidth(1);
                gy->SetMarkerSize(0);
                gy->Draw("[] SAME");

                leg->AddEntry(gs,Form("spot %d",id),"pl");
                ++ig;
            }

            leg->Draw();
            c->SaveAs(Form("%s/grouped/%s_luminosity_vs_phase_group_%02d_order_%d.png",
                           out.c_str(),pref.c_str(),counter,gkv.first));
            delete c;
            ++counter;
        }
    }

    // =====================================================================
    // RATIO FOR EACH ANNEALING PHASE RELATIVE TO BEFORE ANNEALING
    // Both statistical and systematic ratio errors are propagated separately.
    // The systematic ratio error is symmetric because deltaL is symmetric.
    // =====================================================================
    const string before="before_annealing";

    for (const auto& ph:phases) {
        if (ph==before) continue;

        vector<double> x,y,ex,estat,esyst;

        for (auto& skv:rep) {
            int id=skv.first;
            if (!skv.second.count(before) || !skv.second.count(ph)) continue;

            const auto& a=skv.second[before];
            const auto& b=skv.second[ph];
            if (a.luminosity==0.0) continue;

            double R=b.luminosity/a.luminosity;
            double stat=std::sqrt(
                std::pow(b.error/a.luminosity,2) +
                std::pow(b.luminosity*a.error/(a.luminosity*a.luminosity),2)
            );
            double syst=std::sqrt(
                std::pow(b.deltaL/a.luminosity,2) +
                std::pow(b.luminosity*a.deltaL/(a.luminosity*a.luminosity),2)
            );

            x.push_back(id);
            y.push_back(R);
            ex.push_back(0.0);
            estat.push_back(stat);
            esyst.push_back(syst);
        }

        if (x.empty()) continue;

        double ymin=1e99,ymax=-1e99;
        for (size_t i=0;i<y.size();++i) {
            double emax=std::max(estat[i],esyst[i]);
            ymin=std::min(ymin,y[i]-emax);
            ymax=std::max(ymax,y[i]+emax);
        }
        if (ymin>0) ymin*=.8; else ymin*=1.2;
        ymax*=1.25;

        TCanvas* c=new TCanvas(Form("c_ratio_%s",safe_token(ph).c_str()),
                               Form("Ratio %s / before_annealing",ph.c_str()),
                               1100,750);
        c->SetGrid();

        TGraphErrors* gs=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),estat.data());
        gs->SetTitle(Form("L(%s) / L(before_annealing);spot;ratio",ph.c_str()));
        gs->SetMarkerStyle(20);
        gs->SetMarkerSize(1.2);
        gs->SetLineWidth(2);
        gs->Draw("AP");
        gs->GetYaxis()->SetRangeUser(ymin,ymax);

        TGraphErrors* gy=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),esyst.data());
        gy->SetLineColor(kGray+2);
        gy->SetLineWidth(2);
        gy->SetMarkerSize(0);
        gy->Draw("[] SAME");
        gs->Draw("P SAME");

        c->Update();
        double xmin=gs->GetXaxis()->GetXmin();
        double xmax=gs->GetXaxis()->GetXmax();
        TLine* one=new TLine(xmin,1.0,xmax,1.0);
        one->SetLineStyle(2);
        one->SetLineWidth(2);
        one->Draw("SAME");

        TLegend* leg=new TLegend(.64,.74,.90,.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gs,"Statistical uncertainty","lep");
        leg->AddEntry(gy,"Systematic uncertainty","l");
        leg->Draw();

        c->SaveAs(Form("%s/ratios/%s_ratio_%s_over_before.png",
                       out.c_str(),pref.c_str(),safe_token(ph).c_str()));
        delete c;
    }

    // =====================================================================
    // B OR LAMBDA EVOLUTION WITH PHASE: 8 HOTSPOTS PER 4x2 CANVAS
    // =====================================================================
    if (T_const && !has_vfin) {
        std::cerr << "Warning: v_fin is absent. B-vs-phase canvases cannot be produced reliably and will be skipped."
                  << std::endl;
    }

    vector<int> ids;
    for (auto& skv:all_by_spot) ids.push_back(skv.first);
    std::sort(ids.begin(),ids.end());

    int canvas_id=0;
    for (int start=0;start<(int)ids.size();start+=8) {
        TCanvas* c=new TCanvas(Form("c_parameter_phase_%d",canvas_id),
                               T_const?"B vs phase":"lambda vs phase",
                               1800,900);
        c->Divide(4,2);
        bool any=false;

        for (int j=0;j<8 && start+j<(int)ids.size();++j) {
            int id=ids[start+j];
            c->cd(j+1);

            vector<double> x,y,ex,ey;
            for (int i=0;i<(int)phases.size();++i) {
                auto it=all_by_spot[id].find(phases[i]);
                if (it==all_by_spot[id].end()) continue;

                PhaseFitPoint p=fit_phase_parameter(it->second,T_const,has_vfin);
                if (!p.ok) continue;

                x.push_back(i);
                y.push_back(p.val);
                ex.push_back(0.0);
                ey.push_back(p.err);
            }

            if (x.empty()) continue;
            any=true;

            TGraphErrors* gr=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
            gr->SetTitle(Form("Spot %d;%s;%s",id,"phase",T_const?"B":"#lambda"));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(1.0);
            gr->SetLineWidth(2);
            gr->Draw("APL");
            gr->GetXaxis()->SetLimits(-.5,(double)phases.size()-.5);
            gr->GetXaxis()->SetNdivisions((int)phases.size(),false);
        }

        if (any) {
            c->SaveAs(Form("%s/fit_parameter_vs_phase/%s_%s_vs_phase_canvas_%02d.png",
                           out.c_str(),pref.c_str(),T_const?"B":"lambda",canvas_id));
        }
        delete c;
        ++canvas_id;
    }
}

#endif
