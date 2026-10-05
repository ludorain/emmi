#ifndef PHASE_COMMON_H
#define PHASE_COMMON_H

#include "analysis_common.h"
#include "TCanvas.h"
#include "TGraph.h"
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

inline AnalysisRow choose_phase_representative(const vector<AnalysisRow>& rows,bool T_const) {
    AnalysisRow best=rows.front();
    for (const auto& r:rows) {
        if (T_const) { if (r.v>best.v) best=r; }
        else { if (r.T>best.T) best=r; }
    }
    return best;
}

inline PhaseFitPoint fit_phase_parameter(vector<AnalysisRow> rows,bool T_const,bool has_vfin) {
    PhaseFitPoint p;
    rows=filter_detected(rows);
    if (T_const) {
        if (!has_vfin || rows.size()<2) return p;
        std::sort(rows.begin(),rows.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.v_fin<b.v_fin;});
        vector<double>x,y,ex,ey;
        for(const auto&r:rows){x.push_back(r.v_fin);y.push_back(r.luminosity);ex.push_back(0.0);ey.push_back(r.error);}
        TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
        TF1 f(Form("phase_B_tmp_%p",(void*)&gr),"[0]*TMath::Power(x,[1])",0,8);
        f.SetParameters(1,2);
        auto res=gr.Fit(&f,"QRS0");
        if((int)res==0&&finite_number(f.GetParameter(1))){p.ok=true;p.val=f.GetParameter(1);p.err=f.GetParError(1);}
        return p;
    }
    if(rows.size()<3)return p;
    std::sort(rows.begin(),rows.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
    double Tmin=rows.front().T,Tmax=rows.back().T,d=Tmax-Tmin;if(d<=0)d=1.0;
    vector<double>x,y,ex,ey;for(const auto&r:rows){x.push_back(r.T);y.push_back(r.luminosity);ex.push_back(0.0);ey.push_back(r.error);}
    double A0=1.0,l0=0.0;vector<AnalysisRow>pos;for(const auto&r:rows)if(r.luminosity>0)pos.push_back(r);
    if(pos.size()>=2&&std::fabs(pos.back().T-pos.front().T)>1e-12){l0=(std::log(pos.back().luminosity)-std::log(pos.front().luminosity))/(pos.back().T-pos.front().T);A0=std::exp(std::log(pos.front().luminosity)-l0*pos.front().T);}
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
    TF1 f(Form("phase_lambda_tmp_%p",(void*)&gr),"[0]*exp([1]*x)",Tmin-.1*d,Tmax+.1*d);
    f.SetParameters(A0,l0);auto res=gr.Fit(&f,"QRS0");
    if((int)res==0&&finite_number(f.GetParameter(1))){p.ok=true;p.val=f.GetParameter(1);p.err=f.GetParError(1);}
    return p;
}

inline vector<string> phases_present_in_rows(const vector<AnalysisRow>& rows) {
    std::set<string> present; for(const auto&r:rows)present.insert(r.phase);
    vector<string> out;
    for(const auto&p:standard_phases()) if(present.count(p)) out.push_back(p);
    for(const auto&p:present) if(std::find(out.begin(),out.end(),p)==out.end()) out.push_back(p);
    return out;
}

inline const AnalysisRow* find_operating_point(const vector<AnalysisRow>& rows,const AnalysisRow& ref,bool T_const,double tol=1e-9){
    const double xref=T_const?ref.v_fin:ref.T;
    for(const auto&r:rows){double x=T_const?r.v_fin:r.T;if(std::fabs(x-xref)<=tol)return &r;}
    return nullptr;
}

inline void common_radius_rows(const vector<AnalysisRow>& r20,const vector<AnalysisRow>& r16,const vector<AnalysisRow>& r24,
                               bool T_const,vector<AnalysisRow>& o20,vector<AnalysisRow>& o16,vector<AnalysisRow>& o24){
    auto a20=filter_detected(r20),a16=filter_detected(r16),a24=filter_detected(r24);
    for(const auto&r:a20){const AnalysisRow*p16=find_operating_point(a16,r,T_const);const AnalysisRow*p24=find_operating_point(a24,r,T_const);if(p16&&p24){o20.push_back(r);o16.push_back(*p16);o24.push_back(*p24);}}
}

inline void save_canvas_both(TCanvas* c,const string& base){c->SaveAs((base+".png").c_str());c->SaveAs((base+".pdf").c_str());}

struct PhaseInfo {
    string raw;
    bool before = false;
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double hours = std::numeric_limits<double>::quiet_NaN();
};

static PhaseInfo parse_phase_for_spacing(const string& phase)
{
    PhaseInfo p;
    p.raw = phase;

    if (phase == "before_annealing") {
        p.before = true;
        return p;
    }

    const string prefix = "annealing_T=";
    if (phase.find(prefix) != 0) return p;

    const size_t pos_h = phase.find("_h=");
    if (pos_h == string::npos) return p;

    try {
        const string Tstr = phase.substr(prefix.size(), pos_h - prefix.size());
        const string hstr = phase.substr(pos_h + 3);
        p.temperature = stod(Tstr);
        p.hours = stod(hstr);
    }
    catch (...) {
        p.temperature = std::numeric_limits<double>::quiet_NaN();
        p.hours = std::numeric_limits<double>::quiet_NaN();
    }

    return p;
}

static double phase_step_units(const string& phase)
{
    const PhaseInfo p = parse_phase_for_spacing(phase);

    if (p.before) return 0.0;

    if (finite_number(p.hours) && p.hours > 0.0)
        return p.hours / 5.0;

    return 1.0;
}

static vector<double> build_phase_positions(const vector<string>& phases)
{
    vector<double> x(phases.size(), 0.0);
    if (phases.empty()) return x;

    x[0] = 0.0;

    for (size_t i = 1; i < phases.size(); ++i) {
        double step = phase_step_units(phases[i]);
        if (!(step > 0.0) || !finite_number(step)) step = 1.0;
        x[i] = x[i - 1] + step;
    }

    return x;
}

// Etichette asse x su due righe
static string two_line_phase_label(const string& phase)
{
    const PhaseInfo p = parse_phase_for_spacing(phase);

    if (p.before)
        return "#splitline{Before}{ann.}";

    if (finite_number(p.temperature) && finite_number(p.hours)) {

        string Ttext;
        string htext;

        if (std::fabs(p.temperature - std::round(p.temperature)) < 1e-9)
            Ttext = Form("%.0f#circC", p.temperature);
        else
            Ttext = Form("%.1f#circC", p.temperature);

        if (std::fabs(p.hours - std::round(p.hours)) < 1e-9)
            htext = Form("%.0f h", p.hours);
        else
            htext = Form("%.1f h", p.hours);

        return "#splitline{" + Ttext + "}{" + htext + "}";
    }

    return phase;
}


//Riadattamento asse x per step proporzionali all'annealing
static TH1D* make_phase_frame(const char* name,
                              const char* title,
                              const vector<string>& phases,
                              const vector<double>& xpos,
                              double ymin,
                              double ymax)
{
    if (phases.empty() || xpos.empty()) return nullptr;

    const int max_unit = std::max(0, (int)std::llround(xpos.back()));
    const int nbins = max_unit + 1;

    TH1D* frame = new TH1D(name, title, nbins, -0.5, max_unit + 0.5);

    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);

    frame->GetXaxis()->SetLabelSize(0.0);   // nasconde le label automatiche
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.95);

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);

    return frame;
}

static void draw_horizontal_phase_labels(TPad* pad,
                                         const vector<string>& phases,
                                         const vector<double>& xpos,
                                         double xmin,
                                         double xmax)
{
    if (!pad || phases.empty() || xpos.empty() || !(xmax > xmin)) return;

    pad->Update();

    const double left   = pad->GetLeftMargin();
    const double right  = pad->GetRightMargin();
    const double usable = 1.0 - left - right;

    for (size_t i = 0; i < phases.size(); ++i) {

        const double frac  = (xpos[i] - xmin) / (xmax - xmin);
        const double x_ndc = left + frac * usable;

        TLatex* lab = new TLatex();
        lab->SetNDC();
        lab->SetTextFont(42);
        lab->SetTextSize(0.028);   // leggermente più piccolo per evitare overlap
        lab->SetTextAlign(22);
        if (i == 0 || i== 2 || i == 4 || i == 6 || i == 8) {
            lab->DrawLatex(x_ndc - 0.003, 0.15, two_line_phase_label(phases[i]).c_str());
        } else {
            lab->DrawLatex(x_ndc + 0.002, 0.15, two_line_phase_label(phases[i]).c_str());
        }
        //lab->DrawLatex(x_ndc, 0.25, two_line_phase_label(phases[i]).c_str());
    }
}



inline void run_phase_analysis(const char* csvfile,const char* output_dir,const char* prefix,bool T_const,
                               const char* csv_R16="",const char* csv_R24="") {
    gStyle->SetOptStat(0);
    const int MAX_SPOTS_PER_CANVAS=6;
    const int DETECTED_COLOR=kP6Blue;
    const int FORCED_COLOR=kP6Red;
    const vector<int> line_colors={kP10Cyan, kP10Ash, kP10Green, kP10Orange, kP10Brown, kP10Red, kP10Yellow,kP10Violet,kP10Blue, kP8Pink, kBlack, kP6Grape, kP10Gray};

    // =====================================================================
    // Manually selected hotspots
    // =====================================================================
    struct ManualPhaseGroup {
        string title;
        vector<int> ids;
    };


    string out=output_dir,pref=prefix;

    vector<ManualPhaseGroup> MANUAL_GROUPS;
    if (pref == "A1_T=20") {

        MANUAL_GROUPS = {
            {
                "Decreasing and disappearing uniformly, L>1000",
                {13, 64}
            },

            {
                "Decreasing and disappearing uniformly, 400< L < 1000",
                {72, 4, 8, 10}
            },

            {
                "Decreasing and disappearing uniformly, L < 400",
                {15, 16, 26, 39}
            },
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                "Decreasing and disappearing fluctuating, L>1000",
                {87}
            },

            {
                "Decreasing and disappearing fluctuating, 400< L < 1000",
                {58, 35}
            },

            {   //Some examples
                "Decreasing and disappearing fluctuating, L < 400",
                {2, 9, 11, 14, 19, 23, 28, 44, 52, 53, 71, 74, 76}
            },
            // ==============================================================
            {
                "Decreasing uniformly, L>800",
                {18, 48, 54, 67, 85, 91}
            },

            {
                "Decreasing uniformly, L < 800",
                {22, 56, 86, 93, 95}
            },

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                "Decreasing fluctuating,  L>800",
                {0, 3, 6, 33, 57, 59, 69, 77, 88, 92}
            },

            {   //Some examples
                "Decreasing fluctuating,  L<800",
                {1, 5, 12, 17, 31, 32, 34, 70}
            },

            // ==============================================================
            {
                "Peak for high L",
                { 37}
            },

            {
                "Peaks for low L",
                {7, 50}
            },

            // ==============================================================
            {
                "Approximately constant hotspots",
                {94, 96}
            },

             // ==============================================================
            {
                "Increasing hotspots",
                {25, 42, 80, 97}
            },

            // ==============================================================
            {
                "Appeared",
                {78, 83}
            }

        };

    }
    else if (pref == "B1_T=20") {
MANUAL_GROUPS = {
            {
                "Decreasing and disappearing uniformly, L>1000",
                {23, 30, 42}
            },

            {
                "Decreasing and disappearing uniformly, 400< L < 1000",
                {0, 17, 18, 38, 61}
            },

            {   //Some examples
                "Decreasing and disappearing uniformly, L < 400",
                {1, 3, 4, 8, 12, 13, 15, 22, 24, 32}
            },
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                "Decreasing and disappearing fluctuating, L>1000",
                {20, 25, 59}
            },

            {
                "Decreasing and disappearing fluctuating, 400< L < 1000",
                {5, 37, 40, 43, 45, 52, 102}
            },

            {   //Some examples
                "Decreasing and disappearing fluctuating, L < 400",
                {10, 19, 27, 28, 35, 48, 58, 62, 69, 72, 75, 81, 82, 96}
            },
            // ==============================================================
            {
                "Decreasing uniformly, L>800",
                {34, 36, 65, 78, 85, 100, 101, 107, 112}
            },

            {
                "Decreasing uniformly, L < 800",
                {6, 74, 104, 105}
            },

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            {
                "Decreasing fluctuating,  L>800",
                {11, 26, 29, 53, 60, 89, 90, 114, 115, 117}
            },

            {   //Some examples
                "Decreasing fluctuating,  L<800",
                {31, 46, 64, 83, 88, 91, 92, 93, 94, 108, 110}
            },

            // ==============================================================
            {
                "Peaks for high L",
                {9, 14, 16, 70, 87}
            },

            {
                "Peaks for low L",
                {106}
            },

            // ==============================================================
            {
                "Approximately constant hotspots",
                {21, 54, 79, 103}
            },

             // ==============================================================
            {
                "Increasing hotspots",
                {109}
            },

            // ==============================================================
            {
                "Appeared high L - (1)",
                {71, 76}
            }, 

            {
                "Appeared high L - (2)",
                {86}
            }, 

            {
                "Appeared low L",
                {77}
            }



        };
        
    }
    ensure_dir(out);
    ensure_dir(out+"/single_spots");
    ensure_dir(out+"/grouped");
    ensure_dir(out+"/manual_grouped");
    ensure_dir(out+"/ratios");
    ensure_dir(out+"/fit_parameter_vs_phase");

    CsvTable table=read_analysis_csv(csvfile,T_const);if(table.rows.empty())return;
    bool has_vfin=has_col(table,"v_fin");
    vector<string> phases=phases_present_in_rows(table.rows);if(phases.empty())return;
    map<string,int>pidx;for(int i=0;i<(int)phases.size();++i)pidx[phases[i]]=i;


    //Pre-calculate the position of phases on x axis
    vector<double> phase_x = build_phase_positions(phases);
    double phase_xmin = -1.0;
    double phase_xmax = phase_x.empty() ? 0.5 : phase_x.back() + 0.5;

    string condition=analysis_condition_label(pref,T_const);
    double representative=-1e99;
    for(const auto&r:table.rows) representative=std::max(representative,T_const?r.v:r.T);
    string sensor = sensor_from_prefix(pref);

    double v_over_rep = representative;

    if (sensor == "A1")
        v_over_rep = 7.0;
    else if (sensor == "B1")
        v_over_rep = 5.0;

    string full_condition = T_const
        ? Form("%s, v_{over} = %.0f V",
            condition.c_str(),
            v_over_rep)
        : Form("%s, T = %.3g #circC",
            condition.c_str(),
            representative);

    map<int,map<string,vector<AnalysisRow>>> all_by_spot;
    for(const auto&r:table.rows)if(pidx.count(r.phase))all_by_spot[r.spot][r.phase].push_back(r);
    map<int,map<string,AnalysisRow>> rep;
    for(auto&skv:all_by_spot)for(auto&pkv:skv.second)if(!pkv.second.empty())rep[skv.first][pkv.first]=choose_phase_representative(pkv.second,T_const);

    // Optional radius tables used only for systematic errors on B/lambda vs phase.
    map<int,map<string,vector<AnalysisRow>>> by16,by24;
    if(csv_R16&&string(csv_R16).size()&&file_exists(csv_R16)){
        CsvTable t=read_analysis_csv(csv_R16,T_const);for(const auto&r:t.rows)if(pidx.count(r.phase))by16[r.spot][r.phase].push_back(r);
    }
    if(csv_R24&&string(csv_R24).size()&&file_exists(csv_R24)){
        CsvTable t=read_analysis_csv(csv_R24,T_const);for(const auto&r:t.rows)if(pidx.count(r.phase))by24[r.spot][r.phase].push_back(r);
    }

    // =====================================================================
    // SINGLE SPOT: all forced+detected values; detected=filled blue circle,
    // forced=open red square. Canvas deliberately twice as wide.
    // =====================================================================
    for(auto&skv:rep){
        int spot=skv.first;vector<double>xall,yall,xd,yd,exd,ed,xf,yf,exf,ef,xs,ys,exs,esyst;double ymin=1e99,ymax=-1e99;
        for(int i=0;i<(int)phases.size();++i){if(!skv.second.count(phases[i]))continue;const auto&r=skv.second[phases[i]];
            xall.push_back(i);yall.push_back(r.luminosity);xs.push_back(i);ys.push_back(r.luminosity);exs.push_back(0);esyst.push_back(r.deltaL);
            if(r.detected){xd.push_back(i);yd.push_back(r.luminosity);exd.push_back(0);ed.push_back(r.error);}else{xf.push_back(i);yf.push_back(r.luminosity);exf.push_back(0);ef.push_back(r.error);}
            double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));ymin=std::min(ymin,r.luminosity-emax);ymax=std::max(ymax,r.luminosity+emax);
        }
        if(xall.empty())continue;if(ymin>0)ymin*=.95;ymax*=1.08;if(ymax<=ymin)ymax=ymin+1;
        TCanvas*c=new TCanvas(Form("c_phase_single_%d",spot),"",2400,850);c->SetGrid();c->SetBottomMargin(.18);
        TH1D*frame=new TH1D(Form("frame_single_%d",spot),Form("%s - Spot %d: luminosity vs phase",full_condition.c_str(),spot),(int)phases.size(),-1.0,(double)phases.size()-1.0);
        frame->SetMinimum(ymin);frame->SetMaximum(ymax);frame->GetXaxis()->SetTitle("annealing phase");frame->GetYaxis()->SetTitle("luminosity");
        for(int i=0;i<(int)phases.size();++i)frame->GetXaxis()->SetBinLabel(i+1,phase_short_label(phases[i]).c_str());
        frame->GetXaxis()->LabelsOption("h");frame->GetXaxis()->SetLabelSize(.034);frame->GetYaxis()->SetTitleOffset(1.15);frame->Draw();
        TGraph*line=new TGraph((int)xall.size(),xall.data(),yall.data());line->SetLineColor(kP10Gray);line->SetLineWidth(2);line->Draw("L SAME");
        TGraphErrors*gy=new TGraphErrors((int)xs.size(),xs.data(),ys.data(),exs.data(),esyst.data());gy->SetLineColor(kP10Violet);gy->SetLineWidth(2);gy->SetMarkerSize(0);gy->Draw("[] SAME");
        TGraphErrors*gd=nullptr;if(!xd.empty()){gd=new TGraphErrors((int)xd.size(),xd.data(),yd.data(),exd.data(),ed.data());gd->SetLineColor(DETECTED_COLOR);gd->SetMarkerColor(DETECTED_COLOR);gd->SetMarkerStyle(20);gd->SetMarkerSize(1.25);gd->Draw("PE SAME");}
        TGraphErrors*gf=nullptr;if(!xf.empty()){gf=new TGraphErrors((int)xf.size(),xf.data(),yf.data(),exf.data(),ef.data());gf->SetLineColor(FORCED_COLOR);gf->SetMarkerColor(FORCED_COLOR);gf->SetMarkerStyle(25);gf->SetMarkerSize(1.25);gf->Draw("PE SAME");}
        TLegend*leg=new TLegend(.66,.68,.92,.90);leg->SetBorderSize(0);leg->SetFillStyle(0);if(gd)leg->AddEntry(gd,"Detected: filled blue circle","lep");if(gf)leg->AddEntry(gf,"Forced: open red square","lep");leg->AddEntry(gy,"Systematic uncertainty","l");leg->Draw();
        save_canvas_both(c,Form("%s/single_spots/%s_luminosity_vs_phase_spot%d",out.c_str(),pref.c_str(),spot));delete c;
    }

    // =====================================================================
    // GROUPED: max six spots/canvas. Marker encodes detection status, line color
    // encodes spot and never uses red/blue.
    // =====================================================================
    map<int,vector<int>>groups;for(auto&skv:rep){double m=-1;for(auto&pkv:skv.second)m=std::max(m,pkv.second.luminosity);groups[m>0?(int)std::floor(std::log10(m)):-999].push_back(skv.first);}
    int counter=0;
    for(auto&gkv:groups){auto ids=gkv.second;std::sort(ids.begin(),ids.end());for(int start=0;start<(int)ids.size();start+=MAX_SPOTS_PER_CANVAS){int end=std::min(start+MAX_SPOTS_PER_CANVAS,(int)ids.size());vector<int>sub(ids.begin()+start,ids.begin()+end);
        double ymin=1e99,ymax=-1e99;for(int id:sub)for(auto&pkv:rep[id]){const auto&r=pkv.second;double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));ymin=std::min(ymin,r.luminosity-emax);ymax=std::max(ymax,r.luminosity+emax);}if(ymin>0)ymin*=.8;else ymin*=1.2;ymax*=1.25;if(ymax<=ymin)ymax=ymin+1;
        string order=(gkv.first==-999)?"zero/negative luminosity":Form("order 10^{%d}",gkv.first);string title=full_condition+" - luminosity vs phase - "+order;
        TCanvas*c=new TCanvas(Form("c_phase_group_%d",counter),title.c_str(),1500,1100);c->SetGrid();c->SetBottomMargin(.24);
        TH1D*frame=new TH1D(Form("frame_group_%d",counter),title.c_str(),(int)phases.size(),-.5,(double)phases.size()-.5);frame->SetMinimum(ymin);frame->SetMaximum(ymax);frame->GetXaxis()->SetTitle("annealing phase");frame->GetYaxis()->SetTitle("luminosity");
        for(int i=0;i<(int)phases.size();++i)frame->GetXaxis()->SetBinLabel(i+1,phase_short_label(phases[i]).c_str());frame->GetXaxis()->LabelsOption("v");frame->GetXaxis()->SetLabelSize(.038);frame->GetYaxis()->SetTitleOffset(1.25);frame->Draw();
        TLegend*leg=new TLegend(.73,.56,.93,.90);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.027);TGraphErrors*det_ex=nullptr;TGraphErrors*for_ex=nullptr;
        int ig=0;for(int id:sub){vector<double>xall,yall,xd,yd,exd,ed,xf,yf,exf,ef,xs,ys,exs,esyst;for(int i=0;i<(int)phases.size();++i){if(!rep[id].count(phases[i]))continue;const auto&r=rep[id][phases[i]];xall.push_back(i);yall.push_back(r.luminosity);xs.push_back(i);ys.push_back(r.luminosity);exs.push_back(0);esyst.push_back(r.deltaL);if(r.detected){xd.push_back(i);yd.push_back(r.luminosity);exd.push_back(0);ed.push_back(r.error);}else{xf.push_back(i);yf.push_back(r.luminosity);exf.push_back(0);ef.push_back(r.error);}}
            if(xall.empty())continue;int col=line_colors[ig%line_colors.size()];TGraph*line=new TGraph((int)xall.size(),xall.data(),yall.data());line->SetLineColor(col);line->SetLineWidth(3);line->Draw("L SAME");
            TGraphErrors*gy=new TGraphErrors((int)xs.size(),xs.data(),ys.data(),exs.data(),esyst.data());gy->SetLineColor(col);gy->SetLineStyle(2);gy->SetLineWidth(1);gy->SetMarkerSize(0);gy->Draw("[] SAME");
            if(!xd.empty()){TGraphErrors*gd=new TGraphErrors((int)xd.size(),xd.data(),yd.data(),exd.data(),ed.data());gd->SetLineColor(DETECTED_COLOR);gd->SetMarkerColor(DETECTED_COLOR);gd->SetMarkerStyle(20);gd->SetMarkerSize(1.15);gd->Draw("PE SAME");if(!det_ex)det_ex=gd;}
            if(!xf.empty()){TGraphErrors*gf=new TGraphErrors((int)xf.size(),xf.data(),yf.data(),exf.data(),ef.data());gf->SetLineColor(FORCED_COLOR);gf->SetMarkerColor(FORCED_COLOR);gf->SetMarkerStyle(25);gf->SetMarkerSize(1.15);gf->Draw("PE SAME");if(!for_ex)for_ex=gf;}
            leg->AddEntry(line,Form("spot %d",id),"l");++ig;}
        leg->Draw();TLegend*status=new TLegend(.10,.76,.42,.90);status->SetBorderSize(0);status->SetFillStyle(0);status->SetTextSize(.030);if(det_ex)status->AddEntry(det_ex,"Detected measurement","lep");if(for_ex)status->AddEntry(for_ex,"Forced measurement","lep");status->Draw();
        save_canvas_both(c,Form("%s/grouped/%s_luminosity_vs_phase_group_%02d_order_%d",out.c_str(),pref.c_str(),counter,gkv.first));delete c;++counter;
    }}

    // =====================================================================
    // MANUALLY SELECTED GROUPS
    //
    // Each entry of MANUAL_GROUPS produces one canvas containing exactly
    // the hotspot IDs specified by the user.
    // =====================================================================
    int manual_canvas_id = 0;

    for (const auto& group : MANUAL_GROUPS) {

        const vector<int>& ids = group.ids;
        const string& manual_title = group.title;

        if (ids.empty())
        continue;

    // -------------------------------------------------------------
    // Determine y range using only the manually selected hotspots.
    // -------------------------------------------------------------
    double ymin = 1e99;
    double ymax = -1e99;
    bool has_valid_spot = false;

    for (int id : ids) {

        if (!rep.count(id)) {
            std::cerr
                << "WARNING: manually selected hotspot "
                << id
                << " is not present in the dataset."
                << std::endl;
            continue;
        }

        for (const auto& pkv : rep[id]) {

            const auto& r = pkv.second;

            double emax =
                std::max(std::fabs(r.error),
                         std::fabs(r.deltaL));

            ymin = std::min(ymin, r.luminosity - emax);
            ymax = std::max(ymax, r.luminosity + emax);

            has_valid_spot = true;
        }
    }

    if (!has_valid_spot)
        continue;

    if (ymin > 0)
        ymin *= 0.8;
    else
        ymin *= 1.2;

    ymax *= 1.25;

    if (ymax <= ymin)
        ymax = ymin + 1.0;


    // -------------------------------------------------------------
    // Canvas
    // -------------------------------------------------------------
    string title = full_condition + " - " + manual_title;

    TCanvas* c = new TCanvas(
        Form("c_phase_manual_%d", manual_canvas_id),
        title.c_str(),
        2400,
        1200
    );

    c->SetGrid();
    c->SetBottomMargin(.20);
    c->SetTopMargin(0.10);

    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.04);

    // -------------------------------------------------------------
    // Frame
    // -------------------------------------------------------------
    TH1D* frame =
        make_phase_frame(
            Form("frame_manual_%d", manual_canvas_id),
            title.c_str(),
            phases,
            phase_x,
            ymin,
            ymax
        );

    frame->GetXaxis()->SetTitle("Annealing phase");
    frame->GetYaxis()->SetTitle("Luminosity");
    frame->GetYaxis()->SetTitleOffset(1.25);

    frame->Draw();


    // -------------------------------------------------------------
    // Legends
    // -------------------------------------------------------------
    TLegend* leg = new TLegend(.73,.56,.93,.90);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(.027);

    TGraphErrors* det_ex = nullptr;
    TGraphErrors* for_ex = nullptr;


    // -------------------------------------------------------------
    // Draw manually selected hotspots.
    // -------------------------------------------------------------
        int ig = 0;

        for (int id : ids) {

            if (!rep.count(id))
                continue;

            vector<double> xall, yall;
            vector<double> xd, yd, exd, ed;
            vector<double> xf, yf, exf, ef;
            vector<double> xs, ys, exs, esyst;

            for (int i=0; i<(int)phases.size(); ++i) {

                if (!rep[id].count(phases[i]))
                    continue;

                const auto& r = rep[id][phases[i]];
                const double xphase = phase_x[i];

                xall.push_back(xphase);
                yall.push_back(r.luminosity);

                xs.push_back(xphase);
                ys.push_back(r.luminosity);
                exs.push_back(0.0);
                esyst.push_back(r.deltaL);

                if (r.detected) {
                    xd.push_back(xphase);
                    yd.push_back(r.luminosity);
                    exd.push_back(0.0);
                    ed.push_back(r.error);
                }
                else {
                    xf.push_back(xphase);
                    yf.push_back(r.luminosity);
                    exf.push_back(0.0);
                    ef.push_back(r.error);
                }
            }

            if (xall.empty())
                continue;


            int col =
                line_colors[ig % line_colors.size()];


            // Connecting line
            TGraph* line = new TGraph(
                    (int)xall.size(),
                    xall.data(),
                    yall.data());

            line->SetLineColor(col);
            line->SetLineWidth(3);
            line->Draw("L SAME");


            // Systematic uncertainty
            TGraphErrors* gy =
                new TGraphErrors(
                    (int)xs.size(),
                    xs.data(),
                    ys.data(),
                    exs.data(),
                    esyst.data()
                );

            gy->SetLineColor(col);
            gy->SetLineStyle(2);
            gy->SetLineWidth(1);
            gy->SetMarkerSize(0);
            gy->Draw("[] SAME");


            // Detected values
            if (!xd.empty()) {

                TGraphErrors* gd =
                    new TGraphErrors(
                        (int)xd.size(),
                        xd.data(),
                        yd.data(),
                        exd.data(),
                        ed.data()
                    );

                gd->SetLineColor(DETECTED_COLOR);
                gd->SetMarkerColor(DETECTED_COLOR);
                gd->SetMarkerStyle(20);
                gd->SetMarkerSize(1.15);
                gd->Draw("PE SAME");

                if (!det_ex)
                    det_ex = gd;
            }


            // Forced values
            if (!xf.empty()) {

                TGraphErrors* gf =
                    new TGraphErrors(
                        (int)xf.size(),
                        xf.data(),
                        yf.data(),
                        exf.data(),
                        ef.data()
                    );

                gf->SetLineColor(DETECTED_COLOR);
                gf->SetMarkerColor(DETECTED_COLOR);
                gf->SetMarkerStyle(25);
                gf->SetMarkerSize(1.15);
                gf->Draw("PE SAME");

                if (!for_ex)
                    for_ex = gf;
            }


            leg->AddEntry(line,Form("spot %d",id),"l");

            ++ig;
        }


        leg->Draw();

        draw_horizontal_phase_labels(
            (TPad*)gPad,
            phases,
            phase_x,
            phase_xmin,
            phase_xmax
        );
        // Detection-status legend
        //TLegend* status = new TLegend(.10,.76,.42,.90);
        /*
        status->SetBorderSize(0);
        status->SetFillStyle(0);
        status->SetTextSize(.030);

        if (det_ex)
            status->AddEntry(
                det_ex,
                "Detected measurement",
                "lep"
            );

        if (for_ex)
            status->AddEntry(
                for_ex,
                "Forced measurement",
                "lep"
            );

        status->Draw();*/


        // -------------------------------------------------------------
        // Save PNG + PDF.
        // -------------------------------------------------------------
        save_canvas_both(
            c,
            Form(
                "%s/manual_grouped/%s_%s",
                out.c_str(),
                pref.c_str(),
                safe_token(manual_title).c_str()
            )
        );

        delete c;

        ++manual_canvas_id;
    }


    // =====================================================================
    // RATIOS: statistical errors only. Produce automatic-range and fixed-focus
    // (-1,2) versions from identical data.
    // =====================================================================
    const string before="before_annealing";
    for(const auto&ph:phases){if(ph==before)continue;vector<double>x,y,ex,estat;for(auto&skv:rep){int id=skv.first;if(!skv.second.count(before)||!skv.second.count(ph))continue;const auto&a=skv.second[before];const auto&b=skv.second[ph];if(a.luminosity==0)continue;double R=b.luminosity/a.luminosity;double stat=std::sqrt(std::pow(b.error/a.luminosity,2)+std::pow(b.luminosity*a.error/(a.luminosity*a.luminosity),2));if(!finite_number(R)||!finite_number(stat))continue;x.push_back(id);y.push_back(R);ex.push_back(0);estat.push_back(stat);}if(x.empty())continue;
        double dmin=1e99,dmax=-1e99;for(size_t i=0;i<y.size();++i){dmin=std::min(dmin,y[i]-estat[i]);dmax=std::max(dmax,y[i]+estat[i]);}double span=dmax-dmin;if(!(span>0))span=std::max(.2,std::fabs(.5*(dmin+dmax))*.2);double ymin=dmin-.15*span,ymax=dmax+.15*span;if(ymax<=ymin)ymax=ymin+1;
        auto draw_ratio=[&](double lo,double hi,const string&suffix){TCanvas*c=new TCanvas(Form("c_ratio_%s_%s",safe_token(ph).c_str(),suffix.c_str()),"",1400,850);c->SetGrid();TGraphErrors*g=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),estat.data());g->SetTitle(Form("%s - L(%s) / L(before_annealing);global spot ID;ratio",condition.c_str(),phase_short_label(ph).c_str()));g->SetMarkerStyle(20);g->SetMarkerColor(kP6Blue);g->SetLineColor(kP6Blue);g->SetMarkerSize(1.2);g->SetLineWidth(2);g->Draw("AP");g->GetYaxis()->SetRangeUser(lo,hi);c->Update();double xmin=g->GetXaxis()->GetXmin(),xmax=g->GetXaxis()->GetXmax();if(1>=lo&&1<=hi){TLine*one=new TLine(xmin,1,xmax,1);one->SetLineStyle(2);one->SetLineWidth(2);one->SetLineColor(kP10Gray);one->Draw("SAME");}TLegend*leg=new TLegend(.68,.80,.92,.90);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->AddEntry(g,"Statistical uncertainty","lep");leg->Draw();save_canvas_both(c,Form("%s/ratios/%s_ratio_%s_over_before%s",out.c_str(),pref.c_str(),safe_token(ph).c_str(),suffix.c_str()));delete c;};
        draw_ratio(ymin,ymax,"");draw_ratio(-1.0,2.0,"_focus_m1_2");
    }

    // =====================================================================
    // B/lambda VS PHASE: six hotspots per 3x2 canvas, explicit phase labels,
    // statistical + R16/R24 systematic uncertainty.
    // =====================================================================
    vector<int>ids;for(auto&skv:all_by_spot)ids.push_back(skv.first);std::sort(ids.begin(),ids.end());int canvas_id=0;
    for(int start=0;start<(int)ids.size();start+=6){TCanvas*c=new TCanvas(Form("c_parameter_phase_%d",canvas_id),T_const?"B vs phase":"lambda vs phase",2400,1200);c->Divide(3,2);bool any=false;
        for(int j=0;j<6&&start+j<(int)ids.size();++j){int id=ids[start+j];c->cd(j+1);gPad->SetBottomMargin(.28);vector<double>x,y,ex,estat,esyst;
            for(int i=0;i<(int)phases.size();++i){auto it=all_by_spot[id].find(phases[i]);if(it==all_by_spot[id].end())continue;PhaseFitPoint p20=fit_phase_parameter(it->second,T_const,has_vfin);if(!p20.ok)continue;double syst=0;bool hs=false;
                if(by16.count(id)&&by16[id].count(phases[i])&&by24.count(id)&&by24[id].count(phases[i])){vector<AnalysisRow>r20,r16,r24;common_radius_rows(it->second,by16[id][phases[i]],by24[id][phases[i]],T_const,r20,r16,r24);PhaseFitPoint q20=fit_phase_parameter(r20,T_const,has_vfin),q16=fit_phase_parameter(r16,T_const,has_vfin),q24=fit_phase_parameter(r24,T_const,has_vfin);if(q20.ok&&q16.ok&&q24.ok){syst=std::max(std::fabs(q24.val-q20.val),std::fabs(q20.val-q16.val));hs=true;p20=q20;}}
                x.push_back(i);y.push_back(p20.val);ex.push_back(0);estat.push_back(p20.err);esyst.push_back(hs?syst:0.0);}
            if(x.empty())continue;any=true;double ymin=1e99,ymax=-1e99;for(size_t k=0;k<y.size();++k){double e=std::max(estat[k],esyst[k]);ymin=std::min(ymin,y[k]-e);ymax=std::max(ymax,y[k]+e);}double sp=ymax-ymin;if(!(sp>0))sp=std::max(.1,std::fabs(ymax)*.2);ymin-=.15*sp;ymax+=.30*sp;
            TH1D*frame=new TH1D(Form("frame_par_%d_%d",canvas_id,id),Form("%s - Spot %d;%s;%s",condition.c_str(),id,"annealing phase",T_const?"B":"#lambda"),(int)phases.size(),-.5,(double)phases.size()-.5);frame->SetMinimum(ymin);frame->SetMaximum(ymax);for(int i=0;i<(int)phases.size();++i)frame->GetXaxis()->SetBinLabel(i+1,phase_short_label(phases[i]).c_str());frame->GetXaxis()->LabelsOption("v");frame->GetXaxis()->SetLabelSize(.055);frame->Draw();
            TGraphErrors*gs=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),estat.data());gs->SetMarkerStyle(20);gs->SetMarkerColor(kBlack);gs->SetLineColor(kBlack);gs->SetMarkerSize(1);gs->SetLineWidth(2);gs->Draw("P SAME");TGraphErrors*gy=new TGraphErrors((int)x.size(),x.data(),y.data(),ex.data(),esyst.data());gy->SetLineColor(kP6Grape);gy->SetLineWidth(2);gy->SetMarkerSize(0);gy->Draw("[] SAME");gs->Draw("P SAME");TLegend*leg=new TLegend(.12,.72,.58,.89);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.035);leg->AddEntry(gs,"Statistical uncertainty","lep");leg->AddEntry(gy,"Systematic uncertainty (R16/R24)","l");leg->Draw();}
        if(any)save_canvas_both(c,Form("%s/fit_parameter_vs_phase/%s_%s_vs_phase_canvas_%02d",out.c_str(),pref.c_str(),T_const?"B":"lambda",canvas_id));delete c;++canvas_id;
    }
}

#endif
