#include "gold_analysis_common.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

struct GoldBSyst { bool ok=false; double deltaB=0.0; };
struct GoldVFit {
    int spot=-1;
    vector<GoldRow> rows;
    double A=0.0,Aerr=0.0,B=0.0,Berr=0.0,chi2=0.0,prob=0.0,edm=-1.0;
    int ndf=0,status=-999,covstatus=-999;
    bool fitted=false,converged=false,has_syst=false;
    double deltaB=0.0;
};

static map<int,GoldBSyst> read_gold_B_syst(const string& filename) {
    map<int,GoldBSyst> out;
    std::ifstream fin(filename); if (!fin.is_open()) return out;
    string line; if (!std::getline(fin,line)) return out;
    auto h=gold_split_csv(line); map<string,int> c; for (int i=0;i<(int)h.size();++i)c[h[i]]=i;
    while (std::getline(fin,line)) {
        if (gold_trim(line).empty()) continue;
        auto f=gold_split_csv(line); if (f.size()<h.size()) f.resize(h.size(),"");
        try {
            const int spot=(int)std::llround(std::stod(f[c["spot"]]));
            GoldBSyst s;
            if (c.count("deltaB") && !f[c["deltaB"]].empty()) {
                s.deltaB=std::fabs(std::stod(f[c["deltaB"]]));
                s.ok=gold_finite(s.deltaB);
            }
            out[spot]=s;
        } catch (...) {}
    }
    return out;
}

void gold_lum_vs_v_fit(const char* csv_R20,
                       const char* systematic_csv,
                       const char* output_dir) {
    gStyle->SetOptStat(0);
    const string outdir=output_dir;
    gold_ensure_dir(outdir);
    gold_ensure_dir(outdir+"/lum_vs_overvoltage_all_spots");

    GoldTable table=read_gold_csv(csv_R20,true);
    auto selected=gold_select_dataset(table.rows,"T20",true);
    if (selected.empty()) { std::cerr << "No detected rows for dataset T20" << std::endl; return; }
    auto byspot=gold_group_spots(selected);
    auto syst=read_gold_B_syst(systematic_csv);
    map<int,GoldVFit> fits;

    for (auto& kv:byspot) {
        const int spot=kv.first;
        auto rows=kv.second;
        std::sort(rows.begin(),rows.end(),[](const GoldRow&a,const GoldRow&b){return a.v_fin<b.v_fin;});
        GoldVFit info; info.spot=spot; info.rows=rows;
        vector<double> xv,yv,exv,eyv,esyst;
        for (const auto& r:rows) {
            if (!gold_finite(r.v_fin)) continue;
            xv.push_back(r.v_fin); yv.push_back(r.luminosity); exv.push_back(0.0); eyv.push_back(r.error); esyst.push_back(r.deltaL);
        }
        if (xv.size()<2) { fits[spot]=info; continue; }

        TGraphErrors gr((int)xv.size(),xv.data(),yv.data(),exv.data(),eyv.data());
        TF1 fitf(Form("gold_v_fit_spot_%d",spot),"[0]*TMath::Power(x,[1])",0.0,8.0);
        fitf.SetParameters(1.0,2.0); fitf.SetParNames("A","B");
        TFitResultPtr fit=gr.Fit(&fitf,"QRS0");
        info.fitted=true; info.status=(int)fit;
        if (fit.Get()) { info.covstatus=fit->CovMatrixStatus(); info.edm=fit->Edm(); }
        info.A=fitf.GetParameter(0); info.Aerr=fitf.GetParError(0);
        info.B=fitf.GetParameter(1); info.Berr=fitf.GetParError(1);
        info.chi2=fitf.GetChisquare(); info.ndf=fitf.GetNDF(); info.prob=fitf.GetProb();
        info.converged=(info.status==0 && gold_finite(info.B));
        if (syst.count(spot) && syst[spot].ok) { info.has_syst=true; info.deltaB=syst[spot].deltaB; }
        fits[spot]=info;

        TCanvas* c=new TCanvas(Form("gold_v_spot_%d",spot),"",1600,1200);
        c->SetLeftMargin(.12); c->SetRightMargin(.05); c->SetBottomMargin(.13); c->SetTopMargin(.08);
        gr.SetTitle(Form("A1, T = 20 #circC - Spot %d: x = %.2f, y = %.2f;Overvoltage (V);Luminosity",
                         spot,rows.front().x,rows.front().y));
        gr.SetMarkerStyle(20); gr.SetMarkerSize(1.8); gr.SetMarkerColor(kBlack); gr.SetLineColor(kBlack); gr.SetLineWidth(2);
        gr.GetXaxis()->SetTitleSize(.045); gr.GetYaxis()->SetTitleSize(.045); gr.GetYaxis()->SetTitleOffset(1.4);
        double ymin=1e99,ymax=-1e99;
        for (size_t i=0;i<yv.size();++i) {
            double e=std::max(std::fabs(eyv[i]),std::fabs(esyst[i]));
            ymin=std::min(ymin,yv[i]-e); ymax=std::max(ymax,yv[i]+e);
        }
        if (gold_finite(ymin)&&gold_finite(ymax)&&ymax>ymin) {
            double span=ymax-ymin; gr.SetMinimum(ymin-.08*span); gr.SetMaximum(ymax+.18*span);
        }
        gr.Draw("AP"); c->Update();
        gold_draw_syst_brackets(xv,yv,esyst,gr.GetXaxis()->GetXmin(),gr.GetXaxis()->GetXmax(),kBlack,4,.009);
        gr.Draw("P SAME");
        fitf.SetLineColor(kBlack); fitf.SetLineWidth(2); fitf.Draw("SAME");

        TLine* sysproxy=new TLine(0,0,1,0); sysproxy->SetLineColor(kBlack); sysproxy->SetLineWidth(4);
        TLegend* leg=new TLegend(.12,.69,.66,.91); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(.030);
        leg->AddEntry(&gr,"Data (statistical uncertainty)","lep");
        leg->AddEntry(sysproxy,"Luminosity systematic uncertainty","l");
        string flabel=info.has_syst
            ? Form("Fit: L = A V_{over}^{B}, B = %.3f #pm %.3f (stat) #pm %.3f (syst)",info.B,info.Berr,info.deltaB)
            : Form("Fit: L = A V_{over}^{B}, B = %.3f #pm %.3f",info.B,info.Berr);
        leg->AddEntry(&fitf,flabel.c_str(),"l"); leg->Draw();
        gold_save_canvas(c,Form("%s/lum_vs_overvoltage_all_spots/A1_T=20_lum_vs_overvoltage_spot%d",outdir.c_str(),spot));
        delete c;
    }

    vector<double> x,ex,B,Bstat,Bsyst;
    for (const auto& kv:fits) {
        const auto& f=kv.second;
        if (!f.converged || !gold_finite(f.Berr)) continue;
        x.push_back(kv.first); ex.push_back(0.0); B.push_back(f.B); Bstat.push_back(f.Berr); Bsyst.push_back(f.has_syst?f.deltaB:0.0);
    }
    if (!B.empty()) {
        double xmin=*std::min_element(x.begin(),x.end())-1.0;
        double xmax=*std::max_element(x.begin(),x.end())+1.0;
        double ymin=2.0,ymax=2.0;
        for(size_t i=0;i<B.size();++i){double e=std::max(std::fabs(Bstat[i]),std::fabs(Bsyst[i]));ymin=std::min(ymin,B[i]-e);ymax=std::max(ymax,B[i]+e);}
        double span=ymax-ymin; if(!(span>0.0)) span=.4;
        TCanvas* c=new TCanvas("gold_B_vs_spot","B vs spot ID",1800,1000);
        TGraphErrors* gr=new TGraphErrors((int)B.size(),x.data(),B.data(),ex.data(),Bstat.data());
        gr->SetTitle("A1, T = 20 #circC - Power-law exponent B;Global spot ID;B");
        gr->SetMarkerStyle(21); gr->SetMarkerSize(1.2); gr->SetMarkerColor(kBlack); gr->SetLineColor(kBlack); gr->SetLineWidth(2);
        gr->SetMinimum(ymin-.18*span); gr->SetMaximum(ymax+.28*span); gr->Draw("AP"); gr->GetXaxis()->SetLimits(xmin,xmax); c->Update();
        gold_draw_syst_brackets(x,B,Bsyst,xmin,xmax,kBlack,4,.009); gr->Draw("P SAME");
        TLine* line2=new TLine(xmin,2.0,xmax,2.0); line2->SetLineColor(kP6Red); line2->SetLineStyle(2); line2->SetLineWidth(3); line2->Draw("SAME");
        TLine* sysproxy=new TLine(0,0,1,0); sysproxy->SetLineColor(kBlack); sysproxy->SetLineWidth(4);
        TLegend* leg=new TLegend(.12,.72,.48,.91); leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(gr,"B statistical uncertainty","lep");
        leg->AddEntry(sysproxy,"B systematic uncertainty","l");
        leg->AddEntry(line2,"B = 2","l"); leg->Draw();
        gold_save_canvas(c,outdir+"/A1_T=20_B_vs_spot"); delete c;
    }

    std::ofstream fout(outdir+"/A1_T=20_B_fit_results.csv");
    fout << "spot,x,y,A,A_stat_error,B,B_stat_error,deltaB,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";
    for (const auto& kv:fits) {
        const auto& f=kv.second; if (f.rows.empty()) continue;
        fout << kv.first << ',' << f.rows.front().x << ',' << f.rows.front().y << ','
             << f.A << ',' << f.Aerr << ',' << f.B << ',' << f.Berr << ',' << (f.has_syst?f.deltaB:0.0) << ','
             << f.chi2 << ',' << f.ndf << ',' << (f.ndf>0?f.chi2/f.ndf:0.0) << ',' << f.prob << ','
             << f.status << ',' << f.covstatus << ',' << f.edm << ',' << (f.converged?1:0) << '\n';
    }
}
