#include "gold_analysis_common.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

struct GoldLambdaSyst { bool ok=false; double delta=0.0; };
struct GoldTFit {
    int spot=-1;
    vector<GoldRow> rows;
    double A=0.0,Aerr=0.0,lambda=0.0,lambdaerr=0.0,chi2=0.0,prob=0.0,edm=-1.0;
    int ndf=0,status=-999,covstatus=-999;
    bool fitted=false,converged=false,has_syst=false;
    double deltaLambda=0.0;
};

static map<int,GoldLambdaSyst> read_gold_lambda_syst(const string& filename) {
    map<int,GoldLambdaSyst> out;
    std::ifstream fin(filename); if (!fin.is_open()) return out;
    string line; if (!std::getline(fin,line)) return out;
    auto h=gold_split_csv(line); map<string,int> c; for(int i=0;i<(int)h.size();++i)c[h[i]]=i;
    while(std::getline(fin,line)){
        if(gold_trim(line).empty())continue; auto f=gold_split_csv(line); if(f.size()<h.size())f.resize(h.size(),"");
        try{int id=(int)std::llround(std::stod(f[c["spot"]]));GoldLambdaSyst s;if(c.count("deltaLambda")&&!f[c["deltaLambda"]].empty()){s.delta=std::fabs(std::stod(f[c["deltaLambda"]]));s.ok=gold_finite(s.delta);}out[id]=s;}catch(...){}
    }
    return out;
}

void gold_lum_vs_T_fit(const char* csv_R20,
                       const char* dataset_key,
                       const char* systematic_csv,
                       const char* output_dir) {
    gStyle->SetOptStat(0);
    const string key=dataset_key, outdir=output_dir;
    gold_ensure_dir(outdir); gold_ensure_dir(outdir+"/lum_vs_T_all_spots_expfit");
    GoldTable table=read_gold_csv(csv_R20,true);
    auto selected=gold_select_dataset(table.rows,key,true);
    if(selected.empty()){std::cerr<<"No detected rows for dataset "<<key<<std::endl;return;}
    auto byspot=gold_group_spots(selected); auto syst=read_gold_lambda_syst(systematic_csv); map<int,GoldTFit> fits;

    double vlabel=selected.front().v;
    for(auto&kv:byspot){
        int spot=kv.first;auto rows=kv.second;std::sort(rows.begin(),rows.end(),[](const GoldRow&a,const GoldRow&b){return a.T<b.T;});
        GoldTFit info;info.spot=spot;info.rows=rows;int n=(int)rows.size();if(n<3){fits[spot]=info;continue;}
        vector<double>xv(n),yv(n),exv(n,0.0),eyv(n),esyst(n);for(int i=0;i<n;++i){xv[i]=rows[i].T;yv[i]=rows[i].luminosity;eyv[i]=rows[i].error;esyst[i]=rows[i].deltaL;}
        double Tmin=rows.front().T,Tmax=rows.back().T,dT=Tmax-Tmin;if(!(dT>0.0))dT=1.0;double A0,l0;gold_estimate_exp(rows,A0,l0);
        TGraphErrors gr(n,xv.data(),yv.data(),exv.data(),eyv.data());TF1 fitf(Form("gold_T_fit_%s_%d",key.c_str(),spot),"[0]*exp([1]*x)",Tmin-.10*dT,Tmax+.10*dT);
        fitf.SetParameters(A0,l0);fitf.SetParNames("A","#lambda");TFitResultPtr fit=gr.Fit(&fitf,"QRS0");
        info.fitted=true;info.status=(int)fit;if(fit.Get()){info.covstatus=fit->CovMatrixStatus();info.edm=fit->Edm();}
        info.A=fitf.GetParameter(0);info.Aerr=fitf.GetParError(0);info.lambda=fitf.GetParameter(1);info.lambdaerr=fitf.GetParError(1);info.chi2=fitf.GetChisquare();info.ndf=fitf.GetNDF();info.prob=fitf.GetProb();info.converged=(info.status==0&&gold_finite(info.lambda));
        if(syst.count(spot)&&syst[spot].ok){info.has_syst=true;info.deltaLambda=syst[spot].delta;}fits[spot]=info;

        TCanvas*c=new TCanvas(Form("gold_T_%s_%d",key.c_str(),spot),"",1400,950);c->SetLeftMargin(.12);c->SetBottomMargin(.12);
        gr.SetTitle(Form("A1, v = %.1f V - Spot %d: x = %.2f, y = %.2f;T (#circC);Luminosity",vlabel,spot,rows.front().x,rows.front().y));
        gr.SetMarkerStyle(20);gr.SetMarkerSize(1.4);gr.SetMarkerColor(kBlack);gr.SetLineColor(kBlack);gr.SetLineWidth(2);gr.GetYaxis()->SetTitleOffset(1.35);
        double ymin=1e99,ymax=-1e99;for(int i=0;i<n;++i){double e=std::max(std::fabs(eyv[i]),std::fabs(esyst[i]));ymin=std::min(ymin,yv[i]-e);ymax=std::max(ymax,yv[i]+e);}if(ymin>0)ymin*=.85;else ymin-=.08*(ymax-ymin);ymax*=1.20;gr.SetMinimum(ymin);gr.SetMaximum(ymax);gr.Draw("AP");c->Update();
        gold_draw_syst_brackets(xv,yv,esyst,gr.GetXaxis()->GetXmin(),gr.GetXaxis()->GetXmax(),kBlack,4,.0045);gr.Draw("P SAME");fitf.SetLineColor(kBlack);fitf.SetLineWidth(2);fitf.Draw("SAME");
        TLine*sysproxy=new TLine(0,0,1,0);sysproxy->SetLineColor(kBlack);sysproxy->SetLineWidth(4);TLegend*leg=new TLegend(.14,.66,.68,.89);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.030);leg->AddEntry(&gr,"Data (statistical uncertainty)","lep");leg->AddEntry(sysproxy,"Luminosity systematic uncertainty","l");
        string flabel=info.has_syst?Form("Fit: L = A e^{#lambda T}, #lambda = %.5f #pm %.5f (stat) #pm %.5f (syst)",info.lambda,info.lambdaerr,info.deltaLambda):Form("Fit: L = A e^{#lambda T}, #lambda = %.5f #pm %.5f",info.lambda,info.lambdaerr);leg->AddEntry(&fitf,flabel.c_str(),"l");leg->Draw();
        gold_save_canvas(c,Form("%s/lum_vs_T_all_spots_expfit/A1_%s_lum_vs_T_spot%d",outdir.c_str(),key.c_str(),spot));delete c;
    }

    vector<double>x,ex,L,Lstat,Lsyst;for(const auto&kv:fits){const auto&f=kv.second;if(!f.converged||!gold_finite(f.lambdaerr))continue;x.push_back(kv.first);ex.push_back(0);L.push_back(f.lambda);Lstat.push_back(f.lambdaerr);Lsyst.push_back(f.has_syst?f.deltaLambda:0);}
    if(!L.empty()){
        double xmin=*std::min_element(x.begin(),x.end())-1.0,xmax=*std::max_element(x.begin(),x.end())+1.0,lmin=1e99,lmax=-1e99;for(size_t i=0;i<L.size();++i){double e=std::max(std::fabs(Lstat[i]),std::fabs(Lsyst[i]));lmin=std::min(lmin,L[i]-e);lmax=std::max(lmax,L[i]+e);}double sp=lmax-lmin;if(!(sp>0))sp=std::max(.01,std::fabs(lmax)*.2);
        TCanvas*c=new TCanvas(Form("gold_lambda_%s",key.c_str()),"lambda vs spot ID",1800,1000);TGraphErrors*gr=new TGraphErrors((int)L.size(),x.data(),L.data(),ex.data(),Lstat.data());gr->SetTitle(Form("A1, v = %.1f V - Exponential parameter #lambda;Global spot ID;#lambda",vlabel));gr->SetMarkerStyle(21);gr->SetMarkerSize(1.1);gr->SetMarkerColor(kBlack);gr->SetLineColor(kBlack);gr->SetLineWidth(2);gr->SetMinimum(lmin-.18*sp);gr->SetMaximum(lmax+.55*sp);gr->Draw("AP");gr->GetXaxis()->SetLimits(xmin,xmax);c->Update();
        gold_draw_syst_brackets(x,L,Lsyst,xmin,xmax,kP6Grape,4,.009);gr->Draw("P SAME");
        TF1*fc=new TF1(Form("gold_lambda_const_%s",key.c_str()),"[0]",xmin,xmax);fc->SetParNames("#lambda_{const}");fc->SetLineColor(kP10Gray);fc->SetLineWidth(2);bool fitted=false;if(L.size()>=2){gr->Fit(fc,"QRS");fc->Draw("SAME");fitted=true;}
        TLine*sysproxy=new TLine(0,0,1,0);sysproxy->SetLineColor(kP6Grape);sysproxy->SetLineWidth(4);TLegend*leg=new TLegend(.10,.68,.58,.91);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->AddEntry(gr,"#lambda statistical uncertainty","lep");leg->AddEntry(sysproxy,"#lambda systematic uncertainty","l");if(fitted)leg->AddEntry(fc,Form("Stat-only constant fit: #lambda = %.5f #pm %.5f",fc->GetParameter(0),fc->GetParError(0)),"l");leg->Draw();
        gold_save_canvas(c,outdir+"/A1_"+key+"_lambda_vs_spot");delete c;
    }

    std::ofstream fout(outdir+"/A1_"+key+"_lambda_fit_results.csv");fout<<"spot,x,y,v,A,A_stat_error,lambda,lambda_stat_error,deltaLambda,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";for(const auto&kv:fits){const auto&f=kv.second;if(f.rows.empty())continue;fout<<kv.first<<','<<f.rows.front().x<<','<<f.rows.front().y<<','<<f.rows.front().v<<','<<f.A<<','<<f.Aerr<<','<<f.lambda<<','<<f.lambdaerr<<','<<(f.has_syst?f.deltaLambda:0)<<','<<f.chi2<<','<<f.ndf<<','<<(f.ndf>0?f.chi2/f.ndf:0)<<','<<f.prob<<','<<f.status<<','<<f.covstatus<<','<<f.edm<<','<<(f.converged?1:0)<<'\n';}
}
