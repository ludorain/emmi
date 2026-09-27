#include "gold_analysis_common.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

struct GoldSummaryLambda { double lambda=0,stat=0,syst=0; bool ok=false; };

static map<int,GoldSummaryLambda> read_gold_lambda_results(const string& filename){
    map<int,GoldSummaryLambda>out;std::ifstream fin(filename);if(!fin.is_open())return out;string line;if(!std::getline(fin,line))return out;auto h=gold_split_csv(line);map<string,int>c;for(int i=0;i<(int)h.size();++i)c[h[i]]=i;
    while(std::getline(fin,line)){if(gold_trim(line).empty())continue;auto f=gold_split_csv(line);if(f.size()<h.size())f.resize(h.size(),"");try{int id=(int)std::llround(std::stod(f[c["spot"]]));GoldSummaryLambda q;q.lambda=std::stod(f[c["lambda"]]);q.stat=std::fabs(std::stod(f[c["lambda_stat_error"]]));q.syst=c.count("deltaLambda")&&!f[c["deltaLambda"]].empty()?std::fabs(std::stod(f[c["deltaLambda"]])):0.0;bool conv=!c.count("converged")||gold_parse_bool(f[c["converged"]]);q.ok=conv&&gold_finite(q.lambda)&&gold_finite(q.stat);out[id]=q;}catch(...){}}
    return out;
}

struct OverlaySeries { vector<GoldRow> rows; TGraphErrors* gr=nullptr; TF1* fit=nullptr; bool fit_ok=false; double lambda=0,lambdaerr=0; };
static OverlaySeries make_overlay_T(vector<GoldRow> rows,int spot,const string&tag,int color,int marker){
    OverlaySeries s;s.rows=rows;if(rows.empty())return s;std::sort(s.rows.begin(),s.rows.end(),[](const GoldRow&a,const GoldRow&b){return a.T<b.T;});int n=s.rows.size();vector<double>x(n),y(n),ex(n,0),ey(n);for(int i=0;i<n;++i){x[i]=s.rows[i].T;y[i]=s.rows[i].luminosity;ey[i]=s.rows[i].error;}s.gr=new TGraphErrors(n,x.data(),y.data(),ex.data(),ey.data());s.gr->SetMarkerStyle(marker);s.gr->SetMarkerSize(1.35);s.gr->SetMarkerColor(color);s.gr->SetLineColor(color);s.gr->SetLineWidth(2);
    if(n>=3){double Tmin=s.rows.front().T,Tmax=s.rows.back().T,dT=Tmax-Tmin;if(!(dT>0))dT=1;double A0,l0;gold_estimate_exp(s.rows,A0,l0);s.fit=new TF1(Form("gold_overlay_%s_%d",tag.c_str(),spot),"[0]*exp([1]*x)",Tmin-.10*dT,Tmax+.10*dT);s.fit->SetParameters(A0,l0);s.fit->SetParNames("A","#lambda");s.fit->SetLineColor(color);s.fit->SetLineWidth(2);TFitResultPtr r=s.gr->Fit(s.fit,"QRS0");s.lambda=s.fit->GetParameter(1);s.lambdaerr=s.fit->GetParError(1);s.fit_ok=((int)r==0&&gold_finite(s.lambda));}
    return s;
}

void gold_compare_v7_v9(const char* csv_R20,
                        const char* lambda_v7_csv,
                        const char* lambda_v9_csv,
                        const char* output_dir){
    gStyle->SetOptStat(0);
    string out=output_dir;
    gold_ensure_dir(out);
    gold_ensure_dir(out+"/lum_vs_T_v7_v9_all_spots");
    GoldTable table=read_gold_csv(csv_R20,true);
    auto m7=gold_group_spots(gold_select_dataset(table.rows,"v7",true));
    auto m9=gold_group_spots(gold_select_dataset(table.rows,"v9",true));
    std::set<int>ids;
    for(auto&kv:m7)ids.insert(kv.first);
    for(auto&kv:m9)ids.insert(kv.first);
    for(int id:ids){OverlaySeries s7=make_overlay_T(m7.count(id)?m7[id]:vector<GoldRow>{},id,"v7",kP6Red,20);OverlaySeries s9=make_overlay_T(m9.count(id)?m9[id]:vector<GoldRow>{},id,"v9",kP6Blue,21);if(!s7.gr&&!s9.gr)continue;TCanvas*c=new TCanvas(Form("gold_compare_T_%d",id),"",1700,1200);c->SetLeftMargin(.12);c->SetBottomMargin(.12);
        double ymin=1e99,ymax=-1e99,Tmin=1e99,Tmax=-1e99;auto scan=[&](const OverlaySeries&s){for(const auto&r:s.rows){double e=std::max(std::fabs(r.error),std::fabs(r.deltaL));ymin=std::min(ymin,r.luminosity-e);ymax=std::max(ymax,r.luminosity+e);Tmin=std::min(Tmin,r.T);Tmax=std::max(Tmax,r.T);}};scan(s7);scan(s9);double ysp=ymax-ymin;if(!(ysp>0))ysp=std::max(1.0,std::fabs(ymax)*.2);double xsp=Tmax-Tmin;if(!(xsp>0))xsp=1.0;
        TGraphErrors*base=s7.gr?s7.gr:s9.gr;base->SetTitle(Form("A1 - Spot %d: luminosity vs temperature at v = 7 V and v = 9 V;T (#circC);Luminosity",id));base->SetMinimum(ymin-.10*ysp);base->SetMaximum(ymax+.20*ysp);base->Draw("AP");base->GetXaxis()->SetLimits(Tmin-.08*xsp,Tmax+.08*xsp);c->Update();
        auto draw_series=[&](OverlaySeries&s,int color){if(!s.gr)return;vector<double>x,y,d;for(const auto&r:s.rows){x.push_back(r.T);y.push_back(r.luminosity);d.push_back(r.deltaL);}gold_draw_syst_brackets(x,y,d,base->GetXaxis()->GetXmin(),base->GetXaxis()->GetXmax(),color,3,.0045);s.gr->Draw("P SAME");if(s.fit)s.fit->Draw("SAME");};draw_series(s7,kP6Red);draw_series(s9,kP6Blue);
        TLine*redproxy=new TLine(0,0,1,0);redproxy->SetLineColor(kP6Red);redproxy->SetLineWidth(3);TLine*blueproxy=new TLine(0,0,1,0);blueproxy->SetLineColor(kP6Blue);blueproxy->SetLineWidth(3);TLegend*leg=new TLegend(.12,.63,.68,.91);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.029);if(s7.gr)leg->AddEntry(s7.gr,"v = 7 V data (statistical uncertainty)","lep");if(s9.gr)leg->AddEntry(s9.gr,"v = 9 V data (statistical uncertainty)","lep");if(s7.gr)leg->AddEntry(redproxy,"v = 7 V luminosity systematic uncertainty","l");if(s9.gr)leg->AddEntry(blueproxy,"v = 9 V luminosity systematic uncertainty","l");if(s7.fit&&s7.fit_ok)leg->AddEntry(s7.fit,Form("v = 7 V fit: #lambda = %.5f #pm %.5f (stat)",s7.lambda,s7.lambdaerr),"l");if(s9.fit&&s9.fit_ok)leg->AddEntry(s9.fit,Form("v = 9 V fit: #lambda = %.5f #pm %.5f (stat)",s9.lambda,s9.lambdaerr),"l");leg->Draw();gold_save_canvas(c,Form("%s/lum_vs_T_v7_v9_all_spots/A1_lum_vs_T_v7_v9_spot%d",out.c_str(),id));delete c;
    }

    auto l7=read_gold_lambda_results(lambda_v7_csv);auto l9=read_gold_lambda_results(lambda_v9_csv);vector<double>x7,y7,e7,s7,x9,y9,e9,s9,ex7,ex9;for(auto&kv:l7)if(kv.second.ok){x7.push_back(kv.first);y7.push_back(kv.second.lambda);e7.push_back(kv.second.stat);s7.push_back(kv.second.syst);ex7.push_back(0);}for(auto&kv:l9)if(kv.second.ok){x9.push_back(kv.first);y9.push_back(kv.second.lambda);e9.push_back(kv.second.stat);s9.push_back(kv.second.syst);ex9.push_back(0);}if(x7.empty()&&x9.empty())return;
    double xmin=1e99,xmax=-1e99,ymin=1e99,ymax=-1e99;auto scanpar=[&](const vector<double>&x,const vector<double>&y,const vector<double>&estat,const vector<double>&esyst){for(size_t i=0;i<x.size();++i){xmin=std::min(xmin,x[i]);xmax=std::max(xmax,x[i]);double e=std::max(std::fabs(estat[i]),std::fabs(esyst[i]));ymin=std::min(ymin,y[i]-e);ymax=std::max(ymax,y[i]+e);}};scanpar(x7,y7,e7,s7);scanpar(x9,y9,e9,s9);xmin-=1;xmax+=1;double sp=ymax-ymin;if(!(sp>0))sp=std::max(.01,std::fabs(ymax)*.2);
    TCanvas*c=new TCanvas("gold_lambda_compare","lambda comparison",1800,1200);TGraphErrors*g7=x7.empty()?nullptr:new TGraphErrors(x7.size(),x7.data(),y7.data(),ex7.data(),e7.data());TGraphErrors*g9=x9.empty()?nullptr:new TGraphErrors(x9.size(),x9.data(),y9.data(),ex9.data(),e9.data());TGraphErrors*base=g7?g7:g9;base->SetTitle("A1 - Exponential parameter #lambda at v = 7 V and v = 9 V;Global spot ID;#lambda");base->SetMinimum(ymin-.18*sp);base->SetMaximum(ymax+.55*sp);base->Draw("AP");base->GetXaxis()->SetLimits(xmin,xmax);c->Update();
    TF1*f7=nullptr,*f9=nullptr;if(g7){g7->SetMarkerStyle(20);g7->SetMarkerSize(1.15);g7->SetMarkerColor(kP6Red);g7->SetLineColor(kP6Red);g7->SetLineWidth(2);gold_draw_syst_brackets(x7,y7,s7,xmin,xmax,kP6Red,3,.009);g7->Draw("P SAME");if(x7.size()>=2){f7=new TF1("gold_lambda_v7_const","[0]",xmin,xmax);f7->SetLineColor(kP6Red);f7->SetLineStyle(2);f7->SetLineWidth(2);g7->Fit(f7,"QRS");f7->Draw("SAME");}}
    if(g9){g9->SetMarkerStyle(21);g9->SetMarkerSize(1.15);g9->SetMarkerColor(kP6Blue);g9->SetLineColor(kP6Blue);g9->SetLineWidth(2);gold_draw_syst_brackets(x9,y9,s9,xmin,xmax,kP6Blue,3,.009);g9->Draw("P SAME");if(x9.size()>=2){f9=new TF1("gold_lambda_v9_const","[0]",xmin,xmax);f9->SetLineColor(kP6Blue);f9->SetLineStyle(2);f9->SetLineWidth(2);g9->Fit(f9,"QRS");f9->Draw("SAME");}}
    TLine*rp=new TLine(0,0,1,0);rp->SetLineColor(kP6Red);rp->SetLineWidth(3);TLine*bp=new TLine(0,0,1,0);bp->SetLineColor(kP6Blue);bp->SetLineWidth(3);TLegend*leg=new TLegend(.10,.64,.62,.91);leg->SetBorderSize(0);leg->SetFillStyle(0);if(g7)leg->AddEntry(g7,"v = 7 V: statistical uncertainty","lep");if(g9)leg->AddEntry(g9,"v = 9 V: statistical uncertainty","lep");if(g7)leg->AddEntry(rp,"v = 7 V: systematic uncertainty","l");if(g9)leg->AddEntry(bp,"v = 9 V: systematic uncertainty","l");if(f7)leg->AddEntry(f7,Form("v = 7 V stat-only constant fit: %.5f #pm %.5f",f7->GetParameter(0),f7->GetParError(0)),"l");if(f9)leg->AddEntry(f9,Form("v = 9 V stat-only constant fit: %.5f #pm %.5f",f9->GetParameter(0),f9->GetParError(0)),"l");leg->Draw();gold_save_canvas(c,out+"/A1_lambda_vs_spot_v7_v9");delete c;
}
