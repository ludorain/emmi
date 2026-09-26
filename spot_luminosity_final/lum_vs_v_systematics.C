#include "analysis_common.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

struct VFitResult {
    int spot=-1;
    double B=std::numeric_limits<double>::quiet_NaN();
    double Berr=std::numeric_limits<double>::quiet_NaN();
    double chi2=std::numeric_limits<double>::quiet_NaN();
    double chi2ndf=std::numeric_limits<double>::quiet_NaN();
    double edm=std::numeric_limits<double>::quiet_NaN();
    int ndf=0,status=-999,covstatus=-999;
    bool fitted=false,converged=false;
};

static const AnalysisRow* find_v_point(const vector<AnalysisRow>& rows,double x,double tol=1e-9) {
    for (const auto& r:rows) if (std::fabs(r.v_fin-x)<=tol) return &r;
    return nullptr;
}

static VFitResult fit_v_rows(const vector<AnalysisRow>& rows,const string& tag,int spot) {
    VFitResult fr; fr.spot=spot;
    if (rows.size()<2) return fr;
    vector<double>x,y,ex,ey;
    for(const auto&r:rows){x.push_back(r.v_fin);y.push_back(r.luminosity);ex.push_back(0.0);ey.push_back(r.error);}
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
    TF1 f(Form("f_sys_v_%s_%d",tag.c_str(),spot),"[0]*TMath::Power(x,[1])",0.0,8.0);
    f.SetParameters(1.0,2.0); f.SetParNames("A","B");
    TFitResultPtr fit=gr.Fit(&f,"QRS0");
    fr.fitted=true; fr.status=(int)fit; fr.covstatus=fit->CovMatrixStatus(); fr.edm=fit->Edm();
    fr.B=f.GetParameter(1); fr.Berr=f.GetParError(1); fr.chi2=f.GetChisquare(); fr.ndf=f.GetNDF();
    fr.chi2ndf=(fr.ndf>0)?fr.chi2/fr.ndf:0.0;
    fr.converged=(fr.status==0 && finite_number(fr.B));
    return fr;
}

// IMPORTANT: all three radii are read here.  For every hotspot the three fits
// use exactly the SAME set of genuinely-detected overvoltage points.  This
// prevents a missing/different operating point at one radius from becoming a
// fake systematic shift of B.
void lum_vs_v_systematics(const char* csv_R20,
                          const char* csv_R16,
                          const char* csv_R24,
                          const char* phase,
                          const char* output_csv) {
    string ph=phase;
    CsvTable t20=read_analysis_csv(csv_R20,true);
    CsvTable t16=read_analysis_csv(csv_R16,true);
    CsvTable t24=read_analysis_csv(csv_R24,true);
    auto a20=filter_phase_detected(t20.rows,ph);
    auto a16=filter_phase_detected(t16.rows,ph);
    auto a24=filter_phase_detected(t24.rows,ph);

    map<int,vector<AnalysisRow>> m20,m16,m24;
    for(const auto&r:a20)m20[r.spot].push_back(r);
    for(const auto&r:a16)m16[r.spot].push_back(r);
    for(const auto&r:a24)m24[r.spot].push_back(r);

    std::ofstream fout(output_csv);
    if(!fout.is_open()){std::cerr<<"Error: cannot create "<<output_csv<<std::endl;return;}
    fout<<"spot,phase,n_common_points,B_R20,Berr_R20,fit_status_R20,B_R16,Berr_R16,fit_status_R16,B_R24,Berr_R24,fit_status_R24,deltaB,converged_all\n";

    for(auto&kv:m20){
        int id=kv.first;
        if(!m16.count(id)||!m24.count(id)) continue;
        vector<AnalysisRow> r20,r16,r24;
        auto base=kv.second;
        std::sort(base.begin(),base.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.v_fin<b.v_fin;});
        for(const auto&r:base){
            const AnalysisRow* p16=find_v_point(m16[id],r.v_fin);
            const AnalysisRow* p24=find_v_point(m24[id],r.v_fin);
            if(!p16||!p24) continue;
            r20.push_back(r); r16.push_back(*p16); r24.push_back(*p24);
        }
        if(ppoint.is_open()) {
            for(size_t ip=0;ip<r20.size();++ip){
                double dL=std::max(std::fabs(r24[ip].luminosity-r20[ip].luminosity),
                                   std::fabs(r20[ip].luminosity-r16[ip].luminosity));
                ppoint<<id<<','<<ph<<','<<r20[ip].T<<','<<r20[ip].v<<','<<r20[ip].v_fin<<','
                      <<r16[ip].luminosity<<','<<r20[ip].luminosity<<','<<r24[ip].luminosity<<','<<dL<<'\n';
            }
        }
        if(r20.size()<2) continue;
        auto f20=fit_v_rows(r20,"R20",id);
        auto f16=fit_v_rows(r16,"R16",id);
        auto f24=fit_v_rows(r24,"R24",id);
        bool ok=f20.converged&&f16.converged&&f24.converged;
        double dB=ok?std::max(std::fabs(f24.B-f20.B),std::fabs(f20.B-f16.B)):
                      std::numeric_limits<double>::quiet_NaN();
        fout<<id<<','<<ph<<','<<r20.size()<<','
            <<f20.B<<','<<f20.Berr<<','<<f20.status<<','
            <<f16.B<<','<<f16.Berr<<','<<f16.status<<','
            <<f24.B<<','<<f24.Berr<<','<<f24.status<<','
            <<dB<<','<<(ok?1:0)<<'\n';
    }
    std::cout<<"Saved common-point power-law systematic fits to "<<output_csv<<std::endl;
}
