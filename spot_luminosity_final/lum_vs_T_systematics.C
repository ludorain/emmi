#include "analysis_common.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

struct TFitResultSys {
    int spot=-1;
    double lambda=std::numeric_limits<double>::quiet_NaN();
    double lambdaerr=std::numeric_limits<double>::quiet_NaN();
    double chi2=std::numeric_limits<double>::quiet_NaN();
    double chi2ndf=std::numeric_limits<double>::quiet_NaN();
    double edm=std::numeric_limits<double>::quiet_NaN();
    int ndf=0,status=-999,covstatus=-999;
    bool fitted=false,converged=false;
};

static void estimate_exp_parameters_sys(const vector<AnalysisRow>& rows,double& A0,double& lambda0){
    vector<AnalysisRow> pos; for(const auto&r:rows)if(r.luminosity>0)pos.push_back(r);
    std::sort(pos.begin(),pos.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
    if(pos.size()>=2&&std::fabs(pos.back().T-pos.front().T)>1e-12){
        lambda0=(std::log(pos.back().luminosity)-std::log(pos.front().luminosity))/(pos.back().T-pos.front().T);
        A0=std::exp(std::log(pos.front().luminosity)-lambda0*pos.front().T);
    }else if(pos.size()==1){A0=pos[0].luminosity;lambda0=0.0;}else{A0=1.0;lambda0=0.0;}
}

static const AnalysisRow* find_T_point(const vector<AnalysisRow>& rows,double x,double tol=1e-9){
    for(const auto&r:rows)if(std::fabs(r.T-x)<=tol)return &r;
    return nullptr;
}

static TFitResultSys fit_T_rows(const vector<AnalysisRow>& rows,const string& tag,int spot){
    TFitResultSys fr; fr.spot=spot; if(rows.size()<3)return fr;
    auto rr=rows; std::sort(rr.begin(),rr.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
    double Tmin=rr.front().T,Tmax=rr.back().T,dT=Tmax-Tmin;if(dT<=0)dT=1.0;
    vector<double>x,y,ex,ey;for(const auto&r:rr){x.push_back(r.T);y.push_back(r.luminosity);ex.push_back(0.0);ey.push_back(r.error);}
    double A0,l0;estimate_exp_parameters_sys(rr,A0,l0);
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
    TF1 f(Form("f_sys_T_%s_%d",tag.c_str(),spot),"[0]*exp([1]*x)",Tmin-.10*dT,Tmax+.10*dT);
    f.SetParNames("A","lambda");f.SetParameters(A0,l0);
    TFitResultPtr fit=gr.Fit(&f,"QRS0");
    fr.fitted=true;fr.status=(int)fit;fr.covstatus=fit->CovMatrixStatus();fr.edm=fit->Edm();
    fr.lambda=f.GetParameter(1);fr.lambdaerr=f.GetParError(1);fr.chi2=f.GetChisquare();fr.ndf=f.GetNDF();
    fr.chi2ndf=(fr.ndf>0)?fr.chi2/fr.ndf:0.0;
    fr.converged=(fr.status==0&&finite_number(fr.lambda));return fr;
}

// R16/R20/R24 are fitted on the exact same genuinely-detected temperature
// points.  This is essential for lambda: comparing independent fits with a
// different T coverage can create a very large artificial systematic shift.
void lum_vs_T_systematics(const char* csv_R20,
                          const char* csv_R16,
                          const char* csv_R24,
                          const char* phase,
                          const char* output_csv){
    string ph=phase;
    CsvTable t20=read_analysis_csv(csv_R20,false),t16=read_analysis_csv(csv_R16,false),t24=read_analysis_csv(csv_R24,false);
    auto a20=filter_phase_detected(t20.rows,ph),a16=filter_phase_detected(t16.rows,ph),a24=filter_phase_detected(t24.rows,ph);
    map<int,vector<AnalysisRow>>m20,m16,m24;for(const auto&r:a20)m20[r.spot].push_back(r);for(const auto&r:a16)m16[r.spot].push_back(r);for(const auto&r:a24)m24[r.spot].push_back(r);
    std::ofstream fout(output_csv);if(!fout.is_open()){std::cerr<<"Error: cannot create "<<output_csv<<std::endl;return;}
    fout<<"spot,phase,n_common_points,lambda_R20,lambdaerr_R20,fit_status_R20,lambda_R16,lambdaerr_R16,fit_status_R16,lambda_R24,lambdaerr_R24,fit_status_R24,deltaLambda,converged_all\n";
    for(auto&kv:m20){
        int id=kv.first;if(!m16.count(id)||!m24.count(id))continue;
        vector<AnalysisRow>r20,r16,r24;auto base=kv.second;std::sort(base.begin(),base.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
        for(const auto&r:base){const AnalysisRow*p16=find_T_point(m16[id],r.T);const AnalysisRow*p24=find_T_point(m24[id],r.T);if(!p16||!p24)continue;r20.push_back(r);r16.push_back(*p16);r24.push_back(*p24);}
        if(r20.size()<3)continue;
        auto f20=fit_T_rows(r20,"R20",id),f16=fit_T_rows(r16,"R16",id),f24=fit_T_rows(r24,"R24",id);
        bool ok=f20.converged&&f16.converged&&f24.converged;
        double dl=ok?std::max(std::fabs(f24.lambda-f20.lambda),std::fabs(f20.lambda-f16.lambda)):std::numeric_limits<double>::quiet_NaN();
        fout<<id<<','<<ph<<','<<r20.size()<<','
            <<f20.lambda<<','<<f20.lambdaerr<<','<<f20.status<<','
            <<f16.lambda<<','<<f16.lambdaerr<<','<<f16.status<<','
            <<f24.lambda<<','<<f24.lambdaerr<<','<<f24.status<<','
            <<dl<<','<<(ok?1:0)<<'\n';
    }
    std::cout<<"Saved common-point exponential systematic fits to "<<output_csv<<std::endl;
}
