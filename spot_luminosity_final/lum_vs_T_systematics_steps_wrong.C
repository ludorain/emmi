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

// Weighted linear regression of ln(L) versus T-Tref.
// For sigma_L>0 the weight is 1/sigma_lnL^2=(L/sigma_L)^2.
static void estimate_exp_parameters_sys(const vector<AnalysisRow>& rows,
                                        double Tref,
                                        double& Aref0,
                                        double& lambda0) {
    double S=0.0,Sx=0.0,Sy=0.0,Sxx=0.0,Sxy=0.0;
    int nused=0;
    for (const auto& r:rows) {
        if (!(r.luminosity>0.0) || !finite_number(r.luminosity) || !finite_number(r.T)) continue;
        const double xx=r.T-Tref;
        const double yy=std::log(r.luminosity);
        double w=1.0;
        if (finite_number(r.error) && r.error>0.0) {
            const double sigma_log=r.error/r.luminosity;
            if (sigma_log>0.0 && finite_number(sigma_log)) w=1.0/(sigma_log*sigma_log);
        }
        S+=w; Sx+=w*xx; Sy+=w*yy; Sxx+=w*xx*xx; Sxy+=w*xx*yy;
        ++nused;
    }
    const double D=S*Sxx-Sx*Sx;
    if (nused>=2 && S>0.0 && std::fabs(D)>1e-20) {
        lambda0=(S*Sxy-Sx*Sy)/D;
        const double intercept=(Sy-lambda0*Sx)/S;
        Aref0=std::exp(intercept);
    } else {
        Aref0=1.0;
        lambda0=0.0;
    }
    if (!(Aref0>0.0) || !finite_number(Aref0)) Aref0=1.0;
    if (!finite_number(lambda0)) lambda0=0.0;
}

static TFitResultSys fit_T_rows(const vector<AnalysisRow>& rows,const string& tag,int spot){
    TFitResultSys fr; fr.spot=spot; if(rows.size()<3)return fr;
    auto rr=rows; std::sort(rr.begin(),rr.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
    double Tmin=rr.front().T,Tmax=rr.back().T,dT=Tmax-Tmin;if(dT<=0)dT=1.0;
    vector<double>x,y,ex,ey;for(const auto&r:rr){x.push_back(r.T);y.push_back(r.luminosity);ex.push_back(0.0);ey.push_back(r.error);}
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());

    // Step 1: centered fit at 19 C, initialized by the weighted log-linear prefit.
    double A19,l0; estimate_exp_parameters_sys(rr,19.0,A19,l0);
    TF1 f19(Form("f_sys_T19_%s_%d",tag.c_str(),spot),"[0]*exp([1]*(x-19.0))",Tmin-.10*dT,Tmax+.10*dT);
    f19.SetParNames("A_{ref,19}","lambda");
    f19.SetParameters(A19,l0);
    f19.SetParLimits(0,1e-300,1e300);
    TFitResultPtr fit19=gr.Fit(&f19,"QRSN");

    double l1=((int)fit19==0 && finite_number(f19.GetParameter(1))) ? f19.GetParameter(1) : l0;
    double A19fit=((int)fit19==0 && f19.GetParameter(0)>0.0 && finite_number(f19.GetParameter(0))) ? f19.GetParameter(0) : A19;

    // Step 2: same model centered at 21 C.  The previous fit is exactly
    // reparameterized before minimization: A21=A19*exp(lambda*(21-19)).
    double A21=A19fit*std::exp(l1*(21.0-19.0));
    if (!(A21>0.0) || !finite_number(A21)) A21=A19;
    TF1 f21(Form("f_sys_T21_%s_%d",tag.c_str(),spot),"[0]*exp([1]*(x-21.0))",Tmin-.10*dT,Tmax+.10*dT);
    f21.SetParNames("A_{ref,21}","lambda");
    f21.SetParameters(A21,l1);
    f21.SetParLimits(0,1e-300,1e300);
    TFitResultPtr fit=gr.Fit(&f21,"QRSN");

    fr.fitted=true;fr.status=(int)fit;fr.covstatus=fit->CovMatrixStatus();fr.edm=fit->Edm();
    fr.lambda=f21.GetParameter(1);fr.lambdaerr=f21.GetParError(1);fr.chi2=f21.GetChisquare();fr.ndf=f21.GetNDF();
    fr.chi2ndf=(fr.ndf>0)?fr.chi2/fr.ndf:0.0;
    fr.converged=(fr.status==0&&finite_number(fr.lambda));return fr;
}

static const AnalysisRow* find_T_point(const vector<AnalysisRow>& rows,double x,double tol=1e-9){
    for(const auto&r:rows)if(std::fabs(r.T-x)<=tol)return &r;
    return nullptr;
}

// R16/R20/R24 are fitted on the exact same genuinely-detected temperature
// points. This is essential for lambda: comparing independent fits with a
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
