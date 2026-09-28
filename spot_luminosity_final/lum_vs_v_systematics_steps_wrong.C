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

// Weighted linear regression ln(L)=ln(A)+B ln(V).
static void estimate_power_parameters_sys(const vector<AnalysisRow>& rows,double& A0,double& B0) {
    double S=0.0,Sx=0.0,Sy=0.0,Sxx=0.0,Sxy=0.0;
    int nused=0;
    for (const auto& r:rows) {
        if (!(r.v_fin>0.0) || !(r.luminosity>0.0) || !finite_number(r.v_fin) || !finite_number(r.luminosity)) continue;
        const double xx=std::log(r.v_fin);
        const double yy=std::log(r.luminosity);
        double w=1.0;
        if (finite_number(r.error) && r.error>0.0) {
            const double sigma_log=r.error/r.luminosity;
            if (sigma_log>0.0 && finite_number(sigma_log)) w=1.0/(sigma_log*sigma_log);
        }
        S+=w;Sx+=w*xx;Sy+=w*yy;Sxx+=w*xx*xx;Sxy+=w*xx*yy;++nused;
    }
    const double D=S*Sxx-Sx*Sx;
    if (nused>=2 && S>0.0 && std::fabs(D)>1e-20) {
        B0=(S*Sxy-Sx*Sy)/D;
        const double intercept=(Sy-B0*Sx)/S;
        A0=std::exp(intercept);
    } else { A0=1.0; B0=2.0; }
    if (!(A0>0.0) || !finite_number(A0)) A0=1.0;
    if (!finite_number(B0)) B0=2.0;
}

// Extract the two central distinct overvoltage values directly from the data.
static void central_vrefs_sys(const vector<AnalysisRow>& rows,double& vref1,double& vref2) {
    vector<double> values;
    for (const auto& r:rows) if (r.v_fin>0.0 && finite_number(r.v_fin)) values.push_back(r.v_fin);
    std::sort(values.begin(),values.end());
    vector<double> u;
    for (double v:values) if (u.empty() || std::fabs(v-u.back())>1e-9) u.push_back(v);
    if (u.empty()) { vref1=1.0; vref2=1.0; return; }
    if (u.size()==1) { vref1=vref2=u[0]; return; }
    if (u.size()%2==0) { vref1=u[u.size()/2-1]; vref2=u[u.size()/2]; }
    else {
        size_t m=u.size()/2;
        vref1=u[m];
        vref2=(m+1<u.size()) ? u[m+1] : u[m-1];
    }
}

static VFitResult fit_v_rows(const vector<AnalysisRow>& rows,const string& tag,int spot) {
    VFitResult fr; fr.spot=spot;
    if (rows.size()<2) return fr;
    vector<double>x,y,ex,ey;
    for(const auto&r:rows){x.push_back(r.v_fin);y.push_back(r.luminosity);ex.push_back(0.0);ey.push_back(r.error);}
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());

    double A0,B0; estimate_power_parameters_sys(rows,A0,B0);
    double vref1,vref2; central_vrefs_sys(rows,vref1,vref2);

    // Step 1: centered nonlinear fit at the first central overvoltage.
    double Aref1=A0*std::pow(vref1,B0);
    if (!(Aref1>0.0) || !finite_number(Aref1)) Aref1=1.0;
    TF1 f1(Form("f_sys_v1_%s_%d",tag.c_str(),spot),
           Form("[0]*TMath::Power(x/%.17g,[1])",vref1),0.0,8.0);
    f1.SetParNames("A_ref","B"); f1.SetParameters(Aref1,B0); f1.SetParLimits(0,1e-300,1e300);
    TFitResultPtr fit1=gr.Fit(&f1,"QRSN");
    double B1=((int)fit1==0 && finite_number(f1.GetParameter(1))) ? f1.GetParameter(1) : B0;
    double A1=((int)fit1==0 && f1.GetParameter(0)>0.0 && finite_number(f1.GetParameter(0))) ? f1.GetParameter(0) : Aref1;

    // Step 2: exact reparameterization at the second central overvoltage.
    double Aref2=A1*std::pow(vref2/vref1,B1);
    if (!(Aref2>0.0) || !finite_number(Aref2)) Aref2=Aref1;
    TF1 f2(Form("f_sys_v2_%s_%d",tag.c_str(),spot),
           Form("[0]*TMath::Power(x/%.17g,[1])",vref2),0.0,8.0);
    f2.SetParNames("A_ref","B"); f2.SetParameters(Aref2,B1); f2.SetParLimits(0,1e-300,1e300);
    TFitResultPtr fit=gr.Fit(&f2,"QRSN");

    fr.fitted=true; fr.status=(int)fit; fr.covstatus=fit->CovMatrixStatus(); fr.edm=fit->Edm();
    fr.B=f2.GetParameter(1); fr.Berr=f2.GetParError(1); fr.chi2=f2.GetChisquare(); fr.ndf=f2.GetNDF();
    fr.chi2ndf=(fr.ndf>0)?fr.chi2/fr.ndf:0.0;
    fr.converged=(fr.status==0 && finite_number(fr.B));
    return fr;
}

// IMPORTANT: all three radii are read here. For every hotspot the three fits
// use exactly the SAME set of genuinely-detected overvoltage points.
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
