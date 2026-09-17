#include "analysis_common.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

struct TFitResultSys {
    int spot=-1;
    double x=0.0,y=0.0;
    double A=0.0,Aerr=0.0,lambda=0.0,lambdaerr=0.0;
    double chi2=0.0,chi2ndf=0.0,prob=0.0,edm=-1.0;
    int ndf=0,status=-999,covstatus=-999;
    bool fitted=false,converged=false;
};

static void estimate_exp_parameters_sys(const vector<AnalysisRow>& rows,double& A0,double& lambda0) {
    vector<AnalysisRow> pos;
    for (const auto& r: rows) if (r.luminosity>0) pos.push_back(r);
    std::sort(pos.begin(),pos.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
    if (pos.size()>=2 && std::fabs(pos.back().T-pos.front().T)>1e-12) {
        lambda0=(std::log(pos.back().luminosity)-std::log(pos.front().luminosity))/(pos.back().T-pos.front().T);
        A0=std::exp(std::log(pos.front().luminosity)-lambda0*pos.front().T);
    } else if (pos.size()==1) { A0=pos[0].luminosity; lambda0=0.0; }
    else { A0=1.0; lambda0=0.0; }
}

static map<int,TFitResultSys> fit_T_file(const string& filename,const string& phase) {
    CsvTable t=read_analysis_csv(filename,false);
    auto rows=filter_phase(t.rows,phase);
    map<int,vector<AnalysisRow>> by_spot;
    for (auto&r:rows) by_spot[r.spot].push_back(r);

    map<int,TFitResultSys> out;
    for (auto&kv:by_spot) {
        int spot=kv.first; auto rr=kv.second;
        std::sort(rr.begin(),rr.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});
        TFitResultSys fr; fr.spot=spot;
        if(!rr.empty()){fr.x=rr[0].x;fr.y=rr[0].y;}
        if(rr.size()<3){out[spot]=fr;continue;}

        double Tmin=rr.front().T,Tmax=rr.back().T,dT=Tmax-Tmin; if(dT<=0)dT=1.0;
        double fitmin=Tmin-0.10*dT,fitmax=Tmax+0.10*dT;
        vector<double>xv,yv,exv,eyv;
        for(auto&r:rr){xv.push_back(r.T);yv.push_back(r.luminosity);exv.push_back(0.0);eyv.push_back(r.error);}
        TGraphErrors gr((int)xv.size(),xv.data(),yv.data(),exv.data(),eyv.data());
        double A0,lambda0; estimate_exp_parameters_sys(rr,A0,lambda0);
        TF1 f(Form("f_sys_T_%d_%p",spot,(void*)&gr),"[0]*exp([1]*x)",fitmin,fitmax);
        f.SetParNames("A","lambda"); f.SetParameters(A0,lambda0);
        TFitResultPtr fit=gr.Fit(&f,"QRS0");
        fr.fitted=true; fr.status=(int)fit; fr.covstatus=fit->CovMatrixStatus(); fr.edm=fit->Edm();
        fr.A=f.GetParameter(0); fr.Aerr=f.GetParError(0); fr.lambda=f.GetParameter(1); fr.lambdaerr=f.GetParError(1);
        fr.chi2=f.GetChisquare();fr.ndf=f.GetNDF();fr.chi2ndf=(fr.ndf>0)?fr.chi2/fr.ndf:0.0;fr.prob=f.GetProb();
        fr.converged=(fr.status==0 && finite_number(fr.lambda));
        out[spot]=fr;
    }
    return out;
}

// Fits R=16 and R=24 with the same exponential model used by the nominal T analysis.
void lum_vs_T_systematics(const char* csv_R16,
                          const char* csv_R24,
                          const char* phase,
                          const char* output_csv) {
    string ph=phase;
    auto r16=fit_T_file(csv_R16,ph), r24=fit_T_file(csv_R24,ph);
    std::set<int> ids; for(auto&kv:r16)ids.insert(kv.first);for(auto&kv:r24)ids.insert(kv.first);
    std::ofstream fout(output_csv);
    if(!fout.is_open()){std::cerr<<"Error: cannot create "<<output_csv<<std::endl;return;}
    fout<<"spot,phase,lambda_R16,lambdaerr_R16,chi2_R16,ndf_R16,chi2ndf_R16,fit_status_R16,covmatrix_status_R16,edm_R16,converged_R16,"
          "lambda_R24,lambdaerr_R24,chi2_R24,ndf_R24,chi2ndf_R24,fit_status_R24,covmatrix_status_R24,edm_R24,converged_R24\n";
    const double NaN=std::numeric_limits<double>::quiet_NaN();
    for(int id:ids){
        TFitResultSys a,b; bool ha=r16.count(id),hb=r24.count(id);
        if(ha)a=r16[id];else{a.spot=id;a.lambda=a.lambdaerr=a.chi2=a.chi2ndf=a.edm=NaN;}
        if(hb)b=r24[id];else{b.spot=id;b.lambda=b.lambdaerr=b.chi2=b.chi2ndf=b.edm=NaN;}
        fout<<id<<','<<ph<<','
            <<a.lambda<<','<<a.lambdaerr<<','<<a.chi2<<','<<a.ndf<<','<<a.chi2ndf<<','<<a.status<<','<<a.covstatus<<','<<a.edm<<','<<(a.converged?1:0)<<','
            <<b.lambda<<','<<b.lambdaerr<<','<<b.chi2<<','<<b.ndf<<','<<b.chi2ndf<<','<<b.status<<','<<b.covstatus<<','<<b.edm<<','<<(b.converged?1:0)<<'\n';
    }
    std::cout<<"Saved exponential systematic-fit values to "<<output_csv<<std::endl;
}
