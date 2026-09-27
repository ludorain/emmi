#include "gold_analysis_common.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

struct GoldTFitResult {
    double lambda=std::numeric_limits<double>::quiet_NaN();
    double lambdaerr=std::numeric_limits<double>::quiet_NaN();
    int status=-999;
    bool converged=false;
};

static GoldTFitResult gold_fit_T_rows(const vector<GoldRow>& rows,const string& tag,int spot) {
    GoldTFitResult fr;
    if (rows.size()<3) return fr;
    auto rr=rows;
    std::sort(rr.begin(),rr.end(),[](const GoldRow&a,const GoldRow&b){return a.T<b.T;});
    const double Tmin=rr.front().T, Tmax=rr.back().T;
    double dT=Tmax-Tmin; if (!(dT>0.0)) dT=1.0;
    vector<double> x,y,ex,ey;
    for (const auto& r:rr) { x.push_back(r.T); y.push_back(r.luminosity); ex.push_back(0.0); ey.push_back(r.error); }
    double A0,l0; gold_estimate_exp(rr,A0,l0);
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
    TF1 f(Form("gold_sys_T_%s_%d",tag.c_str(),spot),"[0]*exp([1]*x)",Tmin-.10*dT,Tmax+.10*dT);
    f.SetParameters(A0,l0); f.SetParNames("A","#lambda");
    TFitResultPtr fit=gr.Fit(&f,"QRS0");
    fr.status=(int)fit;
    fr.lambda=f.GetParameter(1);
    fr.lambdaerr=f.GetParError(1);
    fr.converged=(fr.status==0 && gold_finite(fr.lambda));
    return fr;
}

void gold_T_systematics(const char* csv_R20,
                        const char* csv_R16,
                        const char* csv_R24,
                        const char* dataset_key,
                        const char* output_csv) {
    const string key=dataset_key;
    GoldTable t20=read_gold_csv(csv_R20,false);
    GoldTable t16=read_gold_csv(csv_R16,false);
    GoldTable t24=read_gold_csv(csv_R24,false);

    auto m20=gold_group_spots(gold_select_dataset(t20.rows,key,true));
    auto m16=gold_group_spots(gold_select_dataset(t16.rows,key,true));
    auto m24=gold_group_spots(gold_select_dataset(t24.rows,key,true));

    std::ofstream out(output_csv);
    if (!out.is_open()) { std::cerr << "Error: cannot create " << output_csv << std::endl; return; }
    out << "spot,dataset_key,n_common_points,lambda_R20,lambdaerr_R20,fit_status_R20,"
           "lambda_R16,lambdaerr_R16,fit_status_R16,lambda_R24,lambdaerr_R24,fit_status_R24,"
           "deltaLambda,converged_all\n";

    for (auto& kv:m20) {
        const int spot=kv.first;
        if (!m16.count(spot) || !m24.count(spot)) continue;
        auto base=kv.second;
        std::sort(base.begin(),base.end(),[](const GoldRow&a,const GoldRow&b){return a.T<b.T;});
        vector<GoldRow> r20,r16,r24;
        for (const auto& r:base) {
            const GoldRow* p16=gold_find_T(m16[spot],r.T);
            const GoldRow* p24=gold_find_T(m24[spot],r.T);
            if (!p16 || !p24) continue;
            r20.push_back(r); r16.push_back(*p16); r24.push_back(*p24);
        }
        if (r20.size()<3) continue;

        auto f20=gold_fit_T_rows(r20,"R20",spot);
        auto f16=gold_fit_T_rows(r16,"R16",spot);
        auto f24=gold_fit_T_rows(r24,"R24",spot);
        const bool ok=f20.converged&&f16.converged&&f24.converged;
        const double dl=ok ? std::max(std::fabs(f24.lambda-f20.lambda),std::fabs(f20.lambda-f16.lambda))
                           : std::numeric_limits<double>::quiet_NaN();
        out << spot << ',' << key << ',' << r20.size() << ','
            << f20.lambda << ',' << f20.lambdaerr << ',' << f20.status << ','
            << f16.lambda << ',' << f16.lambdaerr << ',' << f16.status << ','
            << f24.lambda << ',' << f24.lambdaerr << ',' << f24.status << ','
            << dl << ',' << (ok?1:0) << '\n';
    }
    std::cout << "Saved lambda systematic fits for " << key << " to " << output_csv << std::endl;
}
