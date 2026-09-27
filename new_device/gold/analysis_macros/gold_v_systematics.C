#include "gold_analysis_common.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

struct GoldVFitResult {
    double B=std::numeric_limits<double>::quiet_NaN();
    double Berr=std::numeric_limits<double>::quiet_NaN();
    int status=-999;
    bool converged=false;
};

static GoldVFitResult gold_fit_v_rows(const vector<GoldRow>& rows,const string& tag,int spot) {
    GoldVFitResult fr;
    if (rows.size()<2) return fr;
    vector<double> x,y,ex,ey;
    for (const auto& r:rows) {
        if (!gold_finite(r.v_fin)) continue;
        x.push_back(r.v_fin); y.push_back(r.luminosity); ex.push_back(0.0); ey.push_back(r.error);
    }
    if (x.size()<2) return fr;
    TGraphErrors gr((int)x.size(),x.data(),y.data(),ex.data(),ey.data());
    TF1 f(Form("gold_sys_v_%s_%d",tag.c_str(),spot),"[0]*TMath::Power(x,[1])",0.0,8.0);
    f.SetParameters(1.0,2.0); f.SetParNames("A","B");
    TFitResultPtr fit=gr.Fit(&f,"QRS0");
    fr.status=(int)fit;
    fr.B=f.GetParameter(1);
    fr.Berr=f.GetParError(1);
    fr.converged=(fr.status==0 && gold_finite(fr.B));
    return fr;
}

void gold_v_systematics(const char* csv_R20,
                        const char* csv_R16,
                        const char* csv_R24,
                        const char* output_csv) {
    const string key="T20";
    GoldTable t20=read_gold_csv(csv_R20,false);
    GoldTable t16=read_gold_csv(csv_R16,false);
    GoldTable t24=read_gold_csv(csv_R24,false);

    auto m20=gold_group_spots(gold_select_dataset(t20.rows,key,true));
    auto m16=gold_group_spots(gold_select_dataset(t16.rows,key,true));
    auto m24=gold_group_spots(gold_select_dataset(t24.rows,key,true));

    std::ofstream out(output_csv);
    if (!out.is_open()) { std::cerr << "Error: cannot create " << output_csv << std::endl; return; }
    out << "spot,n_common_points,B_R20,Berr_R20,fit_status_R20,"
           "B_R16,Berr_R16,fit_status_R16,B_R24,Berr_R24,fit_status_R24,deltaB,converged_all\n";

    for (auto& kv:m20) {
        const int spot=kv.first;
        if (!m16.count(spot) || !m24.count(spot)) continue;
        auto base=kv.second;
        std::sort(base.begin(),base.end(),[](const GoldRow&a,const GoldRow&b){return a.v_fin<b.v_fin;});
        vector<GoldRow> r20,r16,r24;
        for (const auto& r:base) {
            const GoldRow* p16=gold_find_vfin(m16[spot],r.v_fin);
            const GoldRow* p24=gold_find_vfin(m24[spot],r.v_fin);
            if (!p16 || !p24) continue;
            r20.push_back(r); r16.push_back(*p16); r24.push_back(*p24);
        }
        if (r20.size()<2) continue;

        auto f20=gold_fit_v_rows(r20,"R20",spot);
        auto f16=gold_fit_v_rows(r16,"R16",spot);
        auto f24=gold_fit_v_rows(r24,"R24",spot);
        const bool ok=f20.converged&&f16.converged&&f24.converged;
        const double dB=ok ? std::max(std::fabs(f24.B-f20.B),std::fabs(f20.B-f16.B))
                           : std::numeric_limits<double>::quiet_NaN();
        out << spot << ',' << r20.size() << ','
            << f20.B << ',' << f20.Berr << ',' << f20.status << ','
            << f16.B << ',' << f16.Berr << ',' << f16.status << ','
            << f24.B << ',' << f24.Berr << ',' << f24.status << ','
            << dB << ',' << (ok?1:0) << '\n';
    }
    std::cout << "Saved B systematic fits to " << output_csv << std::endl;
}
