#include "analysis_common.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

struct VFitResult {
    int spot=-1;
    double x=0.0,y=0.0;
    double A=0.0,Aerr=0.0,B=0.0,Berr=0.0;
    double chi2=0.0,chi2ndf=0.0,prob=0.0,edm=-1.0;
    int ndf=0,status=-999,covstatus=-999;
    bool fitted=false,converged=false;
};

static map<int,VFitResult> fit_v_file(const string& filename, const string& phase) {
    CsvTable t = read_analysis_csv(filename, true);
    vector<AnalysisRow> rows = filter_phase(t.rows, phase);
    map<int, vector<AnalysisRow>> by_spot;
    for (const auto& r : rows) by_spot[r.spot].push_back(r);

    map<int,VFitResult> out;
    for (auto& kv : by_spot) {
        int spot = kv.first;
        auto rr = kv.second;
        std::sort(rr.begin(), rr.end(), [](const AnalysisRow& a,const AnalysisRow& b){return a.v_fin<b.v_fin;});
        VFitResult fr; fr.spot=spot;
        if (!rr.empty()) { fr.x=rr[0].x; fr.y=rr[0].y; }
        if (rr.size()<2) { out[spot]=fr; continue; }

        vector<double> xv, yv, exv, eyv;
        for (const auto& r: rr) {
            xv.push_back(r.v_fin); yv.push_back(r.luminosity);
            exv.push_back(0.0); eyv.push_back(r.error);
        }

        TGraphErrors gr((int)xv.size(),xv.data(),yv.data(),exv.data(),eyv.data());
        TF1 f(Form("f_sys_v_%d_%p",spot,(void*)&gr), "[0]*TMath::Power(x,[1])", 0.0, 8.0);
        f.SetParameters(1.0,2.0);
        f.SetParNames("A","B");

        TFitResultPtr fit = gr.Fit(&f,"QRS0");
        fr.fitted=true;
        fr.status=(int)fit;
        fr.covstatus=fit->CovMatrixStatus();
        fr.edm=fit->Edm();
        fr.A=f.GetParameter(0); fr.Aerr=f.GetParError(0);
        fr.B=f.GetParameter(1); fr.Berr=f.GetParError(1);
        fr.chi2=f.GetChisquare(); fr.ndf=f.GetNDF();
        fr.chi2ndf=(fr.ndf>0)?fr.chi2/fr.ndf:0.0;
        fr.prob=f.GetProb();
        fr.converged=(fr.status==0 && finite_number(fr.B));
        out[spot]=fr;
    }
    return out;
}

// Fits the R=16 and R=24 datasets with exactly the same power-law model
// used by the nominal overvoltage analysis. The output is recreated at each call.
void lum_vs_v_systematics(const char* csv_R16,
                          const char* csv_R24,
                          const char* phase,
                          const char* output_csv) {
    string ph = phase;
    auto r16 = fit_v_file(csv_R16, ph);
    auto r24 = fit_v_file(csv_R24, ph);

    std::set<int> ids;
    for (auto& kv:r16) ids.insert(kv.first);
    for (auto& kv:r24) ids.insert(kv.first);

    std::ofstream fout(output_csv);
    if (!fout.is_open()) { std::cerr << "Error: cannot create " << output_csv << std::endl; return; }
    fout << "spot,phase,B_R16,Berr_R16,chi2_R16,ndf_R16,chi2ndf_R16,fit_status_R16,covmatrix_status_R16,edm_R16,converged_R16,"
            "B_R24,Berr_R24,chi2_R24,ndf_R24,chi2ndf_R24,fit_status_R24,covmatrix_status_R24,edm_R24,converged_R24\n";

    const double NaN = std::numeric_limits<double>::quiet_NaN();
    for (int id: ids) {
        VFitResult a,b;
        bool ha=r16.count(id), hb=r24.count(id);
        if (ha) a=r16[id]; else {a.spot=id;a.B=a.Berr=a.chi2=a.chi2ndf=a.edm=NaN;}
        if (hb) b=r24[id]; else {b.spot=id;b.B=b.Berr=b.chi2=b.chi2ndf=b.edm=NaN;}
        fout << id << ',' << ph << ','
             << a.B << ',' << a.Berr << ',' << a.chi2 << ',' << a.ndf << ',' << a.chi2ndf << ',' << a.status << ',' << a.covstatus << ',' << a.edm << ',' << (a.converged?1:0) << ','
             << b.B << ',' << b.Berr << ',' << b.chi2 << ',' << b.ndf << ',' << b.chi2ndf << ',' << b.status << ',' << b.covstatus << ',' << b.edm << ',' << (b.converged?1:0) << '\n';
    }
    std::cout << "Saved power-law systematic-fit values to " << output_csv << std::endl;
}
