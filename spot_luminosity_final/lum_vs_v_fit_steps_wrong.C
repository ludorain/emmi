#include "analysis_common.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TMath.h"

struct BSystematicPair {
    bool has16=false, has24=false;
    double B16=0.0, B24=0.0;
    bool conv16=false, conv24=false;
    bool has_delta=false;
    double deltaB=0.0;
};

struct VNominalFit {
    int spot=-1;
    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TF1* func=nullptr;

    double A=0.0, Aerr=0.0;
    double B=0.0, Berr=0.0;
    double chi2=0.0, chi2ndf=0.0, prob=0.0, edm=-1.0;
    int ndf=0, status=-999, covstatus=-999;

    bool fit_done=false, converged=false;
    bool has_param_syst=false;
    double deltaB=0.0;
};

struct VRadiusOverlayFit {
    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TF1* func=nullptr;
    double B=0.0, Berr=0.0;
    int status=-999;
    bool fit_done=false, converged=false;
};

// -----------------------------------------------------------------------------
// Read the R=16 and R=24 fit values produced by lum_vs_v_systematics.C.
// ROOT fit status == 0 is the convergence criterion; ndf is intentionally not
// required because a two-point/two-parameter fit can converge with ndf == 0.
// -----------------------------------------------------------------------------
static map<int,BSystematicPair> read_B_systematics(const string& filename) {
    map<int,BSystematicPair> out;
    std::ifstream fin(filename);
    if (!fin.is_open()) return out;

    string line;
    if (!std::getline(fin,line)) return out;
    auto header=split_csv_simple(line);
    map<string,int> col;
    for (int i=0;i<(int)header.size();++i) col[header[i]]=i;

    while (std::getline(fin,line)) {
        if (trim_copy(line).empty()) continue;
        auto f=split_csv_simple(line);
        if (f.size()<header.size()) f.resize(header.size(),"");

        try {
            int id=(int)std::llround(std::stod(f[col["spot"]]));
            BSystematicPair s;

            if (col.count("B_R16") && !f[col["B_R16"]].empty()) {
                s.B16=std::stod(f[col["B_R16"]]);
                s.has16=finite_number(s.B16);
            }
            if (col.count("B_R24") && !f[col["B_R24"]].empty()) {
                s.B24=std::stod(f[col["B_R24"]]);
                s.has24=finite_number(s.B24);
            }

            if (col.count("fit_status_R16") && !f[col["fit_status_R16"]].empty())
                s.conv16=(std::stoi(f[col["fit_status_R16"]])==0 && s.has16);
            else if (col.count("converged_R16"))
                s.conv16=parse_bool_safe(f[col["converged_R16"]]);

            if (col.count("fit_status_R24") && !f[col["fit_status_R24"]].empty())
                s.conv24=(std::stoi(f[col["fit_status_R24"]])==0 && s.has24);
            else if (col.count("converged_R24"))
                s.conv24=parse_bool_safe(f[col["converged_R24"]]);

            if (col.count("deltaB") && !f[col["deltaB"]].empty()) {
                s.deltaB=std::fabs(std::stod(f[col["deltaB"]]));
                s.has_delta=finite_number(s.deltaB);
            }
            out[id]=s;
        } catch (...) {
            // Skip malformed systematic rows without aborting the full analysis.
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
// Draw a symmetric systematic uncertainty as two green horizontal brackets.
// The bracket width is defined as a small fraction of the visible x range, so it
// remains only slightly wider than the experimental marker for any axis scale.
// This is deliberately different from a conventional statistical error bar.
// -----------------------------------------------------------------------------
static void draw_horizontal_syst_brackets(const vector<double>& x,
                                           const vector<double>& y,
                                           const vector<double>& delta,
                                           double xmin,
                                           double xmax,
                                           int color=kBlack,
                                           int line_width=4) {
    if (x.empty() || y.size()!=x.size() || delta.size()!=x.size()) return;

    double span=xmax-xmin;
    if (!(span>0.0)) span=1.0;

    // Total cap width = 0.9% of the visible x range.
    // On the 1800-pixel canvas this is slightly wider than MarkerSize(1.5).
    const double half_width=0.009*span;

    for (size_t i=0;i<x.size();++i) {
        if (!finite_number(delta[i]) || delta[i]<=0.0) continue;

        const double ytop=y[i]+delta[i];
        const double ybottom=y[i]-delta[i];
        const double hook=0.20*delta[i];

        // Upper horizontal bracket with short downward hooks.
        TLine* top=new TLine(x[i]-half_width,ytop,x[i]+half_width,ytop);
        top->SetLineColor(color); top->SetLineWidth(line_width); top->Draw("SAME");
        TLine* top_l=new TLine(x[i]-half_width,ytop,x[i]-half_width,ytop-hook);
        top_l->SetLineColor(color); top_l->SetLineWidth(line_width); top_l->Draw("SAME");
        TLine* top_r=new TLine(x[i]+half_width,ytop,x[i]+half_width,ytop-hook);
        top_r->SetLineColor(color); top_r->SetLineWidth(line_width); top_r->Draw("SAME");

        // Lower horizontal bracket with short upward hooks.
        TLine* bottom=new TLine(x[i]-half_width,ybottom,x[i]+half_width,ybottom);
        bottom->SetLineColor(color); bottom->SetLineWidth(line_width); bottom->Draw("SAME");
        TLine* bottom_l=new TLine(x[i]-half_width,ybottom,x[i]-half_width,ybottom+hook);
        bottom_l->SetLineColor(color); bottom_l->SetLineWidth(line_width); bottom_l->Draw("SAME");
        TLine* bottom_r=new TLine(x[i]+half_width,ybottom,x[i]+half_width,ybottom+hook);
        bottom_r->SetLineColor(color); bottom_r->SetLineWidth(line_width); bottom_r->Draw("SAME");
    }
}

// Weighted linear prefit ln(L)=ln(A)+B ln(V), using all positive points.
static void estimate_power_parameters_nom(const vector<AnalysisRow>& rows,
                                          double A_fallback,double B_fallback,
                                          double& A0,double& B0) {
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
    } else { A0=A_fallback; B0=B_fallback; }
    if (!(A0>0.0) || !finite_number(A0)) A0=A_fallback;
    if (!finite_number(B0)) B0=B_fallback;
}

// The two reference overvoltages are read from the actual v_fin values.
static void central_vrefs_nom(const vector<AnalysisRow>& rows,double& vref1,double& vref2) {
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

// Two centered nonlinear fits.  The final function uses the second central
// overvoltage.  A_ref is positive; A0,B0 from the log-linear prefit are only
// starting values and the fit itself is performed in the original L(V) scale.
static TFitResultPtr fit_power_centered_nom(TGraphErrors* graph,
                                           const vector<AnalysisRow>& rows,
                                           const string& final_name,
                                           double fit_xmin,double fit_xmax,
                                           double A_fallback,double B_fallback,
                                           const char* final_options,
                                           TF1*& final_func,
                                           double& final_vref) {
    double A0,B0;
    estimate_power_parameters_nom(rows,A_fallback,B_fallback,A0,B0);
    double vref1,vref2;
    central_vrefs_nom(rows,vref1,vref2);

    double Aref1=A0*std::pow(vref1,B0);
    if (!(Aref1>0.0) || !finite_number(Aref1)) Aref1=A_fallback;
    TF1 f1((final_name+"_prefit1").c_str(),
           Form("[0]*TMath::Power(x/%.17g,[1])",vref1),fit_xmin,fit_xmax);
    f1.SetParNames("A_ref","B");
    f1.SetParameters(Aref1,B0);
    f1.SetParLimits(0,1e-300,1e300);
    TFitResultPtr fit1=graph->Fit(&f1,"QRSN");

    const double B1=((int)fit1==0 && finite_number(f1.GetParameter(1))) ? f1.GetParameter(1) : B0;
    const double A1=((int)fit1==0 && f1.GetParameter(0)>0.0 && finite_number(f1.GetParameter(0))) ? f1.GetParameter(0) : Aref1;
    double Aref2=A1*std::pow(vref2/vref1,B1);
    if (!(Aref2>0.0) || !finite_number(Aref2)) Aref2=Aref1;

    final_vref=vref2;
    final_func=new TF1(final_name.c_str(),
                       Form("[0]*TMath::Power(x/%.17g,[1])",vref2),fit_xmin,fit_xmax);
    final_func->SetParNames("A_ref","B");
    final_func->SetParameters(Aref2,B1);
    final_func->SetParLimits(0,1e-300,1e300);
    return graph->Fit(final_func,final_options);
}

// Convert the final centered amplitude back to the historical coefficient A
// in L=A*V^B, preserving the existing CSV/output meaning.
static void convert_power_Aref_to_A(const TFitResultPtr& fit,
                                    double Aref,double ArefErr,
                                    double B,double Berr,double vref,
                                    double& A,double& Aerr) {
    const double scale=std::pow(vref,-B);
    A=Aref*scale;
    double cov=0.0;
    if ((int)fit>=0) cov=fit->CovMatrix(0,1);
    const double dA_dAref=scale;
    const double dA_dB=-std::log(vref)*A;
    double var=dA_dAref*dA_dAref*ArefErr*ArefErr
              +dA_dB*dA_dB*Berr*Berr
              +2.0*dA_dAref*dA_dB*cov;
    if (var<0.0 && std::fabs(var)<1e-12*std::max(1.0,A*A)) var=0.0;
    Aerr=(var>=0.0 && finite_number(var)) ? std::sqrt(var) : std::fabs(scale*ArefErr);
}

struct BSimultaneousResult {
    bool ok=false;
    double B=std::numeric_limits<double>::quiet_NaN();
    double Berr=std::numeric_limits<double>::quiet_NaN();
    double chi2=std::numeric_limits<double>::quiet_NaN();
    double prob=0.0;
    double vref=1.0;
    int ndf=0,nspots=0,npoints=0;
};

static double profile_common_B_chi2(const map<int,vector<AnalysisRow>>& spots,
                                    const std::set<int>& accepted,
                                    double B,double vref,
                                    int* npoints_out=nullptr,
                                    int* nspots_out=nullptr) {
    double chi2=0.0;
    int npoints=0,nspots=0;
    for (int id:accepted) {
        auto it=spots.find(id);
        if (it==spots.end()) continue;
        const auto& rows=it->second;
        double num=0.0,den=0.0;
        int nv=0;
        for (const auto& r:rows) {
            if (!(r.v_fin>0.0) || !finite_number(r.v_fin) || !finite_number(r.luminosity)) continue;
            const double sigma=(finite_number(r.error) && r.error>0.0) ? r.error : 1.0;
            const double g=std::pow(r.v_fin/vref,B);
            if (!finite_number(g)) return 1e300;
            const double w=1.0/(sigma*sigma);
            num+=w*r.luminosity*g;
            den+=w*g*g;
            ++nv;
        }
        if (nv<2 || !(den>0.0)) continue;
        double Aref=num/den;
        if (!(Aref>0.0) || !finite_number(Aref)) Aref=1e-300;
        for (const auto& r:rows) {
            if (!(r.v_fin>0.0) || !finite_number(r.v_fin) || !finite_number(r.luminosity)) continue;
            const double sigma=(finite_number(r.error) && r.error>0.0) ? r.error : 1.0;
            const double model=Aref*std::pow(r.v_fin/vref,B);
            if (!finite_number(model)) return 1e300;
            const double pull=(r.luminosity-model)/sigma;
            chi2+=pull*pull;
            ++npoints;
        }
        ++nspots;
    }
    if (npoints_out) *npoints_out=npoints;
    if (nspots_out) *nspots_out=nspots;
    return chi2;
}

static double golden_min_B(const map<int,vector<AnalysisRow>>& spots,
                           const std::set<int>& accepted,double vref,
                           double lo,double hi) {
    const double gr=0.6180339887498948482;
    double c=hi-gr*(hi-lo),d=lo+gr*(hi-lo);
    double fc=profile_common_B_chi2(spots,accepted,c,vref),fd=profile_common_B_chi2(spots,accepted,d,vref);
    for (int i=0;i<180;++i) {
        if (fc<fd) { hi=d; d=c; fd=fc; c=hi-gr*(hi-lo); fc=profile_common_B_chi2(spots,accepted,c,vref); }
        else       { lo=c; c=d; fc=fd; d=lo+gr*(hi-lo); fd=profile_common_B_chi2(spots,accepted,d,vref); }
    }
    return 0.5*(lo+hi);
}

static double profile_B_crossing(const map<int,vector<AnalysisRow>>& spots,
                                 const std::set<int>& accepted,double vref,
                                 double best,double chi2min,
                                 int direction,double first_step) {
    const double target=chi2min+1.0;
    double inner=best,outer=best+direction*first_step;
    double fouter=profile_common_B_chi2(spots,accepted,outer,vref);
    for (int i=0;i<60 && fouter<target;++i) {
        inner=outer;
        first_step*=1.7;
        outer=best+direction*first_step;
        fouter=profile_common_B_chi2(spots,accepted,outer,vref);
    }
    if (!(fouter>=target) || !finite_number(fouter)) return std::numeric_limits<double>::quiet_NaN();
    double a=std::min(inner,outer),b=std::max(inner,outer);
    for (int i=0;i<100;++i) {
        const double m=0.5*(a+b);
        const double fm=profile_common_B_chi2(spots,accepted,m,vref);
        if (direction<0) { if (fm>=target) a=m; else b=m; }
        else             { if (fm>=target) b=m; else a=m; }
    }
    return 0.5*(a+b);
}

static BSimultaneousResult simultaneous_B_fit(const map<int,vector<AnalysisRow>>& spots,
                                              const std::set<int>& accepted,
                                              const vector<double>& individual_B,
                                              double vref) {
    BSimultaneousResult out; out.vref=vref;
    if (accepted.empty() || individual_B.empty() || !(vref>0.0)) return out;
    double mn=*std::min_element(individual_B.begin(),individual_B.end());
    double mx=*std::max_element(individual_B.begin(),individual_B.end());
    const double span=mx-mn;
    const double margin=std::max(2.0,2.0*span+0.5);
    double lo=mn-margin,hi=mx+margin,best=0.0;
    for (int k=0;k<5;++k) {
        best=golden_min_B(spots,accepted,vref,lo,hi);
        const double edge=0.02*(hi-lo);
        if (best-lo<edge) { const double w=hi-lo; lo-=w; continue; }
        if (hi-best<edge) { const double w=hi-lo; hi+=w; continue; }
        break;
    }
    int npoints=0,nspots=0;
    const double chi2=profile_common_B_chi2(spots,accepted,best,vref,&npoints,&nspots);
    if (!finite_number(chi2) || nspots<1 || npoints<=nspots) return out;
    const double step=std::max(1e-3,0.01*(hi-lo));
    const double left=profile_B_crossing(spots,accepted,vref,best,chi2,-1,step);
    const double right=profile_B_crossing(spots,accepted,vref,best,chi2,+1,step);
    double err=std::numeric_limits<double>::quiet_NaN();
    if (finite_number(left) && finite_number(right)) err=0.5*(right-left);
    else if (finite_number(left)) err=best-left;
    else if (finite_number(right)) err=right-best;
    out.ok=true; out.B=best; out.Berr=err; out.chi2=chi2;
    out.npoints=npoints; out.nspots=nspots; out.ndf=npoints-nspots-1;
    out.prob=(out.ndf>0)?TMath::Prob(chi2,out.ndf):0.0;
    return out;
}

// -----------------------------------------------------------------------------
// Write an analysis CSV using the new symmetric-systematic convention.
// -----------------------------------------------------------------------------
static void write_augmented_v_csv(const CsvTable& original,
                                  const string& phase,
                                  const map<int,VNominalFit>& fits,
                                  const string& output_csv) {
    std::ofstream fout(output_csv);
    if (!fout.is_open()) return;

    vector<int> kept_columns;
    bool input_has_deltaL=false;
    for (int i=0;i<(int)original.header.size();++i) {
        const string& name=original.header[i];
        if (name=="deltaL_plus" || name=="deltaL_minus") continue;
        if (name=="deltaL") input_has_deltaL=true;
        kept_columns.push_back(i);
    }

    bool first=true;
    for (int idx:kept_columns) {
        if (!first) fout << ',';
        fout << original.header[idx];
        first=false;
    }
    if (!input_has_deltaL) fout << ",deltaL";
    fout << ",B,B_stat_error,deltaB,A,A_stat_error,fit_chi2,fit_ndf,fit_chi2ndf,fit_status,covmatrix_status,edm,fit_converged\n";

    for (const auto& r:original.rows) {
        if (r.phase!=phase) continue;

        first=true;
        for (int idx:kept_columns) {
            if (!first) fout << ',';
            if (idx<(int)r.raw_fields.size()) fout << r.raw_fields[idx];
            first=false;
        }
        if (!input_has_deltaL) fout << ',' << r.deltaL;

        auto it=fits.find(r.spot);
        if (it==fits.end()) {
            fout << ",,,,,,,,,,,,0\n";
            continue;
        }

        const auto& f=it->second;
        fout << ',' << f.B
             << ',' << f.Berr
             << ',' << (f.has_param_syst ? f.deltaB : 0.0)
             << ',' << f.A
             << ',' << f.Aerr
             << ',' << f.chi2
             << ',' << f.ndf
             << ',' << f.chi2ndf
             << ',' << f.status
             << ',' << f.covstatus
             << ',' << f.edm
             << ',' << (f.converged ? 1 : 0)
             << '\n';
    }
}


// -----------------------------------------------------------------------------
// Build the optional R=16/R=24 overlays. The same power-law model and starting
// parameters used for the nominal R=20 fit are used here. "N" prevents ROOT
// from attaching extra fit-statistics boxes to the comparison graphs; the fit
// functions are drawn manually as dashed lines.
// -----------------------------------------------------------------------------
static map<int,VRadiusOverlayFit> build_v_radius_overlays(
    const string& filename,
    const string& phase,
    int color,
    int marker,
    const string& tag,
    double fit_xmin,
    double fit_xmax,
    double A_init,
    double B_init) {

    map<int,VRadiusOverlayFit> out;
    if (filename.empty()) return out;

    CsvTable table=read_analysis_csv(filename,true);
    auto selected=filter_phase_detected(table.rows,phase);
    if (selected.empty()) {
        std::cerr << "Warning: no comparison rows for phase " << phase
                  << " in " << filename << std::endl;
        return out;
    }

    map<int,vector<AnalysisRow>> spots;
    for (const auto& r:selected) spots[r.spot].push_back(r);

    for (auto& entry:spots) {
        int spot=entry.first;
        auto rows=entry.second;
        std::sort(rows.begin(),rows.end(),
                  [](const AnalysisRow& a,const AnalysisRow& b){return a.v_fin<b.v_fin;});

        VRadiusOverlayFit info;
        info.rows=rows;

        int n=(int)rows.size();
        if (n<2) {
            out[spot]=info;
            continue;
        }

        vector<double> xv(n),yv(n),exv(n,0.0),eyv(n);
        for (int i=0;i<n;++i) {
            xv[i]=rows[i].v_fin;
            yv[i]=rows[i].luminosity;
            eyv[i]=rows[i].error;
        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
        info.graph->SetName(Form("gr_%s_spot_%d",tag.c_str(),spot));
        info.graph->SetMarkerColor(color);
        info.graph->SetLineColor(color);
        info.graph->SetMarkerStyle(marker);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetLineWidth(2);

        double vref_final=1.0;
        TFitResultPtr fit=fit_power_centered_nom(info.graph,rows,
                                                 Form("fit_%s_spot_%d",tag.c_str(),spot),
                                                 fit_xmin,fit_xmax,A_init,B_init,
                                                 "QRSN",info.func,vref_final);
        info.func->SetLineColor(color);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(2);
        info.fit_done=true;
        info.status=(int)fit;
        info.B=info.func->GetParameter(1);
        info.Berr=info.func->GetParError(1);
        info.converged=(info.status==0 && finite_number(info.B));

        out[spot]=info;
    }

    return out;
}

void lum_vs_v_fit(const char* all_phases_csv,
                  const char* phase,
                  const char* systematic_values_csv,
                  const char* output_dir,
                  const char* prefix,
                  const char* r16_all_phases_csv = "",
                  const char* r24_all_phases_csv = "") {
    gStyle->SetOptFit(111);

    const double FIT_XMIN=0.0;
    const double FIT_XMAX=8.0;
    const double A_INIT=1.0;
    const double B_INIT=2.0;

    string ph=phase, outdir=output_dir, pref=prefix;
    const string r16_file=(r16_all_phases_csv ? r16_all_phases_csv : "");
    const string r24_file=(r24_all_phases_csv ? r24_all_phases_csv : "");
    const bool comparison_mode=(!r16_file.empty() || !r24_file.empty());
    ensure_dir(outdir);

    CsvTable table=read_analysis_csv(all_phases_csv,true);
    auto selected=filter_phase_detected(table.rows,ph);
    if (selected.empty()) {
        std::cerr << "No rows for phase " << ph << std::endl;
        return;
    }

    auto sys=read_B_systematics(systematic_values_csv);

    map<int,vector<AnalysisRow>> spots;
    for (const auto& r:selected) spots[r.spot].push_back(r);

    map<int,VNominalFit> fits;

    // =====================================================================
    // FIT STAGE: each hotspot is fitted once and the result is reused below.
    // =====================================================================
    for (auto& entry:spots) {
        int spot=entry.first;
        auto rows=entry.second;
        std::sort(rows.begin(),rows.end(),
                  [](const AnalysisRow& a,const AnalysisRow& b){return a.v_fin<b.v_fin;});

        VNominalFit info;
        info.spot=spot;
        info.rows=rows;

        int n=(int)rows.size();
        if (n<2) {
            fits[spot]=info;
            continue;
        }

        vector<double> xv(n), yv(n), exv(n,0.0), eyv(n);
        for (int i=0;i<n;++i) {
            xv[i]=rows[i].v_fin;
            yv[i]=rows[i].luminosity;
            eyv[i]=rows[i].error;
        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
        info.graph->SetName(Form("gr_spot_%d",spot));
        info.graph->SetTitle(Form("%s - Spot %d: x = %.2f, y = %.2f, T = %.1f #circC, %s;Overvoltage (V);Luminosity",
                                  pref.c_str(),spot,rows[0].x,rows[0].y,rows[0].T,ph.c_str()));
        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetLineWidth(2);

        double vref_final=1.0;
        TFitResultPtr fit=fit_power_centered_nom(info.graph,rows,
                                                 Form("f_quad_spot_%d",spot),
                                                 FIT_XMIN,FIT_XMAX,A_INIT,B_INIT,
                                                 "QRS0",info.func,vref_final);
        info.fit_done=true;
        info.status=(int)fit;
        info.covstatus=fit->CovMatrixStatus();
        info.edm=fit->Edm();
        info.B=info.func->GetParameter(1);
        info.Berr=info.func->GetParError(1);
        convert_power_Aref_to_A(fit,info.func->GetParameter(0),info.func->GetParError(0),
                                info.B,info.Berr,vref_final,info.A,info.Aerr);
        info.chi2=info.func->GetChisquare();
        info.ndf=info.func->GetNDF();
        info.chi2ndf=(info.ndf>0) ? info.chi2/info.ndf : 0.0;
        info.prob=info.func->GetProb();

        // A point enters Canvas 9 when the ROOT minimisation converged and B is
        // finite. No chi2/ndf or covariance-status cut is applied here.
        info.converged=(info.status==0 && finite_number(info.B));

        // Symmetric systematic on B:
        // deltaB = max(|B24-B20|, |B20-B16|).
        if (sys.count(spot) && sys[spot].has_delta) {
            // deltaB is computed upstream from R16/R20/R24 fits performed on
            // the exact same common set of detected overvoltage points.
            info.deltaB=sys[spot].deltaB;
            info.has_param_syst=true;
        }

        fits[spot]=info;
    }

    map<int,VRadiusOverlayFit> fits16;
    map<int,VRadiusOverlayFit> fits24;
    if (!r16_file.empty()) {
        fits16=build_v_radius_overlays(r16_file,ph,kP6Blue,24,"R16",
                                       FIT_XMIN,FIT_XMAX,A_INIT,B_INIT);
    }
    if (!r24_file.empty()) {
        fits24=build_v_radius_overlays(r24_file,ph,kP6Red,25,"R24",
                                       FIT_XMIN,FIT_XMAX,A_INIT,B_INIT);
    }

    // =====================================================================
    // SINGLE-HOTSPOT CANVASES
    // =====================================================================
    string spotdir=outdir+"/lum_vs_v_all_spots";
    ensure_dir(spotdir);

    for (auto& entry:fits) {
        auto& info=entry.second;
        int id=entry.first;
        if (!info.fit_done || !info.graph || !info.func) continue;

        TCanvas* c=new TCanvas(Form("c_v_spot_%d",id),
                               Form("Spot %d - Luminosity vs overvoltage",id),
                               1600,1200);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.13);
        c->SetTopMargin(0.08);

        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.8);
        info.graph->SetMarkerColor(kBlack);
        info.graph->SetLineColor(kBlack);
        info.graph->SetLineWidth(2);
        info.func->SetLineColor(kBlack);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(1);
        info.graph->GetXaxis()->SetTitleSize(0.045);
        info.graph->GetYaxis()->SetTitleSize(0.045);
        info.graph->GetYaxis()->SetTitleOffset(1.4);

        // Make sure R=20 (stat+syst) and the optional R=16/R=24
        // statistical error bars all fit inside the visible y range.
        double ymin=1e99, ymax=-1e99;
        for (const auto& r:info.rows) {
            const double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));
            ymin=std::min(ymin,r.luminosity-emax);
            ymax=std::max(ymax,r.luminosity+emax);
        }

        auto expand_overlay_range=[&](const map<int,VRadiusOverlayFit>& overlays) {
            auto it=overlays.find(id);
            if (it==overlays.end()) return;
            for (const auto& r:it->second.rows) {
                ymin=std::min(ymin,r.luminosity-std::fabs(r.error));
                ymax=std::max(ymax,r.luminosity+std::fabs(r.error));
            }
        };
        expand_overlay_range(fits16);
        expand_overlay_range(fits24);

        if (finite_number(ymin) && finite_number(ymax) && ymax>ymin) {
            double span=ymax-ymin;
            info.graph->SetMinimum(ymin-0.08*span);
            info.graph->SetMaximum(ymax+0.15*span);
        }

        info.graph->Draw("AP");
        c->Update();

        vector<double> sx, sy, sd;
        for (const auto& r:info.rows) {
            sx.push_back(r.v_fin);
            sy.push_back(r.luminosity);
            sd.push_back(r.deltaL);
        }
        draw_horizontal_syst_brackets(sx,sy,sd,
                                      info.graph->GetXaxis()->GetXmin(),
                                      info.graph->GetXaxis()->GetXmax());

        // Experimental points and fit remain on top of the systematic brackets.
        info.graph->Draw("P SAME");
        info.func->Draw("SAME");

        VRadiusOverlayFit* ov16=nullptr;
        VRadiusOverlayFit* ov24=nullptr;

        auto it16=fits16.find(id);
        if (it16!=fits16.end() && it16->second.graph) {
            ov16=&it16->second;
            ov16->graph->Draw("PE SAME");
            if (ov16->fit_done && ov16->func) ov16->func->Draw("SAME");
        }

        auto it24=fits24.find(id);
        if (it24!=fits24.end() && it24->second.graph) {
            ov24=&it24->second;
            ov24->graph->Draw("PE SAME");
            if (ov24->fit_done && ov24->func) ov24->func->Draw("SAME");
        }

        TLine* syst_proxy=new TLine(0,0,1,0);
        syst_proxy->SetLineColor(kBlack);
        syst_proxy->SetLineWidth(4);

        TLegend* leg=nullptr;
        if (!comparison_mode) {
            // Original 5analysis legend.
            leg=new TLegend(0.12,0.73,0.48,0.91);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(info.graph,"Data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"Systematic uncertainty","l");
            leg->AddEntry(info.func,"Fit: Lum = A_{ref}(V/V_{ref})^{B}","l");
        } else {
            // Radius-comparison legend: R=20 remains the nominal dataset.
            leg=new TLegend(0.12,0.67,0.66,0.91);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.030);

            leg->AddEntry(info.graph,"R = 20 px data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"R = 20 px luminosity systematic uncertainty","l");

            string nominal_fit_label;
            if (info.has_param_syst) {
                nominal_fit_label=Form("R = 20 px fit: B = %.3f #pm %.3f (stat) #pm %.3f (syst)",
                                       info.B,info.Berr,info.deltaB);
            } else {
                nominal_fit_label=Form("R = 20 px fit: B = %.3f #pm %.3f",
                                       info.B,info.Berr);
            }
            leg->AddEntry(info.func,nominal_fit_label.c_str(),"l");

            if (ov16) {
                leg->AddEntry(ov16->graph,"R = 16 px data (statistical uncertainty)","lep");
                if (ov16->fit_done && ov16->func) {
                    string label16=Form("R = 16 px fit: B = %.3f #pm %.3f",
                                        ov16->B,ov16->Berr);
                    leg->AddEntry(ov16->func,label16.c_str(),"l");
                }
            }

            if (ov24) {
                leg->AddEntry(ov24->graph,"R = 24 px data (statistical uncertainty)","lep");
                if (ov24->fit_done && ov24->func) {
                    string label24=Form("R = 24 px fit: B = %.3f #pm %.3f",
                                        ov24->B,ov24->Berr);
                    leg->AddEntry(ov24->func,label24.c_str(),"l");
                }
            }
        }
        leg->Draw();

        c->Update();
        
        TPaveStats* stats=(TPaveStats*)info.graph->FindObject("stats");
        if (!stats) {
            // Fallback: explicitly create a TPaveStats so fit information is
            // always visible even when ROOT does not attach an automatic box
            // after a quiet/zero-draw fit.
            stats=new TPaveStats(0.12,comparison_mode ? 0.43 : 0.51,
                                 comparison_mode ? 0.66 : 0.48,
                                 comparison_mode ? 0.65 : 0.71,"brNDC");
            stats->SetName(Form("fit_stats_spot_%d",id));
            stats->AddText(Form("A = %.4g #pm %.3g",info.A,info.Aerr));
            stats->AddText(Form("B = %.4g #pm %.3g",info.B,info.Berr));
            stats->AddText(Form("#chi^{2}/ndf = %.3g / %d",info.chi2,info.ndf));
        } else {
            stats->SetX1NDC(0.12);
            stats->SetX2NDC(comparison_mode ? 0.66 : 0.48);
            stats->SetY1NDC(comparison_mode ? 0.43 : 0.51);
            stats->SetY2NDC(comparison_mode ? 0.65 : 0.71);
        }
        stats->SetTextSize(0.020);
        stats->SetBorderSize(1);
        stats->SetFillStyle(0);
        stats->Draw();

        c->Modified();
        c->Update();
        c->SaveAs(Form("%s/%s_lum_vs_v_spot%d.png",spotdir.c_str(),pref.c_str(),id));
        c->SaveAs(Form("%s/%s_lum_vs_v_spot%d.pdf",spotdir.c_str(),pref.c_str(),id));
        delete c;
    }

    if (comparison_mode) {
        // The comparison pipeline is intended only to redraw the individual
        // luminosity curves with R=16/20/24 overlays. The nominal B-vs-spot
        // and CSV products remain the responsibility of 5analysis.sh.
        return;
    }

    // =====================================================================
    // CANVAS 9: B vs global hotspot ID.
    // The displayed individual B values retain the original chi2/ndf<4 cut.
    // The horizontal common B is obtained from a simultaneous fit to the
    // original L(V) points with B common and one free A_ref,i per hotspot.
    // =====================================================================
    vector<double> x, ex, B, Bstat, Bsyst;
    std::set<int> accepted_spots;
    for (auto& kv:fits) {
        const auto& f=kv.second;
        if (!f.fit_done ||
            !f.converged ||
            !finite_number(f.B) ||
            !finite_number(f.Berr) ||
            f.ndf <= 0 ||
            f.chi2ndf >= 4) {
            continue;
        }
        x.push_back((double)kv.first);
        ex.push_back(0.0);
        B.push_back(f.B);
        Bstat.push_back(f.Berr);
        Bsyst.push_back(f.has_param_syst ? f.deltaB : 0.0);
        accepted_spots.insert(kv.first);
    }

    double global_vref1=1.0,global_vref2=1.0;
    central_vrefs_nom(selected,global_vref1,global_vref2);
    BSimultaneousResult sim_B=simultaneous_B_fit(spots,accepted_spots,B,global_vref2);
    if (sim_B.ok) {
        string sim_csv=outdir+"/"+pref+"_B_simultaneous_fit.csv";
        std::ofstream fsim(sim_csv);
        fsim << "B_common,B_stat_error,chi2,ndf,chi2ndf,prob,n_spots,n_points,V_ref\n";
        fsim << sim_B.B << ',' << sim_B.Berr << ',' << sim_B.chi2 << ',' << sim_B.ndf << ','
             << (sim_B.ndf>0?sim_B.chi2/sim_B.ndf:0.0) << ',' << sim_B.prob << ','
             << sim_B.nspots << ',' << sim_B.npoints << ',' << sim_B.vref << '\n';
    }

    if (!B.empty()) {
        const string sensor=sensor_from_prefix(pref);
        double xmin=*std::min_element(x.begin(),x.end())-1.0;
        double xmax=*std::max_element(x.begin(),x.end())+1.0;
        double ymin=1e99,ymax=-1e99;
        for(size_t i=0;i<B.size();++i){const double emax=std::max(std::fabs(Bstat[i]),std::fabs(Bsyst[i]));ymin=std::min(ymin,B[i]-emax);ymax=std::max(ymax,B[i]+emax);}
        double yspan=ymax-ymin;if(!(yspan>0.0))yspan=std::max(0.2,std::fabs(ymax)*0.2);

        auto draw_B_canvas=[&](const string& tag,bool focus){
            TCanvas* c9=new TCanvas(Form("c9_%s",tag.c_str()),"B vs spot",1800,1000);
            TGraphErrors* grB=new TGraphErrors((int)B.size(),x.data(),B.data(),ex.data(),Bstat.data());
            grB->SetTitle(Form("Sensor %s - Power-law exponent B - %s;Global spot ID;B",sensor.c_str(),ph.c_str()));
            grB->SetMarkerStyle(21);
            grB->SetMarkerSize(1.2);
            grB->SetMarkerColor(kBlack);
            grB->SetLineColor(kBlack);
            grB->SetLineWidth(2);
            grB->SetMinimum(focus?1.0:ymin-0.15*yspan);
            grB->SetMaximum(focus?3.0:ymax+0.35*yspan);
            grB->Draw("AP");
            grB->GetXaxis()->SetLimits(xmin,xmax);

            draw_horizontal_syst_brackets(x,B,Bsyst,xmin,xmax,kP6Grape,4);grB->Draw("P SAME");
            TF1* line_B=new TF1(Form("line_B_%s",tag.c_str()),"2",xmin,xmax);
            line_B->SetLineColor(kP6Red);
            line_B->SetLineStyle(2);
            line_B->SetLineWidth(2);
            line_B->Draw("SAME");

            TF1* fit_B=nullptr;
            if (sim_B.ok) {
                fit_B=new TF1(Form("fit_B_common_%s",tag.c_str()),"[0]",xmin,xmax);
                fit_B->SetParameter(0,sim_B.B);
                fit_B->SetLineColor(kP10Gray);
                fit_B->SetLineWidth(2);
                fit_B->Draw("SAME");
            }
            c9->Update();

            TPaveStats* stats=(TPaveStats*)grB->FindObject("stats");
            if(!stats){
                stats=new TPaveStats(.68,.72,.94,.90,"brNDC");
                stats->SetName(Form("B_common_stats_%s",tag.c_str()));
                if(sim_B.ok){
                    stats->AddText(Form("B_{common} = %.4g #pm %.3g",sim_B.B,sim_B.Berr));
                    stats->AddText(Form("#chi^{2}/ndf = %.3g / %d",sim_B.chi2,sim_B.ndf));
                }
            } else {
                stats->SetX1NDC(.68);
                stats->SetX2NDC(.94);
                stats->SetY1NDC(.72);
                stats->SetY2NDC(.90);
            }

            stats->SetTextSize(.022);
            stats->SetFillStyle(0);
            stats->Draw();
            TLine* syst_proxy=new TLine(0,0,1,0);syst_proxy->SetLineColor(kP6Grape);syst_proxy->SetLineWidth(4);
            TLegend* leg=new TLegend(.12,.70,.48,.91);leg->SetBorderSize(0);leg->SetFillStyle(0);
            leg->AddEntry(grB,"B statistical uncertainty","lep");
            leg->AddEntry(syst_proxy,"B systematic uncertainty","l");
            if(fit_B)leg->AddEntry(fit_B,"Simultaneous fit: common B, free A_{i}","l");
            leg->AddEntry(line_B,"Line: B = 2","l");leg->Draw();
            c9->Modified();c9->Update();
            string base=outdir+"/"+pref+"_B_vs_spot"+(focus?"_focus_1_3":"");c9->SaveAs((base+".png").c_str());c9->SaveAs((base+".pdf").c_str());delete c9;
        };
        draw_B_canvas("full",false);
        draw_B_canvas("focus",true);
    } else {
        std::cerr << "Warning: no ROOT-converged B fit is available for Canvas 9 in phase "
                  << ph << std::endl;
    }

    // =====================================================================
    // CSV OUTPUTS
    // =====================================================================
    string result_csv=outdir+"/"+pref+"_B_fit_results.csv";
    std::ofstream fout(result_csv);
    fout << "spot,x,y,phase,A,A_stat_error,B,B_stat_error,deltaB,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";

    for (auto& kv:fits) {
        const auto& f=kv.second;
        double xx=f.rows.empty()?0.0:f.rows[0].x;
        double yy=f.rows.empty()?0.0:f.rows[0].y;
        fout << kv.first << ',' << xx << ',' << yy << ',' << ph
             << ',' << f.A << ',' << f.Aerr
             << ',' << f.B << ',' << f.Berr
             << ',' << (f.has_param_syst ? f.deltaB : 0.0)
             << ',' << f.chi2 << ',' << f.ndf << ',' << f.chi2ndf
             << ',' << f.prob << ',' << f.status << ',' << f.covstatus
             << ',' << f.edm << ',' << (f.converged?1:0) << '\n';
    }
    fout.close();

    write_augmented_v_csv(table,ph,fits,outdir+"/"+pref+"_analysis.csv");
}
