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
#include "TLine.h"
#include "TMath.h"

struct LambdaSystematicPair {
    bool has16=false,has24=false;
    double l16=0.0,l24=0.0;
    bool conv16=false,conv24=false;
    bool has_delta=false;
    double deltaLambda=0.0;
};

struct TNominalFit {
    int spot=-1;
    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TGraphErrors* syst_graph=nullptr;
    TF1* func=nullptr;
    double A=0.0,Aerr=0.0,lambda=0.0,lambdaerr=0.0;
    double chi2=0.0,chi2ndf=0.0,prob=0.0,edm=-1.0;
    int ndf=0,status=-999,covstatus=-999;
    bool fit_done=false,converged=false;
    bool has_param_syst=false;
    double deltaLambda=0.0;
};

struct TRadiusOverlayFit {
    vector<AnalysisRow> rows;
    TGraphErrors* graph=nullptr;
    TF1* func=nullptr;
    double lambda=0.0,lambdaerr=0.0;
    int status=-999;
    bool fit_done=false,converged=false;
};

// -----------------------------------------------------------------------------
// Draw a symmetric systematic uncertainty as two green horizontal brackets.
// The horizontal cap is intentionally only slightly wider than the data marker,
// so statistical and systematic uncertainties remain visually distinct.
// -----------------------------------------------------------------------------
static void draw_horizontal_syst_brackets_T(const vector<double>& x,
                                             const vector<double>& y,
                                             const vector<double>& delta,
                                             double xmin,
                                             double xmax,
                                             int color=kBlack,
                                             int line_width=4) {
    if (x.empty() || y.size()!=x.size() || delta.size()!=x.size()) return;

    double span=xmax-xmin;
    if (!(span>0.0)) span=1.0;

    const double half_width=0.0045*span;

    for (size_t i=0;i<x.size();++i) {
        if (!finite_number(delta[i]) || delta[i]<=0.0) continue;

        const double ytop=y[i]+delta[i];
        const double ybottom=y[i]-delta[i];
        const double hook=0.10*delta[i];

        TLine* top=new TLine(x[i]-half_width,ytop,x[i]+half_width,ytop);
        top->SetLineColor(color); top->SetLineWidth(line_width); top->Draw("SAME");
        TLine* top_l=new TLine(x[i]-half_width,ytop,x[i]-half_width,ytop-hook);
        top_l->SetLineColor(color); top_l->SetLineWidth(line_width); top_l->Draw("SAME");
        TLine* top_r=new TLine(x[i]+half_width,ytop,x[i]+half_width,ytop-hook);
        top_r->SetLineColor(color); top_r->SetLineWidth(line_width); top_r->Draw("SAME");

        TLine* bottom=new TLine(x[i]-half_width,ybottom,x[i]+half_width,ybottom);
        bottom->SetLineColor(color); bottom->SetLineWidth(line_width); bottom->Draw("SAME");
        TLine* bottom_l=new TLine(x[i]-half_width,ybottom,x[i]-half_width,ybottom+hook);
        bottom_l->SetLineColor(color); bottom_l->SetLineWidth(line_width); bottom_l->Draw("SAME");
        TLine* bottom_r=new TLine(x[i]+half_width,ybottom,x[i]+half_width,ybottom+hook);
        bottom_r->SetLineColor(color); bottom_r->SetLineWidth(line_width); bottom_r->Draw("SAME");
    }
}

static void estimate_exp_parameters_nom(const vector<AnalysisRow>& rows,
                                        double Tref,
                                        double& Aref0,
                                        double& lambda0) {
    // Weighted regression of ln(L) versus T-Tref using all positive points.
    // sigma_lnL = sigma_L/L, hence w = (L/sigma_L)^2 when sigma_L is available.
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

// Two centered nonlinear steps: first Tref=19 C, then Tref=21 C.
// The final TF1 is always the 21 C parameterization and A_ref is constrained positive.
static TFitResultPtr fit_exp_centered_nom(TGraphErrors* graph,
                                         const vector<AnalysisRow>& rows,
                                         const string& final_name,
                                         double fitmin,
                                         double fitmax,
                                         const char* final_options,
                                         TF1*& final_func) {
    double A19,l0;
    estimate_exp_parameters_nom(rows,19.0,A19,l0);

    TF1 f19((final_name+"_prefit19").c_str(),"[0]*exp([1]*(x-19.0))",fitmin,fitmax);
    f19.SetParNames("A_{ref,19}","#lambda");
    f19.SetParameters(A19,l0);
    f19.SetParLimits(0,1e-300,1e300);
    TFitResultPtr fit19=graph->Fit(&f19,"QRSN");

    const double l1=((int)fit19==0 && finite_number(f19.GetParameter(1))) ? f19.GetParameter(1) : l0;
    const double A19fit=((int)fit19==0 && f19.GetParameter(0)>0.0 && finite_number(f19.GetParameter(0))) ? f19.GetParameter(0) : A19;
    double A21=A19fit*std::exp(l1*(21.0-19.0));
    if (!(A21>0.0) || !finite_number(A21)) A21=A19;

    final_func=new TF1(final_name.c_str(),"[0]*exp([1]*(x-21.0))",fitmin,fitmax);
    final_func->SetParNames("A_{ref,21}","#lambda");
    final_func->SetParameters(A21,l1);
    final_func->SetParLimits(0,1e-300,1e300);
    return graph->Fit(final_func,final_options);
}

// Keep the historical output parameter A of L=A exp(lambda T), even though
// the minimization is performed with the better-conditioned A_ref at Tref=21 C.
static void convert_exp_Aref_to_A(const TFitResultPtr& fit,
                                  double Aref,double ArefErr,
                                  double lambda,double lambdaErr,
                                  double& A,double& Aerr) {
    const double Tref=21.0;
    const double scale=std::exp(-lambda*Tref);
    A=Aref*scale;
    double cov=0.0;
    if ((int)fit>=0) cov=fit->CovMatrix(0,1);
    const double dA_dAref=scale;
    const double dA_dlambda=-Tref*A;
    double var=dA_dAref*dA_dAref*ArefErr*ArefErr
              +dA_dlambda*dA_dlambda*lambdaErr*lambdaErr
              +2.0*dA_dAref*dA_dlambda*cov;
    if (var<0.0 && std::fabs(var)<1e-12*std::max(1.0,A*A)) var=0.0;
    Aerr=(var>=0.0 && finite_number(var)) ? std::sqrt(var) : std::fabs(scale*ArefErr);
}

struct TSimultaneousResult {
    bool ok=false;
    double lambda=std::numeric_limits<double>::quiet_NaN();
    double lambdaerr=std::numeric_limits<double>::quiet_NaN();
    double chi2=std::numeric_limits<double>::quiet_NaN();
    double prob=0.0;
    int ndf=0,nspots=0,npoints=0;
};

// Profile chi2 for a common lambda. For each trial lambda the best A_ref,i is
// obtained analytically, so this is exactly the simultaneous weighted least-
// squares fit with one common lambda and one free amplitude per hotspot.
static double profile_common_lambda_chi2(const map<int,vector<AnalysisRow>>& spots,
                                         const std::set<int>& accepted,
                                         double lambda,
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
            if (!finite_number(r.T) || !finite_number(r.luminosity)) continue;
            const double sigma=(finite_number(r.error) && r.error>0.0) ? r.error : 1.0;
            const double expo=lambda*(r.T-21.0);
            if (expo>700.0 || expo<-700.0) return 1e300;
            const double g=std::exp(expo);
            const double w=1.0/(sigma*sigma);
            num+=w*r.luminosity*g;
            den+=w*g*g;
            ++nv;
        }
        if (nv<3 || !(den>0.0)) continue;
        double Aref=num/den;
        if (!(Aref>0.0) || !finite_number(Aref)) Aref=1e-300;
        for (const auto& r:rows) {
            if (!finite_number(r.T) || !finite_number(r.luminosity)) continue;
            const double sigma=(finite_number(r.error) && r.error>0.0) ? r.error : 1.0;
            const double expo=lambda*(r.T-21.0);
            if (expo>700.0 || expo<-700.0) return 1e300;
            const double model=Aref*std::exp(expo);
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

static double golden_min_lambda(const map<int,vector<AnalysisRow>>& spots,
                                const std::set<int>& accepted,
                                double lo,double hi) {
    const double gr=0.6180339887498948482;
    double c=hi-gr*(hi-lo),d=lo+gr*(hi-lo);
    double fc=profile_common_lambda_chi2(spots,accepted,c),fd=profile_common_lambda_chi2(spots,accepted,d);
    for (int i=0;i<180;++i) {
        if (fc<fd) { hi=d; d=c; fd=fc; c=hi-gr*(hi-lo); fc=profile_common_lambda_chi2(spots,accepted,c); }
        else       { lo=c; c=d; fc=fd; d=lo+gr*(hi-lo); fd=profile_common_lambda_chi2(spots,accepted,d); }
    }
    return 0.5*(lo+hi);
}

static double profile_lambda_crossing(const map<int,vector<AnalysisRow>>& spots,
                                      const std::set<int>& accepted,
                                      double best,double chi2min,
                                      int direction,double first_step) {
    const double target=chi2min+1.0;
    double inner=best,outer=best+direction*first_step;
    double fouter=profile_common_lambda_chi2(spots,accepted,outer);
    for (int i=0;i<60 && fouter<target;++i) {
        inner=outer;
        first_step*=1.7;
        outer=best+direction*first_step;
        fouter=profile_common_lambda_chi2(spots,accepted,outer);
    }
    if (!(fouter>=target) || !finite_number(fouter)) return std::numeric_limits<double>::quiet_NaN();
    double a=std::min(inner,outer),b=std::max(inner,outer);
    for (int i=0;i<100;++i) {
        double m=0.5*(a+b);
        double fm=profile_common_lambda_chi2(spots,accepted,m);
        if (direction<0) { if (fm>=target) a=m; else b=m; }
        else             { if (fm>=target) b=m; else a=m; }
    }
    return 0.5*(a+b);
}

static TSimultaneousResult simultaneous_lambda_fit(const map<int,vector<AnalysisRow>>& spots,
                                                   const std::set<int>& accepted,
                                                   const vector<double>& individual_lambda) {
    TSimultaneousResult out;
    if (accepted.empty() || individual_lambda.empty()) return out;
    double mn=*std::min_element(individual_lambda.begin(),individual_lambda.end());
    double mx=*std::max_element(individual_lambda.begin(),individual_lambda.end());
    double span=mx-mn;
    double margin=std::max(0.25,2.0*span+0.05);
    double lo=mn-margin,hi=mx+margin,best=0.0;
    for (int k=0;k<5;++k) {
        best=golden_min_lambda(spots,accepted,lo,hi);
        const double edge=0.02*(hi-lo);
        if (best-lo<edge) { const double w=hi-lo; lo-=w; continue; }
        if (hi-best<edge) { const double w=hi-lo; hi+=w; continue; }
        break;
    }
    int npoints=0,nspots=0;
    const double chi2=profile_common_lambda_chi2(spots,accepted,best,&npoints,&nspots);
    if (!finite_number(chi2) || nspots<1 || npoints<=nspots) return out;
    const double step=std::max(1e-4,0.01*(hi-lo));
    const double left=profile_lambda_crossing(spots,accepted,best,chi2,-1,step);
    const double right=profile_lambda_crossing(spots,accepted,best,chi2,+1,step);
    double err=std::numeric_limits<double>::quiet_NaN();
    if (finite_number(left) && finite_number(right)) err=0.5*(right-left);
    else if (finite_number(left)) err=best-left;
    else if (finite_number(right)) err=right-best;
    out.ok=true; out.lambda=best; out.lambdaerr=err; out.chi2=chi2;
    out.npoints=npoints; out.nspots=nspots; out.ndf=npoints-nspots-1;
    out.prob=(out.ndf>0)?TMath::Prob(chi2,out.ndf):0.0;
    return out;
}

// -----------------------------------------------------------------------------
// Build optional R=16/R=24 temperature overlays. They use the same exponential
// model and initial-parameter estimate as the nominal R=20 analysis. "N" keeps
// ROOT from creating extra fit-statistics boxes for the comparison datasets.
// -----------------------------------------------------------------------------
static map<int,TRadiusOverlayFit> build_T_radius_overlays(
    const string& filename,
    const string& phase,
    int color,
    int marker,
    const string& tag) {

    map<int,TRadiusOverlayFit> out;
    if (filename.empty()) return out;

    CsvTable table=read_analysis_csv(filename,false);
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
                  [](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});

        TRadiusOverlayFit info;
        info.rows=rows;

        int n=(int)rows.size();
        if (n<3) {
            out[spot]=info;
            continue;
        }

        double Tmin=rows.front().T,Tmax=rows.back().T,dT=Tmax-Tmin;
        if (dT<=0) dT=1.0;
        double fitmin=Tmin-0.10*dT,fitmax=Tmax+0.10*dT;

        vector<double> xv(n),yv(n),exv(n,0.0),eyv(n);
        for (int i=0;i<n;++i) {
            xv[i]=rows[i].T;
            yv[i]=rows[i].luminosity;
            eyv[i]=rows[i].error;
        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
        info.graph->SetName(Form("gr_%s_spot_%d",tag.c_str(),spot));
        info.graph->SetMarkerColor(color);
        info.graph->SetLineColor(color);
        info.graph->SetMarkerStyle(marker);
        info.graph->SetMarkerSize(1.25);
        info.graph->SetLineWidth(2);

        TFitResultPtr fit=fit_exp_centered_nom(info.graph,rows,
                                               Form("fit_%s_spot_%d",tag.c_str(),spot),
                                               fitmin,fitmax,"QRSN",info.func);
        info.func->SetLineColor(color);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(2);
        info.fit_done=true;
        info.status=(int)fit;
        info.lambda=info.func->GetParameter(1);
        info.lambdaerr=info.func->GetParError(1);
        info.converged=(info.status==0 && finite_number(info.lambda));

        out[spot]=info;
    }

    return out;
}

static map<int,LambdaSystematicPair> read_lambda_systematics(const string& filename) {
    map<int,LambdaSystematicPair> out;
    std::ifstream fin(filename);
    if (!fin.is_open()) return out;

    string line;
    if (!std::getline(fin,line)) return out;
    auto h=split_csv_simple(line);
    map<string,int> c;
    for (int i=0;i<(int)h.size();++i) c[h[i]]=i;

    while (std::getline(fin,line)) {
        if (trim_copy(line).empty()) continue;
        auto f=split_csv_simple(line);
        if (f.size()<h.size()) f.resize(h.size(),"");
        try {
            int id=(int)std::llround(std::stod(f[c["spot"]]));
            LambdaSystematicPair s;

            if (c.count("lambda_R16") && !f[c["lambda_R16"]].empty()) {
                s.l16=std::stod(f[c["lambda_R16"]]);
                s.has16=finite_number(s.l16);
            }
            if (c.count("lambda_R24") && !f[c["lambda_R24"]].empty()) {
                s.l24=std::stod(f[c["lambda_R24"]]);
                s.has24=finite_number(s.l24);
            }

            if (c.count("fit_status_R16") && !f[c["fit_status_R16"]].empty())
                s.conv16=(std::stoi(f[c["fit_status_R16"]])==0 && s.has16);
            else if (c.count("converged_R16"))
                s.conv16=parse_bool_safe(f[c["converged_R16"]]);

            if (c.count("fit_status_R24") && !f[c["fit_status_R24"]].empty())
                s.conv24=(std::stoi(f[c["fit_status_R24"]])==0 && s.has24);
            else if (c.count("converged_R24"))
                s.conv24=parse_bool_safe(f[c["converged_R24"]]);

            if (c.count("deltaLambda") && !f[c["deltaLambda"]].empty()) {
                s.deltaLambda=std::fabs(std::stod(f[c["deltaLambda"]]));
                s.has_delta=finite_number(s.deltaLambda);
            }
            out[id]=s;
        } catch (...) {
            // Skip malformed rows and continue with the remaining hotspots.
        }
    }
    return out;
}

// Write a phase CSV using only the new symmetric-systematic columns.
static void write_augmented_T_csv(const CsvTable& original,
                                  const string& phase,
                                  const map<int,TNominalFit>& fits,
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
    fout << ",lambda,lambda_stat_error,deltaLambda,A,A_stat_error,fit_chi2,fit_ndf,fit_chi2ndf,fit_status,covmatrix_status,edm,fit_converged\n";

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
        fout << ',' << f.lambda
             << ',' << f.lambdaerr
             << ',' << (f.has_param_syst?f.deltaLambda:0.0)
             << ',' << f.A
             << ',' << f.Aerr
             << ',' << f.chi2
             << ',' << f.ndf
             << ',' << f.chi2ndf
             << ',' << f.status
             << ',' << f.covstatus
             << ',' << f.edm
             << ',' << (f.converged?1:0)
             << '\n';
    }
}

void lum_vs_T_fit(const char* all_phases_csv,
                  const char* phase,
                  const char* systematic_values_csv,
                  const char* output_dir,
                  const char* prefix,
                  const char* r16_all_phases_csv = "",
                  const char* r24_all_phases_csv = "") {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);

    string ph=phase,outdir=output_dir,pref=prefix;
    const string r16_file=(r16_all_phases_csv ? r16_all_phases_csv : "");
    const string r24_file=(r24_all_phases_csv ? r24_all_phases_csv : "");
    const bool comparison_mode=(!r16_file.empty() || !r24_file.empty());
    ensure_dir(outdir);

    CsvTable table=read_analysis_csv(all_phases_csv,false);
    auto selected=filter_phase_detected(table.rows,ph);
    if (selected.empty()) {
        std::cerr << "No rows for phase " << ph << std::endl;
        return;
    }

    auto sys=read_lambda_systematics(systematic_values_csv);
    map<int,vector<AnalysisRow>> spots;
    for (const auto& r:selected) spots[r.spot].push_back(r);
    map<int,TNominalFit> fits;

    for (auto& entry:spots) {
        int spot=entry.first;
        auto rows=entry.second;
        std::sort(rows.begin(),rows.end(),[](const AnalysisRow&a,const AnalysisRow&b){return a.T<b.T;});

        TNominalFit info;
        info.spot=spot;
        info.rows=rows;
        int n=(int)rows.size();
        if (n<3) {
            fits[spot]=info;
            continue;
        }

        double Tmin=rows.front().T,Tmax=rows.back().T,dT=Tmax-Tmin;
        if (dT<=0) dT=1.0;
        double fitmin=Tmin-0.10*dT,fitmax=Tmax+0.10*dT;

        vector<double> xv(n),yv(n),exv(n,0.0),eyv(n),esyst(n);
        for (int i=0;i<n;++i) {
            xv[i]=rows[i].T;
            yv[i]=rows[i].luminosity;
            eyv[i]=rows[i].error;
            esyst[i]=rows[i].deltaL;
        }

        info.graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),eyv.data());
        info.graph->SetTitle(Form("%s - Spot %d: x = %.2f, y = %.2f, v = %.1f V, %s;T (#circC);Luminosity",
                                  pref.c_str(),spot,rows[0].x,rows[0].y,rows[0].v,ph.c_str()));
        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.4);
        info.graph->SetLineWidth(2);

        info.syst_graph=new TGraphErrors(n,xv.data(),yv.data(),exv.data(),esyst.data());
        info.syst_graph->SetLineColor(kBlack);
        info.syst_graph->SetLineWidth(4);
        info.syst_graph->SetMarkerSize(0);

        TFitResultPtr fit=fit_exp_centered_nom(info.graph,rows,
                                               Form("fit_exp_spot_%d",spot),
                                               fitmin,fitmax,"QRS0",info.func);
        info.func->SetLineColor(kBlack);
        info.func->SetLineWidth(2);

        info.fit_done=true;
        info.status=(int)fit;
        info.covstatus=fit->CovMatrixStatus();
        info.edm=fit->Edm();
        info.lambda=info.func->GetParameter(1);
        info.lambdaerr=info.func->GetParError(1);
        convert_exp_Aref_to_A(fit,info.func->GetParameter(0),info.func->GetParError(0),
                              info.lambda,info.lambdaerr,info.A,info.Aerr);
        info.chi2=info.func->GetChisquare();
        info.ndf=info.func->GetNDF();
        info.chi2ndf=(info.ndf>0)?info.chi2/info.ndf:0.0;
        info.prob=info.func->GetProb();
        info.converged=(info.status==0 && finite_number(info.lambda));

        // Symmetric systematic on lambda:
        // deltaLambda=max(|lambda24-lambda20|, |lambda20-lambda16|).
        if (sys.count(spot) && sys[spot].has_delta) {
            // Computed upstream from R16/R20/R24 fits using exactly the same
            // common set of genuinely-detected temperature points.
            info.deltaLambda=sys[spot].deltaLambda;
            info.has_param_syst=true;
        }

        fits[spot]=info;
    }

    map<int,TRadiusOverlayFit> fits16;
    map<int,TRadiusOverlayFit> fits24;
    if (!r16_file.empty()) {
        fits16=build_T_radius_overlays(r16_file,ph,kP6Blue,24,"R16");
    }
    if (!r24_file.empty()) {
        fits24=build_T_radius_overlays(r24_file,ph,kP6Red,25,"R24");
    }

    // Individual hotspot canvases. The original T-fit aesthetics are preserved.
    string spotdir=outdir+"/lum_vs_T_all_spots_expfit";
    ensure_dir(spotdir);
    for (auto& kv:fits) {
        auto& info=kv.second;
        int spot=kv.first;
        if (!info.fit_done || !info.graph) continue;

        TCanvas* c=new TCanvas(Form("c_T_spot_%d",spot),
                               Form("Spot %d - Luminosity vs T",spot),
                               1200,800);

        info.graph->SetMarkerColor(kBlack);
        info.graph->SetLineColor(kBlack);
        info.func->SetLineColor(kBlack);
        info.func->SetLineWidth(2);
        info.func->SetLineStyle(1);

        double ymin=1e99,ymax=-1e99;
        for (const auto& r:info.rows) {
            const double emax=std::max(std::fabs(r.error),std::fabs(r.deltaL));
            ymin=std::min(ymin,r.luminosity-emax);
            ymax=std::max(ymax,r.luminosity+emax);
        }

        auto expand_overlay_range=[&](const map<int,TRadiusOverlayFit>& overlays) {
            auto it=overlays.find(spot);
            if (it==overlays.end()) return;
            for (const auto& r:it->second.rows) {
                ymin=std::min(ymin,r.luminosity-std::fabs(r.error));
                ymax=std::max(ymax,r.luminosity+std::fabs(r.error));
            }
        };
        expand_overlay_range(fits16);
        expand_overlay_range(fits24);

        if (ymin>0) ymin*=0.80; else ymin=0.0;
        ymax*=1.25;

        info.graph->SetMinimum(ymin);
        info.graph->SetMaximum(ymax);
        info.graph->Draw("AP");
        c->Update();

        vector<double> sx,sy,sd;
        for (const auto& r:info.rows) {
            sx.push_back(r.T);
            sy.push_back(r.luminosity);
            sd.push_back(r.deltaL);
        }
        draw_horizontal_syst_brackets_T(sx,sy,sd,
                                        info.graph->GetXaxis()->GetXmin(),
                                        info.graph->GetXaxis()->GetXmax());

        // Redraw the experimental markers and the fit above the systematic brackets.
        info.graph->Draw("P SAME");
        info.func->Draw("SAME");

        TRadiusOverlayFit* ov16=nullptr;
        TRadiusOverlayFit* ov24=nullptr;

        auto it16=fits16.find(spot);
        if (it16!=fits16.end() && it16->second.graph) {
            ov16=&it16->second;
            ov16->graph->Draw("PE SAME");
            if (ov16->fit_done && ov16->func) ov16->func->Draw("SAME");
        }

        auto it24=fits24.find(spot);
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
            leg=new TLegend(0.14,0.68,0.60,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(info.graph,"Data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"Systematic uncertainty","l");
            leg->AddEntry(info.func,"Fit: Lum = A_{ref} e^{#lambda(T-21)}","l");
        } else {
            leg=new TLegend(0.14,0.48,0.66,0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.030);

            leg->AddEntry(info.graph,"R = 20 px data (statistical uncertainty)","lep");
            leg->AddEntry(syst_proxy,"R = 20 px luminosity systematic uncertainty","l");

            string nominal_fit_label;
            if (info.has_param_syst) {
                nominal_fit_label=Form("R = 20 px fit: #lambda = %.5f #pm %.5f (stat) #pm %.5f (syst)",
                                       info.lambda,info.lambdaerr,info.deltaLambda);
            } else {
                nominal_fit_label=Form("R = 20 px fit: #lambda = %.5f #pm %.5f",
                                       info.lambda,info.lambdaerr);
            }
            leg->AddEntry(info.func,nominal_fit_label.c_str(),"l");

            if (ov16) {
                leg->AddEntry(ov16->graph,"R = 16 px data (statistical uncertainty)","lep");
                if (ov16->fit_done && ov16->func) {
                    string label16=Form("R = 16 px fit: #lambda = %.5f #pm %.5f",
                                        ov16->lambda,ov16->lambdaerr);
                    leg->AddEntry(ov16->func,label16.c_str(),"l");
                }
            }

            if (ov24) {
                leg->AddEntry(ov24->graph,"R = 24 px data (statistical uncertainty)","lep");
                if (ov24->fit_done && ov24->func) {
                    string label24=Form("R = 24 px fit: #lambda = %.5f #pm %.5f",
                                        ov24->lambda,ov24->lambdaerr);
                    leg->AddEntry(ov24->func,label24.c_str(),"l");
                }
            }
        }
        leg->Draw();

        c->SaveAs(Form("%s/%s_lum_vs_T_spot%d.png",spotdir.c_str(),pref.c_str(),spot));
        c->SaveAs(Form("%s/%s_lum_vs_T_spot%d.pdf",spotdir.c_str(),pref.c_str(),spot));
        delete c;
    }

    if (comparison_mode) {
        // Only the individual R=16/20/24 comparison canvases are generated in
        // comparison mode. Nominal lambda-vs-spot and CSV products remain in
        // the standard 5analysis.sh workflow.
        return;
    }

    // lambda vs spot: keep the original quality selection for the displayed
    // individual values. The horizontal result is NOT a fit of these lambda_i:
    // it is a simultaneous fit of the original L(T) points with common lambda
    // and one free amplitude A_ref,i per accepted hotspot.
    vector<double> x,ex,L,Lstat,Lsyst;
    std::set<int> accepted_spots;
    for (auto& kv:fits) {
        const auto& f=kv.second;
        if (!f.fit_done ||
            !finite_number(f.lambda) ||
            !finite_number(f.lambdaerr) ||
            f.ndf <= 0 ||
            f.chi2ndf >= 4) {
            continue;
        }
        x.push_back(kv.first);
        ex.push_back(0.0);
        L.push_back(f.lambda);
        Lstat.push_back(f.lambdaerr);
        Lsyst.push_back(f.has_param_syst?f.deltaLambda:0.0);
        accepted_spots.insert(kv.first);
    }

    TSimultaneousResult sim_lambda=simultaneous_lambda_fit(spots,accepted_spots,L);
    if (sim_lambda.ok) {
        string sim_csv=outdir+"/"+pref+"_lambda_simultaneous_fit.csv";
        std::ofstream fsim(sim_csv);
        fsim << "lambda_common,lambda_stat_error,chi2,ndf,chi2ndf,prob,n_spots,n_points,T_ref\n";
        fsim << sim_lambda.lambda << ',' << sim_lambda.lambdaerr << ','
             << sim_lambda.chi2 << ',' << sim_lambda.ndf << ','
             << (sim_lambda.ndf>0?sim_lambda.chi2/sim_lambda.ndf:0.0) << ','
             << sim_lambda.prob << ',' << sim_lambda.nspots << ',' << sim_lambda.npoints << ",21\n";
    }

    if (!L.empty()) {
        TCanvas* c5=new TCanvas("c5_lambda_vs_spot_constfit","lambda vs spot ID",1800,1000);
        TGraphErrors* gr=new TGraphErrors((int)L.size(),x.data(),L.data(),ex.data(),Lstat.data());
        gr->SetTitle(Form("%s - Exponential fit parameter #lambda - %s;Global spot ID;#lambda",pref.c_str(),ph.c_str()));
        gr->SetMarkerStyle(21);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);
        gr->Draw("AP");

        double xmin=*std::min_element(x.begin(),x.end())-1.0;
        double xmax=*std::max_element(x.begin(),x.end())+1.0;
        gr->GetXaxis()->SetLimits(xmin,xmax);

        double lmin=*std::min_element(L.begin(),L.end());
        double lmax=*std::max_element(L.begin(),L.end());
        double span=lmax-lmin;
        if (span<=0) span=std::max(0.1,std::fabs(lmax)*0.2);
        gr->GetYaxis()->SetRangeUser(lmin-0.25*span,lmax+0.90*span);

        c5->Update();
        draw_horizontal_syst_brackets_T(x,L,Lsyst,xmin,xmax,kBlack,4);
        gr->Draw("P SAME");

        TF1* fc=nullptr;
        if (sim_lambda.ok) {
            fc=new TF1("fit_lambda_common","[0]",xmin,xmax);
            fc->SetParNames("lambda_{common}");
            fc->SetParameter(0,sim_lambda.lambda);
            fc->SetLineColor(kBlack);
            fc->SetLineWidth(2);
            fc->Draw("SAME");
        }

        TLine* syst_proxy_lambda=new TLine(0,0,1,0);
        syst_proxy_lambda->SetLineColor(kBlack);
        syst_proxy_lambda->SetLineWidth(4);

        TLegend* leg=new TLegend(0.10,0.68,0.60,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr,"Individual #lambda (statistical uncertainty)","lep");
        leg->AddEntry(syst_proxy_lambda,"Systematic uncertainty (R=16/24)","l");
        if (fc) {
            leg->AddEntry(fc,Form("Simultaneous fit: #lambda = %.5f #pm %.5f, #chi^{2}/ndf = %.1f/%d",
                                  sim_lambda.lambda,sim_lambda.lambdaerr,sim_lambda.chi2,sim_lambda.ndf),"l");
        }
        leg->Draw();

        c5->SaveAs(Form("%s/%s_lambda_vs_spot.png",outdir.c_str(),pref.c_str()));
        c5->SaveAs(Form("%s/%s_lambda_vs_spot.pdf",outdir.c_str(),pref.c_str()));
        delete c5;
    }

    string result_csv=outdir+"/"+pref+"_lambda_fit_results.csv";
    std::ofstream fout(result_csv);
    fout << "spot,x,y,phase,A,A_stat_error,lambda,lambda_stat_error,deltaLambda,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";
    for (auto& kv:fits) {
        const auto& f=kv.second;
        double xx=f.rows.empty()?0.0:f.rows[0].x;
        double yy=f.rows.empty()?0.0:f.rows[0].y;
        fout << kv.first << ',' << xx << ',' << yy << ',' << ph
             << ',' << f.A << ',' << f.Aerr
             << ',' << f.lambda << ',' << f.lambdaerr
             << ',' << (f.has_param_syst?f.deltaLambda:0.0)
             << ',' << f.chi2 << ',' << f.ndf << ',' << f.chi2ndf
             << ',' << f.prob << ',' << f.status << ',' << f.covstatus
             << ',' << f.edm << ',' << (f.converged?1:0) << '\n';
    }
    fout.close();

    write_augmented_T_csv(table,ph,fits,outdir+"/"+pref+"_analysis.csv");
}
