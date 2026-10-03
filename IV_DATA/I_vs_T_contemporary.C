//root -l 'I_vs_T_contemporary.C("A1_I_vs_T.csv","analysis","A1")'


#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TAxis.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using std::map;
using std::set;
using std::string;
using std::vector;


// ============================================================================
// Structures & Helpers
// ============================================================================

namespace {

struct IVRow {

    string phase;

    double temperature_C = 0.0;
    double current = 0.0;
    double current_error = 0.0;
};


// Comparatore mancante aggiunto qui:
bool IVRow_less(const IVRow& a, const IVRow& b) {
    return a.temperature_C < b.temperature_C;
}


struct IVFitResult {

    string phase;

    vector<IVRow> rows;

    TGraphErrors* graph = nullptr;
    TF1* func = nullptr;

    double A = std::numeric_limits<double>::quiet_NaN();
    double Aerr = std::numeric_limits<double>::quiet_NaN();

    double lambda = std::numeric_limits<double>::quiet_NaN();
    double lambdaerr = std::numeric_limits<double>::quiet_NaN();

    double chi2 = std::numeric_limits<double>::quiet_NaN();
    double chi2ndf = std::numeric_limits<double>::quiet_NaN();
    double prob = std::numeric_limits<double>::quiet_NaN();
    double edm = std::numeric_limits<double>::quiet_NaN();

    int ndf = 0;
    int status = -999;
    int covstatus = -999;

    bool fit_done = false;
    bool converged = false;
};


struct PhaseKey {

    int group = 2;

    double T = 1e99;
    double h = 1e99;

    string name;
};


// ============================================================================
// Generic helpers
// ============================================================================

string trim_copy(const string& s)
{
    const size_t first =
        s.find_first_not_of(" \t\r\n");

    if (first == string::npos)
        return "";

    const size_t last =
        s.find_last_not_of(" \t\r\n");

    return s.substr(
        first,
        last-first+1
    );
}


vector<string> split_csv(const string& line)
{
    vector<string> fields;

    std::stringstream ss(line);

    string field;

    while (std::getline(ss,field,',')) {

        fields.push_back(
            trim_copy(field)
        );
    }

    return fields;
}


int find_col(
    const vector<string>& header,
    const string& name
)
{
    for (int i=0;
         i<(int)header.size();
         ++i) {

        if (trim_copy(header[i]) == name)
            return i;
    }

    return -1;
}


bool finite_number(double x)
{
    return std::isfinite(x);
}


// ============================================================================
// Phase handling
// ============================================================================

PhaseKey phase_key(const string& phase)
{
    PhaseKey k;

    k.name = phase;


    // Accept both possible names for the first phase
    if (phase == "bef_ann" ||
        phase == "before_annealing") {

        k.group = 0;

        k.T = -1e99;
        k.h = -1e99;

        return k;
    }


    double T = 0.0;
    double h = 0.0;


    if (std::sscanf(
            phase.c_str(),
            "annealing_T=%lf_h=%lf",
            &T,
            &h
        ) == 2) {

        k.group = 1;

        k.T = T;
        k.h = h;
    }


    return k;
}


bool phase_less(
    const string& a,
    const string& b
)
{
    const PhaseKey ka =
        phase_key(a);

    const PhaseKey kb =
        phase_key(b);


    if (ka.group != kb.group)
        return ka.group < kb.group;

    if (ka.T != kb.T)
        return ka.T < kb.T;

    if (ka.h != kb.h)
        return ka.h < kb.h;

    return ka.name < kb.name;
}


string pretty_phase(const string& phase)
{
    if (phase == "bef_ann" ||
        phase == "before_annealing")
        return "Before annealing";


    double T = 0.0;
    double h = 0.0;


    if (std::sscanf(
            phase.c_str(),
            "annealing_T=%lf_h=%lf",
            &T,
            &h
        ) == 2) {

        return Form(
            "%g #circC, %g h",
            T,
            h
        );
    }


    return phase;
}


// ============================================================================
// Plot style
// ============================================================================

int phase_color(size_t i)
{
    static const int colors[] = {

        kBlack,
        kRed+1,
        kOrange+7,
        kGreen+2,
        kCyan+2,
        kAzure+1,
        kBlue+1,
        kViolet+1,
        kMagenta+1,
        kGray+2,
        kPink+7,
        kTeal+3
    };


    return colors[
        i %
        (sizeof(colors)/sizeof(colors[0]))
    ];
}


int phase_marker(size_t i)
{
    static const int markers[] = {

        20,21,22,23,24,25,
        26,27,28,30,33,34
    };


    return markers[
        i %
        (sizeof(markers)/sizeof(markers[0]))
    ];
}


// ============================================================================
// Read IV CSV
// ============================================================================

vector<IVRow> read_IV_csv(
    const string& filename
)
{
    vector<IVRow> rows;

    std::ifstream fin(filename);

    if (!fin.is_open()) {

        std::cerr
            << "ERROR: cannot open CSV: "
            << filename
            << std::endl;

        return rows;
    }

    string line;

    if (!std::getline(fin,line)) {

        std::cerr
            << "ERROR: empty CSV: "
            << filename
            << std::endl;

        return rows;
    }

    const vector<string> header =
        split_csv(line);

    const int c_phase =
        find_col(header,"phase");

    const int c_T =
        find_col(header,"temperature_C");

    const int c_I =
        find_col(header,"current");

    const int c_Ierr =
        find_col(header,"current_error");


    if (c_phase < 0 ||
        c_T < 0 ||
        c_I < 0 ||
        c_Ierr < 0) {

        std::cerr
            << "ERROR: required columns not found."
            << std::endl;

        return rows;
    }

    int line_number = 1;

    while (std::getline(fin,line)) {

        ++line_number;

        if (trim_copy(line).empty())
            continue;

        vector<string> f =
            split_csv(line);

        const int max_col =
            std::max(
                std::max(c_phase,c_T),
                std::max(c_I,c_Ierr)
            );

        if ((int)f.size() <= max_col) {

            std::cerr
                << "WARNING: malformed row "
                << line_number
                << std::endl;

            continue;
        }

        try {

            IVRow r;

            r.phase =
                trim_copy(
                    f[c_phase]
                );

            r.temperature_C =
                std::stod(
                    f[c_T]
                );

            r.current =
                std::fabs(
                    std::stod(
                        f[c_I]
                    )
                );

            r.current_error =
                std::fabs(
                    std::stod(
                        f[c_Ierr]
                    )
                );

            if (r.phase.empty())
                continue;

            if (!finite_number(r.temperature_C) ||
                !finite_number(r.current) ||
                !finite_number(r.current_error))
                continue;

            rows.push_back(r);
        }

        catch (...) {

            std::cerr
                << "WARNING: cannot parse row "
                << line_number
                << std::endl;
        }
    }

    return rows;
}


// ============================================================================
// LOGARITHMIC PREFIT
// ============================================================================

void estimate_exp_parameters(
    const vector<IVRow>& rows,
    double& A0,
    double& lambda0
)
{
    vector<IVRow> positive;

    for (const auto& r : rows) {

        if (r.current > 0.0 &&
            finite_number(r.current) &&
            finite_number(r.temperature_C)) {

            positive.push_back(r);
        }
    }

    std::sort(
        positive.begin(),
        positive.end(),
        IVRow_less
    );

    if (positive.size() >= 2 &&
        std::fabs(
            positive.back().temperature_C -
            positive.front().temperature_C
        ) > 1e-12) {

        TGraph gr_log;

        for (size_t i=0;
             i<positive.size();
             ++i) {

            gr_log.SetPoint(
                i,
                positive[i].temperature_C,
                std::log(
                    positive[i].current
                )
            );
        }

        static int prefit_counter = 0;

        const string fname =
            Form(
                "f_log_prefit_%d",
                prefit_counter++
            );

        TF1 f_lin(
            fname.c_str(),
            "pol1",
            positive.front().temperature_C,
            positive.back().temperature_C
        );

        TFitResultPtr fit_lin =
            gr_log.Fit(
                &f_lin,
                "Q0SN"
            );

        const int status =
            (int)fit_lin;

        if (status == 0 &&
            finite_number(
                f_lin.GetParameter(0)
            ) &&
            finite_number(
                f_lin.GetParameter(1)
            )) {

            lambda0 =
                f_lin.GetParameter(1);

            A0 =
                std::exp(
                    f_lin.GetParameter(0)
                );

            return;
        }

        lambda0 =
            (
                std::log(
                    positive.back().current
                )
                -
                std::log(
                    positive.front().current
                )
            )
            /
            (
                positive.back().temperature_C
                -
                positive.front().temperature_C
            );

        A0 =
            std::exp(
                std::log(
                    positive.front().current
                )
                -
                lambda0 *
                positive.front().temperature_C
            );

        return;
    }

    if (positive.size() == 1) {

        A0 =
            positive[0].current;

        lambda0 =
            0.0;

        return;
    }

    A0 =
        1e-9;

    lambda0 =
        0.06;
}

} // namespace


// ============================================================================
// Main function
// ============================================================================

void I_vs_T_contemporary(
    const char* input_csv,
    const char* output_dir,
    const char* sensor
)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const string infile = input_csv ? input_csv : "";
    const string outdir = output_dir ? output_dir : "";
    const string sens   = sensor ? sensor : "";

    if (infile.empty() || outdir.empty() || sens.empty()) {
        std::cerr << "ERROR: empty required argument." << std::endl;
        return;
    }

    gSystem->mkdir(outdir.c_str(), true);

    vector<IVRow> rows = read_IV_csv(infile);

    if (rows.empty()) {
        std::cerr << "ERROR: no valid rows found in " << infile << std::endl;
        return;
    }

    set<string> phase_set;
    for (const auto& r : rows)
        phase_set.insert(r.phase);

    vector<string> phases(phase_set.begin(), phase_set.end());
    std::sort(phases.begin(), phases.end(), phase_less);

    std::cout << "\nSensor " << sens << ": found " << phases.size() << " phases" << std::endl;

    map<string,vector<IVRow>> data;
    for (const auto& r : rows)
        data[r.phase].push_back(r);

    for (auto& kv : data) {
       std::sort(kv.second.begin(), kv.second.end(), IVRow_less);
    }

    map<string,IVFitResult> fits;

    for (size_t ip=0; ip<phases.size(); ++ip) {

        const string& ph = phases[ip];

        IVFitResult info;
        info.phase = ph;
        info.rows = data[ph];

        const int n = (int)info.rows.size();

        if (n == 0) {
            fits[ph] = info;
            continue;
        }

        vector<double> xv(n);
        vector<double> yv(n);
        vector<double> exv(n,0.0);
        vector<double> eyv(n);

        for (int i=0; i<n; ++i) {
            xv[i]  = info.rows[i].temperature_C;
            yv[i]  = info.rows[i].current;
            eyv[i] = info.rows[i].current_error;
        }

        info.graph = new TGraphErrors(n, xv.data(), yv.data(), exv.data(), eyv.data());
        info.graph->SetName(Form("gr_%s_I_vs_T_phase_%zu", sens.c_str(), ip));

        const int color  = phase_color(ip);
        const int marker = phase_marker(ip);

        info.graph->SetMarkerColor(color);
        info.graph->SetLineColor(color);
        info.graph->SetMarkerStyle(marker);
        info.graph->SetMarkerSize(1.35);
        info.graph->SetLineWidth(2);

        if (n < 3) {
            std::cerr << "WARNING: phase " << ph << " has only " << n << " points. Fit skipped." << std::endl;
            fits[ph] = info;
            continue;
        }

        const double Tmin = info.rows.front().temperature_C;
        const double Tmax = info.rows.back().temperature_C;

        double dT = Tmax-Tmin;
        if (!(dT > 0.0)) dT = 1.0;

        const double fitmin = Tmin - 0.10*dT;
        const double fitmax = Tmax + 0.10*dT;

        double A0 = 0.0;
        double lambda0 = 0.0;

        estimate_exp_parameters(info.rows, A0, lambda0);

        std::cout << "\n" << sens << " | " << ph << std::endl;
        std::cout << "  logarithmic prefit:\n    A0      = " << A0 << "\n    lambda0 = " << lambda0 << std::endl;

        info.func = new TF1(Form("fit_%s_I_vs_T_phase_%zu", sens.c_str(), ip), "[0]*exp([1]*x)", fitmin, fitmax);
        info.func->SetParNames("A", "#lambda");
        info.func->SetParameters(A0, lambda0);
        info.func->SetLineColor(color);
        info.func->SetLineWidth(3);
        info.func->SetLineStyle(1);

        TFitResultPtr fit = info.graph->Fit(info.func, "QRS0N");

        info.fit_done = true;
        info.status   = (int)fit;
        info.A        = info.func->GetParameter(0);
        info.Aerr     = info.func->GetParError(0);
        info.lambda   = info.func->GetParameter(1);
        info.lambdaerr= info.func->GetParError(1);
        info.chi2     = info.func->GetChisquare();
        info.ndf      = info.func->GetNDF();
        info.chi2ndf  = (info.ndf > 0) ? info.chi2/info.ndf : 0.0;
        info.prob     = info.func->GetProb();

        if (fit.Get()) {
            info.covstatus = fit->CovMatrixStatus();
            info.edm       = fit->Edm();
        }

        info.converged = (info.status == 0) && finite_number(info.A) && finite_number(info.lambda);

        std::cout << "  final exponential fit:\n"
                  << "    A        = " << info.A << " +/- " << info.Aerr << "\n"
                  << "    lambda   = " << info.lambda << " +/- " << info.lambdaerr << " C^-1\n"
                  << "    chi2/ndf = " << info.chi2ndf << "\n"
                  << "    status   = " << info.status << std::endl;

        fits[ph] = info;
    }

    double xmin =  1e99, xmax = -1e99;
    double ymin =  1e99, ymax = -1e99;

    for (const auto& r : rows) {
        xmin = std::min(xmin, r.temperature_C);
        xmax = std::max(xmax, r.temperature_C);
        ymin = std::min(ymin, r.current - std::fabs(r.current_error));
        ymax = std::max(ymax, r.current + std::fabs(r.current_error));
    }

    for (const auto& kv : fits) {
        const IVFitResult& f = kv.second;
        if (!f.fit_done || !f.func || f.rows.empty()) continue;

        const double y1 = f.func->Eval(f.func->GetXmin());
        const double y2 = f.func->Eval(f.func->GetXmax());

        if (finite_number(y1)) { ymin = std::min(ymin, y1); ymax = std::max(ymax, y1); }
        if (finite_number(y2)) { ymin = std::min(ymin, y2); ymax = std::max(ymax, y2); }
    }

    double xspan = xmax-xmin;
    if (!(xspan > 0.0)) xspan = 1.0;

    const double frame_xmin = xmin - 0.10*xspan;
    const double frame_xmax = xmax + 0.10*xspan;

    double yspan = ymax-ymin;
    if (!(yspan > 0.0)) yspan = std::max(1e-12, std::fabs(ymax)*0.20);

    double frame_ymin = ymin - 0.10*yspan;
    double frame_ymax = ymax + 0.55*yspan;

    if (frame_ymin < 0.0) frame_ymin = 0.0;

    TCanvas* c = new TCanvas(Form("c_%s_I_vs_T", sens.c_str()), Form("%s - Current vs temperature", sens.c_str()), 1900, 1000);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.11);
    c->SetTopMargin(0.08);

    TH1* frame = c->DrawFrame(
        frame_xmin, frame_ymin, frame_xmax, frame_ymax,
        Form("Sensor %s - Dark Current vs temperature for all annealing phases;T (#circ C);Dark Current (A)", sens.c_str())
    );

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleOffset(1.20);

    TLegend* leg = new TLegend(0.12, 0.70, 0.88, 0.91);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.020);
    leg->SetNColumns(3);

    for (size_t ip=0; ip<phases.size(); ++ip) {
        const string& ph = phases[ip];
        auto it = fits.find(ph);
        if (it == fits.end()) continue;

        IVFitResult& info = it->second;
        if (!info.graph) continue;

        info.graph->Draw("PZ SAME");

        if (info.fit_done && info.func) {
            info.func->Draw("SAME");
        }

        string label = pretty_phase(ph);

        if (info.fit_done && info.converged && finite_number(info.lambda) && finite_number(info.lambdaerr)) {
            label += Form("  (#lambda = %.5f #pm %.5f #circ C^{-1})", info.lambda, info.lambdaerr);
        } else {
            label += "  (fit unavailable)";
        }

        leg->AddEntry(info.graph, label.c_str(), "pe");
    }

    leg->AddEntry((TObject*)0, "Solid curves: I = A e^{#lambda T}", "");
    leg->AddEntry((TObject*)0, "Vertical bars: statistical uncertainty", "");
    leg->Draw();

    c->Modified();
    c->Update();

    c->SaveAs(Form("%s/%s_I_vs_T.png", outdir.c_str(), sens.c_str()));
    c->SaveAs(Form("%s/%s_I_vs_T.pdf", outdir.c_str(), sens.c_str()));

    const string result_csv = outdir + "/" + sens + "_I_vs_T_fit_results.csv";
    std::ofstream fout(result_csv);

    fout << "phase,A,A_stat_error,lambda,lambda_stat_error,chi2,ndf,chi2ndf,prob,fit_status,covmatrix_status,edm,converged\n";
    fout << std::setprecision(12);

    for (const string& ph : phases) {
        auto it = fits.find(ph);
        if (it == fits.end()) continue;

        const IVFitResult& f = it->second;

        fout << ph << ',' << f.A << ',' << f.Aerr << ',' << f.lambda << ',' << f.lambdaerr << ','
             << f.chi2 << ',' << f.ndf << ',' << f.chi2ndf << ',' << f.prob << ','
             << f.status << ',' << f.covstatus << ',' << f.edm << ',' << (f.converged ? 1 : 0) << '\n';
    }

    fout.close();

    std::cout << "\n=============================================\n"
              << "I vs T analysis completed for " << sens << "\n"
              << "Canvas saved to " << outdir << "\n"
              << "=============================================" << std::endl;
}