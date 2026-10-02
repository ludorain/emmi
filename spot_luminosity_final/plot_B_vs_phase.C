// Input summary CSV expected columns:
//   phase,B_mean,B_mean_stat_error,B_rms,n_spots
//
// Example:
// root -l 'plot_B_vs_phase.C("analysis/A1_T=20/A1_T=20_B_vs_phase.csv","A1","A1_B_vs_annealing")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <limits>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TString.h"
#include "TLatex.h"
#include "TPad.h"

using namespace std;

static vector<string> split_csv_B_phase(const string& line)
{
    vector<string> fields;
    string field;
    stringstream ss(line);
    while (getline(ss, field, ',')) fields.push_back(field);
    return fields;
}

static bool finite_value(double x) { return std::isfinite(x); }

static string dirname_of(const string& path)
{
    const size_t pos = path.find_last_of("/\\");
    if (pos == string::npos) return ".";
    if (pos == 0) return "/";
    return path.substr(0, pos);
}

static string basename_of(const string& path)
{
    const size_t pos = path.find_last_of("/\\");
    if (pos == string::npos) return path;
    return path.substr(pos + 1);
}

static string condition_from_summary_filename(const string& csv_file)
{
    string name = basename_of(csv_file);
    const string suffix = "_B_vs_phase.csv";
    if (name.size() >= suffix.size() &&
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) {
        return name.substr(0, name.size() - suffix.size());
    }
    return "";
}

struct PhaseInfo {
    string raw;
    bool before = false;
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double hours = std::numeric_limits<double>::quiet_NaN();
};

static PhaseInfo parse_phase(const string& phase)
{
    PhaseInfo p;
    p.raw = phase;
    if (phase == "before_annealing") {
        p.before = true;
        return p;
    }

    const string prefix = "annealing_T=";
    if (phase.find(prefix) != 0) return p;

    const size_t pos_h = phase.find("_h=");
    if (pos_h == string::npos) return p;

    try {
        const string Tstr = phase.substr(prefix.size(), pos_h - prefix.size());
        const string hstr = phase.substr(pos_h + 3);
        p.temperature = stod(Tstr);
        p.hours = stod(hstr);
    } catch (...) {
        p.temperature = std::numeric_limits<double>::quiet_NaN();
        p.hours = std::numeric_limits<double>::quiet_NaN();
    }
    return p;
}

static double phase_step_units(const string& phase)
{
    const PhaseInfo p = parse_phase(phase);
    if (p.before) return 0.0;
    if (finite_value(p.hours) && p.hours > 0.0) return p.hours / 5.0;
    return 1.0;
}

static vector<double> build_phase_positions(const vector<string>& phases)
{
    vector<double> x(phases.size(), 0.0);
    if (phases.empty()) return x;

    x[0] = 0.0;
    for (size_t i = 1; i < phases.size(); ++i) {
        double step = phase_step_units(phases[i]);
        if (!(step > 0.0) || !finite_value(step)) step = 1.0;
        x[i] = x[i - 1] + step;
    }
    return x;
}

static string two_line_phase_label(const string& phase)
{
    const PhaseInfo p = parse_phase(phase);
    if (p.before) return "#splitline{Before}{annealing}";

    if (finite_value(p.temperature) && finite_value(p.hours)) {
        const double Tround = std::round(p.temperature);
        const double hround = std::round(p.hours);

        string Ttext;
        string htext;

        if (std::fabs(p.temperature - Tround) < 1e-9)
            Ttext = Form("%.0f^{#circ}C", p.temperature);
        else
            Ttext = Form("%.1f^{#circ}C", p.temperature);

        if (std::fabs(p.hours - hround) < 1e-9)
            htext = Form("%.0f h", p.hours);
        else
            htext = Form("%.1f h", p.hours);

        return "#splitline{" + Ttext + "}{" + htext + "}";
    }

    return phase;
}

static TH1D* make_phase_frame(const char* name,
                              const char* title,
                              const vector<string>& phases,
                              const vector<double>& xpos,
                              double ymin,
                              double ymax)
{
    if (phases.empty() || xpos.empty()) return nullptr;

    const int max_unit = std::max(0, (int)std::llround(xpos.back()));
    const int nbins = max_unit + 1;

    TH1D* frame = new TH1D(name, title, nbins, -0.5, max_unit + 0.5);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);

    // Hide ROOT's categorical bin labels. They are drawn manually below
    // with TLatex, which guarantees that they remain horizontal.
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleOffset(1.85);

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.15);

    return frame;
}


// Draw phase labels manually in NDC coordinates. This avoids ROOT's
// automatic rotation of alphanumeric bin labels and keeps every label
// horizontal, including the two-line temperature/time labels.
static void draw_horizontal_phase_labels(TPad* pad,
                                         const vector<string>& phases,
                                         const vector<double>& xpos,
                                         double xmin,
                                         double xmax)
{
    if (!pad || phases.empty() || xpos.empty() || !(xmax > xmin)) return;

    pad->Update();

    const double left  = pad->GetLeftMargin();
    const double right = pad->GetRightMargin();
    const double usable = 1.0 - left - right;

    for (size_t i = 0; i < phases.size(); ++i) {
        const double frac = (xpos[i] - xmin) / (xmax - xmin);
        const double x_ndc = left + frac * usable;

        TLatex* lab = new TLatex();
        lab->SetNDC();
        lab->SetTextFont(42);
        lab->SetTextSize(0.031);
        lab->SetTextAlign(22);
        lab->DrawLatex(x_ndc, 0.135,
                       two_line_phase_label(phases[i]).c_str());
    }
}

struct PhaseSummary {
    string phase;
    double B_mean = std::numeric_limits<double>::quiet_NaN();
    double B_mean_stat_error = std::numeric_limits<double>::quiet_NaN();
    double B_rms = std::numeric_limits<double>::quiet_NaN();
    int n_spots = 0;
};

static bool read_B_summary(const string& csv_file,
                                vector<PhaseSummary>& data)
{
    ifstream fin(csv_file);
    if (!fin.is_open()) {
        cerr << "ERROR: cannot open summary CSV: " << csv_file << endl;
        return false;
    }

    string line;
    if (!getline(fin, line)) {
        cerr << "ERROR: empty summary CSV: " << csv_file << endl;
        return false;
    }

    const vector<string> header = split_csv_B_phase(line);
    map<string, int> col;
    for (int i = 0; i < (int)header.size(); ++i) col[header[i]] = i;

    const vector<string> required = {
        "phase", "B_mean", "B_mean_stat_error", "B_rms", "n_spots"
    };

    for (const auto& name : required) {
        if (!col.count(name)) {
            cerr << "ERROR: summary CSV is missing column '" << name << "'." << endl;
            return false;
        }
    }

    while (getline(fin, line)) {
        if (line.empty()) continue;

        vector<string> fields = split_csv_B_phase(line);
        if (fields.size() < header.size()) fields.resize(header.size(), "");

        try {
            PhaseSummary row;
            row.phase = fields[col["phase"]];
            row.B_mean = stod(fields[col["B_mean"]]);
            row.B_mean_stat_error = stod(fields[col["B_mean_stat_error"]]);
            row.B_rms = stod(fields[col["B_rms"]]);
            row.n_spots = stoi(fields[col["n_spots"]]);

            if (!finite_value(row.B_mean) ||
                !finite_value(row.B_mean_stat_error) ||
                !finite_value(row.B_rms)) continue;

            data.push_back(row);
        } catch (...) {
            cerr << "WARNING: malformed summary row skipped:\n" << line << endl;
        }
    }

    fin.close();

    if (data.empty()) {
        cerr << "ERROR: no valid rows found in summary CSV." << endl;
        return false;
    }

    return true;
}

static bool read_individual_B_values(const string& csv_file,
                                     vector<double>& B_values)
{
    ifstream fin(csv_file);
    if (!fin.is_open()) {
        cerr << "WARNING: cannot open individual-fit CSV: " << csv_file << endl;
        return false;
    }

    string line;
    if (!getline(fin, line)) return false;

    const vector<string> header = split_csv_B_phase(line);
    map<string, int> col;
    for (int i = 0; i < (int)header.size(); ++i) col[header[i]] = i;

    const vector<string> required = {
        "B", "B_stat_error", "ndf", "chi2ndf", "converged"
    };

    for (const auto& name : required) {
        if (!col.count(name)) {
            cerr << "WARNING: " << csv_file
                 << " is missing column '" << name << "'." << endl;
            return false;
        }
    }

    while (getline(fin, line)) {
        if (line.empty()) continue;

        vector<string> fields = split_csv_B_phase(line);
        if (fields.size() < header.size()) fields.resize(header.size(), "");

        try {
            const double B = stod(fields[col["B"]]);
            const double B_err = stod(fields[col["B_stat_error"]]);
            const int ndf = stoi(fields[col["ndf"]]);
            const double chi2ndf = stod(fields[col["chi2ndf"]]);
            const int converged = stoi(fields[col["converged"]]);

            if (!converged ||
                !finite_value(B) ||
                !finite_value(B_err) ||
                ndf <= 0 ||
                chi2ndf >= 3.5) {
                continue;
            }

            B_values.push_back(B);
        } catch (...) {
            // Skip malformed rows.
        }
    }

    fin.close();
    return !B_values.empty();
}

void plot_B_vs_phase(
    const char* csv_file,
    const char* sensor,
    const char* output_name = "B_vs_phase")
{
    gStyle->SetOptStat(0);

    vector<PhaseSummary> summary;
    if (!read_B_summary(csv_file, summary)) return;

    vector<string> phases;
    vector<double> B_mean;
    vector<double> B_mean_err;
    vector<double> B_rms;
    vector<int> n_spots;

    for (const auto& row : summary) {
        phases.push_back(row.phase);
        B_mean.push_back(row.B_mean);
        B_mean_err.push_back(row.B_mean_stat_error);
        B_rms.push_back(row.B_rms);
        n_spots.push_back(row.n_spots);
    }

    const int n = (int)phases.size();

    // Non-uniform annealing-phase positions:
    // 5 h -> +1, 25 h -> +5.
    vector<double> x = build_phase_positions(phases);
    vector<double> ex(n, 0.0);

    // =====================================================================
    // CANVAS 1: mean B vs annealing phase
    // =====================================================================

    double ymin1 = 1e99;
    double ymax1 = -1e99;

    for (int i = 0; i < n; ++i) {
        const double extent = std::max(B_mean_err[i], B_rms[i]);
        ymin1 = std::min(ymin1, B_mean[i] - extent);
        ymax1 = std::max(ymax1, B_mean[i] + extent);
    }

    double yrange1 = ymax1 - ymin1;
    if (!(yrange1 > 0.0)) yrange1 = std::max(0.01, std::fabs(ymax1) * 0.10);
    ymin1 -= 0.25 * yrange1;
    ymax1 += 0.35 * yrange1;

    TCanvas* c1 = new TCanvas(
        "c_B_mean_vs_phase",
        "Mean B vs annealing phase",
        2100,
        950
    );

    c1->SetBottomMargin(0.20);
    c1->SetLeftMargin(0.10);
    c1->SetRightMargin(0.04);
    c1->SetTopMargin(0.10);

    string title1 = Form(
        "Sensor %s - Mean fit parameter B vs annealing phases;Annealing phase;Mean B",
        sensor
    );

    TH1D* frame1 = make_phase_frame(
        "frame_B_mean_phase",
        title1.c_str(),
        phases,
        x,
        ymin1,
        ymax1
    );

    if (!frame1) return;
    frame1->Draw();

    // RMS of the hotspot-to-hotspot B distribution. Draw it with
    // ROOT's bracket option, exactly as a systematic uncertainty.
    TGraphErrors* gr_rms = new TGraphErrors(
        n,
        x.data(),
        B_mean.data(),
        ex.data(),
        B_rms.data()
    );

    gr_rms->SetMarkerSize(0);
    gr_rms->SetLineColor(kBlack);
    gr_rms->SetLineWidth(3);
    gr_rms->Draw("[] SAME");

    // Statistical uncertainty on the mean: vertical bar only.
    TGraphErrors* gr_mean = new TGraphErrors(
        n,
        x.data(),
        B_mean.data(),
        ex.data(),
        B_mean_err.data()
    );

    gr_mean->SetMarkerStyle(20);
    gr_mean->SetMarkerSize(1.5);
    gr_mean->SetLineWidth(2);
    gr_mean->Draw("PZ SAME");

    TLegend* leg1 = new TLegend(0.13, 0.70, 0.48, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.034);
    leg1->AddEntry(gr_mean, "Mean B", "p");
    leg1->AddEntry((TObject*)nullptr,
                   "Vertical bars: statistical uncertainty on mean",
                   "");
    leg1->AddEntry((TObject*)nullptr,
                   "Brackets: RMS",
                   "");
    leg1->Draw();

    c1->Update();
    draw_horizontal_phase_labels(
        c1, phases, x,
        frame1->GetXaxis()->GetXmin(),
        frame1->GetXaxis()->GetXmax()
    );
    c1->Update();
    c1->SaveAs(Form("%s.png", output_name));
    c1->SaveAs(Form("%s.pdf", output_name));

    // =====================================================================
    // CANVAS 2: all individual B values + highlighted phase mean
    // =====================================================================

    const string summary_path = csv_file;
    const string analysis_base = dirname_of(summary_path);
    const string condition = condition_from_summary_filename(summary_path);

    if (condition.empty()) {
        cerr << "ERROR: cannot infer the analysis condition from summary filename.\n"
             << "Expected a filename like A1_T=20_B_vs_phase.csv" << endl;
        return;
    }

    vector<double> all_x;
    vector<double> all_B;
    vector<vector<double>> phase_B_values(n);

    for (int i = 0; i < n; ++i) {
        const string individual_csv =
            analysis_base + "/" + phases[i] + "/" +
            condition + "_" + phases[i] + "_B_fit_results.csv";

        vector<double> vals;
        if (!read_individual_B_values(individual_csv, vals)) {
            cerr << "WARNING: no individual B values found for phase "
                 << phases[i] << endl;
            continue;
        }

        phase_B_values[i] = vals;

        if ((int)vals.size() != n_spots[i]) {
            cerr << "WARNING: phase " << phases[i]
                 << " contains " << vals.size()
                 << " selected individual B values, while summary CSV reports n_spots="
                 << n_spots[i] << "." << endl;
        }

        const double jitter_half_width = 0.28;

        for (size_t j = 0; j < vals.size(); ++j) {
            const double u = std::fmod((double)(j + 1) * 0.6180339887498949, 1.0);
            const double jitter = (2.0 * u - 1.0) * jitter_half_width;

            all_x.push_back(x[i] + jitter);
            all_B.push_back(vals[j]);
        }
    }

    if (all_B.empty()) {
        cerr << "WARNING: second plot was not created because no individual B values could be read." << endl;
        return;
    }

    double ymin2 = 1e99;
    double ymax2 = -1e99;

    for (double value : all_B) {
        ymin2 = std::min(ymin2, value);
        ymax2 = std::max(ymax2, value);
    }

    for (int i = 0; i < n; ++i) {
        const double extent = std::max(B_mean_err[i], B_rms[i]);
        ymin2 = std::min(ymin2, B_mean[i] - extent);
        ymax2 = std::max(ymax2, B_mean[i] + extent);
    }

    double yrange2 = ymax2 - ymin2;
    if (!(yrange2 > 0.0)) yrange2 = std::max(0.01, std::fabs(ymax2) * 0.10);
    ymin2 -= 0.08 * yrange2;
    ymax2 += 0.15 * yrange2;

    TCanvas* c2 = new TCanvas(
        "c_B_all_points_vs_phase",
        "Individual B values vs annealing phase",
        2200,
        1000
    );

    c2->SetBottomMargin(0.20);
    c2->SetLeftMargin(0.10);
    c2->SetRightMargin(0.04);
    c2->SetTopMargin(0.10);

    string title2 = Form(
        "Sensor %s - Individual fit parameter B values vs annealing phases;Annealing phase;B",
        sensor
    );

    TH1D* frame2 = make_phase_frame(
        "frame_B_all_points_phase",
        title2.c_str(),
        phases,
        x,
        ymin2,
        ymax2
    );

    if (!frame2) return;
    frame2->Draw();

    TGraph* gr_individual = new TGraph(
        (int)all_B.size(),
        all_x.data(),
        all_B.data()
    );

    gr_individual->SetMarkerStyle(20);
    gr_individual->SetMarkerSize(0.75);
    gr_individual->SetMarkerColorAlpha(kP6Blue, 0.45);
    gr_individual->Draw("P SAME");

    // RMS around the phase mean, shown as ROOT brackets.
    TGraphErrors* gr_rms_overlay = new TGraphErrors(
        n,
        x.data(),
        B_mean.data(),
        ex.data(),
        B_rms.data()
    );

    gr_rms_overlay->SetMarkerSize(0);
    gr_rms_overlay->SetLineColor(kBlack);
    gr_rms_overlay->SetLineWidth(3);
    gr_rms_overlay->Draw("[] SAME");

    // Highlighted mean with statistical uncertainty only.
    TGraphErrors* gr_mean_overlay = new TGraphErrors(
        n,
        x.data(),
        B_mean.data(),
        ex.data(),
        B_mean_err.data()
    );

    gr_mean_overlay->SetMarkerStyle(21);
    gr_mean_overlay->SetMarkerSize(1.7);
    gr_mean_overlay->SetMarkerColor(kBlack);
    gr_mean_overlay->SetLineColor(kBlack);
    gr_mean_overlay->SetLineWidth(3);
    gr_mean_overlay->Draw("PZ SAME");

    TLegend* leg2 = new TLegend(0.13, 0.69, 0.50, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.033);
    leg2->AddEntry(gr_individual, "Individual hotspot B", "p");
    leg2->AddEntry(gr_mean_overlay, "Phase mean B", "p");
    leg2->AddEntry((TObject*)nullptr,
                   "Vertical bars: statistical uncertainty on mean",
                   "");
    leg2->AddEntry((TObject*)nullptr,
                   "Brackets: RMS",
                   "");
    leg2->Draw();

    c2->Update();
    draw_horizontal_phase_labels(
        c2, phases, x,
        frame2->GetXaxis()->GetXmin(),
        frame2->GetXaxis()->GetXmax()
    );
    c2->Update();
    c2->SaveAs(Form("%s_all_points.png", output_name));
    c2->SaveAs(Form("%s_all_points.pdf", output_name));

    cout << "Saved plots:" << endl;
    cout << "  " << output_name << ".png" << endl;
    cout << "  " << output_name << ".pdf" << endl;
    cout << "  " << output_name << "_all_points.png" << endl;
    cout << "  " << output_name << "_all_points.pdf" << endl;

    cout << "Phase x positions (5 h = 1 unit, 25 h = 5 units):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "  " << phases[i] << " -> x = " << x[i] << endl;
    }
}
