// ROOT usage (from the directory containing the CSVs):
// root -l 'plot_mean_luminosity_vs_phase.C("A1_T=20_v=5_total_luminosity_sum.csv","B1_T=20_v=3_total_luminosity_sum.csv","A1_B1_mean_luminosity_vs_phase")'
//
// Expected CSV columns: phase,mean_luminosity,error
// Accepts both annealing_150_125 and annealing_T=150_h=125 phase names.

#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

// ============================================================================
// Data and phase parsing
// ============================================================================

struct LuminosityRow {
    string phase;
    double mean_luminosity;
    double error;
};

struct PhaseInfo {
    bool valid = false;
    bool before = false;
    double temperature = numeric_limits<double>::quiet_NaN();
    double hours = numeric_limits<double>::quiet_NaN();
};

static string trim_copy(const string& s)
{
    const auto first = s.find_first_not_of(" \t\n\r\f\v");
    if (first == string::npos) return "";
    const auto last = s.find_last_not_of(" \t\n\r\f\v");
    return s.substr(first, last - first + 1);
}

static vector<string> split_csv_simple(const string& line)
{
    vector<string> result;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ',')) result.push_back(trim_copy(item));
    return result;
}

static PhaseInfo parse_phase(const string& input)
{
    const string name = trim_copy(input);
    PhaseInfo p;

    if (name == "before_annealing" || name == "before_ann" || name == "bef_ann") {
        p.valid = true;
        p.before = true;
        return p;
    }

    double T = 0.0, h = 0.0;
    int consumed = -1;

    // Format 1: annealing_75_5
    bool matched =
        sscanf(name.c_str(), "annealing_%lf_%lf%n", &T, &h, &consumed) == 2 &&
        consumed == (int)name.size();

    if (!matched) {
        consumed = -1;
        // Format 2: annealing_T=75_h=5
        matched =
            sscanf(name.c_str(), "annealing_T=%lf_h=%lf%n", &T, &h, &consumed) == 2 &&
            consumed == (int)name.size();
    }

    if (matched && isfinite(T) && isfinite(h) && h > 0) {
        p.valid = true;
        p.temperature = T;
        p.hours = h;
    }
    return p;
}

// Normalize both naming conventions to the same map key.
static string canonical_phase(const string& input)
{
    const PhaseInfo p = parse_phase(input);
    if (p.before) return "before_annealing";
    if (!p.valid) return trim_copy(input);
    ostringstream ss;
    ss << "annealing_" << p.temperature << "_" << p.hours;
    return ss.str();
}

// ============================================================================
// Annealing positions
// ============================================================================

static vector<double> build_phase_positions(const vector<string>& phases)
{
    vector<double> x(phases.size(), 0.0);

    for (size_t i = 1; i < phases.size(); ++i) {
        const PhaseInfo current  = parse_phase(phases[i]);
        const PhaseInfo previous = parse_phase(phases[i - 1]);
        if (!current.valid || current.before)
            throw runtime_error("Unknown annealing phase: " + phases[i]);

        // Keep the original positions of all the existing phases:
        // 0, 1, 6, 7, 12, 13, 18, 19, 24.
        double step = current.hours / 5.0;

        // Extra-long annealing at the SAME temperature:
        // 150 C, 25 h -> 150 C, 125 h: delta_t = 100 h = 20 units.
        // Thus the new position is x = 24 + 20 = 44.
        if (previous.valid && !previous.before &&
            fabs(current.temperature - previous.temperature) < 1e-9 &&
            current.hours > 25.0 && current.hours > previous.hours) {
            step = (current.hours - previous.hours) / 5.0;
        }

        if (!(step > 0.0) || !isfinite(step))
            throw runtime_error("Invalid phase spacing: " + phases[i]);
        x[i] = x[i - 1] + step;
    }
    return x;
}

// ============================================================================
// CSV reading (independently for A1 and B1)
// ============================================================================

static map<string, pair<double, double>> read_luminosity_csv(const char* path)
{
    ifstream fin(path);
    if (!fin.is_open()) throw runtime_error(string("Cannot open CSV file: ") + path);

    string line;
    if (!getline(fin, line)) throw runtime_error(string("Empty CSV file: ") + path);

    const vector<string> header = split_csv_simple(line);
    map<string, size_t> col;
    for (size_t i = 0; i < header.size(); ++i) col[header[i]] = i;

    for (const string& required : {"phase", "mean_luminosity", "error"}) {
        if (!col.count(required))
            throw runtime_error(string("Missing column '") + required + "' in " + path);
    }

    const size_t max_col = max({col["phase"], col["mean_luminosity"], col["error"]});
    map<string, pair<double, double>> data;

    while (getline(fin, line)) {
        if (trim_copy(line).empty()) continue;
        const vector<string> fields = split_csv_simple(line);
        if (fields.size() <= max_col) {
            cerr << "WARNING: incomplete row in " << path << ": " << line << endl;
            continue;
        }
        try {
            const string phase = canonical_phase(fields[col["phase"]]);
            const double y = stod(fields[col["mean_luminosity"]]);
            const double ey = stod(fields[col["error"]]);
            if (!isfinite(y) || !isfinite(ey) || ey < 0.0) continue;
            if (data.count(phase))
                cerr << "WARNING: repeated phase " << phase << " in " << path
                     << " (keeping last value)" << endl;
            data[phase] = {y, ey};
        } catch (const exception&) {
            cerr << "WARNING: invalid row in " << path << ": " << line << endl;
        }
    }
    if (data.empty()) throw runtime_error(string("No valid data in ") + path);
    return data;
}

// Build one graph. Missing phases remain gaps in the observations;
// positions are ALWAYS shared between the two sensors.
static TGraphErrors* make_sensor_graph(
    const string& name,
    const map<string, pair<double, double>>& data,
    const vector<string>& phases,
    const vector<double>& x,
    int color,
    int marker)
{
    vector<double> xv, yv, ex, ey;
    for (size_t i = 0; i < phases.size(); ++i) {
        const auto it = data.find(phases[i]);
        if (it == data.end()) {
            cerr << "WARNING: " << name << " missing " << phases[i] << endl;
            continue;
        }
        xv.push_back(x[i]);
        yv.push_back(it->second.first);
        ex.push_back(0.0);
        ey.push_back(it->second.second);
    }

    if (xv.empty()) throw runtime_error("No plotted values for " + name);
    auto* gr = new TGraphErrors((int)xv.size(), xv.data(), yv.data(),
                                ex.data(), ey.data());
    gr->SetName(("gr_mean_luminosity_" + name).c_str());
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(1.45);
    gr->SetLineWidth(3);
    return gr;
}

// ============================================================================
// Bottom time labels and top annealing temperature labels
// ============================================================================

static void draw_time_labels(TPad* pad, const vector<string>& phases,
                             const vector<double>& x, double xmin, double xmax)
{
    if (!pad || phases.empty() || !(xmax > xmin)) return;
    pad->Update();
    const double left = pad->GetLeftMargin();
    const double usable = 1.0 - left - pad->GetRightMargin();

    for (size_t i = 0; i < phases.size(); ++i) {
        const PhaseInfo p = parse_phase(phases[i]);
        if (p.before || !p.valid) continue;

        double x_ndc = left + (x[i] - xmin) / (xmax - xmin) * usable;
        if (i == 2 || i == 4 || i == 6 || i == 9){
            x_ndc = x_ndc - 0.01;
        } else if (i == 3 || i == 5 || i == 7 ){
            x_ndc = x_ndc + 0.01;
        }

        const string label = Form("%g h", p.hours);
        auto* lab = new TLatex();
        lab->SetNDC();
        lab->SetTextFont(42);
        lab->SetTextSize(0.032);
        lab->SetTextAlign(22);
        lab->DrawLatex(x_ndc, 0.13, label.c_str());
    }

    auto* title = new TLatex();
    title->SetNDC();
    title->SetTextFont(42);
    title->SetTextSize(0.047);
    title->SetTextAlign(22);
    title->DrawLatex(0.5, 0.08, "Annealing time");
}

static void draw_temperature_blocks(const vector<string>& phases,
                                    const vector<double>& x,
                                    double ymin, double ymax)
{
    map<double, vector<double>> by_temperature;
    for (size_t i = 0; i < phases.size(); ++i) {
        const PhaseInfo p = parse_phase(phases[i]);
        if (p.valid && !p.before) by_temperature[p.temperature].push_back(x[i]);
    }

    const double ytext = ymax - 0.055 * (ymax - ymin);
    TLatex latex;
    latex.SetTextFont(42);
    latex.SetTextAlign(22);
    latex.SetTextSize(0.030);
    latex.DrawLatex(x.front() - 1.0, ytext, "#splitline{Before}{annealing}");
    latex.SetTextSize(0.034);

    for (const auto& kv : by_temperature) {
        const auto mm = minmax_element(kv.second.begin(), kv.second.end());
        const double centre = 0.5 * (*mm.first + *mm.second);
        const string label = Form("T = %g#circC", kv.first);
        latex.DrawLatex(centre, ytext, label.c_str());
    }

    // Draw separators at temperature changes, using actual x positions.
    for (size_t i = 1; i < phases.size(); ++i) {
        const PhaseInfo prev = parse_phase(phases[i - 1]);
        const PhaseInfo next = parse_phase(phases[i]);
        if (next.before || !next.valid) continue;
        if (prev.before || fabs(prev.temperature - next.temperature) > 1e-9) {
            const double separator = 0.5 * (x[i - 1] + x[i]);
            auto* line = new TLine(separator, ymin, separator, ymax);
            line->SetLineColor(kGray + 1);
            line->SetLineStyle(2);
            line->SetLineWidth(1);
            line->Draw("same");
        }
    }
}

// ============================================================================
// Main macro: A1 red, B1 yellow, BOTH on the same canvas
// ============================================================================

void plot_mean_luminosity_vs_phase(
    const char* csv_A1 = "A1_T=20_v=5_total_luminosity_sum.csv",
    const char* csv_B1 = "B1_T=20_v=3_total_luminosity_sum.csv",
    const char* output_prefix = "A1_B1_mean_luminosity_vs_phase")
{
    gStyle->SetOptStat(0);

    const vector<string> phases = {
        "before_annealing",
        "annealing_75_5",  "annealing_75_25",
        "annealing_100_5", "annealing_100_25",
        "annealing_125_5", "annealing_125_25",
        "annealing_150_5", "annealing_150_25",
        "annealing_150_125"
    };
    const vector<double> x = build_phase_positions(phases);

    const auto data_A1 = read_luminosity_csv(csv_A1);
    const auto data_B1 = read_luminosity_csv(csv_B1);

    // Range of both sensors together (including measurement errors).
    double ymin_data = numeric_limits<double>::infinity();
    double ymax_data = -numeric_limits<double>::infinity();
    for (const auto& data : {data_A1, data_B1}) {
        for (const auto& phase : phases) {
            const auto it = data.find(phase);
            if (it == data.end()) continue;
            ymin_data = min(ymin_data, it->second.first - it->second.second);
            ymax_data = max(ymax_data, it->second.first + it->second.second);
        }
    }
    if (!isfinite(ymin_data) || !isfinite(ymax_data))
        throw runtime_error("No valid phases found in either CSV.");

    double yrange = ymax_data - ymin_data;
    if (!(yrange > 0)) yrange = max(1.0, fabs(ymax_data) * 0.10);
    const double ymin = max(0.0, ymin_data - 0.12 * yrange);
    const double ymax = ymax_data + 0.25 * yrange;
    const double xmin = -4.5;
    const double xmax = x.back() + 0.5;
    const int nbins = (int)llround(x.back()) + 1;

    auto* c = new TCanvas("c_mean_luminosity_vs_phase",
                          "mean luminosity vs annealing phase", 2250, 1000);
    c->SetBottomMargin(0.16);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.10);
    c->SetTicks(1, 1);

    const string plot_title =
        "Sensors A1 and B1 - Mean luminosity vs annealing phase, T = 20#circC";
    auto* frame = new TH1D("frame_mean_luminosity_vs_phase",
                           Form("%s;;Mean luminosity (a.u)", plot_title.c_str()),
                           nbins, xmin, xmax);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->SetTitleSize(0.045);
    frame->SetTitleOffset(1.0);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetNdivisions(0);
    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.047);
    frame->GetYaxis()->SetTitleOffset(1.05);
    frame->Draw();

    const double tick_height = 0.018 * (ymax - ymin);
    for (double xi : x) {
        auto* tick = new TLine(xi, ymin, xi, ymin + tick_height);
        tick->SetLineColor(kBlack);
        tick->SetLineWidth(1);
        tick->Draw("same");
    }

    draw_temperature_blocks(phases, x, ymin, ymax);

    const int red_A1    = kP6Red;
    const int yellow_B1 = kP6Yellow;
    auto* gr_A1 = make_sensor_graph("A1", data_A1, phases, x, red_A1, 20);
    auto* gr_B1 = make_sensor_graph("B1", data_B1, phases, x, yellow_B1, 22);
    gr_A1->Draw("PL E1 SAME");
    gr_B1->Draw("PL E1 SAME");

    auto* leg = new TLegend(0.69, 0.58, 0.88, 0.71);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->AddEntry(gr_A1, "A1 (v_{over} = 5 V)", "lep");
    leg->AddEntry(gr_B1, "B1 (v_{over} = 3 V)", "lep");
    leg->Draw();

    c->Modified();
    c->Update();
    draw_time_labels(c, phases, x,
                     frame->GetXaxis()->GetXmin(),
                     frame->GetXaxis()->GetXmax());
    gPad->RedrawAxis();
    c->Modified();
    c->Update();

    c->SaveAs(Form("%s.png", output_prefix));
    c->SaveAs(Form("%s.pdf", output_prefix));

    cout << "\nSensors: A1 (20 C, v_over=5 V) and B1 (20 C, v_over=3 V)\n";
    cout << "CSV A1: " << csv_A1 << "\nCSV B1: " << csv_B1 << "\n";
    cout << "Phase positions (5 h = 1 unit):\n";
    for (size_t i = 0; i < phases.size(); ++i)
        cout << "  " << phases[i] << " -> x = " << x[i] << '\n';
    cout << "150 C: 25 h -> 125 h = +100 h = +20 units\n";
    cout << "Saved: " << output_prefix << ".png and "
         << output_prefix << ".pdf\n\n";
}
