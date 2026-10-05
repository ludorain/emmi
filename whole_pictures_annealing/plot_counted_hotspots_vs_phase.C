// root -l 'plot_counted_hotspots_vs_phase.C("hotspot_phase_summary_gold_error.csv", "hotspots_gold")'

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
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
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct SummaryRow {
    string sensor;
    string phase;
    double counted_hotspots = numeric_limits<double>::quiet_NaN();
    double counted_hotspots_error = numeric_limits<double>::quiet_NaN();
};

struct PhaseInfo {
    string raw;
    bool before = false;
    double temperature = numeric_limits<double>::quiet_NaN();
    double hours = numeric_limits<double>::quiet_NaN();
};

static string trim_copy(const string& s)
{
    const char* ws = " \t\n\r\f\v";
    const auto begin = s.find_first_not_of(ws);
    if (begin == string::npos) return "";
    const auto end = s.find_last_not_of(ws);
    return s.substr(begin, end - begin + 1);
}

static vector<string> split_csv_simple(const string& line)
{
    vector<string> out;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ',')) out.push_back(trim_copy(item));
    return out;
}

static bool finite_value(double x)
{
    return std::isfinite(x);
}

// Supports both naming conventions used in the project:
//   before_annealing
//   annealing_75_5
//   annealing_T=75_h=5
static PhaseInfo parse_phase(const string& phase)
{
    PhaseInfo p;
    p.raw = phase;

    if (phase == "before_annealing") {
        p.before = true;
        return p;
    }

    // New whole-picture naming convention: annealing_75_5
    const string simple_prefix = "annealing_";
    if (phase.find(simple_prefix) == 0 && phase.find("annealing_T=") != 0) {
        const string rest = phase.substr(simple_prefix.size());
        const size_t pos = rest.find('_');
        if (pos != string::npos) {
            try {
                p.temperature = stod(rest.substr(0, pos));
                p.hours = stod(rest.substr(pos + 1));
                return p;
            } catch (...) {
                p.temperature = numeric_limits<double>::quiet_NaN();
                p.hours = numeric_limits<double>::quiet_NaN();
            }
        }
    }

    // Older naming convention: annealing_T=75_h=5
    const string old_prefix = "annealing_T=";
    if (phase.find(old_prefix) == 0) {
        const size_t pos_h = phase.find("_h=");
        if (pos_h != string::npos) {
            try {
                p.temperature = stod(phase.substr(old_prefix.size(), pos_h - old_prefix.size()));
                p.hours = stod(phase.substr(pos_h + 3));
            } catch (...) {
                p.temperature = numeric_limits<double>::quiet_NaN();
                p.hours = numeric_limits<double>::quiet_NaN();
            }
        }
    }

    return p;
}

// spacing logic:
// 5 h -> +1 x-unit
// 25 h -> +5 x-units
static double phase_step_units(const string& phase)
{
    const PhaseInfo p = parse_phase(phase);
    if (p.before) return 0.0;

    if (finite_value(p.hours) && p.hours > 0.0)
        return p.hours / 5.0;

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

    if (p.before)
        return "#splitline{Before}{annealing}";

    if (finite_value(p.temperature) && finite_value(p.hours)) {
        const double Tround = round(p.temperature);
        const double hround = round(p.hours);

        string Ttext;
        string htext;

        if (fabs(p.temperature - Tround) < 1e-9)
            Ttext = Form("%.0f#circC", p.temperature);
        else
            Ttext = Form("%.1f#circC", p.temperature);

        if (fabs(p.hours - hround) < 1e-9)
            htext = Form("%.0f h", p.hours);
        else
            htext = Form("%.1f h", p.hours);

        return "#splitline{" + Ttext + "}{" + htext + "}";
    }

    return phase;
}

static vector<SummaryRow> read_summary_csv(const char* csv_path)
{
    ifstream fin(csv_path);
    if (!fin.is_open())
        throw runtime_error(string("Cannot open CSV file: ") + csv_path);

    string header;
    if (!getline(fin, header))
        throw runtime_error(string("CSV file is empty: ") + csv_path);

    vector<string> header_fields = split_csv_simple(header);
    map<string, int> col;
    for (int i = 0; i < (int)header_fields.size(); ++i)
        col[header_fields[i]] = i;

    const vector<string> required = {
        "sensor", "phase", "counted_hotspots", "counted_hotspots_error"
    };

    for (const auto& name : required) {
        if (!col.count(name))
            throw runtime_error("Missing required column in CSV: " + name);
    }

    vector<SummaryRow> rows;
    string line;

    while (getline(fin, line)) {
        line = trim_copy(line);
        if (line.empty()) continue;

        // Harmless if the CSV was copied from a text block ending lines with '\\'.
        if (!line.empty() && line.back() == '\\') line.pop_back();
        line = trim_copy(line);
        if (line.empty()) continue;

        vector<string> fields = split_csv_simple(line);
        const int max_needed = max({
            col["sensor"], col["phase"], col["counted_hotspots"], col["counted_hotspots_error"]
        });

        if ((int)fields.size() <= max_needed) {
            cerr << "WARNING: malformed row skipped:\n" << line << endl;
            continue;
        }

        try {
            SummaryRow r;
            r.sensor = fields[col["sensor"]];
            r.phase = fields[col["phase"]];
            r.counted_hotspots = stod(fields[col["counted_hotspots"]]);
            r.counted_hotspots_error = stod(fields[col["counted_hotspots_error"]]);

            if (!finite_value(r.counted_hotspots)) continue;
            if (!finite_value(r.counted_hotspots_error)) continue;
            rows.push_back(r);
        } catch (...) {
            cerr << "WARNING: malformed row skipped:\n" << line << endl;
        }
    }

    if (rows.empty())
        throw runtime_error("No valid rows found in CSV.");

    return rows;
}

// Draw only the annealing TIME labels below the x axis.
// The true non-uniform x positions are preserved exactly as built by
// build_phase_positions(): 5 h -> +1 unit, 25 h -> +5 units.
// Temperature labels are drawn separately inside the plot area.
static void draw_time_labels(
    TPad* pad,
    const vector<string>& phases,
    const vector<double>& xpos,
    double xmin,
    double xmax)
{
    if (!pad || phases.empty() || xpos.empty() || !(xmax > xmin)) return;

    pad->Update();

    const double left   = pad->GetLeftMargin();
    const double right  = pad->GetRightMargin();
    const double usable = 1.0 - left - right;

    for (size_t i = 0; i < phases.size(); ++i) {
        const double frac  = (xpos[i] - xmin) / (xmax - xmin);
        const double x_ndc = left + frac * usable;

        string label;
        const PhaseInfo p = parse_phase(phases[i]);

        // The baseline point has no annealing time, so it gets no bottom label.
        // It is identified separately at the top of the plot.
        if (p.before) continue;

        if (finite_value(p.hours)) {
            const double hround = round(p.hours);
            if (fabs(p.hours - hround) < 1e-9)
                label = Form("%.0f h", p.hours);
            else
                label = Form("%.1f h", p.hours);
        } else {
            label = phases[i];
        }

        TLatex* lab = new TLatex();
        lab->SetNDC();
        lab->SetTextFont(42);
        lab->SetTextSize(0.032);
        lab->SetTextAlign(22);
        lab->DrawLatex(x_ndc, 0.13, label.c_str());
    }

    TLatex* title = new TLatex();
    title->SetNDC();
    title->SetTextFont(42);
    title->SetTextSize(0.047);
    title->SetTextAlign(22);
    title->DrawLatex(0.5, 0.035, "Annealing time");
}

// Draw one temperature label centred above each 5 h / 25 h pair and
// vertical lines separating the temperature blocks.
static void draw_temperature_blocks(
    const vector<string>& phases,
    const vector<double>& xpos,
    double ymin,
    double ymax)
{
    if (phases.size() != xpos.size() || phases.empty()) return;

    map<double, vector<double>> temp_positions;
    for (size_t i = 0; i < phases.size(); ++i) {
        const PhaseInfo p = parse_phase(phases[i]);
        if (!p.before && finite_value(p.temperature))
            temp_positions[p.temperature].push_back(xpos[i]);
    }

    const double yrange = ymax - ymin;
    const double ytext  = ymax - 0.055 * yrange;

    TLatex latex;
    latex.SetTextFont(42);
    latex.SetTextSize(0.034);
    latex.SetTextAlign(22);

    // Baseline label: it is not an annealing-time point, so keep it out of
    // the bottom x-axis labels and identify it here instead.
    latex.SetTextSize(0.030);
    latex.DrawLatex(xpos.front()-1, ytext, "#splitline{Before}{annealing}");
    latex.SetTextSize(0.034);

    for (const auto& kv : temp_positions) {
        if (kv.second.empty()) continue;

        double xmin_block = *min_element(kv.second.begin(), kv.second.end());
        double xmax_block = *max_element(kv.second.begin(), kv.second.end());
        const double xcentre = 0.5 * (xmin_block + xmax_block);

        const double Tround = round(kv.first);
        string tlabel;
        if (fabs(kv.first - Tround) < 1e-9)
            tlabel = Form("T = %.0f#circC", kv.first);
        else
            tlabel = Form("T = %.1f#circC", kv.first);

        latex.DrawLatex(xcentre, ytext, tlabel.c_str());
    }

    // Boundaries between baseline / 75 C and between successive temperatures.
    // With x = 0,1,6,7,12,13,18,19,24 these are 0.5, 6.5, 12.5, 18.5.
    const vector<double> separators = {0.5, 6.5, 12.5, 18.5};
    for (double xsep : separators) {
        TLine* line = new TLine(xsep, ymin, xsep, ymax);
        line->SetLineColor(kGray + 1);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("same");
    }
}

void plot_counted_hotspots_vs_phase(
    const char* csv_path = "hotspot_phase_summary_all_sensor_error.csv",
    const char* output_prefix = "counted_hotspots_vs_phase")
{
    gStyle->SetOptStat(0);

    const vector<string> phases = {
        "before_annealing",
        "annealing_75_5",
        "annealing_75_25",
        "annealing_100_5",
        "annealing_100_25",
        "annealing_125_5",
        "annealing_125_25",
        "annealing_150_5",
        "annealing_150_25"
    };

    // IMPORTANT: this is exactly the same cumulative spacing prescription as
    // plot_lambda_vs_phase.C:
    //   before            -> x = 0
    //   75 C, 5 h         -> x = 1
    //   75 C, 25 h        -> x = 6
    //   100 C, 5 h        -> x = 7
    //   100 C, 25 h       -> x = 12
    //   ...
    // Therefore every 5 h step occupies 1 x-unit and every 25 h step 5 units.
    const vector<double> x = build_phase_positions(phases);

    map<string, double> x_phase;
    for (size_t i = 0; i < phases.size(); ++i)
        x_phase[phases[i]] = x[i];

    const vector<string> sensors = {"A1", "A2", "B1", "B2"};

    const map<string, Color_t> sensor_color = {
        {"A1", kP6Red},
        {"A2", kP6Blue},
        {"B1", kP6Yellow},
        {"B2", kP8Green}
    };

    // Marker shapes are intentionally different as a second visual cue.
    const map<string, Style_t> sensor_marker = {
        {"A1", 20},
        {"A2", 21},
        {"B1", 22},
        {"B2", 33}
    };

    const vector<SummaryRow> rows = read_summary_csv(csv_path);

    map<string, map<string, double>> data;
    // Systematic uncertainties
    map<string, map<string, double>> data_error;
    double ymin_data = numeric_limits<double>::infinity();
    double ymax_data = -numeric_limits<double>::infinity();

    for (const auto& r : rows) {
        if (!x_phase.count(r.phase)) continue;
        data[r.sensor][r.phase] = r.counted_hotspots;
        data_error[r.sensor][r.phase] = r.counted_hotspots_error;
        // Include the systematic uncertainty in the y-axis range.
        ymin_data = min(ymin_data, r.counted_hotspots - r.counted_hotspots_error);
        ymax_data = max(ymax_data, r.counted_hotspots + r.counted_hotspots_error);
    }

    if (!finite_value(ymin_data) || !finite_value(ymax_data))
        throw runtime_error("No valid counted_hotspots values found in the CSV.");

    double yrange = ymax_data - ymin_data;
    if (!(yrange > 0.0)) yrange = max(1.0, fabs(ymax_data) * 0.10);

    const double ymin = max(0.0, ymin_data - 0.12 * yrange);
    const double ymax = ymax_data + 0.25 * yrange;

    // Same frame philosophy as plot_lambda_vs_phase.C.
    const int max_unit = max(0, (int)llround(x.back()));
    // Extend the numerical x range to the left of the before-annealing
    // point. The actual phase coordinates are NOT changed, so the relative
    // spacing remains exactly 5 h -> +1 unit and 25 h -> +5 units.
    //
    // This extra empty x-range moves the x=0 "Before annealing" label away
    // from the y-axis labels/title without altering the annealing spacing.
    const double xmin = -2.5;
    const double xmax = max_unit + 0.5;
    const int nbins = max_unit + 1;

    TCanvas* c = new TCanvas(
        "c_counted_hotspots",
        "counted hotspots vs annealing phase",
        2250,
        1000
    );

    c->SetBottomMargin(0.20);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.09);
    c->SetTicks(1,1);

    TH1D* frame = new TH1D(
        "frame_counted_hotspots",
        ";;Counted hotspots",
        nbins,
        xmin,
        xmax
    );

    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);

    // We hide ROOT's own x labels and divisions. The phase labels and ticks
    // are placed manually at the actual non-uniform x coordinates.
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetNdivisions(0);
    frame->GetXaxis()->SetTickLength(0.0);

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.047);
    frame->GetYaxis()->SetTitleOffset(1.05);

    frame->Draw();

    // Small x-axis tick marks exactly at each phase position.
    const double tick_height = 0.018 * (ymax - ymin);
    for (double xpos : x) {
        TLine* tick = new TLine(xpos, ymin, xpos, ymin + tick_height);
        tick->SetLineColor(kBlack);
        tick->SetLineWidth(1);
        tick->Draw("same");
    }

    // Temperature labels and vertical separators are drawn first, so the
    // sensor curves and markers remain visually on top of the separator lines.
    draw_temperature_blocks(phases, x, ymin, ymax);

    TLegend* leg = new TLegend(0.80, 0.60, 0.93, 0.78);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);

    vector<TGraph*> graphs;
    vector<TGraphErrors*> syst_graphs;

    for (const auto& sensor : sensors) {

        vector<double> xv;
        vector<double> yv;

        vector<double> exv;
        vector<double> eyv;


        auto sensor_it = data.find(sensor);

        if (sensor_it == data.end())
            continue;


        auto error_sensor_it = data_error.find(sensor);

        if (error_sensor_it == data_error.end())
            continue;


        for (size_t i = 0; i < phases.size(); ++i) {

            auto phase_it =
                sensor_it->second.find(phases[i]);

            if (phase_it == sensor_it->second.end())
                continue;


            auto error_phase_it =
                error_sensor_it->second.find(phases[i]);

            if (error_phase_it ==
                error_sensor_it->second.end())
                continue;


            xv.push_back(x[i]);

            yv.push_back(
                phase_it->second
            );

            // No uncertainty along x.
            exv.push_back(0.0);

            // Systematic uncertainty along y.
            eyv.push_back(
                error_phase_it->second
            );
        }


        if (xv.empty())
            continue;


        // =====================================================
        // SYSTEMATIC UNCERTAINTY
        //
        // Draw only ROOT brackets.
        // =====================================================
        TGraphErrors* gr_syst =
            new TGraphErrors(
                (int)xv.size(),
                xv.data(),
                yv.data(),
                exv.data(),
                eyv.data()
            );


        gr_syst->SetName(
            Form("gr_syst_%s", sensor.c_str())
        );

        gr_syst->SetLineColor(
            sensor_color.at(sensor)
        );

        gr_syst->SetLineWidth(2);


        // ROOT option []:
        // draw uncertainties as brackets.
        gr_syst->Draw("[] SAME");


        // =====================================================
        // CENTRAL VALUES + CONNECTING LINE
        // =====================================================
        TGraph* gr =
            new TGraph(
                (int)xv.size(),
                xv.data(),
                yv.data()
            );


        gr->SetName(
            Form("gr_%s", sensor.c_str())
        );

        gr->SetLineColor(
            sensor_color.at(sensor)
        );

        gr->SetMarkerColor(
            sensor_color.at(sensor)
        );

        gr->SetMarkerStyle(
            sensor_marker.at(sensor)
        );

        gr->SetMarkerSize(1.45);

        gr->SetLineWidth(3);


        // Draw after systematic brackets so that
        // markers remain clearly visible.
        gr->Draw("PL SAME");


        leg->AddEntry(
            gr,
            sensor.c_str(),
            "lp"
        );


        graphs.push_back(gr);
        syst_graphs.push_back(gr_syst);
    }

    leg->Draw();

    c->Modified();
    c->Update();

    // Only time labels are shown below the x axis.
    draw_time_labels(
        c,
        phases,
        x,
        frame->GetXaxis()->GetXmin(),
        frame->GetXaxis()->GetXmax()
    );

    gPad->RedrawAxis();
    c->Modified();
    c->Update();

    c->SaveAs(Form("%s.png", output_prefix));
    c->SaveAs(Form("%s.pdf", output_prefix));

    cout << "Saved plots:" << endl;
    cout << "  " << output_prefix << ".png" << endl;
    cout << "  " << output_prefix << ".pdf" << endl;

    cout << "Phase x positions (same prescription as plot_lambda_vs_phase.C):" << endl;
    for (size_t i = 0; i < phases.size(); ++i)
        cout << "  " << phases[i] << " -> x = " << x[i] << endl;

    cout << "Spacing rule: 5 h = +1 unit, 25 h = +5 units." << endl;
}
