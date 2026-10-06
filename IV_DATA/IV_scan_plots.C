// ============================================================================
// IV_scan_plots.C
//
// Read:
//   A1_IV_scan.csv
//   B1_IV_scan.csv
//
// Expected CSV structure:
//   phase,temperature_C,voltage,current,current_error
//
// All plots use the signed zero-subtracted current from makeiv-like processing.
// The y axis is logarithmic, therefore only points with current > 0 are drawn.
//
// Produce three sets of IV plots in:
//   IV_plots/1_T21_single_phase/
//   IV_plots/2_all_temperatures/
//   IV_plots/3_T21_all_phases/
//
// Usage from emmi/IV_DATA:
//
//   root -l -q 'IV_scan_plots.C(".")'
//
// or, interactively:
//
//   root -l
//   .L IV_scan_plots.C+
//   IV_scan_plots(".")
//
// ============================================================================

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "Rtypes.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
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

namespace {

// ============================================================================
// Data structures
// ============================================================================

struct IVPoint {
    string phase;
    double temperature_C = 0.0;
    double voltage = 0.0;
    double current = 0.0;
    double current_error = 0.0;
};

struct PhaseKey {
    int group = 2;
    double annealing_T = 1e99;
    double annealing_h = 1e99;
    string name;
};


// ============================================================================
// String helpers
// ============================================================================

string trim_copy(const string& s)
{
    const size_t first = s.find_first_not_of(" \t\r\n");

    if (first == string::npos)
        return "";

    const size_t last = s.find_last_not_of(" \t\r\n");

    return s.substr(first, last-first+1);
}


vector<string> split_csv_simple(const string& line)
{
    vector<string> fields;
    std::stringstream ss(line);
    string field;

    while (std::getline(ss, field, ','))
        fields.push_back(trim_copy(field));

    return fields;
}


int find_col(const vector<string>& header, const string& name)
{
    for (int i=0; i<(int)header.size(); ++i) {
        if (trim_copy(header[i]) == name)
            return i;
    }

    return -1;
}


bool finite_number(double x)
{
    return std::isfinite(x);
}


string filename_safe(string s)
{
    for (char& c : s) {
        const bool ok =
            (c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') ||
            c == '_' || c == '-';

        if (!ok)
            c = '_';
    }

    return s;
}


// ============================================================================
// Phase handling
// ============================================================================

PhaseKey phase_key(const string& phase)
{
    PhaseKey k;
    k.name = phase;

    if (phase == "bef_ann" ||
        phase == "before_annealing") {

        k.group = 0;
        k.annealing_T = -1e99;
        k.annealing_h = -1e99;

        return k;
    }

    double T = 0.0;
    double h = 0.0;

    // Convention produced by IV_scan.sh:
    // annealing_T=75_h=5
    if (std::sscanf(
            phase.c_str(),
            "annealing_T=%lf_h=%lf",
            &T,
            &h
        ) == 2) {

        k.group = 1;
        k.annealing_T = T;
        k.annealing_h = h;

        return k;
    }

    // Also accept:
    // annealing_75_5
    if (std::sscanf(
            phase.c_str(),
            "annealing_%lf_%lf",
            &T,
            &h
        ) == 2) {

        k.group = 1;
        k.annealing_T = T;
        k.annealing_h = h;

        return k;
    }

    // Also accept:
    // ann_75_5
    if (std::sscanf(
            phase.c_str(),
            "ann_%lf_%lf",
            &T,
            &h
        ) == 2) {

        k.group = 1;
        k.annealing_T = T;
        k.annealing_h = h;

        return k;
    }

    return k;
}


bool phase_less(const string& a, const string& b)
{
    const PhaseKey ka = phase_key(a);
    const PhaseKey kb = phase_key(b);

    if (ka.group != kb.group)
        return ka.group < kb.group;

    if (ka.annealing_T != kb.annealing_T)
        return ka.annealing_T < kb.annealing_T;

    if (ka.annealing_h != kb.annealing_h)
        return ka.annealing_h < kb.annealing_h;

    return ka.name < kb.name;
}


string pretty_phase(const string& phase)
{
    if (phase == "bef_ann" ||
        phase == "before_annealing")
        return "Before annealing";

    const PhaseKey k = phase_key(phase);

    if (k.group == 1) {
        return Form(
            "%.0f #circC, %.0f h",
            k.annealing_T,
            k.annealing_h
        );
    }

    return phase;
}


// ============================================================================
// Hard-coded phase colors
//
// IMPORTANT:
// The color depends on the PHASE NAME, not on the position of the phase in
// the input file. Therefore, if a phase is absent, its color is skipped.
//
// Standard annealing sequence:
//   before
//   75 C  5 h
//   75 C 25 h
//   100 C  5 h
//   100 C 25 h
//   125 C  5 h
//   125 C 25 h
//   150 C  5 h
//   150 C 25 h
//
// 150 C, 25 h is explicitly assigned kPink+7 as requested.
// ============================================================================

int phase_color(const string& phase)
{
    const PhaseKey k = phase_key(phase);

    if (k.group == 0)
        return kBlack;

    if (k.group == 1) {

        if (std::fabs(k.annealing_T - 75.0) < 1e-9 &&
            std::fabs(k.annealing_h - 5.0) < 1e-9)
            return kRed+1;

        if (std::fabs(k.annealing_T - 75.0) < 1e-9 &&
            std::fabs(k.annealing_h - 25.0) < 1e-9)
            return kOrange+7;

        if (std::fabs(k.annealing_T - 100.0) < 1e-9 &&
            std::fabs(k.annealing_h - 5.0) < 1e-9)
            return kGreen+2;

        if (std::fabs(k.annealing_T - 100.0) < 1e-9 &&
            std::fabs(k.annealing_h - 25.0) < 1e-9)
            return kCyan+2;

        if (std::fabs(k.annealing_T - 125.0) < 1e-9 &&
            std::fabs(k.annealing_h - 5.0) < 1e-9)
            return kAzure+1;

        if (std::fabs(k.annealing_T - 125.0) < 1e-9 &&
            std::fabs(k.annealing_h - 25.0) < 1e-9)
            return kBlue+1;

        if (std::fabs(k.annealing_T - 150.0) < 1e-9 &&
            std::fabs(k.annealing_h - 5.0) < 1e-9)
            return kViolet+1;

        if (std::fabs(k.annealing_T - 150.0) < 1e-9 &&
            std::fabs(k.annealing_h - 25.0) < 1e-9)
            return kPink+7;
    }

    // Reserved fallback colors for additional/unexpected phases.
    // They keep the same family requested for the annealing plots.
    static int fallback_index = 0;
    static const int fallback_colors[] = {
        kMagenta+1,
        kGray+2
    };

    const int color =
        fallback_colors[
            fallback_index %
            (sizeof(fallback_colors)/sizeof(fallback_colors[0]))
        ];

    ++fallback_index;

    return color;
}


// ============================================================================
// Temperature colors
//
// Temperatures are sorted in increasing order before this function is called.
//
// Requested initial sequence:
//   kP10Violet
//   kP10Blue
//   kP8Pink
//   kP10Green
//   kP10Orange
//
// Further entries use colors from the same accessible ROOT color families.
// ============================================================================

int temperature_color(size_t index)
{
    static const int colors[] = {
        kP10Violet,
        kP10Blue,
        kP8Pink,
        kP10Green,
        kP10Orange,
        kP10Cyan,
        kP10Red,
        kP10Brown,
        kP10Ash,
        kP10Gray
    };

    return colors[
        index %
        (sizeof(colors)/sizeof(colors[0]))
    ];
}


// ============================================================================
// CSV reader
// ============================================================================

vector<IVPoint> read_iv_csv(const string& filename)
{
    vector<IVPoint> rows;

    std::ifstream fin(filename);

    if (!fin.is_open()) {
        std::cerr
            << "ERROR: cannot open "
            << filename
            << std::endl;

        return rows;
    }

    string line;

    if (!std::getline(fin, line)) {
        std::cerr
            << "ERROR: empty CSV: "
            << filename
            << std::endl;

        return rows;
    }

    const vector<string> header =
        split_csv_simple(line);

    const int cphase =
        find_col(header, "phase");

    const int ctemp =
        find_col(header, "temperature_C");

    const int cvoltage =
        find_col(header, "voltage");

    const int ccurrent =
        find_col(header, "current");

    const int cerror =
        find_col(header, "current_error");

    if (cphase < 0 ||
        ctemp < 0 ||
        cvoltage < 0 ||
        ccurrent < 0 ||
        cerror < 0) {

        std::cerr
            << "ERROR: missing required columns in "
            << filename
            << std::endl;

        std::cerr
            << "Required columns:"
            << std::endl
            << "  phase"
            << std::endl
            << "  temperature_C"
            << std::endl
            << "  voltage"
            << std::endl
            << "  current"
            << std::endl
            << "  current_error"
            << std::endl;

        return rows;
    }

    const int maxcol =
        std::max(
            std::max(cphase, ctemp),
            std::max(
                cvoltage,
                std::max(ccurrent, cerror)
            )
        );

    int line_number = 1;

    while (std::getline(fin, line)) {

        ++line_number;

        if (trim_copy(line).empty())
            continue;

        vector<string> fields =
            split_csv_simple(line);

        if ((int)fields.size() <= maxcol) {
            std::cerr
                << "WARNING: malformed row "
                << line_number
                << " in "
                << filename
                << std::endl;

            continue;
        }

        try {

            IVPoint r;

            r.phase =
                trim_copy(fields[cphase]);

            r.temperature_C =
                std::stod(fields[ctemp]);

            r.voltage =
                std::stod(fields[cvoltage]);

            r.current =
                std::stod(fields[ccurrent]);

            r.current_error =
                std::stod(fields[cerror]);

            if (r.phase.empty())
                continue;

            if (!finite_number(r.temperature_C) ||
                !finite_number(r.voltage) ||
                !finite_number(r.current))
                continue;

            // If a single-voltage statistical uncertainty could not be
            // evaluated, keep the point and draw it with zero error.
            if (!finite_number(r.current_error) ||
                r.current_error < 0.0)
                r.current_error = 0.0;

            rows.push_back(r);
        }
        catch (...) {

            std::cerr
                << "WARNING: cannot parse row "
                << line_number
                << " in "
                << filename
                << std::endl;
        }
    }

    return rows;
}


// ============================================================================
// Data helpers
// ============================================================================

bool point_voltage_less(
    const IVPoint& a,
    const IVPoint& b
)
{
    return a.voltage < b.voltage;
}


vector<IVPoint> select_phase_temperature(
    const vector<IVPoint>& rows,
    const string& phase,
    double temperature
)
{
    vector<IVPoint> out;

    for (const auto& r : rows) {

        if (r.phase != phase)
            continue;

        if (std::fabs(
                r.temperature_C - temperature
            ) > 1e-6)
            continue;

        out.push_back(r);
    }

    std::sort(
        out.begin(),
        out.end(),
        point_voltage_less
    );

    return out;
}


set<double> get_temperatures_for_phase(
    const vector<IVPoint>& rows,
    const string& phase
)
{
    set<double> temperatures;

    for (const auto& r : rows) {
        if (r.phase == phase)
            temperatures.insert(r.temperature_C);
    }

    return temperatures;
}


vector<string> get_phases(
    const vector<IVPoint>& rows
)
{
    set<string> phase_set;

    for (const auto& r : rows) {
        if (!r.phase.empty())
            phase_set.insert(r.phase);
    }

    vector<string> phases(
        phase_set.begin(),
        phase_set.end()
    );

    std::sort(
        phases.begin(),
        phases.end(),
        phase_less
    );

    return phases;
}


// ============================================================================
// Graph creation
// ============================================================================

TGraphErrors* make_graph(
    const vector<IVPoint>& points,
    const string& name,
    int color,
    int marker_style = 20,
    double marker_size = 1.10
)
{
    // ------------------------------------------------------------------------
    // IMPORTANT: do NOT reflect negative zero-subtracted currents with fabs().
    // The reference makeiv.C keeps the sign after subtracting the zero level.
    // Since these canvases use log(y), only strictly positive currents can be
    // represented and are therefore retained here.
    // ------------------------------------------------------------------------

    vector<double> x;
    vector<double> y;
    vector<double> ex;
    vector<double> ey;

    x.reserve(points.size());
    y.reserve(points.size());
    ex.reserve(points.size());
    ey.reserve(points.size());

    for (const auto& p : points) {

        if (!(p.current > 0.0) ||
            !finite_number(p.current))
            continue;

        x.push_back(p.voltage);
        y.push_back(p.current);
        ex.push_back(0.0);
        ey.push_back(std::fabs(p.current_error));
    }

    if (x.empty())
        return nullptr;

    const int n =
        (int)x.size();

    TGraphErrors* graph =
        new TGraphErrors(
            n,
            x.data(),
            y.data(),
            ex.data(),
            ey.data()
        );

    graph->SetName(
        name.c_str()
    );

    graph->SetMarkerStyle(
        marker_style
    );

    graph->SetMarkerSize(
        marker_size
    );

    graph->SetMarkerColor(
        color
    );

    graph->SetLineColor(
        color
    );

    graph->SetLineWidth(2);

    return graph;
}


// ============================================================================
// Plot ranges
// ============================================================================

void update_range(
    const vector<IVPoint>& points,
    double& xmin,
    double& xmax,
    double& ymin,
    double& ymax
)
{
    for (const auto& p : points) {

        // Keep the visible x range consistent with the points that can actually
        // be drawn on a logarithmic y axis.
        if (!(p.current > 0.0) ||
            !finite_number(p.current))
            continue;

        xmin = std::min(xmin, p.voltage);
        xmax = std::max(xmax, p.voltage);

        const double y  = p.current;
        const double ey = std::fabs(p.current_error);

        const double lower = y - ey;

        if (lower > 0.0 && finite_number(lower))
            ymin = std::min(ymin, lower);
        else
            ymin = std::min(ymin, y);

        const double upper = y + ey;

        if (upper > 0.0 && finite_number(upper))
            ymax = std::max(ymax, upper);
        else
            ymax = std::max(ymax, y);
    }
}


void expand_range(
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double& frame_xmin,
    double& frame_xmax,
    double& frame_ymin,
    double& frame_ymax
)
{
    // Linear padding on x.
    // If no positive-y point was available, callers should normally have
    // skipped graph creation; these defaults keep the frame numerically safe.
    if (!finite_number(xmin) || !finite_number(xmax) || xmin > xmax) {
        xmin = 0.0;
        xmax = 1.0;
    }

    double xspan = xmax - xmin;

    if (!(xspan > 0.0))
        xspan = std::max(1.0, std::fabs(xmax)*0.10);

    frame_xmin = xmin - 0.04*xspan;
    frame_xmax = xmax + 0.04*xspan;

    // Logarithmic padding on y.
    // This guarantees strictly positive frame boundaries.
    if (!(ymin > 0.0) || !finite_number(ymin))
        ymin = 1e-15;

    if (!(ymax > ymin) || !finite_number(ymax))
        ymax = ymin * 10.0;

    const double log_ymin = std::log10(ymin);
    const double log_ymax = std::log10(ymax);

    double log_span = log_ymax - log_ymin;

    if (!(log_span > 0.0))
        log_span = 1.0;

    frame_ymin = std::pow(10.0, log_ymin - 0.10*log_span);
    frame_ymax = std::pow(10.0, log_ymax + 0.15*log_span);
}


// ============================================================================
// Common canvas/frame style
// ============================================================================

TCanvas* create_canvas(
    const string& name,
    const string& title
)
{
    TCanvas* c =
        new TCanvas(
            name.c_str(),
            title.c_str(),
            1500,
            950
        );

    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.10);

    return c;
}


TH1* create_frame(
    TCanvas* c,
    double xmin,
    double ymin,
    double xmax,
    double ymax,
    const string& title
)
{
    TH1* frame =
        c->DrawFrame(
            xmin,
            ymin,
            xmax,
            ymax,
            title.c_str()
        );

    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);

    frame->GetXaxis()->SetLabelSize(0.038);
    frame->GetYaxis()->SetLabelSize(0.038);

    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetTitleOffset(1.30);

    return frame;
}


void save_canvas(
    TCanvas* c,
    const string& output_base
)
{
    c->Modified();
    c->Update();

    c->SaveAs(
        (output_base + ".png").c_str()
    );

    c->SaveAs(
        (output_base + ".pdf").c_str()
    );
}


// ============================================================================
// SET 1
//
// Individual IV curve for each annealing phase at T = 21 C.
// ============================================================================

int produce_set1(
    const string& sensor,
    const vector<IVPoint>& rows,
    const vector<string>& phases,
    const string& output_dir
)
{
    const double target_temperature =
        21.0;

    int created = 0;

    for (const string& phase : phases) {

        vector<IVPoint> points =
            select_phase_temperature(
                rows,
                phase,
                target_temperature
            );

        if (points.empty())
            continue;

        double xmin =  1e99;
        double xmax = -1e99;
        double ymin =  1e99;
        double ymax = -1e99;

        update_range(
            points,
            xmin,
            xmax,
            ymin,
            ymax
        );

        double fxmin, fxmax, fymin, fymax;

        expand_range(
            xmin,
            xmax,
            ymin,
            ymax,
            fxmin,
            fxmax,
            fymin,
            fymax
        );

        const string phase_file =
            filename_safe(phase);

        TCanvas* c =
            create_canvas(
                Form(
                    "c_set1_%s_%s",
                    sensor.c_str(),
                    phase_file.c_str()
                ),
                Form(
                    "%s - %s - T=21 C",
                    sensor.c_str(),
                    phase_file.c_str()
                )
            );

        c->SetLogy();

        const string plot_title =
            Form(
                "Sensor %s - %s - T = 21 #circC;"
                "Voltage (V);Current (A)",
                sensor.c_str(),
                pretty_phase(phase).c_str()
            );

        create_frame(
            c,
            fxmin,
            fymin,
            fxmax,
            fymax,
            plot_title
        );

        TGraphErrors* graph =
            make_graph(
                points,
                Form(
                    "gr_set1_%s_%s",
                    sensor.c_str(),
                    phase_file.c_str()
                ),
                kBlack,
                20,
                1.15
            );

        if (!graph) {
            delete c;
            continue;
        }

        // P = markers
        // Z = vertical error bars without horizontal end caps
        graph->Draw("PZ SAME");

        const string outbase =
            output_dir
            + "/"
            + sensor
            + "_"
            + phase_file
            + "_T21_IV";

        save_canvas(
            c,
            outbase
        );

        delete c;
        ++created;
    }

    return created;
}


// ============================================================================
// SET 2
//
// Same sensor + same annealing phase, all available temperatures.
//
// Filled markers only, vertical error bars only.
// ============================================================================

int produce_set2(
    const string& sensor,
    const vector<IVPoint>& rows,
    const vector<string>& phases,
    const string& output_dir
)
{
    int created = 0;

    for (const string& phase : phases) {

        const set<double> temp_set =
            get_temperatures_for_phase(
                rows,
                phase
            );

        if (temp_set.empty())
            continue;

        vector<double> temperatures(
            temp_set.begin(),
            temp_set.end()
        );

        double xmin =  1e99;
        double xmax = -1e99;
        double ymin =  1e99;
        double ymax = -1e99;

        map<double,vector<IVPoint>> curves;

        for (double temperature : temperatures) {

            vector<IVPoint> points =
                select_phase_temperature(
                    rows,
                    phase,
                    temperature
                );

            if (points.empty())
                continue;

            curves[temperature] =
                points;

            update_range(
                points,
                xmin,
                xmax,
                ymin,
                ymax
            );
        }

        if (curves.empty())
            continue;

        double fxmin, fxmax, fymin, fymax;

        expand_range(
            xmin,
            xmax,
            ymin,
            ymax,
            fxmin,
            fxmax,
            fymin,
            fymax
        );

        const string phase_file =
            filename_safe(phase);

        TCanvas* c =
            create_canvas(
                Form(
                    "c_set2_%s_%s",
                    sensor.c_str(),
                    phase_file.c_str()
                ),
                Form(
                    "%s - %s - all temperatures",
                    sensor.c_str(),
                    phase_file.c_str()
                )
            );

        c->SetLogy();

        const string plot_title =
            Form(
                "Sensor %s - %s;"
                "Voltage (V);Current (A)",
                sensor.c_str(),
                pretty_phase(phase).c_str()
            );

        create_frame(
            c,
            fxmin,
            fymin,
            fxmax,
            fymax,
            plot_title
        );

        TLegend* leg =
            new TLegend(
                0.15,
                0.66,
                0.39,
                0.88
            );

        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.030);

        size_t itemp = 0;

        for (double temperature : temperatures) {

            auto it =
                curves.find(temperature);

            if (it == curves.end())
                continue;

            const int color =
                temperature_color(itemp);

            TGraphErrors* graph =
                make_graph(
                    it->second,
                    Form(
                        "gr_set2_%s_%s_T%.0f",
                        sensor.c_str(),
                        phase_file.c_str(),
                        temperature
                    ),
                    color,
                    20,
                    1.10
                );

            if (!graph) {
                ++itemp;
                continue;
            }

            graph->Draw("PZ SAME");

            leg->AddEntry(
                graph,
                Form(
                    "T = %.0f #circC",
                    temperature
                ),
                "pe"
            );

            ++itemp;
        }

        leg->Draw();

        const string outbase =
            output_dir
            + "/"
            + sensor
            + "_"
            + phase_file
            + "_all_temperatures_IV";

        save_canvas(
            c,
            outbase
        );

        delete c;
        ++created;
    }

    return created;
}


// ============================================================================
// SET 3
//
// Same sensor + T = 21 C, all annealing phases.
//
// Phase colors are hard-coded through phase_color().
// ============================================================================

int produce_set3(
    const string& sensor,
    const vector<IVPoint>& rows,
    const vector<string>& phases,
    const string& output_dir
)
{
    const double target_temperature =
        21.0;

    double xmin =  1e99;
    double xmax = -1e99;
    double ymin =  1e99;
    double ymax = -1e99;

    map<string,vector<IVPoint>> curves;

    for (const string& phase : phases) {

        vector<IVPoint> points =
            select_phase_temperature(
                rows,
                phase,
                target_temperature
            );

        if (points.empty())
            continue;

        curves[phase] =
            points;

        update_range(
            points,
            xmin,
            xmax,
            ymin,
            ymax
        );
    }

    if (curves.empty())
        return 0;

    double fxmin, fxmax, fymin, fymax;

    expand_range(
        xmin,
        xmax,
        ymin,
        ymax,
        fxmin,
        fxmax,
        fymin,
        fymax
    );

    TCanvas* c =
        create_canvas(
            Form(
                "c_set3_%s_T21_all_phases",
                sensor.c_str()
            ),
            Form(
                "%s - T=21 C - all annealing phases",
                sensor.c_str()
            )
        );

        c->SetLogy();

    const string plot_title =
        Form(
            "Sensor %s - T = 21 #circC - all annealing phases;"
            "Voltage (V);Current (A)",
            sensor.c_str()
        );

    create_frame(
        c,
        fxmin,
        fymin,
        fxmax,
        fymax,
        plot_title
    );

    TLegend* leg =
        new TLegend(
            0.14,
            0.56,
            0.46,
            0.88
        );

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.026);

    for (const string& phase : phases) {

        auto it =
            curves.find(phase);

        if (it == curves.end())
            continue;

        const int color =
            phase_color(phase);

        const string phase_file =
            filename_safe(phase);

        TGraphErrors* graph =
            make_graph(
                it->second,
                Form(
                    "gr_set3_%s_%s_T21",
                    sensor.c_str(),
                    phase_file.c_str()
                ),
                color,
                20,
                1.05
            );

        if (!graph)
            continue;

        graph->Draw("PZ SAME");

        leg->AddEntry(
            graph,
            pretty_phase(phase).c_str(),
            "pe"
        );
    }

    leg->Draw();

    const string outbase =
        output_dir
        + "/"
        + sensor
        + "_T21_all_annealing_phases_IV";

    save_canvas(
        c,
        outbase
    );

    delete c;

    return 1;
}

} // namespace


// ============================================================================
// Main macro
// ============================================================================

void IV_scan_plots(const char* base_dir = ".")
{
    gROOT->SetBatch(kTRUE);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Keep scientific notation compact on the current axis.
    TGaxis::SetMaxDigits(3);

    const string base =
        base_dir ? base_dir : ".";

    const string file_A1 =
        base + "/A1_IV_scan.csv";

    const string file_B1 =
        base + "/B1_IV_scan.csv";

    const string plot_dir =
        base + "/IV_plots";

    const string set1_dir =
        plot_dir + "/1_T21_single_phase";

    const string set2_dir =
        plot_dir + "/2_all_temperatures";

    const string set3_dir =
        plot_dir + "/3_T21_all_phases";

    gSystem->mkdir(
        plot_dir.c_str(),
        true
    );

    gSystem->mkdir(
        set1_dir.c_str(),
        true
    );

    gSystem->mkdir(
        set2_dir.c_str(),
        true
    );

    gSystem->mkdir(
        set3_dir.c_str(),
        true
    );

    struct SensorInput {
        string sensor;
        string filename;
    };

    const vector<SensorInput> inputs = {
        {"A1", file_A1},
        {"B1", file_B1}
    };

    int total_set1 = 0;
    int total_set2 = 0;
    int total_set3 = 0;

    for (const auto& input : inputs) {

        std::cout
            << "\n============================================================"
            << std::endl
            << "Reading "
            << input.sensor
            << ": "
            << input.filename
            << std::endl
            << "============================================================"
            << std::endl;

        const vector<IVPoint> rows =
            read_iv_csv(
                input.filename
            );

        if (rows.empty()) {

            std::cerr
                << "WARNING: no valid rows for sensor "
                << input.sensor
                << ". Skipping."
                << std::endl;

            continue;
        }

        const vector<string> phases =
            get_phases(
                rows
            );

        std::cout
            << "Read "
            << rows.size()
            << " IV points and found "
            << phases.size()
            << " phases."
            << std::endl;

        total_set1 +=
            produce_set1(
                input.sensor,
                rows,
                phases,
                set1_dir
            );

        total_set2 +=
            produce_set2(
                input.sensor,
                rows,
                phases,
                set2_dir
            );

        total_set3 +=
            produce_set3(
                input.sensor,
                rows,
                phases,
                set3_dir
            );
    }

    std::cout
        << "\n============================================================"
        << std::endl
        << "IV plot production completed"
        << std::endl
        << "Set 1 canvases: "
        << total_set1
        << std::endl
        << "Set 2 canvases: "
        << total_set2
        << std::endl
        << "Set 3 canvases: "
        << total_set3
        << std::endl
        << "Output directory: "
        << plot_dir
        << std::endl
        << "============================================================"
        << std::endl;
}
