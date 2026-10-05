//root -l 'plot_mean_luminosity_vs_phase.C("A1_T=20_v=5_total_luminosity_sum.csv","A1","A1_mean_luminosity_vs_phase")'
//root -l 'plot_mean_luminosity_vs_phase.C("B1_T=20_v=3_total_luminosity_sum.csv","B1","B1_mean_luminosity_vs_phase")'


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
// Data structures
// ============================================================================

struct LuminosityRow {

    string phase;

    double mean_luminosity =
        numeric_limits<double>::quiet_NaN();

    double error =
        numeric_limits<double>::quiet_NaN();
};


struct PhaseInfo {

    string raw;

    bool before = false;

    double temperature =
        numeric_limits<double>::quiet_NaN();

    double hours =
        numeric_limits<double>::quiet_NaN();
};


// ============================================================================
// General utilities
// ============================================================================

static string trim_copy(const string& s)
{
    const char* ws = " \t\n\r\f\v";

    const auto begin =
        s.find_first_not_of(ws);

    if (begin == string::npos)
        return "";

    const auto end =
        s.find_last_not_of(ws);

    return s.substr(
        begin,
        end - begin + 1
    );
}


static vector<string> split_csv_simple(
    const string& line)
{
    vector<string> out;

    stringstream ss(line);

    string item;

    while (getline(ss, item, ','))

        out.push_back(
            trim_copy(item)
        );

    return out;
}


static bool finite_value(double x)
{
    return std::isfinite(x);
}


// ============================================================================
// Phase handling
// ============================================================================

static string canonical_phase(
    const string& phase)
{
    string p =
        trim_copy(phase);

    if (p == "bef_ann" ||
        p == "before_ann" ||
        p == "before_annealing")

        return "before_annealing";

    return p;
}


// Supports:
//
// before_annealing
// annealing_T=75_h=5
// annealing_75_5
//
static PhaseInfo parse_phase(
    const string& phase_input)
{
    const string phase =
        canonical_phase(
            phase_input
        );


    PhaseInfo p;

    p.raw = phase;


    // ------------------------------------------------------------------------
    // Before annealing
    // ------------------------------------------------------------------------

    if (phase == "before_annealing") {

        p.before = true;

        return p;
    }


    // ------------------------------------------------------------------------
    // Convention:
    //
    // annealing_75_5
    // ------------------------------------------------------------------------

    const string simple_prefix =
        "annealing_";


    if (phase.find(simple_prefix) == 0 &&
        phase.find("annealing_T=") != 0) {

        const string rest =
            phase.substr(
                simple_prefix.size()
            );


        const size_t pos =
            rest.find('_');


        if (pos != string::npos) {

            try {

                p.temperature =
                    stod(
                        rest.substr(
                            0,
                            pos
                        )
                    );


                p.hours =
                    stod(
                        rest.substr(
                            pos + 1
                        )
                    );


                return p;

            } catch (...) {

                p.temperature =
                    numeric_limits<double>::quiet_NaN();

                p.hours =
                    numeric_limits<double>::quiet_NaN();
            }
        }
    }


    // ------------------------------------------------------------------------
    // Convention:
    //
    // annealing_T=75_h=5
    // ------------------------------------------------------------------------

    const string old_prefix =
        "annealing_T=";


    if (phase.find(old_prefix) == 0) {

        const size_t pos_h =
            phase.find("_h=");


        if (pos_h != string::npos) {

            try {

                p.temperature =
                    stod(
                        phase.substr(
                            old_prefix.size(),
                            pos_h -
                            old_prefix.size()
                        )
                    );


                p.hours =
                    stod(
                        phase.substr(
                            pos_h + 3
                        )
                    );

            } catch (...) {

                p.temperature =
                    numeric_limits<double>::quiet_NaN();

                p.hours =
                    numeric_limits<double>::quiet_NaN();
            }
        }
    }


    return p;
}


// ============================================================================
// Annealing phase spacing
//
// before annealing -> x = 0
//
// 75 C, 5 h       -> x = 1
// 75 C, 25 h      -> x = 6
//
// 100 C, 5 h      -> x = 7
// 100 C, 25 h     -> x = 12
//
// 125 C, 5 h      -> x = 13
// 125 C, 25 h     -> x = 18
//
// 150 C, 5 h      -> x = 19
// 150 C, 25 h     -> x = 24
//
// ============================================================================

static double phase_step_units(
    const string& phase)
{
    const PhaseInfo p =
        parse_phase(
            phase
        );


    if (p.before)
        return 0.0;


    if (finite_value(p.hours) &&
        p.hours > 0.0)

        return p.hours / 5.0;


    return 1.0;
}


static vector<double> build_phase_positions(
    const vector<string>& phases)
{
    vector<double> x(
        phases.size(),
        0.0
    );


    if (phases.empty())
        return x;


    x[0] = 0.0;


    for (size_t i = 1;
         i < phases.size();
         ++i) {

        double step =
            phase_step_units(
                phases[i]
            );


        if (!(step > 0.0) ||
            !finite_value(step))

            step = 1.0;


        x[i] =
            x[i - 1] +
            step;
    }


    return x;
}


// ============================================================================
// CSV reader
// ============================================================================

static vector<LuminosityRow>
read_luminosity_csv(
    const char* csv_path)
{
    ifstream fin(
        csv_path
    );


    if (!fin.is_open())

        throw runtime_error(
            string(
                "Cannot open CSV file: "
            ) +
            csv_path
        );


    string header;


    if (!getline(fin, header))

        throw runtime_error(
            string(
                "CSV file is empty: "
            ) +
            csv_path
        );


    vector<string> header_fields =
        split_csv_simple(
            header
        );


    map<string, int> col;


    for (int i = 0;
         i < (int)header_fields.size();
         ++i)

        col[
            header_fields[i]
        ] = i;


    const vector<string> required = {

        "phase",
        "mean_luminosity",
        "error"

    };


    for (const auto& name :
         required) {

        if (!col.count(name))

            throw runtime_error(
                "Missing required column in CSV: "
                + name
            );
    }


    vector<LuminosityRow> rows;

    string line;


    while (getline(fin, line)) {

        line =
            trim_copy(
                line
            );


        if (line.empty())
            continue;


        vector<string> fields =
            split_csv_simple(
                line
            );


        const int max_needed =
            max({

                col["phase"],
                col["mean_luminosity"],
                col["error"]

            });


        if ((int)fields.size() <=
            max_needed) {

            cerr
                << "WARNING: malformed row skipped:"
                << endl
                << line
                << endl;

            continue;
        }


        try {

            LuminosityRow r;


            r.phase =
                canonical_phase(
                    fields[
                        col["phase"]
                    ]
                );


            r.mean_luminosity =
                stod(
                    fields[
                        col["mean_luminosity"]
                    ]
                );


            r.error =
                stod(
                    fields[
                        col["error"]
                    ]
                );


            if (!finite_value(
                    r.mean_luminosity
                ) ||
                !finite_value(
                    r.error
                ))

                continue;


            rows.push_back(
                r
            );


        } catch (...) {

            cerr
                << "WARNING: malformed row skipped:"
                << endl
                << line
                << endl;
        }
    }


    if (rows.empty())

        throw runtime_error(
            "No valid rows found in CSV."
        );


    return rows;
}


// ============================================================================
// Bottom annealing-time labels
// ============================================================================

static void draw_time_labels(
    TPad* pad,
    const vector<string>& phases,
    const vector<double>& xpos,
    double xmin,
    double xmax)
{
    if (!pad ||
        phases.empty() ||
        xpos.empty() ||
        !(xmax > xmin))

        return;


    pad->Update();


    const double left =
        pad->GetLeftMargin();

    const double right =
        pad->GetRightMargin();

    const double usable =
        1.0 -
        left -
        right;


    for (size_t i = 0;
         i < phases.size();
         ++i) {

        const double frac =
            (
                xpos[i] -
                xmin
            )
            /
            (
                xmax -
                xmin
            );


        const double x_ndc =
            left +
            frac * usable;


        const PhaseInfo p =
            parse_phase(
                phases[i]
            );


        // No time label for before annealing.
        if (p.before)
            continue;


        string label;


        if (finite_value(
                p.hours
            )) {

            const double hround =
                round(
                    p.hours
                );


            if (fabs(
                    p.hours -
                    hround
                ) < 1e-9)

                label =
                    Form(
                        "%.0f h",
                        p.hours
                    );

            else

                label =
                    Form(
                        "%.1f h",
                        p.hours
                    );

        } else {

            label =
                phases[i];
        }


        TLatex* lab =
            new TLatex();


        lab->SetNDC();

        lab->SetTextFont(42);

        lab->SetTextSize(
            0.032
        );

        lab->SetTextAlign(
            22
        );


        // ------------------------------------------------------------
        // Lower than in the original version because the bottom margin
        // has been reduced.
        // ------------------------------------------------------------

        lab->DrawLatex(
            x_ndc,
            0.10,
            label.c_str()
        );
    }


    // ========================================================================
    // X-axis title
    // ========================================================================

    TLatex* title =
        new TLatex();


    title->SetNDC();

    title->SetTextFont(42);

    title->SetTextSize(
        0.047
    );

    title->SetTextAlign(
        22
    );


    title->DrawLatex(
        0.5,
        0.032,
        "Annealing time"
    );
}


// ============================================================================
// Annealing-temperature labels and vertical separators
// ============================================================================

static void draw_temperature_blocks(
    const vector<string>& phases,
    const vector<double>& xpos,
    double ymin,
    double ymax)
{
    if (phases.size() !=
            xpos.size() ||
        phases.empty())

        return;


    map<
        double,
        vector<double>
    > temp_positions;


    for (size_t i = 0;
         i < phases.size();
         ++i) {

        const PhaseInfo p =
            parse_phase(
                phases[i]
            );


        if (!p.before &&
            finite_value(
                p.temperature
            ))

            temp_positions[
                p.temperature
            ].push_back(
                xpos[i]
            );
    }


    const double yrange =
        ymax -
        ymin;


    const double ytext =
        ymax -
        0.055 * yrange;


    TLatex latex;


    latex.SetTextFont(
        42
    );

    latex.SetTextSize(
        0.034
    );

    latex.SetTextAlign(
        22
    );


    // ========================================================================
    // Before annealing
    // ========================================================================

    latex.SetTextSize(
        0.030
    );


    latex.DrawLatex(
        xpos.front() - 1.0,
        ytext,
        "#splitline{Before}{annealing}"
    );


    latex.SetTextSize(
        0.034
    );


    // ========================================================================
    // Annealing temperatures
    // ========================================================================

    for (const auto& kv :
         temp_positions) {

        if (kv.second.empty())
            continue;


        const double xmin_block =
            *min_element(
                kv.second.begin(),
                kv.second.end()
            );


        const double xmax_block =
            *max_element(
                kv.second.begin(),
                kv.second.end()
            );


        const double xcentre =
            0.5 *
            (
                xmin_block +
                xmax_block
            );


        const double Tround =
            round(
                kv.first
            );


        string tlabel;


        if (fabs(
                kv.first -
                Tround
            ) < 1e-9)

            tlabel =
                Form(
                    "T = %.0f#circC",
                    kv.first
                );

        else

            tlabel =
                Form(
                    "T = %.1f#circC",
                    kv.first
                );


        latex.DrawLatex(
            xcentre,
            ytext,
            tlabel.c_str()
        );
    }


    // ========================================================================
    // Vertical separators
    // ========================================================================

    const vector<double>
    separators = {

        0.5,
        6.5,
        12.5,
        18.5

    };


    for (double xsep :
         separators) {

        TLine* line =
            new TLine(
                xsep,
                ymin,
                xsep,
                ymax
            );


        line->SetLineColor(
            kGray + 1
        );

        line->SetLineStyle(
            2
        );

        line->SetLineWidth(
            1
        );


        line->Draw(
            "same"
        );
    }
}


// ============================================================================
// MAIN MACRO
// ============================================================================

void plot_mean_luminosity_vs_phase(
    const char* csv_path = "mean_luminosity_vs_phase.csv",
    const char* sensor = "A1",
    const char* output_prefix = "mean_luminosity_vs_phase")
{
    gStyle->SetOptStat(
        0
    );


    // ========================================================================
    // Measurement conditions
    // ========================================================================

    const double measurement_temperature =
        20.0;


    double overvoltage =
        numeric_limits<double>::quiet_NaN();


    if (string(sensor) == "A1")

        overvoltage =
            5.0;


    else if (string(sensor) == "B1")

        overvoltage =
            3.0;


    else

        throw runtime_error(
            string(
                "Unknown sensor: "
            ) +
            sensor +
            ". Expected A1 or B1."
        );


    // ========================================================================
    // Plot title
    // ========================================================================

    const string plot_title =
        Form(
            "Sensor %s - Mean luminosity vs annealing phase, "
            "T = %.0f#circC, v_{over} = %.0f V",
            sensor,
            measurement_temperature,
            overvoltage
        );


    // ========================================================================
    // Annealing phases
    // ========================================================================

    const vector<string> phases = {

        "before_annealing",

        "annealing_T=75_h=5",
        "annealing_T=75_h=25",

        "annealing_T=100_h=5",
        "annealing_T=100_h=25",

        "annealing_T=125_h=5",
        "annealing_T=125_h=25",

        "annealing_T=150_h=5",
        "annealing_T=150_h=25"

    };


    // ========================================================================
    // Phase positions
    // ========================================================================

    const vector<double> x =
        build_phase_positions(
            phases
        );


    map<string, double>
    x_phase;


    for (size_t i = 0;
         i < phases.size();
         ++i)

        x_phase[
            phases[i]
        ] =
            x[i];


    // ========================================================================
    // Read CSV
    // ========================================================================

    const vector<LuminosityRow> rows =
        read_luminosity_csv(
            csv_path
        );


    map<
        string,
        pair<double, double>
    > data;


    double ymin_data =
        numeric_limits<double>::infinity();


    double ymax_data =
        -numeric_limits<double>::infinity();


    for (const auto& r :
         rows) {

        if (!x_phase.count(
                r.phase
            ))

            continue;


        data[
            r.phase
        ] = {

            r.mean_luminosity,
            r.error

        };


        ymin_data =
            min(
                ymin_data,
                r.mean_luminosity -
                r.error
            );


        ymax_data =
            max(
                ymax_data,
                r.mean_luminosity +
                r.error
            );
    }


    if (!finite_value(
            ymin_data
        ) ||
        !finite_value(
            ymax_data
        ))

        throw runtime_error(
            "No valid mean_luminosity values found in the CSV."
        );


    // ========================================================================
    // Graph vectors
    // ========================================================================

    vector<double> xv;
    vector<double> yv;

    vector<double> ex;
    vector<double> ey;


    for (size_t i = 0;
         i < phases.size();
         ++i) {

        auto it =
            data.find(
                phases[i]
            );


        if (it ==
            data.end())

            continue;


        xv.push_back(
            x[i]
        );


        yv.push_back(
            it->second.first
        );


        ex.push_back(
            0.0
        );


        ey.push_back(
            it->second.second
        );
    }


    if (xv.empty())

        throw runtime_error(
            "No valid points available for plotting."
        );


    // ========================================================================
    // Y range
    // ========================================================================

    double yrange =
        ymax_data -
        ymin_data;


    if (!(yrange > 0.0))

        yrange =
            max(
                1.0,
                fabs(
                    ymax_data
                ) * 0.10
            );


    const double ymin =
        max(
            0.0,
            ymin_data -
            0.12 * yrange
        );


    const double ymax =
        ymax_data +
        0.25 * yrange;


    // ========================================================================
    // X range
    // ========================================================================

    const int max_unit =
        max(
            0,
            (int)llround(
                x.back()
            )
        );


    const double xmin =
        -2.5;


    const double xmax =
        max_unit +
        0.5;


    const int nbins =
        max_unit +
        1;


    // ========================================================================
    // Canvas
    // ========================================================================

    TCanvas* c =
        new TCanvas(
            "c_mean_luminosity_vs_phase",
            "mean luminosity vs annealing phase",
            2250,
            1000
        );


    // ------------------------------------------------------------------------
    // Reduced from 0.20 to 0.16.
    //
    // This extends the plot area downward and removes part of the large
    // empty gap previously present between the graph and the x-axis title.
    // ------------------------------------------------------------------------

    c->SetBottomMargin(
        0.16
    );


    c->SetLeftMargin(
        0.10
    );


    c->SetRightMargin(
        0.04
    );


    c->SetTopMargin(
        0.10
    );


    c->SetTicks(
        1,
        1
    );


    // ========================================================================
    // Frame
    // ========================================================================

    TH1D* frame =
        new TH1D(
            "frame_mean_luminosity_vs_phase",

            Form(
                "%s;;Mean luminosity",
                plot_title.c_str()
            ),

            nbins,
            xmin,
            xmax
        );


    frame->SetMinimum(
        ymin
    );


    frame->SetMaximum(
        ymax
    );


    // ========================================================================
    // Plot title style
    // ========================================================================

    frame->SetTitleSize(
        0.045
    );


    frame->SetTitleOffset(
        1.0
    );


    // ========================================================================
    // X axis
    // ========================================================================

    frame->GetXaxis()
        ->SetLabelSize(
            0.0
        );


    frame->GetXaxis()
        ->SetTitle(
            ""
        );


    frame->GetXaxis()
        ->SetNdivisions(
            0
        );


    frame->GetXaxis()
        ->SetTickLength(
            0.0
        );


    // ========================================================================
    // Y axis
    // ========================================================================

    frame->GetYaxis()
        ->SetLabelSize(
            0.040
        );


    frame->GetYaxis()
        ->SetTitleSize(
            0.047
        );


    frame->GetYaxis()
        ->SetTitleOffset(
            1.05
        );


    frame->Draw();


    // ========================================================================
    // Manual x ticks
    // ========================================================================

    const double tick_height =
        0.018 *
        (
            ymax -
            ymin
        );


    for (double xpos :
         x) {

        TLine* tick =
            new TLine(
                xpos,
                ymin,
                xpos,
                ymin +
                tick_height
            );


        tick->SetLineColor(
            kBlack
        );


        tick->SetLineWidth(
            1
        );


        tick->Draw(
            "same"
        );
    }


    // ========================================================================
    // Temperature labels + separators
    // ========================================================================

    draw_temperature_blocks(
        phases,
        x,
        ymin,
        ymax
    );


    // ========================================================================
    // Mean luminosity graph
    // ========================================================================

    TGraphErrors* gr =
        new TGraphErrors(
            (int)xv.size(),
            xv.data(),
            yv.data(),
            ex.data(),
            ey.data()
        );


    gr->SetName(
        "gr_mean_luminosity"
    );


    gr->SetLineColor(
        kP6Red
    );


    gr->SetMarkerColor(
        kP6Red
    );


    gr->SetMarkerStyle(
        20
    );


    gr->SetMarkerSize(
        1.45
    );


    gr->SetLineWidth(
        3
    );


    gr->Draw(
        "PL E1 SAME"
    );


    // ========================================================================
    // Legend
    // ========================================================================

    TLegend* leg = new TLegend(0.70,0.60,0.83,0.70);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->AddEntry(gr,"mean luminosity","lep");
    leg->Draw();


    // ========================================================================
    // Bottom labels
    // ========================================================================

    c->Modified();
    c->Update();


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


    // ========================================================================
    // Save
    // ========================================================================

    c->SaveAs(
        Form(
            "%s.png",
            output_prefix
        )
    );


    c->SaveAs(
        Form(
            "%s.pdf",
            output_prefix
        )
    );


    // ========================================================================
    // Terminal output
    // ========================================================================

    cout << endl;


    cout
        << "Sensor: "
        << sensor
        << endl;


    cout
        << "Conditions: T = "
        << measurement_temperature
        << " C, v_over = "
        << overvoltage
        << " V"
        << endl;


    cout
        << "Number of plotted phases: "
        << xv.size()
        << endl;


    cout
        << "Saved plots:"
        << endl;


    cout
        << "  "
        << output_prefix
        << ".png"
        << endl;


    cout
        << "  "
        << output_prefix
        << ".pdf"
        << endl;


    cout << endl;


    cout
        << "Phase x positions:"
        << endl;


    for (size_t i = 0;
         i < phases.size();
         ++i)

        cout
            << "  "
            << phases[i]
            << " -> x = "
            << x[i]
            << endl;


    cout
        << "Spacing rule: "
        << "5 h = +1 unit, "
        << "25 h = +5 units."
        << endl;
}