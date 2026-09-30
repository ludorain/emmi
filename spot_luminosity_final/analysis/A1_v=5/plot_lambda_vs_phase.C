//root -l 'plot_lambda_vs_phase.C("A1_v=5_global_lambda_vs_phase.csv","A1","A1_lambda_vs_annealing")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;


// -------------------------------------------------------------
// Simple CSV split
// -------------------------------------------------------------
static vector<string> split_csv(const string& line)
{
    vector<string> fields;
    string field;
    stringstream ss(line);

    while (getline(ss, field, ',')) {
        fields.push_back(field);
    }

    return fields;
}


// -------------------------------------------------------------
// Make annealing phase names nicer for the x axis
//
// before_annealing
// annealing_T=75_h=5
//
// become
//
// Before annealing
// 75 C, 5 h
// -------------------------------------------------------------
static string pretty_phase_name(const string& phase)
{
    if (phase == "before_annealing")
        return "Before annealing";

    const string prefix = "annealing_T=";

    if (phase.find(prefix) == 0) {

        size_t pos_h = phase.find("_h=");

        if (pos_h != string::npos) {

            string T =
                phase.substr(prefix.size(),
                             pos_h - prefix.size());

            string h =
                phase.substr(pos_h + 3);

            return T + "#circC " + h + " h";
        }
    }

    return phase;
}


// -------------------------------------------------------------
// Main macro
// -------------------------------------------------------------
void plot_lambda_vs_phase(
    const char* csv_file,
    const char* sensor,
    const char* output_name = "lambda_vs_phase")
{
    gStyle->SetOptStat(0);

    // ---------------------------------------------------------
    // Open CSV
    // ---------------------------------------------------------

    ifstream fin(csv_file);

    if (!fin.is_open()) {
        cerr << "ERROR: cannot open file "
             << csv_file << endl;
        return;
    }


    // ---------------------------------------------------------
    // Read header
    // ---------------------------------------------------------

    string line;

    if (!getline(fin, line)) {
        cerr << "ERROR: empty CSV file." << endl;
        return;
    }

    vector<string> header = split_csv(line);

    map<string, int> col;

    for (int i = 0; i < (int)header.size(); ++i)
        col[header[i]] = i;


    // Required columns
    if (!col.count("phase") ||
        !col.count("lambda_common") ||
        !col.count("lambda_common_stat_error")) {

        cerr
            << "ERROR: CSV must contain columns:"
            << endl
            << "  phase"
            << endl
            << "  lambda_common"
            << endl
            << "  lambda_common_stat_error"
            << endl;

        return;
    }


    // ---------------------------------------------------------
    // Read data
    // ---------------------------------------------------------

    vector<string> phases;

    vector<double> x;
    vector<double> ex;

    vector<double> lambda;
    vector<double> lambda_err;


    while (getline(fin, line)) {

        if (line.empty())
            continue;

        vector<string> fields = split_csv(line);

        if (fields.size() < header.size())
            continue;

        try {

            string phase =
                fields[col["phase"]];

            double l =
                stod(fields[col["lambda_common"]]);

            double el =
                stod(fields[col["lambda_common_stat_error"]]);


            int i = phases.size();

            phases.push_back(phase);

            x.push_back((double)i);
            ex.push_back(0.0);

            lambda.push_back(l);
            lambda_err.push_back(el);

        }
        catch (...) {

            cerr
                << "WARNING: malformed row skipped:"
                << endl
                << line
                << endl;
        }
    }

    fin.close();


    if (lambda.empty()) {
        cerr << "ERROR: no valid data found." << endl;
        return;
    }


    const int n = lambda.size();


    // ---------------------------------------------------------
    // Canvas
    // ---------------------------------------------------------

    TCanvas* c =
        new TCanvas(
            "c_lambda_vs_phase",
            "Global lambda vs annealing phase",
            1900,
            900
        );

    // More space for horizontal phase labels
    c->SetBottomMargin(0.20);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.10);


    // ---------------------------------------------------------
    // Determine y range
    // ---------------------------------------------------------

    double ymin = 1e99;
    double ymax = -1e99;

    for (int i = 0; i < n; ++i) {

        ymin = std::min(
            ymin,
            lambda[i] - lambda_err[i]
        );

        ymax = std::max(
            ymax,
            lambda[i] + lambda_err[i]
        );
    }

    double yrange = ymax - ymin;

    if (yrange <= 0.0)
        yrange = std::max(0.01, fabs(ymax) * 0.1);

    ymin -= 0.20 * yrange;
    ymax += 0.35 * yrange;


    // ---------------------------------------------------------
    // Frame with categorical x-axis
    // ---------------------------------------------------------

    string title =
        Form(
            "Sensor %s - Simultaneous fit parameter #lambda vs annealing phases",
            sensor
        );

    TH1D* frame =
        new TH1D(
            "frame_lambda_phase",
            Form("%s;Annealing phase;", title.c_str()),
            n,
            -0.5,
            n - 0.5
        );

    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);


    // ---------------------------------------------------------
    // X-axis labels
    // ---------------------------------------------------------

    for (int i = 0; i < n; ++i) {

        string label =
            pretty_phase_name(phases[i]);

        frame->GetXaxis()->SetBinLabel(
            i + 1,
            label.c_str()
        );
    }

    // Horizontal phase labels
    frame->GetXaxis()->LabelsOption("h");

    frame->GetXaxis()->SetLabelSize(0.030);

    // Lower the x-axis title slightly
    frame->GetXaxis()->SetTitleOffset(1.65);
    frame->GetXaxis()->SetTitleSize(0.045);


    // ---------------------------------------------------------
    // Y axis
    //
    // ROOT normally rotates the Y-axis title.
    // Leave the normal title empty and draw an horizontal one
    // manually with TLatex.
    // ---------------------------------------------------------

    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitle("");

    frame->Draw();


    // ---------------------------------------------------------
    // Horizontal Y-axis title
    // ---------------------------------------------------------

    TLatex* ytitle = new TLatex();

    ytitle->SetNDC();
    ytitle->SetTextFont(42);
    ytitle->SetTextSize(0.045);
    ytitle->SetTextAlign(13);

    ytitle->DrawLatex(
        0.015,
        0.92,
        "#lambda_{global}"
    );


    // ---------------------------------------------------------
    // Graph
    // ---------------------------------------------------------

    TGraphErrors* gr =
        new TGraphErrors(
            n,
            x.data(),
            lambda.data(),
            ex.data(),
            lambda_err.data()
        );

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.4);
    gr->SetLineWidth(2);

    gr->Draw("P SAME");


    // ---------------------------------------------------------
    // Legend
    // ---------------------------------------------------------

    TLegend* leg =
        new TLegend(
            0.18,
            0.72,
            0.44,
            0.87
        );

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    leg->AddEntry(
        gr,
        "#lambda_{global} from fit",
        "p"
    );

    leg->AddEntry(
        gr,
        "Statistical uncertainty",
        "e"
    );

    leg->Draw();


    // ---------------------------------------------------------
    // Save
    // ---------------------------------------------------------

    c->Update();

    c->SaveAs(
        Form("%s.png", output_name)
    );

    c->SaveAs(
        Form("%s.pdf", output_name)
    );

    cout
        << "Saved lambda vs annealing phase plot:"
        << endl
        << "  " << output_name << ".png"
        << endl
        << "  " << output_name << ".pdf"
        << endl;
}