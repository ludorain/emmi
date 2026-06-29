//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_annealed_75_25_run=20260530-001924_total_global_ID.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TObject.h"

using namespace std;

struct Row {
    int spot;
    double x, y;
    double luminosity, error;
    double T;
    double v;       // tensione assoluta
    double v_fin;   // overvoltage finale
};

vector<Row> read_csv(const char* filename) {
    vector<Row> data;
    ifstream fin(filename);

    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire " << filename << endl;
        return data;
    }

    string line;
    getline(fin, line); // header

    while (getline(fin, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string field;
        vector<string> f;

        while (getline(ss, field, ',')) {
            f.push_back(field);
        }

        if (f.size() < 8) continue;

        Row r;
        r.spot       = stoi(f[0]);
        r.x          = stod(f[1]);
        r.y          = stod(f[2]);
        r.luminosity = stod(f[3]);
        r.error      = stod(f[4]);
        r.T          = stod(f[5]);
        r.v          = stod(f[6]);
        r.v_fin      = stod(f[7]);

        data.push_back(r);
    }

    return data;
}

void lum_vs_V_fit_irradiated(const char* filename = "data.csv") {

    gStyle->SetOptFit(111);

    vector<Row> data = read_csv(filename);

    map<int, vector<Row>> spots;

    for (auto &r : data) {
        spots[r.spot].push_back(r);
    }

    int nspots = spots.size();
    
    //Parametri display canvas
    int max_spot_display = 8;
    int ncols = 4;
    int nrows = 2;

    //Variables for the fit
    vector<int> spot_ids;
    vector<double> A_vals, A_errs;
    vector<double> B_vals, B_errs;
    vector<double> chi2ndf_vals;

    // ============================================================
    // CANVAS 5: luminosity vs spot at max v_fin
    // ============================================================

    // trova il massimo valore di v_fin presente nel file
    double vmax = -1e9;

    for (const auto& r : data) {
        if (r.v_fin > vmax)
            vmax = r.v_fin;
    }

    vector<double> x_spot;
    vector<double> y_lum;
    vector<double> ey_lum;

    for (auto &entry : spots) {

        int spot = entry.first;
        const vector<Row>& rows = entry.second;

        bool found = false;

        for (const auto& r : rows) {

            // tolleranza per confronti floating point
            if (fabs(r.v_fin - vmax) < 1e-6) {

                x_spot.push_back(spot);
                y_lum.push_back(r.luminosity);
                ey_lum.push_back(r.error);

                found = true;
                break;
            }
        }

        if (!found) {
            cout << "Spot " << spot
                << " non contiene v_fin = "
                << vmax << endl;
        }
    }

    int nmax = x_spot.size();

    vector<double> ex_spot(nmax, 0.0);

    TCanvas* c5 = new TCanvas("c5",
                            Form("Luminosity at v_fin = %.2f V", vmax),
                            1800,
                            600);

    TGraphErrors* grMax =
        new TGraphErrors(nmax,
                        x_spot.data(),
                        y_lum.data(),
                        ex_spot.data(),
                        ey_lum.data());

    grMax->SetTitle(
        Form("Luminosity at maximum overvoltage (v_{fin}=%.2f V) annealing 75#circC 25h;Spot;Luminosity",
            vmax));

    grMax->SetMarkerStyle(20);
    grMax->SetMarkerSize(1.2);
    grMax->SetLineWidth(2);

    grMax->Draw("AP");

    c5->SaveAs("A1_ann_75_25_lum_vs_v_max_vfin.png");

    // ============================================================
    // CANVAS 5.1: luminosity vs spot at max v_fin - focus y range
    // ============================================================

    // trova il massimo valore di v_fin presente nel file
    double vmax_c51 = -1e9;

    for (const auto& r : data) {
        if (r.v_fin > vmax_c51)
            vmax_c51 = r.v_fin;
    }

    vector<double> x_spot_c51;
    vector<double> y_lum_c51;
    vector<double> ey_lum_c51;

    for (auto &entry_c51 : spots) {

        int spot_c51 = entry_c51.first;
        const vector<Row>& rows_c51 = entry_c51.second;

        bool found_c51 = false;

        for (const auto& r_c51 : rows_c51) {

            // tolleranza per confronti floating point
            if (fabs(r_c51.v_fin - vmax_c51) < 1e-6) {

                x_spot_c51.push_back(spot_c51);
                y_lum_c51.push_back(r_c51.luminosity);
                ey_lum_c51.push_back(r_c51.error);

                found_c51 = true;
                break;
            }
        }

        if (!found_c51) {
            cout << "Spot " << spot_c51
                << " non contiene v_fin = "
                << vmax_c51 << endl;
        }
    }

    int nmax_c51 = x_spot_c51.size();

    vector<double> ex_spot_c51(nmax_c51, 0.0);

    TCanvas* c51 = new TCanvas("c51",
                            Form("Luminosity at v_fin = %.2f V - focus", vmax_c51),
                            1800,
                            600);

    TGraphErrors* grMax_c51 =
        new TGraphErrors(nmax_c51,
                        x_spot_c51.data(),
                        y_lum_c51.data(),
                        ex_spot_c51.data(),
                        ey_lum_c51.data());

    grMax_c51->SetTitle(
        Form("Luminosity at maximum overvoltage (v_{fin}=%.2f V) annealing 75#circC 25h - focus;Spot;Luminosity",
            vmax_c51));

    grMax_c51->SetMarkerStyle(20);
    grMax_c51->SetMarkerSize(1.2);
    grMax_c51->SetLineWidth(2);

    // Focus sull'asse y
    grMax_c51->SetMinimum(-5);
    grMax_c51->SetMaximum(1000);

    grMax_c51->Draw("AP");

    c51->SaveAs("A1_ann_75_25_lum_vs_v_max_vfin_focus.png");


    // ============================================================
    // CANVAS 6: luminosity vs v_fin per i primi 8 spot
    // ============================================================

    TCanvas* c6 = new TCanvas("c6", "Luminosity vs overvoltage", 1600, 1000);
    c6->Divide(4, 2);

    int displayed = 0;

    for (auto &entry : spots) {

    int spot = entry.first;
    vector<Row> rows = entry.second;

    sort(rows.begin(), rows.end(),
         [](const Row& a, const Row& b) {
             return a.v_fin < b.v_fin;
         });

    int n = rows.size();

    vector<double> xv(n), yv(n), exv(n), eyv(n);

    for (int i = 0; i < n; i++) {
        xv[i]  = rows[i].v_fin;
        yv[i]  = rows[i].luminosity;
        exv[i] = 0.0;
        eyv[i] = rows[i].error;
    }

    TGraphErrors* gr = new TGraphErrors(n,xv.data(),yv.data(), exv.data(), eyv.data());

    gr->SetTitle(Form("Spot %d: x=%.2f y=%.2f",
                      spot,
                      rows[0].x,
                      rows[0].y));

    gr->GetXaxis()->SetTitle("Overvoltage (V)");
    gr->GetYaxis()->SetTitle("Luminosity");

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);

    TF1* f_quad =
        new TF1(Form("f_quad_%d", spot),
                "[0]*TMath::Power(x,[1])",
                0,
                10);

    f_quad->SetParameters(1.0, 2.0);

    gr->Fit(f_quad, "Q");

    double A       = f_quad->GetParameter(0);
    double A_err   = f_quad->GetParError(0);

    double B       = f_quad->GetParameter(1);
    double B_err   = f_quad->GetParError(1);

    double chi2    = f_quad->GetChisquare();
    double ndf     = f_quad->GetNDF();

    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    // salva SEMPRE i risultati
    spot_ids.push_back(spot);

    A_vals.push_back(A);
    A_errs.push_back(A_err);

    B_vals.push_back(B);
    B_errs.push_back(B_err);

    chi2ndf_vals.push_back(chi2ndf);

    // disegna solo i primi 8 spot
    if (displayed < 8) {
        c6->cd(displayed + 1);
        gr->Draw("AP");
        displayed++;
    }

    }
    c6->SaveAs("A1_ann_75_25_lum_vs_v.png");

    // ============================================================
    // CANVAS 7: luminosity vs V_over for a single spot with fit
    // ============================================================
    TCanvas* c7 = new TCanvas("c7", "Spot 0 - Luminosity vs overvoltage", 800, 600);

    // prendi solo spot 2
    int spot_0 = 2;
    vector<Row> rows_0 = spots[spot_0];

    // ordina per overvoltage
    sort(rows_0.begin(), rows_0.end(),
        [](const Row& a, const Row& b) {
            return a.v_fin < b.v_fin;
        });

    int n = rows_0.size();

    vector<double> x_0(n), y_0(n), ex_0(n), ey_0(n);

    for (int i = 0; i < n; i++) {
        x_0[i]  = rows_0[i].v_fin;
        y_0[i]  = rows_0[i].luminosity;
        ex_0[i] = 0.0;
        ey_0[i] = rows_0[i].error;
    }

    TGraphErrors* gr = new TGraphErrors(n, x_0.data(), y_0.data(), ex_0.data(), ey_0.data());

    gr->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f #circC annealing 75#circC 25h",
                    spot_0, rows_0[0].x, rows_0[0].y, rows_0[0].T));

    gr->GetXaxis()->SetTitle("Overvoltage (V)");
    gr->GetXaxis()->SetTitleSize(0.045);
    //gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetTitle("Luminosity");
    gr->GetYaxis()->SetTitleSize(0.045);
    //gr->GetYaxis()->SetTitleOffset(1.3);

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(2.0);
    gr->SetLineWidth(2);
    
    // fit
    TF1* f_quad = new TF1("f_quad_0", "[0] * TMath::Power(x,[1])", 0, 10);
    f_quad->SetParameters(1.0, 2.0);
    gr->Fit(f_quad, "Q");


    gr->Draw("AP");
    c7->SaveAs("A1_ann_75_25_lum_vs_v_spot0.png");

   

     int n_fit = spot_ids.size();

    vector<double> x(n_fit), ex(n_fit, 0.0);
    for (int i = 0; i < n_fit; i++) x[i] = spot_ids[i];
    // ============================================================
    // CANVAS 8: A vs spot ID
    // ============================================================
        /*
    TCanvas* c8 = new TCanvas("c8", "A vs spot", 800, 600);

    TGraphErrors* grA = new TGraphErrors(n_fit, x.data(), A_vals.data(), ex.data(), A_errs.data());

    grA->SetTitle("Parametro A;Spot;A");
    grA->SetMarkerStyle(20);
    grA->Draw("AP");
    grA->GetXaxis()->SetLimits(-1, nspots + 0.5);


    // Legend
    TLegend* leg8 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg8->SetBorderSize(0);
    leg8->SetFillStyle(0);
    leg8->Draw();

    //c8->SaveAs("canvas_A_vs_spot.png");
        */


    // ============================================================
    // CANVAS 9: B vs spot ID con fit
    // ============================================================
   /*
    TCanvas* c9 = new TCanvas("c9", "B vs spot", 800, 600);

    TGraphErrors* grB = new TGraphErrors(n_fit, x.data(), B_vals.data(), ex.data(), B_errs.data());

    grB->SetTitle("Power law exponent:B;Spot;B");
    grB->SetMarkerStyle(21);
    grB->Draw("AP");
    grB->GetXaxis()->SetLimits(-1, nspots + 0.5);

    // Linear fit B vs spot ID
    TF1* fit_B = new TF1("fit_B", "[0]", 0, n_fit);
    fit_B->SetLineColor(kRed+1);
    grB->Fit(fit_B, "RQ");
    fit_B->Draw("SAME");

    //Fit results
    double chi2_B = fit_B->GetChisquare();
    int ndf_B = fit_B->GetNDF();
    double fit_result_B = fit_B->GetParameter(0);
    double fit_error_B = fit_B->GetParError(0);


    // Legend
    TLegend* leg9 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->AddEntry(grB, "Power law exponent", "lep");
    leg9->AddEntry(fit_B, "Fit: A*x^B", "lep");
    //leg9->AddEntry((TObject*)0, Form("p0 = %.2f #pm %.2f", fit_result_B, fit_error_B), "");
    //leg9->AddEntry(fit_B, Form("fit: #chi^{2}/ndf = %.2f/%d", chi2_B, ndf_B), "l");
    leg9->Draw();

*/

    // ============================================================
    // CANVAS 9: B vs spot ID con linea B=2
    // ============================================================
    TCanvas* c9 = new TCanvas("c9", "B vs spot", 1800, 600);

    TGraphErrors* grB = new TGraphErrors(n_fit, x.data(), B_vals.data(), ex.data(), B_errs.data());

    grB->SetTitle("Power law exponent:B annealing 75#circC 25h;Spot;B");
    grB->SetMarkerStyle(21);
    grB->Draw("AP");
    grB->GetXaxis()->SetLimits(-1, nspots + 0.5);

    // Linea B=2
    TF1* line_B = new TF1("line_B", "2", 0, n_fit);
    line_B->SetLineColor(kRed+1);
    line_B->SetLineStyle(2);
    line_B->Draw("SAME");

    // Linear fit B vs spot ID
    TF1* fit_B = new TF1("fit_B", "[0]", 0, n_fit);
    fit_B->SetLineColor(kMagenta+1);
    grB->Fit(fit_B, "RQ");
    fit_B->Draw("SAME");

    //Fit results
    double chi2_B = fit_B->GetChisquare();
    int ndf_B = fit_B->GetNDF();
    double fit_result_B = fit_B->GetParameter(0);
    double fit_error_B = fit_B->GetParError(0);


    // Legend
    TLegend* leg9 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->AddEntry(grB, "Power law exponent", "lep");
    leg9->AddEntry(fit_B, "Fit: A*x^B", "lep");
    leg9->AddEntry(line_B, "Line: B = 2", "l");
    leg9->AddEntry((TObject*)0, Form("p0 = %.2f #pm %.2f", fit_result_B, fit_error_B), "");
    //leg9->AddEntry(fit_B, Form("fit: #chi^{2}/ndf = %.2f/%d", chi2_B, ndf_B), "l");

    leg9->Draw();
    c9->SaveAs("A1_ann_75_25_lum_vs_v_B.png");

    // ============================================================
    // CANVAS 10: chi square vs spot ID
    // ============================================================

    TCanvas* c10 = new TCanvas("c10", "Chi2/ndf vs spot", 800, 600);

    TGraph* grChi2 = new TGraph(n_fit, x.data(), chi2ndf_vals.data());

    grChi2->SetTitle("#chi^{2}/ndf;Spot;#chi^{2}/ndf");
    grChi2->SetMarkerStyle(22);
    grChi2->Draw("AP");
    grChi2->GetXaxis()->SetLimits(-1, nspots + 0.5);
    TLegend* leg10 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg10->SetBorderSize(0);
    leg10->SetFillStyle(0);
    leg10->AddEntry(grChi2, "Qualit#grave{a} del fit", "p");
    leg10->Draw();

    //c10->SaveAs("A1_ann_75_25_lum_vs_v_chi2.png");
}