//root -l 'lum_vs_V_fit.C("run=202508_scanV.csv")'

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
    double T, v;
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

        if (f.size() < 7) continue;

        Row r;
        r.spot      = stoi(f[0]);
        r.x         = stod(f[1]);
        r.y         = stod(f[2]);
        r.luminosity= stod(f[3]);
        r.error     = stod(f[4]);
        r.T         = stod(f[5]);
        r.v         = stod(f[6]);

        data.push_back(r);
    }

    return data;
}

void lum_vs_V_fit(const char* filename = "data.csv") {

    gStyle->SetOptFit(111);

    vector<Row> data = read_csv(filename);

    map<int, vector<Row>> spots;

    for (auto &r : data) {
        spots[r.spot].push_back(r);
    }

    int nspots = spots.size();

    int ncols = ceil(sqrt(nspots));
    int nrows = ceil((double)nspots / ncols);

    //Variables for the fit
    vector<int> spot_ids;
    vector<double> A_vals, A_errs;
    vector<double> B_vals, B_errs;
    vector<double> chi2ndf_vals;

    // ============================================================
    // CANVAS 6: luminosity vs v per ogni spot
    // ============================================================

    TCanvas* c6 = new TCanvas("c6", "Luminosity vs overvoltage", 1600, 1000);
    c6->Divide(ncols, nrows);

    int pad = 1;

    for (auto &entry : spots) {

        int spot = entry.first;
        vector<Row> rows = entry.second;

        sort(rows.begin(), rows.end(),
             [](const Row& a, const Row& b) {
                 return a.v < b.v;
             });

        int n = rows.size();

        vector<double> x(n), y(n), ex(n), ey(n);

        for (int i = 0; i < n; i++) {
            x[i]  = rows[i].v;
            y[i]  = rows[i].luminosity;
            ex[i] = 0.0;
            ey[i] = rows[i].error;
        }

        c6->cd(pad++);

        TGraphErrors* gr = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());

        gr->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f #circC", spot, rows[0].x, rows[0].y, rows[0].T));
        gr->GetXaxis()->SetTitle("Overvoltage (V)");
        gr->GetYaxis()->SetTitle("Luminosity");

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);

        gr->Draw("AP");

        TF1* f_quad = new TF1(Form("f_quad_%d", spot), "[0] * TMath::Power(x,[1])", 0, 10);
        f_quad->SetParameters(1.0, 2.0);

        gr->Fit(f_quad, "Q");

        double A     = f_quad->GetParameter(0);
        double A_err = f_quad->GetParError(0);

        double B     = f_quad->GetParameter(1);
        double B_err = f_quad->GetParError(1);

        double chi2  = f_quad->GetChisquare();
        double ndf   = f_quad->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2/ndf : 0.0;

        // salva risultati
        spot_ids.push_back(spot);

        A_vals.push_back(A);
        A_errs.push_back(A_err);

        B_vals.push_back(B);
        B_errs.push_back(B_err);

        chi2ndf_vals.push_back(chi2ndf);
        /*
        TLegend* leg = new TLegend(0.15, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        //leg->AddEntry(gr, Form("x = %.2f, y = %.2f", rows[0].x, rows[0].y), "lep");
        //leg->AddEntry((TObject*)0, Form("T = %.1f", rows[0].T), "");
        leg->Draw(); */
    }
    //c6->SaveAs("canvas_luminosity_vs_overvoltage.png");

    // ============================================================
    // CANVAS 7: luminosity vs V_over for a single spot with fit
    // ============================================================
    TCanvas* c7 = new TCanvas("c7", "Spot 2 - Luminosity vs overvoltage", 800, 600);

    // prendi solo spot 2
    int spot_2 = 2;
    vector<Row> rows_2 = spots[spot_2];

    // ordina per overvoltage
    sort(rows_2.begin(), rows_2.end(),
        [](const Row& a, const Row& b) {
            return a.v < b.v;
        });

    int n = rows_2.size();

    vector<double> x_2(n), y_2(n), ex_2(n), ey_2(n);

    for (int i = 0; i < n; i++) {
        x_2[i]  = rows_2[i].v;
        y_2[i]  = rows_2[i].luminosity;
        ex_2[i] = 0.0;
        ey_2[i] = rows_2[i].error;
    }

    TGraphErrors* gr = new TGraphErrors(n, x_2.data(), y_2.data(), ex_2.data(), ey_2.data());

    gr->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f #circC",
                    spot_2, rows_2[0].x, rows_2[0].y, rows_2[0].T));

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
    TF1* f_quad = new TF1("f_quad_2", "[0] * TMath::Power(x,[1])", 0, 10);
    f_quad->SetParameters(1.0, 2.0);
    gr->Fit(f_quad, "Q");


    gr->Draw("AP");

   

   // int n = spot_ids.size();

    vector<double> x(n), ex(n, 0.0);
    for (int i = 0; i < n; i++) x[i] = spot_ids[i];
    // ============================================================
    // CANVAS 8: A vs spot ID
    // ============================================================

   

    TCanvas* c8 = new TCanvas("c8", "A vs spot", 800, 600);

    TGraphErrors* grA = new TGraphErrors(n, x.data(), A_vals.data(), ex.data(), A_errs.data());

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
    TCanvas* c9 = new TCanvas("c9", "B vs spot", 800, 600);

    TGraphErrors* grB = new TGraphErrors(n, x.data(), B_vals.data(), ex.data(), B_errs.data());

    grB->SetTitle("Power law exponent:B;Spot;B");
    grB->SetMarkerStyle(21);
    grB->Draw("AP");
    grB->GetXaxis()->SetLimits(-1, nspots + 0.5);

    // Linear fit B vs spot ID
    TF1* fit_B = new TF1("fit_B", "[0]", 0, n);
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



    // ============================================================
    // CANVAS 9: B vs spot ID con linea B=2
    // ============================================================
    TCanvas* c9 = new TCanvas("c9", "B vs spot", 800, 600);

    TGraphErrors* grB = new TGraphErrors(n, x.data(), B_vals.data(), ex.data(), B_errs.data());

    grB->SetTitle("Power law exponent:B;Spot;B");
    grB->SetMarkerStyle(21);
    grB->Draw("AP");
    grB->GetXaxis()->SetLimits(-1, nspots + 0.5);

    // Linea B=2
    TF1* line_B = new TF1("line_B", "2", 0, n);
    line_B->SetLineColor(kRed+1);
    line_B->SetLineStyle(2);
    line_B->Draw("SAME");

    // Legend
    TLegend* leg9 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg9->SetBorderSize(0);
    leg9->SetFillStyle(0);
    leg9->AddEntry(grB, "Power law exponent", "lep");
    leg9->AddEntry(line_B, "Line: B = 2", "l");
    //leg9->AddEntry((TObject*)0, Form("p0 = %.2f #pm %.2f", fit_result_B, fit_error_B), "");
    //leg9->AddEntry(fit_B, Form("fit: #chi^{2}/ndf = %.2f/%d", chi2_B, ndf_B), "l");
    leg9->Draw();

    // ============================================================
    // CANVAS 10: chi square vs spot ID
    // ============================================================

    TCanvas* c10 = new TCanvas("c10", "Chi2/ndf vs spot", 800, 600);

    TGraph* grChi2 = new TGraph(n, x.data(), chi2ndf_vals.data());

    grChi2->SetTitle("#chi^{2}/ndf;Spot;#chi^{2}/ndf");
    grChi2->SetMarkerStyle(22);
    grChi2->Draw("AP");
    grChi2->GetXaxis()->SetLimits(-1, nspots + 0.5);
    TLegend* leg10 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg10->SetBorderSize(0);
    leg10->SetFillStyle(0);
    leg10->AddEntry(grChi2, "Qualit#grave{a} del fit", "p");
    leg10->Draw();

    //c10->SaveAs("canvas_chi2_vs_spot.png");
}