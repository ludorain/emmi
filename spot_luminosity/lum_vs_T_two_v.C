//root -l 'lum_vs_T_two_v.C("run=202508_2overvoltages.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TF1.h"

using namespace std;

struct CurveData {
    double x = 0;
    double y = 0;
    vector<double> T;
    vector<double> luminosity;
    vector<double> error;
};

double ratio_error(double L1, double E1, double L2, double E2) {
    if (L1 <= 0 || L2 <= 0) return 0.0;
    double R = L1 / L2;
    return R * sqrt(pow(E1 / L1, 2) + pow(E2 / L2, 2));
}

void lum_vs_T_two_v(const char* csvfile = "data.csv") {

    gStyle->SetOptStat(0);

    ifstream fin(csvfile);
    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file " << csvfile << endl;
        return;
    }

    string line;
    getline(fin, line); // header

    map<int, map<double, CurveData>> data;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string token;

        double x, y, luminosity, error, T, v;
        int spot;

        getline(ss, token, ','); x = stod(token);
        getline(ss, token, ','); y = stod(token);
        getline(ss, token, ','); luminosity = stod(token);
        getline(ss, token, ','); error = stod(token);
        getline(ss, token, ','); T = stod(token);
        getline(ss, token, ','); v = stod(token);
        getline(ss, token, ','); spot = stoi(token);

        data[spot][v].x = x;
        data[spot][v].y = y;
        data[spot][v].T.push_back(T);
        data[spot][v].luminosity.push_back(luminosity);
        data[spot][v].error.push_back(error);
    }

    fin.close();

    int nspots = data.size();
    int ncols = ceil(sqrt(nspots));
    int nrows = ceil((double)nspots / ncols);

    int color_v7 = kBlue + 1;
    int color_v9 = kRed + 1;

    // =========================================================
    // CANVAS 1: due curve per ogni spot
    // =========================================================
    TCanvas* c1 = new TCanvas("c1", "Luminosity vs T per spot", 1800, 1000);
    c1->Divide(ncols, nrows);

    // =========================================================
    // CANVAS 2: L(16)/L(22) vs spot per v=7 e v=9
    // =========================================================
    TCanvas* c2 = new TCanvas("c2", "Rapporto L(16)/L(22)", 1200, 800);

    vector<double> spot_v7, ratio_v7, err_ratio_v7, ex_v7;
    vector<double> spot_v9, ratio_v9, err_ratio_v9, ex_v9;

    int ipad = 1;

    for (auto& spot_entry : data) {
        int spot = spot_entry.first;

        c1->cd(ipad);
        ipad++;

        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.13, 0.62, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        double ymax_global = 0.0;

        for (double v_value : {7.0, 9.0}) {

            if (spot_entry.second.find(v_value) == spot_entry.second.end()) continue;

            CurveData& d = spot_entry.second[v_value];
            int n = d.T.size();
            if (n == 0) continue;

            vector<int> idx(n);
            for (int i = 0; i < n; i++) idx[i] = i;

            sort(idx.begin(), idx.end(), [&](int a, int b) {
                return d.T[a] < d.T[b];
            });

            vector<double> T_sorted, L_sorted, E_sorted, ex;
            for (int i : idx) {
                T_sorted.push_back(d.T[i]);
                L_sorted.push_back(d.luminosity[i]);
                E_sorted.push_back(d.error[i]);
                ex.push_back(0.0);
            }

            int color = (fabs(v_value - 7.0) < 1e-6) ? color_v7 : color_v9;

            TGraphErrors* gr = new TGraphErrors(
                n,
                T_sorted.data(),
                L_sorted.data(),
                ex.data(),
                E_sorted.data()
            );

            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.9);
            gr->SetMarkerColor(color);
            gr->SetLineColor(color);
            gr->SetLineWidth(2);

            mg->Add(gr, "PL");

            leg->AddEntry(
                gr,
                Form("v = %.0f, x = %.2f, y = %.2f", v_value, d.x, d.y),
                "lep"
            );

            double ymax = *max_element(L_sorted.begin(), L_sorted.end());
            if (ymax > ymax_global) ymax_global = ymax;

            // =================================================
            // Rapporto L(16)/L(22) per canvas 2
            // =================================================
            double L16 = -1.0, E16 = 0.0;
            double L22 = -1.0, E22 = 0.0;

            for (int i = 0; i < n; i++) {
                if (fabs(T_sorted[i] - 16.0) < 1e-6) {
                    L16 = L_sorted[i];
                    E16 = E_sorted[i];
                }
                if (fabs(T_sorted[i] - 22.0) < 1e-6) {
                    L22 = L_sorted[i];
                    E22 = E_sorted[i];
                }
            }

            if (L16 > 0.0 && L22 > 0.0) {
                double R = L16 / L22;
                double ER = ratio_error(L16, E16, L22, E22);

                if (fabs(v_value - 7.0) < 1e-6) {
                    spot_v7.push_back((double)spot);
                    ratio_v7.push_back(R);
                    err_ratio_v7.push_back(ER);
                    ex_v7.push_back(0.0);
                } else {
                    spot_v9.push_back((double)spot);
                    ratio_v9.push_back(R);
                    err_ratio_v9.push_back(ER);
                    ex_v9.push_back(0.0);
                }
            }
        }

        mg->SetTitle(Form("Spot %d;T;Luminosity", spot));
        mg->Draw("A");

        mg->GetXaxis()->SetLimits(0.0, 25.0);
        mg->SetMinimum(0.0);
        mg->SetMaximum(ymax_global * 1.20);

        leg->Draw();
        gPad->SetGrid();
    }

    // =========================================================
    // Disegno CANVAS 2
    // =========================================================
    c2->cd();

    TMultiGraph* mg_ratio = new TMultiGraph();

    TGraphErrors* gr_ratio_v7 = new TGraphErrors(
        spot_v7.size(),
        spot_v7.data(),
        ratio_v7.data(),
        ex_v7.data(),
        err_ratio_v7.data()
    );

    gr_ratio_v7->SetMarkerStyle(20);
    gr_ratio_v7->SetMarkerSize(1.2);
    gr_ratio_v7->SetMarkerColor(color_v7);
    gr_ratio_v7->SetLineColor(color_v7);
    gr_ratio_v7->SetLineWidth(2);

    TGraphErrors* gr_ratio_v9 = new TGraphErrors(
        spot_v9.size(),
        spot_v9.data(),
        ratio_v9.data(),
        ex_v9.data(),
        err_ratio_v9.data()
    );

    gr_ratio_v9->SetMarkerStyle(21);
    gr_ratio_v9->SetMarkerSize(1.2);
    gr_ratio_v9->SetMarkerColor(color_v9);
    gr_ratio_v9->SetLineColor(color_v9);
    gr_ratio_v9->SetLineWidth(2);

    // Solo punti, niente linea tra i punti
    mg_ratio->Add(gr_ratio_v7, "P");
    mg_ratio->Add(gr_ratio_v9, "P");

    mg_ratio->SetTitle("Rapporto L(16)/L(22) per spot;Spot;L(16)/L(22)");
    mg_ratio->Draw("A");

    mg_ratio->GetXaxis()->SetLimits(-0.5, nspots - 0.5);
    mg_ratio->SetMinimum(0.0);
    mg_ratio->SetMaximum(1.2);

    gPad->Modified();
    gPad->Update();

    // Fit costante per v = 7
    TF1* fit_v7 = new TF1("fit_v7", "[0]", -0.5, nspots - 0.5);
    fit_v7->SetParameter(0, 0.8);
    fit_v7->SetLineColor(color_v7);
    fit_v7->SetLineWidth(2);

    // Fit costante per v = 9
    TF1* fit_v9 = new TF1("fit_v9", "[0]", -0.5, nspots - 0.5);
    fit_v9->SetParameter(0, 0.8);
    fit_v9->SetLineColor(color_v9);
    fit_v9->SetLineWidth(2);

    // Fit con opzione R sul range, e + per non cancellare oggetti già disegnati
    gr_ratio_v7->Fit(fit_v7, "R+");
    gr_ratio_v9->Fit(fit_v9, "R+");

    // Stampa valori fit su terminale
    cout << endl;
    cout << "===== FIT RAPPORTO L(16)/L(22) =====" << endl;
    cout << "v = 7 : " 
        << fit_v7->GetParameter(0) 
        << " +/- " 
        << fit_v7->GetParError(0) 
        << endl;

    cout << "v = 9 : " 
        << fit_v9->GetParameter(0) 
        << " +/- " 
        << fit_v9->GetParError(0) 
        << endl;

    // Legenda con anche i valori dei fit
    TLegend* leg_ratio = new TLegend(0.60, 0.68, 0.88, 0.88);
    leg_ratio->SetBorderSize(0);
    leg_ratio->SetFillStyle(0);

    leg_ratio->AddEntry(gr_ratio_v7, "v = 7", "lep");
    leg_ratio->AddEntry(
        fit_v7,
        Form("Fit v=7: %.3f #pm %.3f", 
            fit_v7->GetParameter(0), 
            fit_v7->GetParError(0)),
        "l"
    );

    leg_ratio->AddEntry(gr_ratio_v9, "v = 9", "lep");
    leg_ratio->AddEntry(
        fit_v9,
        Form("Fit v=9: %.3f #pm %.3f", 
            fit_v9->GetParameter(0), 
            fit_v9->GetParError(0)),
        "l"
    );

    leg_ratio->Draw();

    gPad->SetGrid();
    c1->Update();
    c2->Update();
}