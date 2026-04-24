// root -l 'lum_vs_T.C("run=20250826-035846_v=7.csv")'

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

struct SpotData {
    double x;
    double y;
    vector<double> T;
    vector<double> luminosity;
    vector<double> error;
    double overvoltage;
};

void lum_vs_T(const char* csvfile = "data.csv") {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    ifstream fin(csvfile);
    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file " << csvfile << endl;
        return;
    }

    string line;

    // Salta header
    getline(fin, line);

    map<int, SpotData> spots;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string token;

        double x, y, luminosity, error, T, overvoltage;
        int spot;

        getline(ss, token, ','); x = stod(token);
        getline(ss, token, ','); y = stod(token);
        getline(ss, token, ','); luminosity = stod(token);
        getline(ss, token, ','); error = stod(token);
        getline(ss, token, ','); T = stod(token);
        getline(ss, token, ','); overvoltage = stod(token);
        getline(ss, token, ','); spot = stoi(token);

        spots[spot].x = x;
        spots[spot].y = y;
        spots[spot].T.push_back(T);
        spots[spot].luminosity.push_back(luminosity);
        spots[spot].error.push_back(error);
        spots[spot].overvoltage = overvoltage;
    }

    fin.close();

    int nspots = spots.size();
    if (nspots == 0) {
        cerr << "Nessuno spot trovato." << endl;
        return;
    }

    int ncols = ceil(sqrt(nspots));
    int nrows = ceil((double)nspots / ncols);

    vector<int> colors = {
        kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2,
        kViolet, kPink+7, kAzure+2, kSpring+5, kTeal+3, kBlack
    };

    // =========================================================
    // CANVAS 1: un grafico per spot
    // =========================================================
    TCanvas* c1 = new TCanvas("c1", "Grafici per spot", 1600, 900);
    c1->Divide(ncols, nrows);

    // =========================================================
    // CANVAS 2: tutti i grafici insieme
    // =========================================================
    TCanvas* c2 = new TCanvas("c2", "Tutti gli spot insieme", 1200, 800);
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg_all = new TLegend(0.72, 0.15, 0.90, 0.88);
    leg_all->SetBorderSize(0);
    leg_all->SetFillStyle(0);

    // =========================================================
    // CANVAS 4: rapporto L(16)/L(22) vs spot
    // =========================================================
    TCanvas* c4 = new TCanvas("c4", "Rapporto L(16)/L(22) per spot", 1200, 800);

    vector<double> spot_x;
    vector<double> ratio_16_22;
    vector<double> ratio_16_22_err;
    vector<double> spot_ex;

    int ipad = 1;
    int icolor = 0;

    for (auto& entry : spots) {
        int spot = entry.first;
        SpotData& d = entry.second;

        int n = d.T.size();
        if (n == 0) continue;

        vector<double> ex(n, 0.0);

        // Ordino per T crescente
        vector<int> idx(n);
        for (int i = 0; i < n; i++) idx[i] = i;

        sort(idx.begin(), idx.end(), [&](int a, int b) {
            return d.T[a] < d.T[b];
        });

        vector<double> T_sorted, L_sorted, E_sorted;
        for (int i : idx) {
            T_sorted.push_back(d.T[i]);
            L_sorted.push_back(d.luminosity[i]);
            E_sorted.push_back(d.error[i]);
}

        int color = colors[icolor % colors.size()];
        icolor++;

        // =====================================================
        // CANVAS 1
        // =====================================================
        TGraphErrors* gr = new TGraphErrors(
            n,
            T_sorted.data(),
            L_sorted.data(),
            ex.data(),
            E_sorted.data()
        );

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);

        c1->cd(ipad);
        ipad++;

        gr->SetTitle(Form("Spot %d;T;Luminosity", spot));
        gr->Draw("AP");

        gr->GetXaxis()->SetLimits(0.0, 25.0);
        gr->SetMinimum(0.0);

        double ymax = *max_element(L_sorted.begin(), L_sorted.end());
        double emax = *max_element(E_sorted.begin(), E_sorted.end());
        gr->SetMaximum(ymax + emax + 0.15 * ymax);

        TLegend* leg = new TLegend(0.15, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr, Form("x = %.2f, y = %.2f", d.x, d.y), "lep");

        TString v_string = Form("Overvoltage = %.1f", d.overvoltage);
        leg->AddEntry((TObject*)0, v_string.Data(), "");

        leg->Draw();

        gPad->SetGrid();

        // =====================================================
        // CANVAS 2
        // =====================================================
        TGraphErrors* gr2 = new TGraphErrors(
            n,
            T_sorted.data(),
            L_sorted.data(),
            ex.data(),
            E_sorted.data()
        );

        gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(1.0);
        gr2->SetMarkerColor(color);
        gr2->SetLineColor(color);
        gr2->SetLineWidth(2);

        mg->Add(gr2, "PL");
        leg_all->AddEntry(gr2, Form("Spot %d", spot), "lp");

        // =====================================================
        // CANVAS 4 : L(16)/L(22) per spot
        // =====================================================
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
            double R16 = L16 / L22;
            double ER16 = R16 * sqrt( pow(E16/L16, 2) + pow(E22/L22, 2) );

            spot_x.push_back((double)spot);
            spot_ex.push_back(0.0);
            ratio_16_22.push_back(R16);
            ratio_16_22_err.push_back(ER16);
        }
    }

    // =========================================================
    // Disegno canvas 2
    // =========================================================
    c2->cd();
    mg->SetTitle("Luminosity vs T per tutti gli spot;T;Luminosity");
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(0.0, 25.0);
    mg->SetMinimum(0.0);
    gPad->Modified();
    gPad->Update();
    leg_all->Draw();
    gPad->SetGrid();

    // =========================================================
    // Disegno canvas 4 + fit
    // =========================================================
    c4->cd();

    if (!spot_x.empty()) {
        TGraphErrors* gr_spot_ratio = new TGraphErrors(
            spot_x.size(),
            spot_x.data(),
            ratio_16_22.data(),
            spot_ex.data(),
            ratio_16_22_err.data()
        );

        gr_spot_ratio->SetTitle("Rapporto L(16)/L(22) per spot;Spot;L(16)/L(22)");
        gr_spot_ratio->SetMarkerStyle(20);
        gr_spot_ratio->SetMarkerSize(1.2);
        gr_spot_ratio->SetLineWidth(2);

        gr_spot_ratio->Draw("AP");
        gr_spot_ratio->GetXaxis()->SetLimits(-0.5, nspots - 0.5);
        gr_spot_ratio->SetMinimum(0.0);
        gr_spot_ratio->SetMaximum(1.2);

        TF1* fit_const = new TF1("fit_const", "[0]", -0.5, nspots - 0.5);
        fit_const->SetParameter(0, 0.8);
        fit_const->SetLineColor(kRed+1);
        fit_const->SetLineWidth(2);

        gr_spot_ratio->Fit(fit_const, "R");

        gPad->SetGrid();
    }

    c1->Update();
    c2->Update();
    c4->Update();
}