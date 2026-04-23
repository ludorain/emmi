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
#include "TLatex.h"

using namespace std;

struct SpotData {
    double x;
    double y;
    vector<double> T;
    vector<double> luminosity;
    vector<double> error;
    vector<double> v;
};

void lum_vs_T(const char* csvfile = "data.csv") {

    gStyle->SetOptStat(0);

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

        double x, y, luminosity, error, T, v;
        int spot;

        getline(ss, token, ','); x = stod(token);
        getline(ss, token, ','); y = stod(token);
        getline(ss, token, ','); luminosity = stod(token);
        getline(ss, token, ','); error = stod(token);
        getline(ss, token, ','); T = stod(token);
        getline(ss, token, ','); v = stod(token);
        getline(ss, token, ','); spot = stoi(token);

        spots[spot].x = x;
        spots[spot].y = y;
        spots[spot].T.push_back(T);
        spots[spot].luminosity.push_back(luminosity);
        spots[spot].error.push_back(error);
        spots[spot].v.push_back(v);
    }

    fin.close();

    int nspots = spots.size();
    if (nspots == 0) {
        cerr << "Nessuno spot trovato." << endl;
        return;
    }

    // Calcolo disposizione canvas 1
    int ncols = ceil(sqrt(nspots));
    int nrows = ceil((double)nspots / ncols);

    // Colori ROOT
    vector<int> colors = {
        kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2,
        kViolet, kPink+7, kAzure+2, kSpring+5, kTeal+3, kBlack
    };

    // -----------------------------
    // Canvas 1: un grafico per spot
    // -----------------------------
    TCanvas* c1 = new TCanvas("c1", "Grafici per spot", 1600, 900);
    c1->Divide(ncols, nrows);

    // -----------------------------
    // Canvas 2: tutti insieme
    // -----------------------------
    TCanvas* c2 = new TCanvas("c2", "Tutti gli spot insieme", 1200, 800);

    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg_all = new TLegend(0.72, 0.15, 0.90, 0.88);
    leg_all->SetBorderSize(0);
    leg_all->SetFillStyle(0);

    int ipad = 1;
    int icolor = 0;

    for (auto& entry : spots) {
        int spot = entry.first;
        SpotData& d = entry.second;

        int n = d.T.size();
        if (n == 0) continue;

        vector<double> ex(n, 0.0);

        TGraphErrors* gr = new TGraphErrors(
            n,
            d.T.data(),
            d.luminosity.data(),
            ex.data(),
            d.error.data()
        );

        int color = colors[icolor % colors.size()];
        icolor++;

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color);
        gr->SetLineWidth(2);

        // ==========
        // Canvas 1
        // ==========
        c1->cd(ipad);
        ipad++;

        gr->SetTitle(Form("Spot %d;T;Luminosity", spot));
        gr->Draw("AP");

        // assi da 0
        gr->GetXaxis()->SetLimits(0.0, 25.0);
        gr->SetMinimum(0.0);

        // massimo y un po' più alto del massimo del grafico
        double ymax = *max_element(d.luminosity.begin(), d.luminosity.end());
        double emax = *max_element(d.error.begin(), d.error.end());
        gr->SetMaximum(ymax + emax + 0.15 * ymax);

        gr->GetXaxis()->SetTitleSize(0.05);
        gr->GetYaxis()->SetTitleSize(0.05);
        gr->GetXaxis()->SetLabelSize(0.045);
        gr->GetYaxis()->SetLabelSize(0.045);

        TLegend* leg = new TLegend(0.15, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr, Form("x = %.2f, y = %.2f", d.x, d.y), "lep");

        TString v_string = "v = ";
        for (int i = 0; i < (int)d.v.size(); i++) {
            v_string += Form("%.1f", d.v[i]);
            if (i < (int)d.v.size() - 1) v_string += ", ";
        }
        leg->AddEntry((TObject*)0, v_string.Data(), "");
        leg->Draw();

        gPad->SetGrid();

        // ==========
        // Canvas 2
        // ==========
        TGraphErrors* gr2 = new TGraphErrors(
            n,
            d.T.data(),
            d.luminosity.data(),
            ex.data(),
            d.error.data()
        );

        gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(1.0);
        gr2->SetMarkerColor(color);
        gr2->SetLineColor(color);
        gr2->SetLineWidth(2);

        mg->Add(gr2, "PL");
        leg_all->AddEntry(gr2, Form("Spot %d", spot), "lp");
    }

    c2->cd();
    mg->SetTitle("Luminosity vs T per tutti gli spot;T;Luminosity");
    mg->Draw("A");

    // assi da 0
    mg->GetXaxis()->SetLimits(0.0, 25.0);
    mg->SetMinimum(0.0);

    gPad->Modified();
    gPad->Update();

    leg_all->Draw();
    gPad->SetGrid();

    c1->Update();
    c2->Update();
}