//root -l 'lum_vs_V.C("run=20250826-221930_T=20.csv")'
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

using namespace std;

struct Row {
    double x, y;
    double luminosity, error;
    double T, v, v_fin;
    int spot;
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
        r.x = stod(f[0]);
        r.y = stod(f[1]);
        r.luminosity = stod(f[2]);
        r.error = stod(f[3]);
        r.T = stod(f[4]);
        r.v = stod(f[5]);
        r.spot = stoi(f[6]);
        r.v_fin = stod(f[7]);

        data.push_back(r);
    }

    return data;
}

void lum_vs_V(const char* filename = "data.csv") {

    gStyle->SetOptFit(1111);

    vector<Row> data = read_csv(filename);

    map<int, vector<Row>> spots;

    for (auto &r : data) {
        spots[r.spot].push_back(r);
    }

    int nspots = spots.size();

    int ncols = ceil(sqrt(nspots));
    int nrows = ceil((double)nspots / ncols);

    // ============================================================
    // CANVAS 1: luminosity vs v_fin per ogni spot
    // ============================================================

    TCanvas* c1 = new TCanvas("c1", "Luminosity vs overvoltage", 1600, 1000);
    c1->Divide(ncols, nrows);

    for (auto &entry : spots) {

        int spot = entry.first;
        vector<Row> rows = entry.second;

        sort(rows.begin(), rows.end(),
             [](const Row& a, const Row& b) {
                 return a.v_fin < b.v_fin;
             });

        int n = rows.size();

        vector<double> x(n), y(n), ex(n), ey(n);

        for (int i = 0; i < n; i++) {
            x[i]  = rows[i].v_fin;
            y[i]  = rows[i].luminosity;
            ex[i] = 0.0;
            ey[i] = rows[i].error;
        }

        c1->cd(spot + 1);

        TGraphErrors* gr = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());

        gr->SetTitle(Form("Spot %d", spot));
        gr->GetXaxis()->SetTitle("overvoltage");
        gr->GetYaxis()->SetTitle("luminosity");

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);

        gr->Draw("AP");

        TLegend* leg = new TLegend(0.15, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(gr, Form("x = %.2f, y = %.2f", rows[0].x, rows[0].y), "lep");
        leg->AddEntry((TObject*)0, Form("T = %.1f", rows[0].T), "");
        leg->Draw();
    }

    c1->SaveAs("canvas_luminosity_vs_overvoltage.png");

    // ============================================================
    // CANVAS 2: rapporto L(v_fin=8) / L(v_fin=2)
    // ============================================================

    vector<double> xspot;
    vector<double> ratio;
    vector<double> ex;
    vector<double> eratio;

    for (auto &entry : spots) {

        int spot = entry.first;
        vector<Row> rows = entry.second;

        double L2 = 0, eL2 = 0;
        double L8 = 0, eL8 = 0;

        bool found2 = false;
        bool found8 = false;

        for (auto &r : rows) {
            if (fabs(r.v_fin - 2.0) < 1e-6) {
                L2 = r.luminosity;
                eL2 = r.error;
                found2 = true;
            }

            if (fabs(r.v_fin - 8.0) < 1e-6) {
                L8 = r.luminosity;
                eL8 = r.error;
                found8 = true;
            }
        }

        if (!found2 || !found8 || L2 == 0) {
            cerr << "Spot " << spot
                 << ": impossibile calcolare il rapporto" << endl;
            continue;
        }

        double R = L2 / L8;

        double eR = R * sqrt(pow(eL8 / L8, 2) + pow(eL2 / L2, 2));

        xspot.push_back(spot);
        ratio.push_back(R);
        ex.push_back(0.0);
        eratio.push_back(eR);
    }

    int n = xspot.size();

    TCanvas* c2 = new TCanvas("c2", "Luminosity variation", 1000, 700);

    TGraphErrors* gr_ratio = new TGraphErrors(
        n,
        xspot.data(),
        ratio.data(),
        ex.data(),
        eratio.data()
    );

    gr_ratio->SetTitle("Luminosity variation: L(v_{fin}=2.0) / L(v_{fin}=8.0)");
    gr_ratio->GetXaxis()->SetTitle("spot");
    gr_ratio->GetYaxis()->SetTitle("L(2.0) / L(8.0)");

    gr_ratio->SetMarkerStyle(20);
    gr_ratio->SetMarkerSize(1.2);
    gr_ratio->SetLineWidth(2);

    gr_ratio->Draw("AP");

    TF1* fit = new TF1("fit", "pol0", -0.5, nspots - 0.5);
    gr_ratio->Fit(fit, "R");

    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();

    TLegend* leg2 = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(gr_ratio, "Data", "lep");
    leg2->AddEntry(fit, "Fit costante", "l");
    leg2->AddEntry((TObject*)0, Form("#chi^{2} = %.2f", chi2), "");
    leg2->AddEntry((TObject*)0, Form("ndf = %d", ndf), "");
    leg2->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2 / ndf), "");
    leg2->Draw();

    //c2->SaveAs("canvas_luminosity_ratio_fit.png");
}