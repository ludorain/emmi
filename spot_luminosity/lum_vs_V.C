//root -l -q 'lum_vs_V.C("run=202508_scanV.csv")'

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
    // CANVAS 1: luminosity vs v per ogni spot
    // ============================================================

    TCanvas* c1 = new TCanvas("c1", "Luminosity vs overvoltage", 1600, 1000);
    c1->Divide(ncols, nrows);

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

        c1->cd(pad++);

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
    // CANVAS 2: rapporto L(v=2) / L(v=8)
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
            if (fabs(r.v - 2.0) < 1e-6) {
                L2 = r.luminosity;
                eL2 = r.error;
                found2 = true;
            }

            if (fabs(r.v - 8.0) < 1e-6) {
                L8 = r.luminosity;
                eL8 = r.error;
                found8 = true;
            }
        }

        if (!found2 || !found8 || L2 == 0 || L8 == 0) {
            cerr << "Spot " << spot
                 << ": impossibile calcolare il rapporto" << endl;
            continue;
        }

        double R = L2 / L8;

        double eR = fabs(R) * sqrt(pow(eL2 / L2, 2) + pow(eL8 / L8, 2));

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

    gr_ratio->SetTitle("Luminosity variation: L(v=2) / L(v=8)");
    gr_ratio->GetXaxis()->SetTitle("spot");
    gr_ratio->GetYaxis()->SetTitle("L(2) / L(8)");

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

    if (ndf > 0)
        leg2->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2 / ndf), "");

    leg2->Draw();

    c2->SaveAs("canvas_luminosity_ratio_fit.png");

    // ============================================================
    // CANVAS 3: luminosity vs v per tutti gli spot insieme
    // ============================================================
    double ymin = 1e9;
    double ymax = -1e9;

    for (auto &entry : spots) {
        for (auto &r : entry.second) {
            if (r.luminosity < ymin) ymin = r.luminosity;
            if (r.luminosity > ymax) ymax = r.luminosity;
        }
    }

    // margine visivo
    double margin = 0.1 * (ymax - ymin);
    ymin -= margin;
    ymax += margin;

    TCanvas* c3 = new TCanvas("c3", "All spots: luminosity vs overvoltage", 1000, 700);
    // c3->SetLogy();

    TLegend* leg3 = new TLegend(0.72, 0.15, 0.90, 0.88);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);

    bool first = true;
    int colorIndex = 0;

    vector<int> colors = {
        kBlack, kRed+1, kBlue+1, kGreen+2,
        kMagenta+1, kCyan+2, kOrange+7, kViolet,
        kAzure+1, kSpring+5, kPink+7, kTeal+3
    };

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

        TGraphErrors* gr_all = new TGraphErrors(
            n,
            x.data(),
            y.data(),
            ex.data(),
            ey.data()
        );

        int color = colors[colorIndex % colors.size()];
        colorIndex++;

        gr_all->SetTitle("Luminosity vs overvoltage - all spots");
        gr_all->GetXaxis()->SetTitle("overvoltage");
        gr_all->GetYaxis()->SetTitle("luminosity");

        gr_all->SetMarkerStyle(20);
        gr_all->SetMarkerSize(1.0);
        gr_all->SetMarkerColor(color);
        gr_all->SetLineColor(color);
        gr_all->SetLineWidth(2);

        if (first) {
            gr_all->Draw("APL");
            gr_all->GetXaxis()->SetLimits(0, 10);
            gr_all->GetYaxis()->SetRangeUser(ymin, ymax);
            first = false;
        } else {
            gr_all->Draw("PL SAME");
        }

        leg3->AddEntry(gr_all, Form("Spot %d", spot), "lep");
    }

    leg3->Draw();

    c3->SaveAs("canvas_all_spots_luminosity_vs_overvoltage.png");
}