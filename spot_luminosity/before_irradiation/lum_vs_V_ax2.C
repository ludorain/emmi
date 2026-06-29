//root -l 'lum_vs_V_ax2.C("run=202508_scanV.csv")'

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

void lum_vs_V_ax2(const char* filename = "data.csv") {

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
    vector<double> A_vals, A_errs, chi2_vals;

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

        gr->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f °C", spot, rows[0].x, rows[0].y, rows[0].T));
        gr->GetXaxis()->SetTitle("Overvoltage (V)");
        gr->GetYaxis()->SetTitle("Luminosity");

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);

        gr->Draw("AP");

        TF1* f_quad = new TF1("f_quad", "[0]*x*x", 0, 10);
        f_quad->SetParameter(0, 1.0);

        gr->Fit(f_quad, "Q");  // Q = quiet

        double A     = f_quad->GetParameter(0);
        double A_err = f_quad->GetParError(0);
        double chi2  = f_quad->GetChisquare();

        // salva risultati
        spot_ids.push_back(spot);
        A_vals.push_back(A);
        A_errs.push_back(A_err);
        chi2_vals.push_back(chi2);
        /*
        TLegend* leg = new TLegend(0.15, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        //leg->AddEntry(gr, Form("x = %.2f, y = %.2f", rows[0].x, rows[0].y), "lep");
        //leg->AddEntry((TObject*)0, Form("T = %.1f", rows[0].T), "");
        leg->Draw(); */
    }

    //c1->SaveAs("canvas_luminosity_vs_overvoltage.png");

    
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
        gr_all->GetXaxis()->SetTitle("Overvoltage (V)");
        gr_all->GetYaxis()->SetTitle("Luminosity");

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
        TF1* f_quad = new TF1(Form("f_quad_%d", spot), "[0]*x*x", 0, 10);
        f_quad->SetParameter(0, 1.0);
        gr_all->Fit(f_quad, "Q SAME");
        
        leg3->AddEntry(gr_all, Form("Spot %d", spot), "lep");
    }

    leg3->Draw();

    //c3->SaveAs("canvas_all_spots_luminosity_vs_overvoltage.png");

    // ============================================================
    // CANVAS 4: A vs spot ID
    // ============================================================
    int n = spot_ids.size();

    vector<double> x(n), ex(n, 0.0);

    for (int i = 0; i < n; i++) {
        x[i] = spot_ids[i];
    }

    TCanvas* c4 = new TCanvas("c4", "Fit parameter A vs spot", 800, 600);

    TGraphErrors* grA = new TGraphErrors(n, x.data(), A_vals.data(), ex.data(), A_errs.data());

    grA->SetTitle("Parametro A del fit quadratico;Spot;A");
    grA->SetMarkerStyle(20);
    grA->Draw("AP");

    TLegend* leg4 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->AddEntry(grA, "Fit: A x^{2}", "lep");
    leg4->Draw();

    //c4->SaveAs("canvas_A_vs_spot.png");



    // ============================================================
    // CANVAS 5: chi square vs spot ID
    // ============================================================

    TCanvas* c5 = new TCanvas("c5", "Chi2 vs spot", 800, 600);

    TGraph* grChi2 = new TGraph(n, x.data(), chi2_vals.data());

    grChi2->SetTitle("#chi^{2} del fit;Spot;#chi^{2}");
    grChi2->SetMarkerStyle(21);
    grChi2->Draw("AP");

    TLegend* leg5 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg5->SetBorderSize(0);
    leg5->SetFillStyle(0);
    leg5->AddEntry(grChi2, "Fit quadratico", "p");
    leg5->Draw();

    //c5->SaveAs("canvas_chi2_vs_spot.png");
}