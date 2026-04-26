//root -l 'lum_vs_T_two_v.C("run=202508_2overvoltages_scanTemp.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TString.h"
#include "TF1.h"
#include "TLatex.h"

using namespace std;

struct Point {
    double x, y;
    double luminosity, error;
    double T, v;
};

map<int, vector<Point>> read_csv(const char* filename) {
    map<int, vector<Point>> data;

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
        string val;

        int spot;
        Point p;

        getline(ss, val, ','); spot = stoi(val);
        getline(ss, val, ','); p.x = stod(val);
        getline(ss, val, ','); p.y = stod(val);
        getline(ss, val, ','); p.luminosity = stod(val);
        getline(ss, val, ','); p.error = stod(val);
        getline(ss, val, ','); p.T = stod(val);
        getline(ss, val, ','); p.v = stod(val);

        data[spot].push_back(p);
    }

    return data;
}

TGraphErrors* make_graph(const vector<Point>& points, double vsel, int color) {
    vector<double> T, L, eT, eL;

    for (const auto& p : points) {
        if (fabs(p.v - vsel) < 1e-6) {
            T.push_back(p.T);
            L.push_back(p.luminosity);
            eT.push_back(0.0);
            eL.push_back(p.error);
        }
    }

    TGraphErrors* gr = new TGraphErrors(
        T.size(),
        T.data(),
        L.data(),
        eT.data(),
        eL.data()
    );

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);

    return gr;
}

void divide_canvas(TCanvas* c, int nspots) {
    int nx = ceil(sqrt(nspots));
    int ny = ceil((double)nspots / nx);
    c->Divide(nx, ny);
}

void get_local_range(const vector<Point>& points, double vsel, double& xmax_local, double& ymax_local) {
    xmax_local = 0;
    ymax_local = 0;

    for (const auto& p : points) {
        if (fabs(p.v - vsel) < 1e-6) {
            if (p.T > xmax_local) xmax_local = p.T;
            if (p.luminosity + p.error > ymax_local)
                ymax_local = p.luminosity + p.error;
        }
    }

    xmax_local *= 1.15;
    ymax_local *= 1.20;
}

bool get_lum_at_T(const vector<Point>& points, double Tsel, double vsel, double& L, double& eL) {
    for (const auto& p : points) {
        if (fabs(p.T - Tsel) < 1e-6 && fabs(p.v - vsel) < 1e-6) {
            L = p.luminosity;
            eL = p.error;
            return true;
        }
    }
    return false;
}

void lum_vs_T_two_v(const char* filename = "data.csv") {

    gStyle->SetOptStat(0);

    auto data = read_csv(filename);
    int nspots = data.size();

    double xmax = 0;
    double ymax = 0;

    for (const auto& [spot, points] : data) {
        for (const auto& p : points) {
            if (p.T > xmax) xmax = p.T;
            if (p.luminosity + p.error > ymax) ymax = p.luminosity + p.error;
        }
    }

    xmax *= 1.15;
    ymax *= 1.20;

    vector<int> colors = {
        kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7,
        kCyan+1, kViolet, kSpring+5, kPink+7, kAzure+7,
        kTeal+3, kGray+2, kBlack, kYellow+2
    };

    // ------------------------------------------------------------
    // Canvas 1: v = 7.0, pad separati
    // ------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c1_v7", "Luminosity vs T - v = 7.0", 1400, 900);
    divide_canvas(c1, nspots);

    int ipad = 1;
    for (const auto& [spot, points] : data) {
        c1->cd(ipad++);

        int color = colors[spot % colors.size()];
        TGraphErrors* gr = make_graph(points, 7.0, color);

       double xmax_local, ymax_local;
      get_local_range(points, 7.0, xmax_local, ymax_local);

        gPad->DrawFrame(0, 0, xmax_local, ymax_local,Form("Spot %d;T;luminosity", spot));

        gr->Draw("PL SAME");

        const Point& p0 = points[0];
        TLegend* leg = new TLegend(0.15, 0.70, 0.88, 0.88);
        leg->AddEntry(gr, Form("x = %.2f, y = %.2f, v = 7.0", p0.x, p0.y), "pl");
        leg->Draw();
    }

    // ------------------------------------------------------------
    // Canvas 2: v = 9.0, pad separati
    // ------------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2_v9", "Luminosity vs T - v = 9.0", 1400, 900);
    divide_canvas(c2, nspots);

    ipad = 1;
    for (const auto& [spot, points] : data) {
        c2->cd(ipad++);

        int color = colors[spot % colors.size()];
        TGraphErrors* gr = make_graph(points, 9.0, color);

        double xmax_local, ymax_local;
        get_local_range(points, 9.0, xmax_local, ymax_local);

        gPad->DrawFrame(0, 0, xmax_local, ymax_local, Form("Spot %d;T;luminosity", spot));

        gr->Draw("PL SAME");

        const Point& p0 = points[0];
        TLegend* leg = new TLegend(0.15, 0.70, 0.88, 0.88);
        leg->AddEntry(gr, Form("x = %.2f, y = %.2f, v = 9.0", p0.x, p0.y), "pl");
        leg->Draw();
    }

    // ------------------------------------------------------------
    // Canvas 3: tutti gli spot insieme, v = 7.0
    // ------------------------------------------------------------
    TCanvas* c3 = new TCanvas("c3_all_v7", "All spots - v = 7.0", 1000, 800);

    gPad->DrawFrame(0, 0, xmax, ymax,
                    "All spots - v = 7.0;T;luminosity");

    TLegend* leg3 = new TLegend(0.75, 0.15, 0.95, 0.90);

    for (const auto& [spot, points] : data) {
        int color = colors[spot % colors.size()];
        TGraphErrors* gr = make_graph(points, 7.0, color);
        gr->Draw("PL SAME");
        leg3->AddEntry(gr, Form("Spot %d", spot), "pl");
    }

    leg3->Draw();

    // ------------------------------------------------------------
    // Canvas 4: tutti gli spot insieme, v = 9.0
    // ------------------------------------------------------------
    TCanvas* c4 = new TCanvas("c4_all_v9", "All spots - v = 9.0", 1000, 800);

    gPad->DrawFrame(0, 0, xmax, ymax,
                    "All spots - v = 9.0;T;luminosity");

    TLegend* leg4 = new TLegend(0.75, 0.15, 0.95, 0.90);

    for (const auto& [spot, points] : data) {
        int color = colors[spot % colors.size()];
        TGraphErrors* gr = make_graph(points, 9.0, color);
        gr->Draw("PL SAME");
        leg4->AddEntry(gr, Form("Spot %d", spot), "pl");
    }

    leg4->Draw();

    // ------------------------------------------------------------
    // Canvas 5: per ogni spot, v = 7.0 e v = 9.0 insieme
    // ------------------------------------------------------------
    TCanvas* c5 = new TCanvas("c5_v7_v9", "Luminosity vs T - v = 7.0 and v = 9.0", 1400, 900);
    divide_canvas(c5, nspots);

    ipad = 1;
    for (const auto& [spot, points] : data) {
        c5->cd(ipad++);

        TGraphErrors* gr7 = make_graph(points, 7.0, kBlue+1);
        TGraphErrors* gr9 = make_graph(points, 9.0, kRed+1);

        double xmax7, ymax7, xmax9, ymax9;
        get_local_range(points, 7.0, xmax7, ymax7);
        get_local_range(points, 9.0, xmax9, ymax9);

        double xmax_local = max(xmax7, xmax9);
        double ymax_local = max(ymax7, ymax9);


        const Point& p0 = points[0];
        gPad->DrawFrame(0, 0, xmax_local, ymax_local, Form("Spot %d: x = %.2f, y = %.2f ; T ;luminosity", spot, p0.x, p0.y));

        gr7->Draw("PL SAME");
        gr9->Draw("PL SAME");

        
        TLegend* leg = new TLegend(0.15, 0.70, 0.38, 0.88);
        leg->AddEntry(gr7, Form("v = 7.0"), "pl");
        leg->AddEntry(gr9, Form("v = 9.0"), "pl");
        leg->Draw();
    }
    // ------------------------------------------------------------
    // Canvas 6: Luminosity variation in percentage, v = 7.0 vs v = 9.0
    // ------------------------------------------------------------
    // ------------------------------------------------------------

    TCanvas* c6 = new TCanvas("c6_ratio_T16_T22", "Ratio L(T=16) / L(T=22) vs spot", 1000, 800);

    vector<double> spot_v7, ratio_v7, espot_v7, eratio_v7;
    vector<double> spot_v9, ratio_v9, espot_v9, eratio_v9;

    for (const auto& [spot, points] : data) {
        if (spot == 8) continue;
        if (spot == 5) continue;
        double L16_7, e16_7, L22_7, e22_7;
        double L16_9, e16_9, L22_9, e22_9;

        if (get_lum_at_T(points, 16.0, 7.0, L16_7, e16_7) &&
            get_lum_at_T(points, 22.0, 7.0, L22_7, e22_7) &&
            L22_7 != 0) {

            double r = L16_7 / L22_7;
            double er = r * sqrt(pow(e16_7 / L16_7, 2) +
                                pow(e22_7 / L22_7, 2));

            spot_v7.push_back(spot);
            ratio_v7.push_back(r);
            espot_v7.push_back(0.0);
            eratio_v7.push_back(er);
        }

        if (get_lum_at_T(points, 16.0, 9.0, L16_9, e16_9) &&
            get_lum_at_T(points, 22.0, 9.0, L22_9, e22_9) &&
            L22_9 != 0) {

            double r = L16_9 / L22_9;
            double er = r * sqrt(pow(e16_9 / L16_9, 2) +
                                pow(e22_9 / L22_9, 2));

            spot_v9.push_back(spot);
            ratio_v9.push_back(r);
            espot_v9.push_back(0.0);
            eratio_v9.push_back(er);
        }
    }

    TGraphErrors* gr_ratio_v7 = new TGraphErrors(
        spot_v7.size(),
        spot_v7.data(),
        ratio_v7.data(),
        espot_v7.data(),
        eratio_v7.data()
    );

    TGraphErrors* gr_ratio_v9 = new TGraphErrors(
        spot_v9.size(),
        spot_v9.data(),
        ratio_v9.data(),
        espot_v9.data(),
        eratio_v9.data()
    );

    gr_ratio_v7->SetMarkerStyle(20);
    gr_ratio_v7->SetMarkerColor(kBlue+1);
    gr_ratio_v7->SetLineColor(kBlue+1);

    gr_ratio_v9->SetMarkerStyle(21);
    gr_ratio_v9->SetMarkerColor(kRed+1);
    gr_ratio_v9->SetLineColor(kRed+1);

    double ymax_ratio = 0;
    for (double r : ratio_v7) if (r > ymax_ratio) ymax_ratio = r;
    for (double r : ratio_v9) if (r > ymax_ratio) ymax_ratio = r;
    ymax_ratio *= 1.4;

    gPad->DrawFrame(0, 0, nspots, ymax_ratio,
                    "Ratio L(T=16) / L(T=22);spot;L_{16}/L_{22}");

    gr_ratio_v7->Draw("P SAME");
    gr_ratio_v9->Draw("P SAME");

    // Fit costante
    TF1* fit_v7 = new TF1("fit_v7", "[0]", 0, nspots);
    TF1* fit_v9 = new TF1("fit_v9", "[0]", 0, nspots);

    fit_v7->SetLineColor(kBlue+1);
    fit_v9->SetLineColor(kRed+1);

    gr_ratio_v7->Fit(fit_v7, "RQ");
    gr_ratio_v9->Fit(fit_v9, "RQ");

    fit_v7->Draw("SAME");
    fit_v9->Draw("SAME");

    double chi2_v7 = fit_v7->GetChisquare();
    int ndf_v7 = fit_v7->GetNDF();
    double fit_result_v7 = fit_v7->GetParameter(0);
    double fit_error_v7 = fit_v7->GetParError(0);


    double chi2_v9 = fit_v9->GetChisquare();
    int ndf_v9 = fit_v9->GetNDF();
    double fit_result_v9 = fit_v9->GetParameter(0);
    double fit_error_v9 = fit_v9->GetParError(0);

    TLegend* leg6 = new TLegend(0.15, 0.70, 0.55, 0.88);
    leg6->AddEntry(gr_ratio_v7, "v = 7.0", "p");
    leg6->AddEntry(gr_ratio_v9, "v = 9.0", "p");
    leg6->AddEntry((TObject*)0, Form("p0 = %.2f #pm %.2f", fit_result_v7, fit_error_v7), "");
    leg6->AddEntry((TObject*)0, Form("p0 = %.2f #pm %.2f", fit_result_v9, fit_error_v9), "");
    leg6->AddEntry(fit_v7, Form("fit const v=7: #chi^{2}/ndf = %.2f/%d", chi2_v7, ndf_v7), "l");
    leg6->AddEntry(fit_v9, Form("fit const v=9: #chi^{2}/ndf = %.2f/%d", chi2_v9, ndf_v9), "l");
    leg6->Draw();


   
}