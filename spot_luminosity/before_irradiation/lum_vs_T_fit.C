// root -l 'lum_vs_T_fit.C("run=202508_2overvoltages_scanTemp.csv")'

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
#include "TStyle.h"
#include "TAxis.h"
#include "TString.h"
#include "TF1.h"
#include "TLatex.h"
#include "TMath.h"

using namespace std;

struct Point
{
    double x, y;
    double luminosity, error;
    double T, v;
};

map<int, vector<Point>> read_csv(const char *filename, bool kelvin = false)
{
    map<int, vector<Point>> data;
    ifstream fin(filename);
    if (!fin.is_open())
    {
        cerr << "Errore: impossibile aprire " << filename << endl;
        return data;
    }

    string line;
    getline(fin, line); // skip header

    while (getline(fin, line))
    {
        if (line.empty())
            continue;
        stringstream ss(line);
        string val;
        int spot;
        Point p;

        try
        {
            getline(ss, val, ',');
            spot = stoi(val);
            getline(ss, val, ',');
            p.x = stod(val);
            getline(ss, val, ',');
            p.y = stod(val);
            getline(ss, val, ',');
            p.luminosity = stod(val);
            getline(ss, val, ',');
            p.error = stod(val);
            if (kelvin)
            {
                getline(ss, val, ',');
                p.T = stod(val) + 273.15;
            }
            else
            {
                getline(ss, val, ',');
                p.T = stod(val);
            }
            getline(ss, val, ',');
            p.v = stod(val);
            data[spot].push_back(p);
        }
        catch (...)
        {
            continue;
        }
    }
    return data;
}

TGraphErrors *make_graph(const vector<Point> &points, double vsel, int color)
{
    vector<double> T, L, eT, eL;
    for (const auto &p : points)
    {
        if (fabs(p.v - vsel) < 1e-6)
        {
            T.push_back(p.T);
            L.push_back(p.luminosity);
            eT.push_back(0.0);
            eL.push_back(p.error);
        }
    }
    if (T.empty())
        return nullptr;
    TGraphErrors *gr = new TGraphErrors(T.size(), T.data(), L.data(), eT.data(), eL.data());
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);
    return gr;
}

void divide_canvas(TCanvas *c, int nspots)
{
    int nx = ceil(sqrt(nspots));
    int ny = ceil((double)nspots / nx);
    c->Divide(nx, ny);
}

void get_local_range(const vector<Point> &points, double vsel, double &xmax_local, double &ymax_local)
{
    xmax_local = 0;
    ymax_local = 0;
    for (const auto &p : points)
    {
        if (fabs(p.v - vsel) < 1e-6)
        {
            if (p.T > xmax_local)
                xmax_local = p.T;
            if (p.luminosity + p.error > ymax_local)
                ymax_local = p.luminosity + p.error;
        }
    }
}

void lum_vs_T_fit(const char *filename = "run=202508_2overvoltages_scanTemp.csv", bool kelvin = false, double x_initial = 15., double x_final = 25.)
{
    double power_exp = 0.06;
    if (kelvin)
    {
        x_initial += 273.15;
        x_final += 273.15;
    }
    double T0 = x_initial;
    gStyle->SetOptFit(0);
    auto data = read_csv(filename, kelvin);
    if (data.empty())
        return;

    int nspots = data.size();

    // Vettori per Canvas 9
    vector<double> spot_ids_dbl, ex_dbl;
    vector<double> B_7_vals, B_7_errs, B_9_vals, B_9_errs;

    TCanvas *c5 = new TCanvas("c5", "Luminosity vs T - Multiple Spots", 1400, 900);
    divide_canvas(c5, nspots);

    int ipad = 1;
    for (auto const &[spot, points] : data)
    {
        c5->cd(ipad++);

        TGraphErrors *gr7 = make_graph(points, 7.0, kBlue + 1);
        TGraphErrors *gr9 = make_graph(points, 9.0, kRed + 1);

        double xmax7, ymax7, xmax9, ymax9;
        get_local_range(points, 7.0, xmax7, ymax7);
        get_local_range(points, 9.0, xmax9, ymax9);
        double ymax_local = max(ymax7, ymax9) * 1.2;

        gPad->DrawFrame(x_initial, 0, x_final, ymax_local, Form("Spot %d;T;Luminosity", spot));

        // Fit v=7
        if (gr7)
        {
            TF1 *f7 = new TF1(Form("f7_%d", spot),
                              "[0]*(TMath::Exp([1]*(x)))",
                              x_initial, x_final);
            f7->SetParameters(0.001, power_exp);
            f7->SetLineColor(kBlue + 1);
            gr7->Fit(f7, "RQ");
            gr7->Draw("P SAME");

            B_7_vals.push_back(f7->GetParameter(1));
            B_7_errs.push_back(f7->GetParError(1));
        }

        // Fit v=9
        if (gr9)
        {
            //TF1 *f9 = new TF1(Form("f9_%d", spot),
            //                  "[0]*(TMath::Exp([1]*(x - [2])) - 1)",
            //                  x_initial, x_final);
            //f9->SetParameters(0.01, power_exp, T0);

            TF1 *f9 = new TF1(Form("f9_%d", spot),
                              "[0]*(TMath::Exp([1]*(x)))",
                              x_initial, x_final);
            f9->SetParameters(0.001, power_exp);
            f9->SetLineColor(kRed + 1);
            gr9->Fit(f9, "RQ");
            gr9->Draw("P SAME");

            B_9_vals.push_back(f9->GetParameter(1));
            B_9_errs.push_back(f9->GetParError(1));
        }

        spot_ids_dbl.push_back((double)spot);
        ex_dbl.push_back(0.0);

        TLegend *leg5 = new TLegend(0.6, 0.7, 0.88, 0.88);
        leg5->AddEntry(gr7, "V_{over} = 7.0", "lep");
        leg5->AddEntry(gr9, "V_{over} = 9.0", "lep");
        leg5->AddEntry(B_7_vals.empty() ? nullptr : new TGraphErrors(1, &spot_ids_dbl.back(), &B_7_vals.back(), &ex_dbl.back(), &B_7_errs.back()), Form("Fit V_{over}=7: B = %.5f #pm %.5f", B_7_vals.back(), B_7_errs.back()), "l");
        leg5->AddEntry(B_9_vals.empty() ? nullptr : new TGraphErrors(1, &spot_ids_dbl.back(), &B_9_vals.back(), &ex_dbl.back(), &B_9_errs.back()), Form("Fit V_{over}=9: B = %.5f #pm %.5f", B_9_vals.back(), B_9_errs.back()), "l");

        leg5->Draw();
    }
    // ------------------------------------------------------------
    // --- Canvas 6: Lum vs Temperature for a single spot ---
    // ------------------------------------------------------------

    TCanvas *c6 = new TCanvas("c6", "Luminosity vs Temperature for Spot 3", 800, 600);
    // Prendi solo lo spot 3
    auto it = data.find(3);
    if (it == data.end()) {
            std::cout << "Spot 3 non trovato!" << std::endl;
            return;
        }

    const auto& points = it->second;

    // Grafici per v = 7 e v = 9
    TGraphErrors* gr7 = make_graph(points, 7.0, kBlue+1);
    TGraphErrors* gr9 = make_graph(points, 9.0, kRed+1);

    // Range automatico
    double xmax7, ymax7, xmax9, ymax9;
    get_local_range(points, 7.0, xmax7, ymax7);
    get_local_range(points, 9.0, xmax9, ymax9);

    double ymax = std::max(ymax7, ymax9) * 1.2;

     // Prendi coordinate spot
    const Point& p0 = points[0];

    // Frame unico
    gPad->DrawFrame(15, 0, 23, ymax, Form("Spot 3: x = %.2f, y = %.2f ; T ; luminosity", p0.x, p0.y));


    // Fit v=7
        if (gr7)
        {
            TF1 *f7 = new TF1(Form("f7_%d", 3),
                              "[0]*(TMath::Exp([1]*(x)))",
                              x_initial, x_final);
            f7->SetParameters(0.001, power_exp);
            f7->SetLineColor(kBlue + 1);
            gr7->Fit(f7, "RQ");
            gr7->Draw("P SAME");
        }

        // Fit v=9
        if (gr9)
        {
            TF1 *f9 = new TF1(Form("f9_%d", 3),
                              "[0]*(TMath::Exp([1]*(x)))",
                              x_initial, x_final);
            f9->SetParameters(0.001, power_exp);
            f9->SetLineColor(kRed + 1);
            gr9->Fit(f9, "RQ");
            gr9->Draw("P SAME");
        }
    // Disegno grafici
    //gr7->Draw("P SAME");
    //gr9->Draw("P SAME");

    // Legenda
    TLegend* leg51 = new TLegend(0.15, 0.70, 0.38, 0.88);
    leg51->AddEntry(gr7, "V_over = 7.0V", "pl");
    leg51->AddEntry(gr9, "V_over = 9.0V", "pl");
    leg51->Draw();


    // --- Canvas 9: B vs Spot ---
    TCanvas *c9 = new TCanvas("c9", "B parameter vs Spot ID", 800, 600);
    c9->SetGrid();

    TGraphErrors *grB7 = new TGraphErrors(spot_ids_dbl.size(), spot_ids_dbl.data(), B_7_vals.data(), ex_dbl.data(), B_7_errs.data());
    TGraphErrors *grB9 = new TGraphErrors(spot_ids_dbl.size(), spot_ids_dbl.data(), B_9_vals.data(), ex_dbl.data(), B_9_errs.data());

    grB7->SetMarkerStyle(21);
    grB7->SetMarkerColor(kBlue + 1);
    grB7->SetLineColor(kBlue + 1);
    grB9->SetMarkerStyle(22);
    grB9->SetMarkerColor(kRed + 1);
    grB9->SetLineColor(kRed + 1);

    grB7->SetTitle("Exponent B vs Spot ID;Spot ID;B value");
    grB7->GetXaxis()->SetLimits(-1, nspots + 1);
    grB7->Draw("AP");
    grB9->Draw("P SAME");

    // Fit lineari (costanti)
    TF1 *pol0_7 = new TF1("pol0_7", "pol0", -1, nspots + 1);
    pol0_7->SetLineColor(kBlue + 1);
    grB7->Fit(pol0_7, "RQ");

    TF1 *pol0_9 = new TF1("pol0_9", "pol0", -1, nspots + 1);
    pol0_9->SetLineColor(kRed + 1);
    grB9->Fit(pol0_9, "RQ");

    TLegend *leg9 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg9->AddEntry(grB7, "V_{over} = 7.0", "lep");
    leg9->AddEntry(grB9, "V_{over} = 9.0", "lep");
    leg9->AddEntry(pol0_7, Form("Fit V_{over}=7: B = %.5f #pm %.5f", pol0_7->GetParameter(0), pol0_7->GetParError(0)), "l");
    leg9->AddEntry(pol0_9, Form("Fit V_{over}=9: B = %.5f #pm %.5f", pol0_9->GetParameter(0), pol0_9->GetParError(0)), "l");
    leg9->Draw();
}