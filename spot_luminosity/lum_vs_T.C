#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

struct SpotData {
    vector<double> T;
    vector<double> l_area_const;
    vector<double> l_area_const_err;
};

void lum_vs_T(const char* filename = "data.csv") {
    gStyle->SetOptStat(0);

    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file " << filename << endl;
        return;
    }

    string line;

    // Salta header
    getline(fin, line);

    map<int, SpotData> data_map;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string token;

        int spot;
        double l_area_const, l_area_const_err, T, v;

        getline(ss, token, ',');
        spot = stoi(token);

        getline(ss, token, ',');
        l_area_const = stod(token);

        getline(ss, token, ',');
        l_area_const_err = stod(token);

        getline(ss, token, ',');
        T = stod(token);

        getline(ss, token, ',');
        v = stod(token); // letto ma non usato

        data_map[spot].T.push_back(T);
        data_map[spot].l_area_const.push_back(l_area_const);
        data_map[spot].l_area_const_err.push_back(l_area_const_err);
    }

    fin.close();

    for (const auto& entry : data_map) {
        int spot = entry.first;
        const SpotData& d = entry.second;

        int n = d.T.size();
        if (n == 0) continue;

        vector<double> ex(n, 0.0);

        TCanvas* c = new TCanvas(Form("c_spot_%d", spot),
                                 Form("Spot %d", spot),
                                 800, 600);

        TGraphErrors* gr = new TGraphErrors(
            n,
            const_cast<double*>(d.T.data()),
            const_cast<double*>(d.l_area_const.data()),
            ex.data(),
            const_cast<double*>(d.l_area_const_err.data())
        );

        gr->SetTitle(Form("Spot %d;T;l_area_const", spot));
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.2);
        gr->SetLineWidth(2);
    

        gr->Draw("AP");

        gr->GetXaxis()->SetTitle("Temperature");
        gr->GetYaxis()->SetTitle("Luminosity");
        gr->GetXaxis()->CenterTitle();
        gr->GetYaxis()->CenterTitle();

        TLegend* leg = new TLegend(0.62, 0.75, 0.88, 0.88);
        leg->AddEntry(gr, Form("Spot %d", spot), "lep");
        leg->SetBorderSize(0);
        leg->Draw();

        c->Update();
    }
}