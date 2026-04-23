#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TKey.h"
#include "TEllipse.h"

using namespace std;

// 1. Lettura coordinate dal file .txt
vector<pair<double,double>> read_coordinates(const char* txtfile) {
    vector<pair<double,double>> coords;
    ifstream fin(txtfile);
    if (!fin.is_open()) return coords;
    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        double x, y;
        char comma;
        // Supporta sia "x,y" che "x y"
        if (ss >> x >> comma >> y || (ss.clear(), ss >> x >> y)) {
            coords.push_back({x, y});
        }
    }
    fin.close();
    return coords;
}

// Funzione per estrarre l'istogramma dal file ROOT
TH2F* get_th2f(const char* filename) {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) return nullptr;
    TIter next(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        if (TClass::GetClass(key->GetClassName())->InheritsFrom("TH2F")) 
            return (TH2F*)key->ReadObj();
    }
    return nullptr;
}

// 2 & 3. Calcolo Somma ed Errore
pair<double, double> calculate_luminosity(TH2F* h2, double x0, double y0, double radius) {
    double sum_counts = 0.0;
    double sum_error_sq = 0.0;

    // Definiamo un range di bin attorno al punto per non scansionare tutto l'istogramma
    int binX_min = h2->GetXaxis()->FindBin(x0 - radius);
    int binX_max = h2->GetXaxis()->FindBin(x0 + radius);
    int binY_min = h2->GetYaxis()->FindBin(y0 - radius);
    int binY_max = h2->GetYaxis()->FindBin(y0 + radius);

    for (int ix = binX_min; ix <= binX_max; ix++) {
        for (int iy = binY_min; iy <= binY_max; iy++) {
            if (ix < 1 || ix > h2->GetNbinsX() || iy < 1 || iy > h2->GetNbinsY()) continue;

            double dx = h2->GetXaxis()->GetBinCenter(ix) - x0;
            double dy = h2->GetYaxis()->GetBinCenter(iy) - y0;
            
            if (sqrt(dx*dx + dy*dy) <= radius) {
                sum_counts += h2->GetBinContent(ix, iy);
                double err = h2->GetBinError(ix, iy);
                sum_error_sq += err * err;
            }
        }
    }
    return {sum_counts, sqrt(sum_error_sq)};
}

// Funzione Principale
void spot_luminosity_sum(const char* rootfile, const char* txtfile, double integration_radius = 10.0) {
    TH2F* h2 = get_th2f(rootfile);
    if (!h2) {
        cerr << "Errore: TH2F non trovato." << endl;
        return;
    }

    vector<pair<double,double>> coords = read_coordinates(txtfile);
    if (coords.empty()) {
        cerr << "Errore: nessuna coordinata letta dal file txt." << endl;
        return;
    }

    // 4. Apertura file CSV per il salvataggio
    ofstream csvFile("luminosity_results.csv");
    csvFile << "x,y,luminosity,error" << endl;

    cout << "Analisi in corso su " << coords.size() << " spot..." << endl;

    for (size_t i = 0; i < coords.size(); i++) {
        double x_coord = coords[i].first;
        double y_coord = coords[i].second;

        // Calcolo diretto usando le coordinate del file
        auto results = calculate_luminosity(h2, x_coord, y_coord, integration_radius);

        // Scrittura su CSV
        csvFile << x_coord << "," << y_coord << "," 
                << results.first << "," << results.second << endl;

        cout << "Spot " << i << " [" << x_coord << ", " << y_coord << "] -> Lum: " 
             << results.first << " +/- " << results.second << endl;
    }

    csvFile.close();
    cout << "Completato. Dati salvati in 'luminosity_results.csv'" << endl;
}