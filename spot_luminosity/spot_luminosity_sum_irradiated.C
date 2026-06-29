#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TKey.h"
#include "TEllipse.h"
#include "TClass.h"

using namespace std;

// ------------------------------------------------------------
// Structure for one source
// txt file format:
// x, y, area, radius
// ------------------------------------------------------------
struct SourceInfo {
    double x;
    double y;
    double area;
    double radius;
};

// ------------------------------------------------------------
// Read coordinates and radius from txt file
// Expected format:
// x, y, area, radius
// ------------------------------------------------------------
vector<SourceInfo> read_coordinates(const char* txtfile) {
    vector<SourceInfo> sources;

    ifstream fin(txtfile);

    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file " << txtfile << endl;
        return sources;
    }

    string line;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        // Replace commas with spaces, so both formats work:
        // x, y, area, radius
        // x y area radius
        replace(line.begin(), line.end(), ',', ' ');

        stringstream ss(line);

        SourceInfo src;

        if (ss >> src.x >> src.y >> src.area >> src.radius) {
            sources.push_back(src);
        } else {
            cerr << "Warning: riga non valida ignorata: " << line << endl;
        }
    }

    fin.close();

    return sources;
}

// ------------------------------------------------------------
// Function to extract the TH2F histogram from the ROOT file
// ------------------------------------------------------------
TH2F* get_th2f(const char* filename) {
    TFile* file = TFile::Open(filename, "READ");

    if (!file || file->IsZombie()) {
        cerr << "Errore: impossibile aprire il file ROOT " << filename << endl;
        return nullptr;
    }

    TIter next(file->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        if (TClass::GetClass(key->GetClassName())->InheritsFrom("TH2F")) {
            TH2F* h = (TH2F*)key->ReadObj();
            h->SetDirectory(0);
            file->Close();
            return h;
        }
    }

    file->Close();

    return nullptr;
}

// ------------------------------------------------------------
// Calculate luminosity and error inside a circular region 
// Error as sum of bin errors showed to be overestimated
// Error as sqrt of sum of bin errors squared showed to be underestimated
// ------------------------------------------------------------
/*
pair<double, double> calculate_luminosity(
    TH2F* h2,
    double x0,
    double y0,
    double radius
) {
    double sum_counts = 0.0;
    //double sum_error = 0.0;
    double sum_error_sq = 0.0;

    int binX_min = h2->GetXaxis()->FindBin(x0 - radius);
    int binX_max = h2->GetXaxis()->FindBin(x0 + radius);
    int binY_min = h2->GetYaxis()->FindBin(y0 - radius);
    int binY_max = h2->GetYaxis()->FindBin(y0 + radius);

    for (int ix = binX_min; ix <= binX_max; ix++) {
        for (int iy = binY_min; iy <= binY_max; iy++) {

            if (ix < 1 || ix > h2->GetNbinsX()) continue;
            if (iy < 1 || iy > h2->GetNbinsY()) continue;

            double dx = h2->GetXaxis()->GetBinCenter(ix) - x0;
            double dy = h2->GetYaxis()->GetBinCenter(iy) - y0;

            if (sqrt(dx * dx + dy * dy) <= radius) {
                sum_counts += h2->GetBinContent(ix, iy);

                double err = h2->GetBinError(ix, iy);
                //sum_error += err;
                sum_error_sq += err * err;
            }
        }
    }
    //return {sum_counts, sum_error};
    return {sum_counts, sqrt(sum_error_sq)};
} */

pair<double, double> calculate_luminosity(
    TH2F* h2,
    double x0,
    double y0,
    double radius,
    double rho = 0.1   // correlation coefficient: 0 = quadrature, 1 = linear sum
) {
    double sum_counts = 0.0;

    double sum_error = 0.0;      // sum_i sigma_i
    double sum_error_sq = 0.0;   // sum_i sigma_i^2

    int binX_min = h2->GetXaxis()->FindBin(x0 - radius);
    int binX_max = h2->GetXaxis()->FindBin(x0 + radius);
    int binY_min = h2->GetYaxis()->FindBin(y0 - radius);
    int binY_max = h2->GetYaxis()->FindBin(y0 + radius);

    for (int ix = binX_min; ix <= binX_max; ix++) {
        for (int iy = binY_min; iy <= binY_max; iy++) {

            if (ix < 1 || ix > h2->GetNbinsX()) continue;
            if (iy < 1 || iy > h2->GetNbinsY()) continue;

            double dx = h2->GetXaxis()->GetBinCenter(ix) - x0;
            double dy = h2->GetYaxis()->GetBinCenter(iy) - y0;

            if (sqrt(dx * dx + dy * dy) <= radius) {
                sum_counts += h2->GetBinContent(ix, iy);

                double err = h2->GetBinError(ix, iy);

                sum_error += err;
                sum_error_sq += err * err;
            }
        }
    }

    // Partial-correlation uncertainty:
    // rho = 0 -> independent errors, sum in quadrature, showed to be underestimated
    // rho = 1 -> fully correlated errors, linear sum, showed to be overestimated
    double total_variance =
        (1.0 - rho) * sum_error_sq
        + rho * sum_error * sum_error;

    double total_error = sqrt(total_variance);

    return {sum_counts, total_error};
}

// ------------------------------------------------------------
// Main function
// ------------------------------------------------------------
void spot_luminosity_sum_irradiated(const char* rootfile, const char* txtfile) {

    TH2F* h2 = get_th2f(rootfile);

    if (!h2) {
        cerr << "Errore: TH2F non trovato." << endl;
        return;
    }

    vector<SourceInfo> sources = read_coordinates(txtfile);

    if (sources.empty()) {
        cerr << "Errore: nessuna coordinata letta dal file txt." << endl;
        return;
    }

    ofstream csvFile("luminosity_results.csv");

    csvFile << "x,y,luminosity,error" << endl;

    cout << "Analisi in corso su " << sources.size() << " spot..." << endl;

    for (size_t i = 0; i < sources.size(); i++) {

        double x_coord = sources[i].x;
        double y_coord = sources[i].y;
        double radius = sources[i].radius;

        if (radius <= 0) {
            cerr << "Warning: radius non valido per lo spot " << i
                 << ". Spot ignorato." << endl;
            continue;
        }

        auto results = calculate_luminosity(
            h2,
            x_coord,
            y_coord,
            radius
        );

        csvFile << x_coord << ","
                << y_coord << ","
                << results.first << ","
                << results.second << endl;
    }

    csvFile.close();

    cout << "Completato. Dati salvati in 'luminosity_results.csv'" << endl;
}