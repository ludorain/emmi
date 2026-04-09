#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <cstring>

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TKey.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF2.h"
#include "TLegend.h"

using namespace std;


// ------------------------------------------------------------
// Reading coordinates from txt file
// ------------------------------------------------------------
vector<pair<int,int>> read_coordinates(const char* txtfile) {
    vector<pair<int,int>> coords;
    ifstream fin(txtfile);

    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file txt " << txtfile << endl;
        return coords;
    }

    string line;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        double x, y;
        char comma;

        // legge: numero , numero
        if (ss >> x >> comma >> y) {
            coords.push_back({
                static_cast<int>(std::round(x)),
                static_cast<int>(std::round(y))
            });
        } else {
            cerr << "Warning: riga non valida -> " << line << endl;
        }
    }

    fin.close();
    return coords;
}


// ------------------------------------------------------------
// Find TH2F contained in the ROOT file
// ------------------------------------------------------------
TH2F* get_th2f(const char* filename, const char* histname = "") {

    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Errore: impossibile aprire il file " << filename << endl;
        return nullptr;
    }

    TH2F* hist = nullptr;

    if (strlen(histname) > 0) {
        file->GetObject(histname, hist);
        if (!hist) {
            cerr << "Errore: TH2F \"" << histname << "\" non trovato." << endl;
            return nullptr;
        }
    } else {
        TIter next(file->GetListOfKeys());
        TKey* key;

        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom("TH2F")) {
                hist = (TH2F*)obj;
                break;
            }
        }

        if (!hist) {
            cerr << "Errore: nessun TH2F trovato nel file." << endl;
            return nullptr;
        }
    }

    return hist;
}


// ------------------------------------------------------------
// Crea un TH2F ritagliato attorno a (x_center, y_center)
// dimensione: size_x x size_y pixel
// ------------------------------------------------------------
TH2F* make_crop(TH2F* h2, int x_center, int y_center, int size_x=50, int size_y=50, int index=0) {
    if (!h2) return nullptr;

    int half_x = size_x / 2;
    int half_y = size_y / 2;

    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();

    int x_min = x_center - half_x;
    int x_max = x_center + half_x - 1;
    int y_min = y_center - half_y;
    int y_max = y_center + half_y - 1;

    if (x_min < 1) x_min = 1;
    if (y_min < 1) y_min = 1;
    if (x_max > nx) x_max = nx;
    if (y_max > ny) y_max = ny;

    int new_nx = x_max - x_min + 1;
    int new_ny = y_max - y_min + 1;

    TString hname = Form("spot_%d", index);
    TString htitle = Form("Spot %d centered at (%d,%d)", index, x_center, y_center);

    TH2F* hcrop = new TH2F(
        hname,
        htitle,
        new_nx, x_min, x_max + 1,
        new_ny, y_min, y_max + 1
    );

    for (int ix = 0; ix < new_nx; ix++) {
        for (int iy = 0; iy < new_ny; iy++) {
            int old_xbin = x_min + ix;
            int old_ybin = y_min + iy;
            double content = h2->GetBinContent(old_xbin, old_ybin);
            hcrop->SetBinContent(ix + 1, iy + 1, content);
        }
    }

    return hcrop;
}


// ------------------------------------------------------------
// Gaussiana 2D circolare senza fondo
// par[0] = A
// par[1] = x0
// par[2] = y0
// par[3] = sigma
// ------------------------------------------------------------
Double_t gaus2D_circular(Double_t *x, Double_t *par) {
    Double_t X = x[0];
    Double_t Y = x[1];

    Double_t A     = par[0];
    Double_t x0    = par[1];
    Double_t y0    = par[2];
    Double_t sigma = par[3];

    if (sigma <= 0) return 1e30;

    Double_t dx = X - x0;
    Double_t dy = Y - y0;

    return A * exp(-(dx*dx + dy*dy) / (2.0 * sigma * sigma));
}

// ------------------------------------------------------------
// Trova il massimo locale in una finestrella square_size x square_size
// centrata attorno alle coordinate (x_center, y_center)
// square_size deve essere dispari: 3, 5, 7, ...
// ------------------------------------------------------------
double get_local_max_around_point(TH2F* h2, double x_center, double y_center, int square_size) {
    if (!h2) return 0.0;

    if (square_size < 1) square_size = 1;
    if (square_size % 2 == 0) square_size += 1; // forza dispari

    int half_size = square_size / 2;

    int binx_center = h2->GetXaxis()->FindBin(x_center);
    int biny_center = h2->GetYaxis()->FindBin(y_center);

    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();

    int binx_min = binx_center - half_size;
    int binx_max = binx_center + half_size;
    int biny_min = biny_center - half_size;
    int biny_max = biny_center + half_size;

    // protezione bordi
    if (binx_min < 1) binx_min = 1;
    if (biny_min < 1) biny_min = 1;
    if (binx_max > nx) binx_max = nx;
    if (biny_max > ny) biny_max = ny;

    double local_max = -1e99;

    for (int ix = binx_min; ix <= binx_max; ix++) {
        for (int iy = biny_min; iy <= biny_max; iy++) {
            double val = h2->GetBinContent(ix, iy);
            if (val > local_max) local_max = val;
        }
    }

    return local_max;
}


// ------------------------------------------------------------
// Fit multi-step:
// step 1: libero solo A, x0 y0 sigma fissati
// step 2: rilascio sigma
// step 3: fisso A e sigma, liberi x0 y0
// step 4: rilascio tutti
// ------------------------------------------------------------
int run_multistep_fit(TH2F* hspot, TF2* f2, int x_txt, int y_txt, double sigma_fixed = 2.5) {
    if (!hspot || !f2) return -999;

    double xmin = hspot->GetXaxis()->GetXmin();
    double xmax = hspot->GetXaxis()->GetXmax();
    double ymin = hspot->GetYaxis()->GetXmin();
    double ymax = hspot->GetYaxis()->GetXmax();

    int binx_txt = hspot->GetXaxis()->FindBin(x_txt);
    int biny_txt = hspot->GetYaxis()->FindBin(y_txt);

    double A_init = get_local_max_around_point(hspot, x_txt, y_txt, 5);
    if (A_init <= 0) A_init = 1.0;

    // limiti generali
    f2->SetParLimits(0, 0.0, 1e9);
    f2->SetParLimits(1, xmin, xmax);
    f2->SetParLimits(2, ymin, ymax);
    f2->SetParLimits(3, 0.5, 20.0);

    int status1 = 0, status2 = 0, status3 = 0, status4 = 0;

    // ---------------- STEP 1 ----------------
    // libero solo A
    f2->SetParameter(0, A_init);
    f2->FixParameter(1, x_txt);
    f2->FixParameter(2, y_txt);
    f2->FixParameter(3, sigma_fixed);

    f2->ReleaseParameter(0);
    status1 = hspot->Fit(f2, "RQ0");

    // ---------------- STEP 2 ----------------
    // rilascio sigma
    f2->ReleaseParameter(3);
    status2 = hspot->Fit(f2, "RQ0");

    // salvo i valori trovati
    double A_step2     = f2->GetParameter(0);
    double sigma_step2 = f2->GetParameter(3);

    // ---------------- STEP 3 ----------------
    // fisso A e sigma, liberi x0 e y0
    f2->FixParameter(0, A_step2);
    f2->FixParameter(3, sigma_step2);

    f2->ReleaseParameter(1);
    f2->ReleaseParameter(2);
    status3 = hspot->Fit(f2, "RQ0");

    // ---------------- STEP 4 ----------------
    // rilascio tutti
    f2->ReleaseParameter(0);
    f2->ReleaseParameter(1);
    f2->ReleaseParameter(2);
    f2->ReleaseParameter(3);
    status4 = hspot->Fit(f2, "R0");

    cout << "   Step 1 status = " << status1 << endl;
    cout << "   Step 2 status = " << status2 << endl;
    cout << "   Step 3 status = " << status3 << endl;
    cout << "   Step 4 status = " << status4 << endl;

    return status4; // salvo come status finale
}


// ------------------------------------------------------------
// Funzione principale
// ------------------------------------------------------------
void spot_fit_gaussian_step(const char* rootfile, const char* txtfile) {

    TH2F* h2 = get_th2f(rootfile);
    if (!h2) {
        cerr << "Error: no TH2F found in the ROOT file." << endl;
        return;
    }

    cout << "TH2F found: " << h2->GetName() << endl;

    vector<pair<int,int>> coords = read_coordinates(txtfile);

    int n_spots = coords.size();
    cout << "Number of coordinates read (n_spots) = " << n_spots << endl;

    if (n_spots == 0) {
        cerr << "Error: no valid coordinates found in the txt file." << endl;
        return;
    }

    // 4) creation of cropped TH2F for each coordinate
    vector<TH2F*> spots;
    spots.reserve(n_spots);

    for (int i = 0; i < n_spots; i++) {
        int x = coords[i].first;
        int y = coords[i].second;

        TH2F* hspot = make_crop(h2, x, y, 50, 50, i);
        if (hspot) {
            spots.push_back(hspot);
            cout << "Created spot_" << i << " centered at (" << x << ", " << y << ")" << endl;
        }
    }

 // Draw crops normal
    int ncols = std::ceil(std::sqrt(n_spots));
    int nrows = std::ceil((double)n_spots / ncols);

    TCanvas* c2 = new TCanvas("c2", "Spots", 1400, 900);
    c2->Divide(ncols, nrows);
    for (int i = 0; i < (int)spots.size(); i++) {
        c2->cd(i + 1);
        spots[i]->Draw("COLZ");
    }
    c2->Update();
//Draw crops LEGO
    TCanvas* c3 = new TCanvas("c3", "Spots LEGO drawing option", 1400, 900);
    c3->Divide(ncols, nrows);
    for (int i = 0; i < (int)spots.size(); i++) {
        c3->cd(i + 1);
        spots[i]->Draw("LEGO");
    }
    c3->Update();


    // --------------------------------------------------------
    // 5) Fitting the spots with a 2D Gaussian function
    // --------------------------------------------------------
    vector<TF2*> fit_functions;
    fit_functions.reserve(spots.size());

    vector<int> fit_statuses;
    fit_statuses.reserve(spots.size());

    cout << "\n================ FIT RESULTS ================\n";

    // apro il csv una sola volta
    ofstream fout("fit_results.csv");
    if (fout.is_open()) {
        fout << "spot_index,x_txt,y_txt,A,A_err,x0,x0_err,y0,y0_err,sigma,sigma_err,fit_status\n";
    }

    for (int i = 0; i < (int)spots.size(); i++) {
        TH2F* hspot = spots[i];
        if (!hspot) continue;

        double xmin = hspot->GetXaxis()->GetXmin();
        double xmax = hspot->GetXaxis()->GetXmax();
        double ymin = hspot->GetYaxis()->GetXmin();
        double ymax = hspot->GetYaxis()->GetXmax();

        TString fname = Form("f2_spot_%d", i);
        TF2* f2 = new TF2(fname, gaus2D_circular, xmin, xmax, ymin, ymax, 4);
        f2->SetParNames("A", "x0", "y0", "sigma");
        int fitStatus = run_multistep_fit(hspot, f2, coords[i].first, coords[i].second, 2.5);

        fit_functions.push_back(f2);
        fit_statuses.push_back(fitStatus);

        cout << "\nSpot " << i
             << "  [coord approx input = (" << coords[i].first
             << ", " << coords[i].second << ")]" << endl;
        cout << "Final fit status = " << fitStatus << endl;
        cout << "A      = " << f2->GetParameter(0) << " +/- " << f2->GetParError(0) << endl;
        cout << "x0     = " << f2->GetParameter(1) << " +/- " << f2->GetParError(1) << endl;
        cout << "y0     = " << f2->GetParameter(2) << " +/- " << f2->GetParError(2) << endl;
        cout << "sigma  = " << f2->GetParameter(3) << " +/- " << f2->GetParError(3) << endl;

        if (fout.is_open()) {
            fout << i << ","
                 << coords[i].first << ","
                 << coords[i].second << ","
                 << f2->GetParameter(0) << ","
                 << f2->GetParError(0) << ","
                 << f2->GetParameter(1) << ","
                 << f2->GetParError(1) << ","
                 << f2->GetParameter(2) << ","
                 << f2->GetParError(2) << ","
                 << f2->GetParameter(3) << ","
                 << f2->GetParError(3) << ","
                 << fitStatus << "\n";
        }
    }

    if (fout.is_open()) fout.close();


    // --------------------------------------------------------
    // Draw spots with fit contours
    // --------------------------------------------------------
    TCanvas* c4 = new TCanvas("c4", "Spots with fit", 1400, 900);
    c4->Divide(ncols, nrows);

    for (int i = 0; i < (int)spots.size(); i++) {
        c4->cd(i + 1);
        spots[i]->Draw("COLZ");
        fit_functions[i]->Draw("SAME CONT3");
    }
    c4->Update();


    // --------------------------------------------------------
    // Plot fit parameters
    // --------------------------------------------------------
    TCanvas* c5 = new TCanvas("c5", "Plot of Fit Results", 1400, 900);
    c5->Divide(3, 2);

    const int n_parameters = 4;
    const char* par_names[n_parameters]  = {"A", "x0", "y0", "sigma"};
    const char* par_titles[n_parameters] = {
        "Fit parameter A",
        "Fit parameter x0",
        "Fit parameter y0",
        "Fit parameter sigma"
    };

    vector<TH1F*> parameters_hist;
    parameters_hist.reserve(n_parameters);

    for (int p = 0; p < n_parameters; p++) {
        TString hname  = Form("h_%s", par_names[p]);
        TString htitle = Form("%s;Spot index;%s", par_titles[p], par_names[p]);

        TH1F* hpar = new TH1F(hname, htitle, n_spots, 0.5, n_spots + 0.5);

        for (int i = 0; i < n_spots; i++) {
            TString label = Form("%d", i);
            hpar->GetXaxis()->SetBinLabel(i + 1, label);
        }

        parameters_hist.push_back(hpar);
    }

    for (int i = 0; i < (int)fit_functions.size(); i++) {
        if (fit_statuses[i] != 0) continue;

        for (int p = 0; p < n_parameters; p++) {
            double value = fit_functions[i]->GetParameter(p);
            double error = fit_functions[i]->GetParError(p);

            parameters_hist[p]->SetBinContent(i + 1, value);
            parameters_hist[p]->SetBinError(i + 1, error);
        }
    }

    TH1F* h_x_txt = new TH1F("h_x_txt", "x0 parameter;Spot index;x", n_spots, 0.5, n_spots + 0.5);
    TH1F* h_y_txt = new TH1F("h_y_txt", "y0 parameter;Spot index;y", n_spots, 0.5, n_spots + 0.5);

    for (int i = 0; i < n_spots; i++) {
        h_x_txt->SetBinContent(i + 1, coords[i].first);
        h_y_txt->SetBinContent(i + 1, coords[i].second);

        TString label = Form("%d", i);
        h_x_txt->GetXaxis()->SetBinLabel(i + 1, label);
        h_y_txt->GetXaxis()->SetBinLabel(i + 1, label);
    }

    h_x_txt->SetMarkerStyle(24);
    h_y_txt->SetMarkerStyle(24);
    h_x_txt->SetMarkerSize(0.9);
    h_y_txt->SetMarkerSize(0.9);
    h_x_txt->SetMarkerColor(kGreen+2);
    h_y_txt->SetMarkerColor(kGreen+2);

    for (int p = 0; p < n_parameters; p++) {
        c5->cd(p + 1);
        parameters_hist[p]->SetMarkerStyle(20);
        parameters_hist[p]->SetMarkerSize(0.9);
        parameters_hist[p]->Draw("E1");

        if (p == 1) h_x_txt->Draw("P SAME");
        if (p == 2) h_y_txt->Draw("P SAME");
    }

    TH1F* h_dx = new TH1F("h_dx", "x0(fit) - x(txt);Spot index;#Delta x", n_spots, 0.5, n_spots + 0.5);
    TH1F* h_dy = new TH1F("h_dy", "y0(fit) - y(txt);Spot index;#Delta y", n_spots, 0.5, n_spots + 0.5);

    for (int i = 0; i < (int)fit_functions.size(); i++) {
        if (fit_statuses[i] != 0) continue;

        h_dx->SetBinContent(i + 1, fit_functions[i]->GetParameter(1) - coords[i].first);
        h_dx->SetBinError(i + 1, fit_functions[i]->GetParError(1));

        h_dy->SetBinContent(i + 1, fit_functions[i]->GetParameter(2) - coords[i].second);
        h_dy->SetBinError(i + 1, fit_functions[i]->GetParError(2));
    }

    c5->cd(5);
    h_dx->SetMarkerStyle(20);
    h_dx->SetMarkerSize(0.9);
    h_dx->Draw("E1");

    c5->cd(6);
    h_dy->SetMarkerStyle(20);
    h_dy->SetMarkerSize(0.9);
    h_dy->Draw("E1");

    c5->Update();
}