// ------------------------------------------------------------
// How to run:
// root -l -q 'spot_fit_gaussian_circular.C("run=<root_file>.root","coordinates.txt")'
// ------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TH2F.h"
#include "TKey.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TString.h"
#include <TF2.h>



using namespace std;


// ------------------------------------------------------------
// Reading coordinates from txt file
// da valutare se tenere in double
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

    // Apri file
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
        return nullptr;
    }

    TH2F* hist = nullptr;

    // Caso 1: nome fornito
    if (strlen(histname) > 0) {
        file->GetObject(histname, hist);
        if (!hist) {
            std::cerr << "Errore: TH2F \"" << histname << "\" non trovato." << std::endl;
            return nullptr;
        }
    }
    
    // Caso 2: cerca automaticamente
    else {
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
            std::cerr << "Errore: nessun TH2F trovato nel file." << std::endl;
            return nullptr;
        }
    }

    return hist;
}

// ------------------------------------------------------------
// Crea un TH2F ritagliato attorno a (x_center, y_center)
// dimensione: size_x x size_y pixel
// ------------------------------------------------------------

//prima prova fatta con 50x50, aree più grandi sono peggio
//necessario implementare un modo che trovi l'area corretta, magari dal programma python e salvare nel .txt

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

    // Protezione bordi
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


// Gaussiana 2D circolare + fondo costante
// par[0] = A
// par[1] = x0
// par[2] = y0
// par[3] = sigma
// par[4] = B
// ------------------------------------------------------------
Double_t gaus2D_circular(Double_t *x, Double_t *par) {
    Double_t X = x[0];
    Double_t Y = x[1];

    Double_t A     = par[0];
    Double_t x0    = par[1];
    Double_t y0    = par[2];
    Double_t sigma = par[3];
    Double_t B     = par[4];

    if (sigma <= 0) return 1e30;

    Double_t dx = X - x0;
    Double_t dy = Y - y0;

    return A * exp(-(dx*dx + dy*dy) / (2.0 * sigma * sigma)) + B;
}



// ------------------------------------------------------------
// Funzione principale
// ------------------------------------------------------------
void spot_fit_gaussian_circular(const char* rootfile, const char* txtfile) {
    

    // PRENDO  il th2f dal file
    TH2F* h2 = get_th2f(rootfile);
    if (!h2) {
        cerr << "Error: no TH2F found in the ROOT file." << endl;
        return;
    }

    TCanvas* c1 = new TCanvas("c1", "TH2F display", 800, 600);
    h2->Draw("COLZ");
    c1->Update();

    cout << "TH2F found: " << h2->GetName() << endl;

    // 2) open txt file and read coordinates
    vector<pair<int,int>> coords = read_coordinates(txtfile);

    // 3) number of coordinates
    int n_spots = coords.size();
    cout << "Number of coordinates read (n_spots) = " << n_spots << endl;

    if (n_spots == 0) {
        cerr << "Error: no valid coordinates found in the txt file." << endl;
     //   h2->Close();
        return;
    }

    // 4) creation of cropped TH2F for each coordinate
    vector<TH2F*> spots;
    spots.reserve(n_spots);

    for (int i = 0; i < n_spots; i++) {
        int x = coords[i].first;
        int y = coords[i].second;

        
        TH2F* hspot = make_crop(h2, x, y, 50, 50, i);//again qui da valutare se 50x50 va bene, oppure se è meglio un'altra dimensione o addirittura un modo per trovare la dimensione corretta
        if (hspot) {
            spots.push_back(hspot);
            cout << "Created spot_" << i << " centered at (" << x << ", " << y << ")" << endl;
        }
    }

    // Draw crops in a Canva
    int ncols = std::ceil(std::sqrt(n_spots));
    int nrows = std::ceil((double)n_spots / ncols);

    TCanvas* c2 = new TCanvas("c2", "Spots", 1400, 900);
    c2->Divide(ncols, nrows);

    for (int i = 0; i < (int)spots.size(); i++) {
        c2->cd(i + 1);
        spots[i]->Draw("COLZ");
    }

    c2->Update();

    TCanvas* c3 = new TCanvas("c3", "Spots LEGO drawing option", 1400, 900);
    c3->Divide(ncols, nrows);

    for (int i = 0; i < (int)spots.size(); i++) {
            c3->cd(i + 1);
            spots[i]->Draw("LEGO");
        }
    c3->Update();

    /*
    // Saving crops in a new file root:
    TFile* fout = TFile::Open("spots_output.root", "RECREATE");
    for (auto h : spots) {
        h->Write();
    }
    fout->Close();

    cout << "Salvati i TH2F ritagliati in spots_output.root" << endl;

    */

    //5) Fitting the spots with a 2D Gaussian function
    
    vector<TF2*> fit_functions;
    fit_functions.reserve(spots.size());


    cout << "\n================ FIT RESULTS ================\n";

    for (int i = 0; i < (int)spots.size(); i++) {
        TH2F* hspot = spots[i];
        if (!hspot) continue;

        // range del crop
        double xmin = hspot->GetXaxis()->GetXmin();
        double xmax = hspot->GetXaxis()->GetXmax();
        double ymin = hspot->GetYaxis()->GetXmin();
        double ymax = hspot->GetYaxis()->GetXmax();

        // valori iniziali
        int maxbinx, maxbiny, maxbinz;
        hspot->GetMaximumBin(maxbinx, maxbiny, maxbinz);

        double x0_init = hspot->GetXaxis()->GetBinCenter(coords[i].first);
        double y0_init = hspot->GetYaxis()->GetBinCenter(coords[i].second);
        double maxval  = hspot->GetMaximum();
        double minval  = hspot->GetMinimum();
        double sigma_init = 3.0;   // valore iniziale ragionevole per spot piccoli

        TString fname = Form("f2_spot_%d", i);
        TF2* f2 = new TF2(fname, gaus2D_circular, xmin, xmax, ymin, ymax, 5);
        f2->SetParNames("A", "x0", "y0", "sigma", "B");

        f2->SetParameters(maxval - minval, x0_init, y0_init, sigma_init, minval);

        // limiti utili per stabilizzare il fit
        f2->SetParLimits(0, 0.0, 1e9);         // A > 0
        f2->SetParLimits(1, xmin, xmax);       // x0 nel crop
        f2->SetParLimits(2, ymin, ymax);       // y0 nel crop
        f2->SetParLimits(3, 0.5, 20.0);        // sigma positiva
        f2->SetParLimits(4, -1e9, 1e9);        // fondo libero

        // fit
        int fitStatus = hspot->Fit(f2, "RQ");  
        // R = usa range funzione
        // Q = meno output ROOT a schermo

        fit_functions.push_back(f2);
 

        // stampa risultati
        cout << "\nSpot " << i
             << "  [coord approx input = (" << coords[i].first
             << ", " << coords[i].second << ")]" << endl;
        cout << "Fit status = " << fitStatus << endl;
        cout << "A      = " << f2->GetParameter(0) << " +/- " << f2->GetParError(0) << endl;
        cout << "x0     = " << f2->GetParameter(1) << " +/- " << f2->GetParError(1) << endl;
        cout << "y0     = " << f2->GetParameter(2) << " +/- " << f2->GetParError(2) << endl;
        cout << "sigma  = " << f2->GetParameter(3) << " +/- " << f2->GetParError(3) << endl;
        cout << "B      = " << f2->GetParameter(4) << " +/- " << f2->GetParError(4) << endl;
    }

    TCanvas* c4 = new TCanvas("c4", "Spots with fit", 1400, 900);
    c4->Divide(ncols, nrows);

    for (int i = 0; i < (int)spots.size(); i++) {
        c4->cd(i + 1);
        spots[i]->Draw("COLZ");
        //fit_functions[i]->SetNpx(100);
        //fit_functions[i]->SetNpy(100);
        //fit_functions[i]->SetLineColor(kRed);
        fit_functions[i]->Draw("SAME CONT3");

    c4->Update(); }


}




