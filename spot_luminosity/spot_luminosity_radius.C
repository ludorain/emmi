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
#include "TGraphErrors.h"
#include "TMultiGraph.h"

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
TFitResultPtr run_multistep_fit(TH2F* hspot, TF2* f2, int x_txt, int y_txt, double sigma_fixed = 2.5) {
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
    TFitResultPtr r = hspot->Fit(f2, "RS0");
    status4 = r->Status();


    cout << "   Step 1 status = " << status1 << endl;
    cout << "   Step 2 status = " << status2 << endl;
    cout << "   Step 3 status = " << status3 << endl;
    cout << "   Step 4 status = " << status4 << endl;

    return r; // salvo come status finale
}

// ------------------------------------------------------------
// Calcola la somma dei contenuti dei bin entro un raggio R
// dal centro (x0, y0) ricavato dal fit.
// ------------------------------------------------------------
double luminosity_sum(TH2F* hspot, TF2* ffit, double radius) {
    if (!hspot || !ffit) return 0.0;

    // Recupera il centro dal fit (par 1 = x0, par 2 = y0)
    double x0 = ffit->GetParameter(1);
    double y0 = ffit->GetParameter(2);

    double sum_counts = 0.0;

    // Ottimizzazione: cicliamo solo nel rettangolo che circoscrive il cerchio
    int binX_min = hspot->GetXaxis()->FindBin(x0 - radius);
    int binX_max = hspot->GetXaxis()->FindBin(x0 + radius);
    int binY_min = hspot->GetYaxis()->FindBin(y0 - radius);
    int binY_max = hspot->GetYaxis()->FindBin(y0 + radius);

    // Protezione contro i bordi dell'istogramma
    if (binX_min < 1) binX_min = 1;
    if (binY_min < 1) binY_min = 1;
    if (binX_max > hspot->GetNbinsX()) binX_max = hspot->GetNbinsX();
    if (binY_max > hspot->GetNbinsY()) binY_max = hspot->GetNbinsY();

    for (int ix = binX_min; ix <= binX_max; ix++) {
        for (int iy = binY_min; iy <= binY_max; iy++) {
            // Coordinate del centro del bin corrente
            double bin_x = hspot->GetXaxis()->GetBinCenter(ix);
            double bin_y = hspot->GetYaxis()->GetBinCenter(iy);

            // Calcolo distanza euclidea dal centro del fit
            double dx = bin_x - x0;
            double dy = bin_y - y0;
            double dist = std::sqrt(dx*dx + dy*dy);

            if (dist <= radius) {
                sum_counts += hspot->GetBinContent(ix, iy);
            }
        }
    }

    return sum_counts;
}

// ------------------------------------------------------------
// Calcola la luminosità integrando la funzione di fit bin per bin
// entro un raggio R. Questo metodo rimuove il rumore dei singoli bin
// e considera la forma analitica della Gaussiana.
// ------------------------------------------------------------
double luminosity_analytical_bins(TH2F* hspot, TF2* ffit, double radius) {
    if (!hspot || !ffit) return 0.0;

    double x0 = ffit->GetParameter(1);
    double y0 = ffit->GetParameter(2);

    double total_integral = 0.0;

    // Definiamo i limiti dei bin da scansionare
    int binX_min = hspot->GetXaxis()->FindBin(x0 - radius);
    int binX_max = hspot->GetXaxis()->FindBin(x0 + radius);
    int binY_min = hspot->GetYaxis()->FindBin(y0 - radius);
    int binY_max = hspot->GetYaxis()->FindBin(y0 + radius);

    // Protezione bordi
    binX_min = std::max(1, binX_min);
    binY_min = std::max(1, binY_min);
    binX_max = std::min(hspot->GetNbinsX(), binX_max);
    binY_max = std::min(hspot->GetNbinsY(), binY_max);

    for (int ix = binX_min; ix <= binX_max; ix++) {
        for (int iy = binY_min; iy <= binY_max; iy++) {
            
            // Coordinate del centro e dei bordi del bin
            double x_center = hspot->GetXaxis()->GetBinCenter(ix);
            double y_center = hspot->GetYaxis()->GetBinCenter(iy);
            
            double dx = x_center - x0;
            double dy = y_center - y0;
            double dist = std::sqrt(dx*dx + dy*dy);

            // Se il centro del bin è nel raggio, integriamo la funzione sull'area del bin
            if (dist <= radius) {
                double x_low = hspot->GetXaxis()->GetBinLowEdge(ix);
                double x_up  = hspot->GetXaxis()->GetBinUpEdge(ix);
                double y_low = hspot->GetYaxis()->GetBinLowEdge(iy);
                double y_up  = hspot->GetYaxis()->GetBinUpEdge(iy);

            
                total_integral += ffit->Integral(x_low, x_up, y_low, y_up);
            }
        }
    }
    
    return total_integral;
}

double calculate_analytical_bins_error(TH2F* hspot, TF2* f2, TFitResultPtr r, double radius) {
    double x0 = f2->GetParameter(1);
    double y0 = f2->GetParameter(2);
    double error_sum_2 = 0.0;

    for (int ix = 1; ix <= hspot->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hspot->GetNbinsY(); iy++) {
            double dx = hspot->GetXaxis()->GetBinCenter(ix) - x0;
            double dy = hspot->GetYaxis()->GetBinCenter(iy) - y0;
            
            if (sqrt(dx*dx + dy*dy) <= radius) {
                double xlow = hspot->GetXaxis()->GetBinLowEdge(ix);
                double xup  = hspot->GetXaxis()->GetBinUpEdge(ix);
                double ylow = hspot->GetYaxis()->GetBinLowEdge(iy);
                double yup  = hspot->GetYaxis()->GetBinUpEdge(iy);

                // Definiamo i limiti del bin per l'integrazione
                double x_range[2] = {xlow, xup};
                double y_range[2] = {ylow, yup};

                // Calcola l'errore dell'integrale in questo bin usando la matrice di covarianza
                double bin_err = f2->IntegralError(2, x_range, y_range, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
                error_sum_2 += pow(bin_err, 2);
            }
        }
    }
    return sqrt(error_sum_2);
}

// ------------------------------------------------------------
// Calcola la luminosità come integrale della funzione gaussiana fittata
// ------------------------------------------------------------

double luminosity_analytical_area(TH2F* hspot, TF2* ffit){
    
    double volume = 0.0;
    double A = ffit->GetParameter(0);
    double sigma = ffit->GetParameter(3);

    volume = 2.0 * TMath::Pi() * A *  sigma * sigma;

    return volume;
}

double calculate_analytical_area_error(TF2* f2, TFitResultPtr r) {
    double A = f2->GetParameter(0);
    double sigma = f2->GetParameter(3);
    
    double errA = f2->GetParError(0);
    double errSigma = f2->GetParError(3);
    
    // Recuperiamo la covarianza tra il parametro 0 (A) e il parametro 3 (sigma)
    double covASigma = r->CovMatrix(0, 3);

    double dLdA = 2.0 * M_PI * sigma * sigma;
    double dLdSigma = 4.0 * M_PI * A * sigma;

    double varL = pow(dLdA * errA, 2) + pow(dLdSigma * errSigma, 2) + 2.0 * dLdA * dLdSigma * covASigma;
    
    return sqrt(varL);
}
// ------------------------------------------------------------
// Funzione principale
// ------------------------------------------------------------
void spot_luminosity_radius(const char* rootfile, const char* txtfile) {

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
    //c2->SaveAs("spots.png");

    // --------------------------------------------------------
    // 5) Fitting the spots with a 2D Gaussian function
    // --------------------------------------------------------
    vector<TF2*> fit_functions;
    fit_functions.reserve(spots.size());

    vector<int> fit_statuses;
    fit_statuses.reserve(spots.size());

    cout << "\n================ FIT RESULTS ================\n";

    // apro il csv una sola volta
    /*
    ofstream fout("fit_results.csv");
    if (fout.is_open()) {
        fout << "spot_index,x_txt,y_txt,A,A_err,x0,x0_err,y0,y0_err,sigma,sigma_err,fit_status\n";
    }*/

    std::ofstream csvFile("luminosity_results.csv");
    if (!csvFile.is_open()) {
        std::cerr << "Errore: impossibile creare il file CSV." << std::endl;
        return;
    }

    csvFile << "spot,l_area_const,l_area_const_err\n";


    vector<int> good_spots = {0, 1};

    for (int i : good_spots) {
        TH2F* hspot = spots[i];
        if (!hspot) continue;

        double xmin = hspot->GetXaxis()->GetXmin();
        double xmax = hspot->GetXaxis()->GetXmax();
        double ymin = hspot->GetYaxis()->GetXmin();
        double ymax = hspot->GetYaxis()->GetXmax();

        TString fname = Form("f2_spot_%d", i);
        TF2* f2 = new TF2(fname, gaus2D_circular, xmin, xmax, ymin, ymax, 4);
        f2->SetParNames("A", "x0", "y0", "sigma");
        TFitResultPtr fitResult = run_multistep_fit(hspot, f2, coords[i].first, coords[i].second, 2.5);
        int fitStatus = fitResult->Status();

        fit_functions.push_back(f2);
        fit_statuses.push_back(fitStatus);


        // errore sull'integrale totale della gaussiana: costante per tutti i raggi
        double l_area_const     = luminosity_analytical_area(hspot, f2);
        double l_area_const_err = calculate_analytical_area_error(f2, fitResult);

        csvFile << i << "," << l_area_const << "," << l_area_const_err << "\n";

        /*
        //Valutazione tipologie di luminosità in funzione del raggio 

       vector<double> radii = {10.0, 10.5, 11.0, 11.5, 12.0, 13.0, 14.0, 15.0, 17.5, 20.0};

        TGraphErrors *g_bin  = new TGraphErrors();
        TGraphErrors *g_area = new TGraphErrors();

        

        for (int ir = 0; ir < (int)radii.size(); ir++) {
            double r = radii[ir];

            double l_bin     = luminosity_analytical_bins(hspot, f2, r);
            double l_bin_err = calculate_analytical_bins_error(hspot, f2, fitResult, r);

            g_bin->SetPoint(ir, r, l_bin);
            g_bin->SetPointError(ir, 0.0, l_bin_err);

            g_area->SetPoint(ir, r, l_area_const);
            g_area->SetPointError(ir, 0.0, l_area_const_err);
        }

        // Formattazione grafica
        g_bin->SetMarkerStyle(21);
        g_bin->SetMarkerColor(kBlue);
        g_bin->SetLineColor(kBlue);
         

        g_area->SetMarkerStyle(20);
        g_area->SetMarkerColor(kRed);
        g_area->SetLineColor(kRed);
        g_area->SetLineStyle(2);
        g_area->SetLineWidth(2);

        // Creazione Canvas
        TCanvas *c_comp = new TCanvas(Form("c_comp_%d", i),
                                    Form("Luminosity Comparison Spot %d", i),
                                    800, 600);

        TMultiGraph *mg = new TMultiGraph();
        mg->Add(g_bin, "LP");
        mg->Add(g_area, "LP");

        mg->SetTitle(Form("Growth Curve Spot %d;Radius [pixels];Luminosity [counts]", i));
        mg->Draw("A");

        TLegend *leg = new TLegend(0.55, 0.20, 0.88, 0.40);
        leg->AddEntry(g_bin,  "Analytical Bins Integral", "lp");
        leg->AddEntry(g_area, "Total Gaussian Volume (2#pi A#sigma^{2})", "lp");
        leg->Draw();

        c_comp->Update();
        // c_comp->SaveAs(Form("growth_curve_%d.png", i));
        */

    //salvataggio dei dati nel file csv
        /*
        cout << "\nSpot " << i
             << "  [coord approx input = (" << coords[i].first
             << ", " << coords[i].second << ")]" << endl;
        cout << "x0     = " << f2->GetParameter(1) << " +/- " << f2->GetParError(1) << endl;
        cout << "y0     = " << f2->GetParameter(2) << " +/- " << f2->GetParError(2) << endl;

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
                 << fitStatus << ","
                 << lum_sum << "\n";
        }*/
    }
    csvFile.close();
   // if (fout.is_open()) fout.close();
        


    

}