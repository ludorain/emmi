#include <TFile.h>
#include <TH2.h>
#include <TH2D.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TAxis.h>
#include <TMath.h>
#include <TROOT.h>
#include <TKey.h>
#include <TSystem.h>
#include <TMarker.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

struct Spot {
    double x;
    double y;
    double halfSizeX;
    double halfSizeY;
};

struct FitOutput {
    int index;
    bool success;

    double x_init, y_init;

    double mux, emux;
    double muy, emuy;

    double sigmax, esigmax;
    double sigmay, esigmay;

    double fwhmx, efwhmx;
    double fwhmy, efwhmy;

    double amplitude, eamplitude;
    double background, ebackground;

    double chi2;
    int ndf;

    double xwin_min, xwin_max;
    double ywin_min, ywin_max;
};

const double FWHM_FACTOR = 2.0 * sqrt(2.0 * log(2.0));

vector<Spot> load_spots(const char* filename) {
    vector<Spot> spots;
    ifstream fin(filename);

    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file " << filename << endl;
        return spots;
    }

    string line;
    int lineNumber = 0;

    while (getline(fin, line)) {
        ++lineNumber;

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        istringstream iss(line);
        Spot s;
        s.halfSizeX = -1.0;
        s.halfSizeY = -1.0;

        if (!(iss >> s.x >> s.y)) {
            cerr << "Riga ignorata (" << lineNumber << "): " << line << endl;
            continue;
        }

        if (!(iss >> s.halfSizeX >> s.halfSizeY)) {
            s.halfSizeX = -1.0;
            s.halfSizeY = -1.0;
        }

        spots.push_back(s);
    }

    return spots;
}

TH2* get_histogram(TFile* file, const char* histName = "") {
    if (!file || file->IsZombie()) return nullptr;

    TH2* h2 = nullptr;

    if (histName && strlen(histName) > 0) {
        file->GetObject(histName, h2);
        return h2;
    }

    TIter next(file->GetListOfKeys());
    TObject* keyObj = nullptr;

    while ((keyObj = next())) {
        TKey* key = (TKey*)keyObj;
        TObject* obj = key->ReadObj();
        if (obj && obj->InheritsFrom("TH2")) {
            h2 = (TH2*)obj;
            return h2;
        }
    }

    return nullptr;
}

TH2D* make_local_histogram(TH2* h2, const Spot& spot,
                           double defaultHalfSizeX,
                           double defaultHalfSizeY,
                           double& xmin, double& xmax,
                           double& ymin, double& ymax) {
    TAxis* xax = h2->GetXaxis();
    TAxis* yax = h2->GetYaxis();

    double hx = (spot.halfSizeX > 0 ? spot.halfSizeX : defaultHalfSizeX);
    double hy = (spot.halfSizeY > 0 ? spot.halfSizeY : defaultHalfSizeY);

    xmin = spot.x - hx;
    xmax = spot.x + hx;
    ymin = spot.y - hy;
    ymax = spot.y + hy;

    int binx1 = xax->FindBin(xmin);
    int binx2 = xax->FindBin(xmax);
    int biny1 = yax->FindBin(ymin);
    int biny2 = yax->FindBin(ymax);

    binx1 = max(1, binx1);
    binx2 = min(xax->GetNbins(), binx2);
    biny1 = max(1, biny1);
    biny2 = min(yax->GetNbins(), biny2);

    xmin = xax->GetBinLowEdge(binx1);
    xmax = xax->GetBinUpEdge(binx2);
    ymin = yax->GetBinLowEdge(biny1);
    ymax = yax->GetBinUpEdge(biny2);

    int nx = binx2 - binx1 + 1;
    int ny = biny2 - biny1 + 1;

    TH2D* hlocal = new TH2D(Form("hlocal_tmp_%lld", gSystem->Now()),
                            "Local fit window;X;Y",
                            nx, xmin, xmax,
                            ny, ymin, ymax);

    for (int ix = binx1; ix <= binx2; ++ix) {
        for (int iy = biny1; iy <= biny2; ++iy) {
            double x = xax->GetBinCenter(ix);
            double y = yax->GetBinCenter(iy);
            double z = h2->GetBinContent(ix, iy);

            int localBin = hlocal->FindBin(x, y);
            hlocal->SetBinContent(localBin, z);
        }
    }

    return hlocal;
}

FitOutput fit_one_spot_local(TH2* h2, const Spot& spot, int index,
                             double defaultHalfSizeX = 10.0,
                             double defaultHalfSizeY = 10.0,
                             bool draw = false) {
    FitOutput out;
    out.index = index;
    out.success = false;

    out.x_init = spot.x;
    out.y_init = spot.y;

    out.mux = out.emux = 0.0;
    out.muy = out.emuy = 0.0;
    out.sigmax = out.esigmax = 0.0;
    out.sigmay = out.esigmay = 0.0;
    out.fwhmx = out.efwhmx = 0.0;
    out.fwhmy = out.efwhmy = 0.0;
    out.amplitude = out.eamplitude = 0.0;
    out.background = out.ebackground = 0.0;
    out.chi2 = 0.0;
    out.ndf = 0;
    out.xwin_min = out.xwin_max = 0.0;
    out.ywin_min = out.ywin_max = 0.0;

    double xmin, xmax, ymin, ymax;
    TH2D* hlocal = make_local_histogram(h2, spot,
                                        defaultHalfSizeX, defaultHalfSizeY,
                                        xmin, xmax, ymin, ymax);

    out.xwin_min = xmin;
    out.xwin_max = xmax;
    out.ywin_min = ymin;
    out.ywin_max = ymax;

    int nx = hlocal->GetNbinsX();
    int ny = hlocal->GetNbinsY();

    double zmax = -1e99;
    double zmin =  1e99;
    int maxBinX = 1;
    int maxBinY = 1;

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double z = hlocal->GetBinContent(ix, iy);
            if (z > zmax) {
                zmax = z;
                maxBinX = ix;
                maxBinY = iy;
            }
            if (z < zmin) zmin = z;
        }
    }

    double x0 = hlocal->GetXaxis()->GetBinCenter(maxBinX);
    double y0 = hlocal->GetYaxis()->GetBinCenter(maxBinY);

    double hx = 0.5 * (xmax - xmin);
    double hy = 0.5 * (ymax - ymin);

    double initSigmaX = max(1e-3, hx / 3.0);
    double initSigmaY = max(1e-3, hy / 3.0);

    TF2* f2 = new TF2(Form("f2_gaus_%d", index),
        "[0] + [1]*exp(-0.5*((x-[2])/[3])*((x-[2])/[3]) -0.5*((y-[4])/[5])*((y-[4])/[5]))",
        xmin, xmax, ymin, ymax);

    f2->SetParNames("Bkg", "Amp", "MuX", "SigmaX", "MuY", "SigmaY");
    f2->SetParameters(zmin, zmax - zmin, x0, initSigmaX, y0, initSigmaY);

    f2->SetParLimits(1, 0.0, max(1.0, 10.0 * fabs(zmax - zmin) + 1.0));
    f2->SetParLimits(2, xmin, xmax);
    f2->SetParLimits(3, 1e-3, xmax - xmin);
    f2->SetParLimits(4, ymin, ymax);
    f2->SetParLimits(5, 1e-3, ymax - ymin);

    TFitResultPtr r = hlocal->Fit(f2, "SQNR");

    if ((int)r != 0) {
        if (draw) {
            TCanvas* cfail = new TCanvas(Form("c_fit_fail_%d", index),
                                         Form("Fit fail punto %d", index),
                                         900, 700);
            hlocal->SetTitle(Form("Punto %d - fit non riuscito;X;Y", index));
            hlocal->Draw("COLZ");

            TMarker* minit = new TMarker(spot.x, spot.y, 29);
            minit->SetMarkerSize(2.0);
            minit->SetMarkerColor(kBlue);
            minit->Draw("SAME");

            cfail->Update();
        }

        delete f2;
        return out;
    }

    out.success      = true;
    out.background   = f2->GetParameter(0);
    out.ebackground  = f2->GetParError(0);
    out.amplitude    = f2->GetParameter(1);
    out.eamplitude   = f2->GetParError(1);

    out.mux          = f2->GetParameter(2);
    out.emux         = f2->GetParError(2);
    out.sigmax       = fabs(f2->GetParameter(3));
    out.esigmax      = f2->GetParError(3);

    out.muy          = f2->GetParameter(4);
    out.emuy         = f2->GetParError(4);
    out.sigmay       = fabs(f2->GetParameter(5));
    out.esigmay      = f2->GetParError(5);

    out.fwhmx        = FWHM_FACTOR * out.sigmax;
    out.efwhmx       = FWHM_FACTOR * out.esigmax;
    out.fwhmy        = FWHM_FACTOR * out.sigmay;
    out.efwhmy       = FWHM_FACTOR * out.esigmay;

    out.chi2         = r->Chi2();
    out.ndf          = r->Ndf();

    if (draw) {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);

        TCanvas* c = new TCanvas(Form("c_fit_%d", index),
                                 Form("Fit punto %d", index), 950, 750);

        hlocal->SetTitle(Form("Punto %d;X;Y", index));
        hlocal->Draw("COLZ");

        f2->SetNpx(100);
        f2->SetNpy(100);
        f2->SetLineColor(kRed);
        f2->SetLineWidth(2);
        f2->Draw("CONT3 SAME");

        TMarker* minit = new TMarker(spot.x, spot.y, 29);
        minit->SetMarkerSize(2.0);
        minit->SetMarkerColor(kBlue);
        minit->Draw("SAME");

        TMarker* mfit = new TMarker(out.mux, out.muy, 34);
        mfit->SetMarkerSize(2.0);
        mfit->SetMarkerColor(kGreen + 2);
        mfit->Draw("SAME");

        c->Update();
    }

    return out;
}

void print_results(const vector<FitOutput>& results) {
    cout << fixed << setprecision(4);

    for (size_t i = 0; i < results.size(); ++i) {
        const FitOutput& res = results[i];

        cout << "============================================================" << endl;
        cout << "Punto " << res.index
             << "   coord iniziale = (" << res.x_init << ", " << res.y_init << ")" << endl;

        if (!res.success) {
            cout << "Fit non riuscito." << endl;
            continue;
        }

        cout << "MuX     = " << res.mux    << " +/- " << res.emux    << endl;
        cout << "MuY     = " << res.muy    << " +/- " << res.emuy    << endl;
        cout << "SigmaX  = " << res.sigmax << " +/- " << res.esigmax << endl;
        cout << "SigmaY  = " << res.sigmay << " +/- " << res.esigmay << endl;
        cout << "FWHMX   = " << res.fwhmx  << " +/- " << res.efwhmx  << endl;
        cout << "FWHMY   = " << res.fwhmy  << " +/- " << res.efwhmy  << endl;
        cout << "Amp     = " << res.amplitude  << " +/- " << res.eamplitude  << endl;
        cout << "Bkg     = " << res.background << " +/- " << res.ebackground << endl;
        cout << "Chi2/NDF= " << res.chi2 << " / " << res.ndf;
        if (res.ndf > 0) cout << " = " << res.chi2 / res.ndf;
        cout << endl;
    }

    cout << "============================================================" << endl;
}

void save_results_csv(const vector<FitOutput>& results, const char* csvFile) {
    ofstream fout(csvFile);
    if (!fout.is_open()) {
        cerr << "Errore: impossibile scrivere il file CSV " << csvFile << endl;
        return;
    }

    fout << "index,success,x_init,y_init,"
         << "xwin_min,xwin_max,ywin_min,ywin_max,"
         << "mux,emux,muy,emuy,"
         << "sigmax,esigmax,sigmay,esigmay,"
         << "fwhmx,efwhmx,fwhmy,efwhmy,"
         << "amplitude,eamplitude,background,ebackground,"
         << "chi2,ndf,chi2_ndf\n";

    fout << fixed << setprecision(8);

    for (const auto& res : results) {
        double chi2ndf = (res.ndf > 0) ? res.chi2 / res.ndf : -1.0;

        fout << res.index << ","
             << res.success << ","
             << res.x_init << ","
             << res.y_init << ","
             << res.xwin_min << ","
             << res.xwin_max << ","
             << res.ywin_min << ","
             << res.ywin_max << ","
             << res.mux << ","
             << res.emux << ","
             << res.muy << ","
             << res.emuy << ","
             << res.sigmax << ","
             << res.esigmax << ","
             << res.sigmay << ","
             << res.esigmay << ","
             << res.fwhmx << ","
             << res.efwhmx << ","
             << res.fwhmy << ","
             << res.efwhmy << ","
             << res.amplitude << ","
             << res.eamplitude << ","
             << res.background << ","
             << res.ebackground << ","
             << res.chi2 << ","
             << res.ndf << ","
             << chi2ndf
             << "\n";
    }

    fout.close();
    cout << "Risultati salvati in: " << csvFile << endl;
}

void spot_fit_gaussian_ellisse(const char* rootFile,
                       const char* txtFile,
                       const char* histName = "",
                       int showIndex = -1,
                       double defaultHalfSizeX = 10.0,
                       double defaultHalfSizeY = 10.0,
                       const char* csvFile = "fit_results.csv") {
    TFile* file = TFile::Open(rootFile, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Errore: impossibile aprire il file ROOT " << rootFile << endl;
        return;
    }

    TH2* h2 = get_histogram(file, histName);
    if (!h2) {
        cerr << "Errore: nessun TH2 trovato nel file." << endl;
        file->Close();
        return;
    }

    vector<Spot> spots = load_spots(txtFile);
    if (spots.empty()) {
        cerr << "Errore: nessuna coordinata valida trovata nel file txt." << endl;
        file->Close();
        return;
    }

    if (showIndex != -1 && (showIndex < 1 || showIndex > (int)spots.size())) {
        cerr << "Errore: showIndex fuori range. Deve essere tra 1 e " << spots.size() << endl;
        file->Close();
        return;
    }

    cout << endl;
    cout << "Istogramma usato: " << h2->GetName() << endl;
    cout << "Numero punti da fittare: " << spots.size() << endl;
    cout << "File CSV output: " << csvFile << endl;
    cout << endl;

    vector<FitOutput> results;
    results.reserve(spots.size());

    for (size_t i = 0; i < spots.size(); ++i) {
        bool drawThis = (showIndex == (int)i + 1);

        FitOutput res = fit_one_spot_local(h2, spots[i], i + 1,
                                           defaultHalfSizeX, defaultHalfSizeY,
                                           drawThis);
        results.push_back(res);
    }

    print_results(results);
    save_results_csv(results, csvFile);

    file->Close();
}