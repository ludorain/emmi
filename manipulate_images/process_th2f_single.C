// ============================================================
// Come eseguire:
// .L process_th2f_single.C
// process_th2f("prova.root","contrast_stretching,remove_bias")
// ============================================================

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <cmath>

#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TStyle.h"
#include "TROOT.h"

using namespace std;

// ============================================================
// Utility: split string by comma
// Esempio: "contrast_stretching,remove_bias"
// ============================================================
vector<string> split_operations(const string& ops_string) {
    vector<string> ops;
    stringstream ss(ops_string);
    string item;

    while (getline(ss, item, ',')) {
        // rimuove spazi iniziali/finali
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);

        if (!item.empty())
            ops.push_back(item);
    }
    return ops;
}

// ============================================================
// Cerca un TH2F nel file:
// - se histname è dato, prende quello
// - altrimenti prende il primo TH2F trovato
// ============================================================
TH2F* load_th2f(const char* filename, const char* histname = "") {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Errore: impossibile aprire il file " << filename << endl;
        return nullptr;
    }

    TH2F* hist = nullptr;

    if (strlen(histname) > 0) {
        file->GetObject(histname, hist);
        if (!hist) {
            cerr << "Errore: TH2F \"" << histname << "\" non trovato nel file." << endl;
            file->Close();
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
            file->Close();
            return nullptr;
        }
    }

    // Clono l'istogramma per renderlo indipendente dal file
    TH2F* hclone = (TH2F*)hist->Clone("h_original");
    hclone->SetDirectory(0);

    file->Close();
    return hclone;
}

// ============================================================
// Funzione 1: contrast_stretching
// percentile 1% - 99%, riscalato in [0,1]
// ============================================================
TH2F* contrast_stretching(TH2F* h, const char* newname = "h_contrast") {
    if (!h) {
        std::cerr << "Errore: istogramma nullo." << std::endl;
        return nullptr;
    }
    cout << " --- contrast_stretching" << endl;

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    vector<double> values;
    values.reserve(nx * ny);

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            values.push_back(h->GetBinContent(ix, iy));
        }
    }

    sort(values.begin(), values.end());

    int i1  = (int)(0.01 * (values.size() - 1));
    int i99 = (int)(0.99 * (values.size() - 1));

    double pmin = values[i1];
    double pmax = values[i99];

    TH2F* hout = (TH2F*)h->Clone(newname);
    hout->SetDirectory(0);
    hout->Reset();

    if (pmax <= pmin) {
        cerr << "Errore in contrast_stretching: pmax <= pmin" << endl;
        return hout;
    }

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double val = h->GetBinContent(ix, iy);
            double out = 0.0;

            if (val <= pmin) out = 0.0;
            else if (val >= pmax) out = 1.0;
            else out = (val - pmin) / (pmax - pmin);

            hout->SetBinContent(ix, iy, out);
        }
    }

    return hout;
}

// ============================================================
// Funzione 2: remove_column_bias_median
// sottrae la mediana di ogni colonna
// ============================================================
TH2F* remove__column_bias_median(TH2F* h, const char* newname = "h_remove_bias_median") {
std::cout << " --- remove columns bias" << std::endl;

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    // Clona istogramma per output
    TH2F* h_out = (TH2F*)h2->Clone("h_no_bias");

    for (int ix = 1; ix <= nx; ix++) {

        std::vector<double> column_values;

        // Prendi tutti i valori della colonna ix
        for (int iy = 1; iy <= ny; iy++) {
            column_values.push_back(h2->GetBinContent(ix, iy));
        }

        // Calcolo mediana
        std::sort(column_values.begin(), column_values.end());
        double median;

        int n = column_values.size();
        if (n % 2 == 0)
            median = 0.5 * (column_values[n/2 - 1] + column_values[n/2]);
        else
            median = column_values[n/2];

        // Sottrai la mediana alla colonna
        for (int iy = 1; iy <= ny; iy++) {
            double val = h2->GetBinContent(ix, iy);
            h_out->SetBinContent(ix, iy, val - median);
        }
    }

    return h_out;
}
// ============================================================
// Funzione 3: remove_column_bias_fit
// sottrae il valore di fit polinomiale di ogni colonna, facendo 
// un fit bkg+gaussiana per ogni colonna
// ============================================================



// ============================================================
// Funzione 4: remove_hot_pixels
// sostituisce i pixel sopra mean + nsigma*rms con la media locale 3x3
// ============================================================
TH2F* remove_hot_pixels(TH2F* h, const char* newname = "h_remove_hot_pixels", double nsigma = 5.0) {
    if (!h) return nullptr;

    cout << " --- remove_hot_pixels" << endl;

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    double sum = 0.0;
    double sum2 = 0.0;
    int n = 0;

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double v = h->GetBinContent(ix, iy);
            sum += v;
            sum2 += v * v;
            n++;
        }
    }

    double mean = sum / n;
    double rms = sqrt(sum2 / n - mean * mean);
    double threshold = mean + nsigma * rms;

    TH2F* hout = (TH2F*)h->Clone(newname);
    hout->SetDirectory(0);

    for (int ix = 2; ix < nx; ++ix) {
        for (int iy = 2; iy < ny; ++iy) {
            double val = h->GetBinContent(ix, iy);

            if (val > threshold) {
                double local_sum = 0.0;
                int local_n = 0;

                for (int dx = -1; dx <= 1; ++dx) {
                    for (int dy = -1; dy <= 1; ++dy) {
                        if (dx == 0 && dy == 0) continue;
                        local_sum += h->GetBinContent(ix + dx, iy + dy);
                        local_n++;
                    }
                }

                double local_mean = local_sum / local_n;
                hout->SetBinContent(ix, iy, local_mean);
            }
        }
    }

    return hout;
}

// ============================================================
// Dispatcher: applica una operazione in base alla parola chiave
// ============================================================
TH2F* apply_operation(TH2F* h, const string& op_name, int index) {
    if (!h) return nullptr;

    string hname = "h_" + op_name + "_" + to_string(index);

    if (op_name == "contrast_stretching") {
        return contrast_stretching(h, hname.c_str());
    }
    else if (op_name == "remove_column_bias_median") {
        return remove_column_bias_median(h, hname.c_str());
    }
    else if (op_name == "remove_column_bias_fit") {
        return remove_column_bias_fit(h, hname.c_str());
    }
    else if (op_name == "remove_hot_pixels") {
        return remove_hot_pixels(h, hname.c_str());
    }
    else {
        cerr << "Operazione sconosciuta: " << op_name << endl;
        return nullptr;
    }
}

// ============================================================
// Utility per layout canvas
// ============================================================
void compute_canvas_grid(int npads, int& ncols, int& nrows) {
    ncols = ceil(sqrt(npads));
    nrows = ceil((double)npads / ncols);
}

// ============================================================
// Funzione principale
// filename  = file root
// ops       = stringa tipo "contrast_stretching,remove_bias"
// histname  = opzionale
// ============================================================
void process_th2f(const char* filename,
                  const char* ops = "",
                  const char* histname = "") {
    gStyle->SetOptStat(0);

    TH2F* h_original = load_th2f(filename, histname);
    if (!h_original) return;

    vector<TH2F*> histograms;
    vector<string> titles;

    histograms.push_back(h_original);
    titles.push_back("original");

    vector<string> operations = split_operations(ops);

    for (size_t i = 0; i < operations.size(); ++i) {
        TH2F* hproc = apply_operation(h_original, operations[i], (int)i);
        if (hproc) {
            histograms.push_back(hproc);
            titles.push_back(operations[i]);
        }
    }

    int npads = histograms.size();
    int ncols, nrows;
    compute_canvas_grid(npads, ncols, nrows);

    TCanvas* c1 = new TCanvas("c1", "TH2F processing", 1800, 1080);
    c1->Divide(ncols, nrows);

    for (int i = 0; i < npads; ++i) {
        c1->cd(i + 1);
        histograms[i]->SetTitle(titles[i].c_str());
        histograms[i]->Draw("COLZ");
    }

    c1->Update();
}