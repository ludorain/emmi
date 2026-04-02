#include <iostream>
#include <TMath.h>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TClass.h"

// ============================================================
// Funzione 1: contrast_stretching
// percentile 1% - 99%, riscalato in [0,1]
// ============================================================
TH2F* contrast_stretching(TH2F* hist, const char* newname = "hist_contrast")
{
    std::cout << " --- contrast_stretching" << std::endl;

    if (!hist) {
        std::cerr << "Errore: istogramma nullo." << std::endl;
        return nullptr;
    }

    int nx = hist->GetNbinsX();
    int ny = hist->GetNbinsY();

    // Raccolgo tutti i contenuti dei bin
    std::vector<double> values;
    values.reserve(nx * ny);

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            values.push_back(hist->GetBinContent(ix, iy));
        }
    }

    if (values.empty()) {
        std::cerr << "Errore: nessun bin trovato." << std::endl;
        return nullptr;
    }

    // Ordino i valori per ricavare i percentili
    std::sort(values.begin(), values.end());

    int i1  = static_cast<int>(0.01 * (values.size() - 1));
    int i99 = static_cast<int>(0.99 * (values.size() - 1));

    double pmin = values[i1];
    double pmax = values[i99];

    if (pmax == pmin) {
        std::cerr << "Errore: pmax == pmin, impossibile riscalare." << std::endl;
        return nullptr;
    }

    // Clono l'istogramma di input per mantenere stessa geometria/assi
    TH2F* out = (TH2F*)hist->Clone(newname);
    out->Reset();

    // Riscalatura lineare tra -1 e 1
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double val = hist->GetBinContent(ix, iy);
            double newval;

            if (val <= pmin) {
                newval = -1.0;
            } else if (val >= pmax) {
                newval = 1.0;
            } else {
                newval = -1.0 + 2.0 * (val - pmin) / (pmax - pmin);
            }

            out->SetBinContent(ix, iy, newval);
        }
    }

    return out;
}


// ============================================================
// Funzione principale del codice
// ============================================================

void comparison_contrast(const char* filename1, const char* filename2, const char* histname1 = "", const char* histname2 = "") {

// ============================================================
// Aprire i file e caricare i TH2F
// ============================================================
//Guarda come creare gli istogrammi in maniera più formale (slide Arcelli <3) con i range giusti
/* TH2F *hist_root =new TH2F(hist_root", "Istogramma da Root", 100, 0, 1, 100, 0, 1);
    TH2F *hist_python = new TH2F("hist_python", "Istogramma da Python", 100, 0, 1, 100, 0, 1);
*/

    TH2F* hist_root = nullptr;
    TH2F* hist_python = nullptr;

    //Open first file

    TFile* file1 = TFile::Open(filename1, "READ");
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il primo file " << filename1 << std::endl;
        return;
    }

    // Caso 1.1: nome fornito
    if (strlen(histname1) > 0) {
        file1->GetObject(histname1, hist_root);
        if (!hist_root) {
            std::cerr << "Errore: TH2F \"" << histname1 << "\" non trovato." << std::endl;
            return;
        }
    }      
    // Caso 2.1: cerca automaticamente
    else {
        TIter next(file1->GetListOfKeys());
        TKey* key1;

        while ((key1 = (TKey*)next())) {
            TObject* obj1 = key1->ReadObj();
            if (obj1->InheritsFrom("TH2F")) {
                hist_root = (TH2F*)obj1;
                break;
            }
        }

        if (!hist_root) {
            std::cerr << "Errore: nessun TH2F trovato nel file." << std::endl;
            return;
        }
    }

    //Open second file

        TFile* file2 = TFile::Open(filename2, "READ");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il secondo file " << filename2 << std::endl;
        return;
    }

    //caso 1.2: nome fornito
    if (strlen(histname2) > 0) {
        file2->GetObject(histname2, hist_python);
        if (!hist_python) {
            std::cerr << "Errore: TH2F \"" << histname2 << "\" non trovato." << std::endl;
            return;
        }
    }
    // Caso 2.2: cerca automaticamente
    else {
        TIter next(file2->GetListOfKeys());
        TKey* key2;

        while ((key2 = (TKey*)next())) {
            TObject* obj2 = key2->ReadObj();
            if (obj2->InheritsFrom("TH2F")) {
                hist_python = (TH2F*)obj2;
                break;
            }
        }

        if (!hist_python) {
            std::cerr << "Errore: nessun TH2F trovato nel file." << std::endl;
            return;
        }
    }

    //creating the histograms to be plotted

    TH2F* h1 = contrast_stretching(hist_root, "h1_contrast_root");
    //TH2F* h2 = (TH2F*)hist_python->Clone("h2_contrast_python"); //evitabile, messo solo per chiarezza personale

    TH2F *hist_difference = new TH2F("hist_difference", "Differenza", h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(),
                                            h1->GetNbinsY(), h1->GetYaxis()->GetXmin(), h1->GetYaxis()->GetXmax());
    hist_difference->Add(h1, hist_python, 1.0, -1.0); //h1 - h2                                      


   //Canvas, la finestra grafica (divisa in 4).... 
  TCanvas *ccontrast = new TCanvas("ccontrast", "Contrast stretching ROOT vs Python",10,30,1800,1080);
  ccontrast->Divide(2,2); //divisa in 4 pad
  
   
  //posizionamento su prima pad
  ccontrast->cd(1);
  h1->SetTitle("Istogramma da ROOT (contrast stretching)");
  h1->Draw("COLZ");


  //posizionamento su seconda pad
  ccontrast->cd(2);
  hist_python->SetTitle("Istogramma da Python (contrast stretching)");
  hist_python->Draw("COLZ");
 
  //posizionamento su terza pad
  ccontrast->cd(3);
  hist_difference->SetTitle("Differenza (ROOT - Python) opz COLZ");
  hist_difference->Draw("COLZ");
 
  //posizionamento su quarta pad
  ccontrast->cd(4); 
  TH2F* h4 = (TH2F*)hist_difference->Clone("h4_contrast_difference");
  h4->SetTitle("Differenza (ROOT - Python) opz HISTO");
  h4->Draw("histo");
  
  cout << "*-----------------Siamo arrivati alla fine-----------------------------------------*" <<endl; 

/*  
// Debug: stampa info istogrammi
  std::cout << "h1 ptr     = " << h1
          << "  name = " << h1->GetName()
          << "  title = " << h1->GetTitle() << std::endl;

    
  std::cout << "hist_python ptr = " << hist_python
          << "  name = " << hist_python->GetName()
          << "  title = " << hist_python->GetTitle() << std::endl;

*/
}