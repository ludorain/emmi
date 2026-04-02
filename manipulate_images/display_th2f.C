#include <iostream>
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TClass.h"

void display_th2f(const char* filename, const char* histname = "") {

    // Apri file
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
        return;
    }

    TH2F* hist = nullptr;

    // Caso 1: nome fornito
    if (strlen(histname) > 0) {
        file->GetObject(histname, hist);
        if (!hist) {
            std::cerr << "Errore: TH2F \"" << histname << "\" non trovato." << std::endl;
            return;
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
            return;
        }
    }

    // Canvas
    TCanvas* c1 = new TCanvas("c1", "TH2F display", 800, 600);
    hist->Draw("COLZ");
    c1->Update();
}