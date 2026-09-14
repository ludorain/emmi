// ============================================================
// lum_vs_V_comparison.C
//
// Confronto della luminosità vs overvoltage per il medesimo 
// global hotspot utilizzando quattro scelte di raggio di integrazione:
//   - R best
//   - R = 15 px
//   - R = 20 px
//   - R = 25 px
//
// Modello di fit:
//   Lum = A * V_over^C
// ============================================================

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TAxis.h"

using namespace std;

// Struttura per memorizzare una singola riga di dati dal CSV
struct VRow {
    int spot = -1;
    double x = 0.0;
    double y = 0.0;
    double luminosity = 0.0;
    double error = 0.0;
    double T = 0.0;
    double v = 0.0;
    double v_fin = std::numeric_limits<double>::quiet_NaN();
    double integration_radius = std::numeric_limits<double>::quiet_NaN();
};

// Struttura per memorizzare un intero dataset associato a un determinato raggio
struct VDataset {
    string label;
    string filename;
    int color;
    int marker;
    double nominal_radius;
    map<int, vector<VRow>> spots; // Dati raggruppati per spot ID
};

// Struttura per raccogliere i risultati del fit (rinominato il parametro B in C)
struct VFitSummary {
    TF1* func = nullptr;
    bool ok = false;
    double C = std::numeric_limits<double>::quiet_NaN();
    double C_err = std::numeric_limits<double>::quiet_NaN();
};

// Funzione di utilità per rimuovere spazi bianchi all'inizio e alla fine di una stringa
static string trim_copy(string s) {
    const char* ws = " \t\r\n";
    size_t first = s.find_first_not_of(ws);
    if (first == string::npos) return "";
    size_t last = s.find_last_not_of(ws);
    return s.substr(first, last - first + 1);
}

// Suddivide una riga CSV in base alla virgola
static vector<string> split_csv_simple(const string& line) {
    vector<string> fields;
    string field;
    stringstream ss(line);
    while (getline(ss, field, ',')) fields.push_back(trim_copy(field));
    return fields;
}

// Trova l'indice della colonna a partire dalle intestazioni (permette di avere alias multipli)
static int find_column(const map<string, int>& columns,
                       const vector<string>& aliases,
                       bool required = true) {
    for (const auto& name : aliases) {
        auto it = columns.find(name);
        if (it != columns.end()) return it->second;
    }
    if (required) {
        cerr << "Missing required CSV column. Accepted names:";
        for (const auto& a : aliases) cerr << " " << a;
        cerr << endl;
    }
    return -1;
}

// Cerca di parsare una stringa come numero double (gestendo eccezioni ed errori)
static bool parse_double(const string& text, double& value) {
    try {
        size_t pos = 0;
        value = stod(text, &pos);
        return pos > 0 && std::isfinite(value);
    } catch (...) {
        return false;
    }
}

// Legge e processa un singolo file CSV, mappando ogni riga allo spot corrispondente
static map<int, vector<VRow>> read_v_csv(const string& filename) {
    map<int, vector<VRow>> result;
    ifstream fin(filename);

    if (!fin.is_open()) {
        cerr << "ERROR: cannot open " << filename << endl;
        return result;
    }

    string header_line;
    if (!getline(fin, header_line)) {
        cerr << "ERROR: empty CSV " << filename << endl;
        return result;
    }

    vector<string> header = split_csv_simple(header_line);
    map<string, int> columns;
    for (size_t i = 0; i < header.size(); ++i) columns[header[i]] = int(i);

    int i_spot = find_column(columns, {"spot", "Spot_ID_global", "global_spot_id"});
    int i_x    = find_column(columns, {"x"});
    int i_y    = find_column(columns, {"y"});
    int i_lum  = find_column(columns, {"luminosity"});
    int i_err  = find_column(columns, {"error"});
    int i_T    = find_column(columns, {"T"});
    int i_v    = find_column(columns, {"v"});
    int i_vfin = find_column(columns, {"v_fin", "V_over", "overvoltage"});
    int i_r    = find_column(columns, {"integration_radius", "integration_radius_pix", "radius"}, false);

    if (i_spot < 0 || i_x < 0 || i_y < 0 || i_lum < 0 ||
        i_err < 0 || i_T < 0 || i_v < 0 || i_vfin < 0) {
        return {};
    }

    string line;
    while (getline(fin, line)) {
        if (trim_copy(line).empty()) continue;
        vector<string> f = split_csv_simple(line);

        int max_index = std::max({i_spot, i_x, i_y, i_lum, i_err, i_T, i_v, i_vfin, i_r});
        if ((int)f.size() <= max_index) continue;

        try {
            VRow r;
            r.spot = stoi(f[i_spot]);
            if (!parse_double(f[i_x], r.x)) continue;
            if (!parse_double(f[i_y], r.y)) continue;
            if (!parse_double(f[i_lum], r.luminosity)) continue;
            if (!parse_double(f[i_err], r.error)) continue;
            if (!parse_double(f[i_T], r.T)) continue;
            if (!parse_double(f[i_v], r.v)) continue;
            if (!parse_double(f[i_vfin], r.v_fin)) continue;
            if (i_r >= 0) {
                double tmp;
                if (parse_double(f[i_r], tmp)) r.integration_radius = tmp;
            }
            result[r.spot].push_back(r);
        } catch (...) {
            continue;
        }
    }

    // Ordina i dati di ciascuno spot per overvoltage (v_fin) crescente
    for (auto& kv : result) {
        sort(kv.second.begin(), kv.second.end(),
             [](const VRow& a, const VRow& b) { return a.v_fin < b.v_fin; });
    }

    return result;
}

// Utility per formattare la label "fase"
static string pretty_phase(const string& phase) {
    if (phase == "before_annealing") return "Before annealing";

    const string prefix = "annealing_T=";
    size_t p = phase.find(prefix);
    size_t h = phase.find("_h=");
    if (p == 0 && h != string::npos) {
        string t = phase.substr(prefix.size(), h - prefix.size());
        string hours = phase.substr(h + 3);
        return "Annealing: " + t + " #circC, " + hours + " h";
    }

    return phase;
}

// Utility per rendere il nome file sicuro
static string safe_name(string s) {
    for (char& c : s) {
        if (!(isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.')) c = '_';
    }
    return s;
}

// Trova un raggio d'integrazione valido presente nei dati, se manca usa il fallback del dataset
static double representative_radius(const vector<VRow>& rows, double fallback) {
    for (const auto& r : rows) {
        if (std::isfinite(r.integration_radius)) return r.integration_radius;
    }
    return fallback;
}

//-------------------------------------------------------------------------
// FIT FUNZIONE POTENZA
// Effettua un fit robusto del tipo Lum = A * V_over^C a step.
//-------------------------------------------------------------------------
static VFitSummary make_power_fit(const vector<VRow>& rows,
                                  const string& name,
                                  int color) {
    VFitSummary out;
    vector<VRow> valid;
    for (const auto& r : rows) {
        if (std::isfinite(r.v_fin) && r.v_fin > 0.0 &&
            std::isfinite(r.luminosity) && std::isfinite(r.error)) {
            valid.push_back(r);
        }
    }
    if (valid.size() < 2) return out;

    double xmin = valid.front().v_fin;
    double xmax = valid.back().v_fin;
    if (xmax <= xmin) return out;

    vector<double> xv, yv, exv, eyv;
    for (const auto& r : valid) {
        xv.push_back(r.v_fin);
        yv.push_back(r.luminosity);
        exv.push_back(0.0);       // Niente errore in x
        eyv.push_back(r.error);   // Errore in luminosità fornito
    }

    TGraphErrors tmp(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
    
    // Inizializzazione della funzione e dei nomi parametri
    TF1* f = new TF1(name.c_str(), "[0]*TMath::Power(x,[1])", xmin, xmax);
    f->SetParNames("A", "C");
    
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(1);
    
    // ------- STRATEGIA DI FIT ROBUSTO (Step-by-step) -------
    // Step 1: Fissiamo l'esponente C=2.0 e fittiamo A
    f->FixParameter(1, 2.0);
    // Un buon punto di partenza per A usando l'ultimo punto dei dati:
    f->SetParameter(0, yv.back() / std::pow(xv.back(), 2.0));
    tmp.Fit(f, "Q0"); // Q = Quiet, 0 = Non disegnare
    
    // Step 2: Fissiamo A al valore appena trovato, e rilasciamo C per fittarlo
    double est_A = f->GetParameter(0);
    f->FixParameter(0, est_A);
    f->ReleaseParameter(1);
    f->SetParameter(1, 2.0); // Reimposta il guess di C a 2
    tmp.Fit(f, "Q0");
    
    // Step 3: Rilasciamo entrambi i parametri e procediamo con il fit finale
    f->ReleaseParameter(0);
    tmp.Fit(f, "QRSN"); // R = Rispettare Range, S = Ritorna Risultato, N = Non disegnare
    
    out.func = f;
    out.ok = true;
    out.C = f->GetParameter(1);       // Parametro 1 ora è "C"
    out.C_err = f->GetParError(1);
    return out;
}

// ==========================================================
// FUNZIONE PRINCIPALE DELLA MACRO
// ==========================================================
void lum_vs_V_comparison(
    const char* best_file,
    const char* r15_file,
    const char* r20_file,
    const char* r25_file,
    const char* output_dir,
    const char* sensor,
    const char* fixed_condition,
    const char* phase,
    int selected_spot = -1
) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Definizione dei 4 set di dati con relativi colori e marker
    vector<VDataset> datasets = {
        {"R best",    best_file, kBlack,   20, std::numeric_limits<double>::quiet_NaN(), {}},
        {"R = 15 px", r15_file,  kBlue+1,  21, 15.0, {}},
        {"R = 20 px", r20_file,  kRed+1,   22, 20.0, {}},
        {"R = 25 px", r25_file,  kGreen+2, 23, 25.0, {}}
    };

    // Lettura dei file e identificazione degli Spot
    set<int> all_spots;
    for (auto& d : datasets) {
        d.spots = read_v_csv(d.filename);
        if (d.spots.empty()) {
            cerr << "WARNING: no valid rows read from " << d.filename << endl;
        }
        for (const auto& kv : d.spots) all_spots.insert(kv.first);
    }

    if (all_spots.empty()) {
        cerr << "ERROR: no global hotspot IDs found in the input files." << endl;
        return;
    }

    // Filtraggio: processa tutti gli spot oppure uno specifico se requested
    vector<int> spots_to_plot;
    if (selected_spot >= 0) {
        if (all_spots.count(selected_spot) == 0) {
            cerr << "ERROR: global hotspot ID " << selected_spot
                 << " is not present in any input CSV." << endl;
            return;
        }
        spots_to_plot.push_back(selected_spot);
    } else {
        spots_to_plot.assign(all_spots.begin(), all_spots.end());
    }

    // Creazione della cartella di output
    gSystem->mkdir(output_dir, kTRUE);

    /* 
    --------------- CANVAS FOR COMPARISON OF DIFFERENT RADII ---------------
    */

    for (int spot : spots_to_plot) {
        TCanvas* c = new TCanvas(
            Form("c_V_comparison_spot_%d", spot),
            Form("Spot %d - luminosity vs overvoltage", spot),
            1200, 800);

        c->SetLeftMargin(0.13);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.12);
        c->SetTopMargin(0.10);

        TMultiGraph* mg = new TMultiGraph();
        
        // La legenda spostata e dimensionata in alto a sinistra (Top-Left) 
        // e ingrandita in orizzontale per ospitare raggio e risultati del fit (valore C + errore).
        TLegend* leg = new TLegend(0.15, 0.60, 0.60, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        bool any_graph = false;
        double x_coord = std::numeric_limits<double>::quiet_NaN();
        double y_coord = std::numeric_limits<double>::quiet_NaN();

        vector<TF1*> fit_funcs;
        vector<double> radius_x, radius_ex, fitC_y, fitC_ey;

        for (size_t id = 0; id < datasets.size(); ++id) {
            auto it = datasets[id].spots.find(spot);
            if (it == datasets[id].spots.end() || it->second.empty()) continue;

            const vector<VRow>& rows = it->second;
            if (!std::isfinite(x_coord)) {
                x_coord = rows.front().x;
                y_coord = rows.front().y;
            }

            vector<double> xv, yv, exv, eyv;
            for (const auto& r : rows) {
                if (!std::isfinite(r.v_fin) || !std::isfinite(r.luminosity) ||
                    !std::isfinite(r.error)) continue;
                xv.push_back(r.v_fin);
                yv.push_back(r.luminosity);
                exv.push_back(0.0);
                eyv.push_back(r.error);
            }
            if (xv.empty()) continue;

            TGraphErrors* gr = new TGraphErrors(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
            gr->SetName(Form("gr_V_spot_%d_dataset_%zu", spot, id));
            gr->SetMarkerColor(datasets[id].color);
            gr->SetLineColor(datasets[id].color);
            gr->SetMarkerStyle(datasets[id].marker);
            gr->SetMarkerSize(1.25);
            gr->SetLineWidth(2);

            // Aggiungiamo con "PE" per disegnare le barre d'errore (Point + Errors)
            mg->Add(gr, "PE");

            // Raggio utilizzato (Reale misurato oppure nominale per fallback)
            double current_radius = representative_radius(rows, datasets[id].nominal_radius);

            // Generiamo l'etichetta base
            string legend_label = datasets[id].label;
            if (id == 0 && std::isfinite(current_radius)) {
                legend_label = Form("R best = %.1f px", current_radius);
            }

            // Eseguiamo il fit
            VFitSummary fit = make_power_fit(rows, Form("fit_V_spot_%d_dataset_%zu", spot, id), datasets[id].color);

            // Se il fit ha avuto successo formattiamo la legenda aggiungendo il param. C
            if (fit.ok && fit.func) {
                fit_funcs.push_back(fit.func);
                
                legend_label += Form(" (C = %.2f #pm %.2f)", fit.C, fit.C_err);

                if (std::isfinite(current_radius) &&
                    std::isfinite(fit.C) &&
                    std::isfinite(fit.C_err)) {

                    radius_x.push_back(current_radius);
                    radius_ex.push_back(0.0);
                    fitC_y.push_back(fit.C);
                    fitC_ey.push_back(fit.C_err);
                }
            }
            
            // Aggiungiamo la stringa composta in legenda
            leg->AddEntry(gr, legend_label.c_str(), "lep");
            any_graph = true;
        }

        if (!any_graph) {
            delete leg;
            delete mg;
            delete c;
            continue;
        }

        string title = Form(
            "%s - Global hotspot %d - %s - %s;Overvoltage (V);Luminosity",
            sensor, spot, pretty_phase(phase).c_str(), fixed_condition
        );
        if (std::isfinite(x_coord) && std::isfinite(y_coord)) {
            title = Form(
                "%s - Global hotspot %d (x=%.2f, y=%.2f) - %s - %s;Overvoltage (V);Luminosity",
                sensor, spot, x_coord, y_coord,
                pretty_phase(phase).c_str(), fixed_condition
            );
        }

        mg->SetTitle(title.c_str());
        
        // Disegniamo il TMultiGraph con assi (A), marker (P), e barre di errore (E)
        mg->Draw("APE");
        
        mg->GetXaxis()->SetTitleSize(0.045);
        mg->GetYaxis()->SetTitleSize(0.045);
        mg->GetXaxis()->SetLabelSize(0.040);
        mg->GetYaxis()->SetLabelSize(0.040);
        
        // Traccia le curve fittate
        for (TF1* fit : fit_funcs) fit->Draw("SAME");
        leg->Draw();
        
        c->Modified();
        c->Update();

        string outfile = string(output_dir) + "/" +
            safe_name(Form("%s_%s_%s_lum_vs_V_spot%d_radii_comparison.png", sensor, fixed_condition, phase, spot));
        c->SaveAs(outfile.c_str());
        delete c;

        // Se sono stati estratti i C del fit per ogni raggio, crea il grafico andamento C vs R
        if (!radius_x.empty()) {
            vector<size_t> order(radius_x.size());
            for (size_t i = 0; i < order.size(); ++i) order[i] = i;
            sort(order.begin(), order.end(), [&](size_t a, size_t b) { return radius_x[a] < radius_x[b]; });

            vector<double> rx, rex, cy, cey;
            for (size_t idx : order) {
                rx.push_back(radius_x[idx]);
                rex.push_back(radius_ex[idx]);
                cy.push_back(fitC_y[idx]);
                cey.push_back(fitC_ey[idx]);
            }

            TCanvas* cC = new TCanvas(Form("c_V_fitC_vs_R_spot_%d", spot), Form("Spot %d - fit parameter C vs radius", spot), 900, 700);
            cC->SetLeftMargin(0.13);
            cC->SetRightMargin(0.05);
            cC->SetBottomMargin(0.12);
            cC->SetTopMargin(0.10);

            TGraphErrors* grC = new TGraphErrors(rx.size(), rx.data(), cy.data(), rex.data(), cey.data());
            grC->SetMarkerStyle(20);
            grC->SetMarkerSize(1.3);
            grC->SetLineWidth(2);
            grC->SetTitle(Form("%s - Global hotspot %d - %s - %s;Integration radius (px);Fit parameter C", sensor, spot, pretty_phase(phase).c_str(), fixed_condition));
            
            grC->Draw("APL");
            cC->Modified();
            cC->Update();

            string outfileC = string(output_dir) + "/" + safe_name(Form("%s_%s_%s_fitC_vs_radius_spot%d.png", sensor, fixed_condition, phase, spot));
            cC->SaveAs(outfileC.c_str());
            delete cC;
        }
    }
}