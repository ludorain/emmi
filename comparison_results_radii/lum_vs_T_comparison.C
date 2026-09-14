// ============================================================
// lum_vs_T_comparison.C
//
// Compare luminosity vs temperature for the SAME global hotspot
// using four integration-radius choices:
//   - R best
//   - R = 15 px
//   - R = 20 px
//   - R = 25 px
//
// If selected_spot < 0, one plot is produced for every global
// hotspot ID present in at least one input CSV.
//
// Fit model:
//   Lum(T) = A * exp(B*T)
//
// Additional output:
//   fit parameter B vs integration radius
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
struct TRowComparison {
    int spot = -1;
    double x = 0.0;
    double y = 0.0;
    double luminosity = 0.0;
    double error = 0.0;
    double T = 0.0;
    double v = 0.0;
    double integration_radius = std::numeric_limits<double>::quiet_NaN();
};

// Struttura per rappresentare un set di dati (es. corrispondente a un certo raggio)
struct TDatasetComparison {
    string label;
    string filename;
    int color;
    int marker;
    double nominal_radius;
    map<int, vector<TRowComparison>> spots; // Raggruppa i dati per spot ID
};

// Struttura per salvare i risultati del fit esponenziale
struct TFitSummary {
    TF1* func = nullptr;
    bool ok = false;
    double B = std::numeric_limits<double>::quiet_NaN();
    double B_err = std::numeric_limits<double>::quiet_NaN();
};

// Funzione helper per rimuovere spazi bianchi all'inizio e alla fine di una stringa
static string trim_copy_T(string s) {
    const char* ws = " \t\r\n";
    size_t first = s.find_first_not_of(ws);
    if (first == string::npos) return "";
    size_t last = s.find_last_not_of(ws);
    return s.substr(first, last - first + 1);
}

// Suddivide una linea CSV in campi separati da virgola
static vector<string> split_csv_simple_T(const string& line) {
    vector<string> fields;
    string field;
    stringstream ss(line);
    while (getline(ss, field, ',')) fields.push_back(trim_copy_T(field));
    return fields;
}

// Trova l'indice di una colonna basandosi su una lista di alias accettati
static int find_column_T(const map<string, int>& columns,
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

// Cerca di estrarre un double in modo sicuro
static bool parse_double_T(const string& text, double& value) {
    try {
        size_t pos = 0;
        value = stod(text, &pos);
        return pos > 0 && std::isfinite(value);
    } catch (...) {
        return false;
    }
}

// Legge il CSV in input e mappa le righe corrispondenti per spot
static map<int, vector<TRowComparison>> read_t_csv(const string& filename) {
    map<int, vector<TRowComparison>> result;
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

    vector<string> header = split_csv_simple_T(header_line);
    map<string, int> columns;
    for (size_t i = 0; i < header.size(); ++i) columns[header[i]] = int(i);

    // Identificazione delle colonne richieste (o opzionali)
    int i_spot = find_column_T(columns, {"spot", "Spot_ID_global", "global_spot_id"});
    int i_x    = find_column_T(columns, {"x"});
    int i_y    = find_column_T(columns, {"y"});
    int i_lum  = find_column_T(columns, {"luminosity"});
    int i_err  = find_column_T(columns, {"error"});
    int i_T    = find_column_T(columns, {"T"});
    int i_v    = find_column_T(columns, {"v"});
    int i_r    = find_column_T(columns, {"integration_radius", "integration_radius_pix", "radius"}, false);

    if (i_spot < 0 || i_x < 0 || i_y < 0 || i_lum < 0 ||
        i_err < 0 || i_T < 0 || i_v < 0) {
        return {};
    }

    string line;
    while (getline(fin, line)) {
        if (trim_copy_T(line).empty()) continue;
        vector<string> f = split_csv_simple_T(line);

        int max_index = std::max({i_spot, i_x, i_y, i_lum, i_err, i_T, i_v, i_r});
        if ((int)f.size() <= max_index) continue;

        try {
            TRowComparison r;
            r.spot = stoi(f[i_spot]);
            if (!parse_double_T(f[i_x], r.x)) continue;
            if (!parse_double_T(f[i_y], r.y)) continue;
            if (!parse_double_T(f[i_lum], r.luminosity)) continue;
            if (!parse_double_T(f[i_err], r.error)) continue;
            if (!parse_double_T(f[i_T], r.T)) continue;
            if (!parse_double_T(f[i_v], r.v)) continue;
            if (i_r >= 0) {
                double tmp;
                if (parse_double_T(f[i_r], tmp)) r.integration_radius = tmp;
            }
            result[r.spot].push_back(r);
        } catch (...) {
            continue;
        }
    }

    // Ordina i dati di ogni spot in base alla temperatura (dal più freddo al più caldo)
    for (auto& kv : result) {
        sort(kv.second.begin(), kv.second.end(),
             [](const TRowComparison& a, const TRowComparison& b) { return a.T < b.T; });
    }

    return result;
}

// Pulisce e formatta il nome della fase per i titoli dei grafici
static string pretty_phase_T(const string& phase) {
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

// Rimuove caratteri non ammessi per i nomi di file generati
static string safe_name_T(string s) {
    for (char& c : s) {
        if (!(isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.')) c = '_';
    }
    return s;
}

// Restituisce il primo raggio di integrazione valido trovato
static double representative_radius_T(const vector<TRowComparison>& rows, double fallback) {
    for (const auto& r : rows) {
        if (std::isfinite(r.integration_radius)) return r.integration_radius;
    }
    return fallback;
}

// Stima dei parametri iniziali (guess) per agevolare la convergenza del fit
static void estimate_exp_parameters_T(const vector<TRowComparison>& rows, double& A0, double& B0) {
    vector<TRowComparison> positive;
    for (const auto& r : rows) {
        if (r.luminosity > 0.0 && std::isfinite(r.luminosity)) positive.push_back(r);
    }
    sort(positive.begin(), positive.end(), [](const TRowComparison& a, const TRowComparison& b) { return a.T < b.T; });

    if (positive.size() >= 2) {
        const auto& r1 = positive.front();
        const auto& r2 = positive.back();
        if (fabs(r2.T - r1.T) > 1e-12) {
            B0 = (log(r2.luminosity) - log(r1.luminosity)) / (r2.T - r1.T);
            A0 = exp(log(r1.luminosity) - B0 * r1.T);
            return;
        }
    }

    if (positive.size() == 1) {
        A0 = positive.front().luminosity;
        B0 = 0.0;
    } else {
        A0 = 1.0;
        B0 = 0.0;
    }
}

// Effettua il fit esponenziale sui dati passati in input
static TFitSummary make_exp_fit_T(const vector<TRowComparison>& rows,
                                  const string& name,
                                  int color) {
    TFitSummary out;
    vector<TRowComparison> valid;
    for (const auto& r : rows) {
        if (std::isfinite(r.T) && std::isfinite(r.luminosity) && std::isfinite(r.error)) {
            valid.push_back(r);
        }
    }
    if (valid.size() < 3) return out;

    double xmin = valid.front().T;
    double xmax = valid.back().T;
    if (xmax <= xmin) return out;

    double span = xmax - xmin;
    double fit_min = xmin - 0.10 * span;
    double fit_max = xmax + 0.10 * span;

    vector<double> xv, yv, exv, eyv;
    for (const auto& r : valid) {
        xv.push_back(r.T);
        yv.push_back(r.luminosity);
        exv.push_back(0.0);
        eyv.push_back(r.error);
    }

    double A0 = 1.0, B0 = 0.0;
    estimate_exp_parameters_T(valid, A0, B0);

    TGraphErrors tmp(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
    TF1* f = new TF1(name.c_str(), "[0]*exp([1]*x)", fit_min, fit_max);
    f->SetParNames("A", "B");
    
    // Fit più robusto: fissiamo parametri e limiti iniziali coerenti con la fisica attesa
    f->SetParameters(A0, B0);
    f->SetParLimits(0, 1, 3000); // A (luminosità a T=0) deve essere positiva
    f->SetParLimits(1, -1.0, 2.0); // B ragionevolmente confinato per evitare divergenze numeriche

    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetLineStyle(2);
    
    // Opzioni di fit avanzate:
    // Q = Quiet (silenzioso)
    // R = Range specificato dalla funzione
    // S = Salva il risultato nel TFitResultPtr
    // M = Minos (migliora i fit calcolando meglio minimi ed errori)
    // E = Esegue una migliore stima degli errori se possibile
    // 0 = Non disegnare automaticamente (lo facciamo noi dopo)
    tmp.Fit(f, "QRSME0");

    out.func = f;
    out.ok = true;
    out.B = f->GetParameter(1);
    out.B_err = f->GetParError(1);
    return out;
}

// Entry point della macro
void lum_vs_T_comparison(
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

    // Preparazione dei dataset
    vector<TDatasetComparison> datasets = {
        {"R best",    best_file, kBlack,   20, std::numeric_limits<double>::quiet_NaN(), {}},
        {"R = 15 px", r15_file,  kBlue+1,  21, 15.0, {}},
        {"R = 20 px", r20_file,  kRed+1,   22, 20.0, {}},
        {"R = 25 px", r25_file,  kGreen+2, 23, 25.0, {}}
    };

    // Estrazione dei dati
    set<int> all_spots;
    for (auto& d : datasets) {
        d.spots = read_t_csv(d.filename);
        if (d.spots.empty()) {
            cerr << "WARNING: no valid rows read from " << d.filename << endl;
        }
        for (const auto& kv : d.spots) all_spots.insert(kv.first);
    }

    if (all_spots.empty()) {
        cerr << "ERROR: no global hotspot IDs found in the input files." << endl;
        return;
    }

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

    gSystem->mkdir(output_dir, kTRUE);

    // Ciclo su ciascuno spot per produrre i grafici comparativi
    for (int spot : spots_to_plot) {
        TCanvas* c = new TCanvas(Form("c_T_comparison_spot_%d", spot), Form("Spot %d - luminosity vs temperature", spot), 1000, 750);
        c->SetLeftMargin(0.13);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.12);
        c->SetTopMargin(0.10);

        TMultiGraph* mg = new TMultiGraph();
        
        // Legenda riposizionata: in alto a sinistra per evitare la curva esponenziale a destra
        TLegend* leg = new TLegend(0.15, 0.65, 0.50, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);

        bool any_graph = false;
        double x_coord = std::numeric_limits<double>::quiet_NaN();
        double y_coord = std::numeric_limits<double>::quiet_NaN();

        vector<TF1*> fit_funcs;
        vector<double> radius_x, radius_ex, fitB_y, fitB_ey;

        for (size_t id = 0; id < datasets.size(); ++id) {
            auto it = datasets[id].spots.find(spot);
            if (it == datasets[id].spots.end() || it->second.empty()) continue;

            const vector<TRowComparison>& rows = it->second;
            if (!std::isfinite(x_coord)) {
                x_coord = rows.front().x;
                y_coord = rows.front().y;
            }

            vector<double> xv, yv, exv, eyv;
            for (const auto& r : rows) {
                if (!std::isfinite(r.T) || !std::isfinite(r.luminosity) || !std::isfinite(r.error)) continue;
                xv.push_back(r.T);
                yv.push_back(r.luminosity);
                exv.push_back(0.0);
                eyv.push_back(r.error); // Viene catturato l'errore per disegnarlo
            }
            if (xv.empty()) continue;

            TGraphErrors* gr = new TGraphErrors(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
            gr->SetName(Form("gr_T_spot_%d_dataset_%zu", spot, id));
            gr->SetMarkerColor(datasets[id].color);
            gr->SetLineColor(datasets[id].color);
            gr->SetMarkerStyle(datasets[id].marker);
            gr->SetMarkerSize(1.25);
            gr->SetLineWidth(2);

            // Usa "PE" (Points and Errors) in modo da forzare il disegno delle barre di errore
            mg->Add(gr, "PE");
            
            double current_radius = representative_radius_T(rows, datasets[id].nominal_radius);

            // Richiamo funzione di fit esponenziale migliorato
            TFitSummary fit = make_exp_fit_T(rows, Form("fit_T_spot_%d_dataset_%zu", spot, id), datasets[id].color);

            string legend_label = datasets[id].label;
            if (id == 0 && std::isfinite(current_radius)) {
                legend_label = Form("R best = %.1f px", current_radius);
            }
            
            // Se il fit ha successo, aggiorniamo la label della legenda mostrando il parametro B con errore
            if (fit.ok && fit.func) {
                fit_funcs.push_back(fit.func);
                
                legend_label += Form(" (B = %.4f #pm %.4f)", fit.B, fit.B_err);

                if (std::isfinite(current_radius) &&
                    std::isfinite(fit.B) &&
                    std::isfinite(fit.B_err)) {

                    radius_x.push_back(current_radius);
                    radius_ex.push_back(0.0);
                    fitB_y.push_back(fit.B);
                    fitB_ey.push_back(fit.B_err);
                }
            }
            
            // Aggiungiamo il riquadro personalizzato nella legenda
            leg->AddEntry(gr, legend_label.c_str(), "pe");
            any_graph = true;
        }

        if (!any_graph) {
            delete leg;
            delete mg;
            delete c;
            continue;
        }

        string title = Form("%s - Global hotspot %d - %s - %s;Temperature (#circC);Luminosity", sensor, spot, pretty_phase_T(phase).c_str(), fixed_condition);
        if (std::isfinite(x_coord) && std::isfinite(y_coord)) {
            title = Form("%s - Global hotspot %d (x=%.2f, y=%.2f) - %s - %s;Temperature (#circC);Luminosity", sensor, spot, x_coord, y_coord, pretty_phase_T(phase).c_str(), fixed_condition);
        }

        mg->SetTitle(title.c_str());
        mg->Draw("A"); // "A" disegna gli assi del Multigraph, che poi inietterà in automatico i "PE" dei figli.
        mg->GetXaxis()->SetTitleSize(0.045);
        mg->GetYaxis()->SetTitleSize(0.045);
        mg->GetXaxis()->SetLabelSize(0.040);
        mg->GetYaxis()->SetLabelSize(0.040);
        
        for (TF1* fit : fit_funcs) fit->Draw("SAME");
        leg->Draw();
        c->Modified();
        c->Update();

        string outfile = string(output_dir) + "/" + safe_name_T(Form("%s_%s_%s_lum_vs_T_spot%d_radii_comparison.png", sensor, fixed_condition, phase, spot));
        c->SaveAs(outfile.c_str());
        delete c;

        // Se abbiamo salvato dati dei fit rispetto ai raggi, produciamo il grafico riassuntivo per B vs R
        if (!radius_x.empty()) {
            vector<size_t> order(radius_x.size());
            for (size_t i = 0; i < order.size(); ++i) order[i] = i;
            sort(order.begin(), order.end(), [&](size_t a, size_t b) { return radius_x[a] < radius_x[b]; });

            vector<double> rx, rex, by, bey;
            for (size_t idx : order) {
                rx.push_back(radius_x[idx]);
                rex.push_back(radius_ex[idx]); // l'errore lungo X sarà nullo (0)
                by.push_back(fitB_y[idx]);
                bey.push_back(fitB_ey[idx]);   // L'errore del fit su B
            }

            TCanvas* cB = new TCanvas(Form("c_T_fitB_vs_R_spot_%d", spot), Form("Spot %d - fit parameter B vs radius", spot), 900, 700);
            cB->SetLeftMargin(0.13);
            cB->SetRightMargin(0.05);
            cB->SetBottomMargin(0.12);
            cB->SetTopMargin(0.10);

            TGraphErrors* grB = new TGraphErrors(rx.size(), rx.data(), by.data(), rex.data(), bey.data());
            grB->SetMarkerStyle(20);
            grB->SetMarkerSize(1.3);
            grB->SetLineWidth(2);
            grB->SetTitle(Form("%s - Global hotspot %d - %s - %s;Integration radius (px);Fit parameter B", sensor, spot, pretty_phase_T(phase).c_str(), fixed_condition));
            
            // APLE = Axis + Points + Lines + Errors: garantisce che le barre di errore siano visibili
            grB->Draw("APLE"); 
            cB->Modified();
            cB->Update();

            string outfileB = string(output_dir) + "/" + safe_name_T(Form("%s_%s_%s_fitB_vs_radius_spot%d.png", sensor, fixed_condition, phase, spot));
            cB->SaveAs(outfileB.c_str());
            delete cB;
        }
    }
}