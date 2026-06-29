//root -l 'luminosity_vs_phase_T_const.C("merged_max_vfin_global_ID.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <cctype>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLine.h"
#include "TString.h"

using namespace std;

struct Row {
    int spot;
    double x;
    double y;
    double luminosity;
    double error;
    double T;
    double v;
    string phase;
    bool detected;
};

// ------------------------------------------------------------
// Funzioni di utilita'
// ------------------------------------------------------------

string trim(const string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end   = s.find_last_not_of(" \t\r\n");

    if (start == string::npos) return "";
    return s.substr(start, end - start + 1);
}

vector<string> split_csv_line(const string& line) {
    vector<string> fields;
    string field;
    stringstream ss(line);

    while (getline(ss, field, ',')) {
        fields.push_back(trim(field));
    }

    return fields;
}

bool parse_bool(string s) {
    s = trim(s);

    for (char& c : s) {
        c = std::tolower(static_cast<unsigned char>(c));
    }

    if (s == "true"  || s == "1" || s == "yes") return true;
    if (s == "false" || s == "0" || s == "no")  return false;

    cerr << "Warning: valore detected non riconosciuto: '" << s
         << "'. Lo considero false." << endl;

    return false;
}

string safe_name(string s) {
    for (char& c : s) {
        if (c == '=' || c == '/' || c == ' ' || c == '-' || c == '.') {
            c = '_';
        }
    }
    return s;
}

// ------------------------------------------------------------
// Lettura CSV tramite nomi delle colonne
// CSV atteso:
// spot,x,y,luminosity,error,T,v,phase,detected
// ------------------------------------------------------------

vector<Row> read_csv(const char* filename) {
    vector<Row> data;

    ifstream fin(filename);

    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire il file " << filename << endl;
        return data;
    }

    string line;

    if (!getline(fin, line)) {
        cerr << "Errore: file vuoto." << endl;
        return data;
    }

    vector<string> header = split_csv_line(line);
    map<string, int> col;

    for (int i = 0; i < (int)header.size(); i++) {
        col[header[i]] = i;
    }

    vector<string> required = {
        "spot", "x", "y", "luminosity", "error",
        "T", "v", "phase", "detected"
    };

    for (const auto& name : required) {
        if (col.find(name) == col.end()) {
            cerr << "Errore: manca la colonna '" << name << "' nel CSV." << endl;
            return data;
        }
    }

    while (getline(fin, line)) {
        if (line.empty()) continue;

        vector<string> f = split_csv_line(line);

        if ((int)f.size() < (int)header.size()) {
            cerr << "Warning: riga saltata perche' incompleta:\n"
                 << line << endl;
            continue;
        }

        Row r;

        try {
            // spot puo' essere scritto come 0 oppure 0.0.
            r.spot       = (int)llround(stod(f[col["spot"]]));

            r.x          = stod(f[col["x"]]);
            r.y          = stod(f[col["y"]]);
            r.luminosity = stod(f[col["luminosity"]]);
            r.error      = stod(f[col["error"]]);
            r.T          = stod(f[col["T"]]);
            r.v          = stod(f[col["v"]]);
            r.phase      = f[col["phase"]];
            r.detected   = parse_bool(f[col["detected"]]);
        }
        catch (const std::exception& e) {
            cerr << "Warning: errore nella lettura della riga:\n"
                 << line << endl;
            cerr << "Motivo: " << e.what() << endl;
            continue;
        }

        data.push_back(r);
    }

    cout << "Letti " << data.size()
         << " punti dal file " << filename << endl;

    return data;
}

// ------------------------------------------------------------
// Disegna una canvas ratio vs spot
// ------------------------------------------------------------

void draw_ratio_canvas(const vector<double>& xratio,
                       const vector<double>& yratio,
                       const vector<double>& exratio,
                       const vector<double>& eyratio,
                       const string& canvas_name,
                       const string& canvas_title,
                       const string& png_name,
                       bool save_png,
                       bool fixed_y_range = false,
                       double fixed_ymin = 0.0,
                       double fixed_ymax = 2.0) {

    if (xratio.empty()) return;

    double ymin = 1e99;
    double ymax = -1e99;

    if (fixed_y_range) {
        ymin = fixed_ymin;
        ymax = fixed_ymax;
    } else {
        for (int i = 0; i < (int)yratio.size(); i++) {
            ymin = min(ymin, yratio[i] - eyratio[i]);
            ymax = max(ymax, yratio[i] + eyratio[i]);
        }

        if (ymin > 0) ymin *= 0.8;
        else          ymin *= 1.2;

        ymax *= 1.25;

        if (ymax <= ymin) ymax = ymin + 1.0;
    }

    TCanvas* c_ratio = new TCanvas(
        canvas_name.c_str(),
        canvas_title.c_str(),
        1100,
        750
    );

    c_ratio->SetGrid();

    TGraphErrors* gr_ratio = new TGraphErrors(
        (int)xratio.size(),
        xratio.data(),
        yratio.data(),
        exratio.data(),
        eyratio.data()
    );

    gr_ratio->SetTitle(
        "L(annealing_T=100_h=5) / L(before_annealing);spot;ratio"
    );

    gr_ratio->SetMarkerStyle(20);
    gr_ratio->SetMarkerSize(1.2);
    gr_ratio->SetLineWidth(2);

    gr_ratio->Draw("AP");
    gr_ratio->GetYaxis()->SetRangeUser(ymin, ymax);

    c_ratio->Update();

    double xmin_axis = gr_ratio->GetXaxis()->GetXmin();
    double xmax_axis = gr_ratio->GetXaxis()->GetXmax();

    TLine* line_one = new TLine(
        xmin_axis,
        1.0,
        xmax_axis,
        1.0
    );

    line_one->SetLineStyle(2);
    line_one->SetLineWidth(2);
    line_one->Draw("SAME");

    c_ratio->Update();

    if (save_png) {
        c_ratio->SaveAs(png_name.c_str());
    }
}

// ------------------------------------------------------------
// Macro principale
// ------------------------------------------------------------

void luminosity_vs_phase_T_const(const char* csvfile,
                         bool save_png = true,
                         bool use_only_detected = false) {

    gStyle->SetOptStat(0);

    // Numero massimo di spot rappresentati nella stessa canvas.
    const int MAX_SPOTS_PER_CANVAS = 15;

    // Range della seconda canvas finale, quella di focus sul rapporto.
    const double RATIO_FOCUS_YMIN = -0.1;
    const double RATIO_FOCUS_YMAX = 2.0;

    // Soglia per la significance richiesta.
    const double SIGMA_THRESHOLD = 5.0;

    vector<Row> data = read_csv(csvfile);

    if (data.empty()) {
        cerr << "Nessun dato letto. Interrompo." << endl;
        return;
    }

    // Ordine fisso delle fasi sull'asse x.
    vector<string> phases = {
        "before_annealing",
        "annealing_T=75_h=5",
        "annealing_T=75_h=25",
        "annealing_T=100_h=5"
    };

    map<string, int> phase_index;

    for (int i = 0; i < (int)phases.size(); i++) {
        phase_index[phases[i]] = i;
    }

    // Organizzazione dati:
    // spot -> phase -> Row
    map<int, map<string, Row>> by_spot;

    for (const auto& r : data) {
        if (phase_index.find(r.phase) == phase_index.end()) {
            cerr << "Warning: fase non riconosciuta: "
                 << r.phase << endl;
            continue;
        }

        // Se use_only_detected = true, i punti non rilevati vengono esclusi.
        // Di default e' false, cosi' i punti con luminosity = 0 e detected = False
        // restano visibili nei grafici e nel rapporto.
        if (use_only_detected && !r.detected) continue;

        by_spot[r.spot][r.phase] = r;
    }

    if (by_spot.empty()) {
        cerr << "Nessuno spot disponibile dopo il filtraggio." << endl;
        return;
    }

    // --------------------------------------------------------
    // Raggruppamento degli spot per ordine di grandezza
    // della luminosita' massima dello spot
    // --------------------------------------------------------

    map<int, vector<int>> groups_by_order;

    for (const auto& entry : by_spot) {
        int spot = entry.first;
        const auto& phase_map = entry.second;

        double max_lum = -1.0;

        for (const auto& ph : phases) {
            if (phase_map.find(ph) == phase_map.end()) continue;

            max_lum = max(max_lum, phase_map.at(ph).luminosity);
        }

        if (max_lum <= 0) {
            groups_by_order[-999].push_back(spot);
        } else {
            int order = (int)floor(log10(max_lum));
            groups_by_order[order].push_back(spot);
        }
    }

    // Colori e marker per gli spot.
    vector<int> colors = {
        kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1,
        kOrange+7, kCyan+2, kViolet+1, kAzure+1, kPink+7,
        kSpring+5, kTeal+3, kGray+2, kRed-4, kBlue-4
    };

    vector<int> markers = {
        20, 21, 22, 23, 24,
        25, 26, 27, 28, 29,
        30, 33, 34, 43, 47
    };

    // --------------------------------------------------------
    // Canvas luminosity vs phase
    // --------------------------------------------------------

    int canvas_counter = 0;

    for (auto& group_entry : groups_by_order) {

        int order = group_entry.first;
        vector<int>& spots = group_entry.second;

        sort(spots.begin(), spots.end());

        // Se nello stesso ordine di grandezza ci sono piu' di
        // MAX_SPOTS_PER_CANVAS spot, li divido in sottogruppi.
        for (int start = 0; start < (int)spots.size(); start += MAX_SPOTS_PER_CANVAS) {

            int end = min(start + MAX_SPOTS_PER_CANVAS, (int)spots.size());

            vector<int> subspots(
                spots.begin() + start,
                spots.begin() + end
            );

            double ymin = 1e99;
            double ymax = -1e99;

            for (int spot : subspots) {
                for (const auto& ph : phases) {
                    if (by_spot[spot].find(ph) == by_spot[spot].end()) continue;

                    const Row& r = by_spot[spot][ph];

                    ymin = min(ymin, r.luminosity - r.error);
                    ymax = max(ymax, r.luminosity + r.error);
                }
            }

            if (ymin == 1e99 || ymax == -1e99) continue;

            if (ymin > 0) ymin *= 0.8;
            else          ymin *= 1.2;

            ymax *= 1.25;

            if (ymax <= ymin) ymax = ymin + 1.0;

            string cname = Form("c_phase_group_%d", canvas_counter);

            string ctitle;

            if (order == -999) {
                ctitle = "Luminosity vs phase - zero or negative luminosity";
            } else {
                ctitle = Form("Luminosity vs phase - order 10^{%d}", order);
            }

            TCanvas* c = new TCanvas(
                cname.c_str(),
                ctitle.c_str(),
                1100,
                750
            );

            c->SetGrid();

            TH1D* frame = new TH1D(
                Form("frame_%d", canvas_counter),
                ctitle.c_str(),
                (int)phases.size(),
                -0.5,
                (double)phases.size() - 0.5
            );

            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);

            frame->GetXaxis()->SetTitle("phase");
            frame->GetYaxis()->SetTitle("luminosity");

            for (int i = 0; i < (int)phases.size(); i++) {
                frame->GetXaxis()->SetBinLabel(i + 1, phases[i].c_str());
            }

            frame->GetXaxis()->LabelsOption("v");
            frame->GetXaxis()->SetLabelSize(0.035);
            frame->GetYaxis()->SetTitleOffset(1.25);

            frame->Draw();

            TLegend* leg = new TLegend(0.72, 0.58, 0.90, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            int ig = 0;

            for (int spot : subspots) {

                vector<double> xval;
                vector<double> yval;
                vector<double> exval;
                vector<double> eyval;

                for (int iph = 0; iph < (int)phases.size(); iph++) {
                    const string& ph = phases[iph];

                    if (by_spot[spot].find(ph) == by_spot[spot].end()) continue;

                    const Row& r = by_spot[spot][ph];

                    xval.push_back((double)iph);
                    yval.push_back(r.luminosity);
                    exval.push_back(0.0);
                    eyval.push_back(r.error);
                }

                if (xval.empty()) continue;

                TGraphErrors* gr = new TGraphErrors(
                    (int)xval.size(),
                    xval.data(),
                    yval.data(),
                    exval.data(),
                    eyval.data()
                );

                int col = colors[ig % colors.size()];
                int mar = markers[ig % markers.size()];

                gr->SetLineColor(col);
                gr->SetMarkerColor(col);
                gr->SetMarkerStyle(mar);
                gr->SetMarkerSize(1.1);
                gr->SetLineWidth(2);

                gr->Draw("PL SAME");

                leg->AddEntry(gr, Form("spot %d", spot), "pl");

                ig++;
            }

            leg->Draw();

            c->Update();

            if (save_png) {
                string outname;

                if (order == -999) {
                    outname = Form(
                        "luminosity_vs_phase_group_%02d_zero_or_negative.png",
                        canvas_counter
                    );
                } else {
                    outname = Form(
                        "luminosity_vs_phase_group_%02d_order_%d.png",
                        canvas_counter,
                        order
                    );
                }

                c->SaveAs(outname.c_str());
            }

            canvas_counter++;
        }
    }

    cout << "Create " << canvas_counter
         << " canvas per luminosity vs phase." << endl;


    // --------------------------------------------------------
    // Canvas finali rapporto:
    // L(annealing_T=100_h=5) / L(before_annealing)
    // --------------------------------------------------------

    string phase_before = "before_annealing";
    string phase_num    = "annealing_T=100_h=5";

    vector<double> xratio;
    vector<double> yratio;
    vector<double> exratio;
    vector<double> eyratio;

    for (const auto& entry : by_spot) {

        int spot = entry.first;
        const auto& phase_map = entry.second;

        if (phase_map.find(phase_before) == phase_map.end()) {
            cerr << "Warning: spot " << spot
                 << " non ha la fase " << phase_before << endl;
            continue;
        }

        if (phase_map.find(phase_num) == phase_map.end()) {
            cerr << "Warning: spot " << spot
                 << " non ha la fase " << phase_num << endl;
            continue;
        }

        const Row& r_before = phase_map.at(phase_before);
        const Row& r_num    = phase_map.at(phase_num);

        double L_before = r_before.luminosity;
        double L_num    = r_num.luminosity;

        double e_before = r_before.error;
        double e_num    = r_num.error;

        if (L_before == 0) {
            cerr << "Warning: spot " << spot
                 << " ha luminosity before_annealing = 0. Rapporto saltato."
                 << endl;
            continue;
        }

        // Rapporto richiesto:
        // R = L(annealing_T=100_h=5) / L(before_annealing)
        double ratio = L_num / L_before;

        // Propagazione robusta dell'errore, valida anche se L_num = 0:
        // dR/dL_num    = 1 / L_before
        // dR/dL_before = -L_num / L_before^2
        double eratio = sqrt(
            pow(e_num / L_before, 2) +
            pow((L_num * e_before) / (L_before * L_before), 2)
        );

        xratio.push_back((double)spot);
        yratio.push_back(ratio);
        exratio.push_back(0.0);
        eyratio.push_back(eratio);
    }

    if (!xratio.empty()) {

        // Prima canvas finale: rapporto con range y automatico.
        draw_ratio_canvas(
            xratio,
            yratio,
            exratio,
            eyratio,
            "c_ratio_100_over_before",
            "Ratio annealing_T=100_h=5 / before_annealing",
            "ratio_annealing_T100_h5_over_before.png",
            save_png,
            false
        );

        // Seconda canvas finale: stesso rapporto, focus intorno a 1.
        draw_ratio_canvas(
            xratio,
            yratio,
            exratio,
            eyratio,
            "c_ratio_100_over_before_focus",
            "Ratio annealing_T=100_h=5 / before_annealing - focus",
            "ratio_annealing_T100_h5_over_before_focus.png",
            save_png,
            true,
            RATIO_FOCUS_YMIN,
            RATIO_FOCUS_YMAX
        );

        // ----------------------------------------------------
        // Significance:
        // numero di punti sopra 1 di almeno 5 sigma
        // ----------------------------------------------------

        int n_below_5sigma = 0;

        cout << endl;
        cout << "Points below 1 by at least "
             << SIGMA_THRESHOLD << " sigma:" << endl;

        for (int i = 0; i < (int)yratio.size(); i++) {

            double x  = xratio[i];
            double y  = yratio[i];
            double ey = eyratio[i];

            if (ey <= 0) {
                cout << "  spot " << x
                     << "  ratio = " << y
                     << " +/- " << ey
                     << "  significance non calcolabile: errore nullo o negativo"
                     << endl;
                continue;
            }

            double significance = (1.0 - y) / ey;

            if (significance >= SIGMA_THRESHOLD) {
                n_below_5sigma++;

                cout << "  spot " << x
                     << "  ratio = " << y
                     << " +/- " << ey
                     << "  significance = " << significance
                     << " sigma"
                     << endl;
            }
        }

        cout << "Number of points below 1 by at least "
             << SIGMA_THRESHOLD << " sigma: "
             << n_below_5sigma << " / " << yratio.size()
             << endl;

    } else {
        cerr << "Nessun rapporto calcolabile." << endl;
    }
}
