// ============================================================
// 
// root -l 'lum_vs_T_irradiated.C("A1_v=5_annealed_100_5_run=20260612-001017_total_global_ID.csv")'

// 
//lum_vs_T_irradiated.C
//
// CSV atteso:
// spot,x,y,luminosity,error,T,v
//
// Fit:
// Lum(T) = A * exp(B*T)
//
// Uso:
// root -l
// .x lum_vs_T_irradiated.C("nome_file.csv")
//
// Oppure, scegliendo lo spot della canvas 4:
// .x lum_vs_T_irradiated.C("nome_file.csv", 3)
// ============================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TObject.h"

using namespace std;

struct Row {
    int spot;
    double x, y;
    double luminosity, error;
    double T;
    double v;
};

// ============================================================
// Lettura CSV
// Formato:
// spot,x,y,luminosity,error,T,v
// ============================================================

vector<Row> read_csv(const char* filename) {

    vector<Row> data;
    ifstream fin(filename);

    if (!fin.is_open()) {
        cerr << "Errore: impossibile aprire " << filename << endl;
        return data;
    }

    string line;
    getline(fin, line); // header

    while (getline(fin, line)) {

        if (line.empty()) continue;

        stringstream ss(line);
        string field;
        vector<string> f;

        while (getline(ss, field, ',')) {
            f.push_back(field);
        }

        if (f.size() < 7) continue;

        Row r;

        r.spot       = stoi(f[0]);
        r.x          = stod(f[1]);
        r.y          = stod(f[2]);
        r.luminosity = stod(f[3]);
        r.error      = stod(f[4]);
        r.T          = stod(f[5]);
        r.v          = stod(f[6]);

        data.push_back(r);
    }

    return data;
}

// ============================================================
// Stima iniziale parametri esponenziale
// Lum = A * exp(B*T)
// ============================================================

void estimate_exp_parameters(
    const vector<Row>& rows,
    double& A0,
    double& B0
) {
    vector<Row> positive_rows;

    for (const auto& r : rows) {
        if (r.luminosity > 0) {
            positive_rows.push_back(r);
        }
    }

    sort(
        positive_rows.begin(),
        positive_rows.end(),
        [](const Row& a, const Row& b) {
            return a.T < b.T;
        }
    );

    if (positive_rows.size() >= 2) {

        const Row& r1 = positive_rows.front();
        const Row& r2 = positive_rows.back();

        if (fabs(r2.T - r1.T) > 1e-9) {
            B0 = (log(r2.luminosity) - log(r1.luminosity)) / (r2.T - r1.T);
            A0 = exp(log(r1.luminosity) - B0 * r1.T);
        } else {
            A0 = r1.luminosity;
            B0 = 0.0;
        }

    } else if (positive_rows.size() == 1) {

        A0 = positive_rows[0].luminosity;
        B0 = 0.0;

    } else {

        A0 = 1.0;
        B0 = 0.0;
    }
}

// ============================================================
// Macro principale
// ============================================================

void lum_vs_T_irradiated(const char* filename = "data.csv", int selected_spot = 3) {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);

    vector<Row> data = read_csv(filename);

    if (data.empty()) {
        cerr << "Errore: nessun dato letto dal file " << filename << endl;
        return;
    }

    map<int, vector<Row>> spots;

    for (auto& r : data) {
        spots[r.spot].push_back(r);
    }

    int nspots = spots.size();

    double v_const = data[0].v;

    double Tmin = 1e9;
    double Tmax = -1e9;

    double ymin_global = 1e9;
    double ymax_global = -1e9;

    int spot_min = 1e9;
    int spot_max = -1e9;

    for (const auto& entry : spots) {

        int spot = entry.first;
        const vector<Row>& rows = entry.second;

        if (spot < spot_min) spot_min = spot;
        if (spot > spot_max) spot_max = spot;

        for (const auto& r : rows) {

            if (r.T < Tmin) Tmin = r.T;
            if (r.T > Tmax) Tmax = r.T;

            if (r.luminosity - r.error < ymin_global)
                ymin_global = r.luminosity - r.error;

            if (r.luminosity + r.error > ymax_global)
                ymax_global = r.luminosity + r.error;
        }
    }

    double dT = Tmax - Tmin;
    if (dT <= 0) dT = 1.0;

    double T_fit_min = Tmin - 0.10 * dT;
    double T_fit_max = Tmax + 0.10 * dT;

    vector<int> colors = {
        kRed+1, kBlue+1, kGreen+2, kMagenta+1,
        kOrange+7, kCyan+1, kViolet, kSpring+5,
        kPink+7, kAzure+7, kTeal+3, kGray+2,
        kBlack, kYellow+2
    };

    // ============================================================
    // Vettori per salvare i parametri del fit esponenziale
    // ============================================================

    vector<double> x_spot_fit;
    vector<double> ex_spot_fit;

    vector<double> B_vals;
    vector<double> B_errs;

    vector<double> chi2_vals;
    vector<int> ndf_vals;

    // ============================================================
    // FIT ESPONENZIALE PER TUTTI GLI SPOT
    // Serve per Canvas 5 e Canvas 6
    // ============================================================

    for (auto& entry : spots) {

        int spot = entry.first;
        vector<Row> rows = entry.second;

        sort(
            rows.begin(),
            rows.end(),
            [](const Row& a, const Row& b) {
                return a.T < b.T;
            }
        );

        int n = rows.size();

        vector<double> xT(n), yL(n), exT(n), eyL(n);

        for (int i = 0; i < n; i++) {
            xT[i]  = rows[i].T;
            yL[i]  = rows[i].luminosity;
            exT[i] = 0.0;
            eyL[i] = rows[i].error;
        }

        TGraphErrors* gr_tmp = new TGraphErrors(
            n,
            xT.data(),
            yL.data(),
            exT.data(),
            eyL.data()
        );

        double A0, B0;
        estimate_exp_parameters(rows, A0, B0);

        TF1* fit_exp = new TF1(
            Form("fit_exp_global_spot_%d", spot),
            "[0]*exp([1]*x)",
            T_fit_min,
            T_fit_max
        );

        fit_exp->SetParNames("A", "B");
        fit_exp->SetParameters(A0, B0);

        double B = 0.0;
        double eB = 0.0;
        double chi2 = 0.0;
        int ndf = 0;

        if (n >= 3) {
            gr_tmp->Fit(fit_exp, "RQ0");

            B = fit_exp->GetParameter(1);
            eB = fit_exp->GetParError(1);
            chi2 = fit_exp->GetChisquare();
            ndf  = fit_exp->GetNDF();
        }

        x_spot_fit.push_back(spot);
        ex_spot_fit.push_back(0.0);

        B_vals.push_back(B);
        B_errs.push_back(eB);

        chi2_vals.push_back(chi2);
        ndf_vals.push_back(ndf);
    }

    // ============================================================
    // CANVAS 1:
    // luminosity vs spot at max T
    // ============================================================

    vector<double> x_spot_c1;
    vector<double> y_lum_c1;
    vector<double> ex_spot_c1;
    vector<double> ey_lum_c1;

    for (auto& entry : spots) {

        int spot = entry.first;
        const vector<Row>& rows = entry.second;

        bool found = false;

        for (const auto& r : rows) {

            if (fabs(r.T - Tmax) < 1e-6) {

                x_spot_c1.push_back(spot);
                y_lum_c1.push_back(r.luminosity);
                ex_spot_c1.push_back(0.0);
                ey_lum_c1.push_back(r.error);

                found = true;
                break;
            }
        }

        if (!found) {
            cout << "Spot " << spot
                 << " non contiene T = " << Tmax << endl;
        }
    }

    TCanvas* c1 = new TCanvas(
        "c1_lum_vs_spot_maxT",
        "Luminosity vs spot at max T",
        1800,
        600
    );

    TGraphErrors* gr_c1 = new TGraphErrors(
        x_spot_c1.size(),
        x_spot_c1.data(),
        y_lum_c1.data(),
        ex_spot_c1.data(),
        ey_lum_c1.data()
    );

    gr_c1->SetTitle(
        Form("Luminosity at maximum temperature T = %.1f #circC, v = %.1f - Annealing 100#circC 5h;Spot;Luminosity",
             Tmax, v_const)
    );

    gr_c1->SetMarkerStyle(20);
    gr_c1->SetMarkerSize(1.2);
    gr_c1->SetLineWidth(2);

    gr_c1->Draw("AP");

    c1->SaveAs("A1_ann_100_5_lum_vs_T_maxT.png");

    // ============================================================
    // CANVAS 2:
    // luminosity vs spot at max T - focus y range (-5, 1000)
    // ============================================================

    TCanvas* c2 = new TCanvas(
        "c2_lum_vs_spot_maxT_focus",
        "Luminosity vs spot at max T - focus",
        1800,
        600
    );

    TGraphErrors* gr_c2 = new TGraphErrors(
        x_spot_c1.size(),
        x_spot_c1.data(),
        y_lum_c1.data(),
        ex_spot_c1.data(),
        ey_lum_c1.data()
    );

    gr_c2->SetTitle(
        Form("Luminosity at maximum temperature T = %.1f #circC, v = %.1f Annealing 100#circC 5h - focus;Spot;Luminosity",
             Tmax, v_const)
    );

    gr_c2->SetMarkerStyle(20);
    gr_c2->SetMarkerSize(1.2);
    gr_c2->SetLineWidth(2);

    gr_c2->SetMinimum(-5);
    gr_c2->SetMaximum(1000);

    gr_c2->Draw("AP");

    c2->SaveAs("A1_ann_100_5_lum_vs_T_maxT_focus.png");

    // ============================================================
    // CANVAS 3:
    // luminosity vs T per i primi 8 spot, 8 pad
    // con fit Lum = A exp(BT)
    // ============================================================

    TCanvas* c3 = new TCanvas(
        "c3_lum_vs_T_first8",
        "Luminosity vs T - first 8 spots",
        1600,
        900
    );

    c3->Divide(4, 2);

    int displayed = 0;

    for (auto& entry : spots) {

        if (displayed >= 8) break;

        int spot = entry.first;
        vector<Row> rows = entry.second;

        sort(
            rows.begin(),
            rows.end(),
            [](const Row& a, const Row& b) {
                return a.T < b.T;
            }
        );

        int n = rows.size();

        vector<double> xT(n), yL(n), exT(n), eyL(n);

        double ymin_local = 1e9;
        double ymax_local = -1e9;

        for (int i = 0; i < n; i++) {
            xT[i]  = rows[i].T;
            yL[i]  = rows[i].luminosity;
            exT[i] = 0.0;
            eyL[i] = rows[i].error;

            if (rows[i].luminosity - rows[i].error < ymin_local)
                ymin_local = rows[i].luminosity - rows[i].error;

            if (rows[i].luminosity + rows[i].error > ymax_local)
                ymax_local = rows[i].luminosity + rows[i].error;
        }

        if (ymin_local > 0)
            ymin_local *= 0.80;
        else
            ymin_local = 0.0;

        ymax_local *= 1.25;

        c3->cd(displayed + 1);

        TGraphErrors* gr = new TGraphErrors(
            n,
            xT.data(),
            yL.data(),
            exT.data(),
            eyL.data()
        );

        gr->SetTitle(
            Form("Spot %d: x = %.2f, y = %.2f;T (#circC) B.A;Luminosity",
                 spot, rows[0].x, rows[0].y)
        );

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetLineWidth(2);

        gr->SetMinimum(ymin_local);
        gr->SetMaximum(ymax_local);

        gr->Draw("AP");

        double A0, B0;
        estimate_exp_parameters(rows, A0, B0);

        TF1* fit_exp = new TF1(
            Form("fit_exp_canvas3_spot_%d", spot),
            "[0]*exp([1]*x)",
            T_fit_min,
            T_fit_max
        );

        fit_exp->SetParNames("A", "B");
        fit_exp->SetParameters(A0, B0);
        fit_exp->SetLineColor(kRed + 1);
        fit_exp->SetLineWidth(2);

        if (n >= 3) {
            gr->Fit(fit_exp, "RQ");
        }

        displayed++;
    }

    c3->SaveAs("A1_ann_100_5_lum_vs_T_first8_expfit.png");

    // ============================================================
    // CANVAS 4:
    // luminosity vs T for a single spot with fit
    // ============================================================

    if (spots.find(selected_spot) != spots.end()) {

        TCanvas* c4 = new TCanvas(
            "c4_single_spot_expfit",
            Form("Spot %d - Luminosity vs T", selected_spot),
            900,
            700
        );

        vector<Row> rows = spots[selected_spot];

        sort(
            rows.begin(),
            rows.end(),
            [](const Row& a, const Row& b) {
                return a.T < b.T;
            }
        );

        int n = rows.size();

        vector<double> xT(n), yL(n), exT(n), eyL(n);

        double ymin_local = 1e9;
        double ymax_local = -1e9;

        for (int i = 0; i < n; i++) {
            xT[i]  = rows[i].T;
            yL[i]  = rows[i].luminosity;
            exT[i] = 0.0;
            eyL[i] = rows[i].error;

            if (rows[i].luminosity - rows[i].error < ymin_local)
                ymin_local = rows[i].luminosity - rows[i].error;

            if (rows[i].luminosity + rows[i].error > ymax_local)
                ymax_local = rows[i].luminosity + rows[i].error;
        }

        if (ymin_local > 0)
            ymin_local *= 0.80;
        else
            ymin_local = 0.0;

        ymax_local *= 1.25;

        TGraphErrors* gr = new TGraphErrors(
            n,
            xT.data(),
            yL.data(),
            exT.data(),
            eyL.data()
        );

        gr->SetTitle(
            Form("Spot %d: x = %.2f, y = %.2f, v = %.1f;T (#circC) Annealing 100#circC 5h;Luminosity",
                 selected_spot, rows[0].x, rows[0].y, v_const)
        );

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.4);
        gr->SetLineWidth(2);

        gr->SetMinimum(ymin_local);
        gr->SetMaximum(ymax_local);

        gr->Draw("AP");

        double A0, B0;
        estimate_exp_parameters(rows, A0, B0);

        TF1* fit_exp = new TF1(
            Form("fit_exp_canvas4_spot_%d", selected_spot),
            "[0]*exp([1]*x)",
            T_fit_min,
            T_fit_max
        );

        fit_exp->SetParNames("A", "B");
        fit_exp->SetParameters(A0, B0);
        fit_exp->SetLineColor(kRed + 1);
        fit_exp->SetLineWidth(2);

        double A = 0.0;
        double eA = 0.0;
        double B = 0.0;
        double eB = 0.0;
        double chi2 = 0.0;
        int ndf = 0;

        if (n >= 3) {
            gr->Fit(fit_exp, "RQ");

            A = fit_exp->GetParameter(0);
            eA = fit_exp->GetParError(0);

            B = fit_exp->GetParameter(1);
            eB = fit_exp->GetParError(1);

            chi2 = fit_exp->GetChisquare();
            ndf  = fit_exp->GetNDF();
        }

        TLegend* leg4 = new TLegend(0.14, 0.68, 0.60, 0.88);
        leg4->SetBorderSize(0);
        leg4->SetFillStyle(0);

        leg4->AddEntry(gr, "Data", "lep");

        if (n >= 3) {
            leg4->AddEntry(fit_exp, "Fit: Lum = A e^{BT}", "l");
            leg4->AddEntry((TObject*)0, Form("A = %.3g #pm %.2g", A, eA), "");
            leg4->AddEntry((TObject*)0, Form("B = %.3g #pm %.2g", B, eB), "");
            leg4->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf), "");
        }

        leg4->Draw();

        c4->SaveAs(Form("A1_ann_100_5_lum_vs_T_spot%d.png", selected_spot));

    } else {

        cout << "Spot " << selected_spot << " non presente nel file." << endl;
    }

    // ============================================================
    // CANVAS 5:
    // B vs spot ID, con fit costante B = const
    // Usa solo gli spot con chi2/ndf < 2 nel fit esponenziale
    // ============================================================

    vector<double> x_spot_fit_good;
    vector<double> ex_spot_fit_good;
    vector<double> B_vals_good;
    vector<double> B_errs_good;

    for (size_t i = 0; i < B_vals.size(); i++) {

        if (ndf_vals[i] <= 0) continue;

        double chi2_ndf = chi2_vals[i] / ndf_vals[i];

        if (chi2_ndf < 2.0) {
            x_spot_fit_good.push_back(x_spot_fit[i]);
            ex_spot_fit_good.push_back(ex_spot_fit[i]);
            B_vals_good.push_back(B_vals[i]);
            B_errs_good.push_back(B_errs[i]);
        }
    }

    TCanvas* c5 = new TCanvas("c5_B_vs_spot_constfit", "B vs spot ID",1800,600);

    TGraphErrors* gr_B = new TGraphErrors(
        x_spot_fit_good.size(),
        x_spot_fit_good.data(),
        B_vals_good.data(),
        ex_spot_fit_good.data(),
        B_errs_good.data()
    );

    gr_B->SetTitle("Exponential fit parameter B vs spot ID annealing 100#circC 5h (#chi^{2}/ndf < 2);Spot;B");
    gr_B->SetMarkerStyle(21);
    gr_B->SetMarkerSize(1.0);
    gr_B->SetLineWidth(2);

    gr_B->Draw("AP");
    gr_B->GetXaxis()->SetLimits(spot_min - 1, spot_max + 1);

    double B_min = *min_element(B_vals_good.begin(), B_vals_good.end());
    double B_max = *max_element(B_vals_good.begin(), B_vals_good.end());

    double margin_low  = 0.25 * (B_max - B_min);
    double margin_high = 0.90 * (B_max - B_min);  // più spazio sopra per la legenda

    gr_B->GetYaxis()->SetRangeUser(B_min - margin_low, B_max + margin_high);

    // Fit costante: B = const
    TF1* fit_B_const = new TF1( "fit_B_const","[0]", spot_min - 1,spot_max + 1);

    fit_B_const->SetParNames("B_{const}");
    fit_B_const->SetLineColor(kRed + 1);
    fit_B_const->SetLineWidth(2);

    if (gr_B->GetN() >= 2) {
        gr_B->Fit(fit_B_const, "RQ");
    }

    TLegend* leg5 = new TLegend(0.1, 0.72, 0.52, 0.88);
    leg5->SetBorderSize(0);
    leg5->SetFillStyle(0);

    leg5->AddEntry(gr_B, "B from Lum = A e^{BT}", "lep");

    if (gr_B->GetN() >= 2) {

        double lambda = fit_B_const->GetParameter(0);
        double e_lambda = fit_B_const->GetParError(0);

        // Radioactive-decay convention: N(t) = N0 exp(-lambda*t)
        // Half-time: T_1/2 = ln(2)/lambda
        double half_time = log(2.0) / lambda;
        double e_half_time = log(2.0) * e_lambda / (lambda * lambda);

        leg5->AddEntry(fit_B_const, "Constant fit: B = const", "l");
 
        leg5->AddEntry(
            (TObject*)0,
            Form("T_{1/2} = %.3g #pm %.2g",
                half_time,
                e_half_time),
            ""
        );
    }

    leg5->Draw();

    c5->SaveAs("A1_ann_100_5_lum_vs_T_B.png");
    // ============================================================
    // CANVAS 6:
    // chi square del fit esponenziale vs spot ID
    // ============================================================

    TCanvas* c6 = new TCanvas(
        "c6_chi2_expfit_vs_spot",
        "Chi square exponential fit vs spot ID Annealing 100#circC 5h",
        1000,
        800
    );

    TGraph* gr_chi2 = new TGraph(
        x_spot_fit.size(),
        x_spot_fit.data(),
        chi2_vals.data()
    );

    gr_chi2->SetTitle("#chi^{2} of exponential fit vs spot ID Annealing 100#circC 5h;Spot;#chi^{2}");
    gr_chi2->SetMarkerStyle(22);
    gr_chi2->SetMarkerSize(1.2);
    gr_chi2->SetLineWidth(2);

    gr_chi2->Draw("AP");
    gr_chi2->GetXaxis()->SetLimits(spot_min - 1, spot_max + 1);

    c6->SaveAs("A1_ann_100_5_lum_vs_T_chi2.png");
}