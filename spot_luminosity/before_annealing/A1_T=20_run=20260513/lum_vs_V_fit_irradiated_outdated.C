
//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_run=20260513-032137_total_quadratic.csv")'
//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_run=20260513-032137_total_rho=0.0.csv")'
//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_run=20260513-032137_total_rho=0.1.csv")'
//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_run=20260513-032137_total_rho=0.2.csv")'
//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_run=20260513-032137_total_linear_error.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TObject.h"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TPaveText.h"
#include "TSystem.h"

using namespace std;

struct Row {
    int spot;
    double x, y;
    double luminosity, error;
    double T;
    double v;       // tensione assoluta
    double v_fin;   // overvoltage finale
};

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

        if (f.size() < 8) continue;

        Row r;
        r.spot       = stoi(f[0]);
        r.x          = stod(f[1]);
        r.y          = stod(f[2]);
        r.luminosity = stod(f[3]);
        r.error      = stod(f[4]);
        r.T          = stod(f[5]);
        r.v          = stod(f[6]);
        r.v_fin      = stod(f[7]);

        data.push_back(r);
    }

    return data;
}

void lum_vs_V_fit_irradiated_outdated(const char* filename = "data.csv") {

    gStyle->SetOptFit(111);

    vector<Row> data = read_csv(filename);

    map<int, vector<Row>> spots;

    for (auto &r : data) {
        spots[r.spot].push_back(r);
    }

    int nspots = spots.size();
    
    //Parametri display canvas
    int max_spot_display = 8;
    int ncols = 4;
    int nrows = 2;

    //Variables for the fit
    vector<int> spot_ids;
    vector<double> A_vals, A_errs;
    vector<double> B_vals, B_errs;
    vector<double> chi2ndf_vals;

    // ============================================================
    // CANVAS 5: luminosity vs spot at max v_fin
    // ============================================================

    // trova il massimo valore di v_fin presente nel file
    double vmax = -1e9;

    for (const auto& r : data) {
        if (r.v_fin > vmax)
            vmax = r.v_fin;
    }

    vector<double> x_spot;
    vector<double> y_lum;
    vector<double> ey_lum;

    for (auto &entry : spots) {

        int spot = entry.first;
        const vector<Row>& rows = entry.second;

        bool found = false;

        for (const auto& r : rows) {

            // tolleranza per confronti floating point
            if (fabs(r.v_fin - vmax) < 1e-6) {

                x_spot.push_back(spot);
                y_lum.push_back(r.luminosity);
                ey_lum.push_back(r.error);

                found = true;
                break;
            }
        }

        if (!found) {
            cout << "Spot " << spot
                << " non contiene v_fin = "
                << vmax << endl;
        }
    }

    int nmax = x_spot.size();

    vector<double> ex_spot(nmax, 0.0);

    TCanvas* c5 = new TCanvas("c5",
                            Form("Luminosity at v_fin = %.2f V", vmax),
                            1800,
                            600);

    TGraphErrors* grMax =
        new TGraphErrors(nmax,
                        x_spot.data(),
                        y_lum.data(),
                        ex_spot.data(),
                        ey_lum.data());

    grMax->SetTitle(
        Form("Luminosity at maximum overvoltage (v_{fin}=%.2f V) before annealing;Spot;Luminosity",
            vmax));

    grMax->SetMarkerStyle(20);
    grMax->SetMarkerSize(1.2);
    grMax->SetLineWidth(2);

    grMax->Draw("AP");

    c5->SaveAs("A1_bef_ann_lum_vs_v_max_vfin.png");

    // ============================================================
    // CANVAS 5.1: luminosity vs spot at max v_fin - focus y range
    // ============================================================

    // trova il massimo valore di v_fin presente nel file
    double vmax_c51 = -1e9;

    for (const auto& r : data) {
        if (r.v_fin > vmax_c51)
            vmax_c51 = r.v_fin;
    }

    vector<double> x_spot_c51;
    vector<double> y_lum_c51;
    vector<double> ey_lum_c51;

    for (auto &entry_c51 : spots) {

        int spot_c51 = entry_c51.first;
        const vector<Row>& rows_c51 = entry_c51.second;

        bool found_c51 = false;

        for (const auto& r_c51 : rows_c51) {

            // tolleranza per confronti floating point
            if (fabs(r_c51.v_fin - vmax_c51) < 1e-6) {

                x_spot_c51.push_back(spot_c51);
                y_lum_c51.push_back(r_c51.luminosity);
                ey_lum_c51.push_back(r_c51.error);

                found_c51 = true;
                break;
            }
        }

        if (!found_c51) {
            cout << "Spot " << spot_c51
                << " non contiene v_fin = "
                << vmax_c51 << endl;
        }
    }

    int nmax_c51 = x_spot_c51.size();

    vector<double> ex_spot_c51(nmax_c51, 0.0);

    TCanvas* c51 = new TCanvas("c51",
                            Form("Luminosity at v_fin = %.2f V - focus", vmax_c51),
                            1800,
                            600);

    TGraphErrors* grMax_c51 =
        new TGraphErrors(nmax_c51,
                        x_spot_c51.data(),
                        y_lum_c51.data(),
                        ex_spot_c51.data(),
                        ey_lum_c51.data());

    grMax_c51->SetTitle(
        Form("Luminosity at maximum overvoltage (v_{fin}=%.2f V) before annealing - focus;Spot;Luminosity",
            vmax_c51));

    grMax_c51->SetMarkerStyle(20);
    grMax_c51->SetMarkerSize(1.2);
    grMax_c51->SetLineWidth(2);

    // Focus sull'asse y
    grMax_c51->SetMinimum(-5);
    grMax_c51->SetMaximum(1000);

    grMax_c51->Draw("AP");

    c51->SaveAs("A1_bef_ann_lum_vs_v_max_vfin_focus.png");


    // ============================================================
    // CANVAS 6: luminosity vs v_fin per i primi 8 spot
    // ============================================================

    TCanvas* c6 = new TCanvas("c6", "Luminosity vs overvoltage", 1600, 1000);
    c6->Divide(4, 2);

    int displayed = 0;

    for (auto &entry : spots) {

    int spot = entry.first;
    vector<Row> rows = entry.second;

    sort(rows.begin(), rows.end(),
         [](const Row& a, const Row& b) {
             return a.v_fin < b.v_fin;
         });

    int n = rows.size();

    vector<double> xv(n), yv(n), exv(n), eyv(n);

    for (int i = 0; i < n; i++) {
        xv[i]  = rows[i].v_fin;
        yv[i]  = rows[i].luminosity;
        exv[i] = 0.0;
        eyv[i] = rows[i].error;
    }

    TGraphErrors* gr = new TGraphErrors(n,xv.data(),yv.data(), exv.data(), eyv.data());

    gr->SetTitle(Form("Spot %d: x=%.2f y=%.2f",
                      spot,
                      rows[0].x,
                      rows[0].y));

    gr->GetXaxis()->SetTitle("Overvoltage (V)");
    gr->GetYaxis()->SetTitle("Luminosity");

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);

    TF1* f_quad = new TF1(Form("f_quad_%d", spot), "[0]*TMath::Power(x,[1])", 0, 8);

    f_quad->SetParameters(1.0, 2.0);
    //f_quad->SetParLimits(0, 0, 100);
    f_quad->SetParLimits(1, 1.5, 2.5);

    gr->Fit(f_quad, "Q");

    double A       = f_quad->GetParameter(0);
    double A_err   = f_quad->GetParError(0);

    double B       = f_quad->GetParameter(1);
    double B_err   = f_quad->GetParError(1);

    double chi2    = f_quad->GetChisquare();
    double ndf     = f_quad->GetNDF();

    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    // salva SEMPRE i risultati
    spot_ids.push_back(spot);

    A_vals.push_back(A);
    A_errs.push_back(A_err);

    B_vals.push_back(B);
    B_errs.push_back(B_err);

    chi2ndf_vals.push_back(chi2ndf);

    // disegna solo i primi 8 spot
    if (displayed < 8) {
        c6->cd(displayed + 1);
        gr->Draw("AP");
        displayed++;
    }

    }
    c6->SaveAs("A1_bef_ann_lum_vs_v.png");
/*
    // ============================================================
    // CANVAS 7: luminosity vs V_over for a single spot with fit
    // ============================================================
    TCanvas* c7 = new TCanvas("c7", "Spot 2 - Luminosity vs overvoltage", 800, 600);

    // prendi solo spot 2
    int spot_2 = 2;
    vector<Row> rows_2 = spots[spot_2];

    // ordina per overvoltage
    sort(rows_2.begin(), rows_2.end(),
        [](const Row& a, const Row& b) {
            return a.v_fin < b.v_fin;
        });

    int n = rows_2.size();

    vector<double> x_2(n), y_2(n), ex_2(n), ey_2(n);

    for (int i = 0; i < n; i++) {
        x_2[i]  = rows_2[i].v_fin;
        y_2[i]  = rows_2[i].luminosity;
        ex_2[i] = 0.0;
        ey_2[i] = rows_2[i].error;
    }

    TGraphErrors* gr = new TGraphErrors(n, x_2.data(), y_2.data(), ex_2.data(), ey_2.data());

    gr->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f #circC before annealing",
                    spot_2, rows_2[0].x, rows_2[0].y, rows_2[0].T));

    gr->GetXaxis()->SetTitle("Overvoltage (V)");
    gr->GetXaxis()->SetTitleSize(0.045);
    //gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetTitle("Luminosity");
    gr->GetYaxis()->SetTitleSize(0.045);
    //gr->GetYaxis()->SetTitleOffset(1.3);

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(2.0);
    gr->SetLineWidth(2);
    
    // fit
    TF1* f_quad = new TF1("f_quad_2", "[0] * TMath::Power(x,[1])", 0, 10);
    f_quad->SetParameters(1.0, 2.0);
    gr->Fit(f_quad, "Q");


    gr->Draw("AP");
    c7->SaveAs("A1_bef_ann_lum_vs_v_spot2.png");

   

*/


    int n_fit = spot_ids.size();

    vector<double> x(n_fit), ex(n_fit, 0.0);
    for (int i = 0; i < n_fit; i++) x[i] = spot_ids[i];

    // ============================================================
    // CANVAS 7: luminosity vs V_over for all spots with fit
    // ============================================================

    // Cartella di output
    gSystem->mkdir("A1_bef_ann_lum_vs_v_all_spots", kTRUE);

    for (const auto& entry : spots) {

        int spot_id = entry.first;
        vector<Row> rows_spot = entry.second;

        if (rows_spot.size() < 2) {
            cout << "Skipping spot " << spot_id
                << " because it has fewer than 2 points." << endl;
            continue;
        }

        // ordina per overvoltage
        sort(rows_spot.begin(), rows_spot.end(),
            [](const Row& a, const Row& b) {
                return a.v_fin < b.v_fin;
            });

        int n_spot = rows_spot.size();

        vector<double> x_v(n_spot), y_lum(n_spot), ex_v(n_spot), ey_lum(n_spot);

        for (int i = 0; i < n_spot; i++) {
            x_v[i]    = rows_spot[i].v_fin;
            y_lum[i]  = rows_spot[i].luminosity;
            ex_v[i]   = 0.0;
            ey_lum[i] = rows_spot[i].error;
        }

        TCanvas* c7_spot = new TCanvas(
            Form("c7_spot_%d", spot_id),
            Form("Spot %d - Luminosity vs overvoltage", spot_id),
            800, 600
        );

        TGraphErrors* gr_spot = new TGraphErrors(
            n_spot,
            x_v.data(),
            y_lum.data(),
            ex_v.data(),
            ey_lum.data()
        );

        gr_spot->SetTitle(Form(
            "Spot %d: x = %.2f, y = %.2f, T = %.1f #circC before annealing",
            spot_id,
            rows_spot[0].x,
            rows_spot[0].y,
            rows_spot[0].T
        ));

        gr_spot->GetXaxis()->SetTitle("Overvoltage (V)");
        gr_spot->GetXaxis()->SetTitleSize(0.045);

        gr_spot->GetYaxis()->SetTitle("Luminosity");
        gr_spot->GetYaxis()->SetTitleSize(0.045);

        gr_spot->SetMarkerStyle(20);
        gr_spot->SetMarkerSize(2.0);
        gr_spot->SetLineWidth(2);

        // Fit: Lum = A * V_over^B
        TF1* f_quad_spot = new TF1(
            Form("f_quad_spot_%d", spot_id),
            "[0] * TMath::Power(x,[1])",
            0,
            10
        );

        f_quad_spot->SetParameters(1.0, 2.0);

        gr_spot->Draw("AP");

        // S = salva risultato fit
        // Q = quiet mode
        // R = usa il range definito nella TF1
        TFitResultPtr fit_result = gr_spot->Fit(f_quad_spot, "QRS");

        int fit_status = fit_result->Status();
        int cov_status = fit_result->CovMatrixStatus();
        double edm     = fit_result->Edm();

        double chi2 = f_quad_spot->GetChisquare();
        int ndf     = f_quad_spot->GetNDF();
        double prob = f_quad_spot->GetProb();

        // Disegna esplicitamente il fit
        f_quad_spot->Draw("SAME");

        // Box con informazioni sul fit
        TPaveText* fit_info = new TPaveText(0.12, 0.58, 0.45, 0.88, "NDC");
        fit_info->SetFillColor(0);
        fit_info->SetBorderSize(1);
        fit_info->SetTextAlign(12);
        fit_info->SetTextSize(0.030);

        fit_info->AddText(Form("A = %.3g #pm %.3g",
                            f_quad_spot->GetParameter(0),
                            f_quad_spot->GetParError(0)));

        fit_info->AddText(Form("B = %.3g #pm %.3g",
                            f_quad_spot->GetParameter(1),
                            f_quad_spot->GetParError(1)));

        if (ndf > 0) {
            fit_info->AddText(Form("#chi^{2}/ndf = %.2f/%d = %.2f",
                                chi2, ndf, chi2 / ndf));
        } else {
            fit_info->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf));
        }

        fit_info->AddText(Form("p-value = %.3g", prob));
        fit_info->AddText(Form("EDM = %.3g", edm));
        fit_info->AddText(Form("Fit status = %d", fit_status));
        fit_info->AddText(Form("CovMatrix status = %d", cov_status));

        fit_info->Draw("SAME");

        // Stampa anche su terminale
        cout << "----------------------------------------" << endl;
        cout << "Spot " << spot_id << endl;
        cout << "Fit status       = " << fit_status << endl;
        cout << "CovMatrix status = " << cov_status << endl;
        cout << "EDM              = " << edm << endl;
        cout << "Chi2             = " << chi2 << endl;
        cout << "NDF              = " << ndf << endl;
        if (ndf > 0) {
            cout << "Chi2/NDF         = " << chi2 / ndf << endl;
        }
        cout << "p-value          = " << prob << endl;
        cout << "A                = " << f_quad_spot->GetParameter(0)
            << " +/- " << f_quad_spot->GetParError(0) << endl;
        cout << "B                = " << f_quad_spot->GetParameter(1)
            << " +/- " << f_quad_spot->GetParError(1) << endl;

        c7_spot->SaveAs(Form(
            "A1_bef_ann_lum_vs_v_all_spots/A1_bef_ann_lum_vs_v_spot%d.png",
            spot_id
        ));

        delete c7_spot;
    }




    // ============================================================
    // CANVAS 9: B vs spot ID con selezione su chi2/ndf
    // ============================================================

    TCanvas* c9 = new TCanvas("c9", "B vs spot", 1800, 600);

    // ------------------------------------------------------------
    // Selezione dei soli fit "buoni": chi2/ndf < 5
    // ------------------------------------------------------------

    vector<double> x_B_good;
    vector<double> ex_B_good;
    vector<double> B_good;
    vector<double> B_err_good;

    for (int i = 0; i < n_fit; i++) {

        if (chi2ndf_vals[i] < 2.0) {

            x_B_good.push_back(x[i]);
            ex_B_good.push_back(0.0);

            B_good.push_back(B_vals[i]);
            B_err_good.push_back(B_errs[i]);
        }
    }

    int n_B_good = B_good.size();

    cout << "Canvas 9: numero di fit con chi2/ndf < 5 = "
        << n_B_good << " su " << n_fit << endl;


    // ------------------------------------------------------------
    // Controllo: se non ci sono abbastanza punti, evita crash
    // ------------------------------------------------------------

    if (n_B_good == 0) {

        cout << "ATTENZIONE: nessun punto con chi2/ndf < 5. "
            << "Canvas 9 non disegnata." << endl;

    } else {

        TGraphErrors* grB = new TGraphErrors(
            n_B_good,
            x_B_good.data(),
            B_good.data(),
            ex_B_good.data(),
            B_err_good.data()
        );

        grB->SetTitle("Power law exponent: B - before annealing;Spot;B");
        grB->SetMarkerStyle(21);
        grB->SetMarkerSize(1.2);
        grB->Draw("AP");

        grB->GetXaxis()->SetLimits(-1, nspots + 0.5);

        // ------------------------------------------------------------
        // Linea B = 2
        // ------------------------------------------------------------

        TF1* line_B = new TF1("line_B", "2", -1, nspots + 0.5);
        line_B->SetLineColor(kRed + 1);
        line_B->SetLineStyle(2);
        line_B->SetLineWidth(2);
        line_B->Draw("SAME");


        // ------------------------------------------------------------
        // Fit costante B = const
        // ------------------------------------------------------------

        TF1* fit_B = new TF1("fit_B", "[0]", -1, nspots + 0.5);
        fit_B->SetLineColor(kMagenta + 1);
        fit_B->SetLineWidth(2);

        if (n_B_good >= 2) {
            grB->Fit(fit_B, "RQ");
            fit_B->Draw("SAME");
        } else {
            cout << "ATTENZIONE: c'e' solo un punto buono. "
                << "Fit costante non eseguito." << endl;
        }


        // ------------------------------------------------------------
        // Risultati fit
        // ------------------------------------------------------------

        double chi2_B = 0.0;
        int ndf_B = 0;
        double fit_result_B = 0.0;
        double fit_error_B = 0.0;

        if (n_B_good >= 2) {
            chi2_B = fit_B->GetChisquare();
            ndf_B = fit_B->GetNDF();
            fit_result_B = fit_B->GetParameter(0);
            fit_error_B = fit_B->GetParError(0);
        }


        // ------------------------------------------------------------
        // Legenda a sinistra: oggetti grafici
        // ------------------------------------------------------------

        TLegend* leg9_left = new TLegend(0.12, 0.72, 0.45, 0.88);
        leg9_left->SetBorderSize(0);
        leg9_left->SetFillStyle(0);

        leg9_left->AddEntry(grB,"Power law exponent from fit: A#upointV_{over}^{B}","lep");

        if (n_B_good >= 2) {
            leg9_left->AddEntry(fit_B, "Constant fit: B = const", "l");
        }

        leg9_left->AddEntry(line_B, "Line: B = 2", "l");

        leg9_left->Draw();


        // ------------------------------------------------------------
        // Legenda a destra: risultati numerici
        // ------------------------------------------------------------

        if (n_B_good >= 2) {

            TLegend* leg9_right = new TLegend(0.62, 0.74, 0.88, 0.88);
            leg9_right->SetBorderSize(0);
            leg9_right->SetFillStyle(0);

            //leg9_right->AddEntry((TObject*)0, Form("B = %.2f #pm %.2f", fit_result_B, fit_error_B),"");

            // leg9_right->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f/%d", chi2_B, ndf_B),"");

            leg9_right->Draw();
        }

        c9->SaveAs("A1_bef_ann_lum_vs_v_B.png");
    }

    // ============================================================
    // CANVAS 10: chi square vs spot ID
    // ============================================================

    TCanvas* c10 = new TCanvas("c10", "Chi2/ndf vs spot", 800, 600);

    TGraph* grChi2 = new TGraph(n_fit, x.data(), chi2ndf_vals.data());

    grChi2->SetTitle("#chi^{2}/ndf;Spot;#chi^{2}/ndf");
    grChi2->SetMarkerStyle(22);
    grChi2->Draw("AP");
    grChi2->GetXaxis()->SetLimits(-1, nspots + 0.5);
    TLegend* leg10 = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg10->SetBorderSize(0);
    leg10->SetFillStyle(0);
    leg10->AddEntry(grChi2, "Qualit#grave{a} del fit", "p");
    leg10->Draw();

    c10->SaveAs("A1_bef_ann_lum_vs_v_chi2.png");
}