
//root -l 'lum_vs_V_fit_irradiated.C("A1_T=20_annealed_75_5_run=20260520-085312_total_global_ID.csv")'


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TObject.h"
#include "TAxis.h"
#include "TPad.h"

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

struct FitInfo {
    int spot = -1;
    vector<Row> rows;

    TGraphErrors* graph = nullptr;
    TF1* func = nullptr;

    bool fit_done = false;

    double A = 0.0;
    double A_err = 0.0;
    double B = 0.0;
    double B_err = 0.0;

    double chi2 = 0.0;
    int ndf = 0;
    double chi2ndf = 0.0;
    double prob = 0.0;

    int fit_status = -999;
    int covmatrix_status = -999;
    double edm = -1.0;
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

TPaveText* make_fit_box(const FitInfo& info,
                        double x1 = 0.12,
                        double y1 = 0.58,
                        double x2 = 0.45,
                        double y2 = 0.88) {

    TPaveText* fit_box = new TPaveText(x1, y1, x2, y2, "NDC");
    fit_box->SetFillColor(0);
    fit_box->SetBorderSize(1);
    fit_box->SetTextAlign(12);
    fit_box->SetTextSize(0.030);

    fit_box->AddText(Form("A = %.3g #pm %.3g", info.A, info.A_err));
    fit_box->AddText(Form("B = %.3g #pm %.3g", info.B, info.B_err));

    if (info.ndf > 0) {
        fit_box->AddText(Form("#chi^{2}/ndf = %.2f/%d = %.2f",
                              info.chi2,
                              info.ndf,
                              info.chi2ndf));
    } else {
        fit_box->AddText(Form("#chi^{2}/ndf = %.2f/%d",
                              info.chi2,
                              info.ndf));
    }

    fit_box->AddText(Form("p-value = %.3g", info.prob));
    fit_box->AddText(Form("EDM = %.3g", info.edm));
    fit_box->AddText(Form("Fit status = %d", info.fit_status));
    fit_box->AddText(Form("CovMatrix status = %d", info.covmatrix_status));

    return fit_box;
}

void print_fit_summary(const FitInfo& info) {
    cout << "----------------------------------------" << endl;
    cout << "Spot " << info.spot << endl;
    cout << "Fit status       = " << info.fit_status << endl;
    cout << "CovMatrix status = " << info.covmatrix_status << endl;
    cout << "EDM              = " << info.edm << endl;
    cout << "Chi2             = " << info.chi2 << endl;
    cout << "NDF              = " << info.ndf << endl;
    if (info.ndf > 0) {
        cout << "Chi2/NDF         = " << info.chi2ndf << endl;
    }
    cout << "p-value          = " << info.prob << endl;
    cout << "A                = " << info.A
         << " +/- " << info.A_err << endl;
    cout << "B                = " << info.B
         << " +/- " << info.B_err << endl;
}

void lum_vs_V_fit_irradiated(const char* filename = "data.csv") {

    gStyle->SetOptFit(111);

    // ------------------------------------------------------------
    // Global configuration
    // ------------------------------------------------------------
    const double FIT_XMIN = 0.0;
    const double FIT_XMAX = 8.0;

    const double B_INIT = 2.0;
    const double B_MIN  = 1.5;
    const double B_MAX  = 2.5;

    const double A_INIT = 1.0;

    const double CHI2NDF_CUT_FOR_B_PLOT = 2.0;

    const int max_spot_display = 8;
    const int ncols = 4;
    const int nrows = 2;

    // ------------------------------------------------------------
    // Read input CSV
    // ------------------------------------------------------------
    vector<Row> data = read_csv(filename);

    if (data.empty()) {
        cerr << "Errore: nessun dato letto da " << filename << endl;
        return;
    }

    map<int, vector<Row>> spots;

    for (const auto& r : data) {
        spots[r.spot].push_back(r);
    }

    int nspots = spots.size();

    // ============================================================
    // SINGLE FIT STAGE
    // Every spot is fitted only here.
    // All canvases below reuse these same fit results.
    // ============================================================

    map<int, FitInfo> fit_infos;

    vector<int> spot_ids;
    vector<double> A_vals, A_errs;
    vector<double> B_vals, B_errs;
    vector<double> chi2ndf_vals;

    cout << "\nPerforming one fit per spot only..." << endl;

    for (auto& entry : spots) {

        int spot = entry.first;
        vector<Row> rows = entry.second;

        sort(rows.begin(), rows.end(),
             [](const Row& a, const Row& b) {
                 return a.v_fin < b.v_fin;
             });

        FitInfo info;
        info.spot = spot;
        info.rows = rows;

        int n = rows.size();

        if (n < 2) {
            cout << "Skipping spot " << spot
                 << " because it has fewer than 2 points." << endl;
            fit_infos[spot] = info;
            continue;
        }

        vector<double> xv(n), yv(n), exv(n), eyv(n);

        for (int i = 0; i < n; i++) {
            xv[i]  = rows[i].v_fin;
            yv[i]  = rows[i].luminosity;
            exv[i] = 0.0;
            eyv[i] = rows[i].error;
        }

        info.graph = new TGraphErrors(n,
                                      xv.data(),
                                      yv.data(),
                                      exv.data(),
                                      eyv.data());

        info.graph->SetName(Form("gr_spot_%d", spot));
        info.graph->SetTitle(Form("Spot %d: x = %.2f, y = %.2f, T = %.1f #circC - annealing 75#circC 5h;Overvoltage (V);Luminosity",
                                  spot,
                                  rows[0].x,
                                  rows[0].y,
                                  rows[0].T));

        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(1.0);
        info.graph->SetLineWidth(2);

        // Fit: Lum = A * V_over^B
        // This is the only power-law fit performed for each spot.
        info.func = new TF1(Form("f_quad_spot_%d", spot),
                            "[0]*TMath::Power(x,[1])",
                            FIT_XMIN,
                            FIT_XMAX);

        info.func->SetParameters(A_INIT, B_INIT);
        info.func->SetParNames("A", "B");
        //info.func->SetParLimits(1, B_MIN, B_MAX);

        // S = save fit result
        // Q = quiet mode
        // R = use function range
        // 0 = do not draw during this fit stage
        TFitResultPtr fit_result = info.graph->Fit(info.func, "QRS0");

        info.fit_done = true;

        info.fit_status = int(fit_result);
        info.covmatrix_status = fit_result->CovMatrixStatus();
        info.edm = fit_result->Edm();

        info.A = info.func->GetParameter(0);
        info.A_err = info.func->GetParError(0);

        info.B = info.func->GetParameter(1);
        info.B_err = info.func->GetParError(1);

        info.chi2 = info.func->GetChisquare();
        info.ndf = info.func->GetNDF();
        info.chi2ndf = (info.ndf > 0) ? info.chi2 / info.ndf : 0.0;
        info.prob = info.func->GetProb();

        print_fit_summary(info);

        spot_ids.push_back(spot);
        A_vals.push_back(info.A);
        A_errs.push_back(info.A_err);
        B_vals.push_back(info.B);
        B_errs.push_back(info.B_err);
        chi2ndf_vals.push_back(info.chi2ndf);

        fit_infos[spot] = info;
    }

    int n_fit = spot_ids.size();

    if (n_fit == 0) {
        cerr << "Errore: nessun fit eseguito." << endl;
        return;
    }

    vector<double> x(n_fit), ex(n_fit, 0.0);
    for (int i = 0; i < n_fit; i++) {
        x[i] = spot_ids[i];
    }

    // ============================================================
    // CANVAS 5: luminosity vs spot at max v_fin
    // ============================================================

    double vmax = -1e9;

    for (const auto& r : data) {
        if (r.v_fin > vmax) {
            vmax = r.v_fin;
        }
    }

    vector<double> x_spot;
    vector<double> y_lum;
    vector<double> ey_lum;

    for (const auto& entry : spots) {

        int spot = entry.first;
        const vector<Row>& rows = entry.second;

        bool found = false;

        for (const auto& r : rows) {

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

    TGraphErrors* grMax = new TGraphErrors(nmax,
                                           x_spot.data(),
                                           y_lum.data(),
                                           ex_spot.data(),
                                           ey_lum.data());

    grMax->SetTitle(Form("Luminosity at maximum overvoltage (v_{fin}=%.2f V) - annealing 75#circC 5h;Spot;Luminosity",
                         vmax));

    grMax->SetMarkerStyle(20);
    grMax->SetMarkerSize(1.2);
    grMax->SetLineWidth(2);
    grMax->Draw("AP");

    c5->SaveAs("A1_ann_75_5_lum_vs_v_max_vfin.png");

    // ============================================================
    // CANVAS 5.1: luminosity vs spot at max v_fin - focus y range
    // ============================================================

    TCanvas* c51 = new TCanvas("c51",
                               Form("Luminosity at v_fin = %.2f V - focus", vmax),
                               1800,
                               600);

    TGraphErrors* grMax_c51 = new TGraphErrors(nmax,
                                               x_spot.data(),
                                               y_lum.data(),
                                               ex_spot.data(),
                                               ey_lum.data());

    grMax_c51->SetTitle(Form("Luminosity at maximum overvoltage (v_{fin}=%.2f V) - annealing 75#circC 5h - focus;Spot;Luminosity",
                             vmax));

    grMax_c51->SetMarkerStyle(20);
    grMax_c51->SetMarkerSize(1.2);
    grMax_c51->SetLineWidth(2);
    grMax_c51->SetMinimum(-5);
    grMax_c51->SetMaximum(1000);
    grMax_c51->Draw("AP");

    c51->SaveAs("A1_ann_75_5_lum_vs_v_max_vfin_focus.png");

    // ============================================================
    // CANVAS 6: luminosity vs v_fin per i primi 8 spot
    // Uses the stored fit functions. No new fit is performed here.
    // ============================================================

    TCanvas* c6 = new TCanvas("c6", "Luminosity vs overvoltage", 1600, 1000);
    c6->Divide(ncols, nrows);

    int displayed = 0;

    for (auto& entry : fit_infos) {

        if (displayed >= max_spot_display) break;

        FitInfo& info = entry.second;

        if (!info.fit_done || info.graph == nullptr || info.func == nullptr) continue;

        c6->cd(displayed + 1);

        info.graph->SetMarkerSize(1.0);
        info.graph->Draw("AP");
        info.func->Draw("SAME");

        displayed++;
    }

    c6->SaveAs("A1_ann_75_5_lum_vs_v.png");

    // ============================================================
    // CANVAS 7: luminosity vs V_over for all spots with fit
    // Uses the stored fit functions. No new fit is performed here.
    // ============================================================

    gSystem->mkdir("A1_ann_75_5_lum_vs_v_all_spots", kTRUE);

    for (auto& entry : fit_infos) {

        int spot_id = entry.first;
        FitInfo& info = entry.second;

        if (!info.fit_done || info.graph == nullptr || info.func == nullptr) {
            cout << "Skipping Canvas 7 for spot " << spot_id
                 << " because no valid fit is available." << endl;
            continue;
        }

        TCanvas* c7_spot = new TCanvas(Form("c7_spot_%d", spot_id),
                                       Form("Spot %d - Luminosity vs overvoltage", spot_id),
                                       800,
                                       600);

        info.graph->SetMarkerStyle(20);
        info.graph->SetMarkerSize(2.0);
        info.graph->SetLineWidth(2);
        info.graph->GetXaxis()->SetTitleSize(0.045);
        info.graph->GetYaxis()->SetTitleSize(0.045);

        info.graph->Draw("AP");
        info.func->Draw("SAME");

        TPaveText* fit_info_box = make_fit_box(info, 0.12, 0.58, 0.45, 0.88);
        fit_info_box->Draw("SAME");

        c7_spot->SaveAs(Form("A1_ann_75_5_lum_vs_v_all_spots/A1_ann_75_5_lum_vs_v_spot%d.png",
                             spot_id));

        delete c7_spot;
    }

    // ============================================================
    // CANVAS 9: B vs spot ID con selezione su chi2/ndf
    // Uses the stored B values from the single fit stage.
    // The constant fit B = const is a different fit, performed only here.
    // ============================================================

    TCanvas* c9 = new TCanvas("c9", "B vs spot", 1800, 600);

    vector<double> x_B_good;
    vector<double> ex_B_good;
    vector<double> B_good;
    vector<double> B_err_good;

    for (int i = 0; i < n_fit; i++) {

        if (chi2ndf_vals[i] < CHI2NDF_CUT_FOR_B_PLOT) {

            x_B_good.push_back(x[i]);
            ex_B_good.push_back(0.0);
            B_good.push_back(B_vals[i]);
            B_err_good.push_back(B_errs[i]);
        }
    }

    int n_B_good = B_good.size();

    cout << "Canvas 9: numero di fit con chi2/ndf < "
         << CHI2NDF_CUT_FOR_B_PLOT
         << " = " << n_B_good << " su " << n_fit << endl;

    if (n_B_good == 0) {

        cout << "ATTENZIONE: nessun punto con chi2/ndf < "
             << CHI2NDF_CUT_FOR_B_PLOT
             << ". Canvas 9 non disegnata." << endl;

    } else {

        TGraphErrors* grB = new TGraphErrors(n_B_good,
                                             x_B_good.data(),
                                             B_good.data(),
                                             ex_B_good.data(),
                                             B_err_good.data());

        grB->SetTitle("Power law exponent: B - annealing 75#circC 5h;Spot;B");
        grB->SetMarkerStyle(21);
        grB->SetMarkerSize(1.2);
        grB->Draw("AP");
        grB->GetXaxis()->SetLimits(-1, nspots + 0.5);

        // Linea B = 2
        TF1* line_B = new TF1("line_B", "2", -1, nspots + 0.5);
        line_B->SetLineColor(kRed + 1);
        line_B->SetLineStyle(2);
        line_B->SetLineWidth(2);
        line_B->Draw("SAME");

        // Fit costante B = const
        TF1* fit_B = new TF1("fit_B", "[0]", -1, nspots + 0.5);
        fit_B->SetLineColor(kMagenta + 1);
        fit_B->SetLineWidth(2);

        double chi2_B = 0.0;
        int ndf_B = 0;
        double fit_result_B = 0.0;
        double fit_error_B = 0.0;

        if (n_B_good >= 2) {
            grB->Fit(fit_B, "RQ");
            fit_B->Draw("SAME");

            chi2_B = fit_B->GetChisquare();
            ndf_B = fit_B->GetNDF();
            fit_result_B = fit_B->GetParameter(0);
            fit_error_B = fit_B->GetParError(0);
        } else {
            cout << "ATTENZIONE: c'e' solo un punto buono. "
                 << "Fit costante non eseguito." << endl;
        }

        TLegend* leg9_left = new TLegend(0.12, 0.72, 0.45, 0.88);
        leg9_left->SetBorderSize(0);
        leg9_left->SetFillStyle(0);
        leg9_left->AddEntry(grB, "Power law exponent from fit: A#upointV_{over}^{B}", "lep");

        if (n_B_good >= 2) {
            leg9_left->AddEntry(fit_B, "Constant fit: B = const", "l");
        }

        leg9_left->AddEntry(line_B, "Line: B = 2", "l");
        leg9_left->Draw();

        if (n_B_good >= 2) {

            TLegend* leg9_right = new TLegend(0.62, 0.74, 0.88, 0.88);
            leg9_right->SetBorderSize(0);
            leg9_right->SetFillStyle(0);
            //leg9_right->AddEntry((TObject*)0,Form("B_{const} = %.3f #pm %.3f", fit_result_B, fit_error_B), "");
            //leg9_right->AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.2f/%d", chi2_B, ndf_B),"");
            leg9_right->Draw();
        }

        c9->SaveAs("A1_ann_75_5_lum_vs_v_B.png");
    }

    // ============================================================
    // CANVAS 10: chi square vs spot ID
    // Uses chi2/ndf from the single fit stage.
    // ============================================================

    TCanvas* c10 = new TCanvas("c10", "Chi2/ndf vs spot", 800, 600);

    TGraph* grChi2 = new TGraph(n_fit,
                                x.data(),
                                chi2ndf_vals.data());

    grChi2->SetTitle("#chi^{2}/ndf;Spot;#chi^{2}/ndf");
    grChi2->SetMarkerStyle(22);
    grChi2->SetMarkerSize(1.2);
    grChi2->Draw("AP");
    grChi2->GetXaxis()->SetLimits(-1, nspots + 0.5);

    TLegend* leg10 = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg10->SetBorderSize(0);
    leg10->SetFillStyle(0);
    leg10->AddEntry(grChi2, "Qualit#grave{a} del fit", "p");
    leg10->Draw();

    c10->SaveAs("A1_ann_75_5_lum_vs_v_chi2.png");
}
