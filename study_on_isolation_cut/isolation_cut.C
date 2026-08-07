//
// Input CSV format:
// spot,x,y,luminosity,error,T,v,v_fin
//
// Run with:
// root -l 'isolation_cut.C'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TH1F.h>
#include <TMath.h>

namespace {

// ============================================================
// CONFIGURATION
// ============================================================

const std::string CSV_FILE_1 =
    "DATA/bef_ann_A1_T=20_isolation/luminosity/A1_T=20_isolation_Rchange.csv";

const std::string CSV_FILE_2 =
    "DATA/bef_ann_A1_T=20_no_isolation/luminosity/A1_T=20_no_isolation_Rchange.csv";

const std::string LABEL_FILE_1 = "Isolated defects";
const std::string LABEL_FILE_2 = "Non-isolated defects";

const std::string OUTPUT_DIR = "Results";

// Histogram range and binning
constexpr double LUMINOSITY_MIN = 0.0;
constexpr double LUMINOSITY_MAX = 5500.0;
constexpr double BIN_WIDTH      = 100.0;

constexpr int N_BINS = static_cast<int>((LUMINOSITY_MAX - LUMINOSITY_MIN) / BIN_WIDTH);

constexpr double LOG_LUMINOSITY_MIN = 1.0;
constexpr double LOG_LUMINOSITY_MAX = 10000.0;
constexpr int LOG_N_BINS = 30;

// ============================================================

std::string Trim(const std::string& input)
{
    const auto first = std::find_if_not(
        input.begin(),
        input.end(),
        [](unsigned char c) {
            return std::isspace(c);
        }
    );

    const auto last = std::find_if_not(
        input.rbegin(),
        input.rend(),
        [](unsigned char c) {
            return std::isspace(c);
        }
    ).base();

    if (first >= last) {
        return "";
    }

    return std::string(first, last);
}


std::vector<std::string> SplitCSVLine(const std::string& line)
{
    std::vector<std::string> fields;

    std::stringstream stream(line);
    std::string field;

    while (std::getline(stream, field, ',')) {
        fields.push_back(Trim(field));
    }

    return fields;
}


int FindColumn(const std::vector<std::string>& header,const std::string& columnName)
{
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (Trim(header[i]) == columnName) {
            return static_cast<int>(i);
        }
    }

    return -1;
}


// Read the luminosity column and v_fin for filtering.
std::vector<double> ReadLuminosities(const std::string& filename)
{
    std::ifstream input(filename);

    if (!input.is_open()) {
        throw std::runtime_error(
            "Cannot open CSV file: " + filename
        );
    }

    std::string line;

    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "CSV file is empty: " + filename
        );
    }

    // Remove a possible UTF-8 BOM
    if (line.size() >= 3 &&
        static_cast<unsigned char>(line[0]) == 0xEF &&
        static_cast<unsigned char>(line[1]) == 0xBB &&
        static_cast<unsigned char>(line[2]) == 0xBF) {

        line.erase(0, 3);
    }

    const std::vector<std::string> header = SplitCSVLine(line);

    const int luminosityColumn = FindColumn(header, "luminosity");
    const int vFinColumn = FindColumn(header, "v_fin");

    if (luminosityColumn < 0 || vFinColumn < 0) {
        throw std::runtime_error("Required columns not found in: " + filename);
    }

    std::vector<double> luminosities;

    int lineNumber = 1;
    int outsideRange = 0;

    while (std::getline(input, line)) {
        ++lineNumber;

        line = Trim(line);

        if (line.empty()) {
            continue;
        }

        const std::vector<std::string> fields =
            SplitCSVLine(line);

        const int maximumRequiredColumn =
        std::max(luminosityColumn, vFinColumn);

        if (static_cast<int>(fields.size()) <= maximumRequiredColumn) {

            std::cerr
                << "Warning: incomplete line "
                << lineNumber
                << " in "
                << filename
                << std::endl;

            continue;
        }

        try {

            const double vFin = std::stod(fields[vFinColumn]);

            // Keep only rows with v_fin = 5
            if (std::abs(vFin - 5.0) > 1e-6) {
                    continue;
                }

            const double luminosity =
                std::stod(fields[luminosityColumn]);

            if (!std::isfinite(luminosity)) {
                std::cerr
                    << "Warning: non-finite luminosity at line "
                    << lineNumber
                    << " in "
                    << filename
                    << std::endl;

                continue;
            }

            if (luminosity < LUMINOSITY_MIN ||
                luminosity >= LUMINOSITY_MAX) {

                ++outsideRange;
                continue;
            }

            luminosities.push_back(luminosity);
        }
        catch (const std::exception&) {
            std::cerr
                << "Warning: invalid luminosity at line "
                << lineNumber
                << " in "
                << filename
                << std::endl;
        }
    }

    std::cout
        << "Read "
        << luminosities.size()
        << " luminosity values from "
        << filename
        << std::endl;

    if (outsideRange > 0) {
        std::cout
            << "Ignored "
            << outsideRange
            << " values outside the range ["
            << LUMINOSITY_MIN
            << ", "
            << LUMINOSITY_MAX
            << ")."
            << std::endl;
    }

    return luminosities;
}

// Create and fill the histogram.

TH1D* CreateHistogram(const char* histogramName, const std::vector<double>& luminosities)
{
    TH1D* histogram = new TH1D(histogramName,"", N_BINS, LUMINOSITY_MIN, LUMINOSITY_MAX);

    histogram->SetDirectory(nullptr);
    histogram->Sumw2();

    for (const double luminosity : luminosities) {
        histogram->Fill(luminosity);
    }

    histogram->GetXaxis()->SetTitle("Luminosity");

    histogram->GetYaxis()->SetTitle(
        Form(
            "Number of spots / %.0f luminosity units",
            BIN_WIDTH
        )
    );

    histogram->GetXaxis()->SetTitleSize(0.045);
    histogram->GetYaxis()->SetTitleSize(0.045);

    histogram->GetXaxis()->SetLabelSize(0.040);
    histogram->GetYaxis()->SetLabelSize(0.040);

    histogram->GetXaxis()->SetTitleOffset(1.25);
    histogram->GetYaxis()->SetTitleOffset(1.25);

    histogram->SetLineWidth(2);
    histogram->SetMarkerStyle(20);
    histogram->SetMarkerSize(0.8);

    return histogram;
}


double GetHistogramMedian(TH1D* histogram)
{
    double probability = 0.5;
    double median = 0.0;

    histogram->GetQuantiles(1, &median, &probability);

    return median;
}


TCanvas* CreateIndividualCanvas(const char* canvasName, TH1D* histogram, const std::string& datasetLabel, Color_t histogramColor)
{
    TCanvas* canvas = new TCanvas(canvasName, canvasName, 1800, 1200);

    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.05);
    canvas->SetBottomMargin(0.12);
    canvas->SetTopMargin(0.08);

    canvas->SetTicks(1, 1);

    histogram->SetTitle(Form(
            "%s;Luminosity;"
            "Number of spots / %.0f luminosity units",
            datasetLabel.c_str(),
            BIN_WIDTH
        )
    );

    histogram->SetFillColorAlpha(histogramColor, 0.35);
    histogram->SetFillStyle(1001);
    histogram->SetLineColor(histogramColor);
    histogram->SetMarkerColor(histogramColor);

    histogram->SetMinimum(0.0);

    double maximumY = histogram->GetMaximum();

    if (maximumY <= 0.0) {
        maximumY = 1.0;
    }

    histogram->SetMaximum(1.35 * maximumY);

    // E1 displays the Poisson counting uncertainties.
    histogram->Draw("HIST E1");

    const double mean =
        histogram->GetMean();

    const double median =
        GetHistogramMedian(histogram);

    const double numberOfEntries =
        histogram->GetEntries();

    TLegend* legend = new TLegend(
        0.60,
        0.67,
        0.88,
        0.88
    );

    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);

    legend->AddEntry(histogram, datasetLabel.c_str(), "lep");

    legend->AddEntry(static_cast<TObject*>(nullptr), Form("Entries = %.0f", numberOfEntries),"");

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Mean = %.2f", mean),
        ""
    );

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Median = %.2f", median),
        ""
    );

    legend->Draw();

    canvas->Modified();
    canvas->Update();

    return canvas;
}

// Create normalized copies for the ratio calculation.
TH1D* CreateNormalizedCopy( TH1D* source, const char* name)
{
    TH1D* normalized = static_cast<TH1D*>(source->Clone(name));

    normalized->SetDirectory(nullptr);

    const double integral = normalized->Integral();

    if (integral > 0.0) {
        normalized->Scale(1.0 / integral);
    }

    return normalized;
}


void HistoUtils_BinLogX(TH1* histogram)
{
    TAxis* axis = histogram->GetXaxis();

    const Int_t bins = axis->GetNbins();
    const Axis_t from = axis->GetXmin();
    const Axis_t to = axis->GetXmax();
    const Axis_t width = (to - from) / bins;

    Axis_t* newBins = new Axis_t[bins + 1];

    for (Int_t i = 0; i <= bins; ++i) {
        newBins[i] = TMath::Power(
            10.0,
            from + i * width
        );
    }

    axis->Set(bins, newBins);

    delete[] newBins;
}
}

// ============================================================
// MAIN MACRO
// ============================================================

void isolation_cut()
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.045);

    gSystem->mkdir(OUTPUT_DIR.c_str(), kTRUE);

    try {
        // Read luminosity values from the two CSV files
        const std::vector<double> luminosities1 =
            ReadLuminosities(CSV_FILE_1);

        const std::vector<double> luminosities2 =
            ReadLuminosities(CSV_FILE_2);

        if (luminosities1.empty()) {
            throw std::runtime_error(
                "The first CSV file contains no valid luminosity "
                "values in the selected range."
            );
        }

        if (luminosities2.empty()) {
            throw std::runtime_error(
                "The second CSV file contains no valid luminosity "
                "values in the selected range."
            );
        }

        // Create the two original histograms
        TH1D* histogram1 = CreateHistogram("histogram_file1",luminosities1);

        TH1D* histogram2 = CreateHistogram("histogram_file2",luminosities2);

        // ====================================================
        // Canvas 1: first CSV
        // ====================================================
        /*
        TCanvas* canvas1 = CreateIndividualCanvas(
            "canvas_file1",
            histogram1,
            LABEL_FILE_1,
            kBlue + 1
        );

        // ====================================================
        // Canvas 2: second CSV
        // ====================================================

        TCanvas* canvas2 = CreateIndividualCanvas(
            "canvas_file2",
            histogram2,
            LABEL_FILE_2,
            kGreen + 2
        );
        */
        
        // ====================================================
        // Canvas 4: non-normalized histogram comparison
        // ====================================================

        TCanvas* canvas4 = new TCanvas(
            "canvas_non_normalized_comparison",
            "Non-normalized luminosity distributions",
            1200,
            800
        );

        canvas4->SetLeftMargin(0.12);
        canvas4->SetRightMargin(0.05);
        canvas4->SetBottomMargin(0.12);
        canvas4->SetTopMargin(0.08);
        canvas4->SetTicks(1, 1);

        // First histogram: blue
        histogram1->SetLineColor(kBlue + 1);
        histogram1->SetMarkerColor(kBlue + 1);
        histogram1->SetFillColorAlpha(kBlue + 1, 0.35);
        histogram1->SetFillStyle(1001);
        histogram1->SetLineWidth(2);

        // Second histogram: green
        histogram2->SetLineColor(kGreen + 2);
        histogram2->SetMarkerColor(kGreen + 2);
        histogram2->SetFillColorAlpha(kGreen + 2, 0.35);
        histogram2->SetFillStyle(1001);
        histogram2->SetLineWidth(2);

        histogram1->SetTitle(
            Form(
                "Non-normalized luminosity distributions;"
                "Luminosity;"
                "Number of spots"
            )
        );

        // Set a common y-axis range
        const double maximumCanvas4 = std::max(
            histogram1->GetMaximum(),
            histogram2->GetMaximum()
        );

        histogram1->SetMinimum(0.0);
        histogram1->SetMaximum(1.1 * maximumCanvas4);

        // Draw the two non-normalized histograms
        histogram1->Draw("HIST E1");
        histogram2->Draw("HIST E1 SAME");

        TLegend* legend4 = new TLegend(
            0.62,
            0.72,
            0.88,
            0.87
        );

        legend4->SetBorderSize(0);
        legend4->SetFillStyle(0);
        legend4->SetTextSize(0.035);

        legend4->AddEntry(
            histogram1,
            Form("%s, N = %.0f",
                LABEL_FILE_1.c_str(),
                histogram1->GetEntries()),
            "lf"
        );

        legend4->AddEntry(
            histogram2,
            Form("%s, N = %.0f",
                LABEL_FILE_2.c_str(),
                histogram2->GetEntries()),
            "lf"
        );

        legend4->Draw();

        canvas4->Modified();
        canvas4->Update();


        // ====================================================
        // Logarithmically binned histograms
        // ====================================================


        TH1D* histogramLog1 = new TH1D(
            "histogram_log_file1",
            "",
            LOG_N_BINS,
            std::log10(LOG_LUMINOSITY_MIN),
            std::log10(LOG_LUMINOSITY_MAX)
        );

        TH1D* histogramLog2 = new TH1D(
            "histogram_log_file2",
            "",
            LOG_N_BINS,
            std::log10(LOG_LUMINOSITY_MIN),
            std::log10(LOG_LUMINOSITY_MAX)
        );

        histogramLog1->SetDirectory(nullptr);
        histogramLog2->SetDirectory(nullptr);

        // Calling HistoUtils_BinLogX to set logarithmic binning on the x-axis
        HistoUtils_BinLogX(histogramLog1);
        HistoUtils_BinLogX(histogramLog2);

        histogramLog1->Sumw2();
        histogramLog2->Sumw2();

        // Filling the histograms
        for (const double luminosity : luminosities1) {
            if (luminosity >= LOG_LUMINOSITY_MIN &&
                luminosity < LOG_LUMINOSITY_MAX) {

                histogramLog1->Fill(luminosity);
            }
        }

        for (const double luminosity : luminosities2) {
            if (luminosity >= LOG_LUMINOSITY_MIN &&
                luminosity < LOG_LUMINOSITY_MAX) {

                histogramLog2->Fill(luminosity);
            }
        }

        // Setting the color
        histogramLog1->SetLineColor(kBlue + 1);
        histogramLog1->SetFillColorAlpha(kBlue + 1, 0.35);
        histogramLog1->SetMarkerColor(kBlue + 1);
        histogramLog1->SetLineWidth(2);
        histogramLog1->SetMarkerStyle(20);
        histogramLog1->SetMarkerSize(0.75);

        histogramLog2->SetLineColor(kGreen + 2);
        histogramLog2->SetFillColorAlpha(kGreen + 2, 0.35);
        histogramLog2->SetMarkerColor(kGreen + 2);
        histogramLog2->SetLineWidth(2);
        histogramLog2->SetMarkerStyle(21);
        histogramLog2->SetMarkerSize(0.75);

        histogramLog1->GetXaxis()->SetTitle("Luminosity");
        histogramLog1->GetXaxis()->SetTitleOffset(1.4);
        histogramLog1->GetYaxis()->SetTitle("Number of spots");

        histogramLog2->GetXaxis()->SetTitle("Luminosity");
        histogramLog2->GetYaxis()->SetTitle("Number of spots");


        // ====================================================
        // Canvas 5: normalized logarithmically binned histograms
        // ====================================================

        TCanvas* canvas5 = new TCanvas("canvas_log_normalized","Normalized logarithmically binned histograms",
            1400,
            800
        );

        canvas5->SetLeftMargin(0.12);
        canvas5->SetRightMargin(0.05);
        canvas5->SetBottomMargin(0.10);
        canvas5->SetTopMargin(0.08);
        canvas5->SetTicks(1, 1);
        canvas5->SetLogx();

        //
        TH1D* histogramLog1Norm = CreateNormalizedCopy(histogramLog1, "histogramLog1_norm");
        TH1D* histogramLog2Norm = CreateNormalizedCopy(histogramLog2, "histogramLog2_norm");

        histogramLog1Norm->SetTitle(
            "Normalized luminosity distributions with logarithmic binning;"
            "Luminosity;"
            "Normalized number of spots");

        const double maximumLogNormalized = std::max(
            histogramLog1Norm->GetMaximum(),
            histogramLog2Norm->GetMaximum());

        histogramLog1Norm->SetMinimum(0.0);
        histogramLog1Norm->GetYaxis()->SetRangeUser(0.0, 1.4 * maximumLogNormalized);
        histogramLog1Norm->GetXaxis()->SetRangeUser(5.0, std::pow(10.0, 5.0));
        histogramLog2Norm->GetXaxis()->SetRangeUser(5.0, std::pow(10.0, 5.0));
        histogramLog1Norm->GetXaxis()->SetTitleOffset(1.4);

        // Disegniamo con "E1" senza "HIST" per evitare sovrapposizioni di riempimento
        histogramLog1Norm->Draw("HIST E1");
        histogramLog2Norm->Draw("HIST E1 SAME");

        TLegend* legend5 = new TLegend(0.62, 0.72, 0.88, 0.87);
        legend5->SetBorderSize(0);
        legend5->SetFillStyle(0);
        legend5->SetTextSize(0.035);

        legend5->AddEntry(histogramLog1Norm, LABEL_FILE_1.c_str(), "lep");
        legend5->AddEntry(histogramLog2Norm, LABEL_FILE_2.c_str(), "lep");
        legend5->Draw();

        canvas5->Modified();
        canvas5->Update();

        // ====================================================
        // Canvas 6: non-normalized logarithmically binned histograms
        // ====================================================

        TCanvas* canvas6 = new TCanvas(
            "canvas_log_non_normalized",
            "Non-normalized logarithmically binned histograms",
            1400,
            800);

        canvas6->SetLeftMargin(0.12);
        canvas6->SetRightMargin(0.05);
        canvas6->SetBottomMargin(0.12);
        canvas6->SetTopMargin(0.08);
        canvas6->SetTicks(1, 1);
        canvas6->SetLogx();

        histogramLog1->SetTitle(
            "Non-normalized luminosity distributions with logarithmic binning;"
            "Luminosity;"
            "Number of spots");
        
        histogramLog1->SetFillStyle(1001);
        histogramLog2->SetFillStyle(1001);

        const double maximumLogNonNormalized = std::max(
            histogramLog1->GetMaximum(),
            histogramLog2->GetMaximum());

        histogramLog1->SetMinimum(0.0);
        histogramLog1->GetYaxis()->SetRangeUser(0.0, 1.4 * maximumLogNonNormalized);
        histogramLog1->GetXaxis()->SetRangeUser(5.0, std::pow(10.0, 6.0));
        histogramLog2->GetYaxis()->SetRangeUser(0.0, 1.4 * maximumLogNonNormalized);
        histogramLog2->GetXaxis()->SetRangeUser(5.0, std::pow(10.0, 6.0));


        // Draw the bins with filled area
        histogramLog1->Draw("HIST");
        histogramLog2->Draw("HIST SAME");

        //Draw also the error bars
        histogramLog1->Draw("E1 SAME");
        histogramLog2->Draw("E1 SAME");


        TLegend* legend6 = new TLegend(0.62, 0.72, 0.88, 0.87);
        legend6->SetBorderSize(0);
        legend6->SetFillStyle(0);
        legend6->SetTextSize(0.035);

        legend6->AddEntry(
            histogramLog1,
            Form("%s, N = %.0f", LABEL_FILE_1.c_str(), histogramLog1->GetEntries()),
            "lep");

        legend6->AddEntry(
            histogramLog2,
            Form("%s, N = %.0f", LABEL_FILE_2.c_str(), histogramLog2->GetEntries()),
            "lep");

        legend6->Draw();

        canvas6->Modified();
        canvas6->Update();

         TCanvas* canvas = new TCanvas(
        "canvas_normalized_comparison",
        "Normalized luminosity distributions and ratio",
        1800,
        1200);
        
        // ========================================================
        // Canvas 7: normalized logarithmically binned histograms and ratio
        // ========================================================

        // ========================================================
        // Upper pad: normalized histograms
        // ========================================================

        TPad* upperPad = new TPad("upper_pad", "Normalized logarithmically binned histograms",
            0.0,
            0.30,
            1.0,
            1.0);

        upperPad->SetLeftMargin(0.12);
        upperPad->SetRightMargin(0.05);
        upperPad->SetTopMargin(0.08);
        upperPad->SetBottomMargin(0.025);
        upperPad->SetTicks(1, 1);
        upperPad->Draw();

        // ========================================================
        // Lower pad: ratio
        // ========================================================

        TPad* lowerPad = new TPad("lower_pad", "Ratio",
            0.0,
            0.0,
            1.0,
            0.30);

        lowerPad->SetLeftMargin(0.12);
        lowerPad->SetRightMargin(0.05);
        lowerPad->SetTopMargin(0.025);
        lowerPad->SetBottomMargin(0.32);
        lowerPad->SetTicks(1, 1);
        lowerPad->SetGridy();
        lowerPad->Draw();

        // ========================================================
        // Draw the normalized logarithmically binned histograms
        // ========================================================

        upperPad->cd();
        upperPad->SetLogx();

        histogramLog1Norm->SetLineColor(kBlue + 1);
        histogramLog1Norm->SetMarkerColor(kBlue + 1);
        histogramLog1Norm->SetMarkerStyle(20);
        histogramLog1Norm->SetMarkerSize(0.75);
        histogramLog1Norm->SetLineWidth(2);

        histogramLog2Norm->SetLineColor(kGreen + 2);
        histogramLog2Norm->SetMarkerColor(kGreen + 2);
        histogramLog2Norm->SetMarkerStyle(21);
        histogramLog2Norm->SetMarkerSize(0.75);
        histogramLog2Norm->SetLineWidth(2);

        histogramLog1Norm->SetTitle(
            "Normalized luminosity distributions, logarithmically binned;"
            "Luminosity;"
            "Normalized number of spots");

        const double normalizedMaximum = std::max(
            histogramLog1Norm->GetMaximum(),
            histogramLog2Norm->GetMaximum()
        );

        histogramLog1Norm->SetMinimum(0.0);

        histogramLog1Norm->SetMaximum(
            normalizedMaximum > 0.0
                ? 1.1 * normalizedMaximum
                : 1.0
        );

        histogramLog1Norm->GetXaxis()->SetLabelSize(0.0);
        histogramLog1Norm->GetXaxis()->SetTitleSize(0.0);
        histogramLog1Norm->GetYaxis()->SetTitleSize(0.055);
        histogramLog1Norm->GetYaxis()->SetLabelSize(0.050);
        histogramLog1Norm->GetYaxis()->SetTitleOffset(1.05);

        histogramLog1Norm->Draw("HIST E1");
        histogramLog2Norm->Draw("HIST E1 SAME");

        TLegend* legend = new TLegend(
            0.58,
            0.67,
            0.88,
            0.88
        );

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);

        legend->AddEntry(histogramLog1Norm, LABEL_FILE_1.c_str(), "lep");
        legend->AddEntry(histogramLog2Norm, LABEL_FILE_2.c_str(), "lep");
        legend->Draw();

        upperPad->Modified();
        upperPad->Update();

        // ========================================================
        // Ratio of the normalized histograms
        // ========================================================

        lowerPad->cd();
        lowerPad->SetLogx();

        TH1D* histogram1RatioInput = static_cast<TH1D*>(histogramLog1Norm->Clone("histogram1_ratio_input"));
        histogram1RatioInput->SetDirectory(nullptr);
        
        TH1D* histogram2RatioInput = static_cast<TH1D*>(histogramLog2Norm->Clone("histogram2_ratio_input"));
        histogram2RatioInput->SetDirectory(nullptr);

        TH1D* ratio =static_cast<TH1D*>(histogram1RatioInput->Clone("ratio_normalized_histograms"));
        ratio->SetDirectory(nullptr);

        /*
        * Divide calculates:
        *
        *    normalized histogram 1
        *    ----------------------
        *    normalized histogram 2
        *
        * The statistical errors are propagated by ROOT.
        *
        * Bins in which the second histogram is zero are assigned
        * ratio = 0 by TH1::Divide().
        */
        ratio->Divide(histogram2RatioInput);

        ratio->SetTitle("");

        ratio->GetXaxis()->SetTitle("Luminosity");

        ratio->GetYaxis()->SetTitle("Isolated/Non-isolated" );

        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->GetXaxis()->SetTitleOffset(1.05);

        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetLabelSize(0.085);
        ratio->GetYaxis()->SetTitleOffset(0.55);

        ratio->GetYaxis()->SetNdivisions(505);

        ratio->SetLineColor(kBlack);
        ratio->SetMarkerColor(kBlack);
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(0.75);
        ratio->SetLineWidth(2);

        /*
        * Calculate a useful display range from the non-zero,
        * finite ratio bins.
        */
        double maximumRatio = 0.0;

        for (int bin = 1;
            bin <= ratio->GetNbinsX();
            ++bin) {

            const double value =
                ratio->GetBinContent(bin);

            const double error =
                ratio->GetBinError(bin);

            if (std::isfinite(value) &&
                std::isfinite(error) &&
                value > 0.0) {

                maximumRatio = std::max(
                    maximumRatio,
                    value + error
                );
            }
        }

        ratio->SetMinimum(0.0);

        if (maximumRatio > 0.0) {
            ratio->SetMaximum(
                std::max(2.0, 1.25 * maximumRatio)
            );
        }
        else {
            ratio->SetMaximum(2.0);
        }

        // Draw the ratio with propagated y-errors.
        ratio->Draw("E1");

        TLine* unityLine = new TLine(
            0,
            1.0,
            LOG_LUMINOSITY_MAX,
            1.0
        );

        unityLine->SetLineStyle(2);
        unityLine->SetLineWidth(2);
        unityLine->Draw("SAME");

        lowerPad->Modified();
        lowerPad->Update();

        canvas->cd();
        canvas->Modified();
        canvas->Update();
        
        // canvas3->SaveAs(Form("%s/normalized_histograms_ratio.png",OUTPUT_DIR.c_str()));

        canvas4->SaveAs(Form("%s/non_normalized_histograms.png", OUTPUT_DIR.c_str()));

        canvas5->SaveAs(Form("%s/log_normalized_histograms.png", OUTPUT_DIR.c_str()));

        canvas6->SaveAs(Form("%s/log_non_normalized_histograms.png", OUTPUT_DIR.c_str()));
    }

    
    catch (const std::exception& exception) {
        std::cerr
            << "\nError: "
            << exception.what()
            << std::endl;

            }  

    
    
}
