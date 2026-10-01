// ============================================================================
// luminosity_distribution.C
//
// Input files are the compact CSV files produced by
// 6luminosity_distribution.sh, with columns:
//
//     luminosity,error,deltaL
//
// The statistical and systematic luminosity uncertainties are read and
// validated. They are not used as vertical histogram errors: the vertical
// uncertainty of a luminosity-distribution bin is the counting uncertainty.
//
// Actions:
//   individual     -> one non-normalized log-binned distribution
//   new_before     -> normalized A1 new-device vs before-annealing overlay
//   phase_compare  -> normalized before vs after annealing overlay + ratio
// ============================================================================

#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObject.h>
#include <TPad.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ============================================================================
// CONFIGURATION
// ============================================================================

// The reference macro used 30 bins over four decades (1 -> 1e4), i.e.
// 7.5 bins/decade.  Here the same density is extended to 1e6 so that brighter
// irradiated hotspots are not silently lost.  These limits are common to every
// canvas, which makes the plots directly comparable.
constexpr double LOG_LUMINOSITY_MIN = 1.0;
constexpr double LOG_LUMINOSITY_MAX = 1.0e5;
constexpr int    LOG_N_BINS         = 30;

// Display-only x limit: this adds empty graphical space on the right.
// It does NOT alter the 30 physical logarithmic bins above.
constexpr double DISPLAY_LUMINOSITY_MAX = 9.0e5;

struct Measurement {
    double luminosity = 0.0;
    double statError  = 0.0;
    double systError  = 0.0;
};

// ============================================================================
// SMALL UTILITIES
// ============================================================================

std::string Trim(const std::string& input)
{
    const auto first = std::find_if_not(
        input.begin(), input.end(),
        [](unsigned char c) { return std::isspace(c); }
    );

    const auto last = std::find_if_not(
        input.rbegin(), input.rend(),
        [](unsigned char c) { return std::isspace(c); }
    ).base();

    if (first >= last) return "";
    return std::string(first, last);
}

std::vector<std::string> SplitCSVLine(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;

    while (std::getline(ss, field, ',')) {
        fields.push_back(Trim(field));
    }
    return fields;
}

int FindColumn(const std::vector<std::string>& header,
               const std::string& columnName)
{
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (Trim(header[i]) == columnName) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

double MeanLuminosity(const std::vector<Measurement>& data)
{
    if (data.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& m : data)
        sum += m.luminosity;

    return sum / data.size();
}

double MedianLuminosity(const std::vector<Measurement>& data)
{
    if (data.empty()) return 0.0;

    std::vector<double> values;
    values.reserve(data.size());

    for (const auto& m : data)
        values.push_back(m.luminosity);

    std::sort(values.begin(), values.end());

    const size_t n = values.size();

    if (n % 2 == 0)
        return 0.5 * (values[n/2 - 1] + values[n/2]);

    return values[n/2];
}

std::string SafeName(std::string s)
{
    for (char& c : s) {
        const unsigned char uc = static_cast<unsigned char>(c);
        if (!(std::isalnum(uc) || c == '_' || c == '-' || c == '=')) {
            c = '_';
        }
    }
    return s;
}

std::string PhaseLabel(const std::string& phase)
{
    if (phase == "before_annealing") {
        return "Before annealing";
    }

    if (phase == "new_device") {
        return "New device";
    }

    // Expected form: annealing_T=75_h=5
    const std::string prefix = "annealing_T=";
    const std::size_t pT = phase.find(prefix);
    const std::size_t pH = phase.find("_h=");

    if (pT == 0 && pH != std::string::npos && pH > prefix.size()) {
        const std::string temperature =
            phase.substr(prefix.size(), pH - prefix.size());
        const std::string hours = phase.substr(pH + 3);

        // ROOT TLatex-style degree symbol.
        return temperature + " #circC, " + hours + " h";
    }

    return phase;
}

std::vector<double> MakeLogEdges()
{
    std::vector<double> edges(LOG_N_BINS + 1, 0.0);

    const double logMin = std::log10(LOG_LUMINOSITY_MIN);
    const double logMax = std::log10(LOG_LUMINOSITY_MAX);
    const double step   = (logMax - logMin) / LOG_N_BINS;

    for (int i = 0; i <= LOG_N_BINS; ++i) {
        edges[i] = std::pow(10.0, logMin + i * step);
    }

    return edges;
}

std::vector<Measurement> ReadSelectedCSV(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open CSV file: " + filename);
    }

    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("CSV file is empty: " + filename);
    }

    // UTF-8 BOM, if present.
    if (line.size() >= 3 &&
        static_cast<unsigned char>(line[0]) == 0xEF &&
        static_cast<unsigned char>(line[1]) == 0xBB &&
        static_cast<unsigned char>(line[2]) == 0xBF) {
        line.erase(0, 3);
    }

    const std::vector<std::string> header = SplitCSVLine(line);
    const int lumCol   = FindColumn(header, "luminosity");
    const int errCol   = FindColumn(header, "error");
    const int deltaCol = FindColumn(header, "deltaL");

    if (lumCol < 0 || errCol < 0 || deltaCol < 0) {
        throw std::runtime_error(
            "Required columns luminosity,error,deltaL not found in: " + filename
        );
    }

    const int maxRequired = std::max(lumCol, std::max(errCol, deltaCol));

    std::vector<Measurement> data;
    int lineNumber = 1;

    while (std::getline(input, line)) {
        ++lineNumber;
        line = Trim(line);
        if (line.empty()) continue;

        const std::vector<std::string> fields = SplitCSVLine(line);
        if (static_cast<int>(fields.size()) <= maxRequired) {
            std::cerr << "WARNING: incomplete line " << lineNumber
                      << " in " << filename << std::endl;
            continue;
        }

        try {
            Measurement m;
            m.luminosity = std::stod(fields[lumCol]);
            m.statError  = std::stod(fields[errCol]);
            m.systError  = std::stod(fields[deltaCol]);

            if (!std::isfinite(m.luminosity) ||
                !std::isfinite(m.statError) ||
                !std::isfinite(m.systError)) {
                throw std::runtime_error("non-finite value");
            }

            data.push_back(m);
        }
        catch (const std::exception&) {
            std::cerr << "WARNING: invalid numeric data at line "
                      << lineNumber << " in " << filename << std::endl;
        }
    }

    if (data.empty()) {
        throw std::runtime_error("No valid measurements in: " + filename);
    }

    return data;
}

TH1D* CreateLogHistogram(const std::string& name,
                         const std::vector<Measurement>& data)
{
    const std::vector<double> edges = MakeLogEdges();

    TH1D* h = new TH1D(
        name.c_str(), "", LOG_N_BINS, edges.data()
    );
    h->SetDirectory(nullptr);
    h->Sumw2();

    int below = 0;
    int above = 0;
    int nonPositive = 0;

    for (const Measurement& m : data) {
        if (m.luminosity <= 0.0) {
            ++nonPositive;
            continue;
        }
        if (m.luminosity < LOG_LUMINOSITY_MIN) {
            ++below;
            continue;
        }
        if (m.luminosity >= LOG_LUMINOSITY_MAX) {
            ++above;
            continue;
        }
        h->Fill(m.luminosity);
    }

    if (below || above || nonPositive) {
        std::cerr << "WARNING: histogram " << name
                  << " skipped " << nonPositive << " non-positive, "
                  << below << " below-range and "
                  << above << " above-range luminosity value(s)."
                  << std::endl;
    }

    if (h->GetEntries() <= 0) {
        delete h;
        throw std::runtime_error(
            "No entries fall inside the configured logarithmic histogram range."
        );
    }

    h->GetXaxis()->SetTitle("Luminosity");
    h->GetYaxis()->SetTitle("Number of hotspots");
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.040);
    h->GetYaxis()->SetLabelSize(0.040);
    h->GetXaxis()->SetTitleOffset(1.20);
    h->GetYaxis()->SetTitleOffset(1.25);

    return h;
}

TH1D* CreateNormalizedCopy(TH1D* source, const std::string& name)
{
    TH1D* out = static_cast<TH1D*>(source->Clone(name.c_str()));
    out->SetDirectory(nullptr);

    const double integral = out->Integral(1, out->GetNbinsX());
    if (integral <= 0.0) {
        delete out;
        throw std::runtime_error("Cannot normalize an empty histogram.");
    }

    out->Scale(1.0 / integral);
    return out;
}

void StyleHistogram(TH1D* h, Color_t color, Style_t marker)
{
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetFillColorAlpha(color, 0.25);
    h->SetFillStyle(1001);
    h->SetLineWidth(2);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.75);
}

void SaveBoth(TCanvas* canvas,
              const std::string& outputDir,
              const std::string& basename)
{
    const std::string safe = SafeName(basename);
    canvas->SaveAs((outputDir + "/" + safe + ".png").c_str());
    canvas->SaveAs((outputDir + "/" + safe + ".pdf").c_str());
}

// ============================================================================
// INDIVIDUAL DISTRIBUTION
// ============================================================================

void DrawIndividual(const std::string& csv,
                    const std::string& sensor,
                    const std::string& phase,
                    const std::string& outputDir)
{
    const std::vector<Measurement> data = ReadSelectedCSV(csv);

    std::vector<double> luminosities;
    for (const auto& m : data) luminosities.push_back(m.luminosity);

    std::sort(luminosities.begin(), luminosities.end());

    double median = 0.0;
    if (!luminosities.empty()) {
        const size_t n = luminosities.size();
        median = (n % 2 == 0)
            ? 0.5 * (luminosities[n/2 - 1] + luminosities[n/2])
            : luminosities[n/2];
    }

    TH1D* h = CreateLogHistogram("h_individual", data);
    StyleHistogram(h, kP6Blue, 20);

    TCanvas* c = new TCanvas(
        "c_individual", "Luminosity distribution", 1600, 1050
    );
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.09);
    c->SetTicks(1, 1);
    c->SetLogx();

    const std::string phaseLabel = PhaseLabel(phase);
    const std::string title =
        "Sensor " + sensor + " - " + phaseLabel +
        " luminosity distribution";

    const double ymax = std::max(1.0, h->GetMaximum());

    // The frame extends slightly beyond the histogram, but the histogram
    // itself still contains exactly 30 log bins from 1 to 1e4.
    TH1F* frame = c->DrawFrame(
        LOG_LUMINOSITY_MIN,
        0.0,
        DISPLAY_LUMINOSITY_MAX,
        1.35 * ymax
    );

    frame->SetTitle(title.c_str());
    frame->GetXaxis()->SetTitle("Luminosity");
    frame->GetYaxis()->SetTitle("Number of hotspots");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetTitleOffset(1.20);
    frame->GetYaxis()->SetTitleOffset(1.25);

    // HIST keeps the bin structure visible; E1 draws statistical
    // counting errors (sqrt(N)) from Sumw2().
    h->Draw("HIST SAME");
    h->Draw("E1 SAME");

    TLegend* legend = new TLegend(0.63, 0.70, 0.89, 0.87);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.033);
    legend->AddEntry(h, phaseLabel.c_str(), "lep");
    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("N = %.0f", h->GetEntries()),
        ""
    );
    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Mean = %.2f", h->GetMean()),
        ""
    );

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Median = %.2f", median),
        ""
    );

    legend->Draw();

    c->RedrawAxis();
    c->Modified();
    c->Update();

    SaveBoth(
        c,
        outputDir,
        sensor + "_" + phase + "_luminosity_distribution"
    );

    delete c;
    delete h;
}

// ============================================================================
// NORMALIZED COMPARISON
//
// The two samples are normalized independently to unit integral and overlaid.
// No bin-by-bin ratio is drawn: the datasets contain different hotspot
// populations, so such a ratio would not have a direct physical meaning.
// ============================================================================

void DrawNormalizedComparison(const std::string& csvReference,
                              const std::string& csvComparison,
                              const std::string& referenceLabel,
                              const std::string& comparisonLabel,
                              const std::string& title,
                              const std::string& basename,
                              const std::string& outputDir)
{
    const std::vector<Measurement> referenceData = ReadSelectedCSV(csvReference);
    const std::vector<Measurement> comparisonData = ReadSelectedCSV(csvComparison);

    const double meanRef   = MeanLuminosity(referenceData);
    const double medianRef = MedianLuminosity(referenceData);

    const double meanCmp   = MeanLuminosity(comparisonData);
    const double medianCmp = MedianLuminosity(comparisonData);

    TH1D* hRefRaw = CreateLogHistogram("h_reference_raw", referenceData);
    TH1D* hCmpRaw = CreateLogHistogram("h_comparison_raw", comparisonData);

    TH1D* hRef = CreateNormalizedCopy(hRefRaw, "h_reference_norm");
    TH1D* hCmp = CreateNormalizedCopy(hCmpRaw, "h_comparison_norm");

    StyleHistogram(hRef, kP6Blue, 20);
    StyleHistogram(hCmp, kP6Red, 21);

    TCanvas* c = new TCanvas(
        "c_comparison", "Normalized luminosity distributions", 1600, 1050
    );
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.09);
    c->SetTicks(1, 1);
    c->SetLogx();

    const double ymax = std::max(hRef->GetMaximum(), hCmp->GetMaximum());
    const double yDisplayMax = (ymax > 0.0 ? 1.30 * ymax : 1.0);

    // Same physical binning as the individual histograms; only the frame
    // extends to the right to leave empty graphical space.
    TH1F* frame = c->DrawFrame(
        LOG_LUMINOSITY_MIN,
        0.0,
        DISPLAY_LUMINOSITY_MAX,
        yDisplayMax
    );

    frame->SetTitle(title.c_str());
    frame->GetXaxis()->SetTitle("Luminosity");
    frame->GetYaxis()->SetTitle("Normalized number of hotspots");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetXaxis()->SetTitleOffset(1.20);
    frame->GetYaxis()->SetTitleOffset(1.25);

    hRef->Draw("HIST SAME");
    hCmp->Draw("HIST SAME");

    // Scale() propagates Sumw2 errors, so E1 is still the correct
    // statistical counting uncertainty after normalization.
    hRef->Draw("E1 SAME");
    hCmp->Draw("E1 SAME");

    TLegend* legend = new TLegend(0.57, 0.69, 0.89, 0.87);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.033);
    legend->AddEntry(
        hRef,
        Form("%s, N = %.0f", referenceLabel.c_str(), hRefRaw->GetEntries()),
        "lep"
    );

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Mean = %.2f, Median = %.2f", meanRef, medianRef),
        ""
    );

    legend->AddEntry(
        hCmp,
        Form("%s, N = %.0f", comparisonLabel.c_str(), hCmpRaw->GetEntries()),
        "lep"
    );

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Mean = %.2f, Median = %.2f", meanCmp, medianCmp),
        ""
    );

    legend->Draw();

    c->RedrawAxis();
    c->Modified();
    c->Update();

    SaveBoth(c, outputDir, basename);

    delete c;
    delete hRef;
    delete hCmp;
    delete hRefRaw;
    delete hCmpRaw;
}

// ============================================================================
// NORMALIZED BEFORE/AFTER COMPARISON WITH RATIO
//
// Upper pad: independently normalized luminosity distributions.
// Lower pad: normalized after-annealing distribution / normalized
//            before-annealing distribution, bin by bin.
// ============================================================================

void DrawNormalizedComparisonWithRatio(const std::string& csvBefore,
                                       const std::string& csvAfter,
                                       const std::string& afterLabel,
                                       const std::string& title,
                                       const std::string& basename,
                                       const std::string& outputDir)
{
    const std::vector<Measurement> beforeData = ReadSelectedCSV(csvBefore);
    const std::vector<Measurement> afterData  = ReadSelectedCSV(csvAfter);

    const double meanRef   = MeanLuminosity(beforeData);
    const double medianRef = MedianLuminosity(beforeData);

    const double meanCmp   = MeanLuminosity(afterData);
    const double medianCmp = MedianLuminosity(afterData);

    TH1D* hBeforeRaw = CreateLogHistogram("h_before_raw", beforeData);
    TH1D* hAfterRaw  = CreateLogHistogram("h_after_raw",  afterData);

    TH1D* hBefore = CreateNormalizedCopy(hBeforeRaw, "h_before_norm");
    TH1D* hAfter  = CreateNormalizedCopy(hAfterRaw,  "h_after_norm");

    StyleHistogram(hBefore, kP6Blue, 20);
    StyleHistogram(hAfter,  kP6Red, 21);

    // Ratio of the NORMALIZED distributions: after / before.
    TH1D* hRatio = static_cast<TH1D*>(hAfter->Clone("h_after_over_before"));
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hBefore);
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.70);
    hRatio->SetLineWidth(2);
    hRatio->SetFillStyle(0);

    TCanvas* c = new TCanvas(
        "c_phase_comparison", "Normalized luminosity distributions with ratio",
        1600, 1200
    );

    // ------------------------------------------------------------------------
    // Upper pad
    // ------------------------------------------------------------------------
    TPad* padUp = new TPad("padUp", "padUp", 0.0, 0.30, 1.0, 1.0);
    padUp->SetLeftMargin(0.12);
    padUp->SetRightMargin(0.05);
    padUp->SetBottomMargin(0.025);
    padUp->SetTopMargin(0.10);
    padUp->SetTicks(1, 1);
    padUp->SetLogx();
    padUp->Draw();
    padUp->cd();

    const double ymax = std::max(hBefore->GetMaximum(), hAfter->GetMaximum());
    const double yDisplayMax = (ymax > 0.0 ? 1.30 * ymax : 1.0);

    TH1F* frameUp = padUp->DrawFrame(
        LOG_LUMINOSITY_MIN,
        0.0,
        DISPLAY_LUMINOSITY_MAX,
        yDisplayMax
    );

    frameUp->SetTitle(title.c_str());
    frameUp->GetXaxis()->SetTitle("");
    frameUp->GetXaxis()->SetLabelSize(0.0);
    frameUp->GetYaxis()->SetTitle("Normalized number of hotspots");
    frameUp->GetYaxis()->SetTitleSize(0.055);
    frameUp->GetYaxis()->SetLabelSize(0.048);
    frameUp->GetYaxis()->SetTitleOffset(1.00);

    hBefore->Draw("HIST SAME");
    hAfter->Draw("HIST SAME");
    hBefore->Draw("E1 SAME");
    hAfter->Draw("E1 SAME");

    TLegend* legend = new TLegend(0.55, 0.62, 0.89, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.030);
    
    // QUI LA RIGA CORRETTA
    legend->AddEntry(
        hBefore,
        Form("Before annealing, N = %.0f", hBeforeRaw->GetEntries()),
        "lep"
    );

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Mean = %.2f, Median = %.2f", meanRef, medianRef),
        ""
    );

    legend->AddEntry(
        hAfter,
        Form("%s, N = %.0f", afterLabel.c_str(), hAfterRaw->GetEntries()),
        "lep"
    );

    legend->AddEntry(
        static_cast<TObject*>(nullptr),
        Form("Mean = %.2f, Median = %.2f", meanCmp, medianCmp),
        ""
    );

    legend->Draw();

    padUp->RedrawAxis();

    // ------------------------------------------------------------------------
    // Lower pad: after / before
    // ------------------------------------------------------------------------
    c->cd();
    TPad* padDown = new TPad("padDown", "padDown", 0.0, 0.0, 1.0, 0.30);
    padDown->SetLeftMargin(0.12);
    padDown->SetRightMargin(0.05);
    padDown->SetBottomMargin(0.34);
    padDown->SetTopMargin(0.035);
    padDown->SetTicks(1, 1);
    padDown->SetLogx();
    padDown->Draw();
    padDown->cd();

    // Determine a sensible ratio range from finite, non-zero-denominator bins.
    double ratioMax = 0.0;
    for (int ibin = 1; ibin <= hRatio->GetNbinsX(); ++ibin) {
        const double den = hBefore->GetBinContent(ibin);
        const double r   = hRatio->GetBinContent(ibin);
        if (den > 0.0 && std::isfinite(r)) {
            ratioMax = std::max(ratioMax, r + hRatio->GetBinError(ibin));
        }
    }

    double ratioMin = 0.0;
    for (int ibin =1; ibin<=hRatio->GetNbinsX(); ++ibin) {
        const double den = hBefore->GetBinContent(ibin);
        const double r   = hRatio->GetBinContent(ibin);
        if (den > 0.0 && std::isfinite(r)) {
            ratioMin = std::min(ratioMin, r - hRatio->GetBinError(ibin));
        }
    }

    double ratioYMax = 2.0;
    if (ratioMax > 0.0) {
        ratioYMax = std::max(2.0, 1.20 * ratioMax);
    }

    double ratioYMin = -1.0;
    if (ratioMin < 0.0) {
        ratioYMin = std::min(-1.0, 1.20 * ratioMin);
    }

    TH1F* frameDown = padDown->DrawFrame(
        LOG_LUMINOSITY_MIN,
        ratioYMin,
        DISPLAY_LUMINOSITY_MAX,
        ratioYMax
    );

    frameDown->SetTitle("");
    frameDown->GetXaxis()->SetTitle("Luminosity");
    frameDown->GetYaxis()->SetTitle("After / Before");

    frameDown->GetXaxis()->SetTitleSize(0.115);
    frameDown->GetXaxis()->SetLabelSize(0.095);
    frameDown->GetXaxis()->SetTitleOffset(1.15);

    frameDown->GetYaxis()->SetTitleSize(0.095);
    frameDown->GetYaxis()->SetLabelSize(0.080);
    frameDown->GetYaxis()->SetTitleOffset(0.55);
    frameDown->GetYaxis()->SetNdivisions(505);

    hRatio->Draw("E1 SAME");

    TLine* unity = new TLine(
        LOG_LUMINOSITY_MIN,
        1.0,
        DISPLAY_LUMINOSITY_MAX,
        1.0
    );
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->Draw("SAME");

    padDown->RedrawAxis();

    c->cd();
    c->Modified();
    c->Update();

    SaveBoth(c, outputDir, basename);

    delete c;
    delete hRatio;
    delete hBefore;
    delete hAfter;
    delete hBeforeRaw;
    delete hAfterRaw;
}

} // namespace

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================

void luminosity_distribution(const char* action,
                             const char* csv1,
                             const char* csv2 = "",
                             const char* sensor = "",
                             const char* phase = "",
                             const char* outputDir = "luminosity_distributions")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.042);
    gSystem->mkdir(outputDir, kTRUE);

    const std::string actionS = action ? action : "";
    const std::string csv1S   = csv1 ? csv1 : "";
    const std::string csv2S   = csv2 ? csv2 : "";
    const std::string sensorS = sensor ? sensor : "";
    const std::string phaseS  = phase ? phase : "";
    const std::string outS    = outputDir ? outputDir : "luminosity_distributions";

    try {
        if (actionS == "individual") {
            DrawIndividual(csv1S, sensorS, phaseS, outS);
            return;
        }

        if (actionS == "new_before") {
            DrawNormalizedComparison(
                csv1S,
                csv2S,
                "A1 new device",
                "A1 before annealing",
                "Sensor A1 - New device vs before annealing luminosity distributions",
                "A1_new_vs_before_annealing_luminosity_distribution",
                outS
            );
            return;
        }

        if (actionS == "phase_compare") {
            const std::string phaseLabel = PhaseLabel(phaseS);

            DrawNormalizedComparisonWithRatio(
                csv1S,
                csv2S,
                phaseLabel,
                "Sensor " + sensorS + " - Before annealing vs " + phaseLabel +
                    " luminosity distributions",
                sensorS + "_before_annealing_vs_" + phaseS +
                    "_luminosity_distribution",
                outS
            );
            return;
        }

        throw std::runtime_error(
            "Unknown action '" + actionS +
            "'. Use: individual, new_before, or phase_compare."
        );
    }
    catch (const std::exception& ex) {
        std::cerr << "\nERROR in luminosity_distribution.C: "
                  << ex.what() << std::endl;
        gSystem->Exit(1);
    }
}