// File: luminosity_vs_integration_radius.C
// root -l luminosity_vs_integration_radius.C
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TColor.h>

#include <limits>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct DataPoint {
    double integration_radius;
    double luminosity;
    double error;
    double x;
    double y;
};

void luminosity_vs_integration_radius()
{
    // ============================================================
    // USER CONFIGURATION
    // ============================================================
    const string csv_file = "bef_ann_A1_T=20_radius_study.csv";

    const double selected_v_fin = 5.0;

    // Tolerance used when comparing floating-point v_fin values
    const double v_fin_tolerance = 1.0e-6;

    // Output directory
    const string output_directory = "luminosity_vs_radius_plots";

    // ============================================================
    // OPEN INPUT CSV FILE
    // ============================================================
    ifstream input(csv_file);

    if (!input.is_open()) {
        cerr << "Error: cannot open the CSV file: "
             << csv_file << endl;
        return;
    }

    string line;

    // Skip the header:
    // spot,x,y,luminosity,error,T,v,v_fin,integration_radius
    getline(input, line);

    map<int, vector<DataPoint>> data_by_spot;

    // ============================================================
    // READ CSV FILE
    // ============================================================
    while (getline(input, line)) {

        if (line.empty()) {
            continue;
        }

        stringstream ss(line);
        string value;

        int spot;
        double x;
        double y;
        double luminosity;
        double error;
        double T;
        double v;
        double v_fin;
        double integration_radius;

        try {
            getline(ss, value, ',');
            spot = stoi(value);

            getline(ss, value, ',');
            x = stod(value);

            getline(ss, value, ',');
            y = stod(value);

            getline(ss, value, ',');
            luminosity = stod(value);

            getline(ss, value, ',');
            error = stod(value);

            getline(ss, value, ',');
            T = stod(value);

            getline(ss, value, ',');
            v = stod(value);

            getline(ss, value, ',');
            v_fin = stod(value);

            getline(ss, value, ',');
            integration_radius = stod(value);
        }
        catch (const exception& e) {
            cerr << "Warning: invalid CSV line skipped:"
                 << endl
                 << line << endl;
            continue;
        }

        // Select only the requested spot and v_fin
        if (fabs(v_fin - selected_v_fin) < v_fin_tolerance) {
        data_by_spot[spot].push_back({integration_radius, luminosity, error, x, y});
        }
    }

    input.close();
// ============================================================
// CHECK DATA
// ============================================================
if (data_by_spot.empty()) {
    cerr << "Error: no data found for v_fin = "
         << selected_v_fin << endl;
    return;
}

// Create output directory
gSystem->mkdir(output_directory.c_str(), true);


/*
// ============================================================
// CREATE ONE CANVAS FOR EACH SPOT
// ============================================================
for (auto& entry : data_by_spot) {

    int spot_id = entry.first;
    vector<DataPoint>& spot_data = entry.second;

    // Sort by integration radius
    sort(
        spot_data.begin(),
        spot_data.end(),
        [](const DataPoint& a, const DataPoint& b) {
            return a.integration_radius < b.integration_radius;
        }
    );

    int number_of_points =
        static_cast<int>(spot_data.size());

    vector<double> integration_radii(number_of_points);
    vector<double> luminosities(number_of_points);
    vector<double> luminosity_errors(number_of_points);
    vector<double> radius_errors(number_of_points, 0.0);

    for (int i = 0; i < number_of_points; ++i) {
        integration_radii[i] =
            spot_data[i].integration_radius;

        luminosities[i] =
            spot_data[i].luminosity;

        luminosity_errors[i] =
            spot_data[i].error;
    }

    double spot_x = spot_data.front().x;
    double spot_y = spot_data.front().y;

    // ========================================================
    // GRAPH
    // ========================================================
    TGraphErrors* graph = new TGraphErrors(
        number_of_points,
        integration_radii.data(),
        luminosities.data(),
        radius_errors.data(),
        luminosity_errors.data()
    );

    graph->SetName(
        Form(
            "graph_spot_%d_vfin_%.2f",
            spot_id,
            selected_v_fin
        )
    );

    graph->SetTitle(
        Form(
            "Spot %d, coordinates (%.2f, %.2f), v_{fin} = %.2f;"
            "Integration radius [pixel];"
            "Luminosity",
            spot_id,
            spot_x,
            spot_y,
            selected_v_fin
        )
    );

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.4);
    graph->SetLineWidth(2);

    // ========================================================
    // CANVAS
    // ========================================================
    TCanvas* canvas = new TCanvas(
        Form(
            "canvas_spot_%d_vfin_%.2f",
            spot_id,
            selected_v_fin
        ),
        Form(
            "Spot %d - Luminosity versus integration radius",
            spot_id
        ),
        1000,
        800
    );

    canvas->SetLeftMargin(0.13);
    canvas->SetRightMargin(0.05);
    canvas->SetBottomMargin(0.13);
    canvas->SetTopMargin(0.11);
    canvas->SetGrid();

    graph->Draw("APL");

    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);

    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetLabelSize(0.04);

    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()->SetTitleOffset(1.35);

    canvas->Modified();
    canvas->Update();

    // ========================================================
    // SAVE
    // ========================================================
    string output_name = Form(
        "%s/spot_%d_vfin_%.2f_luminosity_vs_radius.png",
        output_directory.c_str(),
        spot_id,
        selected_v_fin
    );

    canvas->SaveAs(output_name.c_str());

    cout << "Saved spot "
         << spot_id
         << " with "
         << number_of_points
         << " points: "
         << output_name
         << endl;

    delete canvas;
    delete graph;
        }

    cout << "Created plots for "
     << data_by_spot.size()
     << " spots at v_fin = "
     << selected_v_fin
     << endl;

*/

// ============================================================
// GROUPED CANVASES: MAXIMUM 10 SPOTS PER CANVAS
// ============================================================

// Information defining a luminosity group
struct LuminosityGroup {
    string title;
    string filename_label;
    double minimum;
    double maximum;
    vector<int> spots;
};

// The intervals are interpreted as:
// minimum <= maximum luminosity < maximum
vector<LuminosityGroup> luminosity_groups = {

    {
        "0 #leq L_{max} < 50",
        "0_50",
        0.0,
        50.0,
        {}
    },

    {
        "50 #leq L_{max} < 200",
        "50_200",
        50.0,
        200.0,
        {}
    },

    {
        "200 #leq L_{max} < 300",
        "200_300",
        200.0,
        300.0,
        {}
    },

    {
        "300 #leq L_{max} < 500",
        "300_500",
        300.0,
        500.0,
        {}
    },

    {
        "500 #leq L_{max} < 1000",
        "500_1000",
        500.0,
        1000.0,
        {}
    },

    {
        "L_{max} #geq 1000",
        "more_than_1000",
        1000.0,
        numeric_limits<double>::infinity(),
        {}
    }
};

// Store the maximum luminosity of each spot
map<int, double> maximum_luminosity_by_spot;

// ============================================================
// ASSIGN EACH SPOT TO A LUMINOSITY GROUP
// ============================================================
for (const auto& entry : data_by_spot) {

    int spot_id = entry.first;
    const vector<DataPoint>& spot_data = entry.second;

    double maximum_luminosity =
        -numeric_limits<double>::infinity();

    for (const DataPoint& point : spot_data) {
        maximum_luminosity = max(
            maximum_luminosity,
            point.luminosity
        );
    }

    maximum_luminosity_by_spot[spot_id] =
        maximum_luminosity;

    bool group_found = false;

    for (LuminosityGroup& group : luminosity_groups) {

        if (
            maximum_luminosity >= group.minimum &&
            maximum_luminosity < group.maximum
        ) {
            group.spots.push_back(spot_id);
            group_found = true;
            break;
        }
    }

    if (!group_found) {
        cout << "Warning: spot "
             << spot_id
             << " has maximum luminosity "
             << maximum_luminosity
             << " and does not belong to any defined group."
             << endl;
    }
}

// ============================================================
// COLORS AND MARKER STYLES
// ============================================================
const int graph_colors[10] = {
    kBlack,
    kRed + 1,
    kBlue + 1,
    kGreen + 2,
    kMagenta + 1,
    kOrange + 7,
    kCyan + 2,
    kViolet + 1,
    kAzure + 2,
    kPink + 7
};

const int marker_styles[10] = {
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    30
};

const size_t spots_per_canvas = 10;

// ============================================================
// CREATE THE GROUPED CANVASES
// ============================================================
for (LuminosityGroup& group : luminosity_groups) {

    if (group.spots.empty()) {
        cout << "No spots found in luminosity group "
             << group.title << endl;
        continue;
    }

    // Sort spots by maximum luminosity
    sort(
        group.spots.begin(),
        group.spots.end(),
        [&maximum_luminosity_by_spot](int spot_a, int spot_b) {
            return maximum_luminosity_by_spot.at(spot_a)
                 < maximum_luminosity_by_spot.at(spot_b);
        }
    );

    const size_t number_of_canvases =
        (group.spots.size() + spots_per_canvas - 1)
        / spots_per_canvas;

    for (
        size_t canvas_index = 0;
        canvas_index < number_of_canvases;
        ++canvas_index
    ) {
        const size_t first_spot_index =
            canvas_index * spots_per_canvas;

        const size_t last_spot_index = min(
            first_spot_index + spots_per_canvas,
            group.spots.size()
        );

        TCanvas* grouped_canvas = new TCanvas(
            Form(
                "grouped_canvas_%s_%zu",
                group.filename_label.c_str(),
                canvas_index + 1
            ),
            Form(
                "Luminosity group %s",
                group.title.c_str()
            ),
            1200,
            1200
        );

        grouped_canvas->SetLeftMargin(0.12);
        grouped_canvas->SetRightMargin(0.06);
        grouped_canvas->SetBottomMargin(0.12);
        grouped_canvas->SetTopMargin(0.11);
        grouped_canvas->SetGrid();

        TMultiGraph* multigraph = new TMultiGraph();

        multigraph->SetTitle(
            Form(
                "%s, v_{fin} = %.2f, canvas %zu;"
                "Integration radius [pixel];"
                "Luminosity",
                group.title.c_str(),
                selected_v_fin,
                canvas_index + 1
            )
        );

        TLegend* legend = new TLegend(
            0.61,
            0.55,
            0.90,
            0.88
        );

        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->SetTextSize(0.025);

        // ====================================================
        // ADD UP TO 10 SPOTS TO THIS CANVAS
        // ====================================================
        for (
            size_t spot_index = first_spot_index;
            spot_index < last_spot_index;
            ++spot_index
        ) {
            int spot_id = group.spots[spot_index];

            vector<DataPoint>& spot_data =
                data_by_spot[spot_id];

            // Sort again for safety
            sort(
                spot_data.begin(),
                spot_data.end(),
                [](const DataPoint& a, const DataPoint& b) {
                    return a.integration_radius
                         < b.integration_radius;
                }
            );

            const int number_of_points =
                static_cast<int>(spot_data.size());

            vector<double> integration_radii(
                number_of_points
            );

            vector<double> luminosities(
                number_of_points
            );

            vector<double> luminosity_errors(
                number_of_points
            );

            vector<double> radius_errors(
                number_of_points,
                0.0
            );

            for (int i = 0; i < number_of_points; ++i) {

                integration_radii[i] =
                    spot_data[i].integration_radius;

                luminosities[i] =
                    spot_data[i].luminosity;

                luminosity_errors[i] =
                    spot_data[i].error;
            }

            TGraphErrors* graph = new TGraphErrors(
                number_of_points,
                integration_radii.data(),
                luminosities.data(),
                radius_errors.data(),
                luminosity_errors.data()
            );

            const int local_graph_index =
                static_cast<int>(
                    spot_index - first_spot_index
                );

            graph->SetName(
                Form(
                    "group_graph_%s_canvas_%zu_spot_%d",
                    group.filename_label.c_str(),
                    canvas_index + 1,
                    spot_id
                )
            );

            graph->SetMarkerStyle(
                marker_styles[local_graph_index]
            );

            graph->SetMarkerColor(
                graph_colors[local_graph_index]
            );

            graph->SetLineColor(
                graph_colors[local_graph_index]
            );

            graph->SetMarkerSize(1.1);
            graph->SetLineWidth(2);

            multigraph->Add(graph, "LP");

            const double spot_x =
                spot_data.front().x;

            const double spot_y =
                spot_data.front().y;

            legend->AddEntry(
                graph,
                Form(
                    "Spot %d: (%.1f, %.1f), L_{max}=%.1f",
                    spot_id,
                    spot_x,
                    spot_y,
                    maximum_luminosity_by_spot[spot_id]
                ),
                "lp"
            );
        }

        multigraph->Draw("A");


        double current_ymax = multigraph->GetYaxis()->GetXmax();
        double current_ymin = multigraph->GetYaxis()->GetXmin();
        multigraph->GetYaxis()->SetRangeUser(current_ymin, 1.70 * current_ymax);


        multigraph->GetXaxis()->SetTitleSize(0.045);
        multigraph->GetYaxis()->SetTitleSize(0.045);

        multigraph->GetXaxis()->SetLabelSize(0.04);
        multigraph->GetYaxis()->SetLabelSize(0.04);

        multigraph->GetXaxis()->SetTitleOffset(1.15);
        multigraph->GetYaxis()->SetTitleOffset(1.25);

        legend->Draw();

        grouped_canvas->Modified();
        grouped_canvas->Update();

        // ====================================================
        // SAVE THE GROUPED CANVAS
        // ====================================================
        string grouped_output_name = Form(
            "%s/group_%s_vfin_%.2f_canvas_%zu.png",
            output_directory.c_str(),
            group.filename_label.c_str(),
            selected_v_fin,
            canvas_index + 1
        );

        grouped_canvas->SaveAs(
            grouped_output_name.c_str()
        );

        cout << "Saved grouped canvas: "
             << grouped_output_name
             << endl;
    }
}
}