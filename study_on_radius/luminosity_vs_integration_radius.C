// File: luminosity_vs_integration_radius.C
//
// Run with:
// root -l luminosity_vs_integration_radius.C

#include <TMultiGraph.h>
#include <TLegend.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>

#include <limits>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


// ============================================================
// STRUCTURES
// ============================================================

struct DataPoint {

    double integration_radius;
    double luminosity;
    double error;

    double x;
    double y;

    double T;
    double v;
};


struct ReferencePoint {

    double x;
    double y;

    double integration_radius;
};


struct DatasetInfo {

    string sensor;
    string fixed_condition;
    string phase;
};


// ============================================================
// STRING UTILITIES
// ============================================================

string trim(const string& input)
{
    const size_t first =
        input.find_first_not_of(" \t\r\n");

    if (first == string::npos)
        return "";

    const size_t last =
        input.find_last_not_of(" \t\r\n");

    return input.substr(
        first,
        last - first + 1
    );
}


vector<string> split_csv_line(
    const string& line
)
{
    vector<string> values;

    stringstream ss(line);
    string value;

    while (getline(ss, value, ',')) {
        values.push_back(trim(value));
    }

    return values;
}


// ============================================================
// GET BASENAME
// ============================================================

string get_basename(
    const string& path
)
{
    const size_t position =
        path.find_last_of("/\\");

    if (position == string::npos)
        return path;

    return path.substr(position + 1);
}


// ============================================================
// PARSE SENSOR, FIXED CONDITION AND PHASE FROM FILENAME
//
// Example:
//
// A1_T=20_before_annealing_all_radii.csv
//
// becomes:
//
// sensor          = A1
// fixed_condition = T=20
// phase           = before_annealing
//
// Also works for:
//
// A1_T=20_annealing_T=75_h=5_all_radii.csv
// ============================================================

bool parse_dataset_info(
    const string& csv_file,
    DatasetInfo& info
)
{
    string filename =
        get_basename(csv_file);

    const string suffix =
        "_all_radii.csv";

    if (
        filename.size() <= suffix.size() ||
        filename.substr(
            filename.size() - suffix.size()
        ) != suffix
    ) {

        cerr
            << "Error: input filename must end with "
            << suffix
            << endl;

        return false;
    }


    // Remove "_all_radii.csv"
    string stem =
        filename.substr(
            0,
            filename.size() - suffix.size()
        );


    // First underscore separates sensor
    const size_t first_separator =
        stem.find('_');

    if (first_separator == string::npos) {

        cerr
            << "Error: cannot extract sensor from "
            << filename
            << endl;

        return false;
    }


    // Second underscore separates fixed condition from phase
    const size_t second_separator =
        stem.find(
            '_',
            first_separator + 1
        );

    if (second_separator == string::npos) {

        cerr
            << "Error: cannot extract condition and phase from "
            << filename
            << endl;

        return false;
    }


    info.sensor =
        stem.substr(
            0,
            first_separator
        );


    info.fixed_condition =
        stem.substr(
            first_separator + 1,
            second_separator
            - first_separator
            - 1
        );


    info.phase =
        stem.substr(
            second_separator + 1
        );


    return true;
}


// ============================================================
// FIND COLUMN IN CSV HEADER
// ============================================================

int find_column(
    const vector<string>& header,
    const vector<string>& possible_names
)
{
    for (
        size_t i = 0;
        i < header.size();
        ++i
    ) {

        for (
            const string& name :
            possible_names
        ) {

            if (header[i] == name) {
                return static_cast<int>(i);
            }
        }
    }

    return -1;
}


// ============================================================
// FIND LUMINOSITY AT REFERENCE RADIUS
//
// If r_reference is present exactly in the scan:
//     use the measured point.
//
// Otherwise:
//     linearly interpolate between the two surrounding points.
//
// This allows the red circle to remain exactly at r_reference.
// ============================================================

bool get_luminosity_at_radius(
    const vector<DataPoint>& data,
    const double reference_radius,
    double& reference_luminosity
)
{
    if (data.empty())
        return false;


    const double tolerance = 1.0e-9;


    // --------------------------------------------------------
    // First look for an exact point
    // --------------------------------------------------------

    for (const DataPoint& point : data) {

        if (
            fabs(
                point.integration_radius
                - reference_radius
            ) < tolerance
        ) {

            reference_luminosity =
                point.luminosity;

            return true;
        }
    }


    // --------------------------------------------------------
    // Then interpolate
    // --------------------------------------------------------

    for (
        size_t i = 0;
        i + 1 < data.size();
        ++i
    ) {

        const DataPoint& left =
            data[i];

        const DataPoint& right =
            data[i + 1];


        if (
            reference_radius
                > left.integration_radius
            &&
            reference_radius
                < right.integration_radius
        ) {

            const double fraction =
                (
                    reference_radius
                    - left.integration_radius
                )
                /
                (
                    right.integration_radius
                    - left.integration_radius
                );


            reference_luminosity =
                left.luminosity
                +
                fraction
                * (
                    right.luminosity
                    - left.luminosity
                );


            return true;
        }
    }


    // Reference radius outside scan interval
    return false;
}


// ============================================================
// MAIN FUNCTION
// ============================================================

void luminosity_vs_integration_radius()
{
    // ========================================================
    // USER CONFIGURATION
    // ========================================================

    const string csv_file =
        "A1_T=20_before_annealing_luminosities/"
        "A1_T=20_before_annealing_all_radii.csv";


    // --------------------------------------------------------
    // Select operating condition
    //
    // Allowed values:
    //
    //     "v"
    //     "T"
    //
    // Example for A1 at T=20:
    //
    //     selected_variable = "v"
    //     selected_value    = 56.3
    //
    // --------------------------------------------------------

    const string selected_variable = "v";

    const double selected_value = 56.3;


    // Floating point comparison tolerance
    const double selection_tolerance =
        1.0e-6;


    // --------------------------------------------------------
    // Output directory
    //
    // Saving structure remains unchanged.
    // --------------------------------------------------------

    const string output_directory =
        "A1_T=20_bef_ann_plots";


    // ========================================================
    // CHECK SELECTED VARIABLE
    // ========================================================

    if (
        selected_variable != "v"
        &&
        selected_variable != "T"
    ) {

        cerr
            << "Error: selected_variable must be "
            << "\"v\" or \"T\"."
            << endl;

        return;
    }


    // ========================================================
    // EXTRACT DATASET INFORMATION FROM INPUT FILENAME
    // ========================================================

    DatasetInfo dataset_info;


    if (
        !parse_dataset_info(
            csv_file,
            dataset_info
        )
    ) {
        return;
    }


    const string dataset_tag =
        dataset_info.sensor
        + "_"
        + dataset_info.fixed_condition
        + "_"
        + dataset_info.phase;


    cout << endl;

    cout
        << "============================================"
        << endl;

    cout
        << "DATASET INFORMATION"
        << endl;

    cout
        << "============================================"
        << endl;

    cout
        << "Sensor:          "
        << dataset_info.sensor
        << endl;

    cout
        << "Fixed condition: "
        << dataset_info.fixed_condition
        << endl;

    cout
        << "Phase:           "
        << dataset_info.phase
        << endl;

    cout
        << "Selection:       "
        << selected_variable
        << " = "
        << selected_value
        << endl;


    // ========================================================
    // AUTOMATICALLY BUILD REFERENCE FILE NAME
    //
    // From:
    //
    // A1_T=20_before_annealing_all_radii.csv
    //
    // to:
    //
    // A1_T=20_before_annealing_reference_global_coordinates.csv
    // ========================================================

    const string all_radii_suffix =
        "_all_radii.csv";


    string reference_csv_file =
        csv_file.substr(
            0,
            csv_file.size()
            - all_radii_suffix.size()
        )
        +
        "_reference_global_coordinates.csv";


    cout
        << "Reference file:  "
        << reference_csv_file
        << endl;

    cout
        << "============================================"
        << endl
        << endl;


    // ========================================================
    // READ REFERENCE COORDINATES AND RADII
    // ========================================================

    ifstream reference_input(
        reference_csv_file
    );


    if (!reference_input.is_open()) {

        cerr
            << "Error: cannot open reference CSV file:"
            << endl
            << reference_csv_file
            << endl;

        return;
    }


    string line;


    // Read reference header
    getline(
        reference_input,
        line
    );


    vector<string> reference_header =
        split_csv_line(line);


    const int ref_id_column =
        find_column(
            reference_header,
            {
                "ID_hotspot_globale",
                "spot"
            }
        );


    const int ref_x_column =
        find_column(
            reference_header,
            {
                "x"
            }
        );


    const int ref_y_column =
        find_column(
            reference_header,
            {
                "y"
            }
        );


    const int ref_radius_column =
        find_column(
            reference_header,
            {
                "integration_radius",
                "radius"
            }
        );


    if (
        ref_id_column < 0
        ||
        ref_x_column < 0
        ||
        ref_y_column < 0
    ) {

        cerr
            << "Error: reference CSV must contain:"
            << endl
            << "  ID_hotspot_globale (or spot)"
            << endl
            << "  x"
            << endl
            << "  y"
            << endl;

        return;
    }


    if (ref_radius_column < 0) {

        cerr
            << endl
            << "ERROR:"
            << endl
            << "The reference CSV does not contain an "
            << "integration-radius column."
            << endl
            << endl
            << "Expected column name:"
            << endl
            << "    integration_radius"
            << endl
            << "or:"
            << endl
            << "    radius"
            << endl;

        return;
    }


    map<int, ReferencePoint>
        reference_by_spot;


    while (
        getline(
            reference_input,
            line
        )
    ) {

        if (line.empty())
            continue;


        vector<string> values =
            split_csv_line(line);


        const int maximum_column =
            max(
                {
                    ref_id_column,
                    ref_x_column,
                    ref_y_column,
                    ref_radius_column
                }
            );


        if (
            static_cast<int>(
                values.size()
            ) <= maximum_column
        ) {

            cerr
                << "Warning: invalid reference line skipped:"
                << endl
                << line
                << endl;

            continue;
        }


        try {

            const int spot =
                stoi(
                    values[
                        ref_id_column
                    ]
                );


            ReferencePoint reference;


            reference.x =
                stod(
                    values[
                        ref_x_column
                    ]
                );


            reference.y =
                stod(
                    values[
                        ref_y_column
                    ]
                );


            reference.integration_radius =
                stod(
                    values[
                        ref_radius_column
                    ]
                );


            reference_by_spot[spot] =
                reference;
        }

        catch (
            const exception& e
        ) {

            cerr
                << "Warning: invalid reference line skipped:"
                << endl
                << line
                << endl;
        }
    }


    reference_input.close();


    cout
        << "Reference hotspots loaded: "
        << reference_by_spot.size()
        << endl;


    // ========================================================
    // OPEN ALL-RADII CSV
    // ========================================================

    ifstream input(csv_file);


    if (!input.is_open()) {

        cerr
            << "Error: cannot open CSV file:"
            << endl
            << csv_file
            << endl;

        return;
    }


    // ========================================================
    // READ HEADER
    // ========================================================

    getline(
        input,
        line
    );


    vector<string> header =
        split_csv_line(line);


    const int id_column =
        find_column(
            header,
            {
                "ID_hotspot_globale",
                "spot"
            }
        );


    const int x_column = find_column(header, {"x"});


    const int y_column = find_column(header,{"y"});


    const int radius_column = find_column(header,{"integration_radius"});


    const int luminosity_column = find_column(header, {"luminosity"});


    const int error_column = find_column(header, {"error"});


    const int T_column = find_column(header,{"T"});


    const int v_column =find_column(header,{"v"});


    if (
        id_column < 0
        ||
        x_column < 0
        ||
        y_column < 0
        ||
        radius_column < 0
        ||
        luminosity_column < 0
        ||
        error_column < 0
        ||
        T_column < 0
        ||
        v_column < 0
    ) {

        cerr
            << "Error: input CSV must contain columns:"
            << endl
            << "ID_hotspot_globale,x,y,"
            << "integration_radius,"
            << "luminosity,error,T,v"
            << endl;

        return;
    }


    // ========================================================
    // READ DATA
    // ========================================================

    map<int, vector<DataPoint>>
        data_by_spot;


    set<double>
        available_selected_values;


    while (
        getline(
            input,
            line
        )
    ) {

        if (line.empty())
            continue;


        vector<string> values =
            split_csv_line(line);


        const int maximum_column =
            max(
                {
                    id_column,
                    x_column,
                    y_column,
                    radius_column,
                    luminosity_column,
                    error_column,
                    T_column,
                    v_column
                }
            );


        if (
            static_cast<int>(
                values.size()
            ) <= maximum_column
        ) {

            cerr
                << "Warning: invalid CSV line skipped:"
                << endl
                << line
                << endl;

            continue;
        }


        try {

            const int spot =
                stoi(
                    values[id_column]
                );


            DataPoint point;


            point.x =
                stod(
                    values[x_column]
                );


            point.y =
                stod(
                    values[y_column]
                );


            point.integration_radius =
                stod(
                    values[radius_column]
                );


            point.luminosity =
                stod(
                    values[luminosity_column]
                );


            point.error =
                stod(
                    values[error_column]
                );


            point.T =
                stod(
                    values[T_column]
                );


            point.v =
                stod(
                    values[v_column]
                );


            // ------------------------------------------------
            // Select operating condition
            // ------------------------------------------------

            const double condition_value =
                (
                    selected_variable == "v"
                )
                ?
                point.v
                :
                point.T;


            available_selected_values.insert(
                condition_value
            );


            if (
                fabs(
                    condition_value
                    - selected_value
                )
                <
                selection_tolerance
            ) {

                data_by_spot[spot]
                    .push_back(point);
            }
        }

        catch (
            const exception& e
        ) {

            cerr
                << "Warning: invalid CSV line skipped:"
                << endl
                << line
                << endl;
        }
    }


    input.close();


    // ========================================================
    // CHECK DATA
    // ========================================================

    if (
        data_by_spot.empty()
    ) {

        cerr
            << endl
            << "Error: no data found for "
            << selected_variable
            << " = "
            << selected_value
            << endl;


        cerr
            << endl
            << "Available values of "
            << selected_variable
            << ":"
            << endl;


        for (
            const double value :
            available_selected_values
        ) {

            cerr
                << "  "
                << value
                << endl;
        }


        return;
    }


    // ========================================================
    // CREATE OUTPUT DIRECTORY
    // ========================================================

    gSystem->mkdir(
        output_directory.c_str(),
        true
    );


    // ========================================================
    // SORT EVERY SPOT BY RADIUS
    // ========================================================

    for (
        auto& entry :
        data_by_spot
    ) {

        sort(
            entry.second.begin(),
            entry.second.end(),

            [](
                const DataPoint& a,
                const DataPoint& b
            ) {

                return
                    a.integration_radius
                    <
                    b.integration_radius;
            }
        );
    }


    // ========================================================
    // CREATE ONE CANVAS FOR EACH SPOT
    // ========================================================

    for (
        auto& entry :
        data_by_spot
    ) {

        const int spot_id =
            entry.first;


        vector<DataPoint>& spot_data =
            entry.second;


        const int number_of_points =
            static_cast<int>(
                spot_data.size()
            );


        vector<double>
            integration_radii(
                number_of_points
            );


        vector<double>
            luminosities(
                number_of_points
            );


        vector<double>
            luminosity_errors(
                number_of_points
            );


        vector<double>
            radius_errors(
                number_of_points,
                0.0
            );


        for (
            int i = 0;
            i < number_of_points;
            ++i
        ) {

            integration_radii[i] =
                spot_data[i]
                    .integration_radius;


            luminosities[i] =
                spot_data[i]
                    .luminosity;


            luminosity_errors[i] =
                spot_data[i]
                    .error;
        }


        const double spot_x =
            spot_data.front().x;


        const double spot_y =
            spot_data.front().y;


        // ====================================================
        // GRAPH
        // ====================================================

        TGraphErrors* graph =
            new TGraphErrors(
                number_of_points,
                integration_radii.data(),
                luminosities.data(),
                radius_errors.data(),
                luminosity_errors.data()
            );


        graph->SetName(
            Form(
                "graph_%s_spot_%d_%s_%.2f",
                dataset_tag.c_str(),
                spot_id,
                selected_variable.c_str(),
                selected_value
            )
        );


        graph->SetTitle(
            Form(
                "%s, %s, %s = %.2f, %s | "
                "Spot %d, coordinates (%.2f, %.2f);"
                "Integration radius [pixel];"
                "Luminosity",
                dataset_info.sensor.c_str(),
                dataset_info.fixed_condition.c_str(),
                selected_variable.c_str(),
                selected_value,
                dataset_info.phase.c_str(),
                spot_id,
                spot_x,
                spot_y
            )
        );


        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.4);
        graph->SetLineWidth(2);


        // ====================================================
        // CANVAS
        // ====================================================

        TCanvas* canvas =
            new TCanvas(
                Form(
                    "canvas_%s_spot_%d_%s_%.2f",
                    dataset_tag.c_str(),
                    spot_id,
                    selected_variable.c_str(),
                    selected_value
                ),
                Form(
                    "%s - Spot %d",
                    dataset_tag.c_str(),
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


        graph->GetXaxis()
            ->SetTitleSize(0.045);

        graph->GetYaxis()
            ->SetTitleSize(0.045);

        graph->GetXaxis()
            ->SetLabelSize(0.04);

        graph->GetYaxis()
            ->SetLabelSize(0.04);

        graph->GetXaxis()
            ->SetTitleOffset(1.2);

        graph->GetYaxis()
            ->SetTitleOffset(1.35);


        // ====================================================
        // REFERENCE RADIUS MARKER
        // ====================================================

        TGraph* reference_marker = nullptr;


        auto reference_it = reference_by_spot.find(spot_id);


        if (
            reference_it
            !=
            reference_by_spot.end()
        ) {

            const double reference_radius =
                reference_it
                    ->second
                    .integration_radius;


            double reference_luminosity;


            if (
                get_luminosity_at_radius(
                    spot_data,
                    reference_radius,
                    reference_luminosity
                )
            ) {

                reference_marker =
                    new TGraph(1);


                reference_marker
                    ->SetPoint(
                        0,
                        reference_radius,
                        reference_luminosity
                    );


                // Black cross
                reference_marker->SetMarkerStyle(47);

                reference_marker->SetMarkerColor(kBlack);

                reference_marker->SetLineColor(kRed);

                reference_marker->SetMarkerSize(2.2);

                reference_marker->Draw("P SAME");
            }

            else {

                cout
                    << "Warning: reference radius "
                    << reference_radius
                    << " for spot "
                    << spot_id
                    << " is outside scanned range."
                    << endl;
            }
        }


        canvas->Modified();
        canvas->Update();


        // ====================================================
        // SAVE
        // ====================================================

        string output_name =
            Form(
                "%s/%s_%s=%.2f_"
                "spot_%d_luminosity_vs_radius.png",
                output_directory.c_str(),
                dataset_tag.c_str(),
                selected_variable.c_str(),
                selected_value,
                spot_id
            );


        canvas->SaveAs(
            output_name.c_str()
        );


        cout
            << "Saved spot "
            << spot_id
            << ": "
            << output_name
            << endl;


        delete canvas;
        delete graph;

        if (
            reference_marker
            != nullptr
        ) {

            delete reference_marker;
        }
    }


    // ========================================================
    // GROUPED CANVASES
    // ========================================================

    struct LuminosityGroup {

        string title;
        string filename_label;

        double minimum;
        double maximum;

        vector<int> spots;
    };


    // Intervals:
    //
    // minimum <= Lmax < maximum

    vector<LuminosityGroup>
        luminosity_groups =
    {

        {
            "0 #leq L_{max} < 50",
            "0_50",
            0.0,
            50.0,
            {}
        },

        {
            "50 #leq L_{max} < 100",
            "50_100",
            50.0,
            100.0,
            {}
        },

        {
            "100 #leq L_{max} < 200",
            "100_200",
            100.0,
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
            numeric_limits<double>
                ::infinity(),
            {}
        }
    };


    // ========================================================
    // MAXIMUM LUMINOSITY BY SPOT
    // ========================================================

    map<int, double>
        maximum_luminosity_by_spot;


    // ========================================================
    // ASSIGN SPOTS TO GROUPS
    // ========================================================

    for (
        const auto& entry :
        data_by_spot
    ) {

        const int spot_id = entry.first;


        const vector<DataPoint>&
            spot_data =
            entry.second;


        double maximum_luminosity =
            -numeric_limits<double>
                ::infinity();


        for (
            const DataPoint& point :
            spot_data
        ) {

            maximum_luminosity =
                max(
                    maximum_luminosity,
                    point.luminosity
                );
        }


        maximum_luminosity_by_spot[
            spot_id
        ] =
            maximum_luminosity;


        bool group_found = false;


        for (
            LuminosityGroup& group :
            luminosity_groups
        ) {

            if (
                maximum_luminosity
                    >= group.minimum
                &&
                maximum_luminosity
                    < group.maximum
            ) {

                group.spots.push_back(
                    spot_id
                );


                group_found = true;

                break;
            }
        }


        if (!group_found) {

            cout
                << "Warning: spot "
                << spot_id
                << " with Lmax = "
                << maximum_luminosity
                << " is not in any group."
                << endl;
        }
    }


    // ========================================================
    // COLORS AND MARKERS
    // ========================================================

    const int graph_colors[6] =
    {
        kRed + 1,
        kBlue + 1,
        kGreen + 2,
        kMagenta + 1,
        kOrange + 7,
        kAzure + 2
    };


    const int marker_styles[6] =
    {
        20,
        21,
        22,
        23,
        24,
        25
    };


    const size_t spots_per_canvas = 6;


    // ========================================================
    // CREATE GROUPED CANVASES
    // ========================================================

    for (
        LuminosityGroup& group :
        luminosity_groups
    ) {

        if (
            group.spots.empty()
        ) {

            cout
                << "No spots found in luminosity group "
                << group.title
                << endl;

            continue;
        }


        // ----------------------------------------------------
        // Sort spots by maximum luminosity
        // ----------------------------------------------------

        sort(
            group.spots.begin(),
            group.spots.end(),

            [&maximum_luminosity_by_spot](
                int spot_a,
                int spot_b
            ) {

                return
                    maximum_luminosity_by_spot
                        .at(spot_a)
                    <
                    maximum_luminosity_by_spot
                        .at(spot_b);
            }
        );


        const size_t number_of_canvases =
            (
                group.spots.size()
                +
                spots_per_canvas
                -
                1
            )
            /
            spots_per_canvas;


        for (
            size_t canvas_index = 0;
            canvas_index
                <
                number_of_canvases;
            ++canvas_index
        ) {

            const size_t first_spot_index =
                canvas_index
                *
                spots_per_canvas;


            const size_t last_spot_index =
                min(
                    first_spot_index
                    +
                    spots_per_canvas,

                    group.spots.size()
                );


            // =================================================
            // CANVAS
            // =================================================

            TCanvas* grouped_canvas =
                new TCanvas(
                    Form(
                        "grouped_canvas_%s_%s_%zu",
                        dataset_tag.c_str(),
                        group.filename_label.c_str(),
                        canvas_index + 1
                    ),
                    Form(
                        "%s - %s",
                        dataset_tag.c_str(),
                        group.title.c_str()
                    ),
                    1800,
                    1200
                );


            grouped_canvas
                ->SetLeftMargin(0.12);

            grouped_canvas
                ->SetRightMargin(0.06);

            grouped_canvas
                ->SetBottomMargin(0.12);

            grouped_canvas
                ->SetTopMargin(0.11);

            grouped_canvas
                ->SetGrid();


            TMultiGraph* multigraph =
                new TMultiGraph();


            multigraph->SetTitle(
                Form(
                    "%s, %s, %s = %.2f, %s | "
                    "%s, canvas %zu;"
                    "Integration radius [pixel];"
                    "Luminosity",
                    dataset_info.sensor.c_str(),
                    dataset_info.fixed_condition.c_str(),
                    selected_variable.c_str(),
                    selected_value,
                    dataset_info.phase.c_str(),
                    group.title.c_str(),
                    canvas_index + 1
                )
            );


            TLegend* legend =
                new TLegend(
                    0.15,
                    0.65,
                    0.92,
                    0.88
                );


            legend->SetBorderSize(1);
            legend->SetFillStyle(1001);
            legend->SetTextSize(0.025);
            legend->SetNColumns(2);


            // Keep objects alive until canvas is saved
            vector<TGraphErrors*>
                grouped_graphs;


            vector<TGraph*>
                reference_markers;


            TGraph* first_reference_marker =
                nullptr;


            // =================================================
            // ADD SPOTS TO GROUPED CANVAS
            // =================================================

            for (
                size_t spot_index =
                    first_spot_index;

                spot_index
                    <
                    last_spot_index;

                ++spot_index
            ) {

                const int spot_id =
                    group.spots[
                        spot_index
                    ];


                vector<DataPoint>&
                    spot_data =
                    data_by_spot[
                        spot_id
                    ];


                const int number_of_points =
                    static_cast<int>(
                        spot_data.size()
                    );


                vector<double>
                    integration_radii(
                        number_of_points
                    );


                vector<double>
                    luminosities(
                        number_of_points
                    );


                vector<double>
                    luminosity_errors(
                        number_of_points
                    );


                vector<double>
                    radius_errors(
                        number_of_points,
                        0.0
                    );


                for (
                    int i = 0;
                    i < number_of_points;
                    ++i
                ) {

                    integration_radii[i] =
                        spot_data[i]
                            .integration_radius;


                    luminosities[i] =
                        spot_data[i]
                            .luminosity;


                    luminosity_errors[i] =
                        spot_data[i]
                            .error;
                }


                TGraphErrors* graph =
                    new TGraphErrors(
                        number_of_points,
                        integration_radii.data(),
                        luminosities.data(),
                        radius_errors.data(),
                        luminosity_errors.data()
                    );


                const int local_graph_index =
                    static_cast<int>(
                        spot_index
                        -
                        first_spot_index
                    );


                graph->SetName(
                    Form(
                        "group_%s_canvas_%zu_spot_%d",
                        group.filename_label.c_str(),
                        canvas_index + 1,
                        spot_id
                    )
                );


                graph->SetMarkerStyle(
                    marker_styles[
                        local_graph_index
                    ]
                );


                graph->SetMarkerColor(
                    graph_colors[
                        local_graph_index
                    ]
                );


                graph->SetLineColor(
                    graph_colors[
                        local_graph_index
                    ]
                );


                graph->SetMarkerSize(1.1);
                graph->SetLineWidth(2);


                multigraph->Add(
                    graph,
                    "LP"
                );


                grouped_graphs.push_back(
                    graph
                );


                const double spot_x =
                    spot_data.front().x;


                const double spot_y =
                    spot_data.front().y;


                // =============================================
                // LOOK UP REFERENCE RADIUS
                // =============================================

                auto reference_it =
                    reference_by_spot.find(
                        spot_id
                    );


                if (
                    reference_it
                    !=
                    reference_by_spot.end()
                ) {

                    const double reference_radius =
                        reference_it
                            ->second
                            .integration_radius;


                    double reference_luminosity;


                    if (
                        get_luminosity_at_radius(
                            spot_data,
                            reference_radius,
                            reference_luminosity
                        )
                    ) {

                        TGraph* marker =
                            new TGraph(1);


                        marker->SetPoint(
                            0,
                            reference_radius,
                            reference_luminosity
                        );


                        // Black cross
                        marker->SetMarkerStyle(35);
                        marker->SetMarkerColor(kBlack);
                        marker->SetLineColor(kRed);
                        marker->SetMarkerSize(2.2);


                        reference_markers
                            .push_back(
                                marker
                            );


                        if (
                            first_reference_marker
                            ==
                            nullptr
                        ) {

                            first_reference_marker =
                                marker;
                        }


                        legend->AddEntry(
                            graph,
                            Form(
                                "Spot %d: "
                                "(%.1f, %.1f), "
                                "L_{max}=%.1f, "
                                "r_{ref}=%.1f",
                                spot_id,
                                spot_x,
                                spot_y,
                                maximum_luminosity_by_spot[
                                    spot_id
                                ],
                                reference_radius
                            ),
                            "lp"
                        );
                    }

                    else {

                        legend->AddEntry(
                            graph,
                            Form(
                                "Spot %d: "
                                "(%.1f, %.1f), "
                                "L_{max}=%.1f",
                                spot_id,
                                spot_x,
                                spot_y,
                                maximum_luminosity_by_spot[
                                    spot_id
                                ]
                            ),
                            "lp"
                        );


                        cout
                            << "Warning: reference radius "
                            << reference_radius
                            << " for spot "
                            << spot_id
                            << " is outside scan."
                            << endl;
                    }
                }

                else {

                    legend->AddEntry(
                        graph,
                        Form(
                            "Spot %d: "
                            "(%.1f, %.1f), "
                            "L_{max}=%.1f",
                            spot_id,
                            spot_x,
                            spot_y,
                            maximum_luminosity_by_spot[
                                spot_id
                            ]
                        ),
                        "lp"
                    );


                    cout
                        << "Warning: no reference radius "
                        << "for spot "
                        << spot_id
                        << endl;
                }
            }


            // =================================================
            // DRAW MULTIGRAPH
            // =================================================

            multigraph->Draw("A");


            double current_ymax =
                multigraph
                    ->GetYaxis()
                    ->GetXmax();


            double current_ymin =
                multigraph
                    ->GetYaxis()
                    ->GetXmin();


            multigraph
                ->GetYaxis()
                ->SetRangeUser(
                    current_ymin,
                    1.70 * current_ymax
                );


            multigraph
                ->GetXaxis()
                ->SetTitleSize(0.045);

            multigraph
                ->GetYaxis()
                ->SetTitleSize(0.045);


            multigraph
                ->GetXaxis()
                ->SetLabelSize(0.04);

            multigraph
                ->GetYaxis()
                ->SetLabelSize(0.04);


            multigraph
                ->GetXaxis()
                ->SetTitleOffset(1.15);

            multigraph
                ->GetYaxis()
                ->SetTitleOffset(1.25);


            // =================================================
            // DRAW RED REFERENCE CIRCLES
            // =================================================

            for (
                TGraph* marker :
                reference_markers
            ) {

                marker->Draw(
                    "P SAME"
                );
            }


            // One explanation in legend
            if (
                first_reference_marker
                !=
                nullptr
            ) {

                legend->AddEntry(
                    first_reference_marker,
                    "Reference integration radius",
                    "p"
                );
            }


            legend->Draw();


            grouped_canvas->Modified();
            grouped_canvas->Update();


            // =================================================
            // SAVE GROUPED CANVAS
            // =================================================

            string grouped_output_name =
                Form(
                    "%s/%s_%s=%.2f_"
                    "group_%s_canvas_%zu.png",
                    output_directory.c_str(),
                    dataset_tag.c_str(),
                    selected_variable.c_str(),
                    selected_value,
                    group.filename_label.c_str(),
                    canvas_index + 1
                );


            grouped_canvas->SaveAs(
                grouped_output_name.c_str()
            );


            cout
                << "Saved grouped canvas: "
                << grouped_output_name
                << endl;


        
        }
    }


    // ========================================================
    // FINAL REPORT
    // ========================================================

    cout
        << endl
        << "============================================"
        << endl;

    cout
        << "Created plots for "
        << data_by_spot.size()
        << " spots"
        << endl;

    cout
        << "Sensor:    "
        << dataset_info.sensor
        << endl;

    cout
        << "Condition: "
        << dataset_info.fixed_condition
        << endl;

    cout
        << "Phase:     "
        << dataset_info.phase
        << endl;

    cout
        << "Selected:  "
        << selected_variable
        << " = "
        << selected_value
        << endl;

    cout
        << "============================================"
        << endl;
}