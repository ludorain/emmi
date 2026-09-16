5analysis pipeline package
==========================

Installation layout expected inside emmi/:

emmi/
  5analysis.sh
  DATA_irradiated_isolated_R=16/
  DATA_irradiated_isolated_R=20/
  DATA_irradiated_isolated_R=24/
  spot_luminosity_final/
    analysis_common.h
    phase_common.h
    lum_vs_v_systematics.C
    lum_vs_T_systematics.C
    lum_vs_v_fit.C
    lum_vs_T_fit.C
    lum_vs_phase_T_const.C
    lum_vs_phase_v_const.C

Usage
-----
./5analysis.sh A1 T analysis
./5analysis.sh A1 v phases
./5analysis.sh A1 analysis      # both T and v
./5analysis.sh A1 T             # analysis + phases
./5analysis.sh A1               # both constants + both commands

No arguments -> error by design.

Analysis outputs
----------------
spot_luminosity_final/anal/<sensor_condition>/<phase>/

For T constant:
  <prefix>_B_values.csv                 R=16/R=24 fit values
  <prefix>_B_fit_results.csv            nominal R=20 B + statistical/systematic errors
  <prefix>_analysis.csv                 selected phase rows with appended fit results
  <prefix>_B_vs_spot.png
  lum_vs_v_all_spots/*.png

For v constant:
  analogous lambda outputs.

Phase outputs
-------------
spot_luminosity_final/phases/<sensor_condition>/
  single_spots/
  grouped/
  ratios/
  fit_parameter_vs_phase/

Important choices
-----------------
1. R=20 is the nominal fit. R=16/R=24 are used only for the fit-parameter systematic.
2. deltaL_plus/deltaL_minus are displayed as luminosity systematic uncertainties and are NOT used as fit weights.
3. Parameter systematics are:
   delta1 = abs(parameter_R24 - parameter_R20)
   delta2 = abs(parameter_R20 - parameter_R16)
   syst_max = max(delta1, delta2)
   syst_min = min(delta1, delta2)
4. Parameter-vs-spot canvases preserve the original chi2/ndf < 2 selection and additionally require a converged fit.
5. Phase evolution at T constant uses the maximum available V point in each phase.
   Phase evolution at v constant uses the maximum available T point in each phase.
6. Missing annealing phases are allowed. If R=16 or R=24 is missing, the nominal R=20 result is still produced but no parameter systematic is assigned for that hotspot.
7. If several run files exist for one sensor/constant/phase, the lexicographically latest timestamped run is selected.
8. B-vs-phase requires v_fin in the all-phases CSV. If v_fin is absent, those canvases are skipped rather than fitting B against absolute bias voltage.

Validation
----------
The Bash script was syntax-checked with `bash -n`. ROOT is not installed in the execution environment used to build this package, so the ROOT macros could not be runtime-compiled here.
