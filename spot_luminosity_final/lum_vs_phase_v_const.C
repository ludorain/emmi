#include "phase_common.h"

// Overvoltage is fixed; each phase is represented at the maximum available T.
// The macro preserves the original phase-plot canvas sizes and drawing style,
// while adding luminosity systematics and the requested ratio/lambda-evolution plots.
void lum_vs_phase_v_const(const char* csvfile,
                          const char* output_dir,
                          const char* prefix) {
    run_phase_analysis(csvfile, output_dir, prefix, false);
}
