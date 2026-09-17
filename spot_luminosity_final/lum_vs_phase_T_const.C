#include "phase_common.h"

// T is fixed; each phase is represented at the maximum available voltage.
// The macro preserves the original phase-plot canvas sizes and drawing style,
// while adding luminosity systematics and the requested ratio/B-evolution plots.
void lum_vs_phase_T_const(const char* csvfile,
                          const char* output_dir,
                          const char* prefix) {
    run_phase_analysis(csvfile, output_dir, prefix, true);
}
