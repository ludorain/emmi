#include "phase_common.h"
void lum_vs_phase_T_const(const char* csvfile,
                          const char* output_dir,
                          const char* prefix,
                          const char* csv_R16="",
                          const char* csv_R24="") {
    run_phase_analysis(csvfile, output_dir, prefix, true, csv_R16, csv_R24);
}
