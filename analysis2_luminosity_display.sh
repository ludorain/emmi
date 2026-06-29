cd spot_luminosity
root -l 'lum_vs_T_two_v.C("run=202508_2overvoltages_scanTemp.csv")'
root -l 'lum_vs_V.C("run=202508_scanV.csv")'


#Device irradiati
cd spot_luminosity
root -l 'lum_vs_V.C("A1_T=20_run=20260513-032137_total.csv")'