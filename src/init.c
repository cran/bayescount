#include <R.h>
#include <Rinternals.h>

void poweranalysispopulation(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals,
                             double *coeffvarrep, double *coeffvarind, double *coeffvargroup, double *lowerl, double *upperl,
                             int *maxiterations, int *precision, double *lcil, double *ucil, int *print,
                             int *nin, int *ntotal);

void poweranalysispopulationfixed(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals,
                                  double *coeffvarrep, double *coeffvarind, double *coeffvargroup, int *maxiterations, int *print,
                                  double *meancounts);

void poweranalysissample(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals,
                         double *coeffvarrep, double *coeffvarind, double *coeffvargroup, double *lowerl, double *upperl,
                         int *maxiterations, int *precision, double *lcil, double *ucil, int *print,
                         int *nin, int *ntotal);

void poweranalysissamplefixed(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals,
                              double *coeffvarrep, double *coeffvarind, double *coeffvargroup, int *maxiterations, int *print,
                              double *meancounts);

void fecrtpowerpopulation(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates,
                          int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep,
                          double *postcoeffvarind, double *postcoeffvargroup, double *lowerl, double *upperl, int *maxiterations,
                          int *precision, double *lcil, double *ucil, int *print, int *nin,
                          int *ntotal);

void fecrtpowersample(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates,
                      int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep,
                      double *postcoeffvarind, double *postcoeffvargroup, double *lowerl, double *upperl, int *maxiterations,
                      int *precision, double *lcil, double *ucil, int *print, int *nin,
                      int *ntotal);

void fecrtpowerpopulationfixed(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates,
                               int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep,
                               double *postcoeffvarind, double *postcoeffvargroup, int *maxiterations, int *print, double *meanreds);

void fecrtpowersamplefixed(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates,
                           int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep,
                           double *postcoeffvarind, double *postcoeffvargroup, int *maxiterations, int *print, double *meanreds);


static const R_CallMethodDef CallEntries[] = {
  {"poweranalysispopulation", (DL_FUNC) &poweranalysispopulation, 17},
  {"poweranalysispopulationfixed", (DL_FUNC) &poweranalysispopulationfixed, 11},
  {"poweranalysissample", (DL_FUNC) &poweranalysissample, 17},
  {"poweranalysissamplefixed", (DL_FUNC) &poweranalysissamplefixed, 11},
  {"fecrtpowerpopulation", (DL_FUNC) &fecrtpowerpopulation, 21},
  {"fecrtpowersample", (DL_FUNC) &fecrtpowersample, 21},
  {"fecrtpowerpopulationfixed", (DL_FUNC) &fecrtpowerpopulationfixed, 15},
  {"fecrtpowersamplefixed", (DL_FUNC) &fecrtpowersamplefixed, 15},
  {NULL, NULL, 0}
};

void R_init_rjags(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


