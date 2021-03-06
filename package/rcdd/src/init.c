/*
d2q.c
lpcdd.c
lpcdd_f.c
mycddio.c
mycddio_f.c
nonred.c
q2d.c
q2q.c
qmatmult.c
qgram.c
qo.c
qoq.c
qsign.c
qsump.c
qminp.c
scdd.c
scdd_f.c
impliedLinearity.c
impliedLinearity_f.c
*/

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "rcdd.h"

static R_CMethodDef cMethods[] = {
    {NULL, NULL, 0, NULL}
};

static R_CallMethodDef callMethods[]  = {
    {"d2q", (DL_FUNC) &d2q, 1},
    {"q2d", (DL_FUNC) &q2d, 1},
    {"q2q", (DL_FUNC) &q2q, 1},
    {"qo", (DL_FUNC) &qo, 2},
    {"qoq", (DL_FUNC) &qoq, 3},
    {"qsign", (DL_FUNC) &qsign, 1},
    {"qsump", (DL_FUNC) &qsump, 2},
    {"qminp", (DL_FUNC) &qminp, 2},
    {"qmatmult", (DL_FUNC) &qmatmult, 2},
    {"qgram", (DL_FUNC) &qgram, 1},
    {"lpcdd", (DL_FUNC) &lpcdd, 4},
    {"lpcdd_f", (DL_FUNC) &lpcdd_f, 4},
    {"scdd", (DL_FUNC) &scdd, 7},
    {"scdd_f", (DL_FUNC) &scdd_f, 7},
    {"allfaces", (DL_FUNC) &allfaces, 1},
    {"allfaces_f", (DL_FUNC) &allfaces_f, 1},
    {"redundant", (DL_FUNC) &redundant, 2},
    {"redundant_f", (DL_FUNC) &redundant_f, 2},
    {"nonred", (DL_FUNC) &nonred, 2},
    {"test_my_subset", (DL_FUNC) &test_my_subset, 3},
    {"all_intersect", (DL_FUNC) &all_intersect, 2},
    {"all_union", (DL_FUNC) &all_union, 2},
    {"impliedLinearity", (DL_FUNC) &impliedLinearity, 2},
    {"impliedLinearity_f", (DL_FUNC) &impliedLinearity_f, 2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_rcdd(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

