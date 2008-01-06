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
qoq.c
qsign.c
qsump.c
scdd.c
scdd_f.c
*/

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "rcdd.h"

static R_CMethodDef cMethods[] = {
    {NULL, NULL, 0, NULL, NULL}
};

static R_CallMethodDef callMethods[]  = {
    {"d2q", (DL_FUNC) &d2q, 1},
    {"q2d", (DL_FUNC) &q2d, 1},
    {"q2q", (DL_FUNC) &q2q, 1},
    {"qoq", (DL_FUNC) &qoq, 3},
    {"qsign", (DL_FUNC) &qsign, 1},
    {"qsump", (DL_FUNC) &qsump, 2},
    {"qmatmult", (DL_FUNC) &qmatmult, 2},
    {"lpcdd", (DL_FUNC) &lpcdd, 4},
    {"lpcdd_f", (DL_FUNC) &lpcdd_f, 4},
    {"scdd", (DL_FUNC) &scdd, 7},
    {"scdd_f", (DL_FUNC) &scdd_f, 7},
    {"redundant", (DL_FUNC) &redundant, 2},
    {"redundant_f", (DL_FUNC) &redundant_f, 2},
    {"nonred", (DL_FUNC) &nonred, 2},
    {NULL, NULL, 0}
};

void R_init_rcdd(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

