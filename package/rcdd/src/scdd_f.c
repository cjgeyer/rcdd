
/*
 *  rcdd an R interface to cddlib
 *  Copyright (C) 2005    Charles J. Geyer
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 59 Temple Place,
 *  Suite 330, Boston, MA  02111-1307  USA.
 */

#include "setoper.h"
#include "cdd_f.h"
#include <Rinternals.h>
#include "mycddio_f.h"
#include <string.h>
#include "rcdd.h"

SEXP scdd_f(SEXP m, SEXP h, SEXP roworder, SEXP adjacency,
    SEXP inputadjacency, SEXP incidence, SEXP inputincidence)
{
    int i, j, k;

    if (! isMatrix(m))
        error("'m' must be matrix");
    if (! isLogical(h))
        error("'h' must be logical");
    if (! isString(roworder))
        error("'roworder' must be character");
    if (! isLogical(adjacency))
        error("'adjacency' must be logical");
    if (! isLogical(inputadjacency))
        error("'inputadjacency' must be logical");
    if (! isLogical(incidence))
        error("'incidence' must be logical");
    if (! isLogical(inputincidence))
        error("'inputincidence' must be logical");

    if (LENGTH(h) != 1)
        error("'h' must be scalar");
    if (LENGTH(roworder) != 1)
        error("'roworder' must be scalar");
    if (LENGTH(adjacency) != 1)
        error("'adjacency' must be scalar");
    if (LENGTH(inputadjacency) != 1)
        error("'inputadjacency' must be scalar");
    if (LENGTH(incidence) != 1)
        error("'incidence' must be scalar");
    if (LENGTH(inputincidence) != 1)
        error("'inputincidence' must be scalar");

    if (! isReal(m))
        error("'m' must be double");

    SEXP m_dim;
    PROTECT(m_dim = getAttrib(m, R_DimSymbol));
    int nrow = INTEGER(m_dim)[0];
    int ncol = INTEGER(m_dim)[1];
    UNPROTECT(1);

#ifdef BLATHER
    printf("nrow = %d\n", nrow);
    printf("ncol = %d\n", ncol);
#endif /* BLATHER */

    if ((! LOGICAL(h)[0]) && nrow <= 0)
        error("no rows in 'm', not allowed for V-representation");
    if (ncol <= 2)
        error("no cols in m[ , - c(1, 2)]");

    for (i = 0; i < nrow * ncol; i++)
        if (! R_finite(REAL(m)[i]))
            error("'m' not finite-valued");

    for (i = 0; i < nrow; i++) {
        double foo = REAL(m)[i];
        if (! (foo == 0.0 || foo == 1.0))
            error("column one of 'm' not zero-or-one valued");
    }
    if (! LOGICAL(h)[0])
        for (i = nrow; i < 2 * nrow; i++) {
            double foo = REAL(m)[i];
            if (! (foo == 0.0 || foo == 1.0))
                error("column two of 'm' not zero-or-one valued");
        }

    ddf_set_global_constants();

    myfloat value;
    ddf_init(value);

    ddf_MatrixPtr mf = ddf_CreateMatrix(nrow, ncol - 1);
    /* note our matrix has one more column than Fukuda's */

    /* representation */
    if(LOGICAL(h)[0])
        mf->representation = ddf_Inequality;
    else
        mf->representation = ddf_Generator;

    mf->numbtype = ddf_Real;

    /* linearity */
    for (i = 0; i < nrow; i++) {
        double foo = REAL(m)[i];
        if (foo == 1.0)
            set_addelem(mf->linset, i + 1);
        /* note conversion from zero-origin to one-origin indexing */
    }

    /* matrix */
    for (j = 1, k = nrow; j < ncol; j++)
        for (i = 0; i < nrow; i++, k++) {
            ddf_set_d(value, REAL(m)[k]);
            ddf_set(mf->matrix[i][j - 1], value);
            /* note our matrix has one more column than Fukuda's */
        }

    ddf_RowOrderType strategy = ddf_LexMin;
    const char *row_str = CHAR(STRING_ELT(roworder, 0));
    if(strcmp(row_str, "maxindex") == 0)
        strategy = ddf_MaxIndex;
    else if(strcmp(row_str, "minindex") == 0)
        strategy = ddf_MinIndex;
    else if(strcmp(row_str, "mincutoff") == 0)
        strategy = ddf_MinCutoff;
    else if(strcmp(row_str, "maxcutoff") == 0)
        strategy = ddf_MaxCutoff;
    else if(strcmp(row_str, "mixcutoff") == 0)
        strategy = ddf_MixCutoff;
    else if(strcmp(row_str, "lexmin") == 0)
        strategy = ddf_LexMin;
    else if(strcmp(row_str, "lexmax") == 0)
        strategy = ddf_LexMax;
    else if(strcmp(row_str, "randomrow") == 0)
        strategy = ddf_RandomRow;
    else
        error("roworder not recognized");

    ddf_ErrorType err = ddf_NoError;
    ddf_PolyhedraPtr poly = ddf_DDMatrix2Poly2(mf, strategy, &err);

    if (poly->child != NULL && poly->child->CompStatus == ddf_InProgress) {
        ddf_FreeMatrix(mf);
        ddf_FreePolyhedra(poly);
        ddf_clear(value);
        ddf_free_global_constants();
        error("Computation failed, floating-point arithmetic problem\n");
    }

    if (err != ddf_NoError) {
        rrf_WriteErrorMessages(err);
        ddf_FreeMatrix(mf);
        ddf_FreePolyhedra(poly);
        ddf_clear(value);
        ddf_free_global_constants();
        error("failed");
    }

    ddf_MatrixPtr aout = NULL;
    if (poly->representation == ddf_Inequality)
        aout = ddf_CopyGenerators(poly);
    else if (poly->representation == ddf_Generator)
        aout = ddf_CopyInequalities(poly);
    else
        error("Cannot happen!  poly->representation no good\n");
    if (aout == NULL)
        error("Cannot happen!  aout no good\n");

    int mrow = aout->rowsize;
    int mcol = aout->colsize;

    if (mcol + 1 != ncol)
        error("Cannot happen!  computed matrix has wrong number of columns");

#ifdef BLATHER
    printf("mrow = %d\n", mrow);
    printf("mcol = %d\n", mcol);
#endif /* BLATHER */

    SEXP bar;
    PROTECT(bar = allocMatrix(REALSXP, mrow, ncol));

    /* linearity output */
    for (i = 0; i < mrow; i++)
        if (set_member(i + 1, aout->linset))
            REAL(bar)[i] = 1.0;
        else
            REAL(bar)[i] = 0.0;
    /* note conversion from zero-origin to one-origin indexing */

    /* matrix output */
    for (j = 1, k = mrow; j < ncol; j++)
        for (i = 0; i < mrow; i++, k++) {
            double ax = ddf_get_d(aout->matrix[i][j - 1]);
            /* note our matrix has one more column than Fukuda's */
            REAL(bar)[k] = ax;
        }

    int nresult = 1;

    SEXP baz_adj = NULL;
    if (LOGICAL(adjacency)[0]) {
        ddf_SetFamilyPtr sout = ddf_CopyAdjacency(poly);
        PROTECT(baz_adj = rrf_WriteSetFamily(sout));
        ddf_FreeSetFamily(sout);
        nresult++;
    }

    SEXP baz_inp_adj = NULL;
    if (LOGICAL(inputadjacency)[0]) {
        ddf_SetFamilyPtr sout = ddf_CopyInputAdjacency(poly);
        PROTECT(baz_inp_adj = rrf_WriteSetFamily(sout));
        ddf_FreeSetFamily(sout);
        nresult++;
    }

    SEXP baz_inc = NULL;
    if (LOGICAL(incidence)[0]) {
        ddf_SetFamilyPtr sout = ddf_CopyIncidence(poly);
        PROTECT(baz_inc = rrf_WriteSetFamily(sout));
        ddf_FreeSetFamily(sout);
        nresult++;
    }

    SEXP baz_inp_inc = NULL;
    if (LOGICAL(inputincidence)[0]) {
        ddf_SetFamilyPtr sout = ddf_CopyInputIncidence(poly);
        PROTECT(baz_inp_inc = rrf_WriteSetFamily(sout));
        ddf_FreeSetFamily(sout);
        nresult++;
    }

    SEXP result, resultnames;
    PROTECT(result = allocVector(VECSXP, nresult));
    PROTECT(resultnames = allocVector(STRSXP, nresult));

    SET_STRING_ELT(resultnames, 0, mkChar("output"));
    SET_VECTOR_ELT(result, 0, bar);

    int iresult = 1;

    if (baz_adj) {
        SET_STRING_ELT(resultnames, iresult, mkChar("adjacency"));
        SET_VECTOR_ELT(result, iresult, baz_adj);
        iresult++;
    }
    if (baz_inp_adj) {
        SET_STRING_ELT(resultnames, iresult, mkChar("inputadjacency"));
        SET_VECTOR_ELT(result, iresult, baz_inp_adj);
        iresult++;
    }
    if (baz_inc) {
        SET_STRING_ELT(resultnames, iresult, mkChar("incidence"));
        SET_VECTOR_ELT(result, iresult, baz_inc);
        iresult++;
    }
    if (baz_inp_inc) {
        SET_STRING_ELT(resultnames, iresult, mkChar("inputincidence"));
        SET_VECTOR_ELT(result, iresult, baz_inp_inc);
        iresult++;
    }
    namesgets(result, resultnames);

    if (aout->objective != ddf_LPnone)
        error("Cannot happen! aout->objective != ddf_LPnone\n");

    ddf_FreeMatrix(aout);
    ddf_FreeMatrix(mf);
    ddf_FreePolyhedra(poly);
    ddf_clear(value);
    ddf_free_global_constants();

    UNPROTECT(2 + nresult);
    return result;
}

