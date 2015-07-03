
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
#include "cdd.h"
#include <R.h>
#include <Rinternals.h>
#include "mycddio.h"
#include <string.h>
#include "rcdd.h"

SEXP scdd(SEXP m, SEXP h, SEXP roworder, SEXP adjacency,
    SEXP inputadjacency, SEXP incidence, SEXP inputincidence)
{
    int i, j, k;

    GetRNGstate();
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

    if (! isString(m))
        error("'m' must be character");

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

    for (i = 0; i < nrow; i++) {
        char *foo = (char *) CHAR(STRING_ELT(m, i));
        if (strlen(foo) != 1)
            error("column one of 'm' not zero-or-one valued");
        if (! (foo[0] == '0' || foo[0] == '1'))
            error("column one of 'm' not zero-or-one valued");
    }
    if (! LOGICAL(h)[0])
        for (i = nrow; i < 2 * nrow; i++) {
            char *foo = (char *) CHAR(STRING_ELT(m, i));
            if (strlen(foo) != 1)
                error("column two of 'm' not zero-or-one valued");
            if (! (foo[0] == '0' || foo[0] == '1'))
                error("column two of 'm' not zero-or-one valued");
        }

    dd_set_global_constants();

    /* note actual type of "value" is mpq_t (defined in cddmp.h) */
    mytype value;
    dd_init(value);

    dd_MatrixPtr mf = dd_CreateMatrix(nrow, ncol - 1);
    /* note our matrix has one more column than Fukuda's */

    /* representation */
    if(LOGICAL(h)[0])
        mf->representation = dd_Inequality;
    else
        mf->representation = dd_Generator;

    mf->numbtype = dd_Rational;

    /* linearity */
    for (i = 0; i < nrow; i++) {
        char *foo = (char *) CHAR(STRING_ELT(m, i));
        if (foo[0] == '1')
            set_addelem(mf->linset, i + 1);
        /* note conversion from zero-origin to one-origin indexing */
    }

    /* matrix */
    for (j = 1, k = nrow; j < ncol; j++)
        for (i = 0; i < nrow; i++, k++) {
            char *rat_str = (char *) CHAR(STRING_ELT(m, k));
            if (mpq_set_str(value, rat_str, 10) == -1)
                error("error converting string to GMP rational");
            mpq_canonicalize(value);
            dd_set(mf->matrix[i][j - 1], value);
            /* note our matrix has one more column than Fukuda's */
        }

    dd_RowOrderType strategy = dd_LexMin;
    char *row_str = (char *) CHAR(STRING_ELT(roworder, 0));
    if(strcmp(row_str, "maxindex") == 0)
        strategy = dd_MaxIndex;
    else if(strcmp(row_str, "minindex") == 0)
        strategy = dd_MinIndex;
    else if(strcmp(row_str, "mincutoff") == 0)
        strategy = dd_MinCutoff;
    else if(strcmp(row_str, "maxcutoff") == 0)
        strategy = dd_MaxCutoff;
    else if(strcmp(row_str, "mixcutoff") == 0)
        strategy = dd_MixCutoff;
    else if(strcmp(row_str, "lexmin") == 0)
        strategy = dd_LexMin;
    else if(strcmp(row_str, "lexmax") == 0)
        strategy = dd_LexMax;
    else if(strcmp(row_str, "randomrow") == 0)
        strategy = dd_RandomRow;
    else
        error("roworder not recognized");

    dd_ErrorType err = dd_NoError;
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly2(mf, strategy, &err);

    if (poly->child != NULL && poly->child->CompStatus == dd_InProgress) {
        dd_FreeMatrix(mf);
        dd_FreePolyhedra(poly);
        dd_clear(value);
        dd_free_global_constants();
        error("Computation failed, rational arithmetic problem\n");
    }

    if (err != dd_NoError) {
        rr_WriteErrorMessages(err);
        dd_FreeMatrix(mf);
        dd_FreePolyhedra(poly);
        dd_clear(value);
        dd_free_global_constants();
        error("failed");
    }

    dd_MatrixPtr aout = NULL;
    if (poly->representation == dd_Inequality)
        aout = dd_CopyGenerators(poly);
    else if (poly->representation == dd_Generator)
        aout = dd_CopyInequalities(poly);
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
    PROTECT(bar = allocMatrix(STRSXP, mrow, ncol));

    /* linearity output */
    for (i = 0; i < mrow; i++)
        if (set_member(i + 1, aout->linset))
            SET_STRING_ELT(bar, i, mkChar("1"));
        else
            SET_STRING_ELT(bar, i, mkChar("0"));
    /* note conversion from zero-origin to one-origin indexing */

    /* matrix output */
    for (j = 1, k = mrow; j < ncol; j++)
        for (i = 0; i < mrow; i++, k++) {
            dd_set(value, aout->matrix[i][j - 1]);
            /* note our matrix has one more column than Fukuda's */
            char *zstr = NULL;
            zstr = mpq_get_str(zstr, 10, value);
            SET_STRING_ELT(bar, k, mkChar(zstr));
            free(zstr);
        }

    int nresult = 1;

    SEXP baz_adj = NULL;
    if (LOGICAL(adjacency)[0]) {
        dd_SetFamilyPtr sout = dd_CopyAdjacency(poly);
        PROTECT(baz_adj = rr_WriteSetFamily(sout));
        dd_FreeSetFamily(sout);
        nresult++;
    }

    SEXP baz_inp_adj = NULL;
    if (LOGICAL(inputadjacency)[0]) {
        dd_SetFamilyPtr sout = dd_CopyInputAdjacency(poly);
        PROTECT(baz_inp_adj = rr_WriteSetFamily(sout));
        dd_FreeSetFamily(sout);
        nresult++;
    }

    SEXP baz_inc = NULL;
    if (LOGICAL(incidence)[0]) {
        dd_SetFamilyPtr sout = dd_CopyIncidence(poly);
        PROTECT(baz_inc = rr_WriteSetFamily(sout));
        dd_FreeSetFamily(sout);
        nresult++;
    }

    SEXP baz_inp_inc = NULL;
    if (LOGICAL(inputincidence)[0]) {
        dd_SetFamilyPtr sout = dd_CopyInputIncidence(poly);
        PROTECT(baz_inp_inc = rr_WriteSetFamily(sout));
        dd_FreeSetFamily(sout);
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

    if (aout->objective != dd_LPnone)
        error("Cannot happen!  aout->objective != dd_LPnone\n");

    dd_FreeMatrix(aout);
    dd_FreeMatrix(mf);
    dd_FreePolyhedra(poly);
    dd_clear(value);
    dd_free_global_constants();

    UNPROTECT(2 + nresult);
    PutRNGstate();
    return result;
}

