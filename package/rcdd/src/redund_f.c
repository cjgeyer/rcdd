
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
#include "mycddio_f.h"
#include <string.h>
#include "rcdd.h"
#ifdef WOOF
#include <stdio.h>
#endif /* WOOF */

SEXP redundant_f(SEXP m, SEXP h)
{
    GetRNGstate();
    if (! isReal(m))
        error("'m' must be double");
    if (! isMatrix(m))
        error("'m' must be matrix");
    if (! isLogical(h))
        error("'h' must be logical");
    if (LENGTH(h) != 1)
        error("'h' must be scalar");

    SEXP m_dim;
    PROTECT(m_dim = getAttrib(m, R_DimSymbol));
    int nrow = INTEGER(m_dim)[0];
    int ncol = INTEGER(m_dim)[1];
    UNPROTECT(1);

#ifdef WOOF
    printf("nrow = %d\n", nrow);
    printf("ncol = %d\n", ncol);
#endif /* WOOF */

    if (nrow < 2)
        error("less than 2 rows, cannot be redundant");
    if (ncol <= 2)
        error("no cols in m[ , - c(1, 2)]");

   for (int i = 0; i < nrow; i++) {
        double foo = REAL(m)[i];
        if (! (foo == 0.0 || foo == 1.0))
            error("column one of 'm' not zero-or-one valued");
    }
    if (! LOGICAL(h)[0])
        for (int i = nrow; i < 2 * nrow; i++) {
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
    for (int i = 0; i < nrow; i++) {
        double foo = REAL(m)[i];
        if (foo == 1.0)
            set_addelem(mf->linset, i + 1);
        /* note conversion from zero-origin to one-origin indexing */
    }

    /* matrix */
    for (int j = 1, k = nrow; j < ncol; j++)
        for (int i = 0; i < nrow; i++, k++) {
            ddf_set_d(value, REAL(m)[k]);
            ddf_set(mf->matrix[i][j - 1], value);
            /* note our matrix has one more column than Fukuda's */
        }

    ddf_rowset impl_linset, redset;
    ddf_rowindex newpos;
    ddf_ErrorType err = ddf_NoError;

    ddf_MatrixCanonicalize(&mf, &impl_linset, &redset, &newpos, &err);

    if (err != ddf_NoError) {
        rrf_WriteErrorMessages(err);
        ddf_FreeMatrix(mf);
        ddf_clear(value);
        ddf_free_global_constants();
        error("failed");
    }

    int mrow = mf->rowsize;
    int mcol = mf->colsize;

    if (mcol + 1 != ncol)
        error("Cannot happen!  computed matrix has wrong number of columns");

#ifdef WOOF
    printf("mrow = %d\n", mrow);
    printf("mcol = %d\n", mcol);
#endif /* WOOF */

    SEXP bar;
    PROTECT(bar = allocMatrix(REALSXP, mrow, ncol));

    /* linearity output */
    for (int i = 0; i < mrow; i++)
        if (set_member(i + 1, mf->linset))
            REAL(bar)[i] = 1.0;
        else
            REAL(bar)[i] = 0.0;
    /* note conversion from zero-origin to one-origin indexing */

    /* matrix output */
    for (int j = 1, k = mrow; j < ncol; j++)
        for (int i = 0; i < mrow; i++, k++) {
            double ax = ddf_get_d(mf->matrix[i][j - 1]);
            /* note our matrix has one more column than Fukuda's */
            REAL(bar)[k] = ax;
        }

    if (mf->representation == ddf_Inequality) {
        SEXP attr_name, attr_value;
        PROTECT(attr_name = ScalarString(mkChar("representation")));
        PROTECT(attr_value = ScalarString(mkChar("H")));
        setAttrib(bar, attr_name, attr_value);
        UNPROTECT(2);
    }
    if (mf->representation == ddf_Generator) {
        SEXP attr_name, attr_value;
        PROTECT(attr_name = ScalarString(mkChar("representation")));
        PROTECT(attr_value = ScalarString(mkChar("V")));
        setAttrib(bar, attr_name, attr_value);
        UNPROTECT(2);
    }

    int impl_size = set_card(impl_linset);
    int red_size = set_card(redset);

    int nresult = 1;
    int iresult = 1;

    SEXP baz = NULL;
    if (impl_size > 0) {
        PROTECT(baz = rrf_set_fwrite(impl_linset));
        nresult++;
    }

    SEXP qux = NULL;
    if (red_size > 0) {
        PROTECT(qux = rrf_set_fwrite(redset));
        nresult++;
    }

    SEXP fred = NULL;
    {
        PROTECT(fred = allocVector(INTSXP, nrow));
        for (int i = 1; i <= nrow; i++)
            INTEGER(fred)[i - 1] = newpos[i];
        nresult++;
    }

#ifdef WOOF
    fprintf(stderr, "impl_size = %d\n", impl_size);
    fprintf(stderr, "red_size = %d\n", red_size);
    fprintf(stderr, "nresult = %d\n", nresult);
    if (baz)
        fprintf(stderr, "LENGTH(baz) = %d\n", LENGTH(baz));
    if (qux)
        fprintf(stderr, "LENGTH(qux) = %d\n", LENGTH(qux));
#endif /* WOOF */

    SEXP result, resultnames;
    PROTECT(result = allocVector(VECSXP, nresult));
    PROTECT(resultnames = allocVector(STRSXP, nresult));

    SET_STRING_ELT(resultnames, 0, mkChar("output"));
    SET_VECTOR_ELT(result, 0, bar);
    if (baz) {
        SET_STRING_ELT(resultnames, iresult, mkChar("implied.linearity"));
        SET_VECTOR_ELT(result, iresult, baz);
        iresult++;
    }
    if (qux) {
        SET_STRING_ELT(resultnames, iresult, mkChar("redundant"));
        SET_VECTOR_ELT(result, iresult, qux);
        iresult++;
    }
    {
        SET_STRING_ELT(resultnames, iresult, mkChar("new.position"));
        SET_VECTOR_ELT(result, iresult, fred);
        iresult++;
    }
    namesgets(result, resultnames);

    set_free(redset);
    set_free(impl_linset);
    free(newpos);
    ddf_FreeMatrix(mf);
    ddf_clear(value);
    ddf_free_global_constants();

    UNPROTECT(nresult + 2);
    PutRNGstate();
    return result;
}

