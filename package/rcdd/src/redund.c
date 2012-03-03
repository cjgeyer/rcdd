
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
#include <Rinternals.h>
#include "mycddio.h"
#include <string.h>
#include "rcdd.h"
#ifdef WOOF
#include <stdio.h>
#endif /* WOOF */

SEXP redundant(SEXP m, SEXP h)
{
    if (! isString(m))
        error("'m' must be character");
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
        char *foo = (char *) CHAR(STRING_ELT(m, i));
        if (strlen(foo) != 1)
            error("column one of 'm' not zero-or-one valued");
        if (! (foo[0] == '0' || foo[0] == '1'))
            error("column one of 'm' not zero-or-one valued");
    }
    if (! LOGICAL(h)[0])
        for (int i = nrow; i < 2 * nrow; i++) {
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
    for (int i = 0; i < nrow; i++) {
        char *foo = (char *) CHAR(STRING_ELT(m, i));
        if (foo[0] == '1')
            set_addelem(mf->linset, i + 1);
        /* note conversion from zero-origin to one-origin indexing */
    }

    /* matrix */
    for (int j = 1, k = nrow; j < ncol; j++)
        for (int i = 0; i < nrow; i++, k++) {
            char *rat_str = (char *) CHAR(STRING_ELT(m, k));
            if (mpq_set_str(value, rat_str, 10) == -1)
                error("error converting string to GMP rational");
            mpq_canonicalize(value);
            dd_set(mf->matrix[i][j - 1], value);
            /* note our matrix has one more column than Fukuda's */
        }

    dd_rowset impl_linset, redset;
    dd_rowindex newpos;
    dd_ErrorType err = dd_NoError;

    dd_MatrixCanonicalize(&mf, &impl_linset, &redset, &newpos, &err);

    if (err != dd_NoError) {
        rr_WriteErrorMessages(err);
        dd_FreeMatrix(mf);
        dd_clear(value);
        dd_free_global_constants();
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
    PROTECT(bar = allocMatrix(STRSXP, mrow, ncol));

    /* linearity output */
    for (int i = 0; i < mrow; i++)
        if (set_member(i + 1, mf->linset))
            SET_STRING_ELT(bar, i, mkChar("1"));
        else
            SET_STRING_ELT(bar, i, mkChar("0"));
    /* note conversion from zero-origin to one-origin indexing */

    /* matrix output */
    for (int j = 1, k = mrow; j < ncol; j++)
        for (int i = 0; i < mrow; i++, k++) {
            dd_set(value, mf->matrix[i][j - 1]);
            /* note our matrix has one more column than Fukuda's */
            char *zstr = NULL;
            zstr = mpq_get_str(zstr, 10, value);
            SET_STRING_ELT(bar, k, mkChar(zstr));
            free(zstr);
        }

    if (mf->representation == dd_Inequality) {
        SEXP attr_name, attr_value;
        PROTECT(attr_name = ScalarString(mkChar("representation")));
        PROTECT(attr_value = ScalarString(mkChar("H")));
        setAttrib(bar, attr_name, attr_value);
        UNPROTECT(2);
    }
    if (mf->representation == dd_Generator) {
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
        PROTECT(baz = rr_set_fwrite(impl_linset));
        nresult++;
    }

    SEXP qux = NULL;
    if (red_size > 0) {
        PROTECT(qux = rr_set_fwrite(redset));
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
    dd_FreeMatrix(mf);
    dd_clear(value);
    dd_free_global_constants();

    UNPROTECT(nresult + 2);
    return result;
}

