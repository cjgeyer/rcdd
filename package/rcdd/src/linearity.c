
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

SEXP impliedLinearity(SEXP m, SEXP h)
{
    if (! isMatrix(m))
        error("'m' must be matrix");
    if (! isLogical(h))
        error("'h' must be logical");

    if (LENGTH(h) != 1)
        error("'h' must be scalar");

    if (! isString(m))
        error("'m' must be character");

    SEXP m_dim;
    PROTECT(m_dim = getAttrib(m, R_DimSymbol));
    int nrow = INTEGER(m_dim)[0];
    int ncol = INTEGER(m_dim)[1];
    UNPROTECT(1);

    if (nrow <= 1)
        error("no use if only one row");
    if (ncol <= 3)
        error("no use if only one col");

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

    dd_ErrorType err = dd_NoError;
    dd_rowset out = dd_ImplicitLinearityRows(mf, &err);

    if (err != dd_NoError) {
        rr_WriteErrorMessages(err);
        dd_FreeMatrix(mf);
        set_free(out);
        dd_clear(value);
        dd_free_global_constants();
        error("failed");
    }

    SEXP foo;
    PROTECT(foo = rr_set_fwrite(out));

    dd_FreeMatrix(mf);
    set_free(out);
    dd_clear(value);
    dd_free_global_constants();

    UNPROTECT(1);
    return foo;
}

