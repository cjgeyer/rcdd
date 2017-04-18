
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
#include <R.h>
#include <Rinternals.h>
#include "mycddio_f.h"
#include <string.h>
#include "rcdd.h"

SEXP impliedLinearity_f(SEXP m, SEXP h)
{
    GetRNGstate();
    if (! isMatrix(m))
        error("'m' must be matrix");
    if (! isLogical(h))
        error("'h' must be logical");

    if (LENGTH(h) != 1)
        error("'h' must be scalar");

    if (! isReal(m))
        error("'m' must be double");

    SEXP m_dim;
    PROTECT(m_dim = getAttrib(m, R_DimSymbol));
    int nrow = INTEGER(m_dim)[0];
    int ncol = INTEGER(m_dim)[1];
    UNPROTECT(1);

    if (nrow <= 1)
        error("no use if only one row");
    if (ncol <= 3)
        error("no use if only one col");

    for (int i = 0; i < nrow * ncol; i++)
        if (! R_finite(REAL(m)[i]))
            error("'m' not finite-valued");

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

    ddf_ErrorType err = ddf_NoError;
    ddf_rowset out = ddf_ImplicitLinearityRows(mf, &err);

    if (err != ddf_NoError) {
        rrf_WriteErrorMessages(err);
        ddf_FreeMatrix(mf);
        set_free(out);
        ddf_clear(value);
        ddf_free_global_constants();
        error("failed");
    }

    SEXP foo;
    PROTECT(foo = rrf_set_fwrite(out));

    ddf_FreeMatrix(mf);
    set_free(out);
    ddf_clear(value);
    ddf_free_global_constants();

    PutRNGstate();

    UNPROTECT(1);
    return foo;
}

