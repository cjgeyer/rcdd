
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

#include <Rinternals.h>
#include <gmp.h>
#include <stdlib.h>
#include "rcdd.h"

SEXP qmatmult(SEXP foo, SEXP bar)
{
    if ((! isString(foo)) || (! isString(bar)))
        error("arguments must be character");
    if ((! isMatrix(foo)) || (! isMatrix(bar)))
        error("arguments must be matrices");

    SEXP dim;
    PROTECT(dim = getAttrib(foo, R_DimSymbol));
    int nrow_foo = INTEGER(dim)[0];
    int ncol_foo = INTEGER(dim)[1];
    UNPROTECT(1);

    PROTECT(dim = getAttrib(bar, R_DimSymbol));
    int nrow_bar = INTEGER(dim)[0];
    int ncol_bar = INTEGER(dim)[1];
    UNPROTECT(1);

    if (nrow_foo <= 0)
        error("row dimension of 1st arg must be positive");
    if (ncol_foo <= 0)
        error("col dimension of 1st arg must be positive");
    if (nrow_bar <= 0)
        error("row dimension of 2nd arg must be positive");
    if (ncol_bar <= 0)
        error("col dimension of 2nd arg must be positive");
    if (ncol_foo != nrow_bar)
        error("col dimension of 1st arg must match row dimension of 2nd arg");

    SEXP baz;
    PROTECT(baz = allocMatrix(STRSXP, nrow_foo, ncol_bar));

    mpq_t value1, value2, value3;
    mpq_init(value1);
    mpq_init(value2);
    mpq_init(value3);

    for (int i = 0; i < nrow_foo; ++i) {
        for (int j = 0; j < ncol_bar; ++j) {

            mpq_set_si(value3, 0, 1);

            for (int k = 0; k < ncol_foo; ++k) {

                const char *foo_ik = CHAR(STRING_ELT(foo, i + nrow_foo * k));
                const char *bar_kj = CHAR(STRING_ELT(bar, k + nrow_bar * j));

                if (mpq_set_str(value1, foo_ik, 10) == -1) {
                    mpq_clear(value1);
                    mpq_clear(value2);
                    mpq_clear(value3);
                    error("error converting string to GMP rational");
                }
                mpq_canonicalize(value1);

                if (mpq_set_str(value2, bar_kj, 10) == -1) {
                    mpq_clear(value1);
                    mpq_clear(value2);
                    mpq_clear(value3);
                    error("error converting string to GMP rational");
                }
                mpq_canonicalize(value2);

                mpq_mul(value2, value1, value2);
                mpq_add(value3, value3, value2);

            }

            char *baz_ij = mpq_get_str(NULL, 10, value3);
            SET_STRING_ELT(baz, i + nrow_foo * j, mkChar(baz_ij));
            free(baz_ij);

        }
    }

    mpq_clear(value1);
    mpq_clear(value2);
    mpq_clear(value3);
    UNPROTECT(1);
    return(baz);
}
