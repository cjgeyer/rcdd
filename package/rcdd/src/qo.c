
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

SEXP qo(SEXP foo, SEXP op)
{
    if ((! isString(foo)))
        error("argument must be character");
    int n = LENGTH(foo);

    if (! isInteger(op))
        error("'op' must be integer");
    if (LENGTH(op) != 1)
        error("'op' must be scalar");
    int the_op = INTEGER(op)[0];
    if (the_op <= 0 || the_op > 3)
        error("'op' not recognized, must be 1 (negation), 2 (absolute value), 3 (inversion)");

    SEXP baz;
    PROTECT(baz = duplicate(foo));

    mpq_t value;
    mpq_init(value);

    for (int k = 0; k < n; k++) {

        const char *zstr = CHAR(STRING_ELT(foo, k));
        if (mpq_set_str(value, zstr, 10) == -1) {
            mpq_clear(value);
            error("error converting string to GMP rational");
        }
        mpq_canonicalize(value);

        switch (the_op) {
            case 1:
                mpq_neg(value, value);
                break;
            case 2:
                mpq_abs(value, value);
                break;
            case 3:
                if (mpq_sgn(value) == 0) {
                    mpq_clear(value);
                    error("rational divide by zero");
                }
                mpq_inv(value, value);
                break;
        }

        char *zstr2 = mpq_get_str(NULL, 10, value);
        SET_STRING_ELT(baz, k, mkChar(zstr2));
        free(zstr2);
    }

    mpq_clear(value);
    UNPROTECT(1);
    return(baz);
}
