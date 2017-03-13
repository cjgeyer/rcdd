
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

SEXP qsump(SEXP foo, SEXP op)
{
    if ((! isString(foo)))
        error("argument must be character");
    int n = LENGTH(foo);

    if (! isInteger(op))
        error("'op' must be integer");
    if (LENGTH(op) != 1)
        error("'op' must be scalar");
    int the_op = INTEGER(op)[0];
    if (the_op <= 0 || the_op > 2)
        error("'op' not recognized, must be 1 (+), 2 (*)");

    mpq_t value1, value3;
    mpq_init(value1);
    mpq_init(value3);
    if (the_op == 2)
        mpq_set_si(value3, 1, 1);

    int k;
    for (k = 0; k < n; k++) {

        const char *zstr1 = CHAR(STRING_ELT(foo, k));
        if (mpq_set_str(value1, zstr1, 10) == -1) {
            mpq_clear(value1);
            mpq_clear(value3);
            error("error converting string to GMP rational");
        }
        mpq_canonicalize(value1);

        switch (the_op) {
            case 1:
                mpq_add(value3, value1, value3);
                break;
            case 2:
                mpq_mul(value3, value1, value3);
                break;
        }
    }

    SEXP baz;
    char *zstr = mpq_get_str(NULL, 10, value3);
    PROTECT(baz = ScalarString(mkChar(zstr)));
    free(zstr);
    mpq_clear(value1);
    mpq_clear(value3);
    UNPROTECT(1);
    return(baz);
}
