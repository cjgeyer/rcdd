
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

SEXP q2q(SEXP foo)
{
    if (! isString(foo))
        error("argument must be character");
    int n = LENGTH(foo);

    SEXP bar;
    PROTECT(bar = duplicate(foo));

    mpq_t value;
    mpq_init(value);

    int k;
    for (k = 0; k < n; k++) {
        char *zstr = (char *) CHAR(STRING_ELT(foo, k));
        if (mpq_set_str(value, zstr, 10) == -1)
            error("error converting string to GMP rational");
        mpq_canonicalize(value);
        zstr = NULL;
        zstr = mpq_get_str(zstr, 10, value);
        SET_STRING_ELT(bar, k, mkChar(zstr));
        free(zstr);
    }

    mpq_clear(value);
    UNPROTECT(1);
    return(bar);
}
