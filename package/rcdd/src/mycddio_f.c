
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

SEXP rrf_WriteSetFamily(ddf_SetFamilyPtr f) {

    if (f == NULL)
        error("WriteSetFamily: requested family is empty");

    long n = f->famsize;
    SEXP foo;
    PROTECT(foo = allocVector(VECSXP, n));
    long i;
    for (i = 0; i < n; i++)
        SET_VECTOR_ELT(foo, i, rrf_set_fwrite(f->set[i]));
    UNPROTECT(1);
    return foo;
}

SEXP rrf_set_fwrite(set_type set) {

    unsigned long elem;
    long k;
    long card = 0;

    for (elem = 1; elem <= set[0]; elem++)
        if (set_member(elem, set))
            card++;

    SEXP bar;
    PROTECT(bar = allocVector(INTSXP, card));

    for (elem = 1, k = 0; elem <= set[0]; elem++)
        if (set_member(elem, set)) {
            if (k < card)
                INTEGER(bar)[k++] = elem;
            else
                error("Cannot happen!  failure writing set.");
        }

    UNPROTECT(1);
    return bar;
}
