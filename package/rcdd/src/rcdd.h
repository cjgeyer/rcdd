
/*
 *  cdd an R interface to cddlib
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

#ifndef RCDD_RCDD_H
#define RCDD_RCDD_H

#include <Rinternals.h>

SEXP d2q(SEXP foo);
SEXP q2d(SEXP foo);
SEXP q2q(SEXP foo);
SEXP qo(SEXP foo, SEXP op);
SEXP qoq(SEXP foo, SEXP bar, SEXP op);
SEXP qsign(SEXP foo);
SEXP qsump(SEXP foo, SEXP op);
SEXP qminp(SEXP foo, SEXP op);
SEXP qmatmult(SEXP foo, SEXP bar);

SEXP lpcdd(SEXP hrep, SEXP objfun, SEXP minimize, SEXP solver);
SEXP lpcdd_f(SEXP hrep, SEXP objfun, SEXP minimize, SEXP solver);

SEXP scdd(SEXP m, SEXP h, SEXP roworder, SEXP adjacency,
    SEXP inputadjacency, SEXP incidence, SEXP inputincidence);
SEXP scdd_f(SEXP m, SEXP h, SEXP roworder, SEXP adjacency,
    SEXP inputadjacency, SEXP incidence, SEXP inputincidence);

SEXP redundant(SEXP m, SEXP h);
SEXP redundant_f(SEXP m, SEXP h);

SEXP allfaces(SEXP m);
SEXP allfaces_f(SEXP m);

SEXP nonred(SEXP sets, SEXP pow2);
SEXP test_my_subset(SEXP set1, SEXP set2, SEXP pow2);
SEXP all_intersect(SEXP sets, SEXP pow2);
SEXP all_union(SEXP sets, SEXP pow2);

SEXP impliedLinearity(SEXP m, SEXP h);
SEXP impliedLinearity_f(SEXP m, SEXP h);

#endif /* RCDD_RCDD_H */

