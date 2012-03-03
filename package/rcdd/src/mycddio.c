
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

SEXP rr_WriteSetFamily(dd_SetFamilyPtr f) {

    if (f == NULL)
        error("WriteSetFamily: requested family is empty");

    long n = f->famsize;
    SEXP foo;
    PROTECT(foo = allocVector(VECSXP, n));
    long i;
    for (i = 0; i < n; i++)
        SET_VECTOR_ELT(foo, i, rr_set_fwrite(f->set[i]));
    UNPROTECT(1);
    return foo;
}

SEXP rr_set_fwrite(set_type set) {

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
                error("Cannot happen!  failure writing set");
        }

    UNPROTECT(1);
    return bar;
}

void rr_WriteErrorMessages(dd_ErrorType Error)
{
  switch (Error) {

  case dd_DimensionTooLarge:
    REprintf("*Input Error: Input matrix is too large:\n");
    REprintf("*Please increase MMAX and/or NMAX in the source code and recompile.\n");
    break;

  case dd_IFileNotFound:
    REprintf("*Input Error: Specified input file does not exist.\n");
    break;

  case dd_OFileNotOpen:
    REprintf("*Output Error: Specified output file cannot be opened.\n");
    break;

  case dd_NegativeMatrixSize:
    REprintf("*Input Error: Input matrix has a negative size:\n");
    REprintf("*Please check rowsize or colsize.\n");
    break;

  case dd_ImproperInputFormat:
    REprintf("*Input Error: Input format is not correct.\n");
    REprintf("*Format:\n");
    REprintf(" begin\n");
    REprintf("   m   n  NumberType(real, rational or integer)\n");
    REprintf("   b  -A\n");
    REprintf(" end\n");
    break;

  case dd_EmptyVrepresentation:
    REprintf("*Input Error: V-representation is empty:\n");
    REprintf("*cddlib does not accept this trivial case for which output can be any inconsistent system.\n");
    break;

  case dd_EmptyHrepresentation:
    REprintf("*Input Error: H-representation is empty.\n");
    break;

  case dd_EmptyRepresentation:
    REprintf("*Input Error: Representation is empty.\n");
    break;

  case dd_NoLPObjective:
    REprintf("*LP Error: No LP objective (max or min) is set.\n");
    break;

  case dd_NoRealNumberSupport:
    REprintf("*LP Error: The binary (with GMP Rational) does not support Real number input.\n");
    REprintf("         : Use a binary compiled without -DGMPRATIONAL option.\n");
    break;

 case dd_NotAvailForH:
    REprintf("*Error: A function is called with H-rep which does not support an H-representation.\n");
    break;

 case dd_NotAvailForV:
    REprintf("*Error: A function is called with V-rep which does not support an V-representation.\n");
    break;

 case dd_CannotHandleLinearity:
    REprintf("*Error: The function called cannot handle linearity.\n");
    break;

 case dd_RowIndexOutOfRange:
    REprintf("*Error: Specified row index is out of range\n");
    break;

 case dd_ColIndexOutOfRange:
    REprintf("*Error: Specified column index is out of range\n");
    break;

 case dd_LPCycling:
    REprintf("*Error: Possibly an LP cycling occurs.  Use the Criss-Cross method.\n");
    break;
    
 case dd_NumericallyInconsistent:
    REprintf("*Error: Numerical inconsistency is found.  Use the GMP exact arithmetic.\n");
    break;
    
  case dd_NoError:
    REprintf("*No Error found.\n");
    break;
  }
}

