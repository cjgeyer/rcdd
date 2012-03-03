
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

#ifndef RCDD_MYCDDIO_F_H
#define RCDD_MYCDDIO_F_H

#ifndef  _CDDTYPES_HF
#include "cddtypes_f.h"
#endif  /* _CDDTYPES_HF */
#ifndef  __SETOPER_H
#include "setoper.h"
#endif  /* __SETOPER_H */

#include <Rinternals.h>

SEXP rrf_WriteSetFamily(ddf_SetFamilyPtr f);
SEXP rrf_set_fwrite(set_type set);
void rrf_WriteErrorMessages(ddf_ErrorType Error);

#endif /* RCDD_MYCDDIO_F_H */

