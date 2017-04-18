
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
#include <R.h>
#include <Rinternals.h>
#include "mycddio_f.h"
#include <string.h>
#ifdef MOO
#include <stdio.h>
#endif /* MOO */
#include "rcdd.h"

static SEXP FaceEnum(ddf_MatrixPtr M);

SEXP allfaces_f(SEXP hrep)
{
    GetRNGstate();
    if (! isMatrix(hrep))
        error("'hrep' must be matrix");
    if (! isReal(hrep))
        error("'hrep' must be double");

    SEXP hrep_dim;
    PROTECT(hrep_dim = getAttrib(hrep, R_DimSymbol));
    int nrow = INTEGER(hrep_dim)[0];
    int ncol = INTEGER(hrep_dim)[1];
    UNPROTECT(1);

    if (nrow <= 0)
        error("no rows in 'hrep'");
    if (ncol <= 3)
        error("three or fewer cols in hrep");

    for (int i = 0; i < nrow * ncol; ++i)
        if (! R_finite(REAL(hrep)[i]))
            error("'hrep' not finite-valued");

    for (int i = 0; i < nrow; ++i) {
        double foo = REAL(hrep)[i];
        if (! (foo == 0.0 || foo == 1.0))
            error("column one of 'hrep' not zero-or-one valued");
    }

    ddf_set_global_constants();

    myfloat value;
    ddf_init(value);

    ddf_MatrixPtr mf = ddf_CreateMatrix(nrow, ncol - 1);
    /* note our matrix has one more column than Fukuda's */

    mf->representation = ddf_Inequality;
    mf->numbtype = ddf_Real;

    /* linearity */
    for (int i = 0; i < nrow; ++i) {
        double foo = REAL(hrep)[i];
        if (foo == 1.0)
            set_addelem(mf->linset, i + 1);
        /* note conversion from zero-origin to one-origin indexing */
    }

    /* matrix */
    for (int j = 1, k = nrow; j < ncol; ++j)
        for (int i = 0; i < nrow; ++i, ++k) {
            ddf_set_d(value, REAL(hrep)[k]);
            ddf_set(mf->matrix[i][j - 1], value);
            /* note our matrix has one more column than Fukuda's */
        }

    SEXP result;
    PROTECT(result = FaceEnum(mf));

    ddf_FreeMatrix(mf);
    ddf_clear(value);
    ddf_free_global_constants();

    if (result == R_NilValue)
        error("failed");

    PutRNGstate();

    UNPROTECT(1);
    return result;
}

static ddf_ErrorType FaceEnumHelper(ddf_MatrixPtr M, ddf_rowset R, ddf_rowset S);

static SEXP dimlist, riplist, activelist;

static PROTECT_INDEX dimidx, ripidx, activeidx;

static SEXP FaceEnum(ddf_MatrixPtr M)
{
    PROTECT_WITH_INDEX(dimlist = R_NilValue, &dimidx);
    PROTECT_WITH_INDEX(riplist = R_NilValue, &ripidx);
    PROTECT_WITH_INDEX(activelist = R_NilValue, &activeidx);

    ddf_rowset R, S;
    set_initialize(&R, M->rowsize);
    set_initialize(&S, M->rowsize);

    ddf_ErrorType err = FaceEnumHelper(M, R, S);

    set_free(R);
    set_free(S);

    if (err != ddf_NoError) {
#ifdef MOO
        switch (err) {
            case ddf_DimensionTooLarge:
                fprintf(stderr, "err = ddf_DimensionTooLarge\n");
                break;
            case ddf_ImproperInputFormat:
                fprintf(stderr, "err = ddf_ImproperInputFormat\n");
                break;
            case ddf_NegativeMatrixSize:
                fprintf(stderr, "err = ddf_NegativeMatrixSize\n");
                break;
            case ddf_EmptyVrepresentation:
                fprintf(stderr, "err = ddf_EmptyVrepresentation\n");
                break;
            case ddf_EmptyHrepresentation:
                fprintf(stderr, "err = ddf_EmptyHrepresentation\n");
                break;
            case ddf_EmptyRepresentation:
                fprintf(stderr, "err = ddf_EmptyRepresentation\n");
                break;
            case ddf_IFileNotFound:
                fprintf(stderr, "err = ddf_IFileNotFound\n");
                break;
            case ddf_OFileNotOpen:
                fprintf(stderr, "err = ddf_OFileNotOpen\n");
                break;
            case ddf_NoLPObjective:
                fprintf(stderr, "err = ddf_NoLPObjective\n");
                break;
            case ddf_NoRealNumberSupport:
                fprintf(stderr, "err = ddf_NoRealNumberSupport\n");
                break;
            case ddf_NotAvailForH:
                fprintf(stderr, "err = ddf_NotAvailForH\n");
                break;
            case ddf_NotAvailForV:
                fprintf(stderr, "err = ddf_NotAvailForV\n");
                break;
            case ddf_CannotHandleLinearity:
                fprintf(stderr, "err = ddf_CannotHandleLinearity\n");
                break;
            case ddf_RowIndexOutOfRange:
                fprintf(stderr, "err = ddf_RowIndexOutOfRange\n");
                break;
            case ddf_ColIndexOutOfRange:
                fprintf(stderr, "err = ddf_ColIndexOutOfRange\n");
                break;
            case ddf_LPCycling:
                fprintf(stderr, "err = ddf_LPCycling\n");
                break;
            case ddf_NumericallyInconsistent:
                fprintf(stderr, "err = ddf_NumericallyInconsistent\n");
                break;
            case ddf_NoError:
                fprintf(stderr, "err = ddf_NoError\n");
                break;
            default:
                fprintf(stderr, "err bogus, WTF????\n");
        }
#endif /* MOO */
        rrf_WriteErrorMessages(err);
        UNPROTECT(3);
        return R_NilValue;
    }

    SEXP result;
    SEXP resultnames;
    PROTECT(result = allocVector(VECSXP, 3));
    PROTECT(resultnames = allocVector(STRSXP, 3));

    SET_STRING_ELT(resultnames, 0, mkChar("dimension"));
    SET_STRING_ELT(resultnames, 1, mkChar("active.set"));
    SET_STRING_ELT(resultnames, 2, mkChar("relative.interior.point"));
    namesgets(result, resultnames);

    SET_VECTOR_ELT(result, 0, PairToVectorList(dimlist));
    SET_VECTOR_ELT(result, 1, PairToVectorList(activelist));
    SET_VECTOR_ELT(result, 2, PairToVectorList(riplist));

    UNPROTECT(5);
    return result;
}

static ddf_ErrorType FaceEnumHelper(ddf_MatrixPtr M, ddf_rowset R, ddf_rowset S)
{
    ddf_ErrorType err;
    ddf_rowset LL, ImL, RR, SS, Lbasis;
    ddf_rowrange iprev = 0;
    ddf_colrange dim;
    ddf_LPSolutionPtr lps = NULL;

    set_initialize(&LL, M->rowsize);
    set_initialize(&RR, M->rowsize);
    set_initialize(&SS, M->rowsize);
    set_copy(LL, M->linset);
    set_copy(RR, R);
    set_copy(SS, S);

    myfloat value;
    ddf_init(value);

    err = ddf_NoError;
    ddf_boolean foo = ddf_ExistsRestrictedFace(M, R, S, &err);
    if (err != ddf_NoError) {
#ifdef MOO
        fprintf(stderr, "err from ddf_ExistsRestrictedFace\n");
        fprintf(stderr, "err = %d\n", err);
#endif /* MOO */
        set_free(LL);
        set_free(RR);
        set_free(SS);
        ddf_clear(value);
        return err;
    }

    if (foo) {

        set_uni(M->linset, M->linset, R);

        err = ddf_NoError;
        ddf_FindRelativeInterior(M, &ImL, &Lbasis, &lps, &err);
        if (err != ddf_NoError) {
#ifdef MOO
            fprintf(stderr, "err from ddf_FindRelativeInterior\n");
            fprintf(stderr, "err = %d\n", err);
#endif /* MOO */
            set_free(LL);
            set_free(RR);
            set_free(SS);
            ddf_clear(value);
            return err;
        }

        dim = M->colsize - set_card(Lbasis) - 1;
        set_uni(M->linset, M->linset, ImL);

        SEXP mydim, myactive, myrip;
        PROTECT(mydim = ScalarInteger(dim));
        PROTECT(myactive = rrf_set_fwrite(M->linset));
        int myd = (lps->d) - 2;
        PROTECT(myrip = allocVector(REALSXP, myd));
        for (int j = 1; j <= myd; j++) {
            double ax = ddf_get_d(lps->sol[j]);
            REAL(myrip)[j - 1] = ax;
        }
        REPROTECT(dimlist = CONS(mydim, dimlist), dimidx);
        REPROTECT(riplist = CONS(myrip, riplist), ripidx);
        REPROTECT(activelist = CONS(myactive, activelist), activeidx);
        UNPROTECT(3);

        ddf_FreeLPSolution(lps);
        set_free(ImL);
        set_free(Lbasis);
   
        if (dim > 0) {
            for (int i = 1; i <= M->rowsize; i++) {
                if ((! set_member(i, M->linset)) && (! set_member(i, S))) {
                    set_addelem(RR, i);
                    if (iprev) {
                        set_delelem(RR, iprev);
                        set_delelem(M->linset, iprev);
                        set_addelem(SS, iprev);
                    }
                    iprev = i;

                    err = FaceEnumHelper(M, RR, SS);
                    if (err != ddf_NoError) {
#ifdef MOO
                        fprintf(stderr, "err from FaceEnumHelper\n");
                        fprintf(stderr, "err = %d\n", err);
#endif /* MOO */
                        set_copy(M->linset, LL);
                        set_free(LL);
                        set_free(RR);
                        set_free(SS);
                        ddf_clear(value);
                        return err;
                    }
                }
            }
        }
    }

    set_copy(M->linset, LL);
    set_free(LL);
    set_free(RR);
    set_free(SS);
    ddf_clear(value);
    return ddf_NoError;
}

