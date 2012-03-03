
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
#include <string.h>
#ifdef MOO
#include <stdio.h>
#endif /* MOO */
#include "rcdd.h"

static SEXP FaceEnum(dd_MatrixPtr M);

SEXP allfaces(SEXP hrep)
{
    if (! isMatrix(hrep))
        error("'hrep' must be matrix");
    if (! isString(hrep))
        error("'hrep' must be character");

    SEXP hrep_dim;
    PROTECT(hrep_dim = getAttrib(hrep, R_DimSymbol));
    int nrow = INTEGER(hrep_dim)[0];
    int ncol = INTEGER(hrep_dim)[1];
    UNPROTECT(1);

    if (nrow <= 0)
        error("no rows in 'hrep'");
    if (ncol <= 3)
        error("three or fewer cols in hrep");

    for (int i = 0; i < nrow; ++i) {
        char *foo = (char *) CHAR(STRING_ELT(hrep, i));
        if (strlen(foo) != 1)
            error("column one of 'hrep' not zero-or-one valued");
        if (! (foo[0] == '0' || foo[0] == '1'))
            error("column one of 'hrep' not zero-or-one valued");
    }

    dd_set_global_constants();

    /* note actual type of "value" is mpq_t (defined in cddmp.h) */
    mytype value;
    dd_init(value);

    dd_MatrixPtr mf = dd_CreateMatrix(nrow, ncol - 1);
    /* note our matrix has one more column than Fukuda's */

    mf->representation = dd_Inequality;
    mf->numbtype = dd_Rational;

    /* linearity */
    for (int i = 0; i < nrow; ++i) {
        char *foo = (char *) CHAR(STRING_ELT(hrep, i));
        if (foo[0] == '1')
            set_addelem(mf->linset, i + 1);
        /* note conversion from zero-origin to one-origin indexing */
    }

    /* matrix */
    for (int j = 1, k = nrow; j < ncol; ++j)
        for (int i = 0; i < nrow; ++i, ++k) {
            char *rat_str = (char *) CHAR(STRING_ELT(hrep, k));
            if (mpq_set_str(value, rat_str, 10) == -1)
                error("error converting string to GMP rational");
            mpq_canonicalize(value);
            dd_set(mf->matrix[i][j - 1], value);
            /* note our matrix has one more column than Fukuda's */
        }

    SEXP result;
    PROTECT(result = FaceEnum(mf));

    dd_FreeMatrix(mf);
    dd_clear(value);
    dd_free_global_constants();

    UNPROTECT(1);
    if (result == R_NilValue)
        error("failed");
    return result;
}

static dd_ErrorType FaceEnumHelper(dd_MatrixPtr M, dd_rowset R, dd_rowset S);

static SEXP dimlist, riplist, activelist;

static PROTECT_INDEX dimidx, ripidx, activeidx;

static SEXP FaceEnum(dd_MatrixPtr M)
{
    PROTECT_WITH_INDEX(dimlist = R_NilValue, &dimidx);
    PROTECT_WITH_INDEX(riplist = R_NilValue, &ripidx);
    PROTECT_WITH_INDEX(activelist = R_NilValue, &activeidx);

    dd_rowset R, S;
    set_initialize(&R, M->rowsize);
    set_initialize(&S, M->rowsize);

    dd_ErrorType err = FaceEnumHelper(M, R, S);

    set_free(R);
    set_free(S);

    if (err != dd_NoError) {
#ifdef MOO
        switch (err) {
            case dd_DimensionTooLarge:
                fprintf(stderr, "err = dd_DimensionTooLarge\n");
                break;
            case dd_ImproperInputFormat:
                fprintf(stderr, "err = dd_ImproperInputFormat\n");
                break;
            case dd_NegativeMatrixSize:
                fprintf(stderr, "err = dd_NegativeMatrixSize\n");
                break;
            case dd_EmptyVrepresentation:
                fprintf(stderr, "err = dd_EmptyVrepresentation\n");
                break;
            case dd_EmptyHrepresentation:
                fprintf(stderr, "err = dd_EmptyHrepresentation\n");
                break;
            case dd_EmptyRepresentation:
                fprintf(stderr, "err = dd_EmptyRepresentation\n");
                break;
            case dd_IFileNotFound:
                fprintf(stderr, "err = dd_IFileNotFound\n");
                break;
            case dd_OFileNotOpen:
                fprintf(stderr, "err = dd_OFileNotOpen\n");
                break;
            case dd_NoLPObjective:
                fprintf(stderr, "err = dd_NoLPObjective\n");
                break;
            case dd_NoRealNumberSupport:
                fprintf(stderr, "err = dd_NoRealNumberSupport\n");
                break;
            case dd_NotAvailForH:
                fprintf(stderr, "err = dd_NotAvailForH\n");
                break;
            case dd_NotAvailForV:
                fprintf(stderr, "err = dd_NotAvailForV\n");
                break;
            case dd_CannotHandleLinearity:
                fprintf(stderr, "err = dd_CannotHandleLinearity\n");
                break;
            case dd_RowIndexOutOfRange:
                fprintf(stderr, "err = dd_RowIndexOutOfRange\n");
                break;
            case dd_ColIndexOutOfRange:
                fprintf(stderr, "err = dd_ColIndexOutOfRange\n");
                break;
            case dd_LPCycling:
                fprintf(stderr, "err = dd_LPCycling\n");
                break;
            case dd_NumericallyInconsistent:
                fprintf(stderr, "err = dd_NumericallyInconsistent\n");
                break;
            case dd_NoError:
                fprintf(stderr, "err = dd_NoError\n");
                break;
            default:
                fprintf(stderr, "err bogus, WTF????\n");
        }
#endif /* MOO */
        rr_WriteErrorMessages(err);
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

static dd_ErrorType FaceEnumHelper(dd_MatrixPtr M, dd_rowset R, dd_rowset S)
{
    dd_ErrorType err;
    dd_rowset LL, ImL, RR, SS, Lbasis;
    dd_rowrange iprev = 0;
    dd_colrange dim;
    dd_LPSolutionPtr lps = NULL;

    set_initialize(&LL, M->rowsize);
    set_initialize(&RR, M->rowsize);
    set_initialize(&SS, M->rowsize);
    set_copy(LL, M->linset);
    set_copy(RR, R);
    set_copy(SS, S);

    /* note actual type of "value" is mpq_t (defined in cddmp.h) */
    mytype value;
    dd_init(value);

    err = dd_NoError;
    dd_boolean foo = dd_ExistsRestrictedFace(M, R, S, &err);
    if (err != dd_NoError) {
#ifdef MOO
        fprintf(stderr, "err from dd_ExistsRestrictedFace\n");
        fprintf(stderr, "err = %d\n", err);
#endif /* MOO */
        set_free(LL);
        set_free(RR);
        set_free(SS);
        dd_clear(value);
        return err;
    }

    if (foo) {

        set_uni(M->linset, M->linset, R);

        err = dd_NoError;
        dd_FindRelativeInterior(M, &ImL, &Lbasis, &lps, &err);
        if (err != dd_NoError) {
#ifdef MOO
            fprintf(stderr, "err from dd_FindRelativeInterior\n");
            fprintf(stderr, "err = %d\n", err);
#endif /* MOO */
            set_free(LL);
            set_free(RR);
            set_free(SS);
            dd_clear(value);
            return err;
        }

        dim = M->colsize - set_card(Lbasis) - 1;
        set_uni(M->linset, M->linset, ImL);

        SEXP mydim, myactive, myrip;
        PROTECT(mydim = ScalarInteger(dim));
        PROTECT(myactive = rr_set_fwrite(M->linset));
        int myd = (lps->d) - 2;
        PROTECT(myrip = allocVector(STRSXP, myd));
        for (int j = 1; j <= myd; j++) {
            dd_set(value, lps->sol[j]);
            char *zstr = NULL;
            zstr = mpq_get_str(zstr, 10, value);
            SET_STRING_ELT(myrip, j - 1, mkChar(zstr));
            free(zstr);
        }
        REPROTECT(dimlist = CONS(mydim, dimlist), dimidx);
        REPROTECT(riplist = CONS(myrip, riplist), ripidx);
        REPROTECT(activelist = CONS(myactive, activelist), activeidx);
        UNPROTECT(3);

        dd_FreeLPSolution(lps);
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
                    if (err != dd_NoError) {
#ifdef MOO
                        fprintf(stderr, "err from FaceEnumHelper\n");
                        fprintf(stderr, "err = %d\n", err);
#endif /* MOO */
                        set_copy(M->linset, LL);
                        set_free(LL);
                        set_free(RR);
                        set_free(SS);
                        dd_clear(value);
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
    dd_clear(value);
    return dd_NoError;
}

