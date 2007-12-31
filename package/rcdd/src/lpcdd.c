
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
#ifdef WOOF
#include <stdio.h>
#endif /* WOOF */
#include "rcdd.h"

SEXP lpcdd(SEXP hrep, SEXP objfun, SEXP minimize, SEXP solver)
{
    if (! isMatrix(hrep))
        error("'hrep' must be matrix");
    if (! isString(hrep))
        error("'hrep' must be character");
    if (! isString(objfun))
        error("'objfun' must be character");
    if (! isLogical(minimize))
        error("'minimize' must be logical");
    if (! isString(solver))
        error("'solver' must be character");

    if (LENGTH(minimize) != 1)
        error("'minimize' must be scalar");
    if (LENGTH(solver) != 1)
        error("'solver' must be scalar");

    SEXP hrep_dim;
    PROTECT(hrep_dim = getAttrib(hrep, R_DimSymbol));
    int nrow = INTEGER(hrep_dim)[0];
    int ncol = INTEGER(hrep_dim)[1];
    UNPROTECT(1);

    if (nrow <= 0)
        error("no rows in 'hrep'");
    if (ncol <= 2)
        error("no cols in hrep[ , - c(1, 2)]");
    if (LENGTH(objfun) != ncol - 1)
        error("length(objfun) != ncol(hrep) - 1");

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

    /* rowvec */
    for (int j = 0; j < ncol - 1; ++j) {
        char *rat_str = (char *) CHAR(STRING_ELT(objfun, j));
        if (mpq_set_str(value, rat_str, 10) == -1)
            error("error converting string to GMP rational");
        mpq_canonicalize(value);
        dd_set(mf->rowvec[j], value);
    }

    if(LOGICAL(minimize)[0])
        mf->objective = dd_LPmin;
    else
        mf->objective = dd_LPmax;

    dd_ErrorType err = dd_NoError;
    dd_LPPtr lp = dd_Matrix2LP(mf, &err);

    if (err != dd_NoError) {
        dd_WriteErrorMessages(stderr, err);
        dd_FreeLPData(lp);
        dd_FreeMatrix(mf);
        dd_clear(value);
        dd_free_global_constants();
        error("failed");
    }

    dd_LPSolverType sol = dd_DualSimplex;
    char *sol_str = (char *) CHAR(STRING_ELT(solver, 0));
    if(strcmp(sol_str, "DualSimplex") == 0)
        sol = dd_DualSimplex;
    else if(strcmp(sol_str, "CrissCross") == 0)
        sol = dd_CrissCross;
    else
        error("solver not recognized");

#ifdef WOOF
    dd_WriteLP(stderr, lp);
#endif /* WOOF */


    dd_LPSolve(lp, sol, &err); 

    if (err != dd_NoError) {
        dd_WriteErrorMessages(stderr, err);
        dd_FreeLPData(lp);
        dd_FreeMatrix(mf);
        dd_clear(value);
        dd_free_global_constants();
        error("failed");
    }

#ifdef WOOF
    switch (lp->LPS) {
        case dd_LPSundecided:
            fprintf(stderr, "    dd_LPSundecided\n");
            break;
        case dd_Optimal:
            fprintf(stderr, "    dd_Optimal\n");
            break;
        case dd_Inconsistent:
            fprintf(stderr, "    dd_Inconsistent\n");
            break;
        case dd_DualInconsistent:
            fprintf(stderr, "    dd_DualInconsistent\n");
            break;
        case dd_StrucInconsistent:
            fprintf(stderr, "    dd_StrucInconsistent\n");
            break;
        case dd_StrucDualInconsistent:
            fprintf(stderr, "    dd_StrucDualInconsistent\n");
            break;
        case dd_Unbounded:
            fprintf(stderr, "    dd_Unbounded\n");
            break;
        case dd_DualUnbounded:
            fprintf(stderr, "    dd_DualUnbounded\n");
            break;
        default:
            fprintf(stderr, "    default  WTF???\n");
    }
#endif /* WOOF */

    SEXP result;
    if (lp->LPS == dd_Optimal) {
        SEXP resultnames;
        PROTECT(result = allocVector(VECSXP, 4));
        PROTECT(resultnames = allocVector(STRSXP, 4));

        SET_STRING_ELT(resultnames, 0, mkChar("solution.type"));
        SET_STRING_ELT(resultnames, 1, mkChar("primal.solution"));
        SET_STRING_ELT(resultnames, 2, mkChar("dual.solution"));
        SET_STRING_ELT(resultnames, 3, mkChar("optimal.value"));
        namesgets(result, resultnames);

        SEXP foo, bar, baz;
        SET_VECTOR_ELT(result, 0, ScalarString(mkChar("Optimal")));
        PROTECT(foo = allocVector(STRSXP, ncol - 2));
        PROTECT(bar = allocVector(STRSXP, nrow));
        PROTECT(baz = allocVector(STRSXP, 1));
        SET_VECTOR_ELT(result, 1, foo);
        SET_VECTOR_ELT(result, 2, bar);
        SET_VECTOR_ELT(result, 3, baz);

        if (lp->d != ncol - 1)
            error("Can't happen.  Dimension changed.");
        for (int j = 1; j < ncol - 1; j++) {
            dd_set(value, lp->sol[j]);
            char *zstr = NULL;
            zstr = mpq_get_str(zstr, 10, value);
            SET_STRING_ELT(foo, j - 1, mkChar(zstr));
            free(zstr);
        }
        for (int j = 0; j < nrow; j++) {
            SET_STRING_ELT(bar, j, mkChar("0"));
        }
        for (int j = 1; j < ncol - 1; j++) {
            int qux = lp->nbindex[j + 1];
            if (qux > 0) {
                if (! (1 <= qux && qux <= nrow))
                    error("Can't happen.  Dual solution index out of bounds");
                dd_set(value, lp->dsol[j]);
                char *zstr = NULL;
                zstr = mpq_get_str(zstr, 10, value);
                SET_STRING_ELT(bar, qux - 1, mkChar(zstr));
                free(zstr);
            }
        }
        {
            dd_set(value, lp->optvalue);
            char *zstr = NULL;
            zstr = mpq_get_str(zstr, 10, value);
            SET_STRING_ELT(baz, 0, mkChar(zstr));
            free(zstr);
        }

        UNPROTECT(4);
    } else if (lp->LPS == dd_Inconsistent) {
        SEXP resultnames;
        PROTECT(result = allocVector(VECSXP, 2));
        PROTECT(resultnames = allocVector(STRSXP, 2));

        SET_STRING_ELT(resultnames, 0, mkChar("solution.type"));
        SET_STRING_ELT(resultnames, 1, mkChar("dual.direction"));
        namesgets(result, resultnames);

        SEXP foo;
        SET_VECTOR_ELT(result, 0, ScalarString(mkChar("Inconsistent")));
        PROTECT(foo = allocVector(STRSXP, nrow));
        SET_VECTOR_ELT(result, 1, foo);

        if (lp->d != ncol - 1)
            error("Can't happen.  Dimension changed.");
        for (int j = 0; j < nrow; j++) {
            SET_STRING_ELT(foo, j, mkChar("0"));
        }
        for (int j = 1; j < ncol - 1; j++) {
            int qux = lp->nbindex[j + 1];
            if (qux > 0) {
                if (! (1 <= qux && qux <= nrow))
                    error("Can't happen.  Dual solution index out of bounds");
                dd_set(value, lp->dsol[j]);
                char *zstr = NULL;
                zstr = mpq_get_str(zstr, 10, value);
                SET_STRING_ELT(foo, qux - 1, mkChar(zstr));
                free(zstr);
            }
        }
        {
            int qux = lp->re;
            if (! (1 <= qux && qux <= nrow))
                error("Can't happen.  Dual solution index out of bounds");
            SET_STRING_ELT(foo, qux - 1, mkChar("1"));
        }

        UNPROTECT(2);
    } else if (lp->LPS == dd_DualInconsistent ||
        lp->LPS == dd_StrucDualInconsistent) {
        SEXP resultnames;
        PROTECT(result = allocVector(VECSXP, 2));
        PROTECT(resultnames = allocVector(STRSXP, 2));

        SET_STRING_ELT(resultnames, 0, mkChar("solution.type"));
        SET_STRING_ELT(resultnames, 1, mkChar("primal.direction"));
        namesgets(result, resultnames);

        SEXP foo;
        if (lp->LPS == dd_DualInconsistent)
            SET_VECTOR_ELT(result, 0, ScalarString(mkChar("DualInconsistent")));
        else
            SET_VECTOR_ELT(result, 0,
                ScalarString(mkChar("StrucDualInconsistent")));
        PROTECT(foo = allocVector(STRSXP, ncol - 2));
        SET_VECTOR_ELT(result, 1, foo);

        if (lp->d != ncol - 1)
            error("Can't happen.  Dimension changed.");
        for (int j = 1; j < ncol - 1; j++) {
            dd_set(value, lp->sol[j]);
            char *zstr = NULL;
            zstr = mpq_get_str(zstr, 10, value);
            SET_STRING_ELT(foo, j - 1, mkChar(zstr));
            free(zstr);
        }

        UNPROTECT(2);
    } else {
        SEXP resultnames;
        PROTECT(result = allocVector(VECSXP, 1));
        PROTECT(resultnames = allocVector(STRSXP, 1));

        SET_STRING_ELT(resultnames, 0, mkChar("solution.type"));
        namesgets(result, resultnames);

        switch (lp->LPS) {
            case dd_LPSundecided:
                SET_VECTOR_ELT(result, 0,
                    ScalarString(mkChar("LPSundecided")));
                break;
            case dd_StrucInconsistent:
                SET_VECTOR_ELT(result, 0,
                    ScalarString(mkChar("StrucInconsistent")));
                break;
            case dd_Unbounded:
                SET_VECTOR_ELT(result, 0,
                    ScalarString(mkChar("Unbounded")));
                break;
            case dd_DualUnbounded:
                SET_VECTOR_ELT(result, 0,
                    ScalarString(mkChar("DualUnbounded")));
                break;
            default:
                error("unrecognized solution type");
        }

        UNPROTECT(1);
    }

    dd_FreeLPData(lp);
    dd_FreeMatrix(mf);
    dd_clear(value);
    dd_free_global_constants();

    UNPROTECT(1);
    return result;
}

