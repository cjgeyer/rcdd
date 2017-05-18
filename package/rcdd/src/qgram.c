
/* rational gram-schmidt */

#include <Rinternals.h>
#include <gmp.h>
#include <stdlib.h>
#include "rcdd.h"

SEXP qgram(SEXP foo)
{
    if ((! isString(foo)))
        error("argument must be character");
    if ((! isMatrix(foo)))
        error("argument must be matrix");

    SEXP dim;
    PROTECT(dim = getAttrib(foo, R_DimSymbol));
    int nrow_foo = INTEGER(dim)[0];
    int ncol_foo = INTEGER(dim)[1];
    UNPROTECT(1);

    if (nrow_foo <= 0)
        error("row dimension of arg must be positive");
    if (ncol_foo <= 0)
        error("col dimension of arg must be positive");

    SEXP baz;
    PROTECT(baz = allocMatrix(STRSXP, nrow_foo, ncol_foo));

    mpq_t *bar = (mpq_t *) R_alloc(nrow_foo * ncol_foo, sizeof(mpq_t));

    for (int i = 0; i < nrow_foo * ncol_foo; i++)
        mpq_init(bar[i]);
    for (int i = 0; i < nrow_foo * ncol_foo; i++) {
        const char *foo_i = CHAR(STRING_ELT(foo, i));
        if (mpq_set_str(bar[i], foo_i, 10) == -1) {
            for (int j = 0; j < nrow_foo * ncol_foo; j++)
            mpq_clear(bar[j]);
            error("error converting string to GMP rational");
        }
        mpq_canonicalize(bar[i]);
    }

    for (int j = 0; j < ncol_foo; j++) {
        int jbase = nrow_foo * j;
        for (int k = 0; k < j; k++) {
            int kbase = nrow_foo * k;

            /* orthogonalize vector j to vector k */
            mpq_t dot_uu, dot_uv;
            mpq_init(dot_uu);
            mpq_init(dot_uv);

            for (int i = 0; i < nrow_foo; i++) {
                 mpq_t tmp;
                 mpq_init(tmp);
                 mpq_mul(tmp, bar[kbase + i], bar[kbase + i]);
                 mpq_add(dot_uu, dot_uu, tmp);
                 mpq_clear(tmp);
            }

            for (int i = 0; i < nrow_foo; i++) {
                 mpq_t tmp;
                 mpq_init(tmp);
                 mpq_mul(tmp, bar[jbase + i], bar[kbase + i]);
                 mpq_add(dot_uv, dot_uv, tmp);
                 mpq_clear(tmp);
            }

            if (mpq_sgn(dot_uu) != 0)
                 mpq_div(dot_uv, dot_uv, dot_uu);
            /* at this point dot_uv is <u, v> / <u, u> */
            /* where u is vector k and v is vector j */

            for (int i = 0; i < nrow_foo; i++) {
                 mpq_t tmp;
                 mpq_init(tmp);
                 mpq_mul(tmp, dot_uv, bar[kbase + i]);
                 mpq_sub(bar[jbase + i], bar[jbase + i], tmp);
                 mpq_clear(tmp);
            }

            mpq_clear(dot_uu);
            mpq_clear(dot_uv);
        }
    }

    /* L1 normalize each vector (not L2 because that would give irrational) */

    for (int j = 0; j < ncol_foo; j++) {
        int jbase = nrow_foo * j;
        mpq_t norm_u;
        mpq_init(norm_u);
        for (int i = 0; i < nrow_foo; i++) {
            mpq_t tmp;
            mpq_init(tmp);
            mpq_abs(tmp, bar[jbase + i]);
            mpq_add(norm_u, norm_u, tmp);
            mpq_clear(tmp);
        }

        if (mpq_sgn(norm_u) != 0)
            for (int i = 0; i < nrow_foo; i++)
                mpq_div(bar[jbase + i], bar[jbase + i], norm_u);

        mpq_clear(norm_u);
    }

    for (int i = 0; i < nrow_foo * ncol_foo; i++) {
        char *baz_i = mpq_get_str(NULL, 10, bar[i]);
        SET_STRING_ELT(baz, i, mkChar(baz_i));
        free(baz_i);
        mpq_clear(bar[i]);
    }

    UNPROTECT(1);
    return(baz);
}

