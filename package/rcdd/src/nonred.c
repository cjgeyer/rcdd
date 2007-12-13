
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
#include <string.h>

static void hash_setup(int pow2);
static int my_subset(SEXP set1, SEXP set2);

SEXP nonred(SEXP sets, SEXP pow2)
{
    int i, j;

    if (! isNewList(sets))
        error("argument not list");
    if (! isInteger(pow2))
        error("'pow2' not integer");
    if (LENGTH(pow2) > 1)
        error("'pow2' not scalar");

    int n = LENGTH(sets);
    int my_pow2 = INTEGER(pow2)[0];

#ifdef BLEAT
    printf("n = %d\n", n);
    printf("pow2 = %d\n", my_pow2);
#endif /* BLEAT */

    hash_setup(my_pow2);

    for (i = 0; i < n; ++i) {
        SEXP foo = VECTOR_ELT(sets, i);
        if (! isInteger(foo))
            error("argument not list of integer vectors");
        int nfoo = LENGTH(foo);
        for (j = 0; j < nfoo; ++j)
            if (INTEGER(foo)[j] <= 0)
                error("argument not list of positive integer vectors");
    }

    SEXP result;
    PROTECT(result = allocVector(LGLSXP, n));
    for (i = 0; i < n; ++i)
        LOGICAL(result)[i] = TRUE;

    for (i = 0; i < n; ++i) {
        SEXP set1 = VECTOR_ELT(sets, i);
        for (j = 0; j < n; ++j)
            if ((i != j) && LOGICAL(result)[j]) {
                SEXP set2 = VECTOR_ELT(sets, j);
                if (my_subset(set1, set2)) {
                    LOGICAL(result)[i] = FALSE;
                    break;
                }
            }
    }

    UNPROTECT(1);
    return result;
}

static int hashsize;
static unsigned int hashmask;
static int *hashtable;
#define MAX_COLLISION 20

static void hash_setup(int pow2)
{
     hashsize = 1 << pow2;
     hashmask = hashsize - 1;
     hashtable = (int *) R_alloc(hashsize, sizeof(int));
#ifdef BLEAT
    printf("pow2 = %d (decimal)\n", pow2);
    printf("hashsize = %d (decimal) %X (hex)\n", hashsize, hashsize);
    printf("hashmask = %X (hex)\n", hashmask);
    printf("MAX_COLLISION = %d (decimal)\n", MAX_COLLISION);
#endif /* BLEAT */
}

static void hash_clear()
{
    memset(hashtable, -1, hashsize * sizeof(int));
}

static int hash_insert_find(int foo, int insert)
{
    /*
    *  insert (if insert == TRUE) or find (if insert == FALSE)
    *       foo in hashtable
    *  return TRUE if foo was already there
    *  return FALSE otherwise
    */

    unsigned int hash = (2654435761u * foo) & hashmask;
    int collisions = 0;
    while (TRUE) {
        if (collisions > MAX_COLLISION)
            error ("too many collisions in hash table, increase table size");
        if (hashtable[hash] == foo)
            return TRUE;
        if (hashtable[hash] < 0)
            break;
        ++hash;
        hash &= hashmask;
        ++collisions;
    }
    if (insert)
        hashtable[hash] = foo;
    return FALSE;
}

static int my_subset(SEXP set1, SEXP set2)
{
    /*
    *  we've already checked that set1 and set2 are nonnegative integer vectors
    *  return TRUE if set1 is subset of set2
    *  return FALSE otherwise
    */

    int i;

    int n1 = LENGTH(set1);
    int n2 = LENGTH(set2);

    if (n1 == 0)
        return TRUE;
    if (n2 == 0)
        return FALSE;

    hash_clear();

    for (i = 0; i < n2; ++i) {
        int fred = INTEGER(set2)[i];
        hash_insert_find(fred, TRUE);
    }

    /* check subset relation */
    for (i = 0; i < n1; ++i) {
        int fred = INTEGER(set1)[i];
        if (!  hash_insert_find(fred, FALSE)) {
            return FALSE;
        }
    }
    return TRUE;
}

SEXP test_my_subset(SEXP set1, SEXP set2, SEXP pow2)
{
    int j;

    if (! isInteger(set1))
        error("'set1' not integer");
    if (! isInteger(set2))
        error("'set2' not integer");
    if (! isInteger(pow2))
        error("'pow2' not integer");
    if (LENGTH(pow2) > 1)
        error("'pow2' not scalar");

    hash_setup(INTEGER(pow2)[0]);

    for (j = 0; j < LENGTH(set1); ++j)
        if (INTEGER(set1)[j] <= 0)
            error("'set1' not positive");
    for (j = 0; j < LENGTH(set2); ++j)
        if (INTEGER(set2)[j] <= 0)
            error("'set2' not positive");

    if (my_subset(set1, set2))
        return Rf_ScalarLogical(TRUE);
    else
        return Rf_ScalarLogical(FALSE);
}

SEXP all_intersect(SEXP sets, SEXP pow2)
{
    int i, j, k, m, mtoo;

    if (! isNewList(sets))
        error("argument not list");
    if (! isInteger(pow2))
        error("'pow2' not integer");
    if (LENGTH(pow2) > 1)
        error("'pow2' not scalar");

    int n = LENGTH(sets);
    int my_pow2 = INTEGER(pow2)[0];

#ifdef BLEAT
    printf("n = %d\n", n);
    printf("pow2 = %d\n", my_pow2);
#endif /* BLEAT */

    hash_setup(my_pow2);

    for (i = 0; i < n; ++i) {
        SEXP foo = VECTOR_ELT(sets, i);
        if (! isInteger(foo))
            error("argument not list of integer vectors");
        int nfoo = LENGTH(foo);
        for (j = 0; j < nfoo; ++j)
            if (INTEGER(foo)[j] <= 0)
                error("argument not list of positive integer vectors");
    }

    SEXP result;
    PROTECT(result = allocVector(VECSXP, n * (n - 1) / 2));

    for (i = 0, k = 0; i < n; ++i) {
        SEXP foo = VECTOR_ELT(sets, i);
        hash_clear();
        int nfoo = LENGTH(foo);
        for (j = 0; j < nfoo; ++j)
            hash_insert_find(INTEGER(foo)[j], TRUE);
        for (j = i + 1; j < n; ++j, ++k) {
            /* do the intersection of i and j and put in result k */
            /* use sign as flag */
            SEXP bar = VECTOR_ELT(sets, j);
            int nbar = LENGTH(bar);
            int ninter = 0;
            for (m = 0; m < nbar; ++m)
                if (hash_insert_find(INTEGER(bar)[m], FALSE)) {
                    INTEGER(bar)[m] *= -1;
                    ++ninter;
                }
            SET_VECTOR_ELT(result, k, allocVector(INTSXP, ninter));
            SEXP qux = VECTOR_ELT(result, k);
            for (m = 0, mtoo = 0; m < nbar; ++m)
                if (INTEGER(bar)[m] < 0) {
                    INTEGER(bar)[m] *= -1;
                    INTEGER(qux)[mtoo] = INTEGER(bar)[m];
                    ++mtoo;
                }
        }
    }

    UNPROTECT(1);
    return result;
}

SEXP all_union(SEXP sets, SEXP pow2)
{
    int i, j, k, m, mtoo;

    if (! isNewList(sets))
        error("argument not list");
    if (! isInteger(pow2))
        error("'pow2' not integer");
    if (LENGTH(pow2) > 1)
        error("'pow2' not scalar");

    int n = LENGTH(sets);
    int my_pow2 = INTEGER(pow2)[0];

#ifdef BLEAT
    printf("n = %d\n", n);
    printf("pow2 = %d\n", my_pow2);
#endif /* BLEAT */

    hash_setup(my_pow2);

    for (i = 0; i < n; ++i) {
        SEXP foo = VECTOR_ELT(sets, i);
        if (! isInteger(foo))
            error("argument not list of integer vectors");
        int nfoo = LENGTH(foo);
        for (j = 0; j < nfoo; ++j)
            if (INTEGER(foo)[j] <= 0)
                error("argument not list of positive integer vectors");
    }

    SEXP result;
    PROTECT(result = allocVector(VECSXP, n * (n - 1) / 2));

    for (i = 0, k = 0; i < n; ++i) {
        SEXP foo = VECTOR_ELT(sets, i);
        hash_clear();
        int nfoo = LENGTH(foo);
        for (j = 0; j < nfoo; ++j)
            hash_insert_find(INTEGER(foo)[j], TRUE);
        for (j = i + 1; j < n; ++j, ++k) {
            /* do the union of i and j and put in result k */
            /* use sign as flag */
            SEXP bar = VECTOR_ELT(sets, j);
            int nbar = LENGTH(bar);
            int ninter = 0;
            for (m = 0; m < nbar; ++m)
                if (hash_insert_find(INTEGER(bar)[m], FALSE)) {
                    INTEGER(bar)[m] *= -1;
                    ++ninter;
                }
            SET_VECTOR_ELT(result, k, allocVector(INTSXP,
                nfoo + nbar - ninter));
            SEXP qux = VECTOR_ELT(result, k);
            for (m = 0; m < nfoo; ++m)
                INTEGER(qux)[m] = INTEGER(foo)[m];
            for (m = 0, mtoo = 0; m < nbar; ++m)
                if (INTEGER(bar)[m] < 0) {
                    INTEGER(bar)[m] *= -1;
                } else {
                    INTEGER(qux)[nfoo + mtoo] = INTEGER(bar)[m];
                    ++mtoo;
                }
        }
    }

    UNPROTECT(1);
    return result;
}

