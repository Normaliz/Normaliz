/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of e-antic

    e-antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>

#include "../e-antic/fmpz_poly_extra.h"

/*

void fmpz_poly_randtest_irreducible(fmpz_poly_t p, flint_rand_t state, slong len, mp_bitcnt_t bits)
{
    slong i;
    fmpz_t c;
    fmpz_mod_poly_t q;
#if __FLINT_RELEASE >= 20700
    fmpz_mod_ctx_t ctx;
#endif

    fmpz_init(c);

    fmpz_randprime(c, state, bits, 0);
#if __FLINT_RELEASE >= 20700
    fmpz_mod_ctx_init(ctx, c);
    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_randtest_irreducible(q, state, len, ctx);

    fmpz_mod_poly_get_fmpz_poly(p, q, ctx);
#else
    fmpz_mod_poly_init(q, c);
    fmpz_mod_poly_randtest_irreducible(q, state, len);

    fmpz_mod_poly_get_fmpz_poly(p, q);
#endif
*/
    /* After lifting, the coefficients belong to {0, ..., c-1}. We now  */
    /* randomly subtract c so that some of them become negative.        */
/*    for (i = 0; i < p->length; i++)
    {
        if (n_randint(state, 3) == 0)
            fmpz_sub(
                fmpz_poly_get_coeff_ptr(p, i),
                fmpz_poly_get_coeff_ptr(p, i),
                c);
    }

#if __FLINT_RELEASE >= 20700
    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_ctx_clear(ctx);
#else
    fmpz_mod_poly_clear(q);
#endif
    fmpz_clear(c);
}
*/
