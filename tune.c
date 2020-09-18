/*
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of the tune-params code repository.

    tune-params is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
check_arb_accuracy(const arb_t x, const char *name)
{
    slong required_bits = 16;
    slong bits = arb_rel_accuracy_bits(x);
    if (bits < required_bits)
    {
        flint_printf("%s has only %ld relative accuracy bits\n",
                name, bits);
        flint_abort();
    }
}

static void
_arb_div_si_si(arb_t res, slong a, slong b, slong prec)
{
    arb_set_si(res, a);
    arb_div_si(res, res, b, prec);
}

static void
_arb_inv_si(arb_t res, slong n, slong prec)
{
    arb_set_si(res, n);
    arb_inv(res, res, prec);
}

static void
_arb_get_lbound_fmpz(fmpz_t z, const arb_t x, slong prec)
{
    arf_t u;
    arf_init(u);
    arb_get_lbound_arf(u, x, prec);
    arf_get_fmpz(z, u, ARF_RND_DOWN);
    arf_clear(u);
}

/* Aggregate parameter values for multieval error calculations. */
typedef struct
{
    slong A;
    slong B;
    slong J;
    slong K;
    slong sigma_grid;
    slong sigma_interp;
    slong prec;
    slong hnum;
    slong hden;
    slong Hnum;
    slong Hden;
    slong Ns_max;
    fmpz T;
    arb_struct t0;
}
agg_struct;
typedef agg_struct agg_t[1];
typedef agg_struct * agg_ptr;
typedef const agg_struct * agg_srcptr;

static void
agg_init(agg_t agg)
{
    arb_init(&agg->t0);
    fmpz_init(&agg->T);
}

static void
agg_clear(agg_t agg)
{
    arb_clear(&agg->t0);
    fmpz_clear(&agg->T);
}

static void
agg_set(agg_t res, const agg_t p)
{
    res->A = p->A;
    res->B = p->B;
    res->J = p->J;
    res->K = p->K;
    res->sigma_grid = p->sigma_grid;
    res->sigma_interp = p->sigma_interp;
    res->prec = p->prec;
    res->hnum = p->hnum;
    res->hden = p->hden;
    res->Hnum = p->Hnum;
    res->Hden = p->Hden;
    res->Ns_max = p->Ns_max;
    arb_set(&res->t0, &p->t0);
}

static void
agg_printf(const agg_t p)
{
    flint_printf("A = %ld\n", p->A);
    flint_printf("B = %ld\n", p->B);
    flint_printf("J = %ld\n", p->J);
    flint_printf("K = %ld\n", p->K);
    flint_printf("sigma_grid = %ld\n", p->sigma_grid);
    flint_printf("sigma_interp = %ld\n", p->sigma_interp);
    flint_printf("prec (for error estimation) = %ld\n", p->prec);
    flint_printf("prec (for zeros computation) = %ld\n",
            p->prec + fmpz_sizeinbase(&p->T, 2));
    flint_printf("h = %ld/%ld = %lf\n",
            p->hnum, p->hden, p->hnum / (double) p->hden);
    flint_printf("H = %ld/%ld = %lf\n",
            p->Hnum, p->Hden, p->Hnum / (double) p->Hden);
    flint_printf("Ns_max = %ld\n", p->Ns_max);
    flint_printf("T = "); fmpz_print(&p->T); flint_printf("\n");
}

/* approximate interpolation error */
static int
_get_interpolation_error(arb_t res, const agg_t p)
{
    arb_t c1, c3, H;

    arb_init(c1);
    arb_init(c3);
    arb_init(H);

    _arb_div_si_si(H, p->Hnum, p->Hden, p->prec);
    acb_dirichlet_platt_i_bound(
            c1, &p->t0, p->A, H, p->sigma_interp, p->prec);
    acb_dirichlet_platt_bound_C3(
            c3, &p->t0, p->A, H, p->Ns_max, p->prec);
    arb_add(res, c1, c3, p->prec);
    check_arb_accuracy(res, "interpolation error");

    arb_clear(c1);
    arb_clear(c3);
    arb_clear(H);

    return 1;
}

/*
 * Returns nonzero on success.
 */
static int
_get_K_truncation_error(arb_t res, const agg_t p)
{
    arb_t h, c, xi;

    arb_init(h);
    arb_init(c);
    arb_init(xi);

    _arb_div_si_si(h, p->hnum, p->hden, p->prec);
    _arb_inv_si(xi, p->B, p->prec);
    arb_mul_2exp_si(xi, xi, -1);
    arb_sqrt_ui(c, (ulong) p->J, p->prec);
    arb_mul_2exp_si(c, c, 1);
    arb_sub_ui(c, c, 1, p->prec);
    acb_dirichlet_platt_lemma_B2(res, p->K, h, xi, p->prec);
    arb_mul(res, res, c, p->prec);
    check_arb_accuracy(res, "K truncation error");

    arb_clear(h);
    arb_clear(c);
    arb_clear(xi);

    return 1;
}

/*
 * Returns nonzero on success.
 */
static int
_get_J_truncation_error(arb_t res, const agg_t p)
{
    arb_t h;
    arb_init(h);
    _arb_div_si_si(h, p->hnum, p->hden, p->prec);
    acb_dirichlet_platt_lemma_B1(res, p->sigma_grid, &p->t0, h, p->J, p->prec);
    check_arb_accuracy(res, "J truncation error");
    arb_clear(h);
    return 1;
}

/* a very rough error estimate of grid evaluation */
static int
_get_sum_of_errors(arb_t res, const agg_t p)
{
    arb_t x, err, ratio, c, xi, h;
    arb_t min_err, max_err, total_err;
    const int v = 0;

    arb_init(x);
    arb_init(err);
    arb_init(ratio);
    arb_init(c);
    arb_init(xi);
    arb_init(h);
    arb_init(min_err);
    arb_init(max_err);
    arb_init(total_err);

    _arb_div_si_si(h, p->hnum, p->hden, p->prec);
    _arb_inv_si(xi, p->B, p->prec);
    arb_mul_2exp_si(xi, xi, -1);

    if (v) flint_printf("(2) Adding Lemma A.5 truncation error... ");
    acb_dirichlet_platt_lemma_A5(err, p->B, h, 0, p->prec);
    arb_add(total_err, total_err, err, p->prec);

    if (v) flint_printf("(4) Adding Lemma A.7 truncation error... ");
    acb_dirichlet_platt_lemma_A7(
            err, p->sigma_grid, &p->t0, h, 0, p->A, p->prec);
    arb_add(total_err, total_err, err, p->prec);

    if (v) flint_printf("Adding Lemma B.1 truncation error for finite J... ");
    acb_dirichlet_platt_lemma_B1(err, p->sigma_grid, &p->t0, h, p->J, p->prec);
    arb_add(total_err, total_err, err, p->prec);
    if (v) {arb_printd(err, 20); flint_printf("\n");}

    if (v) flint_printf("Adding Lemma B.2 truncation error for finite K... ");
    arb_sqrt_ui(c, (ulong) p->J, p->prec);
    arb_mul_2exp_si(c, c, 1);
    arb_sub_ui(c, c, 1, p->prec);
    acb_dirichlet_platt_lemma_B2(err, p->K, h, xi, p->prec);
    arb_mul(err, err, c, p->prec);
    arb_add(total_err, total_err, err, p->prec);
    if (v) {arb_printd(err, 20); flint_printf("\n");}

    if (v) flint_printf("(6) Adding Lemma A.9 truncation error... ");
    acb_dirichlet_platt_lemma_A9(err, p->sigma_grid, &p->t0, h, p->A, p->prec);
    arb_add(total_err, total_err, err, p->prec);
    if (v) {arb_printd(err, 20); flint_printf("\n");}

    if (v) flint_printf("(8) Adding Lemma A.11 truncation error... ");
    acb_dirichlet_platt_lemma_A11(err, &p->t0, h, p->B, p->prec);
    arb_add(total_err, total_err, err, p->prec);
    if (v) {arb_printd(err, 20); flint_printf("\n");}

    arb_set(res, total_err);
    check_arb_accuracy(res, "sum of grid errors");

    arb_clear(x);
    arb_clear(err);
    arb_clear(ratio);
    arb_clear(c);
    arb_clear(xi);
    arb_clear(h);
    arb_clear(min_err);
    arb_clear(max_err);
    arb_clear(total_err);
    return 1;
}


/* An interpolation parameter. */
static slong
_update_interpolation_sigma(const agg_t agg)
{
    slong i, best_sigma;
    arb_t c1, H;
    mag_t m, best_m;
    agg_t p;
    int isfirst = 1;
    int v = 0;

    arb_init(c1);
    arb_init(H);
    mag_init(m);
    mag_init(best_m);
    agg_init(p);

    agg_set(p, agg);
    _arb_div_si_si(H, p->Hnum, p->Hden, p->prec);

    if (agg->sigma_interp % 2 == 0 || agg->sigma_interp < 3)
    {
        flint_printf("invalid initial sigma for interpolation: %wd\n",
                agg->sigma_interp);
        flint_abort();
    }

    for (i = -10; i <= 10; i += 2)
    {
        if (agg->sigma_interp + i < 3)
            continue;
        if (v) flint_printf("p->sigma_interp: %wd\n", p->sigma_interp);
        p->sigma_interp = agg->sigma_interp + i;
        if (p->sigma_interp % 2 == 0 || p->sigma_interp < 3)
        {
            flint_printf("invalid sigma for interpolation: %wd\n",
                    p->sigma_interp);
            flint_abort();
        }
        acb_dirichlet_platt_i_bound(c1,
                &p->t0, p->A, H, p->sigma_interp, p->prec);
        check_arb_accuracy(c1, "i bound in interpolation sigma update");
        arb_get_mag(m, c1);
        if (isfirst || mag_cmp(m, best_m) < 0)
        {
            mag_set(best_m, m);
            best_sigma = p->sigma_interp;
        }
        isfirst = 0;
    }

    arb_clear(c1);
    arb_clear(H);
    mag_clear(m);
    mag_clear(best_m);
    agg_clear(p);

    return best_sigma;
}

/* An interpolation parameter. */
static slong
_update_Hnum(const agg_t agg)
{
    slong i, best_Hnum;
    arb_t cost;
    mag_t m, best_m;
    agg_t p;
    int isfirst = 1;

    arb_init(cost);
    mag_init(m);
    mag_init(best_m);
    agg_init(p);

    agg_set(p, agg);

    for (i = -agg->Hden; i <= agg->Hden; i++)
    {
        if (agg->Hnum + i < 1)
            continue;
        p->Hnum = agg->Hnum + i;
        _get_interpolation_error(cost, p);
        arb_get_mag(m, cost);
        if (isfirst || mag_cmp(m, best_m) < 0)
        {
            mag_set(best_m, m);
            best_Hnum = p->Hnum;
        }
        isfirst = 0;
    }

    arb_clear(cost);
    mag_clear(m);
    mag_clear(best_m);
    agg_clear(p);

    return best_Hnum;
}

/*
 * Find small K whose associated truncation error is negligible
 * compared to the estimated interpolation error.
 */
static slong
_binary_search_K(const agg_t agg, const arb_t interpolation_error,
        slong low, slong high)
{
    int v = 0;
    if (high - low < 2)
    {
        return high;
    }
    else
    {
        arb_t err, midprop;
        agg_t p;
        slong mid;

        arb_init(err);
        arb_init(midprop);
        agg_init(p);

        agg_set(p, agg);

        mid = (high + low) >> 1;
        p->K = mid;
        if (!_get_K_truncation_error(midprop, p))
        {
            flint_printf("failed to compute a bound\n");
            flint_abort();
        }
        arb_div(midprop, midprop, interpolation_error, p->prec);
        if (v)
        {
            flint_printf("error proportion for K=%wd: ", mid);
            arb_printd(midprop, 20); flint_printf("\n");
        }
        arb_mul_2exp_si(midprop, midprop, 3);
        arb_log(midprop, midprop, p->prec);
        if (arb_is_negative(midprop))
        {
            return _binary_search_K(p, interpolation_error, low, mid);
        }
        else
        {
            return _binary_search_K(p, interpolation_error, mid, high);
        }

        arb_clear(err);
        arb_clear(midprop);
        agg_clear(p);
    }
}

/*
 * Find small J whose associated truncation error is negligible
 * compared to the estimated interpolation error.
 */
static slong
_binary_search_J(const agg_t agg, const arb_t interpolation_error,
        slong low, slong high)
{
    int v = 0;
    if (high - low < 2)
    {
        return high;
    }
    else
    {
        arb_t err, midprop;
        agg_t p;
        slong mid;

        arb_init(err);
        arb_init(midprop);
        agg_init(p);

        agg_set(p, agg);

        mid = (high + low) >> 1;
        p->J = mid;
        if (!_get_J_truncation_error(midprop, p))
        {
            flint_printf("failed to compute a bound\n");
            flint_abort();
        }
        arb_div(midprop, midprop, interpolation_error, p->prec);
        if (v)
        {
            flint_printf("error proportion for J=%wd: ", mid);
            arb_printd(midprop, 20); flint_printf("\n");
        }
        arb_mul_2exp_si(midprop, midprop, 3);
        arb_log(midprop, midprop, p->prec);
        if (arb_is_negative(midprop))
        {
            return _binary_search_J(p, interpolation_error, low, mid);
        }
        else
        {
            return _binary_search_J(p, interpolation_error, mid, high);
        }

        arb_clear(err);
        arb_clear(midprop);
        agg_clear(p);
    }
}

static slong
_update_hnum(const agg_t agg)
{
    arb_t err;
    arb_t best_err;
    slong best_hnum;
    slong i, c;
    int isfirst = 1;
    int v = 0;
    agg_t p;

    agg_init(p);
    arb_init(err);
    arb_init(best_err);

    agg_set(p, agg);

    c = 1;
    for (i = -c*agg->hden; i <= c*agg->hden; i++)
    {
        p->hnum = agg->hnum + i;
        if (!_get_sum_of_errors(err, p))
        {
            flint_printf("failed to compute a bound\n");
            flint_abort();
        }
        if (isfirst || arb_lt(err, best_err))
        {
            arb_set(best_err, err);
            best_hnum = p->hnum;
        }
        if (v)
        {
            flint_printf("%wd / %wd\t", p->hnum, p->hden);
            arb_printd(err, 20); flint_printf("\n");
        }
        isfirst = 0;
    }

    agg_clear(p);
    arb_clear(err);
    arb_clear(best_err);

    return best_hnum;
}

static slong
_update_grid_sigma(const agg_t agg)
{
    agg_t p;
    arb_t err, best_err;
    slong i, best_sigma;
    int isfirst = 1;
    int v = 0;
    slong offset[] = {
        -80, -40, -20, -10, -8, -6, -4, -2, 0,
        2, 4, 6, 8, 10, 20, 40, 80};

    agg_init(p);
    arb_init(err);
    arb_init(best_err);

    agg_set(p, agg);

    if (agg->sigma_grid % 2 == 0 || agg->sigma_grid < 3)
    {
        flint_printf("invalid initial sigma for grid: %wd\n", agg->sigma_grid);
        flint_abort();
    }

    for (i = 0; i < 17; i++)
    {
        if (agg->sigma_grid + offset[i] < 3)
            continue;
        p->sigma_grid = agg->sigma_grid + offset[i];
        if (p->sigma_grid % 2 == 0 || p->sigma_grid < 3)
        {
            flint_printf("invalid sigma for grid: %wd\n", p->sigma_grid);
            flint_abort();
        }
        if (!_get_sum_of_errors(err, p))
        {
            flint_printf("failed to compute a bound\n");
            flint_abort();
        }
        if (isfirst || arb_lt(err, best_err))
        {
            best_sigma = p->sigma_grid;
            arb_set(best_err, err);
        }
        if (v)
        {
            flint_printf("%wd\t", p->sigma_grid);
            arb_printd(err, 20); flint_printf("\n");
        }
        isfirst = 0;
    }

    agg_clear(p);
    arb_clear(err);
    arb_clear(best_err);

    return best_sigma;
}

/* all params are output except for n which indicates nth zero */
static void
_get_params(fmpz_t T, slong *A, slong *B, slong *prec,
        slong *sigma_grid, slong *K, slong *J, slong *Jmax,
        slong *hnum, slong *hden,
        slong *sigma_interp, slong *Hnum, slong *Hden, slong *Ns_max,
        const fmpz_t n)
{
    fmpz_t k;
    arb_t g, rho, pi;
    acb_struct z[2];
    slong Abits, Bbits;

    fmpz_init(k);
    arb_init(g);
    arb_init(rho);
    arb_init(pi);
    acb_init(z+0);
    acb_init(z+1);

    /* Default prec. */
    *prec = 64;

    /* Default max number of interpolation points on either side. */
    *Ns_max = 300;

    /* Estimate the height of the nth zero using gram points --
     * it's predicted to fall between g(n-2) and g(n-1). */
    fmpz_sub_ui(k, n, 2);
    acb_dirichlet_gram_point(g, k, NULL, NULL, *prec);

    /*
     * Estimate the local density of Riemann zeta zeros near height g.
     * Currently unused, but maybe in the future this could
     * be used to set A and B.
     */
    arb_const_pi(pi, *prec);
    acb_set_arb(z, g);
    acb_dirichlet_hardy_theta(z, z, NULL, NULL, 2, *prec);
    arb_div(rho, acb_realref(z+1), pi, *prec);

    Abits = 4;
    Bbits = 13;
    *A = 1 << Abits;
    *B = 1 << Bbits;

    /*
     * Let T be the integer at the center of the evaluation grid.
     * The reason to add B/4 instead of B/2 is that only about the
     * middle half of the grid points have enough precision to isolate zeros.
     * If you want to increase the probability of finding the
     * first requested zero at the expense of the number of zeros
     * found, then don't add anything to T.
     */
    _arb_get_lbound_fmpz(T, g, *prec);
    /* fmpz_add_ui(T, T, (*B)/4); */

    /* According to Platt, J should be like the square root of the height.
     * We use a + b*sqrt(T) where a and b are chosen to yield empirically
     * acceptable values of J at a small and a large value of T. */
    fmpz_sqrt(k, T);
    *Jmax = 1000 + 4*fmpz_get_ui(k);

    *J = 19957244503;
    *K = 63;

    *hnum = 702;
    *hden = 4;

    *sigma_grid = 3901;

    *Hnum = 52;
    *Hden = 64;
    *sigma_interp = 19;

    fmpz_clear(k);
    arb_clear(g);
    arb_clear(rho);
    arb_clear(pi);
    acb_clear(z+0);
    acb_clear(z+1);
}

void run(const fmpz_t n, slong k, slong m)
{
    slong Jmax;
    agg_t p, prev;
    arb_t err, interpolation_error;
    const slong itermax = 10;
    slong iter;

    arb_init(err);
    arb_init(interpolation_error);
    agg_init(p);
    agg_init(prev);

    _get_params(
            &p->T, &p->A, &p->B, &p->prec, &p->sigma_grid,
            &p->K, &p->J, &Jmax, &p->hnum, &p->hden,
            &p->sigma_interp, &p->Hnum, &p->Hden, &p->Ns_max, n);
    arb_set_fmpz(&p->t0, &p->T);

    flint_printf("parameter values:\n");
    agg_printf(p);
    flint_printf("\n");

    flint_printf("updating interpolation related parameters...\n");
    for (iter = 1; iter <= itermax; iter++)
    {
        flint_printf("updating H (iter=%wd/%wd)...\n", iter, itermax);
        p->Hnum = _update_Hnum(p);
        _get_interpolation_error(interpolation_error, p);
        flint_printf("\n");
        flint_printf("H = %wd / %wd\n", p->Hnum, p->Hden);
        flint_printf("interpolation error estimate = ");
        arb_printd(interpolation_error, 20); flint_printf("\n");
        flint_printf("n = %wd*10^%wd\n", k, m);

        flint_printf("updating interpolation sigma (iter=%wd/%wd)...\n",
                iter, itermax);
        p->sigma_interp = _update_interpolation_sigma(p);
        _get_interpolation_error(interpolation_error, p);
        flint_printf("\n");
        flint_printf("sigma_interp = %wd\n", p->sigma_interp);
        flint_printf("interpolation error estimate = ");
        arb_printd(interpolation_error, 20); flint_printf("\n");
        flint_printf("n = %wd*10^%wd\n", k, m);

        if (iter > 1 &&
            p->Hnum == prev->Hnum &&
            p->sigma_interp == prev->sigma_interp)
        {
            break;
        }

        agg_set(prev, p);
    }
    flint_printf("\n");

    /* Requires the interpolation error. */
    flint_printf("updating J (iter=%wd/%wd)...\n", 0, itermax);
    p->J = _binary_search_J(p, interpolation_error, 1, Jmax);
    _get_sum_of_errors(err, p);
    flint_printf("\n");
    flint_printf("J = %wd\n", p->J);
    flint_printf("err = "); arb_printd(err, 20); flint_printf("\n");
    flint_printf("n = %wd*10^%wd\n", k, m);

    /* Requires the interpolation error. */
    flint_printf("updating K (iter=%wd/%wd)...\n", 0, itermax);
    p->K = _binary_search_K(p, interpolation_error, 1, 500);
    _get_sum_of_errors(err, p);
    flint_printf("\n");
    flint_printf("K = %wd\n", p->K);
    flint_printf("err = "); arb_printd(err, 20); flint_printf("\n");
    flint_printf("n = %wd*10^%wd\n", k, m);

    flint_printf("updating multieval related parameters...\n");
    for (iter = 1; iter <= itermax; iter++)
    {
        flint_printf("updating h (iter=%wd/%wd)...\n", iter, itermax);
        p->hnum = _update_hnum(p);
        _get_sum_of_errors(err, p);
        flint_printf("\n");
        flint_printf("h = %wd / %wd\n", p->hnum, p->hden);
        flint_printf("err = "); arb_printd(err, 20); flint_printf("\n");
        flint_printf("n = %wd*10^%wd\n", k, m);

        flint_printf("updating grid sigma (iter=%wd/%wd)...\n",
                iter, itermax);
        p->sigma_grid = _update_grid_sigma(p);
        _get_sum_of_errors(err, p);
        flint_printf("\n");
        flint_printf("sigma_grid = %wd\n", p->sigma_grid);
        flint_printf("err = "); arb_printd(err, 20); flint_printf("\n");
        flint_printf("n = %wd*10^%wd\n", k, m);

        /* Requires the interpolation error. */
        flint_printf("updating J (iter=%wd/%wd)...\n", iter, itermax);
        p->J = _binary_search_J(p, interpolation_error, 1, Jmax);
        _get_sum_of_errors(err, p);
        flint_printf("\n");
        flint_printf("J = %wd\n", p->J);
        flint_printf("err = "); arb_printd(err, 20); flint_printf("\n");
        flint_printf("n = %wd*10^%wd\n", k, m);

        /* Requires the interpolation error. */
        flint_printf("updating K (iter=%wd/%wd)...\n", iter, itermax);
        p->K = _binary_search_K(p, interpolation_error, 1, 500);
        _get_sum_of_errors(err, p);
        flint_printf("\n");
        flint_printf("K = %wd\n", p->K);
        flint_printf("err = "); arb_printd(err, 20); flint_printf("\n");
        flint_printf("n = %wd*10^%wd\n", k, m);

        if (iter > 1 &&
            p->J == prev->J &&
            p->K == prev->K &&
            p->hnum == prev->hnum &&
            p->sigma_grid == prev->sigma_grid)
        {
            break;
        }

        agg_set(prev, p);
    }
    flint_printf("\n");

    flint_printf("*** summary ***\n");
    agg_printf(p);

    _get_sum_of_errors(err, p);
    flint_printf("grid error estimate = ");
    arb_printd(err, 20); flint_printf("\n");

    _get_interpolation_error(err, p);
    flint_printf("interpolation error estimate = ");
    arb_printd(err, 20); flint_printf("\n");

    flint_printf("\n");

finish:

    arb_clear(err);
    arb_clear(interpolation_error);
    agg_clear(p);
    agg_clear(prev);
}

int main()
{
    fmpz_t n;
    slong m, k;

    /*
     * Disable the stdout buffer so that no message is lost
     * if the program crashes or is interrupted
     */
    setbuf(stdout, NULL);

    fmpz_init(n);

    /* n = k * 10^m */
    {
        k = 1;
        m = 22;
        flint_printf("trying n = %wd*10^%wd\n", k, m);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, m);
        fmpz_mul_ui(n, n, k);
        run(n, k, m);
    }

    fmpz_clear(n);

    flint_cleanup();
    return 0;
}
