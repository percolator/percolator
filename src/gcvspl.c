#include "gcvspl.h"


/* gcvspl.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/


#ifdef __cplusplus
extern "C" {
#endif

/* Table of constant values */

static doublereal c_b44 = 1e-15;

/* GCVSPL.FOR, 1986-05-12 */

/* *********************************************************************** */

/* SUBROUTINE GCVSPL (REAL*8) */

/* Purpose: */
/* ******* */

/*       Natural B-spline data smoothing subroutine, using the Generali- */
/*       zed Cross-Validation and Mean-Squared Prediction Error Criteria */
/*       of Craven & Wahba (1979). Alternatively, the amount of smoothing */
/*       can be given explicitly, or it can be based on the effective */
/*       number of degrees of freedom in the smoothing process as defined */
/*       by Wahba (1980). The model assumes uncorrelated, additive noise */
/*       and essentially smooth, underlying functions. The noise may be */
/*       non-stationary, and the independent co-ordinates may be spaced */
/*       non-equidistantly. Multiple datasets, with common independent */
/*       variables and weight factors are accomodated. */


/* Calling convention: */
/* ****************** */

/*       CALL GCVSPL ( X, Y, NY, WX, WY, M, N, K, MD, VAL, C, NC, WK, IER ) */

/* Meaning of parameters: */
/* ********************* */

/*       X(N)    ( I )   Independent variables: strictly increasing knot */
/*                       sequence, with X(I-1).lt.X(I), I=2,...,N. */
/*       Y(NY,K) ( I )   Input data to be smoothed (or interpolated). */
/*       NY      ( I )   First dimension of array Y(NY,K), with NY.ge.N. */
/*       WX(N)   ( I )   Weight factor array; WX(I) corresponds with */
/*                       the relative inverse variance of point Y(I,*). */
/*                       If no relative weighting information is */
/*                       available, the WX(I) should be set to ONE. */
/*                       All WX(I).gt.ZERO, I=1,...,N. */
/*       WY(K)   ( I )   Weight factor array; WY(J) corresponds with */
/*                       the relative inverse variance of point Y(*,J). */
/*                       If no relative weighting information is */
/*                       available, the WY(J) should be set to ONE. */
/*                       All WY(J).gt.ZERO, J=1,...,K. */
/*                       NB: The effective weight for point Y(I,J) is */
/*                       equal to WX(I)*WY(J). */
/*       M       ( I )   Half order of the required B-splines (spline */
/*                       degree 2*M-1), with M.gt.0. The values M = */
/*                       1,2,3,4 correspond to linear, cubic, quintic, */
/*                       and heptic splines, respectively. */
/*       N       ( I )   Number of observations per dataset, with N.ge.2*M. */
/*       K       ( I )   Number of datasets, with K.ge.1. */
/*       MD      ( I )   Optimization mode switch: */
/*                       |MD| = 1: Prior given value for p in VAL */
/*                                 (VAL.ge.ZERO). This is the fastest */
/*                                 use of GCVSPL, since no iteration */
/*                                 is performed in p. */
/*                       |MD| = 2: Generalized cross validation. */
/*                       |MD| = 3: True predicted mean-squared error, */
/*                                 with prior given variance in VAL. */
/*                       |MD| = 4: Prior given number of degrees of */
/*                                 freedom in VAL (ZERO.le.VAL.le.N-M). */
/*                        MD  < 0: It is assumed that the contents of */
/*                                 X, W, M, N, and WK have not been */
/*                                 modified since the previous invoca- */
/*                                 tion of GCVSPL. If MD < -1, WK(4) */
/*                                 is used as an initial estimate for */
/*                                 the smoothing parameter p. */
/*                       Other values for |MD|, and inappropriate values */
/*                       for VAL will result in an error condition, or */
/*                       cause a default value for VAL to be selected. */
/*                       After return from MD.ne.1, the same number of */
/*                       degrees of freedom can be obtained, for identical */
/*                       weight factors and knot positions, by selecting */
/*                       |MD|=1, and by copying the value of p from WK(4) */
/*                       into VAL. In this way, no iterative optimization */
/*                       is required when processing other data in Y. */
/*       VAL     ( I )   Mode value, as described above under MD. */
/*       C(NC,K) ( O )   Spline coefficients, to be used in conjunction */
/*                       with function SPLDER. NB: the dimensions of C */
/*                       in GCVSPL and in SPLDER are different! In SPLDER, */
/*                       only a single column of C(N,K) is needed, and the */
/*                       proper column C(1,J), with J=1...K should be used */
/*                       when calling SPLDER. */
/*       NC       ( I )  First dimension of array C(NC,K), NC.ge.N. */
/*       WK(IWK) (I/W/O) Work vector, with length IWK.ge.6*(N*M+1)+N. */
/*                       On normal exit, the first 6 values of WK are */
/*                       assigned as follows: */

/*                       WK(1) = Generalized Cross Validation value */
/*                       WK(2) = Mean Squared Residual. */
/*                       WK(3) = Estimate of the number of degrees of */
/*                               freedom of the residual sum of squares */
/*                               per dataset, with 0.lt.WK(3).lt.N-M. */
/*                       WK(4) = Smoothing parameter p, multiplicative */
/*                               with the splines' derivative constraint. */
/*                       WK(5) = Estimate of the true mean squared error */
/*                               (different formula for |MD| = 3). */
/*                       WK(6) = Gauss-Markov error variance. */

/*                       If WK(4) -->  0 , WK(3) -->  0 , and an inter- */
/*                       polating spline is fitted to the data (p --> 0). */
/*                       A very small value > 0 is used for p, in order */
/*                       to avoid division by zero in the GCV function. */

/*                       If WK(4) --> inf, WK(3) --> N-M, and a least- */
/*                       squares polynomial of order M (degree M-1) is */
/*                       fitted to the data (p --> inf). For numerical */
/*                       reasons, a very high value is used for p. */

/*                       Upon return, the contents of WK can be used for */
/*                       covariance propagation in terms of the matrices */
/*                       B and WE: see the source listings. The variance */
/*                       estimate for dataset J follows as WK(6)/WY(J). */

/*       IER     ( O )   Error parameter: */

/*                       IER = 0:        Normal exit */
/*                       IER = 1:        M.le.0 .or. N.lt.2*M */
/*                       IER = 2:        Knot sequence is not strictly */
/*                                       increasing, or some weight */
/*                                       factor is not positive. */
/*                       IER = 3:        Wrong mode  parameter or value. */

/* Remarks: */
/* ******* */

/*       (1) GCVSPL calculates a natural spline of order 2*M (degree */
/*       2*M-1) which smoothes or interpolates a given set of data */
/*       points, using statistical considerations to determine the */
/*       amount of smoothing required (Craven & Wahba, 1979). If the */
/*       error variance is a priori known, it should be supplied to */
/*       the routine in VAL, for |MD|=3. The degree of smoothing is */
/*       then determined to minimize an unbiased estimate of the true */
/*       mean squared error. On the other hand, if the error variance */
/*       is not known, one may select |MD|=2. The routine then deter- */
/*       mines the degree of smoothing to minimize the generalized */
/*       cross validation function. This is asymptotically the same */
/*       as minimizing the true predicted mean squared error (Craven & */
/*       Wahba, 1979). If the estimates from |MD|=2 or 3 do not appear */
/*       suitable to the user (as apparent from the smoothness of the */
/*       M-th derivative or from the effective number of degrees of */
/*       freedom returned in WK(3) ), the user may select an other */
/*       value for the noise variance if |MD|=3, or a reasonably large */
/*       number of degrees of freedom if |MD|=4. If |MD|=1, the proce- */
/*       dure is non-iterative, and returns a spline for the given */
/*       value of the smoothing parameter p as entered in VAL. */

/*       (2) The number of arithmetic operations and the amount of */
/*       storage required are both proportional to N, so very large */
/*       datasets may be accomodated. The data points do not have */
/*       to be equidistant in the independant variable X or uniformly */
/*       weighted in the dependant variable Y. However, the data */
/*       points in X must be strictly increasing. Multiple dataset */
/*       processing (K.gt.1) is numerically more efficient dan */
/*       separate processing of the individual datasets (K.eq.1). */

/*       (3) If |MD|=3 (a priori known noise variance), any value of */
/*       N.ge.2*M is acceptable. However, it is advisable for N-2*M */
/*       be rather large (at least 20) if |MD|=2 (GCV). */

/*       (4) For |MD| > 1, GCVSPL tries to iteratively minimize the */
/*       selected criterion function. This minimum is unique for |MD| */
/*       = 4, but not necessarily for |MD| = 2 or 3. Consequently, */
/*       local optima rather that the global optimum might be found, */
/*       and some actual findings suggest that local optima might */
/*       yield more meaningful results than the global optimum if N */
/*       is small. Therefore, the user has some control over the */
/*       search procedure. If MD > 1, the iterative search starts */
/*       from a value which yields a number of degrees of freedom */
/*       which is approximately equal to N/2, until the first (local) */
/*       minimum is found via a golden section search procedure */
/*       (Utreras, 1980). If MD < -1, the value for p contained in */
/*       WK(4) is used instead. Thus, if MD = 2 or 3 yield too noisy */
/*       an estimate, the user might try |MD| = 1 or 4, for suitably */
/*       selected values for p or for the number of degrees of */
/*       freedom, and then run GCVSPL with MD = -2 or -3. The con- */
/*       tents of N, M, K, X, WX, WY, and WK are assumed unchanged */
/*       if MD < 0. */

/*       (5) GCVSPL calculates the spline coefficient array C(N,K); */
/*       this array can be used to calculate the spline function */
/*       value and any of its derivatives up to the degree 2*M-1 */
/*       at any argument T within the knot range, using subrou- */
/*       tines SPLDER and SEARCH, and the knot array X(N). Since */
/*       the splines are constrained at their Mth derivative, only */
/*       the lower spline derivatives will tend to be reliable */
/*       estimates of the underlying, true signal derivatives. */

/*       (6) GCVSPL combines elements of subroutine CRVO5 by Utre- */
/*       ras (1980), subroutine SMOOTH by Lyche et al. (1983), and */
/*       subroutine CUBGCV by Hutchinson (1985). The trace of the */
/*       influence matrix is assessed in a similar way as described */
/*       by Hutchinson & de Hoog (1985). The major difference is */
/*       that the present approach utilizes non-symmetrical B-spline */
/*       design matrices as described by Lyche et al. (1983); there- */
/*       fore, the original algorithm by Erisman & Tinney (1975) has */
/*       been used, rather than the symmetrical version adopted by */
/*       Hutchinson & de Hoog. */

/* References: */
/* ********** */

/*       P. Craven & G. Wahba (1979), Smoothing noisy data with */
/*       spline functions. Numerische Mathematik 31, 377-403. */

/*       A.M. Erisman & W.F. Tinney (1975), On computing certain */
/*       elements of the inverse of a sparse matrix. Communications */
/*       of the ACM 18(3), 177-179. */

/*       M.F. Hutchinson & F.R. de Hoog (1985), Smoothing noisy data */
/*       with spline functions. Numerische Mathematik 47(1), 99-106. */

/*       M.F. Hutchinson (1985), Subroutine CUBGCV. CSIRO Division of */
/*       Mathematics and Statistics, P.O. Box 1965, Canberra, ACT 2601, */
/*       Australia. */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori (1983), Fortran */
/*       subroutines for computing smoothing and interpolating natural */
/*       splines. Advances in Engineering Software 5(1), 2-5. */

/*       F. Utreras (1980), Un paquete de programas para ajustar curvas */
/*       mediante funciones spline. Informe Tecnico MA-80-B-209, Depar- */
/*       tamento de Matematicas, Faculdad de Ciencias Fisicas y Matema- */
/*       ticas, Universidad de Chile, Santiago. */

/*       Wahba, G. (1980). Numerical and statistical methods for mildly, */
/*       moderately and severely ill-posed problems with noisy data. */
/*       Technical report nr. 595 (February 1980). Department of Statis- */
/*       tics, University of Madison (WI), U.S.A. */

/* Subprograms required: */
/* ******************** */

/*       BASIS, PREP, SPLC, BANDET, BANSOL, TRINV */

/* *********************************************************************** */

/* Subroutine */ int gcvspl_(doublereal *x, doublereal *y, integer *ny, 
	doublereal *wx, doublereal *wy, integer *m, integer *n, integer *k, 
	integer *md, doublereal *val, doublereal *c__, integer *nc, 
	doublereal *wk, integer *ier)
{
    /* Initialized data */

    static integer m2 = 0;
    static integer nm1 = 0;
    static doublereal el = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal r1, r2, r3, r4;
    static integer ib;
    static doublereal gf1, gf2, gf3, gf4;
    static integer iwe;
    static doublereal err;
    static integer nm2m1, nm2p1;
    extern doublereal splc_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int prep_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal alpha;
    extern /* Subroutine */ int basis_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


    /* Parameter adjustments */
    --wk;
    --wx;
    --x;
    --wy;
    y_dim1 = *ny;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    c_dim1 = *nc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

/* ***  Parameter check and work array initialization */

    *ier = 0;
/* ***  Check on mode parameter */
    if (abs(*md) > 4 || *md == 0 || (abs(*md) == 1 && *val < 0.) || (abs(*md) == 
	    3 && *val < 0.) || (abs(*md) == 4 && (*val < 0. || *val > (
	    doublereal) (*n - *m)))) {
	*ier = 3;
/* Wrong mode value */
	return 0;
    }
/* ***  Check on M and N */
    if (*md > 0) {
	m2 = *m << 1;
	nm1 = *n - 1;
    } else {
	if (m2 != *m << 1 || nm1 != *n - 1) {
	    *ier = 3;
/* M or N modified since previous call */
	    return 0;
	}
    }
    if (*m <= 0 || *n < m2) {
	*ier = 1;
/* M or N invalid */
	return 0;
    }
/* ***  Check on knot sequence and weights */
    if (wx[1] <= 0.) {
	*ier = 2;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (wx[i__] <= 0. || x[i__ - 1] >= x[i__]) {
	    *ier = 2;
	}
	if (*ier != 0) {
	    return 0;
	}
/* L10: */
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	if (wy[j] <= 0.) {
	    *ier = 2;
	}
	if (*ier != 0) {
	    return 0;
	}
/* L15: */
    }

/* ***  Work array parameters (address information for covariance */
/* ***  propagation by means of the matrices STAT, B, and WE). NB: */
/* ***  BWE cannot be used since it is modified by function TRINV. */

    nm2p1 = *n * (m2 + 1);
    nm2m1 = *n * (m2 - 1);
/*     ISTAT = 1            !Statistics array STAT(6) */
/*     IBWE  = ISTAT + 6      !Smoothing matrix BWE( -M:M  ,N) */
    ib = nm2p1 + 7;
/* Design matrix    B  (1-M:M-1,N) */
    iwe = ib + nm2m1;
/*     IWK   = IWE   + NM2P1      !Total work array length N + 6*(N*M+1) */

/* ***  Compute the design matrices B and WE, the ratio */
/* ***  of their L1-norms, and check for iterative mode. */

/* Design matrix    WE ( -M:M  ,N) */
    if (*md > 0) {
	basis_(m, n, &x[1], &wk[ib], &r1, &wk[7]);
	prep_(m, n, &x[1], &wx[1], &wk[iwe], &el);
	el /= r1;
/* L1-norms ratio (SAVEd upon RETURN) */
    }
    if (abs(*md) != 1) {
	goto L20;
    }
/* ***     Prior given value for p */
    r1 = *val;
    goto L100;

/* ***  Iterate to minimize the GCV function (|MD|=2), */
/* ***  the MSE function (|MD|=3), or to obtain the prior */
/* ***  given number of degrees of freedom (|MD|=4). */

L20:
    if (*md < -1) {
	r1 = wk[4];
/* User-determined starting value */
    } else {
	r1 = 1. / el;
/* Default (DOF ~ 0.5) */
    }
    r2 = r1 * 2.;
    gf2 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r2, &
	    c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7])
	    ;
L40:
    gf1 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r1, &
	    c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7])
	    ;
    if (gf1 > gf2) {
	goto L50;
    }
    if (wk[4] <= 0.) {
	goto L100;
    }
/* Interpolation */
    r2 = r1;
    gf2 = gf1;
    r1 /= 2.;
    goto L40;
L50:
    r3 = r2 * 2.;
L60:
    gf3 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r3, &
	    c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7])
	    ;
    if (gf3 > gf2) {
	goto L70;
    }
    if (wk[4] >= 999999999999999.88) {
	goto L100;
    }
/* Least-squares polynomial */
    r2 = r3;
    gf2 = gf3;
    r3 *= 2.;
    goto L60;
L70:
    r2 = r3;
    gf2 = gf3;
    alpha = (r2 - r1) / 1.618033983;
    r4 = r1 + alpha;
    r3 = r2 - alpha;
    gf3 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r3, &
	    c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7])
	    ;
    gf4 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r4, &
	    c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7])
	    ;
L80:
    if (gf3 <= gf4) {
	r2 = r4;
	gf2 = gf4;
	err = (r2 - r1) / (r1 + r2);
	if (err * err + 1. == 1. || err <= 1e-6) {
	    goto L90;
	}
	r4 = r3;
	gf4 = gf3;
	alpha /= 1.618033983;
	r3 = r2 - alpha;
	gf3 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r3, &
		c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &
		wk[7]);
    } else {
	r1 = r3;
	gf1 = gf3;
	err = (r2 - r1) / (r1 + r2);
	if (err * err + 1. == 1. || err <= 1e-6) {
	    goto L90;
	}
	r3 = r4;
	gf3 = gf4;
	alpha /= 1.618033983;
	r4 = r1 + alpha;
	gf4 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r4, &
		c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &
		wk[7]);
    }
    goto L80;
L90:
    r1 = (r1 + r2) * .5;

/* ***  Calculate final spline coefficients */

L100:
    gf1 = splc_(m, n, k, &y[y_offset], ny, &wx[1], &wy[1], md, val, &r1, &
	    c_b44, &c__[c_offset], nc, &wk[1], &wk[ib], &wk[iwe], &el, &wk[7])
	    ;

/* ***  Ready */

    return 0;
} /* gcvspl_ */

/* BASIS.FOR, 1985-06-03 */

/* *********************************************************************** */

/* SUBROUTINE BASIS (REAL*8) */

/* Purpose: */
/* ******* */

/*       Subroutine to assess a B-spline tableau, stored in vectorized */
/*       form. */

/* Calling convention: */
/* ****************** */

/*       CALL BASIS ( M, N, X, B, BL, Q ) */

/* Meaning of parameters: */
/* ********************* */

/*       M               ( I )   Half order of the spline (degree 2*M-1), */
/*                               M > 0. */
/*       N               ( I )   Number of knots, N >= 2*M. */
/*       X(N)            ( I )   Knot sequence, X(I-1) < X(I), I=2,N. */
/*       B(1-M:M-1,N)    ( O )   Output tableau. Element B(J,I) of array */
/*                               B corresponds with element b(i,i+j) of */
/*                               the tableau matrix B. */
/*       BL              ( O )   L1-norm of B. */
/*       Q(1-M:M)        ( W )   Internal work array. */

/* Remark: */
/* ****** */

/*       This subroutine is an adaptation of subroutine BASIS from the */
/*       paper by Lyche et al. (1983). No checking is performed on the */
/*       validity of M and N. If the knot sequence is not strictly in- */
/*       creasing, division by zero may occur. */

/* Reference: */
/* ********* */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines */
/*       for computing smoothing and interpolating natural splines. */
/*       Advances in Engineering Software 5(1983)1, pp. 2-5. */

/* *********************************************************************** */

/* Subroutine */ int basis_(integer *m, integer *n, doublereal *x, doublereal 
	*b, doublereal *bl, doublereal *q)
{
    /* System generated locals */
    integer b_dim1, b_offset, q_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal u, v, y;
    static integer j1, j2, m2, ir, mm1, mp1;
    static doublereal arg;
    static integer nmip1;



    /* Parameter adjustments */
    q_offset = 1 - *m;
    q -= q_offset;
    b_dim1 = *m - 1 - (1 - *m) + 1;
    b_offset = 1 - *m + b_dim1;
    b -= b_offset;
    --x;

    /* Function Body */
    if (*m == 1) {
/* ***         Linear spline */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    b[i__ * b_dim1] = 1.;
/* L3: */
	}
	*bl = 1.;
	return 0;
    }

/* ***  General splines */

    mm1 = *m - 1;
    mp1 = *m + 1;
    m2 = *m << 1;
    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
/* ***     1st row */
	i__2 = *m;
	for (j = -mm1; j <= i__2; ++j) {
	    q[j] = 0.;
/* L5: */
	}
	q[mm1] = 1.;
	if (l != 1 && l != *n) {
	    q[mm1] = 1. / (x[l + 1] - x[l - 1]);
	}
/* ***     Successive rows */
	arg = x[l];
	i__2 = m2;
	for (i__ = 3; i__ <= i__2; ++i__) {
	    ir = mp1 - i__;
	    v = q[ir];
	    if (l < i__) {
/* ***               Left-hand B-splines */
		i__3 = i__;
		for (j = l + 1; j <= i__3; ++j) {
		    u = v;
		    v = q[ir + 1];
		    q[ir] = u + (x[j] - arg) * v;
		    ++ir;
/* L6: */
		}
	    }
/* Computing MAX */
	    i__3 = l - i__ + 1;
	    j1 = max(i__3,1);
/* Computing MIN */
	    i__3 = l - 1, i__4 = *n - i__;
	    j2 = min(i__3,i__4);
	    if (j1 <= j2) {
/* ***               Ordinary B-splines */
		if (i__ < m2) {
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			y = x[i__ + j];
			u = v;
			v = q[ir + 1];
			q[ir] = u + (v - u) * (y - arg) / (y - x[j]);
			++ir;
/* L8: */
		    }
		} else {
		    i__3 = j2;
		    for (j = j1; j <= i__3; ++j) {
			u = v;
			v = q[ir + 1];
			q[ir] = (arg - x[j]) * u + (x[i__ + j] - arg) * v;
			++ir;
/* L10: */
		    }
		}
	    }
	    nmip1 = *n - i__ + 1;
	    if (nmip1 < l) {
/* ***           Right-hand B-splines */
		i__3 = l - 1;
		for (j = nmip1; j <= i__3; ++j) {
		    u = v;
		    v = q[ir + 1];
		    q[ir] = (arg - x[j]) * u + v;
		    ++ir;
/* L12: */
		}
	    }
/* L13: */
	}
	i__2 = mm1;
	for (j = -mm1; j <= i__2; ++j) {
	    b[j + l * b_dim1] = q[j];
/* L14: */
	}
/* L15: */
    }

/* ***  Zero unused parts of B */

    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mm1;
	for (k = i__; k <= i__2; ++k) {
	    b[-k + i__ * b_dim1] = 0.;
	    b[k + (*n + 1 - i__) * b_dim1] = 0.;
/* L16: */
	}
/* L17: */
    }

/* ***  Assess L1-norm of B */

    *bl = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mm1;
	for (k = -mm1; k <= i__2; ++k) {
	    *bl += (d__1 = b[k + i__ * b_dim1], abs(d__1));
/* L18: */
	}
/* L19: */
    }
    *bl /= *n;

/* ***  Ready */

    return 0;
} /* basis_ */

/* PREP.FOR, 1985-07-04 */

/* *********************************************************************** */

/* SUBROUTINE PREP (REAL*8) */

/* Purpose: */
/* ******* */

/*       To compute the matrix WE of weighted divided difference coeffi- */
/*       cients needed to set up a linear system of equations for sol- */
/*       ving B-spline smoothing problems, and its L1-norm EL. The matrix */
/*       WE is stored in vectorized form. */

/* Calling convention: */
/* ****************** */

/*       CALL PREP ( M, N, X, W, WE, EL ) */

/* Meaning of parameters: */
/* ********************* */

/*       M               ( I )   Half order of the B-spline (degree */
/*                               2*M-1), with M > 0. */
/*       N               ( I )   Number of knots, with N >= 2*M. */
/*       X(N)            ( I )   Strictly increasing knot array, with */
/*                               X(I-1) < X(I), I=2,N. */
/*       W(N)            ( I )   Weight matrix (diagonal), with */
/*                               W(I).gt.0.0, I=1,N. */
/*       WE(-M:M,N)      ( O )   Array containing the weighted divided */
/*                               difference terms in vectorized format. */
/*                               Element WE(J,I) of array E corresponds */
/*                               with element e(i,i+j) of the matrix */
/*                               W**-1 * E. */
/*       EL              ( O )   L1-norm of WE. */

/* Remark: */
/* ****** */

/*       This subroutine is an adaptation of subroutine PREP from the paper */
/*       by Lyche et al. (1983). No checking is performed on the validity */
/*       of M and N. Division by zero may occur if the knot sequence is */
/*       not strictly increasing. */

/* Reference: */
/* ********* */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines */
/*       for computing smoothing and interpolating natural splines. */
/*       Advances in Engineering Software 5(1983)1, pp. 2-5. */

/* *********************************************************************** */

/* Subroutine */ int prep_(integer *m, integer *n, doublereal *x, doublereal *
	w, doublereal *we, doublereal *el)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal f;
    static integer i__, j, k, l;
    static doublereal y, f1;
    static integer i1, i2, m2;
    static doublereal ff;
    static integer jj, jm, kl, nm, ku;
    static doublereal wi;
    static integer n2m, mp1, i2m1, inc, i1p1, m2m1, m2p1;



/* ***  Calculate the factor F1 */

/* WE(-M:M,N) */
    /* Parameter adjustments */
    --we;
    --w;
    --x;

    /* Function Body */
    m2 = *m << 1;
    mp1 = *m + 1;
    m2m1 = m2 - 1;
    m2p1 = m2 + 1;
    nm = *n - *m;
    f1 = -1.;
    if (*m != 1) {
	i__1 = *m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    f1 = -f1 * i__;
/* L5: */
	}
	i__1 = m2m1;
	for (i__ = mp1; i__ <= i__1; ++i__) {
	    f1 *= i__;
/* L6: */
	}
    }

/* ***  Columnwise evaluation of the unweighted design matrix E */

    i1 = 1;
    i2 = *m;
    jm = mp1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	inc = m2p1;
	if (j > nm) {
	    f1 = -f1;
	    f = f1;
	} else {
	    if (j < mp1) {
		inc = 1;
		f = f1;
	    } else {
		f = f1 * (x[j + *m] - x[j - *m]);
	    }
	}
	if (j > mp1) {
	    ++i1;
	}
	if (i2 < *n) {
	    ++i2;
	}
	jj = jm;
/* ***     Loop for divided difference coefficients */
	ff = f;
	y = x[i1];
	i1p1 = i1 + 1;
	i__2 = i2;
	for (i__ = i1p1; i__ <= i__2; ++i__) {
	    ff /= y - x[i__];
/* L11: */
	}
	we[jj] = ff;
	jj += m2;
	i2m1 = i2 - 1;
	if (i1p1 <= i2m1) {
	    i__2 = i2m1;
	    for (l = i1p1; l <= i__2; ++l) {
		ff = f;
		y = x[l];
		i__3 = l - 1;
		for (i__ = i1; i__ <= i__3; ++i__) {
		    ff /= y - x[i__];
/* L12: */
		}
		i__3 = i2;
		for (i__ = l + 1; i__ <= i__3; ++i__) {
		    ff /= y - x[i__];
/* L13: */
		}
		we[jj] = ff;
		jj += m2;
/* L14: */
	    }
	}
	ff = f;
	y = x[i2];
	i__2 = i2m1;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ff /= y - x[i__];
/* L16: */
	}
	we[jj] = ff;
	jj += m2;
	jm += inc;
/* L17: */
    }

/* ***  Zero the upper left and lower right corners of E */

    kl = 1;
    n2m = m2p1 * *n + 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ku = kl + *m - i__;
	i__2 = ku;
	for (k = kl; k <= i__2; ++k) {
	    we[k] = 0.;
	    we[n2m - k] = 0.;
/* L18: */
	}
	kl += m2p1;
/* L19: */
    }

/* ***  Weighted matrix WE = W**-1 * E and its L1-norm */

/* L20: */
    jj = 0;
    *el = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__];
	i__2 = m2p1;
	for (j = 1; j <= i__2; ++j) {
	    ++jj;
	    we[jj] /= wi;
	    *el += (d__1 = we[jj], abs(d__1));
/* L21: */
	}
/* L22: */
    }
    *el /= *n;

/* ***  Ready */

    return 0;
} /* prep_ */

/* SPLC.FOR, 1985-12-12 */

/* Author: H.J. Woltring */

/* Organizations: University of Nijmegen, and */
/*                Philips Medical Systems, Eindhoven */
/*                (The Netherlands) */

/* *********************************************************************** */

/* FUNCTION SPLC (REAL*8) */

/* Purpose: */
/* ******* */

/*       To assess the coefficients of a B-spline and various statistical */
/*       parameters, for a given value of the regularization parameter p. */

/* Calling convention: */
/* ****************** */

/*       FV = SPLC ( M, N, K, Y, NY, WX, WY, MODE, VAL, P, EPS, C, NC, */
/*       1           STAT, B, WE, EL, BWE) */

/* Meaning of parameters: */
/* ********************* */

/*       SPLC            ( O )   GCV function value if |MODE|.eq.2, */
/*                               MSE value if |MODE|.eq.3, and absolute */
/*                               difference with the prior given number of */
/*                               degrees of freedom if |MODE|.eq.4. */
/*       M               ( I )   Half order of the B-spline (degree 2*M-1), */
/*                               with M > 0. */
/*       N               ( I )   Number of observations, with N >= 2*M. */
/*       K               ( I )   Number of datasets, with K >= 1. */
/*       Y(NY,K)         ( I )   Observed measurements. */
/*       NY              ( I )   First dimension of Y(NY,K), with NY.ge.N. */
/*       WX(N)           ( I )   Weight factors, corresponding to the */
/*                               relative inverse variance of each measure- */
/*                               ment, with WX(I) > 0.0. */
/*       WY(K)           ( I )   Weight factors, corresponding to the */
/*                               relative inverse variance of each dataset, */
/*                               with WY(J) > 0.0. */
/*       MODE            ( I )   Mode switch, as described in GCVSPL. */
/*       VAL             ( I )   Prior variance if |MODE|.eq.3, and */
/*                               prior number of degrees of freedom if */
/*                               |MODE|.eq.4. For other values of MODE, */
/*                               VAL is not used. */
/*       P               ( I )   Smoothing parameter, with P >= 0.0. If */
/*                               P.eq.0.0, an interpolating spline is */
/*                               calculated. */
/*       EPS             ( I )   Relative rounding tolerance*10.0. EPS is */
/*                               the smallest positive number such that */
/*                               EPS/10.0 + 1.0 .ne. 1.0. */
/*       C(NC,K)         ( O )   Calculated spline coefficient arrays. NB: */
/*                               the dimensions of in GCVSPL and in SPLDER */
/*                               are different! In SPLDER, only a single */
/*                               column of C(N,K) is needed, and the proper */
/*                               column C(1,J), with J=1...K, should be used */
/*                               when calling SPLDER. */
/*       NC              ( I )   First dimension of C(NC,K), with NC.ge.N. */
/*       STAT(6)         ( O )   Statistics array. See the description in */
/*                               subroutine GCVSPL. */
/*       B (1-M:M-1,N)   ( I )   B-spline tableau as evaluated by subroutine */
/*                               BASIS. */
/*       WE( -M:M  ,N)   ( I )   Weighted B-spline tableau (W**-1 * E) as */
/*                               evaluated by subroutine PREP. */
/*       EL              ( I )   L1-norm of the matrix WE as evaluated by */
/*                               subroutine PREP. */
/*       BWE(-M:M,N)     ( O )   Central 2*M+1 bands of the inverted */
/*                               matrix ( B  +  p * W**-1 * E )**-1 */

/* Remarks: */
/* ******* */

/*       This subroutine combines elements of subroutine SPLC0 from the */
/*       paper by Lyche et al. (1983), and of subroutine SPFIT1 by */
/*       Hutchinson (1985). */

/* References: */
/* ********** */

/*       M.F. Hutchinson (1985), Subroutine CUBGCV. CSIRO division of */
/*       Mathematics and Statistics, P.O. Box 1965, Canberra, ACT 2601, */
/*       Australia. */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines */
/*       for computing smoothing and interpolating natural splines. */
/*       Advances in Engineering Software 5(1983)1, pp. 2-5. */

/* *********************************************************************** */

doublereal splc_(integer *m, integer *n, integer *k, doublereal *y, integer *
	ny, doublereal *wx, doublereal *wy, integer *mode, doublereal *val, 
	doublereal *p, doublereal *eps, doublereal *c__, integer *nc, 
	doublereal *stat, doublereal *b, doublereal *we, doublereal *el, 
	doublereal *bwe)
{
    /* System generated locals */
    integer y_dim1, y_offset, c_dim1, c_offset, b_dim1, b_offset, we_dim1, 
	    we_offset, bwe_dim1, bwe_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val=0., d__1;

    /* Local variables */
    static integer i__, j, l;
    static doublereal dp;
    static integer km;
    static doublereal dt;
    static integer kp;
    static doublereal pel, esn, trn;
    extern doublereal trinv_(doublereal *, doublereal *, integer *, integer *)
	    ;
    extern /* Subroutine */ int bandet_(doublereal *, integer *, integer *), 
	    bansol_(doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *);



/* ***  Check on p-value */

    /* Parameter adjustments */
    bwe_dim1 = *m - (-(*m)) + 1;
    bwe_offset = -(*m) + bwe_dim1;
    bwe -= bwe_offset;
    we_dim1 = *m - (-(*m)) + 1;
    we_offset = -(*m) + we_dim1;
    we -= we_offset;
    b_dim1 = *m - 1 - (1 - *m) + 1;
    b_offset = 1 - *m + b_dim1;
    b -= b_offset;
    --wx;
    --wy;
    y_dim1 = *ny;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    c_dim1 = *nc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --stat;

    /* Function Body */
    dp = *p;
    stat[4] = *p;
    pel = *p * *el;
/* ***  Pseudo-interpolation if p is too small */
    if (pel < *eps) {
	dp = *eps / *el;
	stat[4] = 0.;
    }
/* ***  Pseudo least-squares polynomial if p is too large */
    if (pel * *eps > 1.) {
	dp = 1. / (*el * *eps);
	stat[4] = dp;
    }

/* ***  Calculate  BWE  =  B  +  p * W**-1 * E */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = *m, i__3 = i__ - 1;
	km = -min(i__2,i__3);
/* Computing MIN */
	i__2 = *m, i__3 = *n - i__;
	kp = min(i__2,i__3);
	i__2 = kp;
	for (l = km; l <= i__2; ++l) {
	    if (abs(l) == *m) {
		bwe[l + i__ * bwe_dim1] = dp * we[l + i__ * we_dim1];
	    } else {
		bwe[l + i__ * bwe_dim1] = b[l + i__ * b_dim1] + dp * we[l + 
			i__ * we_dim1];
	    }
/* L30: */
	}
/* L40: */
    }

/* ***  Solve BWE * C = Y, and assess TRACE [ B * BWE**-1 ] */

    bandet_(&bwe[bwe_offset], m, n);
    bansol_(&bwe[bwe_offset], &y[y_offset], ny, &c__[c_offset], nc, m, n, k);
    stat[3] = trinv_(&we[we_offset], &bwe[bwe_offset], m, n) * dp;
/* trace * p = res. d.o. */
    trn = stat[3] / *n;

/* ***  Compute mean-squared weighted residual */

    esn = 0.;
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dt = -y[i__ + j * y_dim1];
/* Computing MIN */
	    i__3 = *m - 1, i__4 = i__ - 1;
	    km = -min(i__3,i__4);
/* Computing MIN */
	    i__3 = *m - 1, i__4 = *n - i__;
	    kp = min(i__3,i__4);
	    i__3 = kp;
	    for (l = km; l <= i__3; ++l) {
		dt += b[l + i__ * b_dim1] * c__[i__ + l + j * c_dim1];
/* L50: */
	    }
	    esn += dt * dt * wx[i__] * wy[j];
/* L60: */
	}
/* L70: */
    }
    esn /= *n * *k;

/* ***  Calculate statistics and function value */

    stat[6] = esn / trn;
/* Estimated variance */
    stat[1] = stat[6] / trn;
/* GCV function value */
    stat[2] = esn;
/*     STAT(3) = trace [p*B * BWE**-1] !Estimated residuals' d.o.f. */
/*     STAT(4) = P                     !Normalized smoothing factor */
/* Mean Squared Residual */
    if (abs(*mode) != 3) {
/* ***     Unknown variance: GCV */
	stat[5] = stat[6] - esn;
	if (abs(*mode) == 1) {
	    ret_val = 0.;
	}
	if (abs(*mode) == 2) {
	    ret_val = stat[1];
	}
	if (abs(*mode) == 4) {
	    ret_val = (d__1 = stat[3] - *val, abs(d__1));
	}
    } else {
/* ***     Known variance: estimated mean squared error */
	stat[5] = esn - *val * (trn * 2. - 1.);
	ret_val = stat[5];
    }

    return ret_val;
} /* splc_ */

/* BANDET.FOR, 1985-06-03 */

/* *********************************************************************** */

/* SUBROUTINE BANDET (REAL*8) */

/* Purpose: */
/* ******* */

/*       This subroutine computes the LU decomposition of an N*N matrix */
/*       E. It is assumed that E has M bands above and M bands below the */
/*       diagonal. The decomposition is returned in E. It is assumed that */
/*       E can be decomposed without pivoting. The matrix E is stored in */
/*       vectorized form in the array E(-M:M,N), where element E(J,I) of */
/*       the array E corresponds with element e(i,i+j) of the matrix E. */

/* Calling convention: */
/* ****************** */

/*       CALL BANDET ( E, M, N ) */

/* Meaning of parameters: */
/* ********************* */

/*       E(-M:M,N)       (I/O)   Matrix to be decomposed. */
/*       M, N            ( I )   Matrix dimensioning parameters, */
/*                               M >= 0, N >= 2*M. */

/* Remark: */
/* ****** */

/*       No checking on the validity of the input data is performed. */
/*       If (M.le.0), no action is taken. */

/* *********************************************************************** */

/* Subroutine */ int bandet_(doublereal *e, integer *m, integer *n)
{
    /* System generated locals */
    integer e_dim1, e_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, l;
    static doublereal di, dl;
    static integer mi, km, lm;
    static doublereal du;



    /* Parameter adjustments */
    e_dim1 = *m - (-(*m)) + 1;
    e_offset = -(*m) + e_dim1;
    e -= e_offset;

    /* Function Body */
    if (*m <= 0) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	di = e[i__ * e_dim1];
/* Computing MIN */
	i__2 = *m, i__3 = i__ - 1;
	mi = min(i__2,i__3);
	if (mi >= 1) {
	    i__2 = mi;
	    for (k = 1; k <= i__2; ++k) {
		di -= e[-k + i__ * e_dim1] * e[k + (i__ - k) * e_dim1];
/* L10: */
	    }
	    e[i__ * e_dim1] = di;
	}
/* Computing MIN */
	i__2 = *m, i__3 = *n - i__;
	lm = min(i__2,i__3);
	if (lm >= 1) {
	    i__2 = lm;
	    for (l = 1; l <= i__2; ++l) {
		dl = e[-l + (i__ + l) * e_dim1];
/* Computing MIN */
		i__3 = *m - l, i__4 = i__ - 1;
		km = min(i__3,i__4);
		if (km >= 1) {
		    du = e[l + i__ * e_dim1];
		    i__3 = km;
		    for (k = 1; k <= i__3; ++k) {
			du -= e[-k + i__ * e_dim1] * e[l + k + (i__ - k) * 
				e_dim1];
			dl -= e[-l - k + (l + i__) * e_dim1] * e[k + (i__ - k)
				 * e_dim1];
/* L20: */
		    }
		    e[l + i__ * e_dim1] = du;
		}
		e[-l + (i__ + l) * e_dim1] = dl / di;
/* L30: */
	    }
	}
/* L40: */
    }

/* ***  Ready */

    return 0;
} /* bandet_ */

/* BANSOL.FOR, 1985-12-12 */

/* *********************************************************************** */

/* SUBROUTINE BANSOL (REAL*8) */

/* Purpose: */
/* ******* */

/*       This subroutine solves systems of linear equations given an LU */
/*       decomposition of the design matrix. Such a decomposition is pro- */
/*       vided by subroutine BANDET, in vectorized form. It is assumed */
/*       that the design matrix is not singular. */

/* Calling convention: */
/* ****************** */

/*       CALL BANSOL ( E, Y, NY, C, NC, M, N, K ) */

/* Meaning of parameters: */
/* ********************* */

/*       E(-M:M,N)       ( I )   Input design matrix, in LU-decomposed, */
/*                               vectorized form. Element E(J,I) of the */
/*                               array E corresponds with element */
/*                               e(i,i+j) of the N*N design matrix E. */
/*       Y(NY,K)         ( I )   Right hand side vectors. */
/*       C(NC,K)         ( O )   Solution vectors. */
/*       NY, NC, M, N, K ( I )   Dimensioning parameters, with M >= 0, */
/*                               N > 2*M, and K >= 1. */

/* Remark: */
/* ****** */

/*       This subroutine is an adaptation of subroutine BANSOL from the */
/*       paper by Lyche et al. (1983). No checking is performed on the */
/*       validity of the input parameters and data. Division by zero may */
/*       occur if the system is singular. */

/* Reference: */
/* ********* */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines */
/*       for computing smoothing and interpolating natural splines. */
/*       Advances in Engineering Software 5(1983)1, pp. 2-5. */

/* *********************************************************************** */

/* Subroutine */ int bansol_(doublereal *e, doublereal *y, integer *ny, 
	doublereal *c__, integer *nc, integer *m, integer *n, integer *k)
{
    /* System generated locals */
    integer e_dim1, e_offset, y_dim1, y_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static doublereal d__;
    static integer i__, j, l, mi, nm1;



/* ***  Check on special cases: M=0, M=1, M>1 */

    /* Parameter adjustments */
    e_dim1 = *m - (-(*m)) + 1;
    e_offset = -(*m) + e_dim1;
    e -= e_offset;
    c_dim1 = *nc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    y_dim1 = *ny;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    nm1 = *n - 1;
    if ((i__1 = *m - 1) < 0) {
	goto L10;
    } else if (i__1 == 0) {
	goto L40;
    } else {
	goto L80;
    }

/* ***  M = 0: Diagonal system */

L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
	    c__[i__ + j * c_dim1] = y[i__ + j * y_dim1] / e[i__ * e_dim1];
/* L20: */
	}
/* L30: */
    }
    return 0;

/* ***  M = 1: Tridiagonal system */

L40:
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	c__[j * c_dim1 + 1] = y[j * y_dim1 + 1];
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
/* Forward sweep */
	    c__[i__ + j * c_dim1] = y[i__ + j * y_dim1] - e[i__ * e_dim1 - 1] 
		    * c__[i__ - 1 + j * c_dim1];
/* L50: */
	}
	c__[*n + j * c_dim1] /= e[*n * e_dim1];
	for (i__ = nm1; i__ >= 1; --i__) {
/* Backward sweep */
	    c__[i__ + j * c_dim1] = (c__[i__ + j * c_dim1] - e[i__ * e_dim1 + 
		    1] * c__[i__ + 1 + j * c_dim1]) / e[i__ * e_dim1];
/* L60: */
	}
/* L70: */
    }
    return 0;

/* ***  M > 1: General system */

L80:
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	c__[j * c_dim1 + 1] = y[j * y_dim1 + 1];
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
/* Forward sweep */
/* Computing MIN */
	    i__3 = *m, i__4 = i__ - 1;
	    mi = min(i__3,i__4);
	    d__ = y[i__ + j * y_dim1];
	    i__3 = mi;
	    for (l = 1; l <= i__3; ++l) {
		d__ -= e[-l + i__ * e_dim1] * c__[i__ - l + j * c_dim1];
/* L90: */
	    }
	    c__[i__ + j * c_dim1] = d__;
/* L100: */
	}
	c__[*n + j * c_dim1] /= e[*n * e_dim1];
	for (i__ = nm1; i__ >= 1; --i__) {
/* Backward sweep */
/* Computing MIN */
	    i__2 = *m, i__3 = *n - i__;
	    mi = min(i__2,i__3);
	    d__ = c__[i__ + j * c_dim1];
	    i__2 = mi;
	    for (l = 1; l <= i__2; ++l) {
		d__ -= e[l + i__ * e_dim1] * c__[i__ + l + j * c_dim1];
/* L110: */
	    }
	    c__[i__ + j * c_dim1] = d__ / e[i__ * e_dim1];
/* L120: */
	}
/* L130: */
    }
    return 0;

} /* bansol_ */

/* TRINV.FOR, 1985-06-03 */

/* *********************************************************************** */

/* FUNCTION TRINV (REAL*8) */

/* Purpose: */
/* ******* */

/*       To calculate TRACE [ B * E**-1 ], where B and E are N * N */
/*       matrices with bandwidth 2*M+1, and where E is a regular matrix */
/*       in LU-decomposed form. B and E are stored in vectorized form, */
/*       compatible with subroutines BANDET and BANSOL. */

/* Calling convention: */
/* ****************** */

/*       TRACE = TRINV ( B, E, M, N ) */

/* Meaning of parameters: */
/* ********************* */

/*       B(-M:M,N)       ( I ) Input array for matrix B. Element B(J,I) */
/*                             corresponds with element b(i,i+j) of the */
/*                             matrix B. */
/*       E(-M:M,N)       (I/O) Input array for matrix E. Element E(J,I) */
/*                             corresponds with element e(i,i+j) of the */
/*                             matrix E. This matrix is stored in LU- */
/*                             decomposed form, with L unit lower tri- */
/*                             angular, and U upper triangular. The unit */
/*                             diagonal of L is not stored. Upon return, */
/*                             the array E holds the central 2*M+1 bands */
/*                             of the inverse E**-1, in similar ordering. */
/*       M, N            ( I ) Array and matrix dimensioning parameters */
/*                             (M.gt.0, N.ge.2*M+1). */
/*       TRINV           ( O ) Output function value TRACE [ B * E**-1 ] */

/* Reference: */
/* ********* */

/*       A.M. Erisman & W.F. Tinney, On computing certain elements of the */
/*       inverse of a sparse matrix. Communications of the ACM 18(1975), */
/*       nr. 3, pp. 177-179. */

/* *********************************************************************** */

doublereal trinv_(doublereal *b, doublereal *e, integer *m, integer *n)
{
    /* System generated locals */
    integer b_dim1, b_offset, e_dim1, e_offset, i__1, i__2, i__3;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, k;
    static doublereal dd, dl;
    static integer mi;
    static doublereal du;
    static integer mn, mp;



/* ***  Assess central 2*M+1 bands of E**-1 and store in array E */

    /* Parameter adjustments */
    e_dim1 = *m - (-(*m)) + 1;
    e_offset = -(*m) + e_dim1;
    e -= e_offset;
    b_dim1 = *m - (-(*m)) + 1;
    b_offset = -(*m) + b_dim1;
    b -= b_offset;

    /* Function Body */
    e[*n * e_dim1] = 1. / e[*n * e_dim1];
/* Nth pivot */
    for (i__ = *n - 1; i__ >= 1; --i__) {
/* Computing MIN */
	i__1 = *m, i__2 = *n - i__;
	mi = min(i__1,i__2);
	dd = 1. / e[i__ * e_dim1];
/* ***     Save Ith column of L and Ith row of U, and normalize U row */
/* Ith pivot */
	i__1 = mi;
	for (k = 1; k <= i__1; ++k) {
	    e[k + *n * e_dim1] = e[k + i__ * e_dim1] * dd;
/* Ith row of U (normalized) */
	    e[-k + e_dim1] = e[-k + (k + i__) * e_dim1];
/* Ith column of L */
/* L10: */
	}
	dd += dd;
/* ***     Invert around Ith pivot */
	for (j = mi; j >= 1; --j) {
	    du = 0.;
	    dl = 0.;
	    i__1 = mi;
	    for (k = 1; k <= i__1; ++k) {
		du -= e[k + *n * e_dim1] * e[j - k + (i__ + k) * e_dim1];
		dl -= e[-k + e_dim1] * e[k - j + (i__ + j) * e_dim1];
/* L20: */
	    }
	    e[j + i__ * e_dim1] = du;
	    e[-j + (j + i__) * e_dim1] = dl;
	    dd -= e[j + *n * e_dim1] * dl + e[-j + e_dim1] * du;
/* L30: */
	}
	e[i__ * e_dim1] = dd * .5;
/* L40: */
    }

/* ***  Assess TRACE [ B * E**-1 ] and clear working storage */

    dd = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = *m, i__3 = i__ - 1;
	mn = -min(i__2,i__3);
/* Computing MIN */
	i__2 = *m, i__3 = *n - i__;
	mp = min(i__2,i__3);
	i__2 = mp;
	for (k = mn; k <= i__2; ++k) {
	    dd += b[k + i__ * b_dim1] * e[-k + (k + i__) * e_dim1];
/* L50: */
	}
/* L60: */
    }
    ret_val = dd;
    i__1 = *m;
    for (k = 1; k <= i__1; ++k) {
	e[k + *n * e_dim1] = 0.;
	e[-k + e_dim1] = 0.;
/* L70: */
    }

/* ***  Ready */

    return ret_val;
} /* trinv_ */

/* SPLDER.FOR, 1985-06-11 */

/* *********************************************************************** */

/* FUNCTION SPLDER (REAL*8) */

/* Purpose: */
/* ******* */

/*       To produce the value of the function (IDER.eq.0) or of the */
/*       IDERth derivative (IDER.gt.0) of a 2M-th order B-spline at */
/*       the point T. The spline is described in terms of the half */
/*       order M, the knot sequence X(N), N.ge.2*M, and the spline */
/*       coefficients C(N). */

/* Calling convention: */
/* ****************** */

/*       SVIDER = SPLDER ( IDER, M, N, T, X, C, L, Q ) */

/* Meaning of parameters: */
/* ********************* */

/*       SPLDER  ( O )   Function or derivative value. */
/*       IDER    ( I )   Derivative order required, with 0.le.IDER */
/*                       and IDER.le.2*M. If IDER.eq.0, the function */
/*                       value is returned; otherwise, the IDER-th */
/*                       derivative of the spline is returned. */
/*       M       ( I )   Half order of the spline, with M.gt.0. */
/*       N       ( I )   Number of knots and spline coefficients, */
/*                       with N.ge.2*M. */
/*       T       ( I )   Argument at which the spline or its deri- */
/*                       vative is to be evaluated, with X(1).le.T */
/*                       and T.le.X(N). */
/*       X(N)    ( I )   Strictly increasing knot sequence array, */
/*                       X(I-1).lt.X(I), I=2,...,N. */
/*       C(N)    ( I )   Spline coefficients, as evaluated by */
/*                       subroutine GVCSPL. */
/*       L       (I/O)   L contains an integer such that: */
/*                       X(L).le.T and T.lt.X(L+1) if T is within */
/*                       the range X(1).le.T and T.lt.X(N). If */
/*                       T.lt.X(1), L is set to 0, and if T.ge.X(N), */
/*                       L is set to N. The search for L is facili- */
/*                       tated if L has approximately the right */
/*                       value on entry. */
/*       Q(2*M)  ( W )   Internal work array. */

/* Remark: */
/* ****** */

/*       This subroutine is an adaptation of subroutine SPLDER of */
/*       the paper by Lyche et al. (1983). No checking is performed */
/*       on the validity of the input parameters. */

/* Reference: */
/* ********* */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines */
/*       for computing smoothing and interpolating natural splines. */
/*       Advances in Engineering Software 5(1983)1, pp. 2-5. */

/* *********************************************************************** */

doublereal splder_(integer *ider, integer *m, integer *n, doublereal *t, 
	doublereal *x, doublereal *c__, integer *l, doublereal *q)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer i__, j, k;
    static doublereal z__;
    static integer i1, j1, k1, j2, m2, ii, jj, ki, jl, lk, mi, nk, lm, ml, jm,
	     ir, ju;
    static doublereal tt;
    static integer lk1, mp1, m2m1, jin, nki, npm, lk1i, nki1, lk1i1;
    static doublereal xjki;
    extern /* Subroutine */ int search_(integer *, doublereal *, doublereal *,
	     integer *);



/* ***  Derivatives of IDER.ge.2*M are alway zero */

    /* Parameter adjustments */
    --q;
    --c__;
    --x;

    /* Function Body */
    m2 = *m << 1;
    k = m2 - *ider;
    if (k < 1) {
	ret_val = 0.;
	return ret_val;
    }

/* ***  Search for the interval value L */

    search_(n, &x[1], t, l);

/* ***  Initialize parameters and the 1st row of the B-spline */
/* ***  coefficients tableau */

    tt = *t;
    mp1 = *m + 1;
    npm = *n + *m;
    m2m1 = m2 - 1;
    k1 = k - 1;
    nk = *n - k;
    lk = *l - k;
    lk1 = lk + 1;
    lm = *l - *m;
    jl = *l + 1;
    ju = *l + m2;
    ii = *n - m2;
    ml = -(*l);
    i__1 = ju;
    for (j = jl; j <= i__1; ++j) {
	if (j >= mp1 && j <= npm) {
	    q[j + ml] = c__[j - *m];
	} else {
	    q[j + ml] = 0.;
	}
/* L2: */
    }

/* ***  The following loop computes differences of the B-spline */
/* ***  coefficients. If the value of the spline is required, */
/* ***  differencing is not necessary. */

    if (*ider > 0) {
	jl -= m2;
	ml += m2;
	i__1 = *ider;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ++jl;
	    ++ii;
	    j1 = max(1,jl);
	    j2 = min(*l,ii);
	    mi = m2 - i__;
	    j = j2 + 1;
	    if (j1 <= j2) {
		i__2 = j2;
		for (jin = j1; jin <= i__2; ++jin) {
		    --j;
		    jm = ml + j;
		    q[jm] = (q[jm] - q[jm - 1]) / (x[j + mi] - x[j]);
/* L3: */
		}
	    }
	    if (jl >= 1) {
		goto L6;
	    }
	    i1 = i__ + 1;
	    j = ml + 1;
	    if (i1 <= ml) {
		i__2 = ml;
		for (jin = i1; jin <= i__2; ++jin) {
		    --j;
		    q[j] = -q[j - 1];
/* L5: */
		}
	    }
L6:
	    ;
	}
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    q[j] = q[j + *ider];
/* L7: */
	}
    }

/* ***  Compute lower half of the evaluation tableau */

    if (k1 >= 1) {
/* Tableau ready if IDER.eq.2*M-1 */
	i__1 = k1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    nki = nk + i__;
	    ir = k;
	    jj = *l;
	    ki = k - i__;
	    nki1 = nki + 1;
/* ***        Right-hand B-splines */
	    if (*l >= nki1) {
		i__2 = *l;
		for (j = nki1; j <= i__2; ++j) {
		    q[ir] = q[ir - 1] + (tt - x[jj]) * q[ir];
		    --jj;
		    --ir;
/* L9: */
		}
	    }
/* ***        Middle B-splines */
	    lk1i = lk1 + i__;
	    j1 = max(1,lk1i);
	    j2 = min(*l,nki);
	    if (j1 <= j2) {
		i__2 = j2;
		for (j = j1; j <= i__2; ++j) {
		    xjki = x[jj + ki];
		    z__ = q[ir];
		    q[ir] = z__ + (xjki - tt) * (q[ir - 1] - z__) / (xjki - x[
			    jj]);
		    --ir;
		    --jj;
/* L11: */
		}
	    }
/* ***        Left-hand B-splines */
	    if (lk1i <= 0) {
		jj = ki;
		lk1i1 = 1 - lk1i;
		i__2 = lk1i1;
		for (j = 1; j <= i__2; ++j) {
		    q[ir] += (x[jj] - tt) * q[ir - 1];
		    --jj;
		    --ir;
/* L13: */
		}
	    }
/* L14: */
	}
    }

/* ***  Compute the return value */

    z__ = q[k];
/* ***  Multiply with factorial if IDER.gt.0 */
    if (*ider > 0) {
	i__1 = m2m1;
	for (j = k; j <= i__1; ++j) {
	    z__ *= j;
/* L16: */
	}
    }
    ret_val = z__;

/* ***  Ready */

    return ret_val;
} /* splder_ */

/* SEARCH.FOR, 1985-06-03 */

/* *********************************************************************** */

/* SUBROUTINE SEARCH (REAL*8) */

/* Purpose: */
/* ******* */

/*       Given a strictly increasing knot sequence X(1) < ... < X(N), */
/*       where N >= 1, and a real number T, this subroutine finds the */
/*       value L such that X(L) <= T < X(L+1).  If T < X(1), L = 0; */
/*       if X(N) <= T, L = N. */

/* Calling convention: */
/* ****************** */

/*       CALL SEARCH ( N, X, T, L ) */

/* Meaning of parameters: */
/* ********************* */

/*       N       ( I )   Knot array dimensioning parameter. */
/*       X(N)    ( I )   Stricly increasing knot array. */
/*       T       ( I )   Input argument whose knot interval is to */
/*                       be found. */
/*       L       (I/O)   Knot interval parameter. The search procedure */
/*                       is facilitated if L has approximately the */
/*                       right value on entry. */

/* Remark: */
/* ****** */

/*       This subroutine is an adaptation of subroutine SEARCH from */
/*       the paper by Lyche et al. (1983). No checking is performed */
/*       on the input parameters and data; the algorithm may fail if */
/*       the input sequence is not strictly increasing. */

/* Reference: */
/* ********* */

/*       T. Lyche, L.L. Schumaker, & K. Sepehrnoori, Fortran subroutines */
/*       for computing smoothing and interpolating natural splines. */
/*       Advances in Engineering Software 5(1983)1, pp. 2-5. */

/* *********************************************************************** */

/* Subroutine */ int search_(integer *n, doublereal *x, doublereal *t, 
	integer *l)
{
    static integer il, iu;



    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*t < x[1]) {
/* ***     Out of range to the left */
	*l = 0;
	return 0;
    }
    if (*t >= x[*n]) {
/* ***     Out of range to the right */
	*l = *n;
	return 0;
    }
/* ***  Validate input value of L */
    *l = max(*l,1);
    if (*l >= *n) {
	*l = *n - 1;
    }

/* ***  Often L will be in an interval adjoining the interval found */
/* ***  in a previous call to search */

    if (*t >= x[*l]) {
	goto L5;
    }
    --(*l);
    if (*t >= x[*l]) {
	return 0;
    }

/* ***  Perform bisection */

    il = 1;
L3:
    iu = *l;
L4:
    *l = (il + iu) / 2;
    if (iu - il <= 1) {
	return 0;
    }
    if (*t < x[*l]) {
	goto L3;
    }
    il = *l;
    goto L4;
L5:
    if (*t < x[*l + 1]) {
	return 0;
    }
    ++(*l);
    if (*t < x[*l + 1]) {
	return 0;
    }
    il = *l + 1;
    iu = *n;
    goto L4;

} /* search_ */

#ifdef __cplusplus
	}
#endif
