#include "blas.h"

#ifdef __cplusplus
extern "C" {
#endif

/* int dgemv_(char *trans, long int *m, long int *n,  */
/* 	   double* alpha, double *a, long int *lda,  */
/* 	   double *x, long int *incx,  double *beta,  */
/* 	   double *y, long int *incy) */
int dgemv_(char *trans, int *m, int *n, 
	   double* alpha, double *a, int *lda, 
	   double *x, int *incx,  double *beta, 
	   double *y, int *incy)
{
    /* System generated locals */
    long int a_dim1, a_offset, i1, i2;

    /* Local variables */
    long int i, j, ix, iy, jx, jy, kx, ky, info;
    double temp;
    long int lenx, leny;
    /* extern long int lsame_(char *, char *); */
    /* extern /\* Subroutine *\/ int xerbla_(char *, long int *); */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    /* if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C") */
    if (! ((*trans == 'N') || (*trans == 'n'))&& 
	! ((*trans == 'T') || (*trans == 't'))&& 
	! ((*trans == 'C') || (*trans == 'c'))
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < MAX(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	/* xerbla_("DGEMV ", &info); */
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.0 && *beta == 1.0) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    /* if (lsame_(trans, "N")) { */
    if ((*trans == 'N') || (*trans == 'n')) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i1 = leny;
		for (i = 1; i <= i1; ++i) {
		    y[i] = 0.;
/* L10: */
		}
	    } else {
		i1 = leny;
		for (i = 1; i <= i1; ++i) {
		    y[i] = *beta * y[i];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i1 = leny;
		for (i = 1; i <= i1; ++i) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i1 = leny;
		for (i = 1; i <= i1; ++i) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    /* if (lsame_(trans, "N")) { */
    if ((*trans == 'N') || (*trans == 'n')) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i1 = *n;
	    for (j = 1; j <= i1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    i2 = *m;
		    for (i = 1; i <= i2; ++i) {
			y[i] += temp * a[i + j * a_dim1];
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i1 = *n;
	    for (j = 1; j <= i1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i2 = *m;
		    for (i = 1; i <= i2; ++i) {
			y[iy] += temp * a[i + j * a_dim1];
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i1 = *n;
	    for (j = 1; j <= i1; ++j) {
		temp = 0.;
		i2 = *m;
		for (i = 1; i <= i2; ++i) {
		    temp += a[i + j * a_dim1] * x[i];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i1 = *n;
	    for (j = 1; j <= i1; ++j) {
		temp = 0.;
		ix = kx;
		i2 = *m;
		for (i = 1; i <= i2; ++i) {
		    temp += a[i + j * a_dim1] * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;
} /* dgemv_ */
#ifdef __cplusplus
}
#endif
