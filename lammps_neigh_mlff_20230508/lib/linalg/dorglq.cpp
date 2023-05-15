#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
int dorglq_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau,
            doublereal *work, integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern int dorgl2_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *),
        dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, ftnlen, ftnlen, ftnlen, ftnlen),
        dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer ldwork, lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    nb = ilaenv_(&c__1, (char *)"DORGLQ", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = max(1, *m) * nb;
    work[1] = (doublereal)lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < *m) {
        *info = -2;
    } else if (*k < 0 || *k > *m) {
        *info = -3;
    } else if (*lda < max(1, *m)) {
        *info = -5;
    } else if (*lwork < max(1, *m) && !lquery) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORGLQ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*m <= 0) {
        work[1] = 1.;
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {
        i__1 = 0, i__2 = ilaenv_(&c__3, (char *)"DORGLQ", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
        nx = max(i__1, i__2);
        if (nx < *k) {
            ldwork = *m;
            iws = ldwork * nb;
            if (*lwork < iws) {
                nb = *lwork / ldwork;
                i__1 = 2,
                i__2 = ilaenv_(&c__2, (char *)"DORGLQ", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
                nbmin = max(i__1, i__2);
            }
        }
    }
    if (nb >= nbmin && nb < *k && nx < *k) {
        ki = (*k - nx - 1) / nb * nb;
        i__1 = *k, i__2 = ki + nb;
        kk = min(i__1, i__2);
        i__1 = kk;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = kk + 1; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] = 0.;
            }
        }
    } else {
        kk = 0;
    }
    if (kk < *m) {
        i__1 = *m - kk;
        i__2 = *n - kk;
        i__3 = *k - kk;
        dorgl2_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &tau[kk + 1], &work[1],
                &iinfo);
    }
    if (kk > 0) {
        i__1 = -nb;
        for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
            i__2 = nb, i__3 = *k - i__ + 1;
            ib = min(i__2, i__3);
            if (i__ + ib <= *m) {
                i__2 = *n - i__ + 1;
                dlarft_((char *)"Forward", (char *)"Rowwise", &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__],
                        &work[1], &ldwork, (ftnlen)7, (ftnlen)7);
                i__2 = *m - i__ - ib + 1;
                i__3 = *n - i__ + 1;
                dlarfb_((char *)"Right", (char *)"Transpose", (char *)"Forward", (char *)"Rowwise", &i__2, &i__3, &ib,
                        &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork, &a[i__ + ib + i__ * a_dim1],
                        lda, &work[ib + 1], &ldwork, (ftnlen)5, (ftnlen)9, (ftnlen)7, (ftnlen)7);
            }
            i__2 = *n - i__ + 1;
            dorgl2_(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
            i__2 = i__ - 1;
            for (j = 1; j <= i__2; ++j) {
                i__3 = i__ + ib - 1;
                for (l = i__; l <= i__3; ++l) {
                    a[l + j * a_dim1] = 0.;
                }
            }
        }
    }
    work[1] = (doublereal)iws;
    return 0;
}
#ifdef __cplusplus
}
#endif
