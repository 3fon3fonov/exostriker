/*  -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/***************************************************************************/
/*                        SYMBA5_GETACCH.F */
/***************************************************************************/
/* This subroutine calculates the acceleration on the massive particles */
/* in the HELIOCENTRIC frame. */
/*             Input: */
/*                 nbod        ==>  number of massive bodies (int scalar) */
/*                nbodm       ==>  Location of last massive body(int scalar)*/
/*                 mass        ==>  mass of bodies (real array) */
/*                 j2rp2,j4rp4 ==>  J2*radii_pl^2 and  J4*radii_pl^4 */
/*                                     (real scalars) */
/*                xh,yh,zh    ==>  position in heliocentric coord (real arrays
)*/
/*                 mtiny       ==>  Small mass  (real array) */
/*                ielc           ==>  number of encounters (int scalar) */
/*               ielst          ==>  list of ecnounters (2D integer*2 array)*/
/*             Output: */
/*                axh,ayh,azh ==>  acceleration in helio coord (real arrays)*/

/* Remarks: Based on helio_getacch.f, but does not include the forces of */
/*          an body B on body A, if body B and A are having an encounter. */
/* Author:  Hal Levison */
/* Date:    3/20/97 */
/* Last revision: 5/28/99 */
/* Subroutine */ int symba5_getacch__(nbod, nbodm, mass, j2rp2, j4rp4, xh, yh,
	 zh, axh, ayh, azh, mtiny, ielc, ielst)
integer *nbod, *nbodm;
doublereal *mass, *j2rp2, *j4rp4, *xh, *yh, *zh, *axh, *ayh, *azh, *mtiny;
integer *ielc;
shortint *ielst;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    doublereal aoblx[4096], aobly[4096], aoblz[4096];
    doublereal ir3h[4096], irh[4096];
    integer i__, j, ie;
    doublereal dx, dy, dz, rji2, faci, facj, irij3;
    doublereal massi, xhi, yhi, zhi, axhi, ayhi, azhi;
    doublereal rh[12288], ah[12288];
    extern /* Subroutine */ int getacch_ir3__();
    extern /* Subroutine */ int obl_acc__();
    doublereal *pi, *pj;

/************************************************************************
***/
/*                        SWIFT.INC */
/************************************************************************
***/
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* ...   Version of Swift */
/* you got it baby */
/* ...   Maximum array size */
/*      parameter  (NPLMAX = 21)   ! max number of planets, including the 
Sun*/
/* max number of planets, including the */
/* ...   Size of the test particle integer status flag */
/* max number of test particles */
/* Number of status parameters */
/* Number of status parameters */
/* ...   Size of the test particle integer status flag */
/* include one for @ pl */
/* ...   convergence criteria for danby */
/* io_init_tp assumes NSTAT==NSTATR */
/* ...    loop limits in the Laguerre attempts */
/* ...    A small number */
/* ...    trig stuff */
/*-----------------------------------------------------------------------
--*/
/* ...  Inputs: */
/************************************************************************
***/
/*                        SYMBA5.INC */
/************************************************************************
***/
/* Include file for the SyMBA subroutines */

/* Remarks:  copied from symba.inc */
/* Author:  Hal Levison */
/* Date:    3/20/97 */
/* Last revision: */
/* ...  Maximum number of encounters */
/* ...	scale factor for hill's sphere to take shorter time step */
/* ...   Ratio of shell radii squared */
/* ..    ratio of the number of time steps in the adjoining shells */
/* RSHELL ~ NTENC^(-2/3) */
/* ...  Outputs: */
/* ...  Internals: */
/* --- */
/* ...  Executable code */
/* ...  Initialize things */
    /* Parameter adjustments */
    --azh;
    --ayh;
    --axh;
    --zh;
    --yh;
    --xh;
    --mass;
    ielst -= 3;

    /* Function Body */
    i__1 = *nbod;
    for (i__ = 1; i__ <= i__1; ++i__) {
        pi = &rh[i__ * 3 - 3];
	*pi     = xh[i__];
	*(++pi) = yh[i__];
	*(++pi) = zh[i__];
        pi = &ah[i__ * 3 - 3];
	*pi     = (float)0.;
	*(++pi) = (float)0.;
	*(++pi) = (float)0.;
    }
/* ...  now the third terms */
    i__1 = *nbodm;
    for (i__ = 2; i__ <= i__1; ++i__) {
	massi = mass[i__];
	pi = &rh[i__ * 3 - 3];
	xhi = *pi;
	yhi = *(++pi);
	zhi = *(++pi);
	pi = &ah[i__ * 3 - 3];
	axhi = *pi;
	ayhi = *(++pi);
	azhi = *(++pi);
	i__2 = *nbod;
	for (j = i__ + 1; j <= i__2; ++j) {
	    pj = &rh[j * 3 - 3];
	    dx = *pj     - xhi;
	    dy = *(++pj) - yhi;
	    dz = *(++pj) - zhi;
	    rji2 = dx * dx + dy * dy + dz * dz;
	    irij3 = 1. / (rji2 * sqrt(rji2));
	    faci = massi * irij3;
	    facj = mass[j] * irij3;
	    pj = &ah[j * 3 - 3];
	    *pj     -= faci * dx;
	    *(++pj) -= faci * dy;
	    *(++pj) -= faci * dz;
	    axhi += facj * dx;
	    ayhi += facj * dy;
	    azhi += facj * dz;
	}
	pi = &ah[i__ * 3 - 3];
	*pi     = axhi;
	*(++pi) = ayhi;
	*(++pi) = azhi;
    }
/* ...  Now subtract off anyone in an encounter */
    i__1 = *ielc;
    for (ie = 1; ie <= i__1; ++ie) {
	i__ = ielst[(ie << 1) + 1];
	j = ielst[(ie << 1) + 2];
	dx = rh[j * 3 - 3] - rh[i__ * 3 - 3];
	dy = rh[j * 3 - 2] - rh[i__ * 3 - 2];
	dz = rh[j * 3 - 1] - rh[i__ * 3 - 1];
	rji2 = dx * dx + dy * dy + dz * dz;
	irij3 = 1. / (rji2 * sqrt(rji2));
	faci = mass[i__] * irij3;
	facj = mass[j] * irij3;
	ah[j * 3 - 3] += faci * dx;
	ah[j * 3 - 2] += faci * dy;
	ah[j * 3 - 1] += faci * dz;
	ah[i__ * 3 - 3] -= facj * dx;
	ah[i__ * 3 - 2] -= facj * dy;
	ah[i__ * 3 - 1] -= facj * dz;
    }
/* ...  Now do j2 and j4 stuff */
    if (*j2rp2 != 0.) {
	getacch_ir3__(nbod, &c__2, &xh[1], &yh[1], &zh[1], ir3h, irh);
	obl_acc__(nbod, &mass[1], j2rp2, j4rp4, &xh[1], &yh[1], &zh[1], irh, 
		aoblx, aobly, aoblz);
	i__1 = *nbod;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if (mass[i__] != 0.) {
		ah[i__ * 3 - 3] = ah[i__ * 3 - 3] + aoblx[i__ - 1] - aoblx[0];
		ah[i__ * 3 - 2] = ah[i__ * 3 - 2] + aobly[i__ - 1] - aobly[0];
		ah[i__ * 3 - 1] = ah[i__ * 3 - 1] + aoblz[i__ - 1] - aoblz[0];
	    }
	}
    }
    i__1 = *nbod;
    for (i__ = 1; i__ <= i__1; ++i__) {
	axh[i__] = ah[i__ * 3 - 3];
	ayh[i__] = ah[i__ * 3 - 2];
	azh[i__] = ah[i__ * 3 - 1];
    }
    return 0;
} /* symba5_getacch__ */

