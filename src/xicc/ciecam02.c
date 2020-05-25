//
//  ciecam02.c
//  hargyllcms-1.9.2
//
//  Created by Teunis on 14/11/2019.
//  Copyright © 2019 Teunis Polak. All rights reserved.
//

/*
 *
 * This file is based on cam02.c which is in turn based on cam97s3.c, both by Graeme Gill.
 *
 * Copyright 2004 - 2011 Graeme W. Gill
 * Please refer to COPYRIGHT file for details.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* Note that XYZ values are normalised to 1.0 consistent */
/* with the ICC convention (not 100.0 as assumed by the CIECAM spec.) */
/* Note that all whites are assumed to be normalised (ie. Y = 1.0) */

#include <copyright.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xcam.h"
#include "ciecam02.h"
#include "numlib.h"
#include "cam02.h"

//struct Correlate corr;

static void cam_free(ciecam02 *s);
static int  set_view(ciecam02 *s, double Wxyz[3], double La, double Yb, double Lv);
static int  XYZ_to_cam(ciecam02 *s, double *xyz, double *Jab);
static int  cam_to_XYZ(ciecam02 *s, double *xyz, double *Jab);


/* Create a ciecam02 object */
ciecam02 *new_ciecam02(void) {
	ciecam02 *s;

	if ((s = (ciecam02 *)calloc(1, sizeof(ciecam02))) == NULL) {
		fprintf(stderr,"ciecam02: malloc failed allocating object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del          = cam_free;
	s->set_view     = set_view;
	s->XYZ_to_cam   = XYZ_to_cam;
	s->cam_to_XYZ   = cam_to_XYZ;
    
	return s;
}


/* Delete a ciecam02 object */
static void cam_free(ciecam02 *s) {

	if (s != NULL)
		free(s);
}


/* Set viewing parameters, all white related parameters
 * for ciecam02 object. Return value is always 0 */
static int set_view
(   ciecam02 *s,
    double Wxyz[3],	        /* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
    double La,		        /* Adapting/Surround Luminance cd/m^2 */
    double Yb,		        /* Relative Luminance of Background to reference white (range 0.0 .. 1.0) */
    double Lv	    )	    /* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
                            /* Ignored if Ev is set to anything else than vc_none */
{
	int i;
    double tt;
    
    /* Transfer input parameters to the object */
    s->Wxyz[0] = Wxyz[0];
    s->Wxyz[1] = Wxyz[1];
    s->Wxyz[2] = Wxyz[2];
    s->La = La;
    s->Yb = Yb > 0.005 ? Yb : 0.005;    /* Set minimum to avoid divide by 0.0 */
    s->Lv = Lv;

    
    /* Surround parameters F, c, Nc */
    
//    /* In case surround condition Ev is given */
//    if (Ev != vc_none) {
//        switch(Ev) {
//            case vc_dark:
//                s->C    = 0.525;
//                s->Nc   = 0.8;
//                s->F    = 0.8;
//                break;
//            case vc_dim:
//                s->C    = 0.59;
//                s->Nc   = 0.95;
//                s->F    = 0.9;
//                break;
//            case vc_average:
//            default:
//                s->C    = 0.69;
//                s->Nc   = 1.0;
//                s->F    = 1.0;
//                break;
//        }
//    }
    
    /* In case La and Lv given (and Ev == 0 or vc_none) */
    double r, bf;
    /*                 Dark   Dim   Avg   >Avg */
    double t_C[4]  = { 0.525, 0.59, 0.69, 1.0  };
    double t_Nc[4] = { 0.800, 0.95, 1.00, 1.0  };
    double t_F[4]  = { 0.800, 0.90, 1.00, 1.0  };

    if (La < 1e-10) 		                                                                /* Hmm. */
        La = 1e-10;
    
    /* Force La/Lv to range [0.0 .. 1.0] */
    r = La/Lv;
    if (r < 0.0)
        r = 0.0;
    else if (r > 1.0)
        r = 1.0;
    
    /* Interpolate values for F, c and Nc based on La/Lv */
    if (r < 0.1) {			                                                                /* Dark to Dim */
        i = 0;
        bf = r/0.1;
    } else if (r < 0.2) {	                                                                /* Dim to Average */
        i = 1;
        bf = (r - 0.1)/0.1;
    } else {				                                                                /* Average to above average */
        i = 2;
        bf = (r - 0.2)/0.8;
    }
    s->C  = t_C[i]  * (1.0 - bf) + t_C[i+1]  * bf;
    s->Nc = t_Nc[i] * (1.0 - bf) + t_Nc[i+1] * bf;
    s->F  = t_F[i]  * (1.0 - bf) + t_F[i+1]  * bf;

    
	/* Lightness contrast factor Fl */
		double k = 1.0 / (5.0 * s->La + 1.0);                                               /* eqn (1) */
		s->Fl = 0.2 * pow(k, 4.0) * 5.0 * s->La
		      + 0.1 * pow(1.0 - pow(k, 4.0), 2.0) * pow(5.0 * s->La, 1.0/3.0);              /* eqn (2) */

    
    /* Background induction factor n */
    s->n = s->Yb/ s->Wxyz[1];                                                               /* eqn (3) */

    
	/* Background induction factor Nbb */
	s->Nbb = 0.725 * pow(1.0/s->n, 0.2);                                                    /* eqn (4) */
    s->Ncb = s->Nbb;

    
	/* Base exponential nonlinearity z */
	s->z = 1.48 + pow(s->n , 0.5);                                                          /* eqn (5) */
    
    
    /* XYZ to CAT02 transformed white values rgbW */                                        /* eqn (6), (7) */
    s->rgbW[0] =  0.7328 * s->Wxyz[0] + 0.4296 * s->Wxyz[1] - 0.1624 * s->Wxyz[2];
    s->rgbW[1] = -0.7036 * s->Wxyz[0] + 1.6975 * s->Wxyz[1] + 0.0061 * s->Wxyz[2];
    s->rgbW[2] =  0.0030 * s->Wxyz[0] + 0.0136 * s->Wxyz[1] + 0.9834 * s->Wxyz[2];
    
    
    /* Degree of chromatic adaptation D */
    s->D = s->F * (1.0 - exp((-s->La - 42.0)/92.0)/3.6);                                    /* eqn (8) */
    
    
    /* Intermediate chromatic transform white values Drgb */                                /* eqn (9) */
    s->Drgb[0] = s->D * (s->Wxyz[1]/s->rgbW[0]) + 1.0 - s->D;
    s->Drgb[1] = s->D * (s->Wxyz[1]/s->rgbW[1]) + 1.0 - s->D;
    s->Drgb[2] = s->D * (s->Wxyz[1]/s->rgbW[2]) + 1.0 - s->D;
    
    
    /* Chromaticaly transformed white values rgbcW */
    s->rgbcW[0] = s->Drgb[0] * s->rgbW[0];
    s->rgbcW[1] = s->Drgb[1] * s->rgbW[1];
    s->rgbcW[2] = s->Drgb[2] * s->rgbW[2];
    
    
    /* CAT02 to HPE transformed sample values rgbpW */                                      /* eqn (10), (11), (12) */
    s->rgbpW[0] =  0.7409792 * s->rgbcW[0] + 0.2180250 * s->rgbcW[1] + 0.0410058 * s->rgbcW[2];
    s->rgbpW[1] =  0.2853532 * s->rgbcW[0] + 0.6242014 * s->rgbcW[1] + 0.0904454 * s->rgbcW[2];
    s->rgbpW[2] = -0.0096280 * s->rgbcW[0] - 0.0056980 * s->rgbcW[1] + 1.0153260 * s->rgbcW[2];


	/* Non-linear adapted white values rgbaW */                                             /* eqn (13) */
    for (i = 0; i < 3; i++) {
        tt = pow(s->Fl * s->rgbpW[i], 0.42);
        s->rgbaW[i] = 400.0 * tt / (tt + 27.13) + 0.1;
    }

    
    /* Achromatic response of white Aw */
    s->Aw = (2.0 * s->rgbaW[0] + s->rgbaW[1] + (1.0/20.0) * s->rgbaW[2] - 0.305) * s->Nbb;      /* eqn (20), for white */
    
    
    /* Intermediate parameter for use in eqn (23) nn */
    s->nn = pow(1.64 - pow(0.29, s->n), 0.73);
    
    
    if (printCAM == 1) {
        printf("Scene parameters (CIECAM02 by the book):\n");
        printf("Ref white Wxyz = %f %f %f\n", s->Wxyz[0], s->Wxyz[1], s->Wxyz[2]);
        printf("Relative luminance of background Yb = %f\n", s->Yb);
        printf("Adapting luminance La = %f\n", s->La);

        printf("Internal parameters:\n");
        printf("Surround Impact C = %f\n", s->C);
        printf("Chromatic Induction Nc = %f\n", s->Nc);
        printf("Adaptation Degree F = %f\n", s->F);

        printf("Pre-computed values\n");
        printf("Sharpened cone white rgbW = %f %f %f\n", s->rgbW[0], s->rgbW[1], s->rgbW[2]);
        printf("Degree of chromatic adaptation D = %f\n", s->D);
        printf("Chromatic transform values Drgb = %f %f %f\n", s->Drgb[0], s->Drgb[1], s->Drgb[2]);
        printf("Chromatically transformed white rgbcW = %f %f %f\n", s->rgbcW[0], s->rgbcW[1], s->rgbcW[2]);
        printf("Hunter-P-E cone response white rgbpW = %f %f %f\n", s->rgbpW[0], s->rgbpW[1], s->rgbpW[2]);
        printf("Background induction factor n = %f\n", s->n);
        printf("Lightness contrast factor Fl = %f\n", s->Fl);
        printf("Background brightness induction factor Nbb = %f\n", s->Nbb);
        printf("Base exponential nonlinearity z = %f\n", s->z);
        printf("Post adapted cone response white rgbaW = %f %f %f\n", s->rgbaW[0], s->rgbaW[1], s->rgbaW[2]);
        printf("Achromatic response of white Aw = %f\n", s->Aw);
        printf("\n");
    }
    
    return 0;
}


/* Forward transform. Return value is always 0 */
static int XYZ_to_cam
(   struct _ciecam02 *s,
    double Jab[3],
    double XYZ[3]        )
{
	int i;
	double rgb[3], rgbc[3], rgbp[3], rgba[3];
    double xyz[3];
	double a, b, J, C, h, H, e, A, M, Q, rS, sat;
    double ttA;
    double ttd, tt, t, temp;

    
    /* XYZ to CAT02 transformed sample values rgb */                                        /* eqn (6), (7) */
    rgb[0] =  0.7328 * XYZ[0] + 0.4296 * XYZ[1] - 0.1624 * XYZ[2];
    rgb[1] = -0.7036 * XYZ[0] + 1.6975 * XYZ[1] + 0.0061 * XYZ[2];
    rgb[2] =  0.0030 * XYZ[0] + 0.0136 * XYZ[1] + 0.9834 * XYZ[2];

    
    /* Chromaticaly transformed sample values rgbc */
    rgbc[0] = s->Drgb[0] * rgb[0];
    rgbc[1] = s->Drgb[1] * rgb[1];
    rgbc[2] = s->Drgb[2] * rgb[2];

    
    /* CAT02 to HPE transformed sample values rgbp */                                       /* eqn (10), (11), (12) */
    rgbp[0] =  0.7409792 * rgbc[0] + 0.2180250 * rgbc[1] + 0.0410058 * rgbc[2];
    rgbp[1] =  0.2853532 * rgbc[0] + 0.6242014 * rgbc[1] + 0.0904454 * rgbc[2];
    rgbp[2] = -0.0096280 * rgbc[0] - 0.0056980 * rgbc[1] + 1.0153260 * rgbc[2];


    /* Non-linear adapted sample values rgba */
    for (i = 0; i < 3; i++) {
        if (rgbp[i] < 0) {
            Jab[0] = Jab[1] = Jab[2] = -1.0;
            s->color_a   = -1;
            s->color_b   = -1;
            s->color_rS  =  0;
            s->color_h   =  0;
            s->color_e   =  0;
            s->color_ttd =  0;
            s->color_t   =  0;
            s->color_H   =  0;
            s->color_A   =  0;
            s->color_J   = -1;
            s->color_JJ  =  0;
            s->color_Q   =  0;
            s->color_C   =  0;
            s->color_M   =  0;
            s->color_s   =  0;
            return 1;
        }
        else {
            tt = pow(s->Fl * rgbp[i], 0.42);
            rgba[i] = 400.0 * tt / (tt + 27.13) + 0.1;                                      /* eqn (13) */
        }
    }


    /* Red-green & yellow-blue opponent dimensions a, b */
    a  = rgba[0] - 12.0 * rgba[1]/11.0 + rgba[2]/11.0;                                      /* eqn (14) */
    b  = (1.0/9.0) * (rgba[0] + rgba[1] - 2.0 * rgba[2]);                                   /* eqn (15) */
    rS = sqrt(a * a + b * b);
    
    
    /* Hue angle h */
    h = (180.0/DBL_PI) * atan2(b,a);                                                        /* eqn (17) */
    h = (h < 0.0) ? h + 360.0 : h;
 
    
    /* Eccentricity factor e */
    e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));              /* eqn (18) */
    
    
    /* Intermediate magnitude t */
    ttd = rgba[0] + rgba[1] + (21.0/20.0) * rgba[2];
    t   = e * sqrt(a * a + b * b) / ttd;                                                    /* eqn (16) */
    
    
    /* Hue composition H */                                                                 /* eqn (19) */
    if (h < 20.14) {
        temp = ((h + 122.47)/1.2) + ((20.14 - h)/0.8);
        H = 300 + (100*((h + 122.47)/1.2)) / temp;
    }
    else if (h < 90.0) {
        temp = ((h - 20.14)/0.8) + ((90.00 - h)/0.7);
        H = (100*((h - 20.14)/0.8)) / temp;
    }
    else if (h < 164.25) {
        temp = ((h - 90.00)/0.7) + ((164.25 - h)/1.0);
        H = 100 + ((100*((h - 90.00)/0.7)) / temp);
    }
    else if (h < 237.53) {
        temp = ((h - 164.25)/1.0) + ((237.53 - h)/1.2);
        H = 200 + ((100*((h - 164.25)/1.0)) / temp);
    }
    else {
        temp = ((h - 237.53)/1.2) + ((360 - h + 20.14)/0.8);
        H = 300 + ((100*((h - 237.53)/1.2)) / temp);
    }
  
    
    /* Achromatic response A */
    A = (2.0 * rgba[0] + rgba[1] + (1.0/20.0) * rgba[2] - 0.305) * s->Nbb;                  /* eqn (20) */


    /* Lightness J */
    J = pow(A / s->Aw, s->C * s->z);                                                        /* eqn (21) */
    

    /* Brightness Q */
    Q = (4.0 / s->C) * sqrt(J) * (s->Aw + 4.0) * pow(s->Fl, 0.25);                          /* eqn (22) */

    
    /* Chroma C */
    C = pow(t, 0.9) * sqrt(J) * s->nn;                                                      /* eqn (23) */

    
    /* Colorfulness M */
    M = C * pow(s->Fl, 0.25);                                                               /* eqn (24) */

    
    /* Saturation s */
    sat = 100 * sqrt(M / Q);                                                                /* eqn (25) */
    
    
    /* Jab[] return values, in 0 .. 100 range */
    Jab[0] = J * 100;
    Jab[1] = a * 100;
    Jab[2] = b * 100;

    
    /* Transfer internal sample values */
    for (int i = 0; i < 3; i++) {
        s->color_rgb[i]  = rgb[i];
        s->color_rgbc[i] = rgbc[i];
        s->color_rgbp[i] = rgbp[i];
        s->color_rgba[i] = rgba[i];
    }
    s->color_a   = a;
    s->color_b   = b;
    s->color_rS  = rS;
    s->color_h   = h;
    s->color_e   = e;
    s->color_ttd = ttd;
    s->color_t   = t;
    s->color_H   = H;
    s->color_A   = A;
    s->color_J   = J;
//  s->color_JJ  = JJ;
    s->color_Q   = Q;
    s->color_C   = C;
    s->color_M   = M;
    s->color_s   = sat;
    

    if (printCAM == 1) {
        printf("Processing XYZ->Jab (CIECAM02 by the book):\n");
        printf("XYZ = %f %f %f\n", XYZ[0]*100, XYZ[1]*100, XYZ[2]*100);
        printf("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
        printf("Chromatically transformed sample value rgbc = %f %f %f\n", rgbc[0], rgbc[1], rgbc[2]);
        printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
        printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
        printf("Prelim red green a = %f, b = %f\n", a, b);
        printf("Hue angle h = %f\n", h);
        printf("Hue composition H = %f\n", H);
        printf("Eccentricity factor e = %f\n", e);
        printf("Preliminary magnitude t = %f\n", t);
        printf("Achromatic response A = %f\n", A);
        printf("alt scale factor C / rS = %f\n", C/rS);
        printf("Lightness J = %f\n", J * 100);
        printf("Chroma C = %f\n", C);
        printf("Brightness Q = %f\n",Q);
        printf("Colorfulness M = %f\n",M);
        printf("Saturation s = %f\n",sat);
        printf("\n");
    }

	return 0;
}


static int cam_to_XYZ
(   struct _ciecam02 *s,
    double XYZ[3],
    double Jab[3]       )
{
	int i;
	double xyz[3], rgbp[3], rgba[3];
	double a, b, ja, jb, J, JJ, C, rC, h, e, A, ss;
	double tt, cJ, ttA;
	double k1, k2, k3;		/* (k1 & k3 are different to the fwd k1 & k3) */

//    #ifdef DIAG2        /* Incase in == out */
//    double Jabi[3];
//
//    Jabi[0] = Jab[0];
//    Jabi[1] = Jab[1];
//    Jabi[2] = Jab[2];
//
//    #endif
//
//    TRACE(("\nCIECAM02 Reverse conversion:\n"))
//    TRACE(("Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]))
//
//    JJ = Jab[0] * 0.01;    /* J/100 */
//    ja = Jab[1];
//    jb = Jab[2];
//
//
//    /* Convert Jab to core A, S & h values: */
//
//    /* Compute hue angle */
//    h = (180.0/DBL_PI) * atan2(jb, ja);
//    h = (h < 0.0) ? h + 360.0 : h;
//
//    /* Compute chroma value */
//    C = sqrt(ja * ja + jb * jb);    /* Must be Always +ve, Can be NZ even if J == 0 */
//
//    /* Preliminary Restricted chroma, always +ve and NZ */
//    /* (The exact value isn't important because the numerator dominates as a,b aproach 0) */
//    rC = C;
//    if (rC < DBL_EPSILON)
//        rC = DBL_EPSILON;
//
//    J = JJ;
//
//#ifndef DISABLE_HHKR
//    #ifndef CLASSIC
//        /* Undo Helmholtz-Kohlrausch effect */
//        if (s->hk && J < 1.0) {
//    //        double kk = C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
//            double kk = s->hkscale * HHKR_MUL * C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
//            if (kk > 1e-6)     /* Limit kk to a reasonable range */
//                kk = 1.0/(s->hklimit + 1.0/kk);
//            J = (JJ - kk)/(1.0 - kk);
//            if (J < 0.0)
//                J = JJ - kk;
//            TRACE(("J = %f from JJ = %f, kk = %f\n",J,JJ,kk))
//        }
//    #endif
//#endif /* DISABLE_HHKR */
//
//    /* Achromatic response */
//#ifndef SYMETRICJ        /* Cut to a straight line */
//    if (J >= s->jlimit) {
//        A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
//    } else {    /* In the straight line segment */
//        A = s->lA/s->jlimit * J;
//        TRACE(("Undo Acromatic straight line\n"))
//    }
//#else            /* Symetric */
//    if (J >= 0.0) {
//        A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
//    } else {    /* In the straight line segment */
//        A = -pow(-J, 1.0/(s->C * s->z)) * s->Aw;
//        TRACE(("Undo symetric Acromatic\n"))
//    }
//#endif
//
//    /* Preliminary Acromatic response +ve */
//    ttA = (A/s->Nbb)+0.305;
//
//    if (A > 0.0) {
//        cJ = pow(A/s->Aw, s->C * s->z);
//        if (cJ < s->ssmincj)
//            cJ = s->ssmincj;
//    } else {
//        cJ = s->ssmincj;
//    }
//    TRACE(("C = %f, A = %f from J = %f, cJ = %f\n",C, A,J,cJ))
//
//    /* Eccentricity factor */
//    e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));
//
//    /* ab scale components */
//    k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rC, 1.0/9.0);
//    k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;
//    k3 = s->dcomp[1] * ja + s->dcomp[2] * jb;
//
//    TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, (k1 - k3)/k2))
//
//#ifdef ENABLE_DDL
//
//    /* Limit ratio of k3 to k1 to stop zero or -ve ss */
//    if (k3 > (k1 * s->ddulimit)) {
//        k3 = k1 * s->ddulimit;
//        TRACE(("k3 limited to %f due to k3:k1 ratio, ss = %f\n",k3,(k1 - k3)/k2))
//    }
//
//    /* See if there is going to be a problem in fwd */
//    if (k3 < -k1 * s->ddllimit/(1.0 - s->ddllimit)) {
//        /* Adjust ss to allow for fwd limitd computation */
//        k3 = -k1 * s->ddllimit/(1.0 - s->ddllimit);
//        TRACE(("k3 set to %f to allow for fk3:fk2 fwd limit\n",k3))
//    }
//
//#endif /* ENABLE_DDL */
//
//#ifdef DISABLE_TTD
//
//    ss = k1/k2;
//
//#else /* !DISABLE_TTD */
//
//    /* Compute the ab scale factor */
//    ss = (k1 - k3)/k2;
//
//    ss = 100;
//
//
//#endif /* !DISABLE_TTD */
//
//#ifdef ENABLE_SS
//    if (ss < SSLLIMIT)
//        ss = SSLLIMIT;
//    else if (ss > SSULIMIT)
//        ss = SSULIMIT;
//#endif /* ENABLE_SS */
//
//    /* Unscale a and b */
//    a = ja / ss;
//    b = jb / ss;
//
//    TRACE(("ss = %f, ttA = %f, a = %f, b = %f\n",ss,ttA,a,b))
//
//#ifdef CLASSIC
//    a = ja / 100;
//    b = jb / 100;
//#endif
//
//    /* Solve for post-adapted cone response of sample */
//    /* using inverse matrix on ttA, a, b */
//    rgba[0] = (20.0/61.0) * ttA
//            + ((41.0 * 11.0)/(61.0 * 23.0)) * a
//            + ((288.0 * 1.0)/(61.0 * 23.0)) * b;
//    rgba[1] = (20.0/61.0) * ttA
//            - ((81.0 * 11.0)/(61.0 * 23.0)) * a
//            - ((261.0 * 1.0)/(61.0 * 23.0)) * b;
//    rgba[2] = (20.0/61.0) * ttA
//            - ((20.0 * 11.0)/(61.0 * 23.0)) * a
//            - ((20.0 * 315.0)/(61.0 * 23.0)) * b;
//
//    TRACE(("rgba %f %f %f\n",rgba[0], rgba[1], rgba[2]))
//
//#ifdef DISABLE_MATRIX
//
//    XYZ[0] = rgba[0];
//    XYZ[1] = rgba[1];
//    XYZ[2] = rgba[2];
//
//#else /* !DISABLE_MATRIX */
//
//#ifdef DISABLE_NONLIN
//    rgbp[0] = 27.13/400.0 * rgba[0];
//    rgbp[1] = 27.13/400.0 * rgba[1];
//    rgbp[2] = 27.13/400.0 * rgba[2];
//#else    /* !DISABLE_NONLIN */
//
///* hier komen de verschillen vandaan   alle drie door vergelijking met 0.42 */
//
//        /* Hunt-Pointer_Estevez cone space
//         * (with linear segment at the +ve end) */
//        for (i = 0; i < 3; i++) {
//            if (rgba[i] < s->nldxval) {
//                rgbp[i] = s->nldlimit + (rgba[i] - s->nldxval)/s->nldxslope;
//            } else if (rgba[i] <= s->nluxval) {
//                tt = rgba[i] - 0.1;
//                rgbp[i] = pow((27.13 * tt)/(400.0 - tt), 1.0/0.42)/s->Fl;
//    //            rgbp[i] = pow((27.13 * tt)/(400.0 - tt), 1.0/1.0)/s->Fl;
//            } else {
//                rgbp[i] = s->nlulimit + (rgba[i] - s->nluxval)/s->nluxslope;
//            }
//        }
//
//#endif /* !DISABLE_NONLIN */
//
//    TRACE(("rgbp %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
//
//#ifdef ENABLE_BLUE_ANGLE_FIX
//    ss = rgbp[0] + rgbp[1] + rgbp[2];
//    if (ss < 1e-9)
//        ss = 0.0;
//    else {
//        ss = (rgbp[2]/ss - 1.0/3.0) * 3.0/2.0;
//        if (ss > 0.0)
//            ss = BLUE_BL_MAX * pow(ss, BLUE_BL_POW);
//    }
//    if (ss < 0.0)
//        ss = 0.0;
//    else if (ss > 1.0)
//        ss = 1.0;
//    tt = 0.5 * (rgbp[0] + rgbp[1]);
//    rgbp[0] = (rgbp[0] - ss * tt)/(1.0 - ss);
//    rgbp[1] = (rgbp[1] - ss * tt)/(1.0 - ss);
//
//    TRACE(("rgbp after blue fix %f = %f %f %f\n",ss,rgbp[0], rgbp[1], rgbp[2]))
//#endif
//
//
//#ifdef ENABLE_DECOMPR
//    /* Undo soft limiting */
//    {
//        double tt;            /* Temporary */
//        double wrgb[3];        /* White target */
//
//        /* Make white target white point with same Y value */
//        tt = rgbp[0] * s->icc[1][0] + rgbp[1] * s->icc[1][1] + rgbp[2] * s->icc[1][2];
//        tt = tt > BC_WHMINY ? tt : BC_WHMINY;    /* Limit to minimum Y */
//        icmScale3(wrgb, s->rgbpW, tt/s->Wxyz[1]);    /* White target at same Y */
//        TRACE(("wrgb %f %f %f\n", wrgb[0], wrgb[1], wrgb[2]))
//
//        /* Un-limit b,g,r in turn */
//        for (i = 2; i >= 0; i--) {             double cv;        /* Compression value */
//            double ctv;        /* Compression target value */
//            double cd;        /* Compression displacement needed */
//            double cvec[3];    /* Normalized correction vector */
//            double isec[3];    /* Intersection with plane */
//            double offs;    /* Offset of point from orgin*/
//            double range;    /* Threshold to start compression */
//            double asym;    /* Compression asymtope */
//
//            /* Compute compression direction as vector towards white */
//            /* (We did try correcting in a blend of limit plane normal and white, */
//            /*  but compressing towards white seems to be the best.) */
//            icmSub3(cvec, wrgb, rgbp);                    /* Direction of white target */
//
//            TRACE(("ch %d, rgbp %f %f %f\n", i, rgbp[0], rgbp[1], rgbp[2]))
//            TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))
//
//            if (cvec[i] < 1e-9) {        /* compression direction can't correct this coord */
//                TRACE(("Intersection with limit plane failed\n"))
//                continue;
//            }
//
//            /* Scale compression vector to make it move a unit in normal direction */
//            icmScale3(cvec, cvec, 1.0/cvec[i]);        /* Normalized vector to white */
//            TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))
//
//            /* Compute intersection of correction direction with this limit plane */
//            /* (This corresponds with finding displacement of rgbp by cvec */
//            /*  such that the current coord value = 0) */
//            icmScale3(isec, cvec, -rgbp[i]);        /* (since cvec[i] == 1.0) */
//            icmAdd3(isec, isec, rgbp);
//            TRACE(("isec %f %f %f\n", isec[0], isec[1], isec[2]))
//
//            /* Compute distance from intersection to origin */
//            offs = pow(icmNorm3(isec), 0.85);
//
//            range = s->crange[i] * offs;    /* Scale range by distance to origin */
//            if (range > BC_MAXRANGE)        /* so that it tapers down as we approach it */
//                range = BC_MAXRANGE;        /* and limit maximum */
//
//            /* Aiming above plane when far from origin, */
//            /* but below plane at the origin, so that black isn't affected. */
//            asym = range - 0.2 * (range + (0.01 * s->crange[i]));
//
//            ctv = cv = rgbp[i];        /* Distance above/below limit plane */
//
//            TRACE(("ch %d, offs %f, range %f asym %f, cv %f\n",i, offs,range,asym,cv))
//
//            if (ctv < (range - 1e-12)) {        /* Need to expand */
//
//                if (ctv <= asym) {
//                    cd = BC_LIMIT;
//                    TRACE(("ctv %f < asym %f\n",ctv,asym))
//                } else {
//                    double aa, bb;
//                    aa = 1.0/(range - ctv);
//                    bb = 1.0/(range - asym);
//                    if (aa > (bb + 1e-12))
//                        cv = range - 1.0/(aa - bb);
//                    cd = ctv - cv;                /* Displacement needed */
//                }
//                if (cd > BC_LIMIT)
//                    cd = BC_LIMIT;
//                TRACE(("ch %d cd = %f, scaled cd %f\n",i,cd,cd))
//
//                if (cd > 1e-9) {
//                    icmScale3(cvec, cvec, -cd);            /* Compression vector */
//                    icmAdd3(rgbp, rgbp, cvec);            /* Compress by displacement */
//                    TRACE(("rgbp after decomp. = %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
//                }
//            }
//        }
//    }
//#endif /* ENABLE_COMPR */
//
//    /* Chromaticaly transformed sample value */
//    /* Spectrally sharpened cone responses */
//    /* XYZ values */
//    icmMulBy3x3(xyz, s->icc, rgbp);
//
//    TRACE(("XYZ = %f %f %f\n",xyz[0], xyz[1], xyz[2]))
//
//    /* Subtract flare
//    XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
//    XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
//    XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);
//    */
//    XYZ[0] = xyz[0];
//    XYZ[1] = xyz[1];
//    XYZ[2] = xyz[2];
//
//#endif /* !DISABLE_MATRIX */
//
//    TRACE(("XYZ after flare = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))
//    TRACE(("\n"))
//
//#ifdef DIAG2
//    printf("Processing:\n");
//    printf("Jab = %f %f %f\n", Jabi[0], Jabi[1], Jabi[2]);
//    printf("Chroma C = %f\n", C);
//    printf("Preliminary Saturation ss = %f\n", ss);
//    printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ * 100);
//    printf("Achromatic response A = %f\n", A);
//    printf("Eccentricity factor e = %f\n", e);
//    printf("Hue angle h = %f\n", h);
//    printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
//    printf("Hunundeft-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
//    printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
//    printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
//    printf("\n");
//#endif

    return 0;
}

