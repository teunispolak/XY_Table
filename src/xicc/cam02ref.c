//
//  cam02ref.c
//  hargyllcms-1.9.2
//
//  Created by Teunis on 09/11/2019.
//  Copyright © 2019 Teunis Polak. All rights reserved.
//

/*
 * cam02, unoptimised, untweaked reference version for testing.
 * with optional trace/range error flags.
 *
 * Color Appearance Model.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2004
 * Version: 1.00
 *
 * This file is based on cam97s3.c by Graeme Gill.
 *
 * Copyright 2004, 2007 Graeme W. Gill
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
#include "icc.h"
#include "xcam.h"
#include "numlib.h"
#include "cam02ref.h"
#include "cam02.h"

#define HHKR_MUL 0.25               /* [0.25] - Helmholtz-Kohlrausch strength multiplier */

#undef DIAG                         /* Print internal value diagnostics for each conversion */




static void cam02ref_free(cam02ref *s);
static int cam02ref_set_view(cam02ref *s, ViewingCondition Ev, double Wxyz[3],
                             double Yb, double La, double Lv, double Yf, double Yg, double Gxyz[3],
                             int hk, double hkscale); // double mtaf, double Wxyz2[3]);
static int cam02ref_XYZ_to_cam(cam02ref *s, double *Jab, double *xyz);
static int cam02ref_cam_to_XYZ(cam02ref *s, double XYZ[3], double Jab[3]);

/* Return a viewing condition enumeration from the given Ambient and
 * Adapting/Surround Luminance. */
static ViewingCondition cam02ref_Ambient2VC(
double La,              /* Ambient Luminence (cd/m^2) */
double Lv               /* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
) {
    double r;
    if (fabs(La) < 1e-10)         /* Hmm. */
        r = 1.0;
    else
        r = La/Lv;
    if (r < 0.01)
        return vc_dark;
    if (r < 0.2)
        return vc_dim;
    return vc_average;
}


/* Create a cam02 conversion object */
cam02ref *new_cam02ref(void) {
    cam02ref *s;
    
    if ((s = (cam02ref *)malloc(sizeof(cam02ref))) == NULL) {
        fprintf(stderr,"cam02: malloc failed allocating object\n");
        exit(-1);
    }
    
    /* Initialise methods */
    s->del          = cam02ref_free;
    s->set_view     = cam02ref_set_view;
    s->XYZ_to_cam   = cam02ref_XYZ_to_cam;
    s->cam_to_XYZ   = cam02ref_cam_to_XYZ;
    
    return s;
}

static void cam02ref_free(cam02ref *s) {
    if (s != NULL)
        free(s);
}

static int cam02ref_set_view(
cam02ref *s,
ViewingCondition Ev,    /* Enumerated Viewing Condition */
double Wxyz[3],         /* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
double La,              /* Adapting/Surround Luminance cd/m^2 */
double Yb,              /* Relative Luminence of Background to reference white */
double Lv,              /* Luminence of white in the Viewing/Scene/Image field (cd/m^2) */
                        /* Ignored if Ev is set to other than vc_none */
double Yf,              /* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
double Yg,              /* Glare as a fraction of the adapting/surround (Y range 0.0 .. 1.0) */
double Gxyz[3],         /* The Glare white coordinates (typically the Ambient color) */
int hk,                 /* Flag, NZ to use Helmholtz-Kohlraush effect */
double hkscale          /* HK effect scaling factor */
) {
//double mtaf,          /* Mid tone partial adapation factor from Wxyz to Wxyz2, <= 0.0 if none */
//double Wxyz2[3]       /* Mid tone Adapted White XYZ (Y range 0.0 .. 1.0) */
//) {
    double tt;
    
    if (Ev == vc_none)        /* Compute enumerated viewing condition */
        Ev = cam02ref_Ambient2VC(La, Lv);
    /* Transfer parameters to the object */
    s->Ev = Ev;
    s->Wxyz[0] = Wxyz[0];
    s->Wxyz[1] = Wxyz[1];
    s->Wxyz[2] = Wxyz[2];
    s->Yb = Yb > 0.005 ? Yb : 0.005;    /* Set minimum to avoid divide by 0.0 */
    s->La = La;
    s->Yf = Yf;
    if (Gxyz[0] > 0.0 && Gxyz[1] > 0.0 && Gxyz[2] > 0.0) {
        tt = Wxyz[1]/Gxyz[1];        /* Scale to white ref white */
        s->Gxyz[0] = tt * Gxyz[0];
        s->Gxyz[1] = tt * Gxyz[1];
        s->Gxyz[2] = tt * Gxyz[2];
    } else {
        s->Gxyz[0] = Wxyz[0];
        s->Gxyz[1] = Wxyz[1];
        s->Gxyz[2] = Wxyz[2];
    }
    s->hk = hk;
    s->hkscale = hkscale;
    
    /* Compute the internal parameters by category */
    switch(s->Ev) {
        case vc_dark:
            s->C = 0.525;
            s->Nc = 0.8;
            s->F = 0.8;
            Lv = La/0.033;
            break;
        case vc_dim:
            s->C = 0.59;
            s->Nc = 0.95;
            s->F = 0.9;
            Lv = La/0.1;
            break;
        case vc_average:
        default:    /* average */
            s->C = 0.69;
            s->Nc = 1.0;
            s->F = 1.0;
            Lv = La/0.2;
            break;
        case vc_cut_sheet:
            s->C = 0.41;
            s->Nc = 0.8;
            s->F = 0.8;
            Lv = La/0.02;   // ???
            break;
    }
    s->Lv = Lv;
    
    /* Compute values that only change with viewing parameters */
    
    /* Figure out the Flare contribution to the flareless XYZ input */
    s->Fsxyz[0] = s->Yf * s->Wxyz[0];
    s->Fsxyz[1] = s->Yf * s->Wxyz[1];
    s->Fsxyz[2] = s->Yf * s->Wxyz[2];
    
    /* Add in the Glare contribution from the adapting/surround */
    tt = s->Yg * s->La/s->Lv;
    s->Fsxyz[0] += tt * s->Gxyz[0];
    s->Fsxyz[1] += tt * s->Gxyz[1];
    s->Fsxyz[2] += tt * s->Gxyz[2];
    
    
    /* Next lines copied from cam02.c - seemed missing here, causing issues in forward conversion
     * and probably backward as well, since s->Fsc and s->Fisc remained undefined */
    
    /* Rescale so that the sum of the flare and the input doesn't exceed white */
    s->Fsc = s->Wxyz[1]/(s->Fsxyz[1] + s->Wxyz[1]);                                     // !!!!!!! potential divide by 0
    s->Fsxyz[0] *= s->Fsc;
    s->Fsxyz[1] *= s->Fsc;
    s->Fsxyz[2] *= s->Fsc;
    s->Fisc = 1.0/s->Fsc;
    // end copy
    
    /* Sharpened cone response white values */
    s->rgbW[0] =  0.7328 * s->Wxyz[0] + 0.4296 * s->Wxyz[1] - 0.1624 * s->Wxyz[2];
    s->rgbW[1] = -0.7036 * s->Wxyz[0] + 1.6975 * s->Wxyz[1] + 0.0061 * s->Wxyz[2];
//    s->rgbW[2] =  0.0000 * s->Wxyz[0] + 0.0000 * s->Wxyz[1] + 1.0000 * s->Wxyz[2];
    
    s->rgbW[2] =  0.0030 * s->Wxyz[0] + 0.0136 * s->Wxyz[1] + 0.9834 * s->Wxyz[2];
    
    /* Degree of chromatic adaption */
    s->D = s->F * (1.0 - exp((-s->La - 42.0)/92.0)/3.6);
    
    /* Precompute Chromatic transform values */
    s->Drgb[0] = s->D * (s->Wxyz[1]/s->rgbW[0]) + 1.0 - s->D;
    s->Drgb[1] = s->D * (s->Wxyz[1]/s->rgbW[1]) + 1.0 - s->D;
    s->Drgb[2] = s->D * (s->Wxyz[1]/s->rgbW[2]) + 1.0 - s->D;
    
    /* Chromaticaly transformed white value */
    s->rgbcW[0] = s->Drgb[0] * s->rgbW[0];
    s->rgbcW[1] = s->Drgb[1] * s->rgbW[1];
    s->rgbcW[2] = s->Drgb[2] * s->rgbW[2];
    
    /* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
//    s->rgbpW[0] =  0.7409744840453773 * s->rgbcW[0]
//                +  0.2180245944753982 * s->rgbcW[1]
//                +  0.0410009214792244 * s->rgbcW[2];
//    s->rgbpW[1] =  0.2853532916858801 * s->rgbcW[0]
//                +  0.6242015741188157 * s->rgbcW[1]
//                +  0.0904451341953042 * s->rgbcW[2];
//    s->rgbpW[2] = -0.0096276087384294 * s->rgbcW[0]
//                -  0.0056980312161134 * s->rgbcW[1]
//                +  1.0153256399545427 * s->rgbcW[2];
    
    s->rgbpW[0] =  0.7409792 * s->rgbcW[0]
                +  0.2180250 * s->rgbcW[1]
                +  0.0410058 * s->rgbcW[2];
    s->rgbpW[1] =  0.2853532 * s->rgbcW[0]
                +  0.6242014 * s->rgbcW[1]
                +  0.0904454 * s->rgbcW[2];
    s->rgbpW[2] = -0.0096280 * s->rgbcW[0]
                -  0.0056980 * s->rgbcW[1]
                +  1.0153260 * s->rgbcW[2];
    
    /* Background induction factor */
    s->n = s->Yb/ s->Wxyz[1];
    s->nn = pow(1.64 - pow(0.29, s->n), 0.73);    /* Pre computed value */
    
    /* Lightness contrast factor ?? */
    {
        double k;
        
        k = 1.0 / (5.0 * s->La + 1.0);
        s->Fl = 0.2 * pow(k , 4.0) * 5.0 * s->La
        + 0.1 * pow(1.0 - pow(k , 4.0) , 2.0) * pow(5.0 * s->La , 1.0/3.0);
    }
    
    /* Background and Chromatic brightness induction factors */
    s->Nbb   = 0.725 * pow(1.0/s->n, 0.2);
    s->Ncb   = s->Nbb;
    
    /* Base exponential nonlinearity */
    s->z = 1.48 + pow(s->n , 0.5);
    
    /* Post-adapted cone response of white */
    tt = pow(s->Fl * s->rgbpW[0], 0.42);
    s->rgbaW[0] = (400.1 * tt + 2.713) / (tt + 27.13);
    tt = pow(s->Fl * s->rgbpW[1], 0.42);
    s->rgbaW[1] = (400.1 * tt + 2.713) / (tt + 27.13);
    tt = pow(s->Fl * s->rgbpW[2], 0.42);
    s->rgbaW[2] = (400.1 * tt + 2.713) / (tt + 27.13);
    
    /* Achromatic response of white */
    s->Aw = (2.0 * s->rgbaW[0] + s->rgbaW[1] + (1.0/20.0) * s->rgbaW[2] - 0.305) * s->Nbb;

    if (printCAM == 1)
    {
        //#ifdef DIAG
        printf("Ref. Scene parameters:\n");
        printf("Viewing condition Ev = %d\n",s->Ev);
        printf("Ref white Wxyz = %f %f %f\n", s->Wxyz[0], s->Wxyz[1], s->Wxyz[2]);
        printf("Relative liminance of background Yb = %f\n", s->Yb);
        printf("Adapting liminance La = %f\n", s->La);
        printf("Flare Yf = %f\n", s->Yf);
        printf("Flare color Fxyz = %f %f %f\n", s->Fsxyz[0], s->Fsxyz[1], s->Fsxyz[2]);
    
        printf("Internal parameters:\n");
        printf("Surround Impact C = %f\n", s->C);
        printf("Chromatic Induction Nc = %f\n", s->Nc);
        printf("Adaption Degree F = %f\n", s->F);
    
        printf("Pre-computed values\n");
        printf("Sharpened cone white rgbW = %f %f %f\n", s->rgbW[0], s->rgbW[1], s->rgbW[2]);
        printf("Degree of chromatic adaption D = %f\n", s->D);
        printf("Chromatic transform values Drgb = %f %f %f\n", s->Drgb[0], s->Drgb[1], s->Drgb[2]);
        printf("Chromatically transformed white rgbcW = %f %f %f\n", s->rgbcW[0], s->rgbcW[1], s->rgbcW[2]);
        printf("Hunter-P-E cone response white rgbpW = %f %f %f\n", s->rgbpW[0], s->rgbpW[1], s->rgbpW[2]);
        printf("Background induction factor n = %f\n", s->n);
        printf("Lightness contrast factor Fl = %f\n", s->Fl);
        printf("Background brightness induction factor Nbb = %f\n", s->Nbb);
        printf("Chromatic brightness induction factor Ncb = %f\n", s->Ncb);
        printf("Base exponential nonlinearity z = %f\n", s->z);
        printf("Post adapted cone response white rgbaW = %f %f %f\n", s->rgbaW[0], s->rgbaW[1], s->rgbaW[2]);
        printf("Achromatic response of white Aw = %f\n", s->Aw);
        printf("\n");
        //#endif
    }
    return 0;
}

/* A version of the pow() function that preserves the
 * sign of its first argument. */
static double spow(double x, double y) {
    return x < 0.0 ? -pow(-x,y) : pow(x,y);
}


#define REFTRACE(xxxx) if (s->trace) printf xxxx ;

/* Conversion. Returns NZ and -1, -1, -1 if input is out of range */
static int cam02ref_XYZ_to_cam(
                               cam02ref *s,
                               double Jab[3],
                               double XYZ[3])
{
    int i;
    double xyz[3], rgb[3], rgbp[3], rgba[3], rgbaW[3], rgbc[3], rgbcW[3];
    double a, b, rS, J, C, h, H, e, A, JJ, M, Q, sat;
    double ttd, tt, t, temp;
    
    REFTRACE(("\nReference forward conversion:\n"))
    REFTRACE(("XYZ %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))
    
    /* Add in flare */
    xyz[0] = s->Fsc * XYZ[0] + s->Fsxyz[0];
    xyz[1] = s->Fsc * XYZ[1] + s->Fsxyz[1];
    xyz[2] = s->Fsc * XYZ[2] + s->Fsxyz[2];
    
    REFTRACE(("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]))
    
    /* Spectrally sharpened cone responses */
    rgb[0] =  0.7328 * xyz[0] + 0.4296 * xyz[1] - 0.1624 * xyz[2];
    rgb[1] = -0.7036 * xyz[0] + 1.6975 * xyz[1] + 0.0061 * xyz[2];
//    rgb[2] =  0.0000 * xyz[0] + 0.0000 * xyz[1] + 1.0000 * xyz[2];
    
    rgb[2] =  0.0030 * XYZ[0] + 0.0136 * XYZ[1] + 0.9834 * XYZ[2];
    
    REFTRACE(("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]))
    
    /* Chromaticaly transformed sample value */
    rgbc[0] = s->Drgb[0] * rgb[0];
    rgbc[1] = s->Drgb[1] * rgb[1];
    rgbc[2] = s->Drgb[2] * rgb[2];
    
    REFTRACE(("Chromatically transformed sample value rgbc = %f %f %f\n", rgb[0], rgb[1], rgb[2]))
    
    /* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
//    rgbp[0] =  0.7409744840453773 * rgbc[0]
//            +  0.2180245944753982 * rgbc[1]
//            +  0.0410009214792244 * rgbc[2];
//    rgbp[1] =  0.2853532916858801 * rgbc[0]
//            +  0.6242015741188157 * rgbc[1]
//            +  0.0904451341953042 * rgbc[2];
//    rgbp[2] = -0.0096276087384294 * rgbc[0]
//            -  0.0056980312161134 * rgbc[1]
//            +  1.0153256399545427 * rgbc[2];
    
    rgbp[0] =  0.7409792 * rgbc[0] + 0.2180250 * rgbc[1] + 0.0410058 * rgbc[2];
    rgbp[1] =  0.2853532 * rgbc[0] + 0.6242014 * rgbc[1] + 0.0904454 * rgbc[2];
    rgbp[2] = -0.0096280 * rgbc[0] - 0.0056980 * rgbc[1] + 1.0153260 * rgbc[2];
    
    REFTRACE(("rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]))
    
    /* Post-adapted cone response of sample. */
    /* rgba[] has a minimum value of 0.1 for XYZ[] = 0 and no flare. */
    /* We add a symetric negative compression region */
    for (i = 0; i < 3; i++) {
        if (s->range && (rgbp[i] < 0.0 || rgbp[i] < s->nldlimit)) {
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
        if (rgbp[i] < 0.0) {
            tt = pow(s->Fl * -rgbp[i], 0.42);
            rgba[i] = (2.713 - 397.387 * tt) / (tt + 27.13);
            
        } else {
            tt = pow(s->Fl * rgbp[i], 0.42);
            rgba[i] = (400.1 * tt + 2.713) / (tt + 27.13);
        }
    }
    
    REFTRACE(("rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]))
    
    /* Preliminary red-green & yellow-blue opponent dimensions */
    a     = rgba[0] - 12.0 * rgba[1]/11.0 + rgba[2]/11.0;
    b     = (1.0/9.0) * (rgba[0] + rgba[1] - 2.0 * rgba[2]);
    rS    = sqrt(a * a + b * b);        /* Normalised a, b */
    
    /* Preliminary Saturation */
    /* Note that the minimum values for rgba[] for XYZ = 0 is 0.1 */
    /* Hence magic 0.305 below comes from the following weighting of rgba[] */
    ttd = rgba[0] + rgba[1] + (21.0/20.0) * rgba[2];
    
    /* Achromatic response */
    /* Note that the minimum values of rgba[] for XYZ = 0 is 0.1, */
    /* hence magic 0.305 below comes from the following weighting of rgba[], */
    /* to base A at 0.0 */
    A = (2.0 * rgba[0] + rgba[1] + (1.0/20.0) * rgba[2] - 0.305) * s->Nbb;
    
    REFTRACE(("a = %f, b = %f, ttd = %f, rS = %f, A = %f\n", a, b, ttd, rS, A))
    
    /* Lightness */
    J = pow(A/s->Aw, s->C * s->z);        /* J/100  - keep Sign */
    
    /* Hue angle */
    h = (180.0/DBL_PI) * atan2(b,a);
    h = (h < 0.0) ? h + 360.0 : h;
    
    /* Eccentricity factor */
    e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));
    
//  if (s->range && (J < DBL_EPSILON || J < s->jlimit || ttd < DBL_EPSILON)) {
//      REFTRACE(("J = %f, ttd = %f, exit with error\n", J, ttd))
//      Jab[0] = Jab[1] = Jab[2] = -1.0;
//      return 1;
//  }
    
    /* Hue composition H */
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
    
    t = e * rS / ttd;
    
    /* Chroma */
    C = pow(t, 0.9) * sqrt(J) * s->nn;
    
    /* Colorfulness */
    M = C * pow(s->Fl, 0.25);
    
    REFTRACE(("C = %f\n", C))
    
    /* Helmholtz-Kohlraush effect */
    JJ = J;
    if (s->hk && J < 1.0) {
        double kk = s->hkscale * HHKR_MUL * C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
        if (kk > 0.9)        /* Limit kk to a reasonable range */
            kk = 0.9;
        JJ = J + (1.0 - J) * kk;
        REFTRACE(("JJ = %f from J = %f, kk = %f\n",JJ,J,kk))
        J = JJ;
    }
    
    /* Brightness Q */
    Q = (4.0 / s->C) * sqrt(J) * (s->Aw + 4.0) * pow(s->Fl, 0.25);
    
    /* Saturation s */
    sat = 100 * sqrt(M / Q);
    
    //J *= 100.0;        /* Scale J */
    
    /* Jab[] return values, in 0 .. 100 range */
    Jab[0] = J * 100;
    if (rS >= DBL_EPSILON) {
        Jab[1] = C * a/rS;
        Jab[2] = C * b/rS;
    } else {
        Jab[1] = 0.0;
        Jab[2] = 0.0;
    }
    REFTRACE(("Returning Jab %f %f %f\n", Jab[0],Jab[1],Jab[2]))
    
    
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
    s->color_JJ  = JJ;
    s->color_Q   = Q;
    s->color_C   = C;
    s->color_M   = M;
    s->color_s   = sat;


    
    if (printCAM == 1) {
        //#ifdef DIAG
            printf("Processing XYZ->Jab (CIECAM02, reference):\n");
            printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
            printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
            printf("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
            printf("Chromatically transformed sample value rgbc = %f %f %f\n", rgbc[0], rgbc[1], rgbc[2]);
            printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
            printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
            printf("Prelim red green a = %f, b = %f\n", a, b);
            printf("Hue angle h = %f\n", h);
            printf("Eccentricity factor e = %f\n", e);
            printf("Preliminary magnitude t = %f\n", t);
            printf("Achromatic response A = %f\n", A);
            printf("Lightness J = %f\n", J);
            printf("HK Lightness JJ = %f\n", JJ * 100);
            printf("Prelim Saturation rS = %f\n", rS);
            printf("alt scale factor C / rS = %f\n", C/rS);
            printf("Chroma C = %f\n", C);
            printf("Jab = %f %f %f\n", Jab[0], Jab[1], Jab[2]);
            printf("\n");
        //#endif
    }
    return 0;
}

static int cam02ref_cam_to_XYZ(
                               cam02ref *s,
                               double XYZ[3],
                               double Jab[3]
                               ) {
    int i;
    double xyz[3], rgb[3], rgbp[3], rgba[3], rgbaW[3], rgbc[3], rgbcW[3];
    double ja, jb, aa, ab, a, b, J, C, h, e, A, ss;
    double tt, ttA, tte;
    
    J = Jab[0] * 0.01;    /* J/100 */
    ja = Jab[1];
    jb = Jab[2];
    
    /* Compute hue angle */
    h = (180.0/DBL_PI) * atan2(jb, ja);
    h = (h < 0.0) ? h + 360.0 : h;
    
    /* Compute chroma value */
    C = sqrt(ja * ja + jb * jb);        /* Must be Always +ve */
    
    /* Helmholtz-Kohlraush effect */
    if (s->hk && J < 1.0) {
        double kk = s->hkscale * HHKR_MUL * C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
        if (kk > 0.9)        /* Limit kk to a reasonable range */
            kk = 0.9;
        J = (J - kk)/(1.0 - kk);
    }
    
    /* Eccentricity factor */
    e = (cos(h * DBL_PI/180.0 + 2.0) + 3.8);
    
    /* Achromatic response */
    A = spow(J, 1.0/(s->C * s->z)) * s->Aw;                /* Keep sign of J */
    
    /* Preliminary Saturation - keep +ve */
    tt = fabs(J);
    ss = pow(C/(sqrt(tt) * s->nn), 1.0/0.9);    /* keep +ve */
    
    /* Compute a & b, taking care of numerical problems */
    aa = fabs(ja);
    ab = fabs(jb);
    ttA = (A/s->Nbb)+0.305;                        /* Common factor */
    tte = 12500.0/13.0 * e * s->Nc * s->Ncb;    /* Common factor */
    
    if (aa < 1e-10 && ab < 1e-10) {
        a = ja;
        b = jb;
    } else if (aa > ab) {
        double tanh = jb/ja;
        double sign = (h > 90.0 && h <= 270.0) ? -1.0 : 1.0;
        
        if (ttA < 0.0)
            sign = -sign;
        
        a = (ss * ttA)
        / (sign * sqrt(1.0 + tanh * tanh) * tte + (ss * (11.0/23.0 + (108.0/23.0) * tanh)));
        b = a * tanh;
        
    } else {    /* ab > aa */
        double itanh = ja/jb;
        double sign = (h > 180.0 && h <= 360.0) ? -1.0 : 1.0;
        
        if (ttA < 0.0)
            sign = -sign;
        
        b = (ss * ttA)
        / (sign * sqrt(1.0 + itanh * itanh) * tte + (ss * (108.0/23.0 + (11.0/23.0) * itanh)));
        a = b * itanh;
    }
    
    /* Post-adapted cone response of sample */
    rgba[0] = (20.0/61.0) * ttA
    + ((41.0 * 11.0)/(61.0 * 23.0)) * a
    + ((288.0 * 1.0)/(61.0 * 23.0)) * b;
    rgba[1] = (20.0/61.0) * ttA
    - ((81.0 * 11.0)/(61.0 * 23.0)) * a
    - ((261.0 * 1.0)/(61.0 * 23.0)) * b;
    rgba[2] = (20.0/61.0) * ttA
    - ((20.0 * 11.0)/(61.0 * 23.0)) * a
    - ((20.0 * 315.0)/(61.0 * 23.0)) * b;
    
    /* Hunt-Pointer_Estevez cone space */
    tt = 1.0/s->Fl;
    for (i = 0; i < 3; i++) {
        if (rgba[i] < 0.1) {
            double ta = rgba[i] > -396.387 ? rgba[i] : -396.387;
            rgbp[i] = -tt * pow((2.713 - 27.13 * rgba[i] )/(397.387 + ta), 1.0/0.42);
        } else {
            double ta = rgba[i] < 399.1 ? rgba[i] : 399.1;
            rgbp[i] =  tt * pow((27.13 * rgba[i] -2.713)/(400.1 - ta), 1.0/0.42);
        }
    }
    
    /* Chromaticaly transformed sample value */
    rgbc[0] =  1.5591523979049677 * rgbp[0]
    -  0.5447226796590880 * rgbp[1]
    -  0.0144453097698588 * rgbp[2];
    rgbc[1] = -0.7143267176368630 * rgbp[0]
    +  1.8503099728895096 * rgbp[1]
    -  0.1359761119854705 * rgbp[2];
    rgbc[2] =  0.0107755117023383 * rgbp[0]
    +  0.0052187662221759 * rgbp[1]
    +  0.9840056143203700 * rgbp[2];
    
    /* Spectrally sharpened cone responses */
    rgb[0]  =  rgbc[0]/s->Drgb[0];
    rgb[1]  =  rgbc[1]/s->Drgb[1];
    rgb[2]  =  rgbc[2]/s->Drgb[2];
    
    /* XYZ values */
    xyz[0] =  1.0961238208355140 * rgb[0]
    -  0.2788690002182872 * rgb[1]
    +  0.1827451793827730 * rgb[2];
    xyz[1] =  0.4543690419753590 * rgb[0]
    +  0.4735331543074120 * rgb[1]
    +  0.0720978037172291 * rgb[2];
    xyz[2] = -0.0096276087384294 * rgb[0]
    -  0.0056980312161134 * rgb[1]
    +  1.0153256399545427 * rgb[2];
    
    /* Subtract flare */
    XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
    XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
    XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);
    
#ifdef DIAG
    printf("Processing:\n");
    printf("Jab = %f %f %f\n", Jab[0], Jab[1], Jab[2]);
    printf("Chroma C = %f\n", C);
    printf("Preliminary Saturation ss = %f\n", ss);
    printf("Lightness J = %f\n", J * 100.0);
    printf("Achromatic response A = %f\n", A);
    printf("Eccentricity factor e = %f\n", e);
    printf("Hue angle h = %f\n", h);
    printf("Prelim red green a = %f, b = %f\n", a, b);
    printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
    printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
    printf("Chromatically transformed sample value rgbc = %f %f %f\n", rgbc[0], rgbc[1], rgbc[2]);
    printf("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
    printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
    printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
    printf("\n");
#endif
    return 0;
}


