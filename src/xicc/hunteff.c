//
//  hunteff.c
//  hargyllcms-1.9.2
//
//  Created by Teunis on 07/11/2019.
//  Copyright © 2019 Teunis Polak. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "icc.h"
#include "xcam.h"
#include "numlib.h"
#include "cam02.h"
#include "cam02ref.h"
#include "ciecam02.h"

#undef DIAG                     /* Print internal value diagnostics for each spot test conversion */
                                /* and print diagnostics for excessive errors, nans etc. */
#undef VERBOSE                  /* Print diagnostic values for every conversion */
#define SPOTTEST                /* ** Test known spot colors */
#define USE_HK 1                /* Use Helmholtz-Kohlraush in testing */


#define MAX_SPOT_ERR 2.0
#define MAX_REF_ERR 2.0         /* Maximum permitted error to reference transform in delta Jab */

#ifndef _isnan
#define _isnan(x) ((x) != (x))
#define _finite(x) ((x) == 0.0 || (x) * 1.0000001 != (x))
#endif



/* CIELab to XYZ */
static void
Lab2XYZ(double *in, double *out) {
    double L = in[0], a = in[1], b = in[2];
    double x,y,z,fx,fy,fz;
    
    if (L > 8.0) {
        fy = (L + 16.0)/116.0;
        y = pow(fy,3.0);
    } else {
        y = L/903.2963058;
        fy = 7.787036979 * y + 16.0/116.0;
    }
    
    fx = a/500.0 + fy;
    if (fx > 24.0/116.0)
        x = pow(fx,3.0);
    else
        x = (fx - 16.0/116.0)/7.787036979;
    
    fz = fy - b/200.0;
    if (fz > 24.0/116.0)
        z = pow(fz,3.0);
    else
        z = (fz - 16.0/116.0)/7.787036979;
    
    out[0] = x * 0.96422;
    out[1] = y * 1.00000;
    out[2] = z * 0.82521;
}

/* CIE XYZ to perceptual Lab */
static void
XYZ2Lab(double *in, double *out) {
    double X = in[0], Y = in[1], Z = in[2];
    double x,y,z,fx,fy,fz;
    
    x = X/0.96422;               /* D50 ? */
    y = Y/1.00000;
    z = Z/0.82521;
    
    if (x > 0.008856451586)
        fx = pow(x,1.0/3.0);
    else
        fx = 7.787036979 * x + 16.0/116.0;
    
    if (y > 0.008856451586)
        fy = pow(y,1.0/3.0);
    else
        fy = 7.787036979 * y + 16.0/116.0;
    
    if (z > 0.008856451586)
        fz = pow(z,1.0/3.0);
    else
        fz = 7.787036979 * z + 16.0/116.0;
    
    out[0] = 116.0 * fy - 16.0;
    out[1] = 500.0 * (fx - fy);
    out[2] = 200.0 * (fy - fz);
}

static void
sRGB2XYZ(double *in, double *out) {
    double r, g, b;
    r = in[0] / 255.0;
    g = in[1] / 255.0;
    b = in[2] / 255.0;
    
    /* sRGB curve to linear */
    r = r > 0.04045 ? pow(((r + 0.055) / 1.055), 2.4) : (r / 12.92);
    g = g > 0.04045 ? pow(((g + 0.055) / 1.055), 2.4) : (g / 12.92);
    b = b > 0.04045 ? pow(((b + 0.055) / 1.055), 2.4) : (b / 12.92);
    
    /* Convert linear rgb to XYZ in [0,100] rather than [0,1] */
    out[0] = ( (r * 0.4124) + (g * 0.3576) + (b * 0.1805) ) * 100 ;
    out[1] = ( (r * 0.2126) + (g * 0.7152) + (b * 0.0722) ) * 100 ;
    out[2] = ( (r * 0.0193) + (g * 0.1192) + (b * 0.9505) ) * 100 ;
}

static void
ProPhoto2XYZ(double *in, double *out) {
    double r, g, b;
    
    r = in[0] / 255.0;
    g = in[1] / 255.0;
    b = in[2] / 255.0;
    
    /* ProPhoto gamma curve to linear */
    r = pow (r, 1.8);
    g = pow (g, 1.8);
    b = pow (b, 1.8);

    /* Convert linear rgb to XYZ in [0,100] rather than [0,1] */
    out[0] = ( (r * 0.7976749) + (g * 0.1351917) + (b * 0.0313534) ) * 100;
    out[1] = ( (r * 0.2880402) + (g * 0.7118741) + (b * 0.0000857) ) * 100;
    out[2] = ( (r * 0.0000000) + (g * 0.0000000) + (b * 0.8252100) ) * 100;
}


static void
XYZ2ProPhoto(double *in, double *out) {

    double x = in[0];
    double y = in[1];
    double z = in[2];
    
    /* XYZ [0 .. 1] to linear ProPhoto RGB */
    double r = ( (x *  1.3459433) + (y * -0.2556075) + (z * -0.0511118) );
    double g = ( (x * -0.5445989) + (y *  1.5081673) + (z *  0.0205351) );
    double b = ( (x *  0.0000000) + (y *  0.0000000) + (z *  1.2118128) );
    
    /* Protect from negative values resulting from rounding */
    r = (r < 0) ? 0.0 : r;
    g = (g < 0) ? 0.0 : g;
    b = (b < 0) ? 0.0 : b;
    
    /* ProPhoto linear RGB to gamma RGB */
    r = pow (r, 1/1.8);
    g = pow (g, 1/1.8);
    b = pow (b, 1/1.8);
    
    /* Gamma RGB [0 .. 1] to [0 .. 255] */
    out[0] = r * 255;
    out[1] = g * 255;
    out[2] = b * 255;
}


static void
Jab2JCh(double *in, double *out) {
    double J = in[0], a = in[1], b = in[2];
    double C, h;

    /* Compute chroma value */
    C = sqrt(a * a + b * b);
    
    /* Compute hue angle */
    h = (180.0/3.14159265359) * atan2(b, a);
    h = (h < 0.0) ? h + 360.0 : h;

    out[0] = J;
    out[1] = C;
    out[2] = h;
}


/* Return maximum difference */
double maxdiff(double in1[3], double in2[3]) {
    int i;
    double tt, rv = 0.0;
    
    /* Deal with the possibility that we have nan's */
    for (i = 0; i < 3; i++) {
        tt = fabs(in1[i] - in2[i]);
        if (!_finite(tt))
            return tt;
        if (tt > rv)
            rv = tt;
    }
    return rv;
}

/* Return absolute difference */
double absdiff(double in1[3], double in2[3]) {
    double tt, rv = 0.0;
    tt = in1[0] - in2[0];
    rv += tt * tt;
    tt = in1[1] - in2[1];
    rv += tt * tt;
    tt = in1[2] - in2[2];
    rv += tt * tt;
    return sqrt(rv);
}

/* Return maximum Lab difference of XYZ */
double maxxyzdiff(double i1[3], double i2[3]) {
    int i;
    double tt, rv = 0.0;
    double in1[3], in2[3];
    
    XYZ2Lab(in1, i1);
    XYZ2Lab(in2, i2);
    
    /* Deal with the possibility that we have nan's */
    for (i = 0; i < 3; i++) {
        tt = fabs(in1[i] - in2[i]);
        if (!_finite(tt))
            return tt;
        if (tt > rv)
            rv = tt;
    }
    return rv;
}


int main(void)
{
    cam02    *cam02;
    cam02ref *cam02ref;
    ciecam02 *ciecam02;
    int i, j;

    cam02           = new_cam02();
    cam02ref        = new_cam02ref();
    ciecam02        = new_ciecam02();
    
    cam02ref->range     = 1;
    cam02ref->nldlimit  = cam02->nldlimit;
    cam02ref->jlimit    = cam02->jlimit;

    cam02->trace        = 0;
    cam02ref->trace     = 0;

    double white[3] = { 0.96422, 1.000, 0.82521 };
    double Yb       = 0.20;
    double Lv       = 0.0;
    double Yf       = 0.0;
    double Yg       = 0.0;
    double hk       = 1.0;
    double La;
    double XYZsample[3], Labsample[3];
    double XYZrev[3], RGBrev[3], Labrev[3];
    double Jab[3], Jabcam02ref[3], Jabciecam02[3], JCh[3];
    double QMh[3], JMh[3];


    /* ColorChecker 24 chart values in ProPhoto RGB
     * Source: RGB values derived from averaged spectral data (30 charts), babelcolor.com */
    double RGBsample[39][3] = {
        { 47,  47,  47},                    /* "Black" frame */
        { 82,  68,  55},                    /* Row 1 */
        {160, 136, 114},
        { 95, 102, 134},
        { 75,  86,  56},
        {119, 111, 154},
        {127, 168, 157},
        {166, 118,  53},                    /* Row 2 */
        { 79,  75, 143},
        {142,  84,  80},
        { 67,  50,  81},
        {144, 169,  74},
        {182, 152,  59},
        { 57,  50, 120},                    /* Row 3 */
        { 85, 123,  69},
        {121,  59,  46},
        {201, 190,  67},
        {143,  85, 126},
        { 78, 110, 146},
        {242, 243, 237},                    /* Row 4 */
        {189, 190, 189},
        {144, 145, 144},
        {101, 102, 102},
        { 66,  67,  67},
        { 38,  38,  38},
        {255, 255, 255},                    /* White and primaries, secondaries */
        {255,   0,   0},
        {  0, 255,   0},
        {  0,   0, 255},
        {  0, 255, 255},
        {255,   0, 255},
        {255, 255,   0},
        {  0,   0,   0},
        {  1,   1,   1},                    /* Near black */
        { 44,  38,   2},                    /* Ellipse */
        { 15,   8,   1},                    /* Door */
        { 20,  20,  20},
        { 15,  15,  15},
        { 10,  10,  10},
        {  5,   5,   5}
    };
    
    
    printf(" i      La      R    G    B      L    a    b            X          Y          Z            J          ja          jb");
    printf("          C          M           h           Q          s          a           b          rS           ss       ba_ss           A");
    printf("           JJ     Jab[0]");
    printf("\n");
    
    
    /* Loop for samples */
    for (i = 16; i < 17; i++) {
        
        ProPhoto2XYZ(RGBsample[i], XYZsample);
        XYZsample[0] = XYZsample[0] / 100;
        XYZsample[1] = XYZsample[1] / 100;
        XYZsample[2] = XYZsample[2] / 100;
        
        XYZ2Lab(XYZsample, Labsample);
        
        /* Loop for Adapting / Surround Luminance values */
        for (j = 0; j < 1; j++) {
            La = 5 * pow(10, j);
            La = 3185;
            Lv = La / 0.2;
        
            printCAM  = 0;
            theSource = 1;
            cam02     ->set_view(cam02,    vc_dark, white, La, Yb, Lv, Yf, Yg, white, USE_HK, hk);
            theSource = 0;
            printCAM  = 0;
            cam02ref  ->set_view(cam02ref, vc_average, white, La, Yb, Lv, Yf, Yg, white, USE_HK, hk);
            ciecam02  ->set_view(ciecam02,             white, La, Yb, Lv);

            
            /* Forward transformation */

            /* default Argyll CIECAM02 model */
            printCAM = 0;
            cam02->XYZ_to_cam(cam02, Jab, XYZsample);
            printCAM = 0;

            
            /* CIECAM02 reference model (Argyll) */
            if (cam02ref->XYZ_to_cam(cam02ref, Jabcam02ref, XYZsample))
                /*printf("Reference XYZ2Jab returned error\n") */;
     
            
            /* Billy Bigg's CIECAM02 implementation, (http://scanline.ca/ciecam02/) */
            if (ciecam02->XYZ_to_cam(ciecam02, Jabciecam02, XYZsample))
                /* printf("CIECAM02 XYZ -> Jab returned error\n") */;

            
            printf("%2d   %5d   ",
                   i, (int)La);
            printf("%4d %4d %4d   ",
                   (int)RGBsample[i][0],    (int)RGBsample[i][1],   (int)RGBsample[i][2]);
            printf("%4d %4d %4d   ",
                   (int)round(Labsample[0]), (int)round(Labsample[1]), (int)round(Labsample[2]));
            printf("%10.6f %10.6f %10.6f   ",
                   XYZsample[0],            XYZsample[1],           XYZsample[2]);
            printf("%10.6f %11.6f %11.6f ",
                   cam02->color_J,          cam02->color_ja,        cam02->color_jb);
            printf("%10.6f %10.6f %11.6f ",
                   cam02->color_C,          cam02->color_M,         cam02->color_h);
            printf("%11.6f %10.6f ",
                   cam02->color_Q,          cam02->color_s);
            printf("%10.6f %11.6f %11.6f   ",
                   cam02->color_a,          cam02->color_b,         cam02->color_rS);
            printf("%10.6f %11.6f %11.6f   ",
                   cam02->color_ss,         cam02->color_bass,      cam02->color_A);
            printf("%10.6f %10.6f",
                   cam02->color_JJ,         Jab[0]);

            printf("\n");
            
            /* Convert to JCh */
            //Jab2JCh(Jab, JCh);

            
            /* Inverse transformation */
////            La = 10 * La;
            La = La/5;
            La = 32;
            Lv = La * 5;
            
            printCAM = 0;
            theDest  = 1;
            cam02    ->set_view(cam02,    vc_average, white, La, Yb, Lv, Yf, Yg, white, USE_HK, hk);
            theDest  = 0;
            printCAM = 0;
            
            /* Jab -> XYZ  DIT KLOPT NIET*/
//            Jab[0] = cam02->color_J;
//            Jab[1] = cam02->color_ja;
//            Jab[2] = cam02->color_jb;

            printCAM = 0;
            cam02->trace = 0;
            cam02->cam_to_XYZ(cam02, XYZrev, Jab);
            printCAM = 0;
            cam02->trace = 0;
            
//            /* QMh -> XYZ */
//            QMh[0] = cam02->color_Q;
//            QMh[1] = cam02->color_M;
//            QMh[2] = cam02->color_h;
//
//            cam02->QMh_to_XYZ(cam02, XYZrev, QMh);
 
            
//            /* JMh -> XYZ */
//            JMh[0] = cam02->color_J * 100;
//            JMh[1] = cam02->color_M;
//            JMh[2] = cam02->color_h;
//            cam02->trace = 0;
//            cam02->JMh_to_XYZ(cam02, XYZrev, JMh);
//            cam02->trace = 0;
            

//            /* Protect from small (+ve or -ve) values resulting from rounding */
//            XYZrev[0] = (XYZrev[0] < DBL_EPSILON) ? 0.0 : XYZrev[0];
//            XYZrev[1] = (XYZrev[1] < DBL_EPSILON) ? 0.0 : XYZrev[1];
//            XYZrev[2] = (XYZrev[2] < DBL_EPSILON) ? 0.0 : XYZrev[2];
//
//
//            /* Clip when beyond White */
//            XYZrev[0] = (XYZrev[0] > white[0]) ? white[0] : XYZrev[0];
//            XYZrev[1] = (XYZrev[1] > white[1]) ? white[1] : XYZrev[1];
//            XYZrev[2] = (XYZrev[2] > white[2]) ? white[2] : XYZrev[2];
            
            
            
            XYZ2ProPhoto(XYZrev, RGBrev);
            
            XYZ2Lab(XYZrev, Labrev);


            printf("%2d   %5d   ",
                   i, (int)La);
            printf("%4d %4d %4d   ",
                   (int)round(RGBrev[0]),   (int)round(RGBrev[1]),  (int)round(RGBrev[2]));
            printf("%4d %4d %4d   ",
                   (int)round(Labrev[0]),   (int)round(Labrev[1]),  (int)round(Labrev[2]));
            printf("%10.6f %10.6f %10.6f   ",
                   XYZrev[0],               XYZrev[1],              XYZrev[2]);
            printf("%10.6f %11.6f %11.6f ",
                   cam02->color_J,          cam02->color_ja,        cam02->color_jb);
            printf("%10.6f %10.6f %11.6f ",
                   cam02->color_C,          cam02->color_M,         cam02->color_h);
            printf("%11.6f %10.6f ",
                   cam02->color_Q,          cam02->color_s);
            printf("%10.6f %11.6f %11.6f   ",
                   cam02->color_a,          cam02->color_b,         cam02->color_rS);
            printf("%10.6f %11.6f %11.6f   ",
                   cam02->color_ss,         cam02->color_bass,      cam02->color_A);
            printf("%10.6f %10.6f",
                   cam02->color_JJ,         Jab[0]);
            printf("\n");
            printf("\n");
        }
    }
}
    
