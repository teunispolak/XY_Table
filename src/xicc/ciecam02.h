
//
//  ciecam02.h
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

//struct Correlate {
//    double J, a, b, h;
//    double H, C, Q, M, s;
//};
//extern struct Correlate corr;

struct _ciecam02 {
/* Public: */
	void (*del)(struct _ciecam02 *s);	/* We're done with it */

	int (*set_view)(
		struct _ciecam02 *s,
		double Wxyz[3],	                /* Reference/Adapted White XYZ (Y scale 1.0) */
		double La,		                /* Adapting/Surround Luminance cd/m^2 */
		double Yb,		                /* Luminance of Background relative to reference white (range 0.0 .. 1.0) */
		double Lv		                /* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
                                        /* Ignored if Ev is set */
	);

	/* Forward and backward transforms. Return nz on error */
	int (*XYZ_to_cam)(struct _ciecam02 *s, double *out, double *in);
	int (*cam_to_XYZ)(struct _ciecam02 *s, double *out, double *in);

/* Private: */
	/* Scene parameters */
	double Lv;		                    /* Luminance of white in the Viewing/Image cd/m^2 */
	double La;		                    /* Adapting/Surround Luminance cd/m^2 */
	double Wxyz[3];	                    /* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
	double Yb;		                    /* Relative Luminance of Background to reference white (Y range 0.0 .. 1.0) */

	/* Internal parameters */
	double  C;		                    /* Surround Impact */
	double Nc;		                    /* Chromatic Induction */
	double  F;		                    /* Adaptation Degree */

	/* Pre-computed values */
	double rgbW[3];		                /* Sharpened cone response white values */
	double D;			                /* Degree of chromatic adaption */
	double Drgb[3];		                /* Chromatic transformation value */
	double rgbcW[3];	                /* Chromatically transformed white value */
	double rgbpW[3];	                /* Hunt-Pointer-Estevez cone response space white */
	double n;			                /* Background induction factor */
	double nn;			                /* Precomputed function of n */
	double Fl;			                /* Lightness contrast factor ?? */
	double Nbb;			                /* Background brightness induction factors */
    double Ncb;                         /* Background brightness induction factors */
	double z;			                /* Base exponential nonlinearity */
	double rgbaW[3];	                /* Post-adapted cone response of white */
	double Aw;			                /* Achromatic response of white */
    
    /* Color, i.e. sample values */
    double color_rgb[3];                /* XYZ to CAT02 transformed sample values */
    double color_rgbc[3];               /* Chromatically transformed sample values */
    double color_rgbp[3];               /* CAT02 the HPE transformed sample values */
    double color_rgba[3];               /* Non-linear adapted sample values */
    double color_a;                     /* Red-green opponent dimension */
    double color_b;                     /* Yellow-blue opponent dimension */
    double color_rS;                    /* Preliminary saturation */
    double color_h;                     /* Hue angle */
    double color_e;                     /* Eccentricity factor */
    double color_ttd;                   /* ttd */
    double color_t;                     /* Intermediate magnitude t */
    double color_H;                     /* Hue composition */
    double color_A;                     /* Achromatic response */
    double color_J;                     /* Lightness */
    double color_JJ;                    /* HK Lightess */
    double color_Q;                     /* Brightness */
    double color_C;                     /* Chroma */
    double color_M;                     /* Colorfulness */
    double color_s;                     /* Saturation */
};

typedef struct _ciecam02 ciecam02;


/* Create a ciecam02 conversion class, with default viewing conditions */
ciecam02 *new_ciecam02(void);

