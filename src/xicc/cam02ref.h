/*
 * cam02ref: unoptimised, untweaked reference version for testing.
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


struct _cam02ref {
    /* Public: */
    void (*del)(struct _cam02ref *s);       /* We're done with it */
    
    int (*set_view)(
        struct _cam02ref *s,
        ViewingCondition Ev,                /* Enumerated Viewing Condition */
        double Wxyz[3],                     /* Reference/Adapted White XYZ (Y scale 1.0) */
        double La,                          /* Adapting/Surround Luminance cd/m^2 */
        double Yb,                          /* Luminance of Background relative to reference white (range 0.0 .. 1.0) */
        double Lv,                          /* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
                                            /* Ignored if Ev is set */
        double Yf,                          /* Flare as a fraction of the reference white (range 0.0 .. 1.0) */
        double Yg,                          /* Glare as a fraction of the adapting/surround (range 0.0 .. 1.0) */
        double Gxyz[3],                     /* The Glare white coordinates (ie. the Ambient color) */
                                            /* If <= 0 will Wxyz will be used. */
        int hk,                             /* Flag, NZ to use Helmholtz-Kohlrausch effect */
        double hkscale                      /* HK effect scaling factor */
        );
    
    /* Conversions. Return nz on error */
    int (*XYZ_to_cam)(struct _cam02ref *s, double *out, double *in);
    int (*cam_to_XYZ)(struct _cam02ref *s, double *out, double *in);
    
    /* Private: */
    /* Scene parameters */
    ViewingCondition Ev;                    /* Enumerated Viewing Condition */
    double Lv;                              /* Luminance of white in the Viewing/Image cd/m^2 */
    double La;                              /* Adapting/Surround Luminance cd/m^2 */
    double Wxyz[3];                         /* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
    double Yb;                              /* Relative Luminance of Background to reference white (Y range 0.0 .. 1.0) */
    double Yf;                              /* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
    double Yg;                              /* Glare as a fraction of the adapting/surround (Y range 0.0 .. 1.0) */
    double Gxyz[3];                         /* The Glare white coordinates (typically the Ambient color) */
    
    /* Internal parameters */
    double  C;                              /* Surround Impact */
    double Nc;                              /* Chromatic Induction */
    double  F;                              /* Adaptation Degree */
    
    /* Pre-computed values */
    double cc[3][3];                        /* Forward cone and chromatic transform */
    double icc[3][3];                       /* Reverse cone and chromatic transform */
    double crange[3];                       /* ENABLE_COMPR compression range */
    double Va[3], Vb[3], VttA[3], Vttd[3];  /* rgba vectors */
    double dcomp[3];                        /* Vttd in terms of VttA, Va, Vb */
    double Fsc;                             /* Flare+Glare scale */
    double Fisc;                            /* Inverse Flare+Glare scale */
    double Fsxyz[3];                        /* Scaled Flare+Glare color coordinates */
    double rgbW[3];                         /* Sharpened cone response white values */
    double D;                               /* Degree of chromatic adaption */
    double Drgb[3];                         /* Chromatic transformation value */
    double rgbcW[3];                        /* Chromatically transformed white value */
    double rgbpW[3];                        /* Hunt-Pointer-Estevez cone response space white */
    double n;                               /* Background induction factor */
    double nn;                              /* Precomuted function of n */
    double Fl;                              /* Lightness contrast factor ?? */
    double Nbb;                             /* Background brightness induction factors */
    double Ncb;                             /* Chromatic brightness induction factors */
    double z;                               /* Base exponential nonlinearity */
    double rgbaW[3];                        /* Post-adapted cone response of white */
    double Aw;                              /* Achromatic response of white */
    double nldxval;                         /* Non-linearity output value at lower crossover to linear */
    double nldxslope;                       /* Non-linearity slope at lower crossover to linear */
    double nluxval;                         /* Non-linearity value at upper crossover to linear */
    double nluxslope;                       /* Non-linearity slope at upper crossover to linear */
    double lA;                              /* JLIMIT Limited A */
    
    /* Option flags, code not always enabled */
    int hk;                                 /* Use Helmholtz-Kohlrausch effect */
    int hkscale;                            /* [1.0] Scale HK effect up/down from default */
    int trace;                              /* Trace values through computation */
    int retss;                              /* Return ss rather than Jab */
    int range;                              /* (for cam02ref.h) return on range error */
    
    double nldlimit;                        /* value of NLDLIMIT, sets non-linearity lower limit */
    double nldicept;                        /* value of NLDLICEPT, sets straight line intercept with 0.1 output */
    double nlulimit;                        /* value of NLULIMIT, sets non-linearity upper limit (tangent) */
    double ddllimit;                        /* value of DDLLIMIT, sets fwd k3 to k2 limit  */
    double ddulimit;                        /* value of DDULIMIT, sets bwd k3 to k1 limit */
    double ssmincj;                         /* value of SSMINJ, sets cJ minimum value */
    double jlimit;                          /* value of JLIMIT, sets cutover to straight line for J point */
    double hklimit;                         /* value of HKLIMIT, sets HK factor upper limit */
    
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
    double color_ssc;                   /* ss_classic, C/rS */
};

typedef struct _cam02ref cam02ref;


/* Create a cam02 conversion class, with default viewing conditions */
cam02ref *new_cam02ref(void);

