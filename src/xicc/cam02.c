
/* 
 * cam02
 *
 * Color Appearance Model, based on
 * CIECAM02, "The CIECAM02 Color Appearance Model"
 * by Nathan Moroney, Mark D. Fairchild, Robert W.G. Hunt, Changjun Li,
 * M. Ronnier Luo and Todd Newman, IS&T/SID Tenth Color Imaging
 * Conference, with the addition of the Viewing Flare
 * model described on page 487 of "Digital Color Management",
 * by Edward Giorgianni and Thomas Madden, and the
 * Helmholtz-Kohlrausch effect, using the equation from
 * the Bradford-Hunt 96C model as detailed in Mark Fairchild's
 * book "Color Appearance Models". 
 * The Slight modification to the Hunt-Pointer-Estevez matrix
 * recommended by Kuo, Zeis and Lai in IS&T/SID 14th Color Imaging
 * Conference "Robust CIECAM02 Implementation and Numerical
 * Experiment within an ICC Workflow", page 215, together with
 * their matrix formulation of inversion has been adopted.
 *
 * In addition the algorithm has unique extensions to allow
 * it to operate reasonably over an unbounded domain.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2004
 * Version: 1.00
 *
 * This file is based on cam97s3.c by Graeme Gill.
 *
 * Copyright 2004 - 2011 Graeme W. Gill
 * Please refer to COPYRIGHT file for details.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* Note that XYZ values are normalised to 1.0 consistent */
/* with the ICC convention (not 100.0 as assumed by the CIECAM spec.) */
/* Note that all whites are assumed to be normalised (ie. Y = 1.0) */

/*
	TTBD: Should convert to using Timo Kunkel and Erik Reinhard's simplified
	and improved version of CIECAM02
	[ "A Neurophysiology-Inspired Steady-State Color Appearance Model",
	  ie. "CIECAM02-KR" ? ]

	The rgbp compression has it's problems in terms of perceptual
	uniformity. A color with one component near zero might shift
	all the components to -ve values on inverse conversion - ie.
	a 1 DE shift in Jab becomes a masive DE in XYZ/Lab/perceptual,
	with (say) a dark red becomong black because the blue
	value is small. One way around this is to re-introduce
	a flag to turn off perfect symetry by disabling
	expansion on the reverse conversion.

 */

/*
	Various additions and changes have been made to allow the CAM
	conversions to and from an unbounded range of XYZ and Jab values,
	in a (somewhat) geometrically consistent maner. This is because
	such values arise in the process of gamut mapping, and in scanning through
	the grid of PCS values needed to fill in the A2B table of an
	ICC profile. Such values have no correlation to a real color
	value, but none the less need to be handled without causing an
	exception, in a geometrically consistent and reversible
	fashion.

	The changes are:

	The Hunt-Pointer-Estevez matrix has been increased in precission by
	1 digit [Kuo, Zeise & Lai, IS&T/SID 14th Color Imaging Conference pp. 215-219]

	The sharpened cone response matrix third row has been changed from
	[0.0030, 0.0136, 0.9834] to [0.0000, 0.0000, 1.0000]
	[Susstrunk & Brill, IS&T/SID 14th Color Imaging Conference Late Breaking News sublement]

	To prevent wild Jab values as the rgb' values approach zero, a region close to each
	primary ground plane is compressed. This expands the region of reasonable
	behaviour to encompass XYZ values that are found in practice (ie. ProPhotoRGB).

	To prevent excessive curvature of hue in high saturation blue regions
	due to the non-linearity exagerating the difference between
	r & b values, a heuristic is introduced that blends the r & b
	values, reducing this effect. This degrades the agreement
	with the reference CIECAM implementation by about 1DE in this region,
	but greatly improves the behaviour in mapping and clipping from
	highly saturated blues (ie. ProPhotoRGB).

	The post-adaptation non-linearity response has had a straight
	line negative linear extension added.
	The positive region has a tangential linear extension added, so
	that the range of values is not limited.
	An adaptive threshold for the low level linear extension,
	to avoid coordinates blowing up when one of the cone responses
	is large, while the others become negative.
		
	Re-arrange ss equation into separated effects,
	k1, k2, k3, and then limit these to avoid divide by zero
	and sign change errors in fwd and bwd converson, as well
	as avoiding the blue saturation doubling back for low
	luminance levels.

	To avoid chroma and hue angle information being lost when the
	J value used to scale the chroma is 0, and to ensure
	that J = 0, a,b != 0 values have a valid mapping into
	XYZ space, the J value used to multiply Chroma, is limited
	to be equivalent to not less than A == 0.1.

	The Helmholtz-Kohlrausch effect is crafted to have resonable
	effects over the range of J from 0 to 100, and have more
	moderate effects outside this range. 

*/

#include <copyright.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "icc.h"
#include "xcam.h"
#include "cam02.h"
#include "numlib.h"

#undef CLASSIC                  /* [Undef] Tepo Oct 2019: when defined, revert to classic CIECAM02 model */

#undef ENABLE_COMPR		        /* [Def] Enable XYZ compression  */
                                ///* Nov 2019: temporarily switched off for comparison purposes */
#undef ENABLE_DECOMPR		    /* [Undef] Enable XYZ de-compression  */
#undef ENABLE_BLUE_ANGLE_FIX	/* [Def] Limit maximum blue angle */
                                ///* Nov 2019: temporarily switched off for comparison purposes */
#undef ENABLE_DDL			    /* [Def] Enable k1,k2,k3 overall ss limit values (seems to be the best scheme) */
                                /* Nov 2019: temporarily switched off for comparison purposes */
#undef ENABLE_SS			    /* [Undef] Disable overall ss limit values (not the scheme used) */

#define ENTRACE				    /* [Undef] Enable internal value runtime tracing if s->trace != 0 */
#undef DOTRACE				    /* [Undef] Trace anyway (ie. set s->trace = 1) */
#undef DIAG1				    /* [Undef] Print internal value diagnostics for conditions setup */
#undef DIAG2				    /* [Undef] Print internal value diagnostics for each conversion */
#undef TRACKMINMAX			    /* [Undef] Track min/max DD & SS limits (run with locus cam02test) */
#undef DISABLE_MATRIX		    /* Debug - disable matrix & non-lin, wire XYZ to rgba */
#undef DISABLE_SCR			    /* Debug - disable Sharpened Cone Response matrix */
#undef DISABLE_HPE			    /* Debug - disable just Hunt-Pointer_Estevez matrix */
#undef DISABLE_NONLIN		    /* Debug - wire rgbp to rgba */
#undef DISABLE_TTD			    /* Debug - disable ttd vector 'tilt' */
#undef DISABLE_HHKR			    /* Debug - disable Helmholtz-Kohlrausch */

							    /* We reduce the HK effect from the Hunt equation, in the */
							    /* light of real wold experiments. */
#define HHKR_MUL 1.0		    /* [0.25] - Helmholtz-Kohlrausch strength multiplier */
                                /* Nov 2019: temporarily switched on (i.e. disabled) for comparison purposes */

#ifdef ENABLE_COMPR
    # define BC_WHMINY 0.2		/* [0.2] Compression direction minimum Y value */
    # define BC_RANGE_R 0.01	/* [0.02] Set comp. range as prop. of distance to neutral - red */
    # define BC_RANGE_G 0.01	/* [0.02] Set comp. range as prop. of distance to neutral - green*/
    # define BC_RANGE_B 0.01	/* [0.02] Set comp. range as prop. of distance to neutral - blue */
    # define BC_MAXRANGE 0.13	/* [0.13] Maximum compression range */
    # define BC_LIMIT 0.7		/* [0.7] Correction limit (abs. rgbp distance shift) */
#endif

#ifdef ENABLE_BLUE_ANGLE_FIX
    # define BLUE_BL_MAX 0.9	/* [0.9] Sets ultimate blue angle, higher = less angle */
    # define BLUE_BL_POW 3.5	/* [3.5] Sets rate at which blue angle is limited */
#endif

#define NLDLIMIT 0.00001	    /* [0.00001] Non-linearity minimum lower crossover to straight line */
#define NLDICEPT -0.18		    /* [-0.18] Input intercept of straight line with 0.1 output */

#define NLULIMIT 1e5		    /* Non-linearity upper crossover to straight line */

#ifdef ENABLE_SS			    /* [Undef] */
    # define SSLLIMIT 0.22		/* Overall ab scale lower limit */
    # define SSULIMIT 580.0     /* Overall ab scale upper limit */
#endif

#define SYMETRICJ			    /* [Undef] Use symetric power about zero, else straigt line -ve */

#define DDLLIMIT 0.55		    /* [0.55] ab component -k3:k2 ratio limit (must be < 1.0) */
//#define DDULIMIT 0.9993		/* [0.9993] ab component k3:k1 ratio limit (must be < 1.0) */
#define DDULIMIT 0.34		    /* [0.34] ab component k3:k1 ratio limit (must be < 1.0) */
#define SSMINcJ 0.005		    /* [0.005] ab scale cJ minimum value */
#define JLIMIT 0.005		    /* [0.005] J encoding cutover point straight line (0 - 1.0 range) */
#define HKLIMIT 0.7			    /* [0.7] Maximum Helmholtz-Kohlrausch lift out of 1.0 */

#ifdef TRACKMINMAX
    double minss = 1e60;
    double maxss = -1e60;
    double minlrat = 0.0;
    double maxurat = 0.0;
    #define noslots 103
    double slotsd[noslots];
    double slotsu[noslots];
    double minj = 1e38, maxj = -1e38;
#endif

#if defined(ENTRACE) || defined(DOTRACE)
    #define TRACE(xxxx) if (s->trace) printf xxxx ;
#else
    #define TRACE(xxxx)
#endif


static void cam_free(cam02 *s);
static int set_view(struct _cam02 *s, ViewingCondition Ev, double Wxyz[3],
	                double La, double Yb, double Lv, double Yf, double Yg, double Gxyz[3],
					int hk, double hkscale);
static int XYZ_to_cam(struct _cam02 *s, double *Jab, double *xyz);
static int cam_to_XYZ(struct _cam02 *s, double *xyz, double *Jab);
static int QMh_to_XYZ(struct _cam02 *s, double *xyz, double *QMh);
static int JMh_to_XYZ(struct _cam02 *s, double *xyz, double *JMh);

static double spow(double val, double pp) {
	if (val < 0.0)
		return -pow(-val, pp);
	else
		return pow(val, pp);
}


/* Create a cam02 conversion object, with default viewing conditions */
cam02 *new_cam02(void) {
	cam02 *s;

	if ((s = (cam02 *)calloc(1, sizeof(cam02))) == NULL) {
		fprintf(stderr,"cam02: malloc failed allocating object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del          = cam_free;
	s->set_view     = set_view;
	s->XYZ_to_cam   = XYZ_to_cam;
	s->cam_to_XYZ   = cam_to_XYZ;
    s->QMh_to_XYZ   = QMh_to_XYZ;
    s->JMh_to_XYZ   = JMh_to_XYZ;

	/* Initialise default parameters */
	s->hkscale = 1.0;

	/* Set default range handling limits */
	s->nldlimit     = NLDLIMIT;
	s->nldicept     = NLDICEPT;
	s->nlulimit     = NLULIMIT;
	s->ddllimit     = DDLLIMIT;
	s->ddulimit     = DDULIMIT;
	s->ssmincj      = SSMINcJ;
	s->jlimit       = JLIMIT;
	s->hklimit      = 1.0 / HKLIMIT;
    
    #ifdef DOTRACE
        s->trace = 1;
    #endif

    #ifdef TRACKMINMAX
	{
		int i;
		for (i = 0; i < noslots; i++) {
			slotsd[i] = 0.0;
			slotsu[i] = 0.0;
		}
	}
    #endif
    
	return s;
}


/* Delete a cam02 object */
static void cam_free(cam02 *s) {

    #ifdef TRACKMINMAX
	{
		int i;
		for (i = 0; i < noslots; i++) {
			printf("Slot %d = %f, %f\n",i,slotsd[i], slotsu[i]);
		}
		printf("minj = %f, maxj = %f\n",minj,maxj);
		printf("minss = %f\n",minss);
		printf("maxss = %f\n",maxss);
		printf("minlrat = %f\n",minlrat);
		printf("maxurat = %f\n",maxurat);
	}
    #endif

	if (s != NULL)
		free(s);
}


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


/* Set viewing parameters, all white related parameters and all matrices
 * for cam02 object. Return value is always 0 */
static int set_view(
cam02 *s,
ViewingCondition Ev,	/* Enumerated Viewing Condition */
double Wxyz[3],	        /* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
double La,		        /* Adapting/Surround Luminance cd/m^2 */
double Yb,		        /* Relative Luminance of Background to reference white (range 0.0 .. 1.0) */
double Lv,		        /* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
                        /* Ignored if Ev is set to other than vc_none */
double Yf,		        /* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
double Yg,		        /* Flare as a fraction of the adapting/surround (Y range 0.0 .. 1.0) */
double Gxyz[3],	        /* The Glare white coordinates (typically the Ambient color) */
                        /* If <= 0 will Wxyz will be used. */
int hk,			        /* Flag, NZ to use Helmholtz-Kohlrausch effect */
double hkscale	        /* HK effect scaling factor */
) {
	double tt, t1, t2;
	double tm[3][3];
	int i;

    /* Compute the internal parameters from the */
    /* ratio of La to Lv by interpolation */
    if (Ev == vc_none) {
		int i;
        double r, bf;
        /*                 Dark   dim   avg   >avg */
		double t_C[4]  = { 0.525, 0.59, 0.69, 1.0  };
		double t_Nc[4] = { 0.800, 0.95, 1.00, 1.0  };
		double t_F[4]  = { 0.800, 0.90, 1.00, 1.0  };

		if (La < 1e-10) 		                                                            /* Hmm. */
			La = 1e-10;
		r = La/Lv;
		if (r < 0.0)                                                                        /* force 0.0 to 1.0 range */
			r = 0.0;
		else if (r > 1.0)
			r = 1.0;
        /* interpolate */
		if (r < 0.1) {			                                                            /* Dark to Dim */
			i = 0;
			bf = r/0.1;
		} else if (r < 0.2) {	                                                            /* Dim to Average */
			i = 1;
			bf = (r - 0.1)/0.1;
		} else {				                                                            /* Average to above average */
			i = 2;
			bf = (r - 0.2)/0.8;
		}
		s->C  = t_C[i]  * (1.0 - bf) + t_C[i+1]  * bf;
		s->Nc = t_Nc[i] * (1.0 - bf) + t_Nc[i+1] * bf;
		s->F  = t_F[i]  * (1.0 - bf) + t_F[i+1]  * bf;
	} else {
		/* Compute the internal parameters by category
		 * Fake up Lv according to category */
		switch(Ev) {
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
			default:
				s->C = 0.69;
				s->Nc = 1.0;
				s->F = 1.0;
				Lv = La/0.2; 
				break;
			case vc_cut_sheet:
				s->C = 0.41;
				s->Nc = 0.8;
				s->F = 0.8;
				Lv = La/0.02; 	// ???
				break;
		}
	}

    
	/* Transfer parameters to the object */
	s->Ev = Ev;
	s->Wxyz[0] = Wxyz[0];
	s->Wxyz[1] = Wxyz[1];
	s->Wxyz[2] = Wxyz[2];
	s->La = La;
	s->Yb = Yb > 0.005 ? Yb : 0.005;	/* Set minimum to avoid divide by 0.0 */
	s->Lv = Lv;
	s->Yf = Yf;
	s->Yg = Yg;
	if (Gxyz[0] > 0.0 && Gxyz[1] > 0.0 && Gxyz[2] > 0.0) {
		tt = Wxyz[1]/Gxyz[1];		/* Scale to white ref white */
		s->Gxyz[0] = tt * Gxyz[0];
		s->Gxyz[1] = tt * Gxyz[1];
		s->Gxyz[2] = tt * Gxyz[2];
	} else {
		s->Gxyz[0] = Wxyz[0];
		s->Gxyz[1] = Wxyz[1];
		s->Gxyz[2] = Wxyz[2];
	}
    s->hk = 1; //hk;
    s->hkscale = 1.0; //hkscale;

	/* The rgba vectors */
	/* used in eqn (14) */
    s->Va[0] = 1.0;
	s->Va[1] = -12.0/11.0;
	s->Va[2] = 1.0/11.0;

	/* used in eqn (15) */
    s->Vb[0] = 1.0/9.0;
	s->Vb[1] = 1.0/9.0;
	s->Vb[2] = -2.0/9.0;

	/* used in eqn (20), first part */
    s->VttA[0] = 2.0;
	s->VttA[1] = 1.0;
	s->VttA[2] = 1.0/20.0;

    /* only used for classic CIECAM02: coeffs in denominator of eqn (16)) */
	s->Vttd[0] = 1.0;
	s->Vttd[1] = 1.0;
	s->Vttd[2] = 21.0/20.0;

	/* Vttd in terms of the VttA, Va and Vb vectors */
	s->dcomp[0] = 1.0;
	s->dcomp[1] = -11.0/23.0;
	s->dcomp[2] = -108.0/23.0;

    
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

	/* Rescale so that the sum of the flare and the input doesn't exceed white */
	s->Fsc = s->Wxyz[1]/(s->Fsxyz[1] + s->Wxyz[1]);
	s->Fsxyz[0] *= s->Fsc;
	s->Fsxyz[1] *= s->Fsc;
	s->Fsxyz[2] *= s->Fsc;
	s->Fisc = 1.0/s->Fsc;
    

    /* Sharpened cone response white values */
    #ifndef DISABLE_SCR
        s->rgbW[0] =  0.7328 * s->Wxyz[0] + 0.4296 * s->Wxyz[1] - 0.1624 * s->Wxyz[2];
        s->rgbW[1] = -0.7036 * s->Wxyz[0] + 1.6975 * s->Wxyz[1] + 0.0061 * s->Wxyz[2];
//        s->rgbW[2] =  0.0000 * s->Wxyz[0] + 0.0000 * s->Wxyz[1] + 1.0000 * s->Wxyz[2];
    
        s->rgbW[2] =  0.0030 * s->Wxyz[0] + 0.0136 * s->Wxyz[1] + 0.9834 * s->Wxyz[2];
    #else
        s->rgbW[0] = s->Wxyz[0];
        s->rgbW[1] = s->Wxyz[1];
        s->rgbW[2] = s->Wxyz[2];
    #endif

    
	/* Degree of chromatic adaptation */
	s->D = s->F * (1.0 - exp((-s->La - 42.0)/92.0)/3.6);                                        /* eqn (8) */


    /* Precompute Chromatic transform values
     *
     * Interestingly, "Color Appearance Models" para 16.3 indicates that
     * not Yw (i.e. s->Wxyz[1]) shall be used, but 100. */
	s->Drgb[0] = s->D * (s->Wxyz[1]/s->rgbW[0]) + 1.0 - s->D;
	s->Drgb[1] = s->D * (s->Wxyz[1]/s->rgbW[1]) + 1.0 - s->D;
	s->Drgb[2] = s->D * (s->Wxyz[1]/s->rgbW[2]) + 1.0 - s->D;

    
	/* Chromaticaly transformed white value */
	s->rgbcW[0] = s->Drgb[0] * s->rgbW[0];
	s->rgbcW[1] = s->Drgb[1] * s->rgbW[1];
	s->rgbcW[2] = s->Drgb[2] * s->rgbW[2];
	

    /* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
    #ifndef DISABLE_HPE
//        s->rgbpW[0] =  0.7409744840453773 * s->rgbcW[0]
//                    +  0.2180245944753982 * s->rgbcW[1]
//                    +  0.0410009214792244 * s->rgbcW[2];
//        s->rgbpW[1] =  0.2853532916858801 * s->rgbcW[0]
//                    +  0.6242015741188157 * s->rgbcW[1]
//                    +  0.0904451341953042 * s->rgbcW[2];
//        s->rgbpW[2] = -0.0096276087384294 * s->rgbcW[0]
//                    -  0.0056980312161134 * s->rgbcW[1]
//                    +  1.0153256399545427 * s->rgbcW[2];
    
        s->rgbpW[0] =  0.7409792 * s->rgbcW[0]
                    +  0.2180250 * s->rgbcW[1]
                    +  0.0410058 * s->rgbcW[2];
        s->rgbpW[1] =  0.2853532 * s->rgbcW[0]
                    +  0.6242014 * s->rgbcW[1]
                    +  0.0904454 * s->rgbcW[2];
        s->rgbpW[2] = -0.0096280 * s->rgbcW[0]
                    -  0.0056980 * s->rgbcW[1]
                    +  1.0153260 * s->rgbcW[2];
    
    #else
        s->rgbpW[0] = s->rgbcW[0];
        s->rgbpW[1] = s->rgbcW[1];
        s->rgbpW[2] = s->rgbcW[2];
    #endif


    /* Create combined cone and chromatic transform matrix:
     * Spectrally sharpened cone responses matrix */
    #ifndef DISABLE_SCR
        s->cc[0][0] =  0.7328; s->cc[0][1] = 0.4296; s->cc[0][2] = -0.1624;
        s->cc[1][0] = -0.7036; s->cc[1][1] = 1.6975; s->cc[1][2] =  0.0061;
//        s->cc[2][0] =  0.0000; s->cc[2][1] = 0.0000; s->cc[2][2] =  1.0000;
    
        s->cc[2][0] =  0.0030; s->cc[2][1] = 0.0136; s->cc[2][2] =  0.9834;
    
    #else
        s->cc[0][0] = 1.0; s->cc[0][1] = 0.0; s->cc[0][2] = 0.0;
        s->cc[1][0] = 0.0; s->cc[1][1] = 1.0; s->cc[1][2] = 0.0;
        s->cc[2][0] = 0.0; s->cc[2][1] = 0.0; s->cc[2][2] = 1.0;
    #endif
    
    
	/* Chromaticaly transformed sample values */
	icmSetUnity3x3(tm);
	tm[0][0] = s->Drgb[0];
	tm[1][1] = s->Drgb[1];
	tm[2][2] = s->Drgb[2];
	icmMul3x3(s->cc, tm);
	

    /* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
    #ifndef DISABLE_HPE
//        tm[0][0] =  0.7409744840453773;
//        tm[0][1] =  0.2180245944753982;
//        tm[0][2] =  0.0410009214792244;
//        tm[1][0] =  0.2853532916858801;
//        tm[1][1] =  0.6242015741188157;
//        tm[1][2] =  0.0904451341953042;
//        tm[2][0] = -0.0096276087384294;
//        tm[2][1] = -0.0056980312161134;
//        tm[2][2] =  1.0153256399545427;
    
        tm[0][0] =  0.7409792;
        tm[0][1] =  0.2180250;
        tm[0][2] =  0.0410058;
        tm[1][0] =  0.2853532;
        tm[1][1] =  0.6242014;
        tm[1][2] =  0.0904454;
        tm[2][0] = -0.0096280;
        tm[2][1] = -0.0056980;
        tm[2][2] =  1.0153260;
    #endif
    
    icmMul3x3(s->cc, tm);

    
	/* Create inverse combined cone and chromatic transform matrix: */
	icmInverse3x3(s->icc, s->cc);		/* Hmm. we assume it cannot fail */

    
    #ifdef ENABLE_COMPR
        /* Compression ranges */
        s->crange[0] = BC_RANGE_R;
        s->crange[1] = BC_RANGE_G;
        s->crange[2] = BC_RANGE_B;
    #endif /* ENABLE_COMPR */

    
	/* Background induction factor */
	s->n = s->Yb/ s->Wxyz[1];                                                               /* eqn (3) */
	s->nn = pow(1.64 - pow(0.29, s->n), 0.73);	                                            /* last term of eqn (23), Chroma C */

    
	/* Lightness contrast factor */
	{
		double k;
		k = 1.0 / (5.0 * s->La + 1.0);                                                      /* eqn (1) */
		s->Fl = 0.2 * pow(k , 4.0) * 5.0 * s->La
		      + 0.1 * pow(1.0 - pow(k , 4.0) , 2.0) * pow(5.0 * s->La , 1.0/3.0);           /* eqn (2) */
	}
    

	/* Background and Chromatic brightness induction factors */
	s->Nbb   = 0.725 * pow(1.0/s->n, 0.2);                                                  /* eqn (4) */
	s->Ncb   = s->Nbb;                                                                      /* eqn (4) */

    
	/* Base exponential nonlinearity */
	s->z = 1.48 + pow(s->n , 0.5);                                                          /* eqn (5) */

    
	/* Post-adapted cone response of white */
	tt = pow(s->Fl * s->rgbpW[0], 0.42);
	s->rgbaW[0] = 400.0 * tt / (tt + 27.13) + 0.1;
	tt = pow(s->Fl * s->rgbpW[1], 0.42);
	s->rgbaW[1] = 400.0 * tt / (tt + 27.13) + 0.1;
	tt = pow(s->Fl * s->rgbpW[2], 0.42);
	s->rgbaW[2] = 400.0 * tt / (tt + 27.13) + 0.1;

    
	/* Achromatic response of white */
	s->Aw = (s->VttA[0] * s->rgbaW[0]
	      +  s->VttA[1] * s->rgbaW[1]
	      +  s->VttA[2] * s->rgbaW[2] - 0.305) * s->Nbb;

    
	/* Non-linearity lower crossover output value */
	tt = pow(s->Fl * s->nldlimit, 0.42);
	s->nldxval = 400.0 * tt / (tt + 27.13) + 0.1;

    
	/* Non-linearity lower crossover slope from lower crossover */
	/* to intercept with 0.1 output */
	s->nldxslope = (s->nldxval - 0.1)/(s->nldlimit - s->nldicept);

    
	/* Non-linearity upper crossover value */
	tt = pow(s->Fl * s->nlulimit, 0.42);
	s->nluxval = 400.0 * tt / (tt + 27.13) + 0.1;

    
	/* Non-linearity upper crossover slope, set to asymtope */
	t1 = s->Fl * s->nlulimit;
	t2 = pow(t1, 0.42) + 27.13;
	s->nluxslope = 0.42 * s->Fl * 400.0 * 27.13 / (pow(t1, 0.58) * t2 * t2);


	/* Limited A value at J = JLIMIT */
	s->lA = pow(s->jlimit, 1.0/(s->C * s->z)) * s->Aw;                                      /* via inverse of eqn (21) */

    
    /* Source and destination values for Fl
     * The following is a bit tricky
     *
     * In order to keep M constant between source and destination when building a profile,
     * source Fl must be known. The value is stored as src_Fl
     *
     * When building a profile, the first time the scene parameters are printed is to dump
     * the source viewing conditions to the log file. The second time is to print the
     * destination viewing conditions.
     *
     */

    
    if (theSource == 1 && srcFl == 0.0) {
        srcFl = s->Fl;
        s->src_Fl = s->Fl;
    }
    if (theDest == 1 && dstFl == 0.0) {
        dstFl = s->Fl;
        s->dst_Fl = s->Fl;
    }
  
    
    if (printCAM == 1) {

    //#ifdef DIAG1
        printf("Scene parameters:\n");
        printf("Viewing condition Ev = %d\n",s->Ev);
        printf("Ref white Wxyz = %f %f %f\n", s->Wxyz[0], s->Wxyz[1], s->Wxyz[2]);
        printf("Relative luminance of background Yb = %f\n", s->Yb);
        printf("Adapting luminance La = %f\n", s->La);
        printf("Flare Yf = %f\n", s->Yf);
        printf("Glare Yg = %f\n", s->Yg);
        printf("Glare color Gxyz = %f %f %f\n", s->Gxyz[0], s->Gxyz[1], s->Gxyz[2]);

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
        printf("                            nn = %f\n", s->nn);
        printf("Source lightness contrast factor Fl = %f\n", s->src_Fl);
        printf("Destination lightness contrast factor Fl = %f\n", s->dst_Fl);
        printf("Background brightness induction factor Nbb = %f\n", s->Nbb);
        printf("Chromatic brightness induction factor Ncb = %f\n", s->Ncb);
        printf("Base exponential nonlinearity z = %f\n", s->z);
        printf("Post adapted cone response white rgbaW = %f %f %f\n", s->rgbaW[0], s->rgbaW[1], s->rgbaW[2]);
        printf("Achromatic response of white Aw = %f\n", s->Aw);
        printf("\n");
    }
    //#endif
	
    
    return 0;
}


/* Conversions. Return values are always 0 */
static int XYZ_to_cam(
struct _cam02 *s,
double Jab[3],
double XYZ[3]
) {
	int i;
	double xyz[3], XYZi[3], rgbp[3], rgba[3];
	double a, b, ja, jb, J, JJ, C, h, e, A, ss;
	double ttA, rS, cJ, tt;
	double k1, k2, k3, ss_classic;
    double ba_ss = 0;
    double t, ttd;
    double H, M, Q, sat, temp;
	int clip = 0;

    TRACE(("\nCIECAM02 Forward conversion:\n"))
    TRACE(("XYZ = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))


    /* Store input XYZ sample values for printing at the end */
    //#ifdef DIAG2                                                                          /* normally skipped */
        XYZi[0] = XYZ[0];
        XYZi[1] = XYZ[1];
        XYZi[2] = XYZ[2];
    //#endif


    #ifdef DISABLE_MATRIX                                                                   /* normally skipped */
        rgba[0] = XYZ[0];
        rgba[1] = XYZ[1];
        rgba[2] = XYZ[2];
    #else                                                                                   /* normally executed */
        /* Add in flare */
        xyz[0] = s->Fsc * XYZ[0] + s->Fsxyz[0];
        xyz[1] = s->Fsc * XYZ[1] + s->Fsxyz[1];
        xyz[2] = s->Fsc * XYZ[2] + s->Fsxyz[2];

        TRACE(("XYZ inc flare = %f %f %f\n",xyz[0], xyz[1], xyz[2]))

    
        /* Spectrally sharpened cone response sample values rgbp[].
         * Apply chromatic transform and transform from spectrally sharpened
         * to Hunt-Pointer_Estevez cone space, all in one step. */
        icmMulBy3x3(rgbp, s->cc, xyz);
        TRACE(("rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]))

    
        /* Compress rgbp[] values
         * Avoid out of cam02 gamut behaviour by compressing
         * rgbp[] to prevent it from becoming less than zero. */
        #ifdef ENABLE_COMPR                                                                 /* normally executed */
            {
                double tt;			/* Temporary */
                double wrgb[3];		/* White target */

                /* Make white target white point with same Y value */
                tt = xyz[1] > BC_WHMINY ? xyz[1] : BC_WHMINY;	                            /* Limit to minimum Y */
                icmScale3(wrgb, s->rgbpW, tt/s->Wxyz[1]);	                                /* White target at same Y */
                TRACE(("wrgb %f %f %f\n", wrgb[0], wrgb[1], wrgb[2]))

                /* Compress r,g,b in turn */
                for (i = 0; i < 3; i++) {
                    double cv;		/* Compression value */
                    double ctv;		/* Compression target value */
                    double cd;		/* Compression displacement needed */
                    double cvec[3];	/* Normalized correction vector */
                    double isec[3];	/* Intersection with plane */
                    double offs;	/* Offset of point from orgin*/
                    double range;	/* Threshold to start compression */
                    double asym;	/* Compression asymtope */

                    /* Compute compression direction as vector towards white
                     * (We did try correcting in a blend of limit plane normal and white,
                     *  but compressing towards white seems to be the best.) */
                    icmSub3(cvec, wrgb, rgbp);					                            /* Direction of white target */
                    TRACE(("ch %d, rgbp %f %f %f\n", i, rgbp[0], rgbp[1], rgbp[2]))
                    TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))

                    if (cvec[i] < 1e-9) {		                                            /* compression direction can't correct this coord */
                        TRACE(("Intersection with limit plane failed\n"))
                        continue;
                    }

                    /* Scale compression vector to make it move a unit in normal direction */
                    icmScale3(cvec, cvec, 1.0/cvec[i]);		                                /* Normalized vector to white */
                    TRACE(("ncvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))

                    /* Compute intersection of correction direction with this limit plane
                     * (This corresponds with finding displacement of rgbp by cvec
                     *  such that the current coord value = 0) */
                    icmScale3(isec, cvec, -rgbp[i]);		                                /* (since cvec[i] == 1.0) */
                    icmAdd3(isec, isec, rgbp);
                    TRACE(("isec %f %f %f\n", isec[0], isec[1], isec[2]))

                    /* Compute distance from intersection to origin */
                    offs = pow(icmNorm3(isec), 0.85);
                    range = s->crange[i] * offs;	                                        /* Scale range by distance to origin */
                    if (range > BC_MAXRANGE)		                                        /* so that it tapers down as we approach it */
                        range = BC_MAXRANGE;		                                        /* and limit maximum */

                    /* Aiming above plane when far from origin,
                     * but below plane at the origin, so that black isn't affected. */
                    asym = range - 0.2 * (range + (0.01 * s->crange[i]));
                    ctv = cv = rgbp[i];		                                                /* Distance above/below limit plane */
                    TRACE(("ch %d, offs %f, range %f asym %f, cv %f\n",i, offs,range,asym,cv))

                    if (cv < (range - 1e-12)) {		                                        /* Need to compress */
                        double aa, bb;
                        aa = 1.0/(range - cv);
                        bb = 1.0/(range - asym);
                        ctv = range - 1.0/(aa + bb);
                        cd = ctv - cv;				                                        /* Displacement needed */
                        if (cd > BC_LIMIT)
                            cd = BC_LIMIT;
                        TRACE(("ch %d cd = %f, scaled cd %f\n",i,cd,cd))

                        /* Apply correction */
                        icmScale3(cvec, cvec, cd);			                                /* Compression vector */
                        icmAdd3(rgbp, rgbp, cvec);			                                /* Compress by displacement */
                        TRACE(("rgbp after comp. = %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
                    }
                }
            }
        #endif


        /* Limit maximum blue angle */
        #ifdef ENABLE_BLUE_ANGLE_FIX                                                        /* normally executed */
            ba_ss = rgbp[0] + rgbp[1] + rgbp[2];
            if (ba_ss < 1e-9) {
                ba_ss = 0.0;
            } else {
                ba_ss = (rgbp[2]/ba_ss - 1.0/3.0) * 3.0/2.0;
                if (ba_ss > 0.0)
                    ba_ss = BLUE_BL_MAX * pow(ba_ss, BLUE_BL_POW);
            }
            if (ba_ss < 0.0)
                ba_ss = 0.0;
            else if (ba_ss > 1.0)
                ba_ss = 1.0;
            tt = 0.5 * (rgbp[0] + rgbp[1]);
            rgbp[0] = ba_ss * tt + (1.0 - ba_ss) * rgbp[0];
            rgbp[1] = ba_ss * tt + (1.0 - ba_ss) * rgbp[1];
            TRACE(("rgbp after blue fix ba_ss = %f, rgbp = %f %f %f\n",ba_ss, rgbp[0], rgbp[1], rgbp[2]))
        #endif


        /* Post-adapted cone response sample values rgba[] */
        #ifdef DISABLE_NONLIN                                                               /* normally skipped */
            for (i = 0; i < 3; i++) {
                rgba[i] = 400.0/27.13 * rgbp[i];
            }
        #else	                                                                            /* normally executed */
            /* rgba[] has a minimum value of 0.1 for XYZ[] = 0 and no flare.
             * We add a negative linear region, plus a linear segment at
             * the end of the +ve conversion to allow numerical handling of a
             * very wide range of values. */
            for (i = 0; i < 3; i++) {
                if (rgbp[i] < s->nldlimit) {
                    rgba[i] = s->nldxval + s->nldxslope * (rgbp[i] - s->nldlimit);
                } else {
                    if (rgbp[i] <= s->nlulimit) {
                        tt = pow(s->Fl * rgbp[i], 0.42);
                        rgba[i] = 400.0 * tt / (tt + 27.13) + 0.1;                          /* eqn (13) */
                    } else {
                        rgba[i] = s->nluxval + s->nluxslope * (rgbp[i] - s->nlulimit);
                    }
                }
            }
        #endif
        TRACE(("rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]))
    #endif  /* end of scope for DISABLE_MATRIX */


    /* Preliminary red-green & yellow-blue opponent dimensions
	 * a, b & ttd form an (almost) orthogonal coordinate set.
	 * ttA is in a different direction */
	a = s->Va[0] * rgba[0] + s->Va[1] * rgba[1] + s->Va[2] * rgba[2];                       /* eqn (14) */
	b = s->Vb[0] * rgba[0] + s->Vb[1] * rgba[1] + s->Vb[2] * rgba[2];                       /* eqn (15) */


	/* Restricted Saturation rS to stop division by zero
	 * The exact value isn't important because the numerator dominates
     * as a,b aproach 0 */
	rS = sqrt(a * a + b * b);
	if (rS < DBL_EPSILON)
		rS = DBL_EPSILON;
    TRACE(("a = %f, b = %f, rS = %f\n", a, b, rS))

    
    /* Final hue angle */
    h = (180.0/DBL_PI) * atan2(b,a);                                                        /* eqn (17) */
    h = (h < 0.0) ? h + 360.0 : h;
    
    
    /* Eccentricity factor */
    e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));              /* eqn (18) */
    
    
    /* ttd and t, not in the original cam02.c */
    ttd = s->Vttd[0] * rgba[0] + s->Vttd[1] * rgba[1] + s->Vttd[2] * rgba[2];
    t = e * rS / ttd;                                                                       /* eqn (16) */
  
    
    /* Preliminary Achromatic response ttA and Achromatic response A
     * Note that the minimum values of rgba[] for XYZ = 0 is 0.1,
     * hence magic 0.305 below comes from the following weighting of rgba[],
     * to base A at 0.0
     * Deal with the core rgb to A, S & h conversion: */
    ttA = s->VttA[0] * rgba[0] + s->VttA[1] * rgba[1] + s->VttA[2] * rgba[2];
    A = (ttA - 0.305) * s->Nbb;                                                             /* eqn (20) */
    TRACE(("ttd = %f, t = %f, h = %f, e = %f, ttA = %f, A = %f\n", ttd, t, h, e, ttA, A))


    /* Lightness J, derived directly from Achromatic response A */
    #ifndef SYMETRICJ                                                                       /* normally skipped */
        if (A >= s->lA) {
            J = pow(A/s->Aw, s->C * s->z);
        } else {
            J = s->jlimit/s->lA * A;
            TRACE(("Limited Achromatic response A to straight line\n"))
        }
    #else                                                                                   /* normally executed */
        if (A >= 0.0) {
            J = pow(A/s->Aw, s->C * s->z);		                                            /* eqn (21) */
        } else {
            J = -pow(-A/s->Aw, s->C * s->z);                                                /* eqn (21), if A < 0 */
            TRACE(("Symmetric Achromatic response A\n"))
        }
    #endif

    
	/* cJ, the constrained J (positive, non-zero) */
	if (A > 0.0) {
		cJ = pow(A/s->Aw, s->C * s->z);
		if (cJ < s->ssmincj)
			cJ = s->ssmincj;
	} else {
		cJ = s->ssmincj;
	}
	TRACE(("J = %f, cJ = %f\n",J,cJ))

    
	/* a, b scaling components k1, k2, k3 */
    k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rS, 1.0/9.0);
    k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;

    
    /* Force ss and thus Chroma to remain at their average viewing conditions' values */
    /* oplossen met #ifdef het een of ander? */
    double cavg = 0.69;
    double Ncavg = 1.0;
    double eavg = (12500.0/13.0 * Ncavg * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));
    double Javg = pow(A/s->Aw, cavg * s->z);
    double k1avg = pow(s->nn, 1.0/0.9) * eavg * pow(Javg, 1.0/1.8)/pow(rS, 1.0/9.0);
    double k2avg = pow(Javg, 1.0/(cavg * s->z)) * s->Aw/s->Nbb + 0.305;
    
	k3 = s->dcomp[1] * a + s->dcomp[2] * b;
    double ssavg = pow(k1avg/(k2avg + k3), 0.9);
    
    TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, pow(k1/(k2 + k3), 0.9)))

    
    #ifdef TRACKMINMAX                                                                      /* normally skipped */
	{
		int sno;
		double lrat, urat;

		ss = pow(k1/(k2 + k3), 0.9);
		if (ss < minss)
			minss = ss;
		if (ss > maxss)
			maxss = ss;
		lrat = -k3/k2;
		urat =  k3 * pow(ss, 10.0/9.0) / k1; 
		if (lrat > minlrat)
			minlrat = lrat;
		if (urat > maxurat)
			maxurat = urat;
		/* Record distribution of ss min/max vs. J for
		 * regions outside a,b == 0 */
		sno = (int)(J * 100.0 + 0.5);
		if (sno < 0) {
			sno = 101;
			if (J < minj)
				minj = J;
		} else if (sno > 100) {
			sno = 102;
			if (J > maxj)
				maxj = J;
		}
		if (slotsd[sno] < lrat)
			slotsd[sno] = lrat;
		if (slotsu[sno] < urat)
			slotsu[sno] = urat;
	}
    #endif

    
    /* Limit k3 values */
    #ifdef ENABLE_DDL                                                                       /* normally executed */
        /* Limit ratio of k3 to k2 to stop zero or -ve ss */
        if (k3 < -k2 * s->ddllimit) {
            k3 = -k2 * s->ddllimit;
            TRACE(("k3 limited to %f due to k3:k2 ratio, ss = %f\n",k3,pow(k1/(k2 + k3), 0.9)))
            clip = 1;
        }
        /* See if there is going to be a problem in bwd, and limit k3 if there is */
        if (k3 > (k2 * s->ddulimit/(1.0 - s->ddulimit))) {
            k3 = (k2 * s->ddulimit/(1.0 - s->ddulimit));
            TRACE(("k3 limited to %f to allow for bk3:bk1 bwd limit\n",k3))
            clip = 1;
        }
    #endif

    
    /* ab scaling factor ss */
    #ifdef DISABLE_TTD                                                                      /* normally skipped */
        ss = pow((k1/k2), 0.9);
    #else                                                                                   /* normally executed */
        ss = pow(k1/(k2 + k3), 0.9);
    #endif

    
    /* Limit ab scaling factor ss */
    #ifdef ENABLE_SS                                                                        /* normally skipped */
        if (ss < SSLLIMIT)
            ss = SSLLIMIT;
        else if (ss > SSULIMIT)
            ss = SSULIMIT;
    #endif

    
    /* Scaling
     *
     * In the Argyll cam02.c implementation, ss is used to scale a, b to ja, jb:
     *      ja = a * ss;
     *      jb = b * ss;
     *
     * Next, ja and jb are used to calculate C:
     *      C  = sqrt (ja * ja + jb * jb);
     *
     * which is equivalent to
     *      C  = ss * rS;
     *      ss =  C / rS;
     *
     * In the classic CIECAM02 model [1] C is calculated from
     *      C  = pow(t, 0.9) * sqrt(J) * s->nn;  (eqn 23)
     *
     * Since for classic CIECAM02 a similar relationship exists between C and rS, we can write
     *      ss_classic = C/rS;
     *      ss_classic = pow(t, 0.9) * sqrt(J) * s->nn / rS;
     *
     *
     * For real-world colors, ss and ss_classic appear to be identical. This could be expected
     * because CIECAM02 behaves well for real world colors, so we don't need any scaling for the model
     * to return meaningful results.
     *
     * For non real-world colors (let's loosely define them as colors outside the CIE xy 'horseshoe' gamut),
     * ss and ss_classic are not identical. The blue primary color of ProPhotoRGB can serve as an example.
     * Here, classic CIECAM02 fails and ss is required to avoid runaway behavior.
     *
     *
     *  [1] The CIECAM02 Color Appearance Model
     *  [2] Billy Bigg's CIECAM02 implementation
     *  [3] Spreadsheet
     */
    
    
    /* Scale a, b and calculate C
     *
     * ss replaces the classic calculation of C:
     * C  = pow(t, 0.9) * sqrt(J) * s->nn;                                                 /* eqn (23) */

    
    /* When using viewing condition dark or dim to compensate for
     * display contrast range: keep M constant. Therefor, keep C constant */
    double Qavg;
    Qavg = (4.0 / cavg) * sqrt(Javg) * (s->Aw + 4.0) * pow(s->Fl, 0.25);
    Q    = (4.0 / s->C) * sqrt(cJ)   * (s->Aw + 4.0) * pow(s->Fl, 0.25);
    ss = ss * sqrt(Q/Qavg);
    ja = a * ss;
    jb = b * ss;
    C  = ss * rS;                                                                           /* or sqrt(ja^2 + jb^2) */

//    printf ("cJ = %f, C = %f, Q = %f, ss = %f, Javg = %f Qavg = %f ssavg = %f sqrt = %f\n", cJ, C, Q, ss, Javg, Qavg, ssavg, sqrt(Q/Qavg));
    
    TRACE(("ss = %f, A = %f, J = %f, C = %f, h = %f\n", ss, A, J, C, h))

    
    /* Helmholtz-Kohlrausch effect JJ */
    JJ = J;                                         /* offe... JJ = cJ? */
//    #ifndef DISABLE_HHKR                                                                    /* normally executed */
//        if (s->hk && J < 1.0) {
//            double strength = s->hkscale * HHKR_MUL;
////            double kk = strength * C / 300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0)) / 180.0);
//            double kk = strength * C / 300.0; // * sin(DBL_PI * fabs(0.5 * (h - 90.0)) / 180.0);
//            if (kk > 1e-6)                                                                     /* Limit kk to a reasonable range */
//                kk = 1.0 / (s->hklimit + 1.0 / kk);
//            JJ = J + (1.0 - (J > 0.0 ? J : 0.0)) * kk;
//            TRACE(("JJ = %f from J = %f, kk = %f\n", JJ, J, kk))
////            J = JJ;                                                                         /* would have been easier */
//        }
//    #endif

    
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
    
    
    /* Brightness Q */
    Q = (4.0 / s->C) * sqrt(JJ) * (s->Aw + 4.0) * pow(s->Fl, 0.25);                         /* eqn (22) */
    
    
    /* Colorfulness M */
    M = C * pow(s->Fl, 0.25);                                                               /* eqn (24) */
    
    
    /* Saturation s (well, actually sat, because s is already in use) */
    sat = 100 * sqrt(M / Q);                                                                /* eqn (25) */

    
    /* Jab[] return values, in 0 .. 100 range */
    Jab[0] = JJ  * 100;
    Jab[1] = ja * 1;
    Jab[2] = jb * 1;

    TRACE(("Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]))
    TRACE(("\n"))
 
    
    /* Transfer internal sample values */
    for (int i = 0; i < 3; i++) {
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
    s->color_J   = J * 100;
    s->color_JJ  = JJ* 100;
    s->color_Q   = Q;
    s->color_C   = C;
    s->color_M   = M;
    s->color_s   = sat;
    s->color_ss  = ss;
    s->color_bass = ba_ss;
    s->color_ja  = ja;
    s->color_jb  = jb;

    
    if (printCAM == 1) {
        #ifdef DIAG2                                                                            /* normally skipped */
            printf("Processing XYZ->Jab (Argyll default cam02):\n");
            printf("XYZ = %f %f %f\n", XYZi[0], XYZi[1], XYZi[2]);
            printf("Including flare XYZ = %f %f %f\n", xyz[0]*100, xyz[1]*100, xyz[2]*100);
            printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
            printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
            printf("Prelim red green a = %f, b = %f\n", a, b);
            printf("Hue angle h = %f\n", h);
            printf("Eccentricity factor e = %f\n", e);
            printf("Preliminary magnitude t = %f\n", t);
            printf("Achromatic response A = %f\n", A);
            printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ * 100);
            printf("Prelim Saturation rS = %f\n", rS);
            printf("a, b scale factor ss = %f\n", ss);
            printf("Chroma C = %f\n", C);
            printf("Brightness Q = %f\n",Q);
            printf("Colorfulness M = %f\n",M);
            printf("Saturation s = %f\n",sat);
            printf("Jab = %f %f %f\n", Jab[0], Jab[1], Jab[2]);
            printf("\n");
        #endif
    }

	return 0;
}


static int cam_to_XYZ(
struct _cam02 *s,
double XYZ[3],
double Jab[3]
) {
	int i;
	double xyz[3], Jabi[3], rgbp[3], rgba[3];
	double a, b, ja, jb, rS, J, JJ=0, C=0, rC, h, hr, e, A, ss;
    double Q, M=0,sat;
    double Cb;
    double ba_ss = 0.0;
	double tt, cJ, ttA;
	double k1, k2, k3;
    double kk=0, strength=0;

#define lomp
#ifdef lomp
    double JMh[3], XYZrev[3];
    
    JJ = Jab[0] / 100;
    ja = Jab[1];
    jb = Jab[2];

    
//    /* Keep M constant */
//    ja = ja * pow(srcFl/dstFl, 0.25);
//    jb = jb * pow(srcFl/dstFl, 0.25);

    
    /* Keep M constant , well, uhh, not entirely */
//    ja = ja * sqrt(pow(srcFl/dstFl, 0.25));
//    jb = jb * sqrt(pow(srcFl/dstFl, 0.25));
    
    
    /* Chroma C */
    C = sqrt(ja * ja + jb * jb);

    
    /* Hue angle h */
    h = (180.0/DBL_PI) * atan2(jb, ja);
    h = (h < 0.0) ? h + 360.0 : h;


    /* Colorfulness M */
    M = C * pow(s->Fl, 0.25);
    
    
//    printf("J, ja, jb, C, M, h, srcFl, dstFl %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", J, ja, jb,
//           C, M, h, srcFl, dstFl);
    

    /* Preliminary Restricted chroma, always +ve and NZ */
    /* (The exact value isn't important because the numerator dominates as a,b aproach 0) */
    rC = C;
    if (rC < DBL_EPSILON)
        rC = DBL_EPSILON;

    
    
    
    
    //#ifdef DIAG2
        Jabi[0] = Jab[0] /100;
        Jabi[1] = Jab[1];
        Jabi[2] = Jab[2];
    //#endif

	TRACE(("\nCIECAM02 Reverse conversion:\n"))
	TRACE(("Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]))

    
    /* Correct Chroma for Colorfulness
     * M = C * Fl^0.25
     * When adjusting La, C remains constant but Fl changes, so
     * M shows the effect of different scene lighting.
     *
     * The ratio of (Fl(forward scene) / Fl(reverse scene)) ^0.25 is
     * used to create a La based Chroma effect, keeping M constant */
    
//    Cb = pow((s->color_Flfwd/s->Fl), 0.25);
//    Cb = 1.0;
//    J  = Jab[0] / 100;
//    ja = Jab[1] * Cb;
//    jb = Jab[2] * Cb;
    

	/* Hue angle h */
    h = (180.0/DBL_PI) * atan2(jb, ja);
    h = (h < 0.0) ? h + 360.0 : h;

    
	/* Chroma C */
    C = sqrt(ja * ja + jb * jb);    /* Must be Always +ve, Can be NZ even if J == 0 */

    
//
//    /* Brightness Q */
//    Q = (4.0 / s->C) * sqrt(J) * (s->Aw + 4.0) * pow(s->Fl, 0.25);                          /* eqn (22) */
//
//
//    /* Colorfulness M */
//    M = C * pow(s->Fl, 0.25);                                                               /* eqn (24) */
//
//
//    /* Saturation s (well, actually sat, because s is already in use) */
//    sat = 100 * sqrt(M / Q);                                                                /* eqn (25) */

    
	/* Preliminary Restricted chroma, always +ve and NZ */
	/* (The exact value isn't important because the numerator dominates as a,b aproach 0) */
    rC = C;
    if (rC < DBL_EPSILON)
        rC = DBL_EPSILON;


    /* Undo Helmholtz-Kohlrausch effect */
    J = JJ;
//    #ifndef DISABLE_HHKR
//        if (s->hk && JJ < 1.0) {
//            strength = s->hkscale * HHKR_MUL;
//            kk = strength * rC / 300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0)) / 180.0);   /* rC instead of C */
//            if (kk > 1e-6)                                                                     /* Limit kk to a reasonable range */
//                kk = 1.0/(s->hklimit + 1.0/kk);
//            J = (JJ - kk)/(1.0 - kk);
//            if (J < 0.0)
//                J = JJ - kk;
//            TRACE(("J = %f from JJ = %f, kk = %f\n",J,JJ,kk))
//        }
//    #endif


	/* Achromatic response A */
    #ifndef SYMETRICJ		                                                                /* normally skipped */
        if (J >= s->jlimit) {
            A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
        } else {	                                                                        /* In the straight line segment */
            A = s->lA/s->jlimit * J;
            TRACE(("Undo Achromatic response straight line\n"))
        }
    #else			                                                                        /* normally executed */
        if (J >= 0.0) {
            A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
        } else {	                                                                        /* In the straight line segment */
            A = -pow(-J, 1.0/(s->C * s->z)) * s->Aw;
            TRACE(("Undo symmetric Achromatic response\n"))
        }
    #endif

    
	/* Preliminary Achromatic response ttA */
	ttA = (A / s->Nbb) + 0.305;

	
    /* cJ, the constrained J (positive, non-zero) */
    if (A > 0.0) {
		cJ = pow(A/s->Aw, s->C * s->z);
		if (cJ < s->ssmincj)
			cJ = s->ssmincj;
	} else {
		cJ = s->ssmincj;
	}
	TRACE(("C = %f, A = %f from J = %f, cJ = %f\n", C, A, J, cJ))

    
	/* Eccentricity factor e */
	e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));

    
	/* a, b scaling components k1, k2, k3 */
	k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rC, 1.0/9.0);
	k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;
	k3 = s->dcomp[1] * ja + s->dcomp[2] * jb;
	TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, (k1 - k3)/k2))

    
    /* Limit k3 values */
    #ifdef ENABLE_DDL                                                                       /* normally executed */
        /* Limit ratio of k3 to k1 to stop zero or -ve ss */
        if (k3 > (k1 * s->ddulimit)) {
            k3 = k1 * s->ddulimit;
            TRACE(("k3 limited to %f due to k3:k1 ratio, ss = %f\n",k3,(k1 - k3)/k2))
        }
        /* See if there is going to be a problem in fwd */
        if (k3 < -k1 * s->ddllimit/(1.0 - s->ddllimit)) {
            /* Adjust ss to allow for fwd limitd computation */
            k3 = -k1 * s->ddllimit/(1.0 - s->ddllimit);
            TRACE(("k3 set to %f to allow for fk3:fk2 fwd limit\n",k3))
        }
    #endif


    /* a, b scaling factor ss */
    #ifdef DISABLE_TTD                                                                      /* normally skipped */
        ss = k1/k2;
    #else                                                                                   /* normally executed */
        ss = (k1 - k3)/k2;
    #endif


    /* Limit a, b scaling factor ss */
    #ifdef ENABLE_SS                                                                        /* normally skipped */
        if (ss < SSLLIMIT)
            ss = SSLLIMIT;
        else if (ss > SSULIMIT)
            ss = SSULIMIT;
    #endif
    
	/* Preliminary a and b */
//    double eff = 100 + (ss - 100) * 0.5;
//    ss = 80;
//    a = ja / eff;
//    b = jb / eff;
    a = ja / ss;
    b = jb / ss;
    rS = sqrt (a * a + b * b);
	TRACE(("ss = %f, ttA = %f, a = %f, b = %f\n",ss,ttA,a,b))


    /* Post-adapted cone response sample values rgba[]
     * Solve for post-adapted cone response of sample
	 * using inverse matrix on ttA, a, b */
	rgba[0] = (20.0/61.0) * ttA
	        + ((41.0 * 11.0)/(61.0 * 23.0)) * a
	        + ((288.0 * 1.0)/(61.0 * 23.0)) * b;
	rgba[1] = (20.0/61.0) * ttA
	        - ((81.0 * 11.0)/(61.0 * 23.0)) * a
	        - ((261.0 * 1.0)/(61.0 * 23.0)) * b;
	rgba[2] = (20.0/61.0) * ttA
	        - ((20.0 * 11.0)/(61.0 * 23.0)) * a
	        - ((20.0 * 315.0)/(61.0 * 23.0)) * b;
	TRACE(("rgba %f %f %f\n",rgba[0], rgba[1], rgba[2]))

    
    #ifdef DISABLE_MATRIX                                                                       /* normally skipped */
        XYZ[0] = rgba[0];
        XYZ[1] = rgba[1];
        XYZ[2] = rgba[2];
    #else                                                                                       /* normally executed */

    
        /* HPE cone response sample values rgbp[] */
        #ifdef DISABLE_NONLIN                                                                   /* normally skipped */
            rgbp[0] = 27.13/400.0 * rgba[0];
            rgbp[1] = 27.13/400.0 * rgba[1];
            rgbp[2] = 27.13/400.0 * rgba[2];
        #else	                                                                                /* normally executed */
            for (i = 0; i < 3; i++) {
                if (rgba[i] < s->nldxval) {
                    rgbp[i] = s->nldlimit + (rgba[i] - s->nldxval)/s->nldxslope;
                } else if (rgba[i] <= s->nluxval) {
                    tt = rgba[i] - 0.1;
                    rgbp[i] = pow((27.13 * tt)/(400.0 - tt), 1.0/0.42)/s->Fl;
                } else {
                    rgbp[i] = s->nlulimit + (rgba[i] - s->nluxval)/s->nluxslope;
                }
            }
        #endif
        TRACE(("rgbp %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))

    
        /* Limit maximum blue angle */
        #ifdef ENABLE_BLUE_ANGLE_FIX                                                            /* normally executed */
            ba_ss = rgbp[0] + rgbp[1] + rgbp[2];
            if (ba_ss < 1e-9)
                ba_ss = 0.0;
            else {
                ba_ss = (rgbp[2]/ba_ss - 1.0/3.0) * 3.0/2.0;
                if (ba_ss > 0.0)
                    ba_ss = BLUE_BL_MAX * pow(ba_ss, BLUE_BL_POW);
            }
            if (ba_ss < 0.0)
                ba_ss = 0.0;
            else if (ba_ss > 1.0)
                ba_ss = 1.0;
            tt = 0.5 * (rgbp[0] + rgbp[1]);
            rgbp[0] = (rgbp[0] - ba_ss * tt)/(1.0 - ba_ss);
            rgbp[1] = (rgbp[1] - ba_ss * tt)/(1.0 - ba_ss);
            TRACE(("rgbp after blue fix ba_ss %f = %f %f %f\n",ba_ss, rgbp[0], rgbp[1], rgbp[2]))
        #endif


        /* Undo soft limiting */
        #ifdef ENABLE_DECOMPR
            {
                double tt;			/* Temporary */
                double wrgb[3];		/* White target */

                /* Make white target white point with same Y value */
                tt = rgbp[0] * s->icc[1][0] + rgbp[1] * s->icc[1][1] + rgbp[2] * s->icc[1][2];
                tt = tt > BC_WHMINY ? tt : BC_WHMINY;	                                       /* Limit to minimum Y */
                icmScale3(wrgb, s->rgbpW, tt/s->Wxyz[1]);	                                   /* White target at same Y */
                TRACE(("wrgb %f %f %f\n", wrgb[0], wrgb[1], wrgb[2]))

                /* Un-limit b,g,r in turn */
                for (i = 2; i >= 0; i--) {
                    double cv;		/* Compression value */
                    double ctv;		/* Compression target value */
                    double cd;		/* Compression displacement needed */
                    double cvec[3];	/* Normalized correction vector */
                    double isec[3];	/* Intersection with plane */
                    double offs;	/* Offset of point from orgin*/
                    double range;	/* Threshold to start compression */
                    double asym;	/* Compression asymtope */

                    /* Compute compression direction as vector towards white */
                    /* (We did try correcting in a blend of limit plane normal and white, */
                    /*  but compressing towards white seems to be the best.) */
                    icmSub3(cvec, wrgb, rgbp);					                                /* Direction of white target */
                    TRACE(("ch %d, rgbp %f %f %f\n", i, rgbp[0], rgbp[1], rgbp[2]))
                    TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))

                    if (cvec[i] < 1e-9) {		                                                /* compression direction can't correct this coord */
                        TRACE(("Intersection with limit plane failed\n"))
                        continue;
                    }

                    /* Scale compression vector to make it move a unit in normal direction */
                    icmScale3(cvec, cvec, 1.0/cvec[i]);		                                    /* Normalized vector to white */
                    TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))

                    /* Compute intersection of correction direction with this limit plane */
                    /* (This corresponds with finding displacement of rgbp by cvec */
                    /*  such that the current coord value = 0) */
                    icmScale3(isec, cvec, -rgbp[i]);		                                    /* (since cvec[i] == 1.0) */
                    icmAdd3(isec, isec, rgbp);
                    TRACE(("isec %f %f %f\n", isec[0], isec[1], isec[2]))

                    /* Compute distance from intersection to origin */
                    offs = pow(icmNorm3(isec), 0.85);
                    range = s->crange[i] * offs;	                                            /* Scale range by distance to origin */
                    if (range > BC_MAXRANGE)		                                            /* so that it tapers down as we approach it */
                        range = BC_MAXRANGE;		                                            /* and limit maximum */

                    /* Aiming above plane when far from origin, */
                    /* but below plane at the origin, so that black isn't affected. */
                    asym = range - 0.2 * (range + (0.01 * s->crange[i]));
                    ctv = cv = rgbp[i];		                                                    /* Distance above/below limit plane */
                    TRACE(("ch %d, offs %f, range %f asym %f, cv %f\n",i, offs,range,asym,cv))

                    if (ctv < (range - 1e-12)) {		                                        /* Need to expand */
                        if (ctv <= asym) {
                            cd = BC_LIMIT;
                            TRACE(("ctv %f < asym %f\n",ctv,asym))
                        } else {
                            double aa, bb;
                            aa = 1.0/(range - ctv);
                            bb = 1.0/(range - asym);
                            if (aa > (bb + 1e-12))
                                cv = range - 1.0/(aa - bb);
                            cd = ctv - cv;				                                        /* Displacement needed */
                        }
                        if (cd > BC_LIMIT)
                            cd = BC_LIMIT;
                        TRACE(("ch %d cd = %f, scaled cd %f\n",i,cd,cd))

                        if (cd > 1e-9) {
                            icmScale3(cvec, cvec, -cd);			                                /* Compression vector */
                            icmAdd3(rgbp, rgbp, cvec);			                                /* Compress by displacement */
                            TRACE(("rgbp after decomp. = %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
                        }
                    }
                }
            }
        #endif

    
        /* Chromaticaly transformed sample value */
        /* Spectrally sharpened cone responses */
        /* XYZ values */
        icmMulBy3x3(xyz, s->icc, rgbp);
        TRACE(("XYZ = %f %f %f\n",xyz[0], xyz[1], xyz[2]))

    
        /* Subtract flare */
        XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
        XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
        XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);
        TRACE(("XYZ after flare = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))
        TRACE(("\n"))
    
    #endif  /* end of scope for DISABLE_MATRIX */

    
    /* Brightness Q */
    Q = (4.0 / s->C) * sqrt(J) * (s->Aw + 4.0) * pow(s->Fl, 0.25);                          /* eqn (22) */
    
    
    /* Colorfulness M */
    M = C * pow(s->Fl, 0.25);                                                               /* eqn (24) */
    
    
    /* Saturation s (well, actually sat, because s is already in use) */
    sat = 100 * sqrt(M / Q);                                                                /* eqn (25) */
    
  
//    printf("J, ja, jb, C, M, h, srcFl, dstFl, hk, hkscale %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %2d %10.6f\n", J, ja, jb, C, M, h, srcFl, dstFl, s->hk, s->hkscale);
    
    /* Transfer internal sample values */
    for (int i = 0; i < 3; i++) {
        //      s->color_rgb[i]  = rgb[i];
        //      s->color_rgbc[i] = rgbc[i];
        s->color_rgbp[i] = rgbp[i];
        s->color_rgba[i] = rgba[i];
    }
    s->color_a   = a;
    s->color_b   = b;
    s->color_rS  = rS;
    s->color_h   = h;
    s->color_e   = e;
    s->color_ttd = 0; //ttd;
    s->color_t   = 0; //t;
    s->color_H   = 0; //H;
    s->color_A   = A;
    s->color_J   = J * 100;
    s->color_JJ  = JJ* 100;
    s->color_Q   = Q;
    s->color_C   = C;
    s->color_M   = M;
    s->color_s   = sat;
    s->color_ss  = ss;
    s->color_bass = ba_ss;
    s->color_ja  = ja;
    s->color_jb  = jb;


    
    /* Protect from small (+ve or -ve) values resulting from rounding */
    XYZ[0] = (XYZ[0] < DBL_EPSILON) ? 0.0 : XYZ[0];
    XYZ[1] = (XYZ[1] < DBL_EPSILON) ? 0.0 : XYZ[1];
    XYZ[2] = (XYZ[2] < DBL_EPSILON) ? 0.0 : XYZ[2];
    
    
    /* Clip when beyond White */
    XYZ[0] = (XYZ[0] > s->Wxyz[0]) ? s->Wxyz[0] : XYZ[0];
    XYZ[1] = (XYZ[1] > s->Wxyz[1]) ? s->Wxyz[1] : XYZ[1];
    XYZ[2] = (XYZ[2] > s->Wxyz[2]) ? s->Wxyz[2] : XYZ[2];

    
//    /* Correct XYZ for dynamic range of output
//     *
//     * On Epson SC P800, black is below Lab = 3
//     *
//     * So, subtract 3 from L scale and compress in remaining 100-3= 97 L levels
//     * or the equivalent reflectance of RGB 12,12,12 in ProPhoto: y = 0.004081
//     * and convert back to XYZ */
//    if (C < 120) {
//        double Lab[3] = { 0.0, 0.0, 0.0 };
//        double Lbefore, y, fy;
//
//        XYZ2Lab(XYZ, Lab);
//        Lbefore = Lab[0];
//
//        if (XYZ[1] > 0.004081)
//        {
//            y = (XYZ[1] - 0.004081) / (1.0 - 0.004081);
//            if (y > 0.008856451586)
//                fy = pow(y,1.0/3.0);
//            else
//                fy = 7.787036979 * y + 16.0/116.0;
//
//            Lab[0] = 116.0 * fy - 16.0;
//
//        } else
//            Lab[0] = 0;
//
//        Lab2XYZ(Lab, XYZ);
////        printf("J, ja, jb, C, M, h, XYZ[0], XYZ[1], XYZ[2], L*before, L*, a*, b* %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", J, ja, jb, C, M, h, XYZ[0], XYZ[1], XYZ[2], Lbefore, Lab[0], Lab[1], Lab[2]);
//    }
    
    if (printCAM == 1) {
    #ifdef DIAG2
        printf("Jab -> XYZ, cam02:\n");
        printf("Jab = %f %f %f\n", Jabi[0], Jabi[1], Jabi[2]);
        printf("Chroma C = %f\n", C);
        printf("Preliminary Saturation ss = %f\n", ss);
        printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ * 100);
        printf("Achromatic response A = %f\n", A);
        printf("Eccentricity factor e = %f\n", e);
        printf("Hue angle h = %f\n", h);
        printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
        printf("Hunundeft-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
        printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
        printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
        printf("\n");
    #endif
    }
    
	return 0;
}

    /* ------------------------------------------------------------------------------------------------------------------*/


static int QMh_to_XYZ(
                      struct _cam02 *s,
                      double XYZ[3],
                      double QMh[3]
                      ) {
    int i;
    double xyz[3], QMhi[3], rgbp[3], rgba[3];
    double a, b, ja, jb, J, JJ, C, rS, h, e, A, ss;
    double tt, cJ, ttA;
    double k1, k2, k3;                                                                      /* (k1 & k3 are different from the fwd k1 & k3) */
    double Q, M;
    double ba_ss;
    
    //#ifdef DIAG2
    QMhi[0] = QMh[0];
    QMhi[1] = QMh[1];
    QMhi[2] = QMh[2];
    //#endif
    
    TRACE(("\nCIECAM02  QMh Reverse conversion:\n"))
    TRACE(("QMh %f %f %f\n",QMh[0], QMh[1], QMh[2]))
    
    Q = QMh[0] ;                                                                         /* note: now J rather than JJ */
    M = QMh[1] ;
    h = QMh[2] ;
    
    
    /* Hue angle h */
//    h = (180.0/DBL_PI) * atan2(jb, ja);
//    h = (h < 0.0) ? h + 360.0 : h;
    
    /* Brightness Q */
//    Q = (4.0 / s->C) * sqrt(J) * (s->Aw + 4.0) * pow(s->Fl, 0.25);                          /* eqn (22) */
//
//    Q / (s->Aw + 4.0) * pow(s->Fl, 0.25) = (4.0 / s->C) * sqrt(J)
//
//    s->C * Q / (s->Aw + 4.0) * pow(s->Fl, 0.25) = 4 sqrt(J)
    
    double xx1 = s->C * Q;
    double xx2 = s->Aw + 4;
    double xx3 = pow(s->Fl, 0.25);
    J = 6.25 / 100 * pow(xx1 / (xx2 * xx3), 2);
    
    
    /* Chroma C */
//    C = sqrt(ja * ja + jb * jb);    /* Must be Always +ve, Can be NZ even if J == 0 */
    C = M / pow(s->Fl, 0.25);
    
//
//    /* Preliminary Restricted chroma, always +ve and NZ */
//    /* (The exact value isn't important because the numerator dominates as a,b aproach 0) */
//    rC = C;
//    if (rC < DBL_EPSILON)
//        rC = DBL_EPSILON;
//
//    //J = JJ;
//
//    //    /* Undo Helmholtz-Kohlrausch effect */
//    //    JJ = J;
//    //    #ifndef DISABLE_HHKR
//    //        if (s->hk && J < 1.0) {
//    //            double strength = s->hkscale * HHKR_MUL;
//    //            double kk = strength * rC / 300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0)) / 180.0);   /* rC instead of C */
//    //            if (kk > 1e-6)                                                                     /* Limit kk to a reasonable range */
//    //                kk = 1.0/(s->hklimit + 1.0/kk);
//    //            J = (JJ - kk)/(1.0 - kk);
//    //            if (J < 0.0)
//    //                J = JJ - kk;
//    //            TRACE(("J = %f from JJ = %f, kk = %f\n",J,JJ,kk))
//    //        }
//    //    #endif
//
//
//    /* Achromatic response A */
//#ifndef SYMETRICJ                                                                        /* normally skipped */
//    if (J >= s->jlimit) {
//        A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
//    } else {                                                                            /* In the straight line segment */
//        A = s->lA/s->jlimit * J;
//        TRACE(("Undo Achromatic response straight line\n"))
//    }
//#else                                                                                    /* normally executed */
//    if (J >= 0.0) {
//        A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
//    } else {                                                                            /* In the straight line segment */
//        A = -pow(-J, 1.0/(s->C * s->z)) * s->Aw;
//        TRACE(("Undo symmetric Achromatic response\n"))
//    }
//#endif
//
//
//    /* Preliminary Achromatic response ttA */
//    ttA = (A / s->Nbb) + 0.305;
//
//
//    /* cJ, the constrained J (positive, non-zero) */
//    if (A > 0.0) {
//        cJ = pow(A/s->Aw, s->C * s->z);
//        if (cJ < s->ssmincj)
//            cJ = s->ssmincj;
//    } else {
//        cJ = s->ssmincj;
//    }
//    TRACE(("C = %f, A = %f from J = %f, cJ = %f\n", C, A, J, cJ))
//
//
//    /* Eccentricity factor e */
//    e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));
//
//
//    /* a, b scaling components k1, k2, k3 */
//    k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rC, 1.0/9.0);
//    k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;
//    k3 = s->dcomp[1] * ja + s->dcomp[2] * jb;
//    TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, (k1 - k3)/k2))
//
//
//    /* Limit k3 values */
//#ifdef ENABLE_DDL                                                                       /* normally executed */
//    /* Limit ratio of k3 to k1 to stop zero or -ve ss */
//    if (k3 > (k1 * s->ddulimit)) {
//        k3 = k1 * s->ddulimit;
//        TRACE(("k3 limited to %f due to k3:k1 ratio, ss = %f\n",k3,(k1 - k3)/k2))
//    }
//    /* See if there is going to be a problem in fwd */
//    if (k3 < -k1 * s->ddllimit/(1.0 - s->ddllimit)) {
//        /* Adjust ss to allow for fwd limitd computation */
//        k3 = -k1 * s->ddllimit/(1.0 - s->ddllimit);
//        TRACE(("k3 set to %f to allow for fk3:fk2 fwd limit\n",k3))
//    }
//#endif
//
//
//    /* a, b scaling factor ss */
//#ifdef DISABLE_TTD                                                                      /* normally skipped */
//    ss = k1/k2;
//#else                                                                                   /* normally executed */
//    ss = (k1 - k3)/k2;
//#endif
//
//
//    /* Limit a, b scaling factor ss */
//#ifdef ENABLE_SS                                                                        /* normally skipped */
//    if (ss < SSLLIMIT)
//        ss = SSLLIMIT;
//    else if (ss > SSULIMIT)
//        ss = SSULIMIT;
//#endif
//
//    //    ss = s->color_ss;
//
//    /* Preliminary a and b */
//    //    double eff = 100 + (ss - 100) * 0.5;
//    //    ss = 80;
//    //    a = ja / eff;
//    //    b = jb / eff;
//    a = ja / ss;
//    b = jb / ss;
//    TRACE(("ss = %f, ttA = %f, a = %f, b = %f\n",ss,ttA,a,b))
    

    double t, et;
    double p1, p2, p3, p4, p5;
    double hr, ca, cb;
    
    
    t = pow(C / (sqrt(J) * s->nn), 1.0/0.9);
    et = (1.0 / 4.0) * (cos(((h * DBL_PI) / 180.0) + 2.0) + 3.8);
    
    A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
    
    p1 = ((50000.0 / 13.0) * s->Nc * s->Ncb) * et/t;
    p2 = (A/s->Nbb) + 0.305;
    p3 = 21.0 / 20.0;
    
    hr = (h * DBL_PI) / 180.0;
    
    if (fabs(sin(hr)) >= fabs(cos(hr))) {
        p4 = p1 / sin(hr);
        cb = (p2 * (2.0 + p3) * (460.0 / 1403.0)) /
             (p4 + (2.0 + p3) * (220.0 / 1403.0) *
             (cos(hr) / sin(hr)) - (27.0 / 1403.0) +
             p3 * (6300.0 / 1403.0));
        ca = cb * (cos(hr) / sin(hr));
    }
    else {
        p5 = p1 / cos(hr);
        ca = (p2 * (2.0 + p3) * (460.0 / 1403.0)) /
             (p5 + (2.0 + p3) * (220.0 / 1403.0) -
             ((27.0 / 1403.0) - p3 * (6300.0 / 1403.0)) *
             (sin(hr) / cos(hr)));
        cb = ca * (sin(hr) / cos(hr));
    }

    a = ca;
    b = cb;
    
    rS = sqrt(a * a + b * b);
    ja = a * C/rS;
    jb = b * C/rS;
    
    ttA = p2;
    
    
    /* Post-adapted cone response sample values rgba[]
     * Solve for post-adapted cone response of sample
     * using inverse matrix on ttA, a, b */
    rgba[0] = (  20.0 /  61.0) * ttA
            + (( 41.0 *  11.0) / (61.0 * 23.0)) * a
            + ((288.0 *   1.0) / (61.0 * 23.0)) * b;
    rgba[1] = (  20.0 /  61.0) * ttA
            - (( 81.0 *  11.0) / (61.0 * 23.0)) * a
            - ((261.0 *   1.0) / (61.0 * 23.0)) * b;
    rgba[2] = (  20.0 /  61.0) * ttA
            - (( 20.0 *  11.0) / (61.0 * 23.0)) * a
            - (( 20.0 * 315.0) / (61.0 * 23.0)) * b;
    TRACE(("rgba %f %f %f\n",rgba[0], rgba[1], rgba[2]))
    
    
    #ifdef DISABLE_MATRIX                                                                   /* normally skipped */
    XYZ[0] = rgba[0];
    XYZ[1] = rgba[1];
    XYZ[2] = rgba[2];

    #else                                                                                   /* normally executed */
    /* HPE cone response sample values rgbp[] */
    #ifdef DISABLE_NONLIN                                                                   /* normally skipped */
        rgbp[0] = 27.13/400.0 * rgba[0];
        rgbp[1] = 27.13/400.0 * rgba[1];
        rgbp[2] = 27.13/400.0 * rgba[2];
    #else                                                                                    /* normally executed */
        for (i = 0; i < 3; i++) {
            if (rgba[i] < s->nldxval) {
                rgbp[i] = s->nldlimit + (rgba[i] - s->nldxval)/s->nldxslope;
            } else if (rgba[i] <= s->nluxval) {
                tt = rgba[i] - 0.1;
                rgbp[i] = pow((27.13 * tt)/(400.0 - tt), 1.0/0.42)/s->Fl;
            } else {
                rgbp[i] = s->nlulimit + (rgba[i] - s->nluxval)/s->nluxslope;
            }
        }
    #endif
    TRACE(("rgbp %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
    
    
    /* Limit maximum blue angle */
    #ifdef ENABLE_BLUE_ANGLE_FIX                                                            /* normally executed */
        ba_ss = rgbp[0] + rgbp[1] + rgbp[2];
        if (ba_ss < 1e-9)
            ba_ss = 0.0;
        else {
            ss = (rgbp[2]/ba_ss - 1.0/3.0) * 3.0/2.0;
            if (ba_ss > 0.0)
                ba_ss = BLUE_BL_MAX * pow(ba_ss, BLUE_BL_POW);
        }
        if (ba_ss < 0.0)
            ba_ss = 0.0;
        else if (ba_ss > 1.0)
            ba_ss = 1.0;
        tt = 0.5 * (rgbp[0] + rgbp[1]);
        rgbp[0] = (rgbp[0] - ba_ss * tt)/(1.0 - ba_ss);
        rgbp[1] = (rgbp[1] - ba_ss * tt)/(1.0 - ba_ss);
        TRACE(("rgbp after blue fix ba_ss %f = %f %f %f\n", ba_ss, rgbp[0], rgbp[1], rgbp[2]))
    #endif
    
    
    /* Undo soft limiting */
    #ifdef ENABLE_DECOMPR
        {
            double tt;            /* Temporary */
            double wrgb[3];        /* White target */
            
            /* Make white target white point with same Y value */
            tt = rgbp[0] * s->icc[1][0] + rgbp[1] * s->icc[1][1] + rgbp[2] * s->icc[1][2];
            tt = tt > BC_WHMINY ? tt : BC_WHMINY;                                           /* Limit to minimum Y */
            icmScale3(wrgb, s->rgbpW, tt/s->Wxyz[1]);                                       /* White target at same Y */
            TRACE(("wrgb %f %f %f\n", wrgb[0], wrgb[1], wrgb[2]))
            
            /* Un-limit b,g,r in turn */
            for (i = 2; i >= 0; i--) {
                double cv;        /* Compression value */
                double ctv;        /* Compression target value */
                double cd;        /* Compression displacement needed */
                double cvec[3];    /* Normalized correction vector */
                double isec[3];    /* Intersection with plane */
                double offs;    /* Offset of point from orgin*/
                double range;    /* Threshold to start compression */
                double asym;    /* Compression asymtope */
                
                /* Compute compression direction as vector towards white */
                /* (We did try correcting in a blend of limit plane normal and white, */
                /*  but compressing towards white seems to be the best.) */
                icmSub3(cvec, wrgb, rgbp);                                                    /* Direction of white target */
                TRACE(("ch %d, rgbp %f %f %f\n", i, rgbp[0], rgbp[1], rgbp[2]))
                TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))
                
                if (cvec[i] < 1e-9) {                                                        /* compression direction can't correct this coord */
                    TRACE(("Intersection with limit plane failed\n"))
                    continue;
                }
                
                /* Scale compression vector to make it move a unit in normal direction */
                icmScale3(cvec, cvec, 1.0/cvec[i]);                                            /* Normalized vector to white */
                TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))
                
                /* Compute intersection of correction direction with this limit plane */
                /* (This corresponds with finding displacement of rgbp by cvec */
                /*  such that the current coord value = 0) */
                icmScale3(isec, cvec, -rgbp[i]);                                            /* (since cvec[i] == 1.0) */
                icmAdd3(isec, isec, rgbp);
                TRACE(("isec %f %f %f\n", isec[0], isec[1], isec[2]))
                
                /* Compute distance from intersection to origin */
                offs = pow(icmNorm3(isec), 0.85);
                range = s->crange[i] * offs;                                                /* Scale range by distance to origin */
                if (range > BC_MAXRANGE)                                                    /* so that it tapers down as we approach it */
                    range = BC_MAXRANGE;                                                    /* and limit maximum */
                
                /* Aiming above plane when far from origin, */
                /* but below plane at the origin, so that black isn't affected. */
                asym = range - 0.2 * (range + (0.01 * s->crange[i]));
                ctv = cv = rgbp[i];                                                            /* Distance above/below limit plane */
                TRACE(("ch %d, offs %f, range %f asym %f, cv %f\n",i, offs,range,asym,cv))
                
                if (ctv < (range - 1e-12)) {                                                /* Need to expand */
                    if (ctv <= asym) {
                        cd = BC_LIMIT;
                        TRACE(("ctv %f < asym %f\n",ctv,asym))
                    } else {
                        double aa, bb;
                        aa = 1.0/(range - ctv);
                        bb = 1.0/(range - asym);
                        if (aa > (bb + 1e-12))
                            cv = range - 1.0/(aa - bb);
                        cd = ctv - cv;                                                        /* Displacement needed */
                    }
                    if (cd > BC_LIMIT)
                        cd = BC_LIMIT;
                    TRACE(("ch %d cd = %f, scaled cd %f\n",i,cd,cd))
                    
                    if (cd > 1e-9) {
                        icmScale3(cvec, cvec, -cd);                                            /* Compression vector */
                        icmAdd3(rgbp, rgbp, cvec);                                            /* Compress by displacement */
                        TRACE(("rgbp after decomp. = %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
                    }
                }
            }
        }
    #endif
    
    
    /* Chromaticaly transformed sample value */
    /* Spectrally sharpened cone responses */
    /* XYZ values */
    icmMulBy3x3(xyz, s->icc, rgbp);
    TRACE(("XYZ = %f %f %f\n",xyz[0], xyz[1], xyz[2]))
    
    
    /* Subtract flare */
    XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
    XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
    XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);
    TRACE(("XYZ after flare = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))
    TRACE(("\n"))
    
#endif  /* end of scope for DISABLE_MATRIX */
    
    /* Transfer internal sample values */
    for (int i = 0; i < 3; i++) {
        //      s->color_rgb[i]  = rgb[i];
        //      s->color_rgbc[i] = rgbc[i];
        s->color_rgbp[i] = rgbp[i];
        s->color_rgba[i] = rgba[i];
    }
    s->color_a   = a;
    s->color_b   = b;
    s->color_rS  = rS;
    s->color_h   = h;
    s->color_e   = e;
    s->color_ttd = 0; //ttd;
    s->color_t   = 0; //t;
    s->color_H   = 0; //H;
    s->color_A   = A;
    s->color_J   = J * 100;
    s->color_JJ  = JJ;
    s->color_Q   = Q;
    s->color_C   = C;
    s->color_M   = M;
    s->color_s   = 0; //sat;
    s->color_bass  = ba_ss;
    s->color_ja  = ja;
    s->color_jb  = jb;
    
    
    if (printCAM == 1) {
#ifdef DIAG2
        printf("QMh -> XYZ, cam02:\n");
        printf("QMh = %f %f %f\n", QMhi[0], QMhi[1], QMhi[2]);
        printf("Chroma C = %f\n", C);
        printf("Preliminary Saturation ss = %f\n", ss);
        printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ * 100);
        printf("Achromatic response A = %f\n", A);
        printf("Eccentricity factor e = %f\n", e);
        printf("Hue angle h = %f\n", h);
        printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
        printf("Hunundeft-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
        printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
        printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
        printf("\n");
#endif
    }
    
    return 0;
}



/* ------------------------------------------------------------------------------------------------------------------*/


static int JMh_to_XYZ(
                      struct _cam02 *s,
                      double XYZ[3],
                      double JMh[3]
                      ) {
    int i;
    double xyz[3], JMhi[3], rgbp[3], rgba[3];
    double a, b, ja, jb, J, JJ, C, rC, h, hr, e, A, ss;
    double tt, cJ, ttA;
    double k1, k2, k3;                                                                      /* (k1 & k3 are different from the fwd k1 & k3) */
    double Q, M, sat;
    double ba_ss;
    
    //#ifdef DIAG2
    JMhi[0] = JMh[0];
    JMhi[1] = JMh[1];
    JMhi[2] = JMh[2];
    //#endif
    
    TRACE(("\nCIECAM02  JMh Reverse conversion:\n"))
    TRACE(("JMh %f %f %f\n",JMh[0], JMh[1], JMh[2]))
    
    J = JMh[0] /100 ;
    M = JMh[1] ;
    h = JMh[2] ;

    
    /* Brightness Q */
    Q = (4.0 / s->C) * sqrt(J) * (s->Aw + 4.0) * pow(s->Fl, 0.25);                          /* eqn (22) */
    
    
    /* Chroma C */
    C = M / pow(s->Fl, 0.25);
    

    /* Preliminary Restricted chroma, always +ve and NZ */
    /* (The exact value isn't important because the numerator dominates as a,b aproach 0) */
    rC = C;
    if (rC < DBL_EPSILON)
        rC = DBL_EPSILON;
    //
    //    //J = JJ;
    //
    //    //    /* Undo Helmholtz-Kohlrausch effect */
    //    //    JJ = J;
    //    //    #ifndef DISABLE_HHKR
    //    //        if (s->hk && J < 1.0) {
    //    //            double strength = s->hkscale * HHKR_MUL;
    //    //            double kk = strength * rC / 300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0)) / 180.0);   /* rC instead of C */
    //    //            if (kk > 1e-6)                                                                     /* Limit kk to a reasonable range */
    //    //                kk = 1.0/(s->hklimit + 1.0/kk);
    //    //            J = (JJ - kk)/(1.0 - kk);
    //    //            if (J < 0.0)
    //    //                J = JJ - kk;
    //    //            TRACE(("J = %f from JJ = %f, kk = %f\n",J,JJ,kk))
    //    //        }
    //    //    #endif
    //
    //
    
    hr = (h * DBL_PI) / 180.0;
    ja = rC * cos(hr);
    jb = rC * sin(hr);
    
    
    /* Saturation s (well, actually sat, because s is already in use) */
    sat = 100 * sqrt(M / Q);                                                                /* eqn (25) */
    
        /* Achromatic response A */
    #ifndef SYMETRICJ                                                                        /* normally skipped */
        if (J >= s->jlimit) {
            A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
        } else {                                                                            /* In the straight line segment */
            A = s->lA/s->jlimit * J;
            TRACE(("Undo Achromatic response straight line\n"))
        }
    #else                                                                                    /* normally executed */
        if (J >= 0.0) {
            A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
        } else {                                                                            /* In the straight line segment */
            A = -pow(-J, 1.0/(s->C * s->z)) * s->Aw;
            TRACE(("Undo symmetric Achromatic response\n"))
        }
    #endif
    
    
        /* Preliminary Achromatic response ttA */
        ttA = (A / s->Nbb) + 0.305;
    
    
        /* cJ, the constrained J (positive, non-zero) */
        if (A > 0.0) {
            cJ = pow(A/s->Aw, s->C * s->z);
            if (cJ < s->ssmincj)
                cJ = s->ssmincj;
        } else {
            cJ = s->ssmincj;
        }
        TRACE(("C = %f, A = %f from J = %f, cJ = %f\n", C, A, J, cJ))
    
    
        /* Eccentricity factor e */
        e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));
    
    
        /* a, b scaling components k1, k2, k3 */
        k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rC, 1.0/9.0);
        k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;
        k3 = s->dcomp[1] * ja + s->dcomp[2] * jb;
        TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, (k1 - k3)/k2))
    
    
        /* Limit k3 values */
    #ifdef ENABLE_DDL                                                                       /* normally executed */
        /* Limit ratio of k3 to k1 to stop zero or -ve ss */
        if (k3 > (k1 * s->ddulimit)) {
            k3 = k1 * s->ddulimit;
            TRACE(("k3 limited to %f due to k3:k1 ratio, ss = %f\n",k3,(k1 - k3)/k2))
        }
        /* See if there is going to be a problem in fwd */
        if (k3 < -k1 * s->ddllimit/(1.0 - s->ddllimit)) {
            /* Adjust ss to allow for fwd limitd computation */
            k3 = -k1 * s->ddllimit/(1.0 - s->ddllimit);
            TRACE(("k3 set to %f to allow for fk3:fk2 fwd limit\n",k3))
        }
    #endif
    
    
        /* a, b scaling factor ss */
    #ifdef DISABLE_TTD                                                                      /* normally skipped */
        ss = k1/k2;
    #else                                                                                   /* normally executed */
        ss = (k1 - k3)/k2;
    #endif
    
    
        /* Limit a, b scaling factor ss */
    #ifdef ENABLE_SS                                                                        /* normally skipped */
        if (ss < SSLLIMIT)
            ss = SSLLIMIT;
        else if (ss > SSULIMIT)
            ss = SSULIMIT;
    #endif
    
        //    ss = s->color_ss;
    
        /* Preliminary a and b */
        //    double eff = 100 + (ss - 100) * 0.5;
        //    ss = 80;
        //    a = ja / eff;
        //    b = jb / eff;
        a = ja / ss;
        b = jb / ss;
        TRACE(("ss = %f, ttA = %f, a = %f, b = %f\n",ss,ttA,a,b))
    
    
    double t, et;
    double p1, p2, p3, p4, p5;
    double ca, cb, rS;
    
    t = pow(C / (sqrt(J) * s->nn), 1.0/0.9);
    et = (1.0 / 4.0) * (cos(((h * DBL_PI) / 180.0) + 2.0) + 3.8);
    
    A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
    
    p1 = ((50000.0 / 13.0) * s->Nc * s->Ncb) * et/t;
    p2 = (A/s->Nbb) + 0.305;
    p3 = 21.0 / 20.0;
    
    hr = (h * DBL_PI) / 180.0;
    
    if (fabs(sin(hr)) >= fabs(cos(hr))) {
        p4 = p1 / sin(hr);
        cb = (p2 * (2.0 + p3) * (460.0 / 1403.0)) /
        (p4 + (2.0 + p3) * (220.0 / 1403.0) *
         (cos(hr) / sin(hr)) - (27.0 / 1403.0) +
         p3 * (6300.0 / 1403.0));
        ca = cb * (cos(hr) / sin(hr));
    }
    else {
        p5 = p1 / cos(hr);
        ca = (p2 * (2.0 + p3) * (460.0 / 1403.0)) /
        (p5 + (2.0 + p3) * (220.0 / 1403.0) -
         ((27.0 / 1403.0) - p3 * (6300.0 / 1403.0)) *
         (sin(hr) / cos(hr)));
        cb = ca * (sin(hr) / cos(hr));
    }
    
//    a = ca;
//    b = cb;
    
    rS = sqrt(a * a + b * b);
//    ja = a * C/rS;
//    jb = b * C/rS;
    
//    ja = a * ss;
//    jb = b * ss;
    
//    ttA = p2;
    
    a = ja / ss;
    b = jb / ss;
    
    
    /* Post-adapted cone response sample values rgba[]
     * Solve for post-adapted cone response of sample
     * using inverse matrix on ttA, a, b */
    rgba[0] = (  20.0 /  61.0) * ttA
    + (( 41.0 *  11.0) / (61.0 * 23.0)) * a
    + ((288.0 *   1.0) / (61.0 * 23.0)) * b;
    rgba[1] = (  20.0 /  61.0) * ttA
    - (( 81.0 *  11.0) / (61.0 * 23.0)) * a
    - ((261.0 *   1.0) / (61.0 * 23.0)) * b;
    rgba[2] = (  20.0 /  61.0) * ttA
    - (( 20.0 *  11.0) / (61.0 * 23.0)) * a
    - (( 20.0 * 315.0) / (61.0 * 23.0)) * b;
    TRACE(("rgba %f %f %f\n",rgba[0], rgba[1], rgba[2]))
    
    
#ifdef DISABLE_MATRIX                                                                   /* normally skipped */
    XYZ[0] = rgba[0];
    XYZ[1] = rgba[1];
    XYZ[2] = rgba[2];
    
#else                                                                                   /* normally executed */
    /* HPE cone response sample values rgbp[] */
#ifdef DISABLE_NONLIN                                                                   /* normally skipped */
    rgbp[0] = 27.13/400.0 * rgba[0];
    rgbp[1] = 27.13/400.0 * rgba[1];
    rgbp[2] = 27.13/400.0 * rgba[2];
#else                                                                                    /* normally executed */
    for (i = 0; i < 3; i++) {
        if (rgba[i] < s->nldxval) {
            rgbp[i] = s->nldlimit + (rgba[i] - s->nldxval)/s->nldxslope;
        } else if (rgba[i] <= s->nluxval) {
            tt = rgba[i] - 0.1;
            rgbp[i] = pow((27.13 * tt)/(400.0 - tt), 1.0/0.42)/s->Fl;
        } else {
            rgbp[i] = s->nlulimit + (rgba[i] - s->nluxval)/s->nluxslope;
        }
    }
#endif
    TRACE(("rgbp %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
    
    
    /* Limit maximum blue angle */
#ifdef ENABLE_BLUE_ANGLE_FIX                                                            /* normally executed */
    ba_ss = rgbp[0] + rgbp[1] + rgbp[2];
    if (ba_ss < 1e-9)
        ba_ss = 0.0;
    else {
        ba_ss = (rgbp[2]/ba_ss - 1.0/3.0) * 3.0/2.0;
        if (ba_ss > 0.0)
            ba_ss = BLUE_BL_MAX * pow(ba_ss, BLUE_BL_POW);
    }
    if (ba_ss < 0.0)
        ba_ss = 0.0;
    else if (ba_ss > 1.0)
        ba_ss = 1.0;
    tt = 0.5 * (rgbp[0] + rgbp[1]);
    rgbp[0] = (rgbp[0] - ba_ss * tt)/(1.0 - ba_ss);
    rgbp[1] = (rgbp[1] - ba_ss * tt)/(1.0 - ba_ss);
    TRACE(("rgbp after blue fix ba_ss %f = %f %f %f\n", ba_ss,rgbp[0], rgbp[1], rgbp[2]))
#endif
    
    
    /* Undo soft limiting */
#ifdef ENABLE_DECOMPR
    {
        double tt;            /* Temporary */
        double wrgb[3];        /* White target */
        
        /* Make white target white point with same Y value */
        tt = rgbp[0] * s->icc[1][0] + rgbp[1] * s->icc[1][1] + rgbp[2] * s->icc[1][2];
        tt = tt > BC_WHMINY ? tt : BC_WHMINY;                                           /* Limit to minimum Y */
        icmScale3(wrgb, s->rgbpW, tt/s->Wxyz[1]);                                       /* White target at same Y */
        TRACE(("wrgb %f %f %f\n", wrgb[0], wrgb[1], wrgb[2]))
        
        /* Un-limit b,g,r in turn */
        for (i = 2; i >= 0; i--) {
            double cv;        /* Compression value */
            double ctv;        /* Compression target value */
            double cd;        /* Compression displacement needed */
            double cvec[3];    /* Normalized correction vector */
            double isec[3];    /* Intersection with plane */
            double offs;    /* Offset of point from orgin*/
            double range;    /* Threshold to start compression */
            double asym;    /* Compression asymtope */
            
            /* Compute compression direction as vector towards white */
            /* (We did try correcting in a blend of limit plane normal and white, */
            /*  but compressing towards white seems to be the best.) */
            icmSub3(cvec, wrgb, rgbp);                                                    /* Direction of white target */
            TRACE(("ch %d, rgbp %f %f %f\n", i, rgbp[0], rgbp[1], rgbp[2]))
            TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))
            
            if (cvec[i] < 1e-9) {                                                        /* compression direction can't correct this coord */
                TRACE(("Intersection with limit plane failed\n"))
                continue;
            }
            
            /* Scale compression vector to make it move a unit in normal direction */
            icmScale3(cvec, cvec, 1.0/cvec[i]);                                            /* Normalized vector to white */
            TRACE(("cvec %f %f %f\n", cvec[0], cvec[1], cvec[2]))
            
            /* Compute intersection of correction direction with this limit plane */
            /* (This corresponds with finding displacement of rgbp by cvec */
            /*  such that the current coord value = 0) */
            icmScale3(isec, cvec, -rgbp[i]);                                            /* (since cvec[i] == 1.0) */
            icmAdd3(isec, isec, rgbp);
            TRACE(("isec %f %f %f\n", isec[0], isec[1], isec[2]))
            
            /* Compute distance from intersection to origin */
            offs = pow(icmNorm3(isec), 0.85);
            range = s->crange[i] * offs;                                                /* Scale range by distance to origin */
            if (range > BC_MAXRANGE)                                                    /* so that it tapers down as we approach it */
                range = BC_MAXRANGE;                                                    /* and limit maximum */
            
            /* Aiming above plane when far from origin, */
            /* but below plane at the origin, so that black isn't affected. */
            asym = range - 0.2 * (range + (0.01 * s->crange[i]));
            ctv = cv = rgbp[i];                                                            /* Distance above/below limit plane */
            TRACE(("ch %d, offs %f, range %f asym %f, cv %f\n",i, offs,range,asym,cv))
            
            if (ctv < (range - 1e-12)) {                                                /* Need to expand */
                if (ctv <= asym) {
                    cd = BC_LIMIT;
                    TRACE(("ctv %f < asym %f\n",ctv,asym))
                } else {
                    double aa, bb;
                    aa = 1.0/(range - ctv);
                    bb = 1.0/(range - asym);
                    if (aa > (bb + 1e-12))
                        cv = range - 1.0/(aa - bb);
                    cd = ctv - cv;                                                        /* Displacement needed */
                }
                if (cd > BC_LIMIT)
                    cd = BC_LIMIT;
                TRACE(("ch %d cd = %f, scaled cd %f\n",i,cd,cd))
                
                if (cd > 1e-9) {
                    icmScale3(cvec, cvec, -cd);                                            /* Compression vector */
                    icmAdd3(rgbp, rgbp, cvec);                                            /* Compress by displacement */
                    TRACE(("rgbp after decomp. = %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))
                }
            }
        }
    }
#endif
    
    
    /* Chromaticaly transformed sample value */
    /* Spectrally sharpened cone responses */
    /* XYZ values */
    icmMulBy3x3(xyz, s->icc, rgbp);
    TRACE(("XYZ = %f %f %f\n",xyz[0], xyz[1], xyz[2]))
    
    
    /* Subtract flare */
    XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
    XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
    XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);
    TRACE(("XYZ after flare = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))
    TRACE(("\n"))
    
#endif  /* end of scope for DISABLE_MATRIX */
    
    /* Transfer internal sample values */
    for (int i = 0; i < 3; i++) {
        //      s->color_rgb[i]  = rgb[i];
        //      s->color_rgbc[i] = rgbc[i];
        s->color_rgbp[i] = rgbp[i];
        s->color_rgba[i] = rgba[i];
    }
    s->color_a   = a;
    s->color_b   = b;
    s->color_rS  = rS;
    s->color_h   = h;
    s->color_e   = e;
    s->color_ttd = 0; //ttd;
    s->color_t   = 0; //t;
    s->color_H   = 0; //H;
    s->color_A   = A;
    s->color_J   = J * 100;
    s->color_JJ  = JJ;
    s->color_Q   = Q;
    s->color_C   = C;
    s->color_M   = M;
    s->color_s   = sat;
    s->color_bass  = ba_ss;
    s->color_ja  = ja;
    s->color_jb  = jb;
    
    
    if (printCAM == 1) {
#ifdef DIAG2
        printf("JMh -> XYZ, cam02:\n");
        printf("JMh = %f %f %f\n", JMhi[0], JMhi[1], JMhi[2]);
        printf("Chroma C = %f\n", C);
        printf("Preliminary Saturation ss = %f\n", ss);
        printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ * 100);
        printf("Achromatic response A = %f\n", A);
        printf("Eccentricity factor e = %f\n", e);
        printf("Hue angle h = %f\n", h);
        printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
        printf("Hunundeft-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
        printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
        printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
        printf("\n");
#endif
    }
    
    return 0;
}
