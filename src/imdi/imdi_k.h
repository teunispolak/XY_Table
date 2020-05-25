/* Integer Multi-Dimensional Interpolation */
/* Declarations for all the generated kernel functions */
/* This file is generated by imdi_make */

/* Copyright 2000 - 2007 Graeme W. Gill */
/* All rights reserved. */
/* This material is licensed under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :- */
/* see the License.txt file for licensing details.*/

#include "imdi_k.c"	/* All the kernel code */

struct {
	void (*interp)(imdi *s, void **outp, int ostride, void **inp, int  istride, unsigned int npix);
	void (*gentab)(genspec *g, tabspec *t);
} ktable[891] = {
	{ imdi_k1, imdi_k1_gentab },
	{ imdi_k2, imdi_k2_gentab },
	{ imdi_k3, imdi_k3_gentab },
	{ imdi_k4, imdi_k4_gentab },
	{ imdi_k5, imdi_k5_gentab },
	{ imdi_k6, imdi_k6_gentab },
	{ imdi_k7, imdi_k7_gentab },
	{ imdi_k8, imdi_k8_gentab },
	{ imdi_k9, imdi_k9_gentab },
	{ imdi_k10, imdi_k10_gentab },
	{ imdi_k11, imdi_k11_gentab },
	{ imdi_k12, imdi_k12_gentab },
	{ imdi_k13, imdi_k13_gentab },
	{ imdi_k14, imdi_k14_gentab },
	{ imdi_k15, imdi_k15_gentab },
	{ imdi_k16, imdi_k16_gentab },
	{ imdi_k17, imdi_k17_gentab },
	{ imdi_k18, imdi_k18_gentab },
	{ imdi_k19, imdi_k19_gentab },
	{ imdi_k20, imdi_k20_gentab },
	{ imdi_k21, imdi_k21_gentab },
	{ imdi_k22, imdi_k22_gentab },
	{ imdi_k23, imdi_k23_gentab },
	{ imdi_k24, imdi_k24_gentab },
	{ imdi_k25, imdi_k25_gentab },
	{ imdi_k26, imdi_k26_gentab },
	{ imdi_k27, imdi_k27_gentab },
	{ imdi_k28, imdi_k28_gentab },
	{ imdi_k29, imdi_k29_gentab },
	{ imdi_k30, imdi_k30_gentab },
	{ imdi_k31, imdi_k31_gentab },
	{ imdi_k32, imdi_k32_gentab },
	{ imdi_k33, imdi_k33_gentab },
	{ imdi_k34, imdi_k34_gentab },
	{ imdi_k35, imdi_k35_gentab },
	{ imdi_k36, imdi_k36_gentab },
	{ imdi_k37, imdi_k37_gentab },
	{ imdi_k38, imdi_k38_gentab },
	{ imdi_k39, imdi_k39_gentab },
	{ imdi_k40, imdi_k40_gentab },
	{ imdi_k41, imdi_k41_gentab },
	{ imdi_k42, imdi_k42_gentab },
	{ imdi_k43, imdi_k43_gentab },
	{ imdi_k44, imdi_k44_gentab },
	{ imdi_k45, imdi_k45_gentab },
	{ imdi_k46, imdi_k46_gentab },
	{ imdi_k47, imdi_k47_gentab },
	{ imdi_k48, imdi_k48_gentab },
	{ imdi_k49, imdi_k49_gentab },
	{ imdi_k50, imdi_k50_gentab },
	{ imdi_k51, imdi_k51_gentab },
	{ imdi_k52, imdi_k52_gentab },
	{ imdi_k53, imdi_k53_gentab },
	{ imdi_k54, imdi_k54_gentab },
	{ imdi_k55, imdi_k55_gentab },
	{ imdi_k56, imdi_k56_gentab },
	{ imdi_k57, imdi_k57_gentab },
	{ imdi_k58, imdi_k58_gentab },
	{ imdi_k59, imdi_k59_gentab },
	{ imdi_k60, imdi_k60_gentab },
	{ imdi_k61, imdi_k61_gentab },
	{ imdi_k62, imdi_k62_gentab },
	{ imdi_k63, imdi_k63_gentab },
	{ imdi_k64, imdi_k64_gentab },
	{ imdi_k65, imdi_k65_gentab },
	{ imdi_k66, imdi_k66_gentab },
	{ imdi_k67, imdi_k67_gentab },
	{ imdi_k68, imdi_k68_gentab },
	{ imdi_k69, imdi_k69_gentab },
	{ imdi_k70, imdi_k70_gentab },
	{ imdi_k71, imdi_k71_gentab },
	{ imdi_k72, imdi_k72_gentab },
	{ imdi_k73, imdi_k73_gentab },
	{ imdi_k74, imdi_k74_gentab },
	{ imdi_k75, imdi_k75_gentab },
	{ imdi_k76, imdi_k76_gentab },
	{ imdi_k77, imdi_k77_gentab },
	{ imdi_k78, imdi_k78_gentab },
	{ imdi_k79, imdi_k79_gentab },
	{ imdi_k80, imdi_k80_gentab },
	{ imdi_k81, imdi_k81_gentab },
	{ imdi_k82, imdi_k82_gentab },
	{ imdi_k83, imdi_k83_gentab },
	{ imdi_k84, imdi_k84_gentab },
	{ imdi_k85, imdi_k85_gentab },
	{ imdi_k86, imdi_k86_gentab },
	{ imdi_k87, imdi_k87_gentab },
	{ imdi_k88, imdi_k88_gentab },
	{ imdi_k89, imdi_k89_gentab },
	{ imdi_k90, imdi_k90_gentab },
	{ imdi_k91, imdi_k91_gentab },
	{ imdi_k92, imdi_k92_gentab },
	{ imdi_k93, imdi_k93_gentab },
	{ imdi_k94, imdi_k94_gentab },
	{ imdi_k95, imdi_k95_gentab },
	{ imdi_k96, imdi_k96_gentab },
	{ imdi_k97, imdi_k97_gentab },
	{ imdi_k98, imdi_k98_gentab },
	{ imdi_k99, imdi_k99_gentab },
	{ imdi_k100, imdi_k100_gentab },
	{ imdi_k101, imdi_k101_gentab },
	{ imdi_k102, imdi_k102_gentab },
	{ imdi_k103, imdi_k103_gentab },
	{ imdi_k104, imdi_k104_gentab },
	{ imdi_k105, imdi_k105_gentab },
	{ imdi_k106, imdi_k106_gentab },
	{ imdi_k107, imdi_k107_gentab },
	{ imdi_k108, imdi_k108_gentab },
	{ imdi_k109, imdi_k109_gentab },
	{ imdi_k110, imdi_k110_gentab },
	{ imdi_k111, imdi_k111_gentab },
	{ imdi_k112, imdi_k112_gentab },
	{ imdi_k113, imdi_k113_gentab },
	{ imdi_k114, imdi_k114_gentab },
	{ imdi_k115, imdi_k115_gentab },
	{ imdi_k116, imdi_k116_gentab },
	{ imdi_k117, imdi_k117_gentab },
	{ imdi_k118, imdi_k118_gentab },
	{ imdi_k119, imdi_k119_gentab },
	{ imdi_k120, imdi_k120_gentab },
	{ imdi_k121, imdi_k121_gentab },
	{ imdi_k122, imdi_k122_gentab },
	{ imdi_k123, imdi_k123_gentab },
	{ imdi_k124, imdi_k124_gentab },
	{ imdi_k125, imdi_k125_gentab },
	{ imdi_k126, imdi_k126_gentab },
	{ imdi_k127, imdi_k127_gentab },
	{ imdi_k128, imdi_k128_gentab },
	{ imdi_k129, imdi_k129_gentab },
	{ imdi_k130, imdi_k130_gentab },
	{ imdi_k131, imdi_k131_gentab },
	{ imdi_k132, imdi_k132_gentab },
	{ imdi_k133, imdi_k133_gentab },
	{ imdi_k134, imdi_k134_gentab },
	{ imdi_k135, imdi_k135_gentab },
	{ imdi_k136, imdi_k136_gentab },
	{ imdi_k137, imdi_k137_gentab },
	{ imdi_k138, imdi_k138_gentab },
	{ imdi_k139, imdi_k139_gentab },
	{ imdi_k140, imdi_k140_gentab },
	{ imdi_k141, imdi_k141_gentab },
	{ imdi_k142, imdi_k142_gentab },
	{ imdi_k143, imdi_k143_gentab },
	{ imdi_k144, imdi_k144_gentab },
	{ imdi_k145, imdi_k145_gentab },
	{ imdi_k146, imdi_k146_gentab },
	{ imdi_k147, imdi_k147_gentab },
	{ imdi_k148, imdi_k148_gentab },
	{ imdi_k149, imdi_k149_gentab },
	{ imdi_k150, imdi_k150_gentab },
	{ imdi_k151, imdi_k151_gentab },
	{ imdi_k152, imdi_k152_gentab },
	{ imdi_k153, imdi_k153_gentab },
	{ imdi_k154, imdi_k154_gentab },
	{ imdi_k155, imdi_k155_gentab },
	{ imdi_k156, imdi_k156_gentab },
	{ imdi_k157, imdi_k157_gentab },
	{ imdi_k158, imdi_k158_gentab },
	{ imdi_k159, imdi_k159_gentab },
	{ imdi_k160, imdi_k160_gentab },
	{ imdi_k161, imdi_k161_gentab },
	{ imdi_k162, imdi_k162_gentab },
	{ imdi_k163, imdi_k163_gentab },
	{ imdi_k164, imdi_k164_gentab },
	{ imdi_k165, imdi_k165_gentab },
	{ imdi_k166, imdi_k166_gentab },
	{ imdi_k167, imdi_k167_gentab },
	{ imdi_k168, imdi_k168_gentab },
	{ imdi_k169, imdi_k169_gentab },
	{ imdi_k170, imdi_k170_gentab },
	{ imdi_k171, imdi_k171_gentab },
	{ imdi_k172, imdi_k172_gentab },
	{ imdi_k173, imdi_k173_gentab },
	{ imdi_k174, imdi_k174_gentab },
	{ imdi_k175, imdi_k175_gentab },
	{ imdi_k176, imdi_k176_gentab },
	{ imdi_k177, imdi_k177_gentab },
	{ imdi_k178, imdi_k178_gentab },
	{ imdi_k179, imdi_k179_gentab },
	{ imdi_k180, imdi_k180_gentab },
	{ imdi_k181, imdi_k181_gentab },
	{ imdi_k182, imdi_k182_gentab },
	{ imdi_k183, imdi_k183_gentab },
	{ imdi_k184, imdi_k184_gentab },
	{ imdi_k185, imdi_k185_gentab },
	{ imdi_k186, imdi_k186_gentab },
	{ imdi_k187, imdi_k187_gentab },
	{ imdi_k188, imdi_k188_gentab },
	{ imdi_k189, imdi_k189_gentab },
	{ imdi_k190, imdi_k190_gentab },
	{ imdi_k191, imdi_k191_gentab },
	{ imdi_k192, imdi_k192_gentab },
	{ imdi_k193, imdi_k193_gentab },
	{ imdi_k194, imdi_k194_gentab },
	{ imdi_k195, imdi_k195_gentab },
	{ imdi_k196, imdi_k196_gentab },
	{ imdi_k197, imdi_k197_gentab },
	{ imdi_k198, imdi_k198_gentab },
	{ imdi_k199, imdi_k199_gentab },
	{ imdi_k200, imdi_k200_gentab },
	{ imdi_k201, imdi_k201_gentab },
	{ imdi_k202, imdi_k202_gentab },
	{ imdi_k203, imdi_k203_gentab },
	{ imdi_k204, imdi_k204_gentab },
	{ imdi_k205, imdi_k205_gentab },
	{ imdi_k206, imdi_k206_gentab },
	{ imdi_k207, imdi_k207_gentab },
	{ imdi_k208, imdi_k208_gentab },
	{ imdi_k209, imdi_k209_gentab },
	{ imdi_k210, imdi_k210_gentab },
	{ imdi_k211, imdi_k211_gentab },
	{ imdi_k212, imdi_k212_gentab },
	{ imdi_k213, imdi_k213_gentab },
	{ imdi_k214, imdi_k214_gentab },
	{ imdi_k215, imdi_k215_gentab },
	{ imdi_k216, imdi_k216_gentab },
	{ imdi_k217, imdi_k217_gentab },
	{ imdi_k218, imdi_k218_gentab },
	{ imdi_k219, imdi_k219_gentab },
	{ imdi_k220, imdi_k220_gentab },
	{ imdi_k221, imdi_k221_gentab },
	{ imdi_k222, imdi_k222_gentab },
	{ imdi_k223, imdi_k223_gentab },
	{ imdi_k224, imdi_k224_gentab },
	{ imdi_k225, imdi_k225_gentab },
	{ imdi_k226, imdi_k226_gentab },
	{ imdi_k227, imdi_k227_gentab },
	{ imdi_k228, imdi_k228_gentab },
	{ imdi_k229, imdi_k229_gentab },
	{ imdi_k230, imdi_k230_gentab },
	{ imdi_k231, imdi_k231_gentab },
	{ imdi_k232, imdi_k232_gentab },
	{ imdi_k233, imdi_k233_gentab },
	{ imdi_k234, imdi_k234_gentab },
	{ imdi_k235, imdi_k235_gentab },
	{ imdi_k236, imdi_k236_gentab },
	{ imdi_k237, imdi_k237_gentab },
	{ imdi_k238, imdi_k238_gentab },
	{ imdi_k239, imdi_k239_gentab },
	{ imdi_k240, imdi_k240_gentab },
	{ imdi_k241, imdi_k241_gentab },
	{ imdi_k242, imdi_k242_gentab },
	{ imdi_k243, imdi_k243_gentab },
	{ imdi_k244, imdi_k244_gentab },
	{ imdi_k245, imdi_k245_gentab },
	{ imdi_k246, imdi_k246_gentab },
	{ imdi_k247, imdi_k247_gentab },
	{ imdi_k248, imdi_k248_gentab },
	{ imdi_k249, imdi_k249_gentab },
	{ imdi_k250, imdi_k250_gentab },
	{ imdi_k251, imdi_k251_gentab },
	{ imdi_k252, imdi_k252_gentab },
	{ imdi_k253, imdi_k253_gentab },
	{ imdi_k254, imdi_k254_gentab },
	{ imdi_k255, imdi_k255_gentab },
	{ imdi_k256, imdi_k256_gentab },
	{ imdi_k257, imdi_k257_gentab },
	{ imdi_k258, imdi_k258_gentab },
	{ imdi_k259, imdi_k259_gentab },
	{ imdi_k260, imdi_k260_gentab },
	{ imdi_k261, imdi_k261_gentab },
	{ imdi_k262, imdi_k262_gentab },
	{ imdi_k263, imdi_k263_gentab },
	{ imdi_k264, imdi_k264_gentab },
	{ imdi_k265, imdi_k265_gentab },
	{ imdi_k266, imdi_k266_gentab },
	{ imdi_k267, imdi_k267_gentab },
	{ imdi_k268, imdi_k268_gentab },
	{ imdi_k269, imdi_k269_gentab },
	{ imdi_k270, imdi_k270_gentab },
	{ imdi_k271, imdi_k271_gentab },
	{ imdi_k272, imdi_k272_gentab },
	{ imdi_k273, imdi_k273_gentab },
	{ imdi_k274, imdi_k274_gentab },
	{ imdi_k275, imdi_k275_gentab },
	{ imdi_k276, imdi_k276_gentab },
	{ imdi_k277, imdi_k277_gentab },
	{ imdi_k278, imdi_k278_gentab },
	{ imdi_k279, imdi_k279_gentab },
	{ imdi_k280, imdi_k280_gentab },
	{ imdi_k281, imdi_k281_gentab },
	{ imdi_k282, imdi_k282_gentab },
	{ imdi_k283, imdi_k283_gentab },
	{ imdi_k284, imdi_k284_gentab },
	{ imdi_k285, imdi_k285_gentab },
	{ imdi_k286, imdi_k286_gentab },
	{ imdi_k287, imdi_k287_gentab },
	{ imdi_k288, imdi_k288_gentab },
	{ imdi_k289, imdi_k289_gentab },
	{ imdi_k290, imdi_k290_gentab },
	{ imdi_k291, imdi_k291_gentab },
	{ imdi_k292, imdi_k292_gentab },
	{ imdi_k293, imdi_k293_gentab },
	{ imdi_k294, imdi_k294_gentab },
	{ imdi_k295, imdi_k295_gentab },
	{ imdi_k296, imdi_k296_gentab },
	{ imdi_k297, imdi_k297_gentab },
	{ imdi_k298, imdi_k298_gentab },
	{ imdi_k299, imdi_k299_gentab },
	{ imdi_k300, imdi_k300_gentab },
	{ imdi_k301, imdi_k301_gentab },
	{ imdi_k302, imdi_k302_gentab },
	{ imdi_k303, imdi_k303_gentab },
	{ imdi_k304, imdi_k304_gentab },
	{ imdi_k305, imdi_k305_gentab },
	{ imdi_k306, imdi_k306_gentab },
	{ imdi_k307, imdi_k307_gentab },
	{ imdi_k308, imdi_k308_gentab },
	{ imdi_k309, imdi_k309_gentab },
	{ imdi_k310, imdi_k310_gentab },
	{ imdi_k311, imdi_k311_gentab },
	{ imdi_k312, imdi_k312_gentab },
	{ imdi_k313, imdi_k313_gentab },
	{ imdi_k314, imdi_k314_gentab },
	{ imdi_k315, imdi_k315_gentab },
	{ imdi_k316, imdi_k316_gentab },
	{ imdi_k317, imdi_k317_gentab },
	{ imdi_k318, imdi_k318_gentab },
	{ imdi_k319, imdi_k319_gentab },
	{ imdi_k320, imdi_k320_gentab },
	{ imdi_k321, imdi_k321_gentab },
	{ imdi_k322, imdi_k322_gentab },
	{ imdi_k323, imdi_k323_gentab },
	{ imdi_k324, imdi_k324_gentab },
	{ imdi_k325, imdi_k325_gentab },
	{ imdi_k326, imdi_k326_gentab },
	{ imdi_k327, imdi_k327_gentab },
	{ imdi_k328, imdi_k328_gentab },
	{ imdi_k329, imdi_k329_gentab },
	{ imdi_k330, imdi_k330_gentab },
	{ imdi_k331, imdi_k331_gentab },
	{ imdi_k332, imdi_k332_gentab },
	{ imdi_k333, imdi_k333_gentab },
	{ imdi_k334, imdi_k334_gentab },
	{ imdi_k335, imdi_k335_gentab },
	{ imdi_k336, imdi_k336_gentab },
	{ imdi_k337, imdi_k337_gentab },
	{ imdi_k338, imdi_k338_gentab },
	{ imdi_k339, imdi_k339_gentab },
	{ imdi_k340, imdi_k340_gentab },
	{ imdi_k341, imdi_k341_gentab },
	{ imdi_k342, imdi_k342_gentab },
	{ imdi_k343, imdi_k343_gentab },
	{ imdi_k344, imdi_k344_gentab },
	{ imdi_k345, imdi_k345_gentab },
	{ imdi_k346, imdi_k346_gentab },
	{ imdi_k347, imdi_k347_gentab },
	{ imdi_k348, imdi_k348_gentab },
	{ imdi_k349, imdi_k349_gentab },
	{ imdi_k350, imdi_k350_gentab },
	{ imdi_k351, imdi_k351_gentab },
	{ imdi_k352, imdi_k352_gentab },
	{ imdi_k353, imdi_k353_gentab },
	{ imdi_k354, imdi_k354_gentab },
	{ imdi_k355, imdi_k355_gentab },
	{ imdi_k356, imdi_k356_gentab },
	{ imdi_k357, imdi_k357_gentab },
	{ imdi_k358, imdi_k358_gentab },
	{ imdi_k359, imdi_k359_gentab },
	{ imdi_k360, imdi_k360_gentab },
	{ imdi_k361, imdi_k361_gentab },
	{ imdi_k362, imdi_k362_gentab },
	{ imdi_k363, imdi_k363_gentab },
	{ imdi_k364, imdi_k364_gentab },
	{ imdi_k365, imdi_k365_gentab },
	{ imdi_k366, imdi_k366_gentab },
	{ imdi_k367, imdi_k367_gentab },
	{ imdi_k368, imdi_k368_gentab },
	{ imdi_k369, imdi_k369_gentab },
	{ imdi_k370, imdi_k370_gentab },
	{ imdi_k371, imdi_k371_gentab },
	{ imdi_k372, imdi_k372_gentab },
	{ imdi_k373, imdi_k373_gentab },
	{ imdi_k374, imdi_k374_gentab },
	{ imdi_k375, imdi_k375_gentab },
	{ imdi_k376, imdi_k376_gentab },
	{ imdi_k377, imdi_k377_gentab },
	{ imdi_k378, imdi_k378_gentab },
	{ imdi_k379, imdi_k379_gentab },
	{ imdi_k380, imdi_k380_gentab },
	{ imdi_k381, imdi_k381_gentab },
	{ imdi_k382, imdi_k382_gentab },
	{ imdi_k383, imdi_k383_gentab },
	{ imdi_k384, imdi_k384_gentab },
	{ imdi_k385, imdi_k385_gentab },
	{ imdi_k386, imdi_k386_gentab },
	{ imdi_k387, imdi_k387_gentab },
	{ imdi_k388, imdi_k388_gentab },
	{ imdi_k389, imdi_k389_gentab },
	{ imdi_k390, imdi_k390_gentab },
	{ imdi_k391, imdi_k391_gentab },
	{ imdi_k392, imdi_k392_gentab },
	{ imdi_k393, imdi_k393_gentab },
	{ imdi_k394, imdi_k394_gentab },
	{ imdi_k395, imdi_k395_gentab },
	{ imdi_k396, imdi_k396_gentab },
	{ imdi_k397, imdi_k397_gentab },
	{ imdi_k398, imdi_k398_gentab },
	{ imdi_k399, imdi_k399_gentab },
	{ imdi_k400, imdi_k400_gentab },
	{ imdi_k401, imdi_k401_gentab },
	{ imdi_k402, imdi_k402_gentab },
	{ imdi_k403, imdi_k403_gentab },
	{ imdi_k404, imdi_k404_gentab },
	{ imdi_k405, imdi_k405_gentab },
	{ imdi_k406, imdi_k406_gentab },
	{ imdi_k407, imdi_k407_gentab },
	{ imdi_k408, imdi_k408_gentab },
	{ imdi_k409, imdi_k409_gentab },
	{ imdi_k410, imdi_k410_gentab },
	{ imdi_k411, imdi_k411_gentab },
	{ imdi_k412, imdi_k412_gentab },
	{ imdi_k413, imdi_k413_gentab },
	{ imdi_k414, imdi_k414_gentab },
	{ imdi_k415, imdi_k415_gentab },
	{ imdi_k416, imdi_k416_gentab },
	{ imdi_k417, imdi_k417_gentab },
	{ imdi_k418, imdi_k418_gentab },
	{ imdi_k419, imdi_k419_gentab },
	{ imdi_k420, imdi_k420_gentab },
	{ imdi_k421, imdi_k421_gentab },
	{ imdi_k422, imdi_k422_gentab },
	{ imdi_k423, imdi_k423_gentab },
	{ imdi_k424, imdi_k424_gentab },
	{ imdi_k425, imdi_k425_gentab },
	{ imdi_k426, imdi_k426_gentab },
	{ imdi_k427, imdi_k427_gentab },
	{ imdi_k428, imdi_k428_gentab },
	{ imdi_k429, imdi_k429_gentab },
	{ imdi_k430, imdi_k430_gentab },
	{ imdi_k431, imdi_k431_gentab },
	{ imdi_k432, imdi_k432_gentab },
	{ imdi_k433, imdi_k433_gentab },
	{ imdi_k434, imdi_k434_gentab },
	{ imdi_k435, imdi_k435_gentab },
	{ imdi_k436, imdi_k436_gentab },
	{ imdi_k437, imdi_k437_gentab },
	{ imdi_k438, imdi_k438_gentab },
	{ imdi_k439, imdi_k439_gentab },
	{ imdi_k440, imdi_k440_gentab },
	{ imdi_k441, imdi_k441_gentab },
	{ imdi_k442, imdi_k442_gentab },
	{ imdi_k443, imdi_k443_gentab },
	{ imdi_k444, imdi_k444_gentab },
	{ imdi_k445, imdi_k445_gentab },
	{ imdi_k446, imdi_k446_gentab },
	{ imdi_k447, imdi_k447_gentab },
	{ imdi_k448, imdi_k448_gentab },
	{ imdi_k449, imdi_k449_gentab },
	{ imdi_k450, imdi_k450_gentab },
	{ imdi_k451, imdi_k451_gentab },
	{ imdi_k452, imdi_k452_gentab },
	{ imdi_k453, imdi_k453_gentab },
	{ imdi_k454, imdi_k454_gentab },
	{ imdi_k455, imdi_k455_gentab },
	{ imdi_k456, imdi_k456_gentab },
	{ imdi_k457, imdi_k457_gentab },
	{ imdi_k458, imdi_k458_gentab },
	{ imdi_k459, imdi_k459_gentab },
	{ imdi_k460, imdi_k460_gentab },
	{ imdi_k461, imdi_k461_gentab },
	{ imdi_k462, imdi_k462_gentab },
	{ imdi_k463, imdi_k463_gentab },
	{ imdi_k464, imdi_k464_gentab },
	{ imdi_k465, imdi_k465_gentab },
	{ imdi_k466, imdi_k466_gentab },
	{ imdi_k467, imdi_k467_gentab },
	{ imdi_k468, imdi_k468_gentab },
	{ imdi_k469, imdi_k469_gentab },
	{ imdi_k470, imdi_k470_gentab },
	{ imdi_k471, imdi_k471_gentab },
	{ imdi_k472, imdi_k472_gentab },
	{ imdi_k473, imdi_k473_gentab },
	{ imdi_k474, imdi_k474_gentab },
	{ imdi_k475, imdi_k475_gentab },
	{ imdi_k476, imdi_k476_gentab },
	{ imdi_k477, imdi_k477_gentab },
	{ imdi_k478, imdi_k478_gentab },
	{ imdi_k479, imdi_k479_gentab },
	{ imdi_k480, imdi_k480_gentab },
	{ imdi_k481, imdi_k481_gentab },
	{ imdi_k482, imdi_k482_gentab },
	{ imdi_k483, imdi_k483_gentab },
	{ imdi_k484, imdi_k484_gentab },
	{ imdi_k485, imdi_k485_gentab },
	{ imdi_k486, imdi_k486_gentab },
	{ imdi_k487, imdi_k487_gentab },
	{ imdi_k488, imdi_k488_gentab },
	{ imdi_k489, imdi_k489_gentab },
	{ imdi_k490, imdi_k490_gentab },
	{ imdi_k491, imdi_k491_gentab },
	{ imdi_k492, imdi_k492_gentab },
	{ imdi_k493, imdi_k493_gentab },
	{ imdi_k494, imdi_k494_gentab },
	{ imdi_k495, imdi_k495_gentab },
	{ imdi_k496, imdi_k496_gentab },
	{ imdi_k497, imdi_k497_gentab },
	{ imdi_k498, imdi_k498_gentab },
	{ imdi_k499, imdi_k499_gentab },
	{ imdi_k500, imdi_k500_gentab },
	{ imdi_k501, imdi_k501_gentab },
	{ imdi_k502, imdi_k502_gentab },
	{ imdi_k503, imdi_k503_gentab },
	{ imdi_k504, imdi_k504_gentab },
	{ imdi_k505, imdi_k505_gentab },
	{ imdi_k506, imdi_k506_gentab },
	{ imdi_k507, imdi_k507_gentab },
	{ imdi_k508, imdi_k508_gentab },
	{ imdi_k509, imdi_k509_gentab },
	{ imdi_k510, imdi_k510_gentab },
	{ imdi_k511, imdi_k511_gentab },
	{ imdi_k512, imdi_k512_gentab },
	{ imdi_k513, imdi_k513_gentab },
	{ imdi_k514, imdi_k514_gentab },
	{ imdi_k515, imdi_k515_gentab },
	{ imdi_k516, imdi_k516_gentab },
	{ imdi_k517, imdi_k517_gentab },
	{ imdi_k518, imdi_k518_gentab },
	{ imdi_k519, imdi_k519_gentab },
	{ imdi_k520, imdi_k520_gentab },
	{ imdi_k521, imdi_k521_gentab },
	{ imdi_k522, imdi_k522_gentab },
	{ imdi_k523, imdi_k523_gentab },
	{ imdi_k524, imdi_k524_gentab },
	{ imdi_k525, imdi_k525_gentab },
	{ imdi_k526, imdi_k526_gentab },
	{ imdi_k527, imdi_k527_gentab },
	{ imdi_k528, imdi_k528_gentab },
	{ imdi_k529, imdi_k529_gentab },
	{ imdi_k530, imdi_k530_gentab },
	{ imdi_k531, imdi_k531_gentab },
	{ imdi_k532, imdi_k532_gentab },
	{ imdi_k533, imdi_k533_gentab },
	{ imdi_k534, imdi_k534_gentab },
	{ imdi_k535, imdi_k535_gentab },
	{ imdi_k536, imdi_k536_gentab },
	{ imdi_k537, imdi_k537_gentab },
	{ imdi_k538, imdi_k538_gentab },
	{ imdi_k539, imdi_k539_gentab },
	{ imdi_k540, imdi_k540_gentab },
	{ imdi_k541, imdi_k541_gentab },
	{ imdi_k542, imdi_k542_gentab },
	{ imdi_k543, imdi_k543_gentab },
	{ imdi_k544, imdi_k544_gentab },
	{ imdi_k545, imdi_k545_gentab },
	{ imdi_k546, imdi_k546_gentab },
	{ imdi_k547, imdi_k547_gentab },
	{ imdi_k548, imdi_k548_gentab },
	{ imdi_k549, imdi_k549_gentab },
	{ imdi_k550, imdi_k550_gentab },
	{ imdi_k551, imdi_k551_gentab },
	{ imdi_k552, imdi_k552_gentab },
	{ imdi_k553, imdi_k553_gentab },
	{ imdi_k554, imdi_k554_gentab },
	{ imdi_k555, imdi_k555_gentab },
	{ imdi_k556, imdi_k556_gentab },
	{ imdi_k557, imdi_k557_gentab },
	{ imdi_k558, imdi_k558_gentab },
	{ imdi_k559, imdi_k559_gentab },
	{ imdi_k560, imdi_k560_gentab },
	{ imdi_k561, imdi_k561_gentab },
	{ imdi_k562, imdi_k562_gentab },
	{ imdi_k563, imdi_k563_gentab },
	{ imdi_k564, imdi_k564_gentab },
	{ imdi_k565, imdi_k565_gentab },
	{ imdi_k566, imdi_k566_gentab },
	{ imdi_k567, imdi_k567_gentab },
	{ imdi_k568, imdi_k568_gentab },
	{ imdi_k569, imdi_k569_gentab },
	{ imdi_k570, imdi_k570_gentab },
	{ imdi_k571, imdi_k571_gentab },
	{ imdi_k572, imdi_k572_gentab },
	{ imdi_k573, imdi_k573_gentab },
	{ imdi_k574, imdi_k574_gentab },
	{ imdi_k575, imdi_k575_gentab },
	{ imdi_k576, imdi_k576_gentab },
	{ imdi_k577, imdi_k577_gentab },
	{ imdi_k578, imdi_k578_gentab },
	{ imdi_k579, imdi_k579_gentab },
	{ imdi_k580, imdi_k580_gentab },
	{ imdi_k581, imdi_k581_gentab },
	{ imdi_k582, imdi_k582_gentab },
	{ imdi_k583, imdi_k583_gentab },
	{ imdi_k584, imdi_k584_gentab },
	{ imdi_k585, imdi_k585_gentab },
	{ imdi_k586, imdi_k586_gentab },
	{ imdi_k587, imdi_k587_gentab },
	{ imdi_k588, imdi_k588_gentab },
	{ imdi_k589, imdi_k589_gentab },
	{ imdi_k590, imdi_k590_gentab },
	{ imdi_k591, imdi_k591_gentab },
	{ imdi_k592, imdi_k592_gentab },
	{ imdi_k593, imdi_k593_gentab },
	{ imdi_k594, imdi_k594_gentab },
	{ imdi_k595, imdi_k595_gentab },
	{ imdi_k596, imdi_k596_gentab },
	{ imdi_k597, imdi_k597_gentab },
	{ imdi_k598, imdi_k598_gentab },
	{ imdi_k599, imdi_k599_gentab },
	{ imdi_k600, imdi_k600_gentab },
	{ imdi_k601, imdi_k601_gentab },
	{ imdi_k602, imdi_k602_gentab },
	{ imdi_k603, imdi_k603_gentab },
	{ imdi_k604, imdi_k604_gentab },
	{ imdi_k605, imdi_k605_gentab },
	{ imdi_k606, imdi_k606_gentab },
	{ imdi_k607, imdi_k607_gentab },
	{ imdi_k608, imdi_k608_gentab },
	{ imdi_k609, imdi_k609_gentab },
	{ imdi_k610, imdi_k610_gentab },
	{ imdi_k611, imdi_k611_gentab },
	{ imdi_k612, imdi_k612_gentab },
	{ imdi_k613, imdi_k613_gentab },
	{ imdi_k614, imdi_k614_gentab },
	{ imdi_k615, imdi_k615_gentab },
	{ imdi_k616, imdi_k616_gentab },
	{ imdi_k617, imdi_k617_gentab },
	{ imdi_k618, imdi_k618_gentab },
	{ imdi_k619, imdi_k619_gentab },
	{ imdi_k620, imdi_k620_gentab },
	{ imdi_k621, imdi_k621_gentab },
	{ imdi_k622, imdi_k622_gentab },
	{ imdi_k623, imdi_k623_gentab },
	{ imdi_k624, imdi_k624_gentab },
	{ imdi_k625, imdi_k625_gentab },
	{ imdi_k626, imdi_k626_gentab },
	{ imdi_k627, imdi_k627_gentab },
	{ imdi_k628, imdi_k628_gentab },
	{ imdi_k629, imdi_k629_gentab },
	{ imdi_k630, imdi_k630_gentab },
	{ imdi_k631, imdi_k631_gentab },
	{ imdi_k632, imdi_k632_gentab },
	{ imdi_k633, imdi_k633_gentab },
	{ imdi_k634, imdi_k634_gentab },
	{ imdi_k635, imdi_k635_gentab },
	{ imdi_k636, imdi_k636_gentab },
	{ imdi_k637, imdi_k637_gentab },
	{ imdi_k638, imdi_k638_gentab },
	{ imdi_k639, imdi_k639_gentab },
	{ imdi_k640, imdi_k640_gentab },
	{ imdi_k641, imdi_k641_gentab },
	{ imdi_k642, imdi_k642_gentab },
	{ imdi_k643, imdi_k643_gentab },
	{ imdi_k644, imdi_k644_gentab },
	{ imdi_k645, imdi_k645_gentab },
	{ imdi_k646, imdi_k646_gentab },
	{ imdi_k647, imdi_k647_gentab },
	{ imdi_k648, imdi_k648_gentab },
	{ imdi_k649, imdi_k649_gentab },
	{ imdi_k650, imdi_k650_gentab },
	{ imdi_k651, imdi_k651_gentab },
	{ imdi_k652, imdi_k652_gentab },
	{ imdi_k653, imdi_k653_gentab },
	{ imdi_k654, imdi_k654_gentab },
	{ imdi_k655, imdi_k655_gentab },
	{ imdi_k656, imdi_k656_gentab },
	{ imdi_k657, imdi_k657_gentab },
	{ imdi_k658, imdi_k658_gentab },
	{ imdi_k659, imdi_k659_gentab },
	{ imdi_k660, imdi_k660_gentab },
	{ imdi_k661, imdi_k661_gentab },
	{ imdi_k662, imdi_k662_gentab },
	{ imdi_k663, imdi_k663_gentab },
	{ imdi_k664, imdi_k664_gentab },
	{ imdi_k665, imdi_k665_gentab },
	{ imdi_k666, imdi_k666_gentab },
	{ imdi_k667, imdi_k667_gentab },
	{ imdi_k668, imdi_k668_gentab },
	{ imdi_k669, imdi_k669_gentab },
	{ imdi_k670, imdi_k670_gentab },
	{ imdi_k671, imdi_k671_gentab },
	{ imdi_k672, imdi_k672_gentab },
	{ imdi_k673, imdi_k673_gentab },
	{ imdi_k674, imdi_k674_gentab },
	{ imdi_k675, imdi_k675_gentab },
	{ imdi_k676, imdi_k676_gentab },
	{ imdi_k677, imdi_k677_gentab },
	{ imdi_k678, imdi_k678_gentab },
	{ imdi_k679, imdi_k679_gentab },
	{ imdi_k680, imdi_k680_gentab },
	{ imdi_k681, imdi_k681_gentab },
	{ imdi_k682, imdi_k682_gentab },
	{ imdi_k683, imdi_k683_gentab },
	{ imdi_k684, imdi_k684_gentab },
	{ imdi_k685, imdi_k685_gentab },
	{ imdi_k686, imdi_k686_gentab },
	{ imdi_k687, imdi_k687_gentab },
	{ imdi_k688, imdi_k688_gentab },
	{ imdi_k689, imdi_k689_gentab },
	{ imdi_k690, imdi_k690_gentab },
	{ imdi_k691, imdi_k691_gentab },
	{ imdi_k692, imdi_k692_gentab },
	{ imdi_k693, imdi_k693_gentab },
	{ imdi_k694, imdi_k694_gentab },
	{ imdi_k695, imdi_k695_gentab },
	{ imdi_k696, imdi_k696_gentab },
	{ imdi_k697, imdi_k697_gentab },
	{ imdi_k698, imdi_k698_gentab },
	{ imdi_k699, imdi_k699_gentab },
	{ imdi_k700, imdi_k700_gentab },
	{ imdi_k701, imdi_k701_gentab },
	{ imdi_k702, imdi_k702_gentab },
	{ imdi_k703, imdi_k703_gentab },
	{ imdi_k704, imdi_k704_gentab },
	{ imdi_k705, imdi_k705_gentab },
	{ imdi_k706, imdi_k706_gentab },
	{ imdi_k707, imdi_k707_gentab },
	{ imdi_k708, imdi_k708_gentab },
	{ imdi_k709, imdi_k709_gentab },
	{ imdi_k710, imdi_k710_gentab },
	{ imdi_k711, imdi_k711_gentab },
	{ imdi_k712, imdi_k712_gentab },
	{ imdi_k713, imdi_k713_gentab },
	{ imdi_k714, imdi_k714_gentab },
	{ imdi_k715, imdi_k715_gentab },
	{ imdi_k716, imdi_k716_gentab },
	{ imdi_k717, imdi_k717_gentab },
	{ imdi_k718, imdi_k718_gentab },
	{ imdi_k719, imdi_k719_gentab },
	{ imdi_k720, imdi_k720_gentab },
	{ imdi_k721, imdi_k721_gentab },
	{ imdi_k722, imdi_k722_gentab },
	{ imdi_k723, imdi_k723_gentab },
	{ imdi_k724, imdi_k724_gentab },
	{ imdi_k725, imdi_k725_gentab },
	{ imdi_k726, imdi_k726_gentab },
	{ imdi_k727, imdi_k727_gentab },
	{ imdi_k728, imdi_k728_gentab },
	{ imdi_k729, imdi_k729_gentab },
	{ imdi_k730, imdi_k730_gentab },
	{ imdi_k731, imdi_k731_gentab },
	{ imdi_k732, imdi_k732_gentab },
	{ imdi_k733, imdi_k733_gentab },
	{ imdi_k734, imdi_k734_gentab },
	{ imdi_k735, imdi_k735_gentab },
	{ imdi_k736, imdi_k736_gentab },
	{ imdi_k737, imdi_k737_gentab },
	{ imdi_k738, imdi_k738_gentab },
	{ imdi_k739, imdi_k739_gentab },
	{ imdi_k740, imdi_k740_gentab },
	{ imdi_k741, imdi_k741_gentab },
	{ imdi_k742, imdi_k742_gentab },
	{ imdi_k743, imdi_k743_gentab },
	{ imdi_k744, imdi_k744_gentab },
	{ imdi_k745, imdi_k745_gentab },
	{ imdi_k746, imdi_k746_gentab },
	{ imdi_k747, imdi_k747_gentab },
	{ imdi_k748, imdi_k748_gentab },
	{ imdi_k749, imdi_k749_gentab },
	{ imdi_k750, imdi_k750_gentab },
	{ imdi_k751, imdi_k751_gentab },
	{ imdi_k752, imdi_k752_gentab },
	{ imdi_k753, imdi_k753_gentab },
	{ imdi_k754, imdi_k754_gentab },
	{ imdi_k755, imdi_k755_gentab },
	{ imdi_k756, imdi_k756_gentab },
	{ imdi_k757, imdi_k757_gentab },
	{ imdi_k758, imdi_k758_gentab },
	{ imdi_k759, imdi_k759_gentab },
	{ imdi_k760, imdi_k760_gentab },
	{ imdi_k761, imdi_k761_gentab },
	{ imdi_k762, imdi_k762_gentab },
	{ imdi_k763, imdi_k763_gentab },
	{ imdi_k764, imdi_k764_gentab },
	{ imdi_k765, imdi_k765_gentab },
	{ imdi_k766, imdi_k766_gentab },
	{ imdi_k767, imdi_k767_gentab },
	{ imdi_k768, imdi_k768_gentab },
	{ imdi_k769, imdi_k769_gentab },
	{ imdi_k770, imdi_k770_gentab },
	{ imdi_k771, imdi_k771_gentab },
	{ imdi_k772, imdi_k772_gentab },
	{ imdi_k773, imdi_k773_gentab },
	{ imdi_k774, imdi_k774_gentab },
	{ imdi_k775, imdi_k775_gentab },
	{ imdi_k776, imdi_k776_gentab },
	{ imdi_k777, imdi_k777_gentab },
	{ imdi_k778, imdi_k778_gentab },
	{ imdi_k779, imdi_k779_gentab },
	{ imdi_k780, imdi_k780_gentab },
	{ imdi_k781, imdi_k781_gentab },
	{ imdi_k782, imdi_k782_gentab },
	{ imdi_k783, imdi_k783_gentab },
	{ imdi_k784, imdi_k784_gentab },
	{ imdi_k785, imdi_k785_gentab },
	{ imdi_k786, imdi_k786_gentab },
	{ imdi_k787, imdi_k787_gentab },
	{ imdi_k788, imdi_k788_gentab },
	{ imdi_k789, imdi_k789_gentab },
	{ imdi_k790, imdi_k790_gentab },
	{ imdi_k791, imdi_k791_gentab },
	{ imdi_k792, imdi_k792_gentab },
	{ imdi_k793, imdi_k793_gentab },
	{ imdi_k794, imdi_k794_gentab },
	{ imdi_k795, imdi_k795_gentab },
	{ imdi_k796, imdi_k796_gentab },
	{ imdi_k797, imdi_k797_gentab },
	{ imdi_k798, imdi_k798_gentab },
	{ imdi_k799, imdi_k799_gentab },
	{ imdi_k800, imdi_k800_gentab },
	{ imdi_k801, imdi_k801_gentab },
	{ imdi_k802, imdi_k802_gentab },
	{ imdi_k803, imdi_k803_gentab },
	{ imdi_k804, imdi_k804_gentab },
	{ imdi_k805, imdi_k805_gentab },
	{ imdi_k806, imdi_k806_gentab },
	{ imdi_k807, imdi_k807_gentab },
	{ imdi_k808, imdi_k808_gentab },
	{ imdi_k809, imdi_k809_gentab },
	{ imdi_k810, imdi_k810_gentab },
	{ imdi_k811, imdi_k811_gentab },
	{ imdi_k812, imdi_k812_gentab },
	{ imdi_k813, imdi_k813_gentab },
	{ imdi_k814, imdi_k814_gentab },
	{ imdi_k815, imdi_k815_gentab },
	{ imdi_k816, imdi_k816_gentab },
	{ imdi_k817, imdi_k817_gentab },
	{ imdi_k818, imdi_k818_gentab },
	{ imdi_k819, imdi_k819_gentab },
	{ imdi_k820, imdi_k820_gentab },
	{ imdi_k821, imdi_k821_gentab },
	{ imdi_k822, imdi_k822_gentab },
	{ imdi_k823, imdi_k823_gentab },
	{ imdi_k824, imdi_k824_gentab },
	{ imdi_k825, imdi_k825_gentab },
	{ imdi_k826, imdi_k826_gentab },
	{ imdi_k827, imdi_k827_gentab },
	{ imdi_k828, imdi_k828_gentab },
	{ imdi_k829, imdi_k829_gentab },
	{ imdi_k830, imdi_k830_gentab },
	{ imdi_k831, imdi_k831_gentab },
	{ imdi_k832, imdi_k832_gentab },
	{ imdi_k833, imdi_k833_gentab },
	{ imdi_k834, imdi_k834_gentab },
	{ imdi_k835, imdi_k835_gentab },
	{ imdi_k836, imdi_k836_gentab },
	{ imdi_k837, imdi_k837_gentab },
	{ imdi_k838, imdi_k838_gentab },
	{ imdi_k839, imdi_k839_gentab },
	{ imdi_k840, imdi_k840_gentab },
	{ imdi_k841, imdi_k841_gentab },
	{ imdi_k842, imdi_k842_gentab },
	{ imdi_k843, imdi_k843_gentab },
	{ imdi_k844, imdi_k844_gentab },
	{ imdi_k845, imdi_k845_gentab },
	{ imdi_k846, imdi_k846_gentab },
	{ imdi_k847, imdi_k847_gentab },
	{ imdi_k848, imdi_k848_gentab },
	{ imdi_k849, imdi_k849_gentab },
	{ imdi_k850, imdi_k850_gentab },
	{ imdi_k851, imdi_k851_gentab },
	{ imdi_k852, imdi_k852_gentab },
	{ imdi_k853, imdi_k853_gentab },
	{ imdi_k854, imdi_k854_gentab },
	{ imdi_k855, imdi_k855_gentab },
	{ imdi_k856, imdi_k856_gentab },
	{ imdi_k857, imdi_k857_gentab },
	{ imdi_k858, imdi_k858_gentab },
	{ imdi_k859, imdi_k859_gentab },
	{ imdi_k860, imdi_k860_gentab },
	{ imdi_k861, imdi_k861_gentab },
	{ imdi_k862, imdi_k862_gentab },
	{ imdi_k863, imdi_k863_gentab },
	{ imdi_k864, imdi_k864_gentab },
	{ imdi_k865, imdi_k865_gentab },
	{ imdi_k866, imdi_k866_gentab },
	{ imdi_k867, imdi_k867_gentab },
	{ imdi_k868, imdi_k868_gentab },
	{ imdi_k869, imdi_k869_gentab },
	{ imdi_k870, imdi_k870_gentab },
	{ imdi_k871, imdi_k871_gentab },
	{ imdi_k872, imdi_k872_gentab },
	{ imdi_k873, imdi_k873_gentab },
	{ imdi_k874, imdi_k874_gentab },
	{ imdi_k875, imdi_k875_gentab },
	{ imdi_k876, imdi_k876_gentab },
	{ imdi_k877, imdi_k877_gentab },
	{ imdi_k878, imdi_k878_gentab },
	{ imdi_k879, imdi_k879_gentab },
	{ imdi_k880, imdi_k880_gentab },
	{ imdi_k881, imdi_k881_gentab },
	{ imdi_k882, imdi_k882_gentab },
	{ imdi_k883, imdi_k883_gentab },
	{ imdi_k884, imdi_k884_gentab },
	{ imdi_k885, imdi_k885_gentab },
	{ imdi_k886, imdi_k886_gentab },
	{ imdi_k887, imdi_k887_gentab },
	{ imdi_k888, imdi_k888_gentab },
	{ imdi_k889, imdi_k889_gentab },
	{ imdi_k890, imdi_k890_gentab },
	{ imdi_k891, imdi_k891_gentab }
};

int no_kfuncs = 891;
