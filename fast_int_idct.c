
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "jpeg.h"
#define Y(i,j)                Y[8*i+j]
#define X(i,j)                (output->linear[i*8 + j])

/* This version is IEEE compliant using 16-bit arithmetic. */
/* The number of bits coefficients are scaled up before 2-D IDCT: */
#define S_BITS                 3
/* The number of bits in the fractional part of a fixed point constant: */
#define C_BITS                14
#define SCALE(x,n)        ((x) << (n))



/* This version is vital in passing overall mean error test. */
//#define DESCALE(x, n)        (((x) + (1 << ((n)-1)) - ((x) < 0)) >> (n))
#define DESCALE(x, n)        (((x) + (1 << (S_BITS + 2)) - ((x) < 0)) >> (S_BITS + 3))

#define ADD(x, y)        (x + y)
#define SUB(x, y)        (x - y)
#define CMUL(C, x)        (((C) * (x) + (0x2000)) >> C_BITS)

#define C_BITS2		16
#define CMUL2(C, x)	(((C) * (x) + (1 << (C_BITS2-1))) >> C_BITS2)


/* Butterfly: but(a,b,x,y) = rot(sqrt(2),4,a,b,x,y) */
#define but(a,b,x,y)        { x = SUB(a,b); y = ADD(a,b); }

#define CK1 35468	
#define SK1 85627	
#define ADDK1 121095	
#define SUBK1 50159

#define CK2 54491	
#define SK2 36410	
#define ADDK2 90901	
#define SUBK2 -18081

#define CK3 64277	
#define SK3 12785	
#define ADDK3 77062	
#define SUBK3 -51491

#define SQRT2 92682



/* Inverse 1-D Discrete Cosine Transform.
   Result Y is scaled up by factor sqrt(8).
   Original Loeffler algorithm.
*/

static inline void idct_1d(int *Y)
{
	int z1[8], z2[8], z3[8];
	int t;

	/* Stage 1: */
	but(Y[0], Y[4], z1[1], z1[0]);

	/* rot(sqrt(2), 6, Y[2], Y[6], &z1[2], &z1[3]); */
	t=CMUL2(CK1,ADD(Y[2],Y[6]));
	z1[2]=SUB(t,CMUL2(ADDK1,Y[6]));
	z1[3]=ADD(t,CMUL2(SUBK1,Y[2]));




	but(Y[1], Y[7], z1[4], z1[7]);
	/* z1[5] = CMUL(sqrt(2), Y[3]);
	   z1[6] = CMUL(sqrt(2), Y[5]);
	 */
	z1[5] = CMUL2(SQRT2, Y[3]);
	z1[6] = CMUL2(SQRT2, Y[5]);

	/* Stage 2: */
	but(z1[0], z1[3], z2[3], z2[0]);
	but(z1[1], z1[2], z2[2], z2[1]);
	but(z1[4], z1[6], z2[6], z2[4]);
	but(z1[7], z1[5], z2[5], z2[7]);

	/* Stage 3: */
	z3[0] = z2[0];
	z3[1] = z2[1];
	z3[2] = z2[2];
	z3[3] = z2[3];

	/* rot(1, 3, z2[4], z2[7], &z3[4], &z3[7]); */
	t=CMUL2(CK2,ADD(z2[4],z2[7]));
	z3[4]=SUB(t,CMUL2(ADDK2,z2[7]));	
	z3[7]=ADD(t,CMUL2(SUBK2,z2[4]));




	/* rot(1, 1, z2[5], z2[6], &z3[5], &z3[6]); */
	t=CMUL2(CK3,ADD(z2[5],z2[6]));
	z3[5]=SUB(t,CMUL2(ADDK3,z2[6]));	
	z3[6]=ADD(t,CMUL2(SUBK3,z2[5]));

	/* Final stage 4: */
	but(z3[0], z3[7], Y[7], Y[0]);
	but(z3[1], z3[6], Y[6], Y[1]);
	but(z3[2], z3[5], Y[5], Y[2]);
	but(z3[3], z3[4], Y[4], Y[3]);

}

void IDCT(const FBlock * input, PBlock * output)
{
        int Y[64];

        int k, l;

        /* Pass 1: process rows. */

        for (l = 0; l < 8; l++)
        {
                /* Prescale k-th row: */
                int Yc[8];
                for (k = 0; k < 8; k++) 
                        Yc[k] = SCALE(input->linear[l*8 + k], S_BITS);
                /* 1-D IDCT on k-th row: */
                idct_1d(Yc);
                /* Result Y is scaled up by factor sqrt(8)*2^S_BITS. */
                for (k = 0; k < 8; k++) 
                        Y[8* l + k] = Yc[k];
        }

        /* Pass 2: process columns. */
        for (l = 0; l < 8; l++) 
        {
                int Yc[8];
                for (k = 0; k < 8; k++)
                        Yc[k] = Y[8 * k + l];
                /* 1-D IDCT on l-th column: */
                idct_1d(Yc);
                /* Result is once more scaled up by a factor sqrt(8). */
                for (k = 0; k < 8; k++) 
                {
                        int r = 128 + DESCALE(Yc[k], S_BITS + 3);        /* includes level shift */

                        /* Clip to 8 bits unsigned: */
                        r = r > 0 ? (r < 255 ? r : 255) : 0;
                        X(k, l) = r;
                }
        }
}
