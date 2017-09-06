#include "floatx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*--------------------------------------------------------------------------------
	Return floatx representation (as defined by *def) which
	best approximates value
-------------------------------------------------------------------------------- */
floatx doubleToFloatx(const floatxDef *def, double value) {
	union {
		floatx fx;
		double d;
		//anonymous union adapted from xmp_float that keeps track of value
	} convert;
	long check;					//eventually will be used to keep track of the sign
	convert.d = value;

	check = (convert.fx & 0x8000000000000000)>>63;  	//grabs first bit: sign
	long exp = (convert.fx & 0x7ff0000000000000)>>52;   //grabs next 11 exp bits
	floatx frac = (convert.fx & 0x000fffffffffffff);    //grabs remaining frac bits
	exp -= 1023; 
	exp += (1<<(def->expBits-1))-1; 	//unbias the exponent value & rebias it

	//check if the value is special, i.e. should be interpreted as INFINITY
	if (exp >= ((1<<(def->expBits))-1)) { 
		int i;
		//INF has s=0, exp=1111 1111 1..., and frac=0000 0000 0...
		for(i=0;i<def->expBits;i++) {
			if (i==0) { exp = 1<<((def->expBits-i)-1); continue; }
			exp |= (1<<((def->expBits-i)-1));
		}
		frac = 0;
	} else {
		if (exp<=0) { 						//subnormal case
			int comp = -(def->totBits-def->expBits-1);
			if (exp <= comp) {
				return 0; 					//if the number is too small, return zero
			} else {
				int trunc = 0;
				frac = frac << 12;
				long j;
				for (j=0;j>=exp;j--) {
					frac = frac >> 1;
					if (j==0) {
						frac |= 0x8000000000000000;
					}
				} 				//ensures proper fraction represenation since exp is 0
				frac = frac >> (64 - (def->totBits) + (def->expBits));
				if (frac & 1) { trunc = 1; }
				frac = (frac >> 1) + trunc;
				exp = 0;
			}
		} else {
			int truncate = 0; 				//rounding purposes
			//intentionally shift one less than we need to see if we have to round
			int eval = (52 - (def->totBits) + (def->expBits));
			if (eval>=0) { 				//if eval is negative (unlikely), shift left;
				frac = frac >> eval;	//otherwise shift right
			} else {
				frac = frac << abs(eval);
			}
			if (frac & 1) { truncate = 1; } //is that last bit 1? if so, account for it
			frac = (frac >> 1) + truncate;
		}
	}

	//assign the bits and return
	convert.fx = 0;
	convert.fx |= (check << (def->totBits-1));
	convert.fx |= (exp << (def->totBits-def->expBits-1));
	convert.fx |= frac;

	return convert.fx;
}

/** Return C double with value which best approximates that of floatx fx
 *  (as defined by *def).
 */
double floatxToDouble(const floatxDef *def, floatx fx) {
	union {
		floatx value;
		double d;
		//another union because doubles are difficult to deal with bitwise
	} convert;
	if (fx==0) { return 0.0; }
	convert.value = 0;
	//SIGN BIT
	convert.value |= (fx >> (def->totBits-1)) << 63;
	
	floatx exp = 0;
	if (!convert.value) {
		exp = fx >> (def->totBits-def->expBits-1); //handle negative case. if negative,
	} else { 					   				   //the sign bit will cause the exp
		exp = fx >> (def->totBits-def->expBits-1); //to seem larger than it is, and this
		int a;					   				   //decision structure fixes it 
		floatx compare = 0;
		for (a=0;a<def->expBits;a++) {
			compare |= 1 << a;
		}
		exp = exp & compare;
	}
	floatx bias = (1 << (def->expBits-1)) - 1;
	long UNBIASED = (exp-(bias));			//bias stuff, this is important for the
	floatx newBias = UNBIASED + 1023; 		//newBias variable
	floatx mask_off = 0;
	int k;
	for(k=0;k<(def->totBits-def->expBits-1);k++) {
		mask_off |= 1 << k;					//gets FRAC with a mask (molly percocet)
	}
	floatx FRAC = fx & mask_off;
	if (exp >= ((1 << (def->expBits)) - 1)) {
		if (convert.value) {
			return -INFINITY;               //SPECIAL CASE if exp is too big
		} return INFINITY;
	}
	if (exp == 0) {
		int counter = 0;
		if (FRAC) {
			while(!(FRAC >> (counter+(def->totBits)-(def->expBits)-2))) {
				FRAC = FRAC << 1; 			//DENORMAL CASE if exp is zero but
				newBias -= 1;				//frac>0
				counter += 1;
			}
			FRAC = FRAC << 1;
		}
	} 

	//EXP BITS
	convert.value |= (newBias << 52);

	//FRAC BITS
	FRAC = FRAC << (64 - (def->totBits) + (def->expBits) + 1);
	FRAC = FRAC >> 12;
	convert.value |= FRAC;
	return convert.d;
}
