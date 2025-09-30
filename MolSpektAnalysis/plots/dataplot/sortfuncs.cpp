//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "sortfuncs.h"


int sortUS(int *L1, int *L2)
{
	if (L1 == 0) return -1;
	if (L2 == 0) return 1;
	if (L1[0] > L2[0]) return 1;
	if (L1[0] < L2[0]) return -1;
	if (L1[1] > L2[1]) return 1;
	if (L1[1] < L2[1]) return -1;
	if (L1[2] > L2[2]) return 1;
	if (L1[2] < L2[2]) return -1;
	return 0;
}

int sortLS(int *L1, int *L2)
{
	if (L1 == 0) return -1;
	if (L2 == 0) return 1;
	if (L1[0] > L2[0]) return 1;
	if (L1[0] < L2[0]) return -1;
	if (L1[3] > L2[3]) return 1;
	if (L1[3] < L2[3]) return -1;
	if (L1[4] > L2[4]) return 1;
	if (L1[4] < L2[4]) return -1;
	return 0;
}
