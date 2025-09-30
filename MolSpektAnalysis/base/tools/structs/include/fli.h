//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FLI_H
#define FLI_H


struct Marker;

class ElState;


struct FLI
{
    FLI(int N, int G);
	
	Marker **marker;
    FLI *Next, *First;
    int *vss, *Jss, *vs, *Js, *Iso, N, G;
    double *DD, *SNR, *oc;
	bool *satellite, *overlap;
	ElState **lState, **uState;
};

#endif
