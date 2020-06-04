//
// C++ Interface: Marker
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MARKER_H
#define MARKER_H


struct Marker
{
    const double *Line, *LMin, *RMin;
    bool Marked, Visible, DisplayData, satellite;
    int x1, y1, x2, y2, vss, Jss, Iso, vs, Js, FC;
    double DD, SNR, oc, HFLM, uncertainty;
	bool overlap, sOverlap, b, added;
	QString IsoName, lState, uState;
	ElState *LState, *UState;
	Molecule *Mol;
};

#endif
