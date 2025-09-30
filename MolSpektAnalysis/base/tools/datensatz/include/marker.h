//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MARKER_H
#define MARKER_H


#include <QString>


class ElState;
class Molecule;


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
