//
// C++ Implementation: TableLineSortFunctor
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "tablelinesortfunctor.h"


bool TableLineSortFunctor::operator()(int n, int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    int sn, sm, ln, lm, in, im;
    for (sn = 0, ln = n; ln >= NTLS[sn]; ln -= NTLS[sn++]) ;
    for (sm = 0, lm = m; lm >= NTLS[sm]; lm -= NTLS[sm++]) ;
    if (TLToSort[sn][ln].Iso < TLToSort[sm][lm].Iso) return true;
    if (TLToSort[sn][ln].Iso > TLToSort[sm][lm].Iso) return false;
    if ((in = abs(TLToSort[sn][ln].Jss - TLToSort[sn][ln].Js)) < (im = abs(TLToSort[sm][lm].Jss - TLToSort[sm][lm].Js))) return true;
    if (in > im) return false;
    if (TLToSort[sn][ln].Jss < TLToSort[sm][lm].Jss) return true;
    if (TLToSort[sn][ln].Jss > TLToSort[sm][lm].Jss) return false;
    if (TLToSort[sn][ln].vss < TLToSort[sm][lm].vss) return true;
    if (TLToSort[sn][ln].vss > TLToSort[sm][lm].vss) return false;
    if (TLToSort[sn][ln].WN < TLToSort[sm][lm].WN) return true;
    return false;
}
