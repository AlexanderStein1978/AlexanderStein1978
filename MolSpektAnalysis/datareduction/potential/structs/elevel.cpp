//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "elevel.h"


bool ELevel::isGT(ELevel L)
{    
    if (Iso > L.Iso) return true;
    if (Iso < L.Iso) return false;
    if (J > L.J) return true;
    if (J < L.J) return false;
    if (v > L.v) return true;
    return false;
}

bool ELevel::isGT(int I, int cJ, int cv)
{    
    if (Iso > I) return true;
    if (Iso < I) return false;
    if (J > cJ) return true;
    if (J < cJ) return false;
    if (v > cv) return true;
    return false;
}
