//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef ISOTAB
#define ISOTAB


#include <QString>


struct IsoTab
{
    IsoTab(int NumIso);
    ~IsoTab();
    QString getIsoName(int N);

    inline int getIsoIndex(int mI1, int mI2)
    {
        for (int n=0; n < numIso; ++n) if ((mNumIso1[n] == mI1 && mNumIso2[n] == mI2) || (mNumIso1[n] == mI2 && mNumIso2[n] == mI1)) return n;
        return -1;
    }
	
    QString *chSymb1, *chSymb2, *texName;
    int numIso, *mNumIso1, *mNumIso2, refIso, *JStep;
    double *relNA, *redMass, *relRedMass, *rootRRM, *mIso1, *mIso2;
};

#endif
