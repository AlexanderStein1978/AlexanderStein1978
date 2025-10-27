//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef ISOTAB
#define ISOTAB


#include <QString>

class QPixmap;


struct IsoTab
{
    IsoTab();
    ~IsoTab();
    int getIsoIndex(int mI1, int mI2);
    void createArrays();
	
    QString *chSymb1 = nullptr, *chSymb2 = nullptr, *texName = nullptr, *Isoname = nullptr;
    int numIso = 0, *mNumIso1 = nullptr, *mNumIso2 = nullptr, refIso = 0, *JStep = nullptr;
    double *relNA = nullptr, *redMass = nullptr, *relRedMass = nullptr, *rootRRM = nullptr, *mIso1 = nullptr, *mIso2 = nullptr;
    QPixmap *IsoImage = nullptr;
};

#endif
