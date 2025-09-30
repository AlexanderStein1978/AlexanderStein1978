//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "fitdatasortfunctions.h"
#include "fitdata.h"
#include "tableline.h"

#include <QTableWidget>

#include <cmath>


bool sortIvJFreqF(const QTableWidget * const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    int in = Tab->item(n, 0)->text().toInt(), Jsn = Tab->item(n, 2)->text().toInt(), Jsm = Tab->item(m, 2)->text().toInt();
    int im = Tab->item(m, 0)->text().toInt(), Jssn = Tab->item(n, 4)->text().toInt(), Jssm = Tab->item(m, 4)->text().toInt();
    if (in < im) return true;
    if (in > im) return false;
    in = Tab->item(n, 1)->text().toInt(), im = Tab->item(m, 1)->text().toInt();
    if (in < im) return true;
    if (in > im) return false;
    if ((in = abs(Jssn - Jsn)) < (im = abs(Jsm - Jssm))) return true;
    if (in > im) return false;
    if (Jsn < Jsm) return true;
    if (Jsn > Jsm) return false;
    if ((in = Tab->item(n, 3)->text().toInt()) < (im = Tab->item(m, 3)->text().toInt())) 
        return true;
    if (in > im) return false;
    if (Jssn < Jssm) return true;
    if (Jssn > Jssm) return false;
    double Intn = Tab->item(n, 8)->text().toDouble(), Intm = Tab->item(m, 8)->text().toDouble();
    if (Intn > Intm) return true;
    return false;
}

bool SortvsIvJ(const QTableWidget *const Tab, const int n, const int m)
{
    int in, im;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->item(n, 3)->text().toInt()) > (im = Tab->item(m, 3)->text().toInt()))
        return true;
    if (in < im) return false;
    if ((in = Tab->item(n, 0)->text().toInt()) < (im = Tab->item(m, 0)->text().toInt()))
        return true;
    if (in > im) return false;
    if ((in = Tab->item(n, 1)->text().toInt()) < (im = Tab->item(m, 1)->text().toInt()))
        return true;
    if (in > im) return false;
    if ((in = Tab->item(n, 2)->text().toInt()) < (im = Tab->item(m, 2)->text().toInt()))
        return true;
    if (in > im) return false;
    if ((in = Tab->item(n, 4)->text().toInt()) < (im = Tab->item(m, 4)->text().toInt()))
        return true;
    return false;
}

bool sortIefJFreq(const QTableWidget *const Tab, const int n, const int m)
{
    int in, im, Jn, Jm;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->item(n, 0)->text().toInt()) < (im = Tab->item(m, 0)->text().toInt())) return true;
    if (in > im) return false;
    if ((in = fabs((Jn = Tab->item(n, 2)->text().toInt()) - Tab->item(n, 4)->text().toInt()))  
         < (im = fabs((Jm = Tab->item(m, 2)->text().toInt()) - Tab->item(m, 4)->text().toInt()))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if (Tab->item(n, 8)->text().toDouble() < Tab->item(m, 8)->text().toDouble()) return true;
    return false;
}

bool sortIefJvFreq(const QTableWidget *const Tab, const int n, const int m)
{
    int in, im, Jn, Jm;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->item(n, 0)->text().toInt()) < (im = Tab->item(m, 0)->text().toInt())) return true;
    if (in > im) return false;
    if ((in = fabs((Jn = Tab->item(n, 2)->text().toInt()) - Tab->item(n, 4)->text().toInt()))
          < (im = fabs((Jm = Tab->item(m, 2)->text().toInt()) - Tab->item(m, 4)->text().toInt()))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if ((in = Tab->item(n, fdcv)->text().toInt()) < (im = Tab->item(m, fdcv)->text().toInt())) return true;
    if (in > im) return false;
    if (Tab->item(n, 8)->text().toDouble() < Tab->item(m, 8)->text().toDouble()) return true;
    return false;
}

bool sortIefJFreqv(const QTableWidget *const Tab, const int n, const int m)
{
    int in, im, Jn, Jm;
    double Im, In;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->item(n, 0)->text().toInt()) < (im = Tab->item(m, 0)->text().toInt())) return true;
    if (in > im) return false;
    if ((in = fabs((Jn = Tab->item(n, 2)->text().toInt()) - Tab->item(n, 4)->text().toInt()))  
         < (im = fabs((Jm = Tab->item(m, 2)->text().toInt()) - Tab->item(m, 4)->text().toInt()))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if ((In = Tab->item(n, 8)->text().toDouble()) < (Im = Tab->item(m, 8)->text().toDouble())) return true;
    if (In > Im) return false;
    if ((in = Tab->item(n, 1)->text().toInt()) < (im = Tab->item(m, 1)->text().toInt())) return true;
    if (in > im) return false;
    if ((in = Tab->item(n, 3)->text().toInt()) < (im = Tab->item(m, 3)->text().toInt())) return true;
    return false;
}

bool sortByProg(const QTableWidget *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    QString LTn = Tab->item(n, 5)->text(), LTm = Tab->item(m, 5)->text();
    if (LTn < LTm) return true;
    if (LTn > LTm) return false;
    if (Tab->item(n, 6)->text().toInt() < Tab->item(m, 6)->text().toInt()) return true;
    return false;
}

bool sortbyDeviation(const QTableWidget *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (fabs(Tab->item(n, 10)->text().toDouble()) > fabs(Tab->item(m, 10)->text().toDouble())) return true;
    return false;
}

bool sortbyDevR(const QTableWidget *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (fabs(Tab->item(n, 11)->text().toDouble()) > fabs(Tab->item(m, 11)->text().toDouble())) return true;
    return false;
}

bool sortforTFGS(const QTableWidget *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    QString Bn = Tab->item(n, 3)->text(), Bm = Tab->item(m, 3)->text();
    if (Bn == "TE" && Bm != "TE") return true;
    if (Bn != "TE" && Bm == "TE") return false;
    QString Sn = Tab->item(n, 5)->text(), Sm = Tab->item(m, 5)->text();
    if (Sn < Sm) return true;
    if (Sn > Sm) return false;
    if (Bn == "nA" && Bm != "nA") return true;
    if (Bn != "nA" && Bm == "nA") return false;
    int vn = Bn.toInt(), vm = Bm.toInt();
    if (vn < vm) return true;
    if (vn > vm) return false;
    int Jn = Tab->item(n, 4)->text().toInt(), Jm = Tab->item(m, 4)->text().toInt();
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if (Tab->item(n, 6)->text().toInt() < Tab->item(m, 6)->text().toInt()) return true;
    return false;
}

bool sortByElState(const QTableWidget *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (Tab->item(n, fdcLineElState)->text() < Tab->item(m, fdcLineElState)->text()) return true;
    return false;
}

bool sortForExtractNewOrChanged(const QTableWidget * const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const QString S1(Tab->item(n, fdcLineElState)->text()), S2(Tab->item(m, fdcLineElState)->text());
    if (S1 < S2) return true;
    if (S1 > S2) return false;
    const int Iso1 = Tab->item(n, fdcIso)->text().toInt(), Iso2 = Tab->item(m, fdcIso)->text().toInt();
    if (Iso1 < Iso2) return true;
    if (Iso1 > Iso2) return false;
    const int v1 = Tab->item(n, fdcv)->text().toInt(), v2 = Tab->item(m, fdcv)->text().toInt();
    if (v1 < v2) return true;
    if (v1 > v2) return false;
    const int J1 = Tab->item(n, fdcJ)->text().toInt(), J2 = Tab->item(m, fdcJ)->text().toInt();
    if (J1 < J2) return true;
    if (J1 > J2) return false;
    const bool ef1 = (Tab->item(n, fdcJs)->text().toInt() == J1), ef2 = (Tab->item(m, fdcJs)->text().toInt() == J2);
    if (ef1 && !ef2) return true;
    if (!ef1 && ef2) return false;
    const QString Source1(Tab->item(n, fdcSource)->text()), Source2(Tab->item(m, fdcSource)->text());
    if (Source1 < Source2) return true;
    if (Source1 > Source2) return false;
    const double E1 = Tab->item(n, fdcEnergy)->text().toDouble(), E2 = Tab->item(m, fdcEnergy)->text().toDouble();
    if (E1 < E2) return true;
    return false;
}
