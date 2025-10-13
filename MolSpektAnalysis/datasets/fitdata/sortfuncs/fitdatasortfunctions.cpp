//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "fitdatasortfunctions.h"
#include "fitdata.h"
#include "tableline.h"

#include <cmath>


bool sortIvJFreqF(const FitDataCore * const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    int in = Tab->getIso(n), Jsn = Tab->getJ(n), Jsm = Tab->getJ(m);
    int im = Tab->getIso(m), Jssn = Tab->getJs(n), Jssm = Tab->getJs(m);
    if (in < im) return true;
    if (in > im) return false;
    in = Tab->get_v(n), im = Tab->get_v(m);
    if (in < im) return true;
    if (in > im) return false;
    if ((in = abs(Jssn - Jsn)) < (im = abs(Jsm - Jssm))) return true;
    if (in > im) return false;
    if (Jsn < Jsm) return true;
    if (Jsn > Jsm) return false;
    if ((in = stoi(Tab->get_vs(n))) < (im = stoi(Tab->get_vs(m))))
        return true;
    if (in > im) return false;
    if (Jssn < Jssm) return true;
    if (Jssn > Jssm) return false;
    double Intn = Tab->getEnergy(n), Intm = Tab->getEnergy(m);
    if (Intn > Intm) return true;
    return false;
}

bool SortvsIvJ(const FitDataCore *const Tab, const int n, const int m)
{
    int in, im;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = stoi(Tab->get_vs(n))) > (im = stoi(Tab->get_vs(m))))
        return true;
    if (in < im) return false;
    if ((in = Tab->getIso(n)) < (im = Tab->getIso(m)))
        return true;
    if (in > im) return false;
    if ((in = Tab->get_v(n)) < (im = Tab->get_v(m)))
        return true;
    if (in > im) return false;
    if ((in = Tab->getJ(n)) < (im = Tab->getJ(m)))
        return true;
    if (in > im) return false;
    if ((in = Tab->getJs(n)) < (im = Tab->getJs(m)))
        return true;
    return false;
}

bool sortIefJFreq(const FitDataCore *const Tab, const int n, const int m)
{
    int in, im, Jn, Jm;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->getIso(n)) < (im = Tab->getIso(m))) return true;
    if (in > im) return false;
    if ((in = abs((Jn = Tab->getJ(n)) - Tab->getJs(n)))
         < (im = abs((Jm = Tab->getJ(m)) - Tab->getJs(m)))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if (Tab->getEnergy(n) < Tab->getEnergy(m)) return true;
    return false;
}

bool sortIefJvFreq(const FitDataCore *const Tab, const int n, const int m)
{
    int in, im, Jn, Jm;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->getIso(n)) < (im = Tab->getIso(m))) return true;
    if (in > im) return false;
    if ((in = abs((Jn = Tab->getJ(n)) - Tab->getJs(n)))
          < (im = abs((Jm = Tab->getJ(m)) - Tab->getJs(m)))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if ((in = Tab->get_v(n)) < (im = Tab->get_v(m))) return true;
    if (in > im) return false;
    if (Tab->getEnergy(n) < Tab->getEnergy(m)) return true;
    return false;
}

bool sortIefJFreqv(const FitDataCore *const Tab, const int n, const int m)
{
    int in, im, Jn, Jm;
    double Im, In;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = Tab->getIso(n)) < (im = Tab->getIso(m))) return true;
    if (in > im) return false;
    if ((in = abs((Jn = Tab->getJ(n)) - Tab->getJs(n)))
         < (im = abs((Jm = Tab->getJ(m)) - Tab->getJs(m)))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if ((In = Tab->getEnergy(n)) < (Im = Tab->getEnergy(m))) return true;
    if (In > Im) return false;
    if ((in = Tab->get_v(n)) < (im = Tab->get_v(m))) return true;
    if (in > im) return false;
    if ((in = stoi(Tab->get_vs(n))) < (im = stoi(Tab->get_vs(m)))) return true;
    return false;
}

bool sortByProg(const FitDataCore *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    std::string LTn = Tab->getSource(n), LTm = Tab->getSource(m);
    if (LTn < LTm) return true;
    if (LTn > LTm) return false;
    if (Tab->getProgression(n) < Tab->getProgression(m)) return true;
    return false;
}

bool sortbyDeviation(const FitDataCore *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (fabs(Tab->getObsCalc(n)) > fabs(Tab->getObsCalc(m))) return true;
    return false;
}

bool sortbyDevR(const FitDataCore *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (fabs(Tab->getDevRatio(n)) > fabs(Tab->getDevRatio(m))) return true;
    return false;
}

bool sortforTFGS(const FitDataCore *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    std::string Bn = Tab->get_vs(n), Bm = Tab->get_vs(m);
    if (Bn == "TE" && Bm != "TE") return true;
    if (Bn != "TE" && Bm == "TE") return false;
    std::string Sn = Tab->getSource(n), Sm = Tab->getSource(m);
    if (Sn < Sm) return true;
    if (Sn > Sm) return false;
    if (Bn == "nA" && Bm != "nA") return true;
    if (Bn != "nA" && Bm == "nA") return false;
    int vn = stoi(Bn), vm = stoi(Bm);
    if (vn < vm) return true;
    if (vn > vm) return false;
    int Jn = Tab->getJs(n), Jm = Tab->getJs(m);
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if (Tab->getProgression(n) < Tab->getProgression(m)) return true;
    return false;
}

bool sortByElState(const FitDataCore *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (Tab->getOtherState(n) < Tab->getOtherState(m)) return true;
    return false;
}

bool sortForExtractNewOrChanged(const FitDataCore * const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const std::string S1(Tab->getOtherState(n)), S2(Tab->getOtherState(m));
    if (S1 < S2) return true;
    if (S1 > S2) return false;
    const int Iso1 = Tab->getIso(n), Iso2 = Tab->getIso(m);
    if (Iso1 < Iso2) return true;
    if (Iso1 > Iso2) return false;
    const int v1 = Tab->get_v(n), v2 = Tab->get_v(m);
    if (v1 < v2) return true;
    if (v1 > v2) return false;
    const int J1 = Tab->getJ(n), J2 = Tab->getJ(m);
    if (J1 < J2) return true;
    if (J1 > J2) return false;
    const bool ef1 = (Tab->getJs(n) == J1), ef2 = (Tab->getJs(m) == J2);
    if (ef1 && !ef2) return true;
    if (!ef1 && ef2) return false;
    const std::string Source1(Tab->getSource(n)), Source2(Tab->getSource(m));
    if (Source1 < Source2) return true;
    if (Source1 > Source2) return false;
    const double E1 = Tab->getEnergy(n), E2 = Tab->getEnergy(m);
    if (E1 < E2) return true;
    return false;
}
