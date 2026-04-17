//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "fitdatasortfunctions.h"
#include "fitdata.h"
#include "tableline.h"
#include "utils.h"


bool sortIvJFreqF(const TableViewWindowCore * const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    int in = core->getIso(n), Jsn = core->getJ(n), Jsm = core->getJ(m);
    int im = core->getIso(m), Jssn = core->getJs(n), Jssm = core->getJs(m);
    if (in < im) return true;
    if (in > im) return false;
    in = core->get_v(n), im = core->get_v(m);
    if (in < im) return true;
    if (in > im) return false;
    if ((in = abs(Jssn - Jsn)) < (im = abs(Jsm - Jssm))) return true;
    if (in > im) return false;
    if (Jsn < Jsm) return true;
    if (Jsn > Jsm) return false;
<<<<<<< HEAD
    if ((in = reinterpret_cast<const FitDataCore*>(core)->get_vs(n).toInt()) < (im = reinterpret_cast<const FitDataCore*>(core)->get_vs(m).toInt())) return true;
=======
<<<<<<< HEAD
    if ((in = stdStringToInt(Tab->get_vs(n))) < (im = stdStringToInt(Tab->get_vs(m))))
        return true;
=======
    if ((in = reinterpret_cast<const FitDataCore*>(core)->get_vs(n).toInt()) < (im = reinterpret_cast<const FitDataCore*>(core)->get_vs(m).toInt())) return true;
>>>>>>> linetable
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    if (in > im) return false;
    if (Jssn < Jssm) return true;
    if (Jssn > Jssm) return false;
    double Intn = core->getEnergy(n), Intm = core->getEnergy(m);
    if (Intn > Intm) return true;
    return false;
}

bool SortvsIvJ(const TableViewWindowCore *const core, const int n, const int m)
{
    int in, im;
    if (n==-1) return false;
    if (m==-1) return true;
<<<<<<< HEAD
    if ((in = reinterpret_cast<const FitDataCore*>(core)->get_vs(n).toInt()) < (im = reinterpret_cast<const FitDataCore*>(core)->get_vs(m).toInt())) return true;
=======
<<<<<<< HEAD
    if ((in = stdStringToInt(Tab->get_vs(n))) > (im = stdStringToInt(Tab->get_vs(m))))
        return true;
=======
    if ((in = reinterpret_cast<const FitDataCore*>(core)->get_vs(n).toInt()) < (im = reinterpret_cast<const FitDataCore*>(core)->get_vs(m).toInt())) return true;
>>>>>>> linetable
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    if (in < im) return false;
    if ((in = core->getIso(n)) < (im = core->getIso(m)))
        return true;
    if (in > im) return false;
    if ((in = core->get_v(n)) < (im = core->get_v(m)))
        return true;
    if (in > im) return false;
    if ((in = core->getJ(n)) < (im = core->getJ(m)))
        return true;
    if (in > im) return false;
    if ((in = core->getJs(n)) < (im = core->getJs(m)))
        return true;
    return false;
}

bool sortIefJFreq(const TableViewWindowCore *const core, const int n, const int m)
{
    int in, im, Jn, Jm;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = core->getIso(n)) < (im = core->getIso(m))) return true;
    if (in > im) return false;
    if ((in = abs((Jn = core->getJ(n)) - core->getJs(n)))
         < (im = abs((Jm = core->getJ(m)) - core->getJs(m)))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if (core->getEnergy(n) < core->getEnergy(m)) return true;
    return false;
}

bool sortIefJvFreq(const TableViewWindowCore *const core, const int n, const int m)
{
    int in, im, Jn, Jm;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = core->getIso(n)) < (im = core->getIso(m))) return true;
    if (in > im) return false;
    if ((in = abs((Jn = core->getJ(n)) - core->getJs(n)))
          < (im = abs((Jm = core->getJ(m)) - core->getJs(m)))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if ((in = core->get_v(n)) < (im = core->get_v(m))) return true;
    if (in > im) return false;
    if (core->getEnergy(n) < core->getEnergy(m)) return true;
    return false;
}

bool sortIefJFreqv(const TableViewWindowCore *const core, const int n, const int m)
{
    int in, im, Jn, Jm;
    double Im, In;
    if (n==-1) return false;
    if (m==-1) return true;
    if ((in = core->getIso(n)) < (im = core->getIso(m))) return true;
    if (in > im) return false;
    if ((in = abs((Jn = core->getJ(n)) - core->getJs(n)))
         < (im = abs((Jm = core->getJ(m)) - core->getJs(m)))) return true;
    if (in > im) return false;
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if ((In = core->getEnergy(n)) < (Im = core->getEnergy(m))) return true;
    if (In > Im) return false;
    if ((in = core->get_v(n)) < (im = core->get_v(m))) return true;
    if (in > im) return false;
<<<<<<< HEAD
    if ((in = reinterpret_cast<const FitDataCore*>(core)->get_vs(n).toInt()) < (im = reinterpret_cast<const FitDataCore*>(core)->get_vs(m).toInt())) return true;
=======
<<<<<<< HEAD
    if ((in = stdStringToInt(Tab->get_vs(n))) < (im = stdStringToInt(Tab->get_vs(m)))) return true;
=======
    if ((in = reinterpret_cast<const FitDataCore*>(core)->get_vs(n).toInt()) < (im = reinterpret_cast<const FitDataCore*>(core)->get_vs(m).toInt())) return true;
>>>>>>> linetable
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    return false;
}

bool sortByProg(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    QString LTn = core->getSourceFile(n), LTm = core->getSourceFile(m);
    if (LTn < LTm) return true;
    if (LTn > LTm) return false;
    if (core->getProgression(n) < core->getProgression(m)) return true;
    return false;
}

bool sortbyDeviation(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (fabs(core->getObsCalc(n)) > fabs(core->getObsCalc(m))) return true;
    return false;
}

bool sortbyDevR(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const FitDataCore* ftc = reinterpret_cast<const FitDataCore*>(core);
    if (fabs(ftc->getDevRatio(n)) > fabs(ftc->getDevRatio(m))) return true;
    return false;
}

bool sortforTFGS(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const FitDataCore* ftc = reinterpret_cast<const FitDataCore*>(core);
    QString Bn = ftc->get_vs(n), Bm = ftc->get_vs(m);
    if (Bn == "TE" && Bm != "TE") return true;
    if (Bn != "TE" && Bm == "TE") return false;
    QString Sn = ftc->getSource(n), Sm = ftc->getSource(m);
    if (Sn < Sm) return true;
    if (Sn > Sm) return false;
    if (Bn == "nA" && Bm != "nA") return true;
    if (Bn != "nA" && Bm == "nA") return false;
<<<<<<< HEAD
    int vn = stdStringToInt(Bn), vm = stdStringToInt(Bm);
=======
    int vn = Bn.toInt(), vm = Bm.toInt();
>>>>>>> linetable
    if (vn < vm) return true;
    if (vn > vm) return false;
    int Jn = ftc->getJs(n), Jm = ftc->getJs(m);
    if (Jn < Jm) return true;
    if (Jn > Jm) return false;
    if (ftc->getProgression(n) < ftc->getProgression(m)) return true;
    return false;
}

bool sortByElState(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const FitDataCore* ftc = reinterpret_cast<const FitDataCore*>(core);
    if (ftc->getOtherState(n) < ftc->getOtherState(m)) return true;
    return false;
}

bool sortForExtractNewOrChanged(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const FitDataCore* ftc = reinterpret_cast<const FitDataCore*>(core);
    const QString S1(ftc->getOtherState(n)), S2(ftc->getOtherState(m));
    if (S1 < S2) return true;
    if (S1 > S2) return false;
    const int Iso1 = ftc->getIso(n), Iso2 = ftc->getIso(m);
    if (Iso1 < Iso2) return true;
    if (Iso1 > Iso2) return false;
    const int v1 = ftc->get_v(n), v2 = ftc->get_v(m);
    if (v1 < v2) return true;
    if (v1 > v2) return false;
    const int J1 = ftc->getJ(n), J2 = ftc->getJ(m);
    if (J1 < J2) return true;
    if (J1 > J2) return false;
    const bool ef1 = (ftc->getJs(n) == J1), ef2 = (ftc->getJs(m) == J2);
    if (ef1 && !ef2) return true;
    if (!ef1 && ef2) return false;
    const QString Source1(ftc->getSource(n)), Source2(ftc->getSource(m));
    if (Source1 < Source2) return true;
    if (Source1 > Source2) return false;
    const double E1 = ftc->getEnergy(n), E2 = ftc->getEnergy(m);
    if (E1 < E2) return true;
    return false;
}

bool sortByIsoColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getIso(n) < core->getIso(m)) return true;
    return false;
}

bool sortBy_vColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->get_v(n) < core->get_v(m)) return true;
    return false;
}

bool sortByJColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getJ(n) < core->getJ(m)) return true;
    return false;
}

bool sortBy_vsColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const FitDataCore* ftc = reinterpret_cast<const FitDataCore*>(core);
    if (ftc->get_vs(n) < ftc->get_vs(m)) return true;
    return false;
}

bool sortByJsColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getJs(n) < core->getJs(m)) return true;
    return false;
}

bool sortBySourceColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    const FitDataCore* ftc = reinterpret_cast<const FitDataCore*>(core);
    if (ftc->getSource(n) < ftc->getSource(m)) return true;
    return false;
}

bool sortByProgressionColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getProgression(n) < core->getProgression(m)) return true;
    return false;
}

bool sortByFileColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getSourceFile(n) < core->getSourceFile(m)) return true;
    return false;
}

bool sortByEnergyColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getEnergy(n) < core->getEnergy(m)) return true;
    return false;
}

bool sortByUncertaintyColumn(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getUncertainty(n) < core->getUncertainty(m)) return true;
    return false;
}
