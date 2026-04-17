//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "linetablesortfunctions.h"
#include "linetable.h"
#include "termenergy.h"

#include <QTableWidget>

#include <cmath>


bool isnSPG(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (core->getIso(n) - 1) / 10;
	int im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m))) return true;
	if (in > im) return false;
	int Jsn = ltc->getJs(n), Jsm = ltc->getJs(m);
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	int Jssn = ltc->getJss(n), Jssm = ltc->getJss(m), dn, dm;
	if ((dn = abs(Jssn - Jsn)) < (dm = abs(Jssm - Jsm))) return true;
	if (dn > dm) return false;
	QString Fn = ltc->getSourceFile(n), Fm = ltc->getSourceFile(m);
<<<<<<< HEAD
	in = Fn.lastIndexOf(QRegularExpression("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegularExpression("[\\/]")) + 1;
=======
	in = Fn.lastIndexOf(QRegExp("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegExp("[\\/]")) + 1;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
	int jn = Fn.indexOf('.', in), jm = Fm.indexOf('.', im);
	Fn = Fn.mid(in, (jn >= 0 ? jn : Fn.length()) - in);
	Fm = Fm.mid(im, (jm >= 0 ? jm : Fm.length()) - im);
	//printf("F1: %s, F2: %s\n", Fn.toAscii().data(), Fm.toAscii().data());
	if ((in = QString::compare(Fn, Fm, Qt::CaseInsensitive)) < 0) return true;
	if (in > 0) return false;
	if ((in = ltc->getProgression(n)) < (im = ltc->getProgression(m))) return true;
	if (in > im) return false;
	if ((in = ltc->get_vss(n)) < (im = ltc->get_vss(m))) return true;
	if (in > im) return false;
	if (Jssn < Jssm) return true;
	return false;
}

bool isnSPG(TermEnergy &T1, TermEnergy &T2)
{
	if (T1.Iso == -1) return false;
	if (T2.Iso == -1) return true;
	if (T1.Iso < T2.Iso) return true;
	if (T1.Iso > T2.Iso) return false;
	if (T1.v < T2.v) return true;
	if (T1.v > T2.v) return false;
	if (T1.J < T2.J) return true;
	if (T1.J > T2.J) return false;
	if (!T1.ef && T2.ef) return true;
	return false;
}

bool sortUtIvJ(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	double En = ltc->getUpperEnergy(n), Em = ltc->getUpperEnergy(m);
	if (En < Em) return true;
	if (En > Em) return false;
	int in = (ltc->getIso(n) - 1) / 10;
	int im = (ltc->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m)))	return true;
	if (in > im) return false;
	int Jsn = ltc->getJs(n), Jsm = ltc->getJs(m);
	if (Jsn < Jsm) return true;
	return false;		
}

bool sortIJvP(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (core->getIso(n) - 1) / 10;
	int im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	int Jsn = core->getJs(n), Jsm = core->getJs(m);
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m))) return true;
	if (in > im) return false;
	QString Fn = ltc->getSourceFile(n), Fm = ltc->getSourceFile(m);
	if (Fn < Fm) return true;
	if (Fn > Fm) return false;
	int Jssn = ltc->getJss(n), Jssm = ltc->getJss(m);
	if ((in = abs(Jsn - Jssn)) < (im = abs(Jsm - Jssm))) return true;
	if (in > im) return false;
	if ((in = ltc->getFineStructureQN(n)) < (im = ltc->getFineStructureQN(m))) return true;
	if (in > im) return false;
	if ((in = ltc->get_vss(n)) < (im = ltc->get_vss(m))) return true;
	if (in > im) return false;
	if (Jssn < Jssm) return true;
	return false;
}

bool sortForSPN(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int jn, in = (core->getIso(n) - 1) / 10;
	int jm, im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m)))
		return true;
	if (in > im) return false;
	int Jsn = core->getJs(n), Jsm = core->getJs(m);
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	int Jssn = ltc->getJss(n), Jssm = ltc->getJss(m), dn, dm;
	if ((dn = abs(Jssn - Jsn)) < (dm = abs(Jssm - Jsm))) return true;
	if (dn > dm) return false;
	QString Fn = core->getSourceFile(n), Fm = core->getSourceFile(m);
<<<<<<< HEAD
	in = Fn.lastIndexOf(QRegularExpression("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegularExpression("[\\/]")) + 1;
=======
	in = Fn.lastIndexOf(QRegExp("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegExp("[\\/]")) + 1;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
	jn = Fn.indexOf('.', in);
	jm = Fm.indexOf('.', im);
	Fn = Fn.mid(in, (jn >= 0 ? jn : Fn.length()) - in);
	Fm = Fm.mid(im, (jm >= 0 ? jm : Fm.length()) - im);
	//printf("F1: %s, F2: %s\n", Fn.toAscii().data(), Fm.toAscii().data());
	if ((in = QString::compare(Fn, Fm, Qt::CaseInsensitive)) < 0) return true;
	if (in > 0) return false;
	if ((in = ltc->getFineStructureQN(n)) < (im = ltc->getFineStructureQN(m))) return true;
	if (in > im || core->columnCount() <= LineTableCore::CEUp) return false;
	if (ltc->getUpperEnergy(n) < ltc->getUpperEnergy(m)) return true;
	return false;
}

bool sortIvPJ(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (core->getIso(n) - 1) / 10;
	int im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m))) return true;
	if (in > im) return false;
	int Jsn = ltc->getJs(n), Jsm = ltc->getJs(m);
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	if ((in = ltc->getFineStructureQN(n)) < (im = ltc->getFineStructureQN(m))) return true;
	if (in > im) return false;
	QString Fn = ltc->getSourceFile(n), Fm = ltc->getSourceFile(m);
	if (Fn < Fm) return true;
	if (Fn > Fm) return false;
	if ((in = ltc->get_vss(n)) < (im = ltc->get_vss(m))) return true;
	if (in > im) return false;
	int Jssn = ltc->getJss(n), Jssm = ltc->getJss(m);
	if (Jssn < Jssm) return true;
	return false;
}

bool sortFPInt(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	QString Fn = core->getSourceFile(n), Fm = core->getSourceFile(m);
	if (Fn < Fm) return true;
	if (Fn > Fm) return false;
	int in = (core->getIso(n) - 1) / 10;
	int im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m))) return true;
	if (in > im) return false;
	int Jsn = ltc->getJs(n), Jsm = ltc->getJs(m);
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	if ((in = ltc->getFineStructureQN(n)) < (im = ltc->getFineStructureQN(m))) return true;
	if (in > im) return false;
	double Intn = ltc->getSNR(n), Intm = ltc->getSNR(m);
	if (Intn > Intm) return true;
	return false;
}

bool sortIJvFreq(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (core->getIso(n) - 1) / 10;
	int im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	in = core->getJs(n), im = core->getJs(m);
	if (in < im) return true;
	if (in > im) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m)))
		return true;
	if (in > im) return false;
	if ((in = ltc->getJss(n)) < (im = ltc->getJss(m)))
		return true;
	if (in > im) return false;
	if ((in = ltc->get_vss(n)) < (im = ltc->get_vss(m)))
		return true;
	if (in > im) return false;
	double Intn = ltc->getWaveNumber(n), Intm = ltc->getWaveNumber(m);
	if (Intn > Intm) return true;
	return false;
}

bool sortfRemDoubl(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int jn, in = (core->getIso(n) - 1) / 10;
	int jm, im = (core->getIso(m) - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	in = core->getJs(n), im = core->getJs(m);
	if (in < im) return true;
	if (in > im) return false;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if ((in = ltc->getJss(n)) < (im = ltc->getJss(m))) return true;
	if (in > im) return false;
	if ((in = ltc->get_vss(n)) < (im = ltc->get_vss(m))) return true;
	if (in > im) return false;
	QString Fn = ltc->getSourceFile(n), Fm = ltc->getSourceFile(m);
<<<<<<< HEAD
	in = Fn.lastIndexOf(QRegularExpression("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegularExpression("[\\/]")) + 1;
=======
	in = Fn.lastIndexOf(QRegExp("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegExp("[\\/]")) + 1;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
	jn = Fn.indexOf('.', in);
	jm = Fm.indexOf('.', im);
	Fn = Fn.mid(in, (jn >= 0 ? jn : Fn.length()) - in);
	Fm = Fm.mid(im, (jm >= 0 ? jm : Fm.length()) - im);
	if ((in = QString::compare(Fn, Fm, Qt::CaseInsensitive)) < 0) return true;
	if (in > 0) return false;
	double Intn = ltc->getWaveNumber(n), Intm = ltc->getWaveNumber(m);
	if (Intn > Intm) return true;
	if (Intn < Intm) return false;
	if ((in = ltc->get_vs(n)) < (im = ltc->get_vs(m))) return true;
	return false;
}

bool sortBySpectrum(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	QString S1 = core->getSourceFile(n), S2 = core->getSourceFile(m);
<<<<<<< HEAD
    const int i = S1.lastIndexOf(QRegularExpression("[\\/]"));
=======
    const int i = S1.lastIndexOf(QRegExp("[\\/]"));
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    S1 = S1.right(S1.length() - i - 1);
	S2 = S2.right(S2.length() - S2.lastIndexOf(QRegularExpression("[\\/]")) - 1);
	return (QString::compare(S1, S2, Qt::CaseInsensitive) < 0 ? true : false);
}

bool sortByvs(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if (ltc->get_vs(n) < ltc->get_vs(m)) return true;
	return false;
}

bool sortByFrequency(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if (ltc->getWaveNumber(n) < ltc->getWaveNumber(m)) return true;
	return false;
}

bool sortByProgression(const TableViewWindowCore *const core, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (core->getProgression(n) < core->getProgression(m)) return true;
    return false;
}

bool sortByJs(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->get_vs(n) < ltc->get_vs(m)) return true;
    return false;
}

bool sortBy_vss(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->get_vss(n) < ltc->get_vss(m)) return true;
    return false;
}

bool sortByJss(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getJss(n) < ltc->getJss(m)) return true;
    return false;
}

bool sortByF(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getFineStructureQN(n) < ltc->getFineStructureQN(m)) return true;
    return false;
}

bool sortBy_err(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
    if (core->getUncertainty(n) < core->getUncertainty(m)) return true;
    return false;
}

bool sortByIso(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
    if (core->getIso(n) < core->getIso(m)) return true;
    return false;
}

bool sortByFile(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
    if (core->getSourceFile(n) < core->getSourceFile(m)) return true;
    return false;
}

bool sortBySNR(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getSNR(n) < ltc->getSNR(m)) return true;
    return false;
}

bool sortByDev(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
    if (core->getObsCalc(n) < core->getObsCalc(m)) return true;
    return false;
}

bool sortByComment(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getComment(n) < ltc->getComment(m)) return true;
    return false;
}

bool sortByFCF(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getFCF(n) < ltc->getFCF(m)) return true;
    return false;
}

bool sortByEUp(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getUpperEnergy(n) < ltc->getUpperEnergy(m)) return true;
    return false;
}

bool sortByEav(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getAverageUpperEnergy(n) < ltc->getAverageUpperEnergy(m)) return true;
    return false;
}

bool sortByEUma(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getDiffToAverageUpperEnergy(n) < ltc->getDiffToAverageUpperEnergy(m)) return true;
    return false;
}

bool sortByEdJ(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
	if (ltc->getDiffToNextJ(n) < ltc->getDiffToNextJ(m)) return true;
    return false;
}

bool sortByCalc(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getCalculatedUpperEnergy(n) < ltc->getCalculatedUpperEnergy(m)) return true;
    return false;
}

bool sortByOmC(const TableViewWindowCore *const core, const int n, const int m)
{
	if (n==-1) return false;
    if (m==-1) return true;
	const LineTableCore* ltc = reinterpret_cast<const LineTableCore*>(core);
    if (ltc->getDeviationToCalculatedUpperEnergy(n) < ltc->getDeviationToCalculatedUpperEnergy(m)) return true;
    return false;
}
