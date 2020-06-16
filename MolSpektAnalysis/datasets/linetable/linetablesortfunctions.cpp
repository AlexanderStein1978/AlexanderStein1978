//
// C++ Implementation: sortfunctions
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2016
//
// Copyright: See README file that comes with this source code
//
//


#include "linetablesortfunctions.h"
#include "linetable.h"
#include "termenergy.h"

#include <QTableWidget>

#include <cmath>


bool isnSPG(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	int Jsn = Tab->item(n, CJs)->text().toInt(), Jsm = Tab->item(m, CJs)->text().toInt();
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	int Jssn = Tab->item(n, CJss)->text().toInt(), Jssm = Tab->item(m, CJss)->text().toInt(), dn, dm;
	if ((dn = fabs(Jssn - Jsn)) < (dm = fabs(Jssm - Jsm))) return true;
	if (dn > dm) return false;
	QString Fn = Tab->item(n, CFile)->text(), Fm = Tab->item(m, CFile)->text();
	in = Fn.lastIndexOf(QRegExp("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegExp("[\\/]")) + 1;
	int jn = Fn.indexOf('.', in), jm = Fm.indexOf('.', im);
	Fn = Fn.mid(in, (jn >= 0 ? jn : Fn.length()) - in);
	Fm = Fm.mid(im, (jm >= 0 ? jm : Fm.length()) - im);
	//printf("F1: %s, F2: %s\n", Fn.toAscii().data(), Fm.toAscii().data());
	if ((in = QString::compare(Fn, Fm, Qt::CaseInsensitive)) < 0) return true;
	if (in > 0) return false;
	if ((in = Tab->item(n, CPN)->text().toInt()) < (im = Tab->item(m, CPN)->text().toInt())) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvss)->text().toInt()) < (im = Tab->item(m, Cvss)->text().toInt())) 
		return true;
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

bool sortUtIvJ(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	double En = Tab->item(n, CEUp)->text().toDouble(), Em = Tab->item(m, CEUp)->text().toDouble();
	if (En < Em) return true;
	if (En > Em) return false;
	int in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	int Jsn = Tab->item(n, CJs)->text().toInt(), Jsm = Tab->item(m, CJs)->text().toInt();
	if (Jsn < Jsm) return true;
	return false;		
}

bool sortIJvP(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	int Jsn = Tab->item(n, CJs)->text().toInt(), Jsm = Tab->item(m, CJs)->text().toInt();
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	QString Fn = Tab->item(n, CFile)->text(), Fm = Tab->item(m, CFile)->text();
	if (Fn < Fm) return true;
	if (Fn > Fm) return false;
	int Jssn = Tab->item(n, CJss)->text().toInt(), Jssm = Tab->item(m, CJss)->text().toInt();
	if ((in = fabs(Jsn - Jssn)) < (im = fabs(Jsm - Jssm))) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, CF)->text().toInt()) < (im = Tab->item(m, CF)->text().toDouble())) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvss)->text().toInt()) < (im = Tab->item(m, Cvss)->text().toInt())) 
		return true;
	if (in > im) return false;
	if (Jssn < Jssm) return true;
	return false;
}

bool sortForSPN(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int jn, in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int jm, im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	int Jsn = Tab->item(n, CJs)->text().toInt(), Jsm = Tab->item(m, CJs)->text().toInt();
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	int Jssn = Tab->item(n, CJss)->text().toInt(), Jssm = Tab->item(m, CJss)->text().toInt(), dn, dm;
	if ((dn = abs(Jssn - Jsn)) < (dm = abs(Jssm - Jsm))) return true;
	if (dn > dm) return false;
	QString Fn = Tab->item(n, CFile)->text(), Fm = Tab->item(m, CFile)->text();
	in = Fn.lastIndexOf(QRegExp("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegExp("[\\/]")) + 1;
	jn = Fn.indexOf('.', in);
	jm = Fm.indexOf('.', im);
	Fn = Fn.mid(in, (jn >= 0 ? jn : Fn.length()) - in);
	Fm = Fm.mid(im, (jm >= 0 ? jm : Fm.length()) - im);
	//printf("F1: %s, F2: %s\n", Fn.toAscii().data(), Fm.toAscii().data());
	if ((in = QString::compare(Fn, Fm, Qt::CaseInsensitive)) < 0) return true;
	if (in > 0) return false;
	if ((in = Tab->item(n, CF)->text().toInt()) < (im = Tab->item(m, CF)->text().toInt())) return true;
	if (in > im || Tab->columnCount() <= CEUp) return false;
	if (Tab->item(n, CEUp)->text().toDouble() < Tab->item(m, CEUp)->text().toDouble()) return true;
	return false;
}

bool sortIvPJ(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	int Jsn = Tab->item(n, CJs)->text().toInt(), Jsm = Tab->item(m, CJs)->text().toInt();
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	if ((in = Tab->item(n, CF)->text().toInt()) < (im = Tab->item(m, CF)->text().toDouble())) return true;
	if (in > im) return false;
	QString Fn = Tab->item(n, CFile)->text(), Fm = Tab->item(m, CFile)->text();
	if (Fn < Fm) return true;
	if (Fn > Fm) return false;
	if ((in = Tab->item(n, Cvss)->text().toInt()) < (im = Tab->item(m, Cvss)->text().toInt())) 
		return true;
	if (in > im) return false;
	int Jssn = Tab->item(n, CJss)->text().toInt(), Jssm = Tab->item(m, CJss)->text().toInt();
	if (Jssn < Jssm) return true;
	return false;
}

bool sortFPInt(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	QString Fn = Tab->item(n, CFile)->text(), Fm = Tab->item(m, CFile)->text();
	if (Fn < Fm) return true;
	if (Fn > Fm) return false;
	int in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	int Jsn = Tab->item(n, CJs)->text().toInt(), Jsm = Tab->item(m, CJs)->text().toInt();
	if (Jsn < Jsm) return true;
	if (Jsn > Jsm) return false;
	if ((in = Tab->item(n, CF)->text().toInt()) < (im = Tab->item(m, CF)->text().toDouble())) return true;
	if (in > im) return false;
	double Intn = Tab->item(n, CSNR)->text().toDouble(), Intm = Tab->item(m, CSNR)->text().toDouble();
	if (Intn > Intm) return true;
	return false;
}

bool sortIJvFreq(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	in = Tab->item(n, CJs)->text().toInt(), im = Tab->item(m, CJs)->text().toInt();
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) 
		return true;
	if (in > im) return false;
	if ((in = Tab->item(n, CJss)->text().toInt()) < (im = Tab->item(m, CJss)->text().toInt())) 
		return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvss)->text().toInt()) < (im = Tab->item(m, Cvss)->text().toInt())) 
		return true;
	if (in > im) return false;
	double Intn = Tab->item(n, CWN)->text().toDouble(), Intm = Tab->item(m, CWN)->text().toDouble();
	if (Intn > Intm) return true;
	return false;
}

bool sortfRemDoubl(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	int jn, in = (Tab->item(n, CIso)->text().toInt() - 1) / 10;
	int jm, im = (Tab->item(m, CIso)->text().toInt() - 1) / 10;
	if (in < im) return true;
	if (in > im) return false;
	in = Tab->item(n, CJs)->text().toInt(), im = Tab->item(m, CJs)->text().toInt();
	if (in < im) return true;
	if (in > im) return false;
	if ((in = Tab->item(n, CJss)->text().toInt()) < (im = Tab->item(m, CJss)->text().toInt())) 
		return true;
	if (in > im) return false;
	if ((in = Tab->item(n, Cvss)->text().toInt()) < (im = Tab->item(m, Cvss)->text().toInt())) 
		return true;
	if (in > im) return false;
	QString Fn = Tab->item(n, CFile)->text(), Fm = Tab->item(m, CFile)->text();
	in = Fn.lastIndexOf(QRegExp("[\\/]")) + 1;
	im = Fm.lastIndexOf(QRegExp("[\\/]")) + 1;
	jn = Fn.indexOf('.', in);
	jm = Fm.indexOf('.', im);
	Fn = Fn.mid(in, (jn >= 0 ? jn : Fn.length()) - in);
	Fm = Fm.mid(im, (jm >= 0 ? jm : Fm.length()) - im);
	if ((in = QString::compare(Fn, Fm, Qt::CaseInsensitive)) < 0) return true;
	if (in > 0) return false;
	double Intn = Tab->item(n, CWN)->text().toDouble(), Intm = Tab->item(m, CWN)->text().toDouble();
	if (Intn > Intm) return true;
	if (Intn < Intm) return false;
	if ((in = Tab->item(n, Cvs)->text().toInt()) < (im = Tab->item(m, Cvs)->text().toInt())) return true;
	return false;
}

bool sortBySpectrum(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	QString S1 = Tab->item(n, CFile)->text(), S2 = Tab->item(m, CFile)->text();
    const int i = S1.lastIndexOf(QRegExp("[\\/]"));
    S1 = S1.right(S1.length() - i - 1);
	S2 = S2.right(S2.length() - S2.lastIndexOf(QRegExp("[\\/]")) - 1);
	return (QString::compare(S1, S2, Qt::CaseInsensitive) < 0 ? true : false);
}

bool sortByvs(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	if (Tab->item(n, Cvs)->text().toInt() < Tab->item(m, Cvs)->text().toInt()) return true;
	return false;
}

bool sortByFrequency(const QTableWidget *const Tab, const int n, const int m)
{
	if (n==-1) return false;
	if (m==-1) return true;
	if (Tab->item(n, CWN)->text().toDouble() < Tab->item(m, CWN)->text().toDouble()) return true;
	return false;
}

bool sortByProgression(const QTableWidget *const Tab, const int n, const int m)
{
    if (n==-1) return false;
    if (m==-1) return true;
    if (Tab->item(n, CPN)->text().toInt() < Tab->item(m, CPN)->text().toInt()) return true;
    return false;
}
