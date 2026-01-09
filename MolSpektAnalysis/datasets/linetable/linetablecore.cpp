//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "linetablecore.h"
#include "utils.h"
#include "basedata.h"
#include "termenergy.h"
#include "Spektrum.h"
#include "elstate.h"
#include "molecule.h"
#include "vsolistelement.h"
#include "linetablesortfunctions.h"

#include <QTextStream>
#include <QPixmap>
#include <QPainter>
#include <QMessageBox>


LineTableCore::LineTableCore(LineTable* lt, Molecule* mol, QObject *parent) : QAbstractTableModel(parent), molecule(mol), lTab(lt)
{
	setColumnCount(TableNormCols);
}

LineTableCore::~LineTableCore()
{
}

QStringList LineTableCore::convertColumnVectorToHeaderStrings(const std::vector<TableCols>& headerVector)
{
	QStringList Result;
	for (auto it = headerVector.begin(); it != headerVector.end(); ++it) Result << convertHeaderEnumToHeaderString(*it);
	return Result;
}

QString LineTableCore::convertHeaderEnumToHeaderString(const TableCols enumValue)
{
	switch (enumValue)
	{
		case CPN:
			return "PN";
			break;
		case Cvs:
			return "v'";
			break;
		case CJs:
			return "J'";
			break;
		case Cvss:
			return "v''";
			break;
		case CJss:
			return "J''";
			break;
		case CF:
			return "FC";
			break;
		case CWN:
			return "wave number";
			break;
		case Cerr:
			return "error";
			break;
		case CIso:
			return "isotope";
			break;
		case CFile:
			return "file name";
			break;
		case CSNR:
			return "SNR";
			break;
		case CDev:
			return "Deviation";
			break;
		case CC:
			return "Comment";
			break;
		case CFCF:
			return "calc. rel. int.";
			break;
		case CEUp:
			return "E_UpTerm";
			break;
		case CEav:
			return "E_avarage";
			break;
		case CEUma:
			return "E_Up - E_av";
			break;
		case CEdJ:
			return "E/(J*(J+1))";
			break;
		case CCalc:
			return "Calc";
			break;
		case COmC:
			return "obs-calc";
			break;
		default:
			return "";
			break;
	}
}

std::vector<LineTableCore::TableCols> LineTableCore::convertHeaderStringsToColumnVector(const QStringList& headerStrings)
{
	std::vector<TableCols> Result;
	for (auto it = headerStrings.begin(); it != headerStrings.end(); ++it)
	{
		if (*it == "PN") Result.push_back(CPN);
		else if (*it == "v'") Result.push_back(Cvs);
		else if (*it == "J'") Result.push_back(CJs);
		else if (*it == "v''") Result.push_back(Cvss);
		else if (*it == "J''") Result.push_back(CJss);
		else if (*it == "FC") Result.push_back(CF);
		else if (*it == "wave number") Result.push_back(CWN);
		else if (*it == "error") Result.push_back(Cerr);
		else if (*it == "isotope") Result.push_back(CIso);
		else if (*it == "file name") Result.push_back(CFile);
		else if (*it == "SNR") Result.push_back(CSNR);
		else if (*it == "Deviation") Result.push_back(CDev);
		else if (*it == "Comment") Result.push_back(CC);
	}
	return Result;
}

BaseData * LineTableCore::convertToBaseData(const QStringList& L) const
{
	LineTableBaseData *data = new LineTableBaseData;
	int lc = L.count();
	IsoTab * Iso = nullptr;
	if (nullptr != molecule) Iso = molecule->getIso();
	for (int n=0; n < lc && n < lc; ++n)
	{
		switch(n)
		{
			case Cvs:
				data->vs = L[n].toInt();
				break;
			case CJs:
				data->Js = L[n].toInt();
				break;
			case Cvss:
				data->vss = L[n].toInt();
				break;
			case CJss:
				data->Jss = L[n].toInt();
				break;
			case CF:
				data->F = L[n].toInt();
				break;
			case CWN:
				data->waveNumber = L[n].toDouble();
				break;
			case Cerr:
				data->uncertainty = L[n].toDouble();
				break;
			case CIso:
				data->isotope = L[n].toInt();
				if (nullptr != Iso) data->IsoIcon = &Iso->IsoImage[data->isotope];
				break;
			case CFile:
				data->file = L[n];
				break;
			case CSNR:
				data->SNR = L[n].toDouble();
				break;
			case CDev:
				data->obsMinusCalc = L[n].toDouble();
				break;
			case CC:
				data->Comment = L[n];
				break;
			case CFCF:
				data->FCF = L[n].toDouble();
				break;
			case CEUp:
				data->upperEnergy = L[n].toDouble();
				break;
			case CEav:
				data->averageUpperEnergy = L[n].toDouble();
				break;
			case CEUma:
				data->diffToAverageEnergy = L[n].toDouble();
				break;
			case CEdJ:
				data->energyDiffToNextJ = L[n].toDouble();
				break;
			case CCalc:
				data->calculatedEnergy = L[n].toDouble();
				break;
			case COmC:
				data->diffToCalculatedEnergy = L[n].toDouble();
				break;
			default:
				// not possible
			break;
		}
	}
	return data;
}

LineTableBaseData * LineTableCore::convertStringArrayToLineTableBaseData(const QString *const array, const int length) const
{
	QStringList L;
	for (int l=0; l < length; ++l) L << array[l];
	return reinterpret_cast<LineTableBaseData*>(convertToBaseData(L));
}

void LineTableCore::setColumnCount(const int count)
{
	if (count < mColumns.size())
	{
		beginRemoveColumns(QModelIndex(), count, columnCount());
		mColumns.resize(count);
		endRemoveColumns();
	}
	else if (count > mColumns.size())
	{
		beginInsertColumns(QModelIndex(), columnCount(), count);
		if (mColumns.size() == CPN) mColumns.push_back(CPN);
		if (mColumns.size() == Cvs) mColumns.push_back(Cvs);
		if (mColumns.size() == CJs) mColumns.push_back(CJs);
		if (mColumns.size() == Cvss) mColumns.push_back(Cvss);
		if (mColumns.size() == CJss) mColumns.push_back(CJss);
		int F=0;
		if (mColumns.size() == CF)
		{
			for (auto it = mData.begin(); it != mData.end(); ++it) if (reinterpret_cast<LineTableBaseData*>(*it)->F >= 0)
			{
				F=1;
				break;
			}
			if (0 == F) mColumns.push_back(CF);
		}
		if (mColumns.size() + F == CWN) mColumns.push_back(CWN);
		if (mColumns.size() + F == Cerr) mColumns.push_back(Cerr);
		if (mColumns.size() + F == CIso) mColumns.push_back(CIso);
		if (mColumns.size() + F == CFile) mColumns.push_back(CFile);
		if (mColumns.size() + F == CSNR) mColumns.push_back(CSNR);
		if (mColumns.size() + F == CDev) mColumns.push_back(CDev);
		if (mColumns.size() + F == CC) mColumns.push_back(CC);
		if (mColumns.size() + F == CFCF) mColumns.push_back(CFCF);
		if (mColumns.size() + F == CEUp) mColumns.push_back(CEUp);
		if (mColumns.size() + F == CEav) mColumns.push_back(CEav);
		if (mColumns.size() + F == CEUma) mColumns.push_back(CEUma);
		if (mColumns.size() + F == CEdJ) mColumns.push_back(CEdJ);
		if (mColumns.size() + F == CCalc) mColumns.push_back(CCalc);
		if (mColumns.size() + F == COmC) mColumns.push_back(COmC);
		endInsertColumns();
	}
}

bool LineTableCore::readData(const int numLines, QString& Buffer, const bool FCA, int& MaxPN, const int d)
{
	QTextStream S2(&Buffer, QIODevice::ReadOnly);
	//printf("Ende: N=%d\n", N);
	QString B, Comm;
	QStringList L;
	int s, i;
	bool Success = true, COK;
    for (int n=0; n < numLines; ++n)
    {
		B = S2.readLine();
        L = B.split(QRegExp("\\s+"), Qt::SkipEmptyParts);
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[n]);
		if ((s = L.size()) < 8)
		{
			printf("Fehler beim Lesen, Zeile %d ist zu kurz!\n", n);
			Success = false;
			continue;
		}
		//printf("d=%d, n=%d, s=%d, NC=%d\n", d, n, s, TableNormCols);
		if (FCA)
		{
			if (s < CC)
			{
				if (s==8)
				{
					data->vs = L[1].toInt();
					data->Js = L[2].toInt();
					data->vss = L[3].toInt();
					data->Jss = L[4].toInt();
					data->waveNumber = L[5].toDouble();
					data->uncertainty = L[6].toDouble();
					data->F = L[7].toInt() + 1;
					data->isotope = 10;
					continue;
				}
				else if (s==11)
				{
					data->progressionNumber = L[0].toInt();
					data->vs = L[1].toInt();
					data->Js = L[2].toInt();
					data->vss = L[3].toInt();
					data->Jss = L[4].toInt();
					data->F = L[5].toInt();
					data->waveNumber = L[6].toDouble();
					data->uncertainty = L[7].toDouble();
					data->isotope = L[8].toInt();
					continue;
				}
				else
				{
					printf("Error, line %d is too short!", n);
					Success = false;
					continue;
				}
			}
			if ((i = L[0].toInt()) > MaxPN) MaxPN = i;
			data->progressionNumber = L[0].toInt();
			data->vs = L[1].toInt();
			data->Js = L[2].toInt();
			data->vss = L[3].toInt();
			data->Jss = L[4].toInt();
			data->F = L[5].toInt();
			data->waveNumber = L[6].toDouble();
			data->uncertainty = L[7].toDouble();
			data->isotope = L[8].toInt();
			data->file = L[9];
			for (L[i].toDouble(&COK); !COK && !L[i].isEmpty(); L[++i].toDouble(&COK))
			{
				L[CFile] += ' ' + L[i];
				if (i==s-1)
				{
					++i;
					break;
				}
			}
			if (i > CSNR) data->file = L[CFile];
			if (COK) data->SNR = L[i++].toDouble();
			if (i<s && L[i].toDouble() != 0.0) data->obsMinusCalc = L[i++].toDouble();
			if (s >= TableNormCols && i<s) data->Comment = L[i];
		}
		else
		{
			if (d==0 && (i = L[0].toInt()) > MaxPN)
			{
				MaxPN = i;
				data->progressionNumber = i;
			}
			data->vs = L[1].toInt();
			data->Js = L[2].toInt();
			data->vss = L[3].toInt();
			data->Jss = L[4].toInt();
			data->waveNumber = L[5].toDouble();
			data->uncertainty = L[6].toDouble();
			data->isotope = L[7].toInt();
			Comm = "";
            for (i=8; (i < s ? L[i] == "laser" || L[i] == "?" || L[i] == "weak" || L[i] == "strong"
						|| L[i] == "overlap" : false); i++)
			{
				if (i==8) Comm = L[i];
				else Comm += " " + L[i];
			}
			data->Comment = Comm;
			data->F = -1;
			if ((i < s ? L[i] != "0" : false)) data->file = L[i];
			if (i == 8 && 9 < s)
			{
				if (L[9] == "laser") data->Comment = Comm = L[9];
				else if (L[9] != "0") data->SNR = L[9].toDouble();
				if (s > 10)
				{
					if (L[10] == "laser") data->Comment = Comm = L[10];
					else if (L[10] != "0") data->obsMinusCalc = L[10].toDouble();
					if (s > 11) data->Comment = L[11];
				}
			}
		}
		lTab->setImported();
		//printf("n=%d, N=%d\n", n, N);
        if (molecule != 0)
        {
            QString CurPath = data->file, MolPath = molecule->getFileName();
            data->file = lTab->getAbsolutePath(CurPath, MolPath));
        }
    }
    return Success;
}

void LineTableCore::setRow(const QStringList& L, const int row)
{
	setRow(convertToBaseData(L), row);
}

void LineTableCore::writeData(QTextStream& S)
{
	int b;
	QString Buffer, LFile;
	for (auto it = mData.begin(); it != mData.end(); ++it)
    {
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(*it);
		S << ((b = data->progressionNumber) > 0 ? ("     " + QString::number(b)).right(5) : "    0");
		S << ("    " + QString::number(data->vs)).right(4);
		S << ("    " + QString::number(data->Js)).right(4);
		S << ("    " + QString::number(data->vss)).right(4);
		S << ("    " + QString::number(data->Jss)).right(4);
		S << ' ' << ("   " + QString::number(data->F)).right(3);
		Buffer = "               " + QString::number(data->waveNumber);
		if (Buffer.indexOf('.') == -1) Buffer += ".0000";
		Buffer = Buffer.right(15);
		if (Buffer[0] != ' ') S << " " << Buffer;
		else S << Buffer;
		if ((Buffer = "          " + QString::number(data->uncertainty, 'f', 8)).indexOf('.') == -1) Buffer += ".000";
		S << ' ' << Buffer.right(9);
		S << ("    " + QString::number(data->isotope)).right(4);
		if ((Buffer = data->file).isEmpty()) Buffer = LFile;
		else if (Buffer == "Aktuelle Markierungen")
		{
		    Buffer = lTab->getInfile();
			data->file = Buffer;
		}
		else LFile = Buffer;
        if (molecule != 0) Buffer = lTab->getRelativePath(Buffer, molecule->getFileName());
		S << " " << Buffer;
		if ((Buffer = QString::number(data->SNR, 'f', 6)).isEmpty()) Buffer = "0";
		S << ("      " + Buffer).right(7);
		if ((Buffer = QString::number(data->obsMinusCalc, 'f', 10)).isEmpty()) Buffer = "0";
		S << ("          " + Buffer).right(11);
		S << " " << data->Comment << char(10);
	}
}

void LineTableCore::WriteTFGS(QTextStream& S, vsOListElement *vsOList, const int vsO)
{
	int i, j, n, N = mData.size(), lvs, vs, lJs, Js, Jss, P, lP, PN, lPN, aIso, lIso, ivsO = vsO;
	IsoTab* Iso = molecule->getIso();
	QString Buffer, IBuff, F, lF;
	vsOListElement* CvsOElement = vsOList;
	int *SO = new int[N], *S1 = heapSort(sortIJvP);
	for (n=0; n<N; n++) SO[S1[n]] = n;
	for (i=0; i<N; ++i)
	{
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[SO[i]]);
		n = data->isotope;
		if (10 * (j = n / 10) != n) continue;
		if (--j < 0) continue;
		S << (IBuff = ("    " + QString::number(Iso->mNumIso1[j])).right(5)
				+ ("    " + QString::number(Iso->mNumIso2[j])).right(5));
		//printf("Nach1\n");
		vs = data->vs;
		Js = data->Js;
		Jss = data->Jss;
		P = abs(Js - Jss);
		F = data->file;
		PN = data->progressionNumber;
		aIso = n;
		if (F.left(5) == "laser") F = F.right(F.length() - 6);
		if (lvs != vs || lJs != Js || lIso != aIso)
		{
			while (CvsOElement->next != 0 ? CvsOElement->next->Iso <= aIso && CvsOElement->next->Js <= Js
					&& CvsOElement->vs <= vs : false) CvsOElement = CvsOElement->next;
			if (CvsOElement->Iso == aIso && CvsOElement->Js == Js && CvsOElement->vs == vs)
				CvsOElement->curMaxOffset = ivsO = (CvsOElement->curMaxOffset >= vsO ? CvsOElement->curMaxOffset + 100 : vsO + vs);
			else ivsO = vsO + vs;
			lvs = vs;
			lJs = Js;
			lF = F;
			lP = P;
			lPN = PN;
			lIso = aIso;
		}
		else if (F != lF || P != lP || PN != lPN)
		{
			if (CvsOElement->Iso == aIso && CvsOElement->Js == Js && CvsOElement->vs == vs)
				CvsOElement->curMaxOffset = ivsO = CvsOElement->curMaxOffset + 100;
			else if ((ivsO += 100) >= vsO + 1000)
			{
				vsOListElement *vsOBuff = new vsOListElement;
				vsOBuff->Iso = aIso;
				vsOBuff->vs = vs;
				vsOBuff->Js = Js;
				vsOBuff->curMaxOffset = ivsO;
				vsOBuff->next = CvsOElement->next;
				CvsOElement = CvsOElement->next = vsOBuff;
			}
			lF = F;
			lP = P;
			lPN = PN;
			lIso = aIso;
		}
		Buffer = "    " + QString::number(ivsO);
		S << Buffer.right(5);
		Buffer = "    " + QString::number(Js);
		S << Buffer.right(5);
		S << ("     " + QString::number(data->vss)).right(5);
		S << ("     " + QString::number(data->Jss)).right(5);
		S << IBuff << "    0    0    0    0\n";
		Buffer = "               ";
		Buffer += QString::number(data->waveNumber, 'f', 14);
		if (Buffer.indexOf(".") == -1) Buffer += ".00";
		S << Buffer.right(15);
		S << ("               " + QString::number(data->uncertainty, 'f', 14)).right(15) << "    2\n";
	}
	delete[] S1;
	delete[] SO;
}

void LineTableCore::Assign_vs(double ****const UD, const double AT, const int NumC, const int mvs, const int* const SO)
{
	double ***UE = 0;
	double dB, TE = 0.0;
	int i, nr = mData.size(), vs = 0, I=0, J=0;
	for (i=0; i<nr; ++i)
	{
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[SO[i]]);
		//printf("i=%d, nr=%d, CEav=%d, CCalc=%d, N=%d\n", i, nr, CEav, CCalc, Tab->columnCount());
		dB = data->upperEnergy;
		//printf("Z-\n");
		if (TE != dB)
		{
			//printf("Z0\n");
			TE = dB;
			//printf("Z1\n");
			if ((J = data->Js) != data->Jss || NumC == 1) UE = UD[0];
			else UE = UD[1];
			//printf("Z2\n");
			I = (data->isotope - 1) / 10;
			//printf("vs=%d, mvs=%d, I=%d, J=%d\n", vs, mvs, I, J);
			for (vs=0; (vs <= mvs ? UE[I][vs][J] < TE : false); vs++) ;
			//printf("Z4\n");
			if ((vs > 0 ? (vs <= mvs ? UE[I][vs][J] - TE > TE - UE[I][vs-1][J] : true)
				: false)) vs--;
			//printf("Z5\n");
			if (fabs(UE[I][vs][J] - TE) > AT) vs = -1;
		}
		//printf("Z6\n");
		data->vs = vs;
		//printf("Z7\n");
		if (vs >= 0)
		{
			//printf("Z8, I=%d, vs=%d, J=%d\n", I, vs, J);
			data->calculatedEnergy = UE[I][vs][J];
			//printf("Z9\n");
			data->diffToCalculatedEnergy = TE - UE[I][vs][J];
		}
	}
}

int LineTableCore::columnCount(const QModelIndex &parent) const
{
	if (parent.isValid()) return 0;
	return mColumns.size();
}

void LineTableCore::AssignFC(const int *const LO, const int *const XIT, const int* const EIT, const int ENv, const int ENJ, double ****EData, const int XNv,
							 const int XNJ, double ****XData, const int XNC, const int cTX[], const int cTE[], const int cTL[])
{
	int n, m, i, j, c, PN, RS, NR = mData.size(), nPN, I, vs, Js, NCol = columnCount(), vss[1000], Jss[1000], LI[1000], NI = molecule->getNumIso(), bc;
	double F[1000], BS = 0;
	for (n=0, m=1, PN = mData[LO[0]]->progressionNumber; m <= NR; ++m) if (m < NR && (nPN = mData[LO[m]]->progressionNumber) != PN)
	{
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[LO[n]]);
		I = (data->isotope - 1) / 10;
		vs = data->vs;
		Js = data->Js;
		if (I >= 0 && I < NI && XIT[I] >= 0 && EIT[I] >= 0 && vs >= 0 && vs < ENv && Js >= 0 && Js < ENJ && EData[0][EIT[I]][vs][Js] != 0.0)
		{
			for (i=n, j=0; i<m; ++i)
			{
				LineTableBaseData* data_i = reinterpret_cast<LineTableBaseData*>(mData[LO[i]]);
				vss[j] = data_i->vss;
				Jss[j] = data_i->Jss;
				LI[j] = LO[i];
				if (vss[j] >= 0 && vss[j] < XNv && Jss[j] >= 0 && Jss[j] < XNJ && XData[0][XIT[I]][vss[j]][Jss[j]] != 0.0) F[j++] = data_i->waveNumber;
			}
			if (j>0)
			{
				for (c=0; c < XNC; ++c)
				{
					for (i=0, RS = 0.0; i<j; ++i)
						RS += fabs(XData[cTX[c]][XIT[I]][vss[i]][Jss[i]] + F[i]
								   - EData[cTE[c]][EIT[I]][vs][Js]);
					if (RS < BS || c==0)
					{
						bc = c;
						BS = RS;
					}
				}
				for (i=n; i<m; ++i) reinterpret_cast<LineTableBaseData*>(mData[LO[i]])->F = cTL[bc];
				if (NCol > COmC)
				{
					for (i=0, RS = 0.0; i<j; ++i) RS += (F[i] += XData[cTX[c]][XIT[I]][vss[i]][Jss[i]]);
					RS /= j;
					for (i=0; i<j; ++i)
					{
						LineTableBaseData* data_i = reinterpret_cast<LineTableBaseData*>(mData[LO[i]]);
						data_i->upperEnergy = F[i];
						data_i->diffToAverageEnergy = F[i] - RS;
						data_i->diffToCalculatedEnergy = F[i] - EData[cTE[c]][EIT[I]][vs][Js];
					}
					for (i=n; i<m; ++i)
					{
						LineTableBaseData* data_i = reinterpret_cast<LineTableBaseData*>(mData[LO[i]]);
						data_i->averageUpperEnergy = RS;
						data_i->energyDiffToNextJ = RS / (Js * (Js + 1));
						data_i->calculatedEnergy = EData[cTE[c]][EIT[I]][vs][Js];
					}
				}
			}
			else for (i=n; i<m; ++i) reinterpret_cast<LineTableBaseData*>(mData[i])->F = -1;
		}
		else for (i=n; i<m; i++) reinterpret_cast<LineTableBaseData*>(mData[i])->F = -1;
		n=m;
		PN = nPN;
	}
}

int LineTableCore::columnCount(const QModelIndex &parent) const
{
	if (parent.isValid()) return 0;
	return mColumns.size();
}

QVariant LineTableCore::data(const QModelIndex &index, int role) const
{
	int row = index.row(), column = index.column(), numDigits;
	if (0 > row || row >= mData.size() || 0 > column || mColumns.size() <= column) return QVariant();
	if (role == Qt::DecorationRole)
	{
		if ((column == 0 || mColumns[column] == CIso) && nullptr != molecule)
		{
			if (nullptr == mData[row]->IsoIcon)
			{
				int n = mData[row]->isotope / 10 - 1;
				IsoTab* Iso = molecule->getIso();
				if (0 <= n && n < Iso->numIso) mData[row]->IsoIcon = &Iso->IsoImage[n];
			}
			return *mData[row]->IsoIcon;
		}
		if (row > mData.size() - lTab->getNumNewLines()) return *NewPix;
	}
	if (role == Qt::TextAlignmentRole)
	{
		if (mColumns[column] == CFile || mColumns[column] == CC) return Qt::AlignLeft;
		return Qt::AlignRight;
	}
	if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();
	if (mColumns[column] == CWN || mColumns[column] == Cerr || mColumns[column] == CDev || mColumns[column] == CEUp || mColumns[column] == CEav || mColumns[column] == CEUma
		 || mColumns[column] == CEdJ || mColumns[column] == CCalc || mColumns[column] == COmC) numDigits = getNumDecimalPlaces(mData[row]->uncertainty);
	switch (mColumns[column])
	{
		case CPN:
			return mData[row]->progressionNumber;
		case Cvs:
			return reinterpret_cast<LineTableBaseData*>(mData[row])->vs;
		case CJs:
			return mData[row]->Js;
		case Cvss:
			return reinterpret_cast<LineTableBaseData*>(mData[row])->vss;
		case CJss:
			return reinterpret_cast<LineTableBaseData*>(mData[row])->Jss;
		case CF:
			return reinterpret_cast<LineTableBaseData*>(mData[row])->F;
		case CWN:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->waveNumber, 'f', numDigits);
		case Cerr:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->uncertainty, 'f', numDigits);
		case CIso:
			return mData[row]->isotope;
		case CFile:
			return mData[row]->file;
		case CSNR:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->SNR, 'f', 3);
		case CDev:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->obsMinusCalc, 'f', numDigits);
		case CC:
			return reinterpret_cast<LineTableBaseData*>(mData[row])->Comment;
		case CFCF:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->FCF, 'f', 2);
		case CEUp:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->upperEnergy, 'f', numDigits);
		case CEav:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->averageUpperEnergy, 'f', numDigits);
		case CEUma:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->diffToAverageEnergy, 'f', numDigits);
		case CEdJ:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->energyDiffToNextJ, 'f', numDigits);
		case CCalc:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->calculatedEnergy, 'f', numDigits);
		case COmC:
			return QString::number(reinterpret_cast<LineTableBaseData*>(mData[row])->diffToCalculatedEnergy, 'f', numDigits);
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

bool LineTableCore::setData(const QModelIndex& index, const QVariant& value, int role)
{
	if (role != Qt::EditRole) return false;
	if (index.column() < 0 || index.column() >= mColumns.size() || index.row() < 0 || index.row() > mData.size()) return false;
	switch (mColumns[index.column()])
	{
		case CPN:
			mData[index.row()]->isotope = value.toInt();
			break;
		case Cvs:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->vs = value.toInt();
			break;
		case CJs:
			mData[index.row()]->Js = value.toInt();
			break;
		case Cvss:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->vss = value.toInt();
			break;
		case CJss:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->Jss = value.toInt();
			break;
		case CF:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->F = value.toInt();
			break;
		case CWN:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->waveNumber = value.toDouble();
			break;
		case Cerr:
			mData[index.row()]->uncertainty = value.toDouble();
			break;
		case CIso:
			mData[index.row()]->isotope = value.toInt();
			break;
		case CFile:
			mData[index.row()]->file = value.toString();
			break;
		case CSNR:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->SNR = value.toDouble();
			break;
		case CDev:
			mData[index.row()]->obsMinusCalc = value.toDouble();
			break;
		case CC:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->Comment = value.toString();
			break;
		case CFCF:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->FCF = value.toDouble();
			break;
		case CEUp:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->upperEnergy = value.toDouble();
			break;
		case CEav:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->averageUpperEnergy = value.toDouble();
			break;
		case CEUma:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->diffToAverageEnergy = value.toDouble();
			break;
		case CEdJ:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->energyDiffToNextJ = value.toDouble();
			break;
		case CCalc:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->calculatedEnergy = value.toDouble();
			break;
		case COmC:
			reinterpret_cast<LineTableBaseData*>(mData[index.row()])->diffToCalculatedEnergy = value.toDouble();
			break;
		default:
			return false;
			break;
	}
	QVector<int> roles;
	roles.append(role);
	emit dataChanged(index, index, roles);
	return true;
}

QVariant LineTableCore::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole) return QVariant();
	if (orientation == Qt::Vertical) return section;
	if (0 <= section && mColumns.size() > section) return convertHeaderEnumToHeaderString(mColumns[section]);
	return QVariant();
}

int LineTableCore::getMarkedLine(const Marker& marker, const QString SFN) const
{
	int n, RC = mData.size();
	for (n=0; n < RC; ++n)
	{
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[n]);
		if (data->file.indexOf(SFN) != -1 && marker.Iso == (data->isotope - 1) / 10 && marker.FC == data->F && marker.vs == data->vs && marker.Js == data->Js
			&& marker.vss == data->vss && marker.Jss == data->Jss && fabs(marker.Line[0] - data->waveNumber) <= 1e-3) return n;
	}
}

void LineTableCore::AddMarked(Marker ** Lines, const Marker* const LaserLine, const int AnzahlMarker, int& Offset, int &NpProg, const int MaxPN, const QString& SpektFile)
{
	for (int n=0; n < AnzahlMarker; ++n)
	{
		//printf("i=%d, n=%d, N=%d\n", Offset, n, N);
		if (n>0 && (Lines[n]->vs != Lines[n-1]->vs || Lines[n]->Iso != Lines[n-1]->Iso
				  || Lines[n]->Js != Lines[n-1]->Js)) ++NpProg;
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[Offset + n]);
	    data->progressionNumber = MaxPN + NpProg;
		if (Lines[n]->DisplayData)
		{
			data->vs = Lines[n]->vs;
			data->Js = Lines[n]->Js;
			data->vss = Lines[n]->vss;
			data->Jss = Lines[n]->Jss;
			data->isotope = 10 * (Lines[n]->Iso + 1);
			data->F = Lines[n]->FC;
		}
		else
		{
			data->vs = -1;
			data->Js = -1;
			data->vss = -1;
			data->Jss = -1;
			data->F = -1;
			data->isotope = 10;
		}
		data->waveNumber = Lines[n]->Line[0];
		data->uncertainty = (Lines[n]->DisplayData && Lines[n]->uncertainty > 0.0 ? Lines[n]->uncertainty : 0.005);
    	data->file = SpektFile;
		data->SNR = Lines[n]->SNR;
		data->obsMinusCalc = Lines[n]->DD;
		QString Buffer;
		if (Lines[n]->satellite) Buffer = "satellite";
		else if (Lines[n] == LaserLine) Buffer = "laser";
		else Buffer = "";
		if (Lines[n]->overlap) Buffer += (Buffer.isEmpty() ? "overlap" : ",overlap");
		data->Comment = Buffer;
    }
}

void LineTableCore::ShowUpTerm(const int* const SO, const int MCT, const int* const CT, const int* const IsoT, double ****ELU, const int Mv, const int MJ, double ****UT,
							   const int mvs, const int Jeo, const int* const UIsoT, const float S, const int UNC, const int UMCT, const int* const UCT)
{
	int i, j, v, c, I, n, N = mData.size(), NI = molecule->getNumIso();
	int **Z = CreateInt(N, 6), li[2] = {-1, -1};
	double E, ES, dE, *T = new double[N];
	QString Buffer, SpektFile;
	for (i=0; i<N; ++i)
    {
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[SO[i]]);
		if ((Buffer = mData[SO[i]]->file).isEmpty()) mData[SO[i]]->file = SpektFile;
		else if (Buffer.left(6) == " laser")
		{
		    for (j=6; Buffer[j] == ' '; ++j) ;
	    	SpektFile = Buffer.right(Buffer.length() - j + 1);
		}
		else SpektFile = Buffer;
		j = data->Jss;
		v = data->vss;
		Z[i][0] = data->isotope;
		I = (Z[i][0] - 1) / 10;
		Z[i][1] = data->vs;
		Z[i][2] = data->Js;
		Z[i][3] = abs(Z[i][2] - j);
		Z[i][4] = j;
		Z[i][5] = data->progressionNumber;
		c = data->F;
		c = (c >= 0 && c <= MCT ? (CT[c] >= 0 ? CT[c] : 0) : 0);
		//printf("i=%d, I=%d, v=%d, J=%d\n", i, I, v, j);
		E = ((I >= 0 && I < NI && IsoT[I] >= 0 && v >= 0 && v <= Mv && j >= 0
				&& j <= MJ && ELU[c][IsoT[I]][v][j] != 0.0) ?
					data->waveNumber + ELU[c][IsoT[I]][v][j] : 0.0);
		data->upperEnergy = E;
		T[i] = E;
    }
	//printf("Vor Schleife2\n");
    for (E=ES=0.0, i=n=0; n<=N; n++)
    {
		//printf("E=%f, n=%d, N=%d, i=%d\n", E, n, N, i);
		if (n<N && n>i && (Z[i][0] != Z[n][0] || Z[i][1] != Z[n][1] || Z[i][2] != Z[n][2]
				  || Z[i][3] != Z[n][3] || Z[i][5] != Z[n][5]))
		{
	    	E /= ES;
	    	j = li[Z[n-1][3]];
	    	if (j>0 && Z[j][0] == Z[n-1][0] && Z[j][1] == Z[n-1][1] && Z[j][2] == Z[n-1][2] - 1) dE = (E - T[j]) / (2 * Z[n-1][2]);
	    	else dE = 0.0;
	    	for (j=i; j<n; ++j)
	    	{
				LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(mData[SO[j]]);
				/*if (n==1482)
					printf("Z[%d][0]=%d, Z[%d][1]=%d, Z[%d][2]=%d, Z[%d][4]=%d\n",
						   j, Z[j][0], j, Z[j][1], j, Z[j][2], j, Z[j][4]);*/
				data->averageUpperEnergy = E;
				data->diffToAverageEnergy = T[j] - E;
				if (dE > 0.0) data->energyDiffToNextJ = dE;
				I = (Z[j][0] - 1) / 10;
				if (T[j] != 0.0 && UT != 0 && Z[j][1] >= 0 && Z[j][1] <= mvs
					&& Z[j][2] <= Jeo && UIsoT[I] >= 0)
				{
					//printf("Beginn Neu, j=%d\n", j);
					v = data->vss;

					//printf("I=%d, Z[j][0]=%d, Z[j][1]=%d, Z[j][2]=%d, Z[j][4]=%d, v=%d\n",
						//   I, Z[j][0], Z[j][1], Z[j][2], Z[j][4], v);
					//printf("ELU=%f, eT=%f\n", ELU[I][v][Z[j][4]], eT[I][Z[j][1]][Z[j][2]]);
					if (S==0.0) c = (Z[j][2] != Z[j][4] || UNC == 1 ? 0 : 1);
					else
					{
						c = reinterpret_cast<LineTableBaseData*>(mData[SO[i]])->F;
						c = (c >= 0 && c <= UMCT ? (UCT[c] >= 0 ? UCT[c] : 0) : 0);
					}
					double Calc = UT[c][UIsoT[I]][Z[j][1]][Z[j][2]];
					//printf("E=%f, Calc=%f\n", E, Calc);
					data->calculatedEnergy = Calc;
					data->diffToCalculatedEnergy = E - Calc;
					//printf("Ende Neu\n");
				}
	    	}
	    	//if (n==1482) printf("n=%d, i=%d\n", n, i);
	    	if (n < N)
	    	{
				li[Z[n-1][3]] = n - 1;
				i = n;
				T[n-1] = E;
				E = ES = 0.0;
	    	}
		}
		if (n<N)
		{
			double EF = mData[SO[n]]->uncertainty;
			EF = 1.0 / (EF * EF);
			E += T[n] * EF;
			ES += EF;
		}
    }
    delete[] T;
}

void LineTableCore::ShowUpTermtable(const int *const SO, int &n, QString ** Data)
{
	int i, j=0, k, N=0, NSO = mData.size(), Iso = 0, vs = 0, iJs, Js = 0, Jss = 0, iJss, AD, vss;
	double EAV, SqSum = 0.0, Diff, MinE, MaxE, TE;
	QString File;
	for (i=0; i < NSO; ++i) if (reinterpret_cast<LineTableBaseData*>(SO[i])->averageUpperEnergy != EAV)
    {
		EAV = reinterpret_cast<LineTableBaseData*>(SO[i])->averageUpperEnergy;
		++N;
    }
    Data = CreateQString(NSO, 10);
    for (n=i=0; i <= NSO; ++i)
    {
		LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(SO[i]);
		if (i < NSO || data->averageUpperEnergy != EAV)
		{
	    	if (i > 0)
	    	{
				Data[n][0] = QString::number(Iso);
				Data[n][1] = QString::number(vs);
				Data[n][2] = QString::number(Js);
				if (Jss == Js) Data[n][3] = "f";
				else Data[n][3] = "e";
				Data[n][4] = QString::number(EAV, 'f', 7);
				if (j > 1)
					Data[n][5] = QString::number(sqrt(SqSum / (j * (j - 1))), 'g', 6);
				Data[n][6] = QString::number(j);
				if (Jss != Js) Data[n][7] = QString::number(AD);
				Data[n][8] = QString::number(MinE, 'g', 9);
				Data[n++][9] = QString::number(MaxE, 'g', 9);
	    	}
	    	if (i < NSO)
	    	{
				EAV = data->averageUpperEnergy;
				iJs = Js = data->Js;
				Jss = data->Jss;
				Iso = data->isotope;
				vs = data->vs;
				j = 0;
				SqSum = 0.0;
				MinE = MaxE = EAV;
				AD = 0;
	    	}
		}
		if (i < NSO)
		{
			Diff = EAV - (TE = data->upperEnergy);
			if (TE < MinE) MinE = TE;
			else if (TE > MaxE) MaxE = TE;
			SqSum += Diff * Diff;
	    	++j;
			if (Jss != Js)
			{
				iJss = 2 * iJs - Jss;
				vss = data->vss;
				File = data->file;
				//printf("Jss=%s, vss=%s, vs=%s, Js=%s, File=%s\n",
					 //  Jss.ascii(), vss.ascii(), vs.ascii(), Js.ascii(), File.ascii());
				for (k=i-1; k>=0 && mData[SO[k]]->Js == Js && reinterpret_cast<LineTableBaseData*>(mData[SO[k]])->vs == vs && mData[SO[k]]->isotope == Iso; --k)
				{
					LineTableBaseData* data = reinterpret_cast<LineTableBaseData*>(SO[k]);
					if (data->vss == vss && data->Jss == 2 * Js - iJss && data->file == File) ++AD;
				}
			}
			//printf("k=%d, i=%d, AD=%d\n", k, i, AD);
		}
    }
}

void LineTableCore::RemoveDoubled(const int *const sortArray)
{
	int i, n1, n2, m1, m2, nr = mData.size();
	QString F1, F2, Comment;
	for (i=1, F2 = mData[sortArray[0]]->file; i < nr; ++i)
	{
		LineTableBaseData* data0 = reinterpret_cast<LineTableBaseData*>(mData[sortArray[i-1]]);
		LineTableBaseData* data1 = reinterpret_cast<LineTableBaseData*>(mData[sortArray[i]]);
		F1 = F2;
		F2 = data1->file;
		//printf("I=%d, vs=%d, Js=%d, vss=%d, Jss=%d, WN=%f\n", Tab->item(S[i], CIso)->text().toInt(), Tab->item(S[i], Cvs)->text().toInt(),
			//   Tab->item(S[i], CJs)->text().toInt(), Tab->item(S[i], Cvss)->text().toInt(), Tab->item(S[i], CJss)->text().toInt(),
			  // Tab->item(S[i], CWN)->text().toDouble());
		if (data1->vss == data0->vss && data1->Jss == data0->Jss && data1->Js == data0->Js)
		{
			n1 = F1.lastIndexOf(QRegExp("[\\/]"));
			n2 = F2.lastIndexOf(QRegExp("[\\/]"));
			m1 = ((m1 = F1.lastIndexOf('.')) != -1 ? m1 : F1.length());
			m2 = ((m2 = F2.lastIndexOf('.')) != -1 ? m2 : F2.length());
			if (F1.mid(n1 + 1, m1 - n1 - 1) == F2.mid(n2 + 1, m2 - n2 - 1) && data1->isotope == data0->isotope && fabs(data1->waveNumber - data0->waveNumber) < 1e-4)
			{
				mData[sortArray[i]] = nullptr;
				if (F1 != F2)
				{
					QFile File1(F1), File2(F2);
					if (!File1.exists() && File2.exists()) data0->file = F2;
				}
				Comment = data1->Comment;
				if (data0->Comment.length() < Comment.length()) data0->Comment = Comment;
				if (data0->uncertainty < data1->uncertainty) data0->uncertainty = data1->uncertainty;
			}
		}
	}
}

void LineTableCore::setData(const std::vector<LineTableCore::TableCols>& columns, const std::vector<BaseData *> data)
{
	if (columns.size() > mColumns.size())
	{
		beginInsertColumns(QModelIndex(), 0, columns.size());
		mColumns = columns;
		endInsertColumns();
	}
	else if (columns.size() < mColumns.size())
	{
		beginRemoveColumns(QModelIndex(), 0, mColumns.size());
		mColumns = columns;
		endRemoveColumns();
	}
	else mColumns = columns;
	emit headerDataChanged(Qt::Horizontal, 0, mColumns.size());
	if (data.size() > mData.size())
	{
		beginInsertRows(QModelIndex(), 0, data.size());
		mData = data;
		endInsertRows();
	}
	else if (data.size() < mData.size())
	{
		beginRemoveRows(QModelIndex(), 0, mData.size());
		mData = data;
		endRemoveRows();
	}
	else mData = data;
	EmitDataChanged(0, mData.size(), 0, mColumns.size());
}

void LineTableCore::setUncertainty(const int row, const double Error, const double OvError)
{
	QString Buffer = reinterpret_cast<LineTableBaseData*>(mData[row])->Comment;
	if (Buffer.indexOf("overlap") > -1 || Buffer.indexOf("laser") > -1 || Buffer.indexOf("weak") > -1)
	{
		if (mData[row]->uncertainty < OvError) mData[row]->uncertainty = OvError;
	}
	else if (mData[row]->uncertainty > Error) mData[row]->uncertainty = Error;
}

void set_vs(const int row, const int vs)
{
	reinterpret_cast<LineTableBaseData*>(mData[row])->vs = vs;
}

void set_vss(const int row, const int vss)
{
	reinterpret_cast<LineTableBaseData*>(mData[row])->vss = vss;
}
