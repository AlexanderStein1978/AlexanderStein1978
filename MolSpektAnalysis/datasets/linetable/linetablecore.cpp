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
	mColumns.push_back(CPN);
	mColumns.push_back(Cvs);
	mColumns.push_back(CJs);
	mColumns.push_back(Cvss);
	mColumns.push_back(CJss);
	mColumns.push_back(CWN);
	mColumns.push_back(Cerr);
	mColumns.push_back(CIso);
	mColumns.push_back(CFile);
	mColumns.push_back(CSNR);
	mColumns.push_back(CDev);
	mColumns.push_back(CC);
}

LineTableCore::~LineTableCore()
{
}

QStringList LineTableCore::convertColumnVectorToHeaderStrings(const std::vector<TableCols>& headerVector)
{
	QStringList Result;
	for (auto it = headerVector.begin(); it != headerVector.end(); ++it)
	{
		switch(*it)
		{
			case CPN:
				Result.push_back("PN");
				break;
			case Cvs:
				Result.push_back("v'");
				break;
			case CJs:
				Result.push_back("J'");
				break;
			case Cvss:
				Result.push_back("v''");
				break;
			case CJss:
				Result.push_back("J''");
				break;
			case CF:
				Result.push_back("FC");
				break;
			case CWN:
				Result.push_back("wave number");
				break;
			case Cerr:
				Result.push_back("error");
				break;
			case CIso:
				Result.push_back("isotope");
				break;
			case CFile:
				Result.push_back("file name");
				break;
			case CSNR:
				Result.push_back("SNR");
				break;
			case CDev:
				Result.push_back("Deviation");
				break;
			case CC:
				Result.push_back("Comment");
				break;
			default:
				// do nothing for now
				break;
		}
	}
	return Result;
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
	BaseData *data = new BaseData;
	int lc = L.count();
	IsoTab * Iso = nullptr;
	if (nullptr != molecule) Iso = molecule->getIso();
	for (int n=0; n < lc && n <= fdcLineElState; n++)
	{
		switch(n)
		{
			case fdcIso:
				data->isotope = static_cast<char>(L[n].toInt());
				if (nullptr != Iso) data->IsoIcon = &Iso->IsoImage[data->isotope];
				break;
			case fdcv:
				data->v = static_cast<ushort>(L[n].toInt());
				break;
			case fdcJ:
				data->J = static_cast<ushort>(L[n].toInt());
				break;
			case fdcvs:
				data->vs = L[n].toStdString();
				break;
			case fdcJs:
				data->Js = static_cast<ushort>(L[n].toInt());
				break;
			case fdcSource:
				data->source = L[n].toStdString();
				break;
			case fdcProg:
				data->prog = L[n].toInt();
				break;
			case fdcFile:
				data->file = L[n].toStdString();
				break;
			case fdcEnergy:
				data->energy = L[n].toDouble();
				break;
			case fdcUncert:
				data->uncert = L[n].toDouble();
				break;
			case fdcObsCalc:
				data->obs_calc = L[n].toDouble();
				break;
			case fdcDevR:
				data->devR = L[n].toFloat();
				break;
			case fdcLineElState:
				data->secondState = L[n].toStdString();
				break;
			default:
				// not possible
			break;
		}
	}
	return data;
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
	if (role == Qt::DecorationRole)
	{
		if (column == fdcIso && nullptr != molecule)
		{
			if (nullptr == mData[row]->IsoIcon)
			{
				int n = mData[row]->isotope / 10 - 1;
				IsoTab* Iso = molecule->getIso();
				if (0 <= n && n < Iso->numIso) mData[row]->IsoIcon = &Iso->IsoImage[n];
			}
			return *mData[row]->IsoIcon;
		}
		if (row >= NSources) return *NewPix;
	}
	if (role == Qt::TextAlignmentRole)
	{
		if (column == fdcSource || column == fdcFile || column == fdcLineElState) return Qt::AlignLeft;
		return Qt::AlignRight;
	}
	if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();
	if (column == fdcEnergy || column == fdcUncert || column == fdcObsCalc) numDigits = getNumDecimalPlaces(mData[row]->uncert);
	switch (column)
	{
		case fdcIso:
			return mData[row]->isotope;
		case fdcv:
			return mData[row]->v;
		case fdcJ:
			return mData[row]->J;
		case fdcvs:
			return mData[row]->vs.c_str();
		case fdcJs:
			return mData[row]->Js;
		case fdcSource:
			return mData[row]->source.c_str();
		case fdcProg:
			return mData[row]->prog;
		case fdcFile:
			return mData[row]->file.c_str();
		case fdcEnergy:
			if (!RWError.isEmpty())
			{
				numDigits = -int(floor(log10(mData[row]->uncert)));
				if (numDigits < 4) numDigits = 4;
			}
			return QString::number(mData[row]->energy, 'f', numDigits);
		case fdcUncert:
			return QString::number(mData[row]->uncert, 'f', numDigits);
		case fdcObsCalc:
			return QString::number(mData[row]->obs_calc, 'f', numDigits);
		case fdcDevR:
			return QString::number(mData[row]->devR, 'f', 3);
		case fdcLineElState:
			return mData[row]->secondState.c_str();
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

bool LineTableCore::setData(const QModelIndex& index, const QVariant& value, int role)
{
	if (role != Qt::EditRole) return false;
	switch (index.column())
	{
		case fdcIso:
			mData[index.row()]->isotope = static_cast<ushort>(value.toUInt());
			break;
		case fdcv:
			mData[index.row()]->v = static_cast<ushort>(value.toUInt());
			break;
		case fdcJ:
			mData[index.row()]->J = static_cast<ushort>(value.toUInt());
			break;
		case fdcvs:
			mData[index.row()]->vs = value.toString().toStdString();
			break;
		case fdcJs:
			mData[index.row()]->Js = static_cast<ushort>(value.toUInt());
			break;
		case fdcSource:
			mData[index.row()]->source = value.toString().toStdString();
			break;
		case fdcProg:
			mData[index.row()]->prog = value.toInt();
			break;
		case fdcFile:
			mData[index.row()]->file = value.toString().toStdString();
			break;
		case fdcEnergy:
			mData[index.row()]->energy = value.toDouble();
			break;
		case fdcUncert:
			mData[index.row()]->uncert = value.toDouble();
			break;
		case fdcObsCalc:
			mData[index.row()]->obs_calc = value.toDouble();
			break;
		case fdcDevR:
			mData[index.row()]->devR = value.toFloat();
			break;
		case fdcLineElState:
			mData[index.row()]->secondState = value.toString().toStdString();
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
	TableCols column = mColumns[section];
	switch (column)
	{
		case CPN:
			return "PN";
		case Cvs:
			return "v'";
		case CJs:
			return "J'";
		case Cvss:
			return "v''";
		case CJss:
			return "J''";
		case CF:
			return "FC";
		case CWN:
			return "wave number";
		case Cerr:
			return "error";
		case CIso:
			return "isotope";
		case CFile:
			return "file name";
		case CSNR:
			return "SNR";
		case CDev:
			return "Deviation";
		case CC:
			return "Comment";
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

