//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "fitdatacore.h"
#include "utils.h"
#include "fitdatabasedata.h"
#include "termenergy.h"
#include "Spektrum.h"
#include "elstate.h"
#include "molecule.h"

#include <QTextStream>
#include <QPixmap>
#include <QPainter>
#include <QMessageBox>


FitDataCore::FitDataCore(Molecule* mol, QObject *parent) : TableViewWindowCore(mol, parent, QRegExp("SourceOffsets:|Begin ResidualFit"))
{
}

FitDataCore::~FitDataCore()
{
}

FitDataBaseData * FitDataCore::convertToFitDataCoreBaseData(const QStringList& L) const
{
	FitDataBaseData *data = new FitDataBaseData;
	int lc = L.count();
	IsoTab * Iso = nullptr;
	if (nullptr != molecule) Iso = molecule->getIso();
	for (int n=0; n < lc && n <= fdcLineElState; n++)
	{
		switch(n)
		{
			case fdcIso:
				data->isotope = L[n].toInt();
				if (nullptr != Iso) data->IsoIcon = &Iso->IsoImage[data->isotope];
				break;
			case fdcv:
				data->v = L[n].toInt();
				break;
			case fdcJ:
				data->J = L[n].toInt();
				break;
			case fdcvs:
				data->vs = L[n];
				break;
			case fdcJs:
				data->Js = L[n].toInt();
				break;
			case fdcSource:
				data->source = L[n];
				break;
			case fdcProg:
				data->progressionNumber = L[n].toInt();
				break;
			case fdcFile:
				data->file = L[n];
				break;
			case fdcEnergy:
				data->energy = L[n].toDouble();
				break;
			case fdcUncert:
				data->uncertainty = L[n].toDouble();
				break;
			case fdcObsCalc:
				data->obsMinusCalc = L[n].toDouble();
				break;
			case fdcDevR:
				data->devR = L[n].toDouble();
				break;
			case fdcLineElState:
				data->secondState = L[n];
				break;
			default:
				// not possible
			break;
		}
	}
	return data;
}

QString FitDataCore::readData(QTextStream& S)
{
	QString Buffer;
	QString Spacer = "\t";
	while (S.readLine().left(15) != "Column titles: ") ;
	beginInsertRows(QModelIndex(), -1, -1);
	while(!S.atEnd())
	{
		Buffer = S.readLine();
		if (Buffer.indexOf(mStartSpecialPart) >= 0) return Buffer;
		const QStringList L = Buffer.split(Spacer);
		if (L.size() < 10) continue;
		mData.push_back(convertToBaseData(L));
	}
	endInsertRows();
	beginInsertRows(QModelIndex(), 0, mData.size() - 2);
	endInsertRows();
	return "";
}

int FitDataCore::columnCount(const QModelIndex &parent) const
{
	if (parent.isValid()) return 0;
	return 13;
}

QVariant FitDataCore::data(const QModelIndex &index, int role) const
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
	FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[row]);
	if (column == fdcEnergy || column == fdcUncert || column == fdcObsCalc) numDigits = getNumDecimalPlaces(element->uncertainty);
	switch (column)
	{
		case fdcIso:
			return element->isotope;
		case fdcv:
			return element->v;
		case fdcJ:
			return element->J;
		case fdcvs:
			return element->vs;
		case fdcJs:
			return element->Js;
		case fdcSource:
			return element->source;
		case fdcProg:
			return element->progressionNumber;
		case fdcFile:
			return element->file;
		case fdcEnergy:
			if (!RWError.isEmpty())
			{
				numDigits = -int(floor(log10(element->uncertainty)));
				if (numDigits < 4) numDigits = 4;
			}
			return QString::number(element->energy, 'f', numDigits);
		case fdcUncert:
			return QString::number(element->uncertainty, 'f', numDigits);
		case fdcObsCalc:
			return QString::number(element->obsMinusCalc, 'f', numDigits);
		case fdcDevR:
			return QString::number(element->devR, 'f', 3);
		case fdcLineElState:
			return element->secondState;
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

bool FitDataCore::setData(const QModelIndex& index, const QVariant& value, int role)
{
	if (role != Qt::EditRole) return false;
	FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[index.row()]);
	switch (index.column())
	{
		case fdcIso:
			element->isotope = value.toUInt();
			break;
		case fdcv:
			element->v = value.toUInt();
			break;
		case fdcJ:
			element->J = value.toUInt();
			break;
		case fdcvs:
			element->vs = value.toString();
			break;
		case fdcJs:
			element->Js = value.toUInt();
			break;
		case fdcSource:
			element->source = value.toString();
			break;
		case fdcProg:
			element->progressionNumber = value.toInt();
			break;
		case fdcFile:
			element->file = value.toString();
			break;
		case fdcEnergy:
			element->energy = value.toDouble();
			break;
		case fdcUncert:
			element->uncertainty = value.toDouble();
			break;
		case fdcObsCalc:
			element->obsMinusCalc = value.toDouble();
			break;
		case fdcDevR:
			element->devR = value.toDouble();
			break;
		case fdcLineElState:
			element->secondState = value.toString();
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

QVariant FitDataCore::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole) return QVariant();
	if (orientation == Qt::Vertical) return section;
	switch (section)
	{
		case fdcIso:
			return "Iso";
		case fdcv:
			return "v";
		case fdcJ:
			return "J";
		case fdcvs:
			return "v'";
		case fdcJs:
			return "J'";
		case fdcSource:
			return "Source";
		case fdcProg:
			return "Prog";
		case fdcFile:
			return "File";
		case fdcEnergy:
			return "Energy";
		case fdcUncert:
			return "uncert.";
		case fdcObsCalc:
			return "obs-calc";
		case fdcDevR:
			return "DevR";
		case fdcLineElState:
			if (RWError.isEmpty()) return "Second State";
			return RWError;
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

int FitDataCore::addMarkedLevel(TermEnergy& TE, Spektrum* Source)
{
	int R = mData.size();
	beginInsertRows(QModelIndex(), R-1, R-1);
	FitDataBaseData* element = new FitDataBaseData;
	element->isotope = static_cast<char>(TE.Iso);
	element->v = TE.v;
	element->J = TE.J;
	element->vs = (TE.FC == -1 ? "TE" : QString::number(-1 - TE.FC));
	element->source = Source->getName();
	element->file = Source->getFileName();
	element->energy = TE.E;
	element->uncertainty = TE.err;
	element->obsMinusCalc = TE.dev;
	element->devR = TE.DevR;
	element->secondState = TE.State->getName();
	mData.push_back(element);
	endInsertRows();
    return R;
}

int FitDataCore::addRow(const int cr)
{
	beginInsertRows(QModelIndex(), cr - 1, cr - 1);
	int r=0, nr = mData.size();
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); r < cr; ++r) ++it;
	mData.insert(it, new FitDataBaseData);
	endInsertRows();
	return nr;
}

void FitDataCore::addRow(const QStringList& L)
{
	TableViewWindowCore::addRow(convertToFitDataCoreBaseData(L));
}

void FitDataCore::addData(const int i_numLines, int *const i_Lines, const FitDataCore& data)
{
	int dataSize = mData.size();
	beginInsertRows(QModelIndex(), dataSize, dataSize + data.mData.size() - 2);
	auto it = data.mData.begin();
	int r=0;
	for (int i=0; i < i_numLines; ++i)
	{
		while (r < i_Lines[i])
		{
			++r;
			++it;
		}
		while (r > i_Lines[i])
		{
			--r;
			--it;
		}
		mData.insert(mData.end(), it, data.mData.end());
	}
	endInsertRows();
}

void FitDataCore::set_v(const int row, const int v)
{
	mData[row]->v = v;
	EmitDataChanged(row, fdcv);
}

const QString & FitDataCore::get_vs(const int row) const
{
	return reinterpret_cast<FitDataBaseData*>(mData[row])->vs;
}

void FitDataCore::set_vs(const int row, const QString& vs)
{
	reinterpret_cast<FitDataBaseData*>(mData[row])->vs = vs;
	EmitDataChanged(row, fdcvs);
}

void FitDataCore::setJ(const int row, const int J)
{
	mData[row]->J = J;
	EmitDataChanged(row, fdcJ);
}

void FitDataCore::setJs(const int row, const int Js)
{
	mData[row]->Js = Js;
	EmitDataChanged(row, fdcJs);
}

void FitDataCore::setIso(const int row, const int iso)
{
	mData[row]->isotope = iso;
	EmitDataChanged(row, fdcIso);
}

const QString & FitDataCore::getSource(const int row) const
{
	return reinterpret_cast<FitDataBaseData*>(mData[row])->source;
}

void FitDataCore::setSource(const int row, const QString& source)
{
	reinterpret_cast<FitDataBaseData*>(mData[row])->source = source;
	EmitDataChanged(row, fdcSource);
}

void FitDataCore::setSourceFile(const int row, const QString& filename)
{
	mData[row]->file = filename;
	EmitDataChanged(row, fdcFile);
}

void FitDataCore::setObsCalc(const int row, const double obsCalc)
{
	mData[row]->obsMinusCalc = obsCalc;
	EmitDataChanged(row, fdcObsCalc);
}

const QString & FitDataCore::getOtherState(const int row) const
{
	return reinterpret_cast<FitDataBaseData*>(mData[row])->secondState;
}

void FitDataCore::setUncertainty(const int row, const double uncertainty)
{
	mData[row]->uncertainty = uncertainty;
	EmitDataChanged(row, fdcUncert);
}

void FitDataCore::setEnergy(const int row, const double energy)
{
	mData[row]->energy = energy;
	EmitDataChanged(row, fdcEnergy);
}

void FitDataCore::setProgression(const int row, const int progression)
{
	mData[row]->progressionNumber = progression;
	EmitDataChanged(row, fdcProg);
}

void FitDataCore::setSecondState(const int row, const QString& state)
{
 	reinterpret_cast<FitDataBaseData*>(mData[row])->secondState = state;
	EmitDataChanged(row, fdcLineElState);
}

float FitDataCore::getDevRatio(const int row) const
{
	return reinterpret_cast<FitDataBaseData*>(mData[row])->devR;
}

void FitDataCore::setDevRatio(const int row, const float DevR)
{
	reinterpret_cast<FitDataBaseData*>(mData[row])->devR = DevR;
	EmitDataChanged(row, fdcDevR);
}

void FitDataCore::setRWError(const QString& headerText)
{
	RWError = headerText;
	emit headerDataChanged(Qt::Horizontal, fdcLineElState, fdcLineElState);
}

void FitDataCore::shrinkAllSpectRefs()
{
	int m, N = mData.size();
    QString FileName;
    for (auto it = mData.begin(); it != mData.end(); ++it)
    {
        FileName = (*it)->file;
        if ((m = FileName.lastIndexOf(QRegExp("[\\/]"))) >= 0) (*it)->file = FileName.right(FileName.length() - m - 1);
    }
    EmitDataChanged(0, mData.size(), fdcFile, fdcFile);
}

void FitDataCore::search(const int* const Rows, const int NRows, const int column, const int value, const int smeqla, QModelIndexList& Result) const
{
	int firstColumn = column, lastColumn = column, c;
	if (column == -1)
	{
		firstColumn = 0;
		lastColumn = columnCount() - 1;
	}
	switch (smeqla)
	{
		case 0:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value == element->isotope) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value == element->v) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value == element->J) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value == element->vs.toInt()) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value == element->Js) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value == element->progressionNumber) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value == static_cast<int>(element->energy)) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value == static_cast<int>(element->uncertainty)) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value == static_cast<int>(element->obsMinusCalc)) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value == static_cast<int>(element->devR)) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 1:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value > element->isotope) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value > element->v) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value > element->J) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value > element->vs.toInt()) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value > element->Js) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value > element->progressionNumber) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value > static_cast<int>(element->energy)) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value > static_cast<int>(element->uncertainty)) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value > static_cast<int>(element->obsMinusCalc)) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value > static_cast<int>(element->devR)) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 2:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value < element->isotope) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value < element->v) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value < element->J) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value < element->vs.toInt()) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value < element->Js) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value < element->progressionNumber) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value < static_cast<int>(element->energy)) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value < static_cast<int>(element->uncertainty)) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value < static_cast<int>(element->obsMinusCalc)) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value < static_cast<int>(element->devR)) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 3:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value > abs(element->isotope)) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value > abs(element->v)) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value > abs(element->J)) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value > abs(element->vs.toInt())) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value > abs(element->Js)) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value > abs(element->progressionNumber)) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value > abs(static_cast<int>(element->energy))) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value > abs(static_cast<int>(element->uncertainty))) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value > abs(static_cast<int>(element->obsMinusCalc))) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value > abs(static_cast<int>(element->devR))) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 4:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value > abs(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value > abs(mData[Rows[n]]->v)) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value > abs(mData[Rows[n]]->J)) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value > abs(element->vs.toInt())) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value > abs(mData[Rows[n]]->Js)) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value > abs(mData[Rows[n]]->progressionNumber)) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value > abs(static_cast<int>(element->energy))) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value > abs(static_cast<int>(element->uncertainty))) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value > abs(static_cast<int>(element->obsMinusCalc))) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value > abs(static_cast<int>(element->devR))) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
	}
	if (Result.empty()) QMessageBox::information(nullptr, "MolSpektAnalysis", "A number fulfilling the condition could not be found.");
}

void FitDataCore::search(const int* const Rows, const int NRows, const int column, const double value, const int smeqla, QModelIndexList& Result) const
{
	int firstColumn = column, lastColumn = column, c;
	if (column == -1)
	{
		firstColumn = 0;
		lastColumn = columnCount() - 1;
	}
	switch (smeqla)
	{
		case 0:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value == static_cast<double>(element->isotope)) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value == static_cast<double>(element->v)) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value == static_cast<double>(element->J)) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value == static_cast<double>(element->vs.toInt())) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value == static_cast<double>(element->Js)) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value == static_cast<double>(element->progressionNumber)) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value == element->energy) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value == element->uncertainty) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value == element->obsMinusCalc) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value == element->devR) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 1:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value > static_cast<double>(element->isotope)) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value > static_cast<double>(element->v)) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value > static_cast<double>(element->J)) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value > static_cast<double>(element->vs.toInt())) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value > static_cast<double>(element->Js)) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value > static_cast<double>(element->progressionNumber)) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value > element->energy) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value > element->uncertainty) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value > element->obsMinusCalc) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value > element->devR) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 2:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value < static_cast<double>(element->isotope)) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value < static_cast<double>(element->v)) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value < static_cast<double>(element->J)) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value < static_cast<double>(element->vs.toInt())) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value < static_cast<double>(element->Js)) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value < static_cast<double>(element->progressionNumber)) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value < element->energy) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value < element->uncertainty) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value < element->obsMinusCalc) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value < element->devR) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 3:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value > abs(static_cast<double>(element->isotope))) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value > abs(static_cast<double>(element->v))) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value > abs(static_cast<double>(element->J))) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value > abs(static_cast<double>(element->vs.toInt()))) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value > abs(static_cast<double>(element->Js))) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value > abs(static_cast<double>(element->progressionNumber))) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value > abs(element->energy)) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value > abs(element->uncertainty)) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value > abs(element->obsMinusCalc)) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value > abs(element->devR)) Result << createIndex(n, c);
							break;
						default:
							// not reasonable
							break;
					}
				}
			}
			break;
		case 4:
            for (int n=0; n < NRows; ++n)
			{
				FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[Rows[n]]);
				for (c = firstColumn; c <= lastColumn; ++c)
				{
					switch (c)
					{
						case fdcIso:
							if (value > abs(static_cast<double>(element->isotope))) Result << createIndex(n, c);
							break;
						case fdcv:
							if (value > abs(static_cast<double>(element->v))) Result << createIndex(n, c);
							break;
						case fdcJ:
							if (value > abs(static_cast<double>(element->J))) Result << createIndex(n, c);
							break;
						case fdcvs:
							if (value > abs(static_cast<double>(element->vs.toInt()))) Result << createIndex(n, c);
							break;
						case fdcJs:
							if (value > abs(static_cast<double>(element->Js))) Result << createIndex(n, c);
							break;
						case fdcProg:
							if (value > abs(static_cast<double>(element->progressionNumber))) Result << createIndex(n, c);
							break;
						case fdcEnergy:
							if (value > abs(element->energy)) Result << createIndex(n, c);
							break;
						case fdcUncert:
							if (value > abs(element->uncertainty)) Result << createIndex(n, c);
							break;
						case fdcObsCalc:
							if (value > abs(element->obsMinusCalc)) Result << createIndex(n, c);
							break;
						case fdcDevR:
							if (value > abs(element->devR)) Result << createIndex(n, c);
							break;
					}
				}
			}
			break;
	}
	if (Result.empty()) QMessageBox::information(nullptr, "MolSpektAnalysis", "A number fullfilling the condition could not be found.");
}

void FitDataCore::search(const int* const Rows, const int NRows, const QString& Text, QModelIndexList& Result, const int column,
					 	 const bool completeCell) const
{
	int firstColumn = column, lastColumn = column, c;
	if (column == -1)
	{
		firstColumn = 0;
		lastColumn = columnCount() - 1;
	}
	for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
	{
		QModelIndex index = createIndex(Rows[n], c);
		QString field = data(index).toString();
		if (completeCell ? field == Text : field.indexOf(Text) >= 0) Result << index; 
	}
	if (Result.empty()) QMessageBox::information(nullptr, "MolSpektAnalysis", "The text \'" + Text + "\' could not be found in the table.");
}

QString FitDataCore::cellToString(const int r, const int c) const
{
    QString R;
	FitDataBaseData* element = reinterpret_cast<FitDataBaseData*>(mData[r]);
    switch(c)
    {
        case fdcIso:
            R = QString::number(element->isotope);
            break;
        case fdcv:
            R =  QString::number(element->v);
            break;
        case fdcJ:
            R =  QString::number(element->J);
            break;
        case fdcvs:
            R =  element->vs;
            break;
        case fdcJs:
            R =  QString::number(element->Js);
            break;
        case fdcSource:
            R =  element->source;
            break;
        case fdcProg:
            R =  QString::number(element->progressionNumber);
            break;
        case fdcFile:
            R = element->file;
            break;
        case fdcEnergy:
            R =  QString::number(element->energy, 'f', getNumDecimalPlaces(element->uncertainty));
            break;
        case fdcUncert:
        {
            double u = element->uncertainty;
            R =  QString::number(u, 'f', getNumDecimalPlaces(u));
            break;
        }
        case fdcObsCalc:
            R =  QString::number(element->obsMinusCalc, 'f', getNumDecimalPlaces(element->uncertainty));
            break;
        case fdcDevR:
            R =  QString::number(element->devR, 'f', 3);
            break;
        case fdcLineElState:
            R =  element->secondState;
    }
    return R;
}
