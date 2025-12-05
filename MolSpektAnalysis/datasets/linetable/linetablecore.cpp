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

#include <QTextStream>
#include <QPixmap>
#include <QPainter>
#include <QMessageBox>


LineTableCore::LineTableCore(Molecule* mol, QObject *parent) : QAbstractTableModel(parent), molecule(mol)
{
	NewPix = new QPixmap(10, 10);
    QPainter P(NewPix);
    P.setPen(QColor(255, 0, 0));
    P.setFont(QFont("Arial", 10));
    P.drawText(0, 10, "N");
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
	for (auto it = mData.begin(); it != mData.end(); ++it) delete *it;
	delete NewPix;
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

QString LineTableCore::readData(QTextStream& S)
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

void LineTableCore::setRow(const QStringList& L, const int row)
{
	setRow(convertToBaseData(L), row);
}

void LineTableCore::writeData(QTextStream& S)
{
	QString Spacer = "\t";
	for (auto it = mData.begin(); it != mData.end(); ++it)
	{
		int numDigits = getNumDecimalPlaces((*it)->uncert);
		S << static_cast<int>((*it)->isotope) << Spacer << (*it)->v << Spacer << (*it)->J << Spacer << (*it)->vs.c_str() << Spacer << (*it)->Js << Spacer << (*it)->source.c_str() << Spacer << (*it)->prog
		  << Spacer << (*it)->file.c_str() << Spacer << QString::number((*it)->energy, 'f', numDigits) << Spacer << QString::number((*it)->uncert, 'f', numDigits) << Spacer
		  << QString::number((*it)->obs_calc, 'f', numDigits) << Spacer << QString::number((*it)->devR, 'f', 3) << Spacer << (*it)->secondState.c_str() << '\n';
	}
	NSources = mData.size();
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

int LineTableCore::getMaxJ()
{
	int MJ=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if ((*it)->J > MJ) MJ = (*it)->J;
	return MJ;
}

int LineTableCore::getMaxv()
{
	int Mv=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if ((*it)->v > Mv) Mv = (*it)->v;
	return Mv;
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

int LineTableCore::rowCount(const QModelIndex& parent) const
{
	if (parent.isValid()) return 0;
	return mData.size();
}

void LineTableCore::setRowCount(const int count)
{
	int currentSize = mData.size();
	if (count > currentSize)
	{
		beginInsertRows(QModelIndex(), currentSize, count - 2);
		mData.resize(count);
		endInsertRows();
	}
	else if (count < currentSize)
	{
		beginRemoveRows(QModelIndex(), count, currentSize - 2);
		mData.resize(count);
		endRemoveRows();
	}
}

int LineTableCore::addMarkedLevel(TermEnergy& TE, Spektrum* Source)
{
	int R = mData.size();
	beginInsertRows(QModelIndex(), R-1, R-1);
	BaseData* element = new BaseData;
	element->isotope = static_cast<char>(TE.Iso);
	element->v = TE.v;
	element->J = TE.J;
	element->vs = (TE.FC == -1 ? "TE" : std::to_string(-1 - TE.FC));
	element->source = Source->getName().toStdString();
	element->file = Source->getFileName().toStdString();
	element->energy = TE.E;
	element->uncert = TE.err;
	element->obs_calc = TE.dev;
	element->devR = TE.DevR;
	element->secondState = TE.State->getName().toStdString();
	mData.push_back(element);
	endInsertRows();
    return R;
}

int LineTableCore::addRow(const int cr)
{
	beginInsertRows(QModelIndex(), cr - 1, cr - 1);
	int r=0, nr = mData.size();
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); r < cr; ++r) ++it;
	mData.insert(it, new BaseData);
	endInsertRows();
	return nr;
}

void LineTableCore::addRow(BaseData* const data)
{
	int rc = mData.size();
	beginInsertRows(QModelIndex(), rc - 1, rc - 1);
	mData.push_back(data);
	endInsertRows();
}

void LineTableCore::addRow(const QStringList& L)
{
	addRow(convertToBaseData(L));
}

void LineTableCore::setRow(BaseData *const data, const int row)
{
	delete mData[row];
	mData[row] = data;
	QModelIndex index1 = createIndex(row, fdcIso), index2 = createIndex(row, fdcLineElState);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index1, index2, roles);
}

BaseData * LineTableCore::getRow(const int row) const
{
	return mData[row];
}

void LineTableCore::addData(const int i_numLines, int *const i_Lines, const LineTableCore& data)
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

void LineTableCore::deleteRows(const int *indices, const int numRows)
{
	int min = indices[0], max = min;
	for (int i=1; i < numRows; ++i)
	{
		if (indices[i] < min) min = indices[i];
		else if (indices[i] > max) max = indices[i];
	}
	beginRemoveRows(QModelIndex(), min, max);
	int size = mData.size(), i, j;
	for (i=0; i < numRows; ++i)
	{
		delete mData[indices[i]];
		mData[indices[i]] = nullptr;
	}
	for (i=0; i < size && mData[i] == nullptr; ++i) ;
	for (j=i+1; j < size; j++) if (mData[j] != nullptr) mData[i++] = mData[j];
	mData.resize(i);
	endRemoveRows();
}

void LineTableCore::deleteRow(const int index)
{
	beginRemoveRows(QModelIndex(), index, index);
	int count = 0;
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); it != mData.end() && count < index; ++it) ++count;
	delete *it;
	mData.erase(it);
	endRemoveRows();
}

int LineTableCore::get_vs(const int row) const
{
	return mData[row]->vs;
}

void LineTableCore::set_vs(const int row, const int vs)
{
	mData[row]->vs = vs;
	QModelIndex index = createIndex(row, Cvs);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const int LineTableCore::get_vss(const int row) const
{
	return mData[row]->vss;
}

void LineTableCore::set_vss(const int row, const int vss)
{
	mData[row]->vss = vss;
	QModelIndex index = createIndex(row, Cvss);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int LineTableCore::getJs(const int row) const
{
	return mData[row]->Js;
}

void LineTableCore::setJs(const int row, const int Js)
{
	mData[row]->Js = Js;
	QModelIndex index = createIndex(row, CJs);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int LineTableCore::getJss(const int row) const
{
	return mData[row]->Jss;
}

void LineTableCore::setJss(const int row, const int Jss)
{
	mData[row]->Jss = Jss;
	QModelIndex index = createIndex(row, CJss);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int LineTableCore::getF(const int row) const
{
	return mData[row]->F;
}

void LineTableCore::setF(const int row, const int F)
{
	mData[row]->F = F;
	QModelIndex index = createIndex(row, CF);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int LineTableCore::getIso(const int row) const
{
	return mData[row]->isotope;
}

void LineTableCore::setIso(const int row, const int iso)
{
	mData[row]->isotope = iso;
	QModelIndex index = createIndex(row, CIso);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & LineTableCore::getSource(const int row) const
{
	return mData[row]->source;
}

void LineTableCore::setSource(const int row, const std::string& source)
{
	mData[row]->source = source;
	QModelIndex index = createIndex(row, fdcSource);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & LineTableCore::getFile(const int row) const
{
	return mData[row]->file;
}

void LineTableCore::setFile(const int row, const std::string& filename)
{
	mData[row]->file = filename;
	QModelIndex index = createIndex(row, CFile);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

double LineTableCore::getObsCalc(const int row) const
{
	return mData[row]->obsMinusCalc;
}

void LineTableCore::setObsCalc(const int row, const double obsCalc)
{
	mData[row]->obsMinusCalc = obsCalc;
	QModelIndex index = createIndex(row, CDev);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & LineTableCore::getOtherState(const int row) const
{
	return mData[row]->secondState;
}

double LineTableCore::getUncertainty(const int row) const
{
	return mData[row]->uncert;
}

void LineTableCore::setUncertainty(const int row, const double uncertainty)
{
	mData[row]->uncert = uncertainty;
	QModelIndex index = createIndex(row, fdcUncert);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

double LineTableCore::getWaveNumber(const int row) const
{
	return mData[row]->waveNumber;
}

void LineTableCore::setWaveNumber(const int row, const double waveNumber)
{
	mData[row]->waveNumber = waveNumber;
	QModelIndex index = createIndex(row, CWN);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int LineTableCore::getProgression(const int row) const
{
	return mData[row]->prog;
}

void LineTableCore::setProgression(const int row, const int progression)
{
	mData[row]->prog = progression;
	QModelIndex index = createIndex(row, fdcProg);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

void LineTableCore::setSecondState(const int row, const std::string& state)
{
	mData[row]->secondState = state;
	QModelIndex index = createIndex(row, fdcLineElState);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

float LineTableCore::getDevRatio(const int row) const
{
	return mData[row]->devR;
}

void LineTableCore::setDevRatio(const int row, const float DevR)
{
	mData[row]->devR = DevR;
	QModelIndex index = createIndex(row, fdcDevR);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

void LineTableCore::setIsoIcon(const int row, const QPixmap *const Icon)
{
	mData[row]->IsoIcon = Icon;
	QModelIndex index = createIndex(row, fdcIso);
	QVector<int> roles;
	roles.push_back(Qt::DecorationRole);
	emit dataChanged(index, index, roles);
}

void LineTableCore::setRWError(const QString& headerText)
{
	RWError = headerText;
	emit headerDataChanged(Qt::Horizontal, fdcLineElState, fdcLineElState);
}

void LineTableCore::setMolecule(Molecule* const mol)
{
	if (mol != molecule)
	{
		molecule = mol;
		if (nullptr != mol)
		{
			IsoTab* Iso = mol->getIso();
			for (auto it = mData.begin(); it != mData.end(); ++it) (*it)->IsoIcon = &Iso->IsoImage[(*it)->isotope];
		}
		else for (auto it = mData.begin(); it != mData.end(); ++it) (*it)->IsoIcon = nullptr;
		QModelIndex first = createIndex(0, fdcIso), last = createIndex(mData.size() - 1, fdcIso);
		QVector<int> roles;
		roles.push_back(Qt::DecorationRole);
		emit dataChanged(first, last, roles);
	}
}

void LineTableCore::shrinkAllSpectRefs()
{
	int m, N = mData.size();
    QString FileName;
    for (auto it = mData.begin(); it != mData.end(); ++it)
    {
        FileName = (*it)->file.c_str();
        if ((m = FileName.lastIndexOf(QRegExp("[\\/]"))) >= 0) (*it)->file = FileName.right(FileName.length() - m - 1).toStdString();
    }
}

void LineTableCore::search(const int* const Rows, const int NRows, const int column, const int value, const int smeqla, QModelIndexList& Result) const
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
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value == static_cast<int>(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value == static_cast<int>(mData[Rows[n]]->v)) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value == static_cast<int>(mData[Rows[n]]->J)) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value == stdStringToInt(mData[Rows[n]]->vs), -1) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value == static_cast<int>(mData[Rows[n]]->Js)) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value == mData[Rows[n]]->prog) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value == static_cast<int>(mData[Rows[n]]->energy)) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value == static_cast<int>(mData[Rows[n]]->uncert)) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value == static_cast<int>(mData[Rows[n]]->obs_calc)) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value == static_cast<int>(mData[Rows[n]]->devR)) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 1:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value > static_cast<int>(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value > static_cast<int>(mData[Rows[n]]->v)) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value > static_cast<int>(mData[Rows[n]]->J)) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value > stdStringToInt(mData[Rows[n]]->vs), -1) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value > static_cast<int>(mData[Rows[n]]->Js)) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value > mData[Rows[n]]->prog) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value > static_cast<int>(mData[Rows[n]]->energy)) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value > static_cast<int>(mData[Rows[n]]->uncert)) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value > static_cast<int>(mData[Rows[n]]->obs_calc)) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value > static_cast<int>(mData[Rows[n]]->devR)) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 2:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value < static_cast<int>(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value < static_cast<int>(mData[Rows[n]]->v)) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value < static_cast<int>(mData[Rows[n]]->J)) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value < stdStringToInt(mData[Rows[n]]->vs, -1)) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value < static_cast<int>(mData[Rows[n]]->Js)) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value < mData[Rows[n]]->prog) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value < static_cast<int>(mData[Rows[n]]->energy)) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value < static_cast<int>(mData[Rows[n]]->uncert)) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value < static_cast<int>(mData[Rows[n]]->obs_calc)) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value < static_cast<int>(mData[Rows[n]]->devR)) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 3:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value > abs(static_cast<int>(mData[Rows[n]]->isotope))) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value > abs(static_cast<int>(mData[Rows[n]]->v))) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value > abs(static_cast<int>(mData[Rows[n]]->J))) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value > abs(stdStringToInt(mData[Rows[n]]->vs, -1))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value > abs(static_cast<int>(mData[Rows[n]]->Js))) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value > abs(mData[Rows[n]]->prog)) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value > abs(static_cast<int>(mData[Rows[n]]->energy))) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value > abs(static_cast<int>(mData[Rows[n]]->uncert))) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value > abs(static_cast<int>(mData[Rows[n]]->obs_calc))) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value > abs(static_cast<int>(mData[Rows[n]]->devR))) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 4:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value > abs(static_cast<int>(mData[Rows[n]]->isotope))) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value > abs(static_cast<int>(mData[Rows[n]]->v))) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value > abs(static_cast<int>(mData[Rows[n]]->J))) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value > abs(stdStringToInt(mData[Rows[n]]->vs, -1))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value > abs(static_cast<int>(mData[Rows[n]]->Js))) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value > abs(mData[Rows[n]]->prog)) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value > abs(static_cast<int>(mData[Rows[n]]->energy))) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value > abs(static_cast<int>(mData[Rows[n]]->uncert))) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value > abs(static_cast<int>(mData[Rows[n]]->obs_calc))) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value > abs(static_cast<int>(mData[Rows[n]]->devR))) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
	}
	if (Result.empty()) QMessageBox::information(nullptr, "MolSpektAnalysis", "A number fulfilling the condition could not be found.");
}

void LineTableCore::search(const int* const Rows, const int NRows, const int column, const double value, const int smeqla, QModelIndexList& Result) const
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
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value == static_cast<double>(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value == static_cast<double>(mData[Rows[n]]->v)) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value == static_cast<double>(mData[Rows[n]]->J)) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value == static_cast<double>(stdStringToInt(mData[Rows[n]]->vs, -1))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value == static_cast<double>(mData[Rows[n]]->Js)) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value == mData[Rows[n]]->prog) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value == mData[Rows[n]]->energy) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value == mData[Rows[n]]->uncert) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value == mData[Rows[n]]->obs_calc) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value == static_cast<double>(mData[Rows[n]]->devR)) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 1:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value > static_cast<double>(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value > static_cast<double>(mData[Rows[n]]->v)) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value > static_cast<double>(mData[Rows[n]]->J)) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value > static_cast<double>(stdStringToInt(mData[Rows[n]]->vs, -1))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value > static_cast<double>(mData[Rows[n]]->Js)) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value > mData[Rows[n]]->prog) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value > mData[Rows[n]]->energy) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value > mData[Rows[n]]->uncert) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value > mData[Rows[n]]->obs_calc) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value > static_cast<double>(mData[Rows[n]]->devR)) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 2:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value < static_cast<double>(mData[Rows[n]]->isotope)) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value < static_cast<double>(mData[Rows[n]]->v)) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value < static_cast<double>(mData[Rows[n]]->J)) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value < static_cast<double>(stdStringToInt(mData[Rows[n]]->vs, -1))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value < static_cast<double>(mData[Rows[n]]->Js)) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value < mData[Rows[n]]->prog) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value < mData[Rows[n]]->energy) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value < mData[Rows[n]]->uncert) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value < mData[Rows[n]]->obs_calc) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value < static_cast<double>(mData[Rows[n]]->devR)) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 3:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value > abs(static_cast<double>(mData[Rows[n]]->isotope))) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value > abs(static_cast<double>(mData[Rows[n]]->v))) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value > abs(static_cast<double>(mData[Rows[n]]->J))) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value > abs(static_cast<double>(stdStringToInt(mData[Rows[n]]->vs, -1)))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value > abs(static_cast<double>(mData[Rows[n]]->Js))) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value > abs(static_cast<double>(mData[Rows[n]]->prog))) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value > abs(mData[Rows[n]]->energy)) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value > abs(mData[Rows[n]]->uncert)) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value > abs(mData[Rows[n]]->obs_calc)) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value > abs(static_cast<double>(mData[Rows[n]]->devR))) Result << createIndex(n, c);
						break;
					default:
						// not reasonable
						break;
				}
			}
			break;
		case 4:
            for (int n=0; n < NRows; ++n) for (c = firstColumn; c <= lastColumn; ++c)
			{
				switch (c)
				{
					case fdcIso:
						if (value > abs(static_cast<double>(mData[Rows[n]]->isotope))) Result << createIndex(n, c);
						break;
					case fdcv:
						if (value > abs(static_cast<double>(mData[Rows[n]]->v))) Result << createIndex(n, c);
						break;
					case fdcJ:
						if (value > abs(static_cast<double>(mData[Rows[n]]->J))) Result << createIndex(n, c);
						break;
					case fdcvs:
						if (value > abs(static_cast<double>(stdStringToInt(mData[Rows[n]]->vs, -1)))) Result << createIndex(n, c);
						break;
					case fdcJs:
						if (value > abs(static_cast<double>(mData[Rows[n]]->Js))) Result << createIndex(n, c);
						break;
					case fdcProg:
						if (value > abs(static_cast<double>(mData[Rows[n]]->prog))) Result << createIndex(n, c);
						break;
					case fdcEnergy:
						if (value > abs(mData[Rows[n]]->energy)) Result << createIndex(n, c); 
						break;
					case fdcUncert:
						if (value > abs(mData[Rows[n]]->uncert)) Result << createIndex(n, c);
						break;
					case fdcObsCalc:
						if (value > abs(mData[Rows[n]]->obs_calc)) Result << createIndex(n, c);
						break;
					case fdcDevR:
						if (value > abs(static_cast<double>(mData[Rows[n]]->devR))) Result << createIndex(n, c);
						break;
				}
			}
			break;
	}
	if (Result.empty()) QMessageBox::information(nullptr, "MolSpektAnalysis", "A number fullfilling the condition could not be found.");
}

void LineTableCore::search(const int* const Rows, const int NRows, const QString& Text, QModelIndexList& Result, const int column,
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
