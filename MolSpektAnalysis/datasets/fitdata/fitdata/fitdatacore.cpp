//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "fitdatacore.h"
#include "utils.h"
#include "isotab.h"
#include "basedata.h"
#include "termenergy.h"
#include "Spektrum.h"
#include "elstate.h"

#include <QTextStream>
#include <QPixmap>
#include <QPainter>


FitDataCore::FitDataCore(QObject *parent) : QAbstractTableModel(parent)
{
	NewPix = new QPixmap(10, 10);
    QPainter P(NewPix);
    P.setPen(QColor(255, 0, 0));
    P.setFont(QFont("Arial", 10));
    P.drawText(0, 10, "N");
}

FitDataCore::~FitDataCore()
{
	for (auto it = mData.begin(); it != mData.end(); ++it) delete *it;
	delete NewPix;
}

QString FitDataCore::readData(QTextStream& S)
{
	QString Buffer;
	QString Spacer = "\t";
	while(!S.atEnd())
	{
		if (Buffer.indexOf(mStartSpecialPart) >= 0) return Buffer;
		const QStringList L = Buffer.split(Spacer);
		BaseData *data = new BaseData;
		int lc = L.count();
		for (int n=0; n < lc && n <= fdcLineElState; n++)
		{
			switch(n)
			{
				case fdcIso:
					data->isotope = static_cast<char>(L[n].toInt());
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
		mData.push_back(data);
	}
	return "";
}

void FitDataCore::writeData(QTextStream& S)
{
	QString Spacer = "\t";
	for (auto it = mData.begin(); it != mData.end(); ++it)
	{
		int numDigits = 2 - static_cast<int>(log10((*it)->uncert));
		S << static_cast<int>((*it)->isotope) << Spacer << (*it)->v << Spacer << (*it)->J << Spacer << (*it)->vs.c_str() << Spacer << (*it)->Js << Spacer << (*it)->source.c_str() << Spacer << (*it)->prog
		  << Spacer << (*it)->file.c_str() << Spacer << QString::number((*it)->energy, 'f', numDigits) << Spacer << QString::number((*it)->uncert, 'f', numDigits) << Spacer
		  << QString::number((*it)->obs_calc, 'f', numDigits) << Spacer << QString::number((*it)->devR, 'f', 3) << Spacer << (*it)->secondState.c_str() << '\n';
	}
	NSources = mData.size();
}

int FitDataCore::columnCount(const QModelIndex &parent) const
{
	if (parent.isValid()) return 0;
	return 13;
}

QVariant FitDataCore::data(const QModelIndex &index, int role) const
{
	int row = index.row(), column = index.column(), numDigits;
	if (role == Qt::DecorationRole && row >= NSources) return *NewPix;
	if (role != Qt::DisplayRole) return QVariant();
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

bool FitDataCore::setData(const QModelIndex& index, const QVariant& value, int role)
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

int FitDataCore::getMaxJ()
{
	int MJ=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if ((*it)->J > MJ) MJ = (*it)->J;
	return MJ;
}

int FitDataCore::getMaxv()
{
	int Mv=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if ((*it)->v > Mv) Mv = (*it)->v;
	return Mv;
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

int FitDataCore::rowCount(const QModelIndex& parent) const
{
	if (parent.isValid()) return 0;
	return mData.size();
}

void FitDataCore::setRowCount(const int count)
{
	int currentSize = mData.size();
	if (count > currentSize)
	{
		beginInsertRows(QModelIndex(), currentSize, count - 1);
		mData.resize(count);
		endInsertRows();
	}
	else if (count < currentSize)
	{
		beginRemoveRows(QModelIndex(), count, currentSize - 1);
		mData.resize(count);
		endRemoveRows();
	}
}

int FitDataCore::addMarkedLevel(TermEnergy& TE, Spektrum* Source)
{
	int R = mData.size();
	beginInsertRows(QModelIndex(), R, R+1);
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

int FitDataCore::addRow(const int cr)
{
	beginInsertRows(QModelIndex(), cr, cr);
	int r=0, nr = mData.size();
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); r < cr; ++r) ++it;
	mData.insert(it, new BaseData);
	endInsertRows();
	return nr;
}

void FitDataCore::addRow(BaseData* const data)
{
	int rc = mData.size();
	beginInsertRows(QModelIndex(), rc, rc);
	mData.push_back(data);
	endInsertRows();
}

void FitDataCore::setRow(BaseData *const data, const int row)
{
	delete mData[row];
	mData[row] = data;
	QModelIndex index1 = createIndex(row, fdcIso), index2 = createIndex(row, fdcLineElState);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index1, index2, roles);
}

void FitDataCore::addData(const int i_numLines, int *const i_Lines, const FitDataCore& data)
{
	int dataSize = mData.size();
	beginInsertRows(QModelIndex(), dataSize, dataSize + data.mData.size());
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

void FitDataCore::deleteRows(const int *indices, const int numRows)
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

void FitDataCore::deleteRow(const int index)
{
	beginRemoveRows(QModelIndex(), index, index);
	int count = 0;
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); it != mData.end() && count < index; ++it) ++count;
	delete *it;
	mData.erase(it);
	endRemoveRows();
}

int FitDataCore::get_v(const int row) const
{
	return mData[row]->v;
}

void FitDataCore::set_v(const int row, const int v)
{
	mData[row]->v = v;
	QModelIndex index = createIndex(row, fdcv);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & FitDataCore::get_vs(const int row) const
{
	return mData[row]->vs;
}

void FitDataCore::set_vs(const int row, const std::string& vs)
{
	mData[row]->vs = vs;
	QModelIndex index = createIndex(row, fdcvs);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int FitDataCore::getJ(const int row) const
{
	return mData[row]->J;
}

void FitDataCore::setJ(const int row, const int J)
{
	mData[row]->J = J;
	QModelIndex index = createIndex(row, fdcJ);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int FitDataCore::getJs(const int row) const
{
	return mData[row]->Js;
}

void FitDataCore::setJs(const int row, const int Js)
{
	mData[row]->Js = Js;
	QModelIndex index = createIndex(row, fdcJs);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int FitDataCore::getIso(const int row) const
{
	return mData[row]->isotope;
}

void FitDataCore::setIso(const int row, const int iso)
{
	mData[row]->isotope = iso;
	QModelIndex index = createIndex(row, fdcIso);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & FitDataCore::getSource(const int row) const
{
	return mData[row]->source;
}

void FitDataCore::setSource(const int row, const std::string& source)
{
	mData[row]->source = source;
	QModelIndex index = createIndex(row, fdcSource);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & FitDataCore::getSourceFile(const int row) const
{
	return mData[row]->file;
}

void FitDataCore::setSourceFile(const int row, const std::string& filename)
{
	mData[row]->file = filename;
	QModelIndex index = createIndex(row, fdcFile);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

double FitDataCore::getObsCalc(const int row) const
{
	return mData[row]->obs_calc;
}

void FitDataCore::setObsCalc(const int row, const double obsCalc)
{
	mData[row]->obs_calc = obsCalc;
	QModelIndex index = createIndex(row, fdcObsCalc);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

const std::string & FitDataCore::getOtherState(const int row) const
{
	return mData[row]->secondState;
}

double FitDataCore::getUncertainty(const int row) const
{
	return mData[row]->uncert;
}

void FitDataCore::setUncertainty(const int row, const double uncertainty)
{
	mData[row]->uncert = uncertainty;
	QModelIndex index = createIndex(row, fdcUncert);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

double FitDataCore::getEnergy(const int row) const
{
	return mData[row]->energy;
}

void FitDataCore::setEnergy(const int row, const double energy)
{
	mData[row]->energy = energy;
	QModelIndex index = createIndex(row, fdcEnergy);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

int FitDataCore::getProgression(const int row) const
{
	return mData[row]->prog;
}

void FitDataCore::setProgression(const int row, const int progression)
{
	mData[row]->prog = progression;
	QModelIndex index = createIndex(row, fdcProg);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

void FitDataCore::setSecondState(const int row, const std::string& state)
{
	mData[row]->secondState = state;
	QModelIndex index = createIndex(row, fdcLineElState);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

float FitDataCore::getDevRatio(const int row) const
{
	return mData[row]->devR;
}

void FitDataCore::setDevRatio(const int row, const float DevR)
{
	mData[row]->devR = DevR;
	QModelIndex index = createIndex(row, fdcDevR);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

void FitDataCore::setRWError(const QString& headerText)
{
	RWError = headerText;
	emit headerDataChanged(Qt::Horizontal, fdcLineElState, fdcLineElState);
}
