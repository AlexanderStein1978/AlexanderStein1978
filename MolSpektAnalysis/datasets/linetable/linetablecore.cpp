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

