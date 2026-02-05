//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "tableviewwindowcore.h"
#include "tableviewwindowcoresortfunctor.h"
#include "basedata.h"
#include "isotab.h"
#include "molecule.h"
#include "heapsort.h"

#include <QPixmap>
#include <QPainter>
#include <QTextStream>


TableViewWindowCore::TableViewWindowCore(Molecule* mol, QObject *parent, QRegExp readSpecialPartRegex) : QAbstractTableModel(parent), mStartSpecialPart(readSpecialPartRegex), molecule(mol)
{
    NewPix = new QPixmap(10, 10);
    QPainter P(NewPix);
    P.setPen(QColor(255, 0, 0));
    P.setFont(QFont("Arial", 10));
    P.drawText(0, 10, "N");
}

TableViewWindowCore::~TableViewWindowCore()
{
    for (auto it = mData.begin(); it != mData.end(); ++it) delete *it;
	delete NewPix;
}

QVariant TableViewWindowCore::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();
	if (orientation == Qt::Vertical)
	{
		if (section < verticalHeaders.size()) return verticalHeaders[section];
		else return section;
	}
	else
	{
		if (section < horizontalHeders.size()) return horizontalHeders[section];
		else return section;
	}
}

void TableViewWindowCore::writeData(QTextStream& S)
{
	char Spacer = '\t';
	NSources = mData.size();
	int NC = columnCount();
	for (int r=0; r < NSources; ++r) for (int c = 0; c < NC; ++c) S << cellToString(r, c) << (c < NC - 1 ? Spacer : '\n');
}

void TableViewWindowCore::EmitDataChanged(const int row, const int column)
{
    QModelIndex index = createIndex(row, column);
    QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

void TableViewWindowCore::EmitDataChanged(const int upperRow, const int lowerRow, const int leftColumn, const int rightColumn)
{
    QModelIndex index1 = createIndex(upperRow, leftColumn);
	QModelIndex index2 = createIndex(lowerRow, rightColumn);
    QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index1, index2, roles);
}

void TableViewWindowCore::setRow(const QStringList& L, const int row)
{
	setRow(convertToBaseData(L), row);
}

int TableViewWindowCore::getMaxJ()
{
	int MJ=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if ((*it)->J > MJ) MJ = (*it)->J;
	return MJ;
}

int TableViewWindowCore::getMaxv()
{
	int Mv=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if ((*it)->v > Mv) Mv = (*it)->v;
	return Mv;
}

int TableViewWindowCore::rowCount(const QModelIndex& parent) const
{
	if (parent.isValid()) return 0;
	return mData.size();
}

void TableViewWindowCore::setRowCount(const int count)
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

int TableViewWindowCore::columnCount(const QModelIndex& parent) const
{
	if (parent.isValid()) return 0;
	return 9;
}

void TableViewWindowCore::addRow(BaseData* const data)
{
	int rc = mData.size();
	beginInsertRows(QModelIndex(), rc - 1, rc - 1);
	mData.push_back(data);
	endInsertRows();
}

void TableViewWindowCore::setRow(BaseData *const data, const int row)
{
 	if (nullptr != mData[row]) delete mData[row];
	mData[row] = data;
	QModelIndex index1 = createIndex(row, 0), index2 = createIndex(row, columnCount() - 1);
	QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index1, index2, roles);
}

BaseData * TableViewWindowCore::getRow(const int row) const
{
	return mData[row];
}

void TableViewWindowCore::insertRows(const int startRow, const std::vector<BaseData *>& newRows)
{
    auto it = mData.begin();
    for (int r=0; it != mData.end() && r < startRow; ++it) ++r;
    mData.insert(it, newRows.begin(), newRows.end());
}

void TableViewWindowCore::deleteRows(const int *indices, const int numRows)
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

void TableViewWindowCore::deleteRow(const int index)
{
	beginRemoveRows(QModelIndex(), index, index);
	int count = 0;
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); it != mData.end() && count < index; ++it) ++count;
	delete *it;
	mData.erase(it);
	endRemoveRows();
}

int TableViewWindowCore::get_v(const int row) const
{
	return mData[row]->v;
}

int TableViewWindowCore::getJ(const int row) const
{
	return mData[row]->J;
}

int TableViewWindowCore::getJs(const int row) const
{
	return mData[row]->Js;
}

int TableViewWindowCore::getIso(const int row) const
{
	return mData[row]->isotope;
}

const QString & TableViewWindowCore::getSourceFile(const int row) const
{
	return mData[row]->file;
}

double TableViewWindowCore::getObsCalc(const int row) const
{
	return mData[row]->obsMinusCalc;
}

double TableViewWindowCore::getUncertainty(const int row) const
{
	return mData[row]->uncertainty;
}

double TableViewWindowCore::getEnergy(const int row) const
{
	return mData[row]->energy;
}

int TableViewWindowCore::getProgression(const int row) const
{
	return mData[row]->progressionNumber;
}

void TableViewWindowCore::setIsoIcon(const int row, const QPixmap *const Icon)
{
	mData[row]->IsoIcon = Icon;
	EmitDataChanged(row, 0);
}

void TableViewWindowCore::setMolecule(Molecule* const mol)
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
		QModelIndex first = createIndex(0, 0), last = createIndex(mData.size() - 1, 0);
		QVector<int> roles;
		roles.push_back(Qt::DecorationRole);
		emit dataChanged(first, last, roles);
	}
}

void TableViewWindowCore::RemoveEmptyRows()
{
    int i1, i2, nr = mData.size();
    for (i1 = 0; i1 < nr; ++i1) if (mData[i1] == nullptr) break;
    beginRemoveRows(QModelIndex(), i1, i1);
    for (i2 = i1 + 1; i2 < nr; ++i2) if (mData[i2] != nullptr) mData[i1++] = mData[i2];
    mData.resize(i1);
    endRemoveRows();
}

void TableViewWindowCore::sortTab(int* S2)
{
	int i, P1=0, n, N = mData.size();
	BaseData* AIt[2];
	for (i=0; i<N; ++i) if (S2[i] != i)
	{
		AIt[P1] = mData[i];
		while (S2[i] != i)
		{
			AIt[1-P1] = mData[S2[i]];
			mData[S2[i]] = AIt[P1];
			P1 = 1 - P1;
			n = S2[S2[i]];
			S2[S2[i]] = S2[i];
			S2[i] = n;
		}
		mData[i] = AIt[P1];
	}
	delete[] S2;
}

void TableViewWindowCore::shrinkAllSpectRefs()
{
	int m, N = mData.size();
    QString FileName;
    for (auto it = mData.begin(); it != mData.end(); ++it)
    {
        FileName = (*it)->file;
        if ((m = FileName.lastIndexOf(QRegExp("[\\/]"))) >= 0) (*it)->file = FileName.right(FileName.length() - m - 1);
    }
}

void TableViewWindowCore::setData(const int row, BaseData * const data)
{
	mData[row] = data;
	EmitDataChanged(row, row, 0, columnCount());
}

void TableViewWindowCore::setHorizontalHeader(const QStringList &Labels)
{
	int oldsize = horizontalHeders.size();
	horizontalHeders = Labels;
    emit headerDataChanged(Qt::Horizontal, 0, (Labels.size() > oldsize ? Labels.size() : oldsize));
}

void TableViewWindowCore::setVerticalHeader(const QStringList &Labels)
{
	int oldsize = verticalHeaders.size();
	verticalHeaders = Labels;
    emit headerDataChanged(Qt::Vertical, 0, (Labels.size() > oldsize ? Labels.size() : oldsize));
}

QString TableViewWindowCore::cellToString(const int r, const int c) const
{
	QModelIndex index = createIndex(r, c);
	return data(index, Qt::DisplayRole).toString();
}

void TableViewWindowCore::setUncertainty(const int row, const double uncertainty)
{
	mData[row]->uncertainty = uncertainty;
	EmitDataChanged(row, row, 0, columnCount());
}

void TableViewWindowCore::setProgression(const int row, const int progression)
{
	mData[row]->progressionNumber = progression;
	EmitDataChanged(row, row, 0, columnCount());
}

void TableViewWindowCore::addRow(const QStringList& L)
{
	int r = mData.size();
	beginInsertRows(QModelIndex(), r, r);
	mData.push_back(convertToBaseData(L));
	endInsertRows();
}
