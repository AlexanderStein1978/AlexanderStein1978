//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "tableviewwindowcore.h"
#include "basedata.h"
#include "isotab.h"
#include "molecule.h"

#include <QPixmap>
#include <QPainter>


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

void TableViewWindowCore::addRow(BaseData* const data)
{
	int rc = mData.size();
	beginInsertRows(QModelIndex(), rc - 1, rc - 1);
	mData.push_back(data);
	endInsertRows();
}

void TableViewWindowCore::setRow(BaseData *const data, const int row)
{
	delete mData[row];
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

void TableViewWindowCore::startSearch(int& N, int *& Rows) const
{
    table->getSelectedRows(Rows, N);
    if (nullptr == Rows)
    {
        N = mCore->rowCount();
        Rows = new int[N];
        for (int n=0; n<N; ++n) Rows[n] = n;
    }
}

void TableViewWindowCore::finishSearch(int *const Rows, const QModelIndexList& Result) const
{
    delete[] Rows;
    if (Result.size() > 0)
    {
        QItemSelectionModel* model = new QItemSelectionModel;
        for (auto it = Result.begin(); it != Result.end(); ++it) model->select(*it, QItemSelectionModel::Select);
        table->blockSignals(true);
        table->setSelectionModel(model);
        table->scrollTo(Result[0]);
        table->blockSignals(false);
    }
}

