//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//


#include "tableviewwindowcore.h"
#include "basedata.h"

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

void TableViewWindowCore::EmitDataChanged(QModelIndex& index)
{
    QVector<int> roles;
	roles.push_back(Qt::EditRole);
	emit dataChanged(index, index, roles);
}

void TableViewWindowCore::setRow(const QStringList& L, const int row)
{
	setRow(convertToBaseData(L), row);
}
