//
// C++ Implementation: Table
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2020
//
// Copyright: See README file that comes with this source code
//
//


#include "Table.h"

#include <QKeyEvent>


Table::Table(int NR, int NC, QWidget *parent) : QTableWidget(NR, NC, parent)
{
}

void Table::keyPressEvent(QKeyEvent *E)
{
	//printf("Table::keyPressEvent\n");
	if (E->key() == Qt::Key_Delete && E->modifiers() == Qt::ControlModifier)
	{
		int i, j, r1, r2, N = rowCount(), C = columnCount();
		bool c[N];
		QList<QTableWidgetSelectionRange> R = selectedRanges();
		for (i=0; i<N; i++) c[i] = true;
		for (i=0; i < R.count(); i++) for (j = R[i].topRow(); j <= R[i].bottomRow(); j++) c[j] = false;
		emit RowDeleted(c);
		for (r1=0; (r1 < N ? c[r1] : false); r1++) ;
		for (r2 = r1 + 1; r2 < N; r2++) if (c[r2])
		{
			for (i=0; i<C; i++) setItem(r1, i, takeItem(r2, i));
			r1++;
		}
		setRowCount(r1);
		E->accept();
	}
	else QTableWidget::keyPressEvent(E);
}
