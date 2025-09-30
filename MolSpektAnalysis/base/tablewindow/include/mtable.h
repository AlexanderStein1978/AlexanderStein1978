//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MTABLE_H
#define MTABLE_H


#include <QTableWidget>


class MTable : public QTableView
{
	Q_OBJECT
	
public:
	MTable(QWidget *parent = 0);
	void getSelectedRows(int *Rows, int N);
	
protected slots:
	
	inline void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected)
	{
		QTableView::selectionChanged(selected, deselected);
		emit SelChanged();
	}
	
signals:
	void SelChanged();
};

#endif
