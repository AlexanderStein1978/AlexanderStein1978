//
// C++ Implementation: MTable
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2017
//
// Copyright: See README file that comes with this source code
//
//


#include "mtable.h"


MTable::MTable(QWidget* parent): QTableView(parent)
{

}

void MTable::getSelectedRows(int* Rows, int NR)
{
	QModelIndexList L = selectedIndexes();
	int n, N = L.count(), M=0, r;
	for (n=0; n<N; n++) if ((r = L[n].row()) > M) M=r;
	bool *B = new bool[M];
	for (n=0; n<M; n++) B[n] = false;
	for (n=0; n<N; n++) B[L[n].row()] = true;
	for (NR = n = 0; n<M; n++) if (B[n]) NR++;
	Rows = new int[NR];
	for (n = NR = 0; n<M; n++) if (B[n]) Rows[NR++] = n;
	delete[] B;
}
