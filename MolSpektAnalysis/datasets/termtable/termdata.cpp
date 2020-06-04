//
// C++ Implementation: termtable
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "termdata.h"


TermData::TermData(QObject *parent) : QAbstractTableModel(parent)
{
	Data = 0;
	Iso = 0;
	numComp = numIso = numRows = numv = numJ = numStates = 0;
	Z = 0;
	vStart = 0;
	CompZ = 0;
	MixCoeff = 0;
}

TermData::~TermData()
{
	if (Data != 0) 
	{
		Destroy(Data, numComp, numIso, numv);
		Destroy(vStart, numIso, numComp);
	}
	if (Iso != 0) delete Iso;
	if (Z != 0) delete[] Z;
	if (CompZ != 0) delete[] CompZ;
}

int TermData::columnCount(const QModelIndex &parent) const
{
	if (parent.isValid()) return 0;
	return 6 + numStates;
}

QVariant TermData::data(const QModelIndex &index, int role) const
{
	//printf("TermData::data\n");
	if (!index.isValid() || Data == 0) return QVariant();
	if (role == Qt::TextAlignmentRole) return Qt::AlignRight; 
	if (role != Qt::DisplayRole) return QVariant();
	int r = index.row(), c = index.column(), I, C, v, J, n;
	//printf("r=%d, c=%d, role=%d\n", r, c, role);
	for (I=1; (I < numIso ? r >= vStart[I][0][0] : false); I++) ;
	I--;
	if (c==0) 
	{
		if (Iso != 0) return (QString::number(Iso->mNumIso1[Z[I]]) + *Iso->chSymb1).toLatin1().data();
		else return I;
	}
	if (c==1) 
	{
		if (Iso != 0) return (QString::number(Iso->mNumIso2[Z[I]]) + *Iso->chSymb2).toLatin1().data();
		else return I;
	}
	for (C=1; (C < numComp ? r >= vStart[I][C][0] : false); C++) ;
	C--;
	if (c==2) return CompZ[C];
	for (v=1; (v < numv ? r >= vStart[I][C][v] : false); v++) ;
	v--;
	if (c==3) return v;
	for (n=0, J = r - vStart[I][C][v]; n<=J; n++) if (Data[C][I][v][n] == 0.0) J++;
	//printf("I=%d, numIso=%d, C=%d, numComp=%d, v=%d, numv=%d\n", I, numIso, C, numComp, v, numv);
	//printf("vStart[%d][%d][%d]=%d\n", I, C, v, vStart[I][C][v]);
	//printf("r=%d, numRows=%d, J=%d, numJ=%d\n", r, numRows, J, numJ);
	if (c==4) return J;
	if (c==5) return QString::number(Data[C][I][v][J], 'f', 4);
	return QString::number(MixCoeff[C][I][v][J][c-6], 'f', 5);
}

void TermData::getJE(int* R, int N, int* J, double* E)
{
	int I, C, v, n, lI, lC;
	for (n=I=C=v=lI=lC=0; n<N; n++)
	{
		while (I < numIso ? R[n] >= vStart[I][0][0] : false) I++;
		I--;
		if (I != lI) C=v=0;
		while (C < numComp ? R[n] >= vStart[I][C][0] : false) C++;
		C--;
		if (C != lC) v=0;
		while (v < numv ? R[n] >= vStart[I][C][v] : false) v++;
		v--;
		J[n] = R[n] - vStart[I][C][v];
		E[n] = Data[C][I][v][J[n]];
		lI = I;
		lC = C;
	}
}

void TermData::getRows(int C, int I, int v, int* J, int N, int* R)
{
	int n;
	for (n=0; n<N; n++) R[n] = vStart[I][C][v] + J[n];
}

int* TermData::getCompZ()
{
	return CompZ;
}

double ****TermData::getData()
{
	return Data;
}

void TermData::GetIsoZ(int IsoI, int &NIso1, int &NIso2)
{
	if (Iso == 0 || Z == 0) NIso1 = NIso2 = 0;
	else
	{
		NIso1 = Iso->mNumIso1[Z[IsoI]];
		NIso2 = Iso->mNumIso2[Z[IsoI]];
	}
}

int TermData::getMaxJ()
{
	return numJ - 1;
}

int TermData::getMaxv()
{
	return numv - 1;
}

int TermData::getNumComp()
{
	return numComp;
}

int TermData::getNumIso()
{
	return numIso;
}

QVariant TermData::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole) return QVariant();
	if (orientation == Qt::Vertical) return section;
	switch (section)
	{
		case 0:
			return "Iso 1";
			break;
		case 1:
			return "Iso 2";
			break;
		case 2:
			return "component";
			break;
		case 3:
			return "v";
			break;
		case 4:
			return "J";
			break;
		case 5:
			return "term energy";
			break;
		default:
			return "Mix C S" + QString::number(section - 6);
			break;
	}
	return QVariant();
}

int TermData::rowCount(const QModelIndex& parent) const
{
	if (parent.isValid()) return 0;
	return numRows;
}

void TermData::setData(double ****nData, int nC, int nI, int mv, int mJ, int *CZ, int nMixC, double *****newMixC)
{
	int r, I, C, v, J;
	if (Data != 0)
	{
		beginRemoveRows(QModelIndex(), 0, numRows - 1);
		Destroy(Data, numComp, numIso, numv);
		Destroy(vStart, numIso, numComp);
		if (MixCoeff != 0)
		{
			beginRemoveColumns(QModelIndex(), 6, 5 + numStates);
			Destroy(MixCoeff, numComp, numIso, numv, numJ);
			MixCoeff = 0;
			numStates = 0;
			endRemoveColumns();
		}
		Data = 0;
		vStart = 0;
		numRows = numIso = numComp = numv = numJ = 0;
		delete[] CompZ;
		endRemoveRows();
	}
	if (CZ != 0) CompZ = CZ;
	else
	{
		CompZ = new int[nC];
		for (r=0; r < nC; r++) CompZ[r] = r;
	}
	vStart = CreateInt(nI, nC, mv+1);
	for (r=I=0; I < nI; I++) for (C=0; C < nC; C++) for (v=0; v <= mv; v++) 
	{
		vStart[I][C][v] = r;
		for (J=0; J <= mJ; J++) if (nData[C][I][v][J] != 0.0) r++;
	}
	beginInsertRows(QModelIndex(), 0, r-1);
	numRows = r;
	Data = nData;
	numComp = nC;
	numIso = nI;
	numv = mv + 1;
	numJ = mJ + 1;
	if (nMixC > 0)
	{
		beginInsertColumns(QModelIndex(), 6, 5 + nMixC);
		numStates = nMixC;
		MixCoeff = newMixC;
		endInsertColumns();
	}
	endInsertRows();
}

void TermData::setIso(IsoTab *nIso)
{
	if (Iso != 0) delete Iso;
	Iso = nIso;
}

void TermData::setIsoZ(int *nZ)
{
	if (Z != 0) delete[] Z;
	Z = nZ;
}

IsoTab *TermData::getIso()
{
	return Iso;
}
