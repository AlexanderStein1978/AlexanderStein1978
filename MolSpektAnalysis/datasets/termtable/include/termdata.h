//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TERMDATA_H
#define TERMDATA_H


#include <QAbstractTableModel>

class IsoTab;


class TermData : public QAbstractTableModel
{
	public:
		TermData(QObject *parent = 0);
		~TermData();
		int rowCount(const QModelIndex &parent = QModelIndex()) const;
		int columnCount(const QModelIndex &parent = QModelIndex()) const;
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;
		double ****getData();
		void setData(double ****nData, int numComp, int numIso, int maxv, int maxJ,
					 int *CompZ = 0, int nStates = 0, double *****nMixKoeff = 0);
		void setIso(IsoTab *Iso);
		IsoTab *getIso();
		void setIsoZ(int *Z);
		void GetIsoZ(int IsoI, int &NIso1, int &NIso2);
		int *getCompZ();
		int getNumComp();
		int getNumIso();
		int getMaxJ();
		int getMaxv();
		void getRows(int C, int I, int v, int *J, int N, int *R);
		void getJE(int *R, int N, int *J, double *E);
		
		inline int *getIsoZ()
		{
			return Z;
		}
		
		inline int getnumStates()
		{
			return numStates;
		}
		
		inline double *****getMixCoeff()
		{
			return MixCoeff;
		}
		
	private:
		int numRows, numv, numJ, numIso, numComp, *Z, ***vStart, *CompZ, numStates;
		double ****Data, *****MixCoeff;
		IsoTab *Iso;
};

#endif
