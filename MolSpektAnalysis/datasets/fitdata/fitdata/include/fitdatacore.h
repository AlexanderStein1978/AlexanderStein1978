//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FITDATACORE_H
#define FITDATACORE_H


#include <QAbstractTableModel>
#include <vector>

struct BaseData;

class TermEnergy;
class Spektrum;
class IsoTab;

class FitDataCore : public QAbstractTableModel
{
	public:
		enum FitDataColumn {fdcIso, fdcv, fdcJ, fdcvs, fdcJs, fdcSource, fdcProg, fdcFile, fdcEnergy,
        fdcUncert, fdcObsCalc, fdcDevR, fdcLineElState};

		FitDataCore(QObject *parent = 0);
		~FitDataCore();
		QString readData(QTextStream& S);
		int rowCount(const QModelIndex &parent = QModelIndex()) const override;
		void setRowCount(const int count);
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		std::vector<BaseData*> getData();
		bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
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
		int addMarkedLevel(TermEnergy& TE, Spektrum *Source);
		int addRow(const int cr);
		void setRow(const BaseData& data, const int row);
		void addData(const int i_numLines, int *const i_Lines, const FitDataCore& data);
		void deleteRows(const int *indices, const int numRows);
		int get_v(const int row) const;
		void set_v(const int row, const int v);
		const std::string& get_vs(const int row) const;
		void set_vs(const int row, const std::string& vs);
		int getJ(const int row) const;
		int getJs(const int row) const;
		int getIso(const int row) const;
		const std::string& getSource(const int row) const;
		void setSource(const int row, const std::string& source);
		const std::string& getSourceFile(const int row) const;
		void setSourceFile(const int row, const std::string& filename);
		int getProgression(const int row) const;
		double getEnergy(const int row) const;
		double getUncertainty(const int row) const;
		void setObsCalc(const int row, const double obsCalc);
		void setDevRatio(const int row, const double DevR);
		const std::string& getOtherState(const int row) const;
		void setSecondState(const int row, const std::string& state);

		inline int *getIsoZ() const
		{
			return Z;
		}
		
		inline int getnumStates() const
		{
			return numStates;
		}

		inline void setColumnCount(const int newCount)
		{
			mColumnCount = newCount;
		}

		inline const BaseData* getData(const int row) const
		{
			return mData[row];
		}

	private:
		int *Z, *CompZ, numStates, mColumnCount;
		std::vector<BaseData*> mData;
		IsoTab *Iso;
		const QRegExp mStartSpecialPart = QRegExp("SourceOffsets:|Begin ResidualFit");
};

#endif
