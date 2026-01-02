//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FITDATACORE_H
#define FITDATACORE_H


#include "tableviewwindowcore.h"
#include "fitdatabasedata.h"

#include <vector>

class TermEnergy;
class Spektrum;
class Molecule;

class QTextStream;

class FitDataCore : public TableViewWindowCore
{
	public:
		enum FitDataColumn {fdcIso, fdcv, fdcJ, fdcvs, fdcJs, fdcSource, fdcProg, fdcFile, fdcEnergy,
        fdcUncert, fdcObsCalc, fdcDevR, fdcLineElState};

		FitDataCore(Molecule* mol = nullptr, QObject *parent = 0);
		~FitDataCore();
		QString readData(QTextStream& S) override;
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		std::vector<FitDataBaseData*> getDataAsFitDataBaseData();
		bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
		int addMarkedLevel(TermEnergy& TE, Spektrum *Source);
		void addRow(const QStringList& L) override;
		int addRow(const int cr);
		void setRow(FitDataBaseData* const data, const int row);
		FitDataBaseData* getRowAsFitDataBaseData(const int row) const;
		void addData(const int i_numLines, int *const i_Lines, const FitDataCore& data);
		void setIso(const int row, const int v) override;
		void set_v(const int row, const int v);
		void setJ(const int row, const int J);
		const QString& get_vs(const int row) const;
		void set_vs(const int row, const QString& vs);
		void setJs(const int row, const int J) override;
		const QString& getSource(const int row) const;
		void setSource(const int row, const QString& source);
		void setSourceFile(const int row, const QString& filename);
		int getProgression(const int row) const;
		void setProgression(const int row, const int progression);
		void setEnergy(const int row, const double energy);
		void setUncertainty(const int row, const double uncertainty);
		void setObsCalc(const int row, const double obsCalc);
		float getDevRatio(const int row) const;
		void setDevRatio(const int row, const float DevR);
		const QString& getOtherState(const int row) const;
		void setSecondState(const int row, const QString& state);
		void setRWError(const QString& headerText);
		void shrinkAllSpectRefs();
		void search(const int* const Rows, const int NRows, const int column, const int value, const int smeqla, QModelIndexList& Result) const;
		void search(const int* const Rows, const int NRows, const int column, const double value, const int smeqla, QModelIndexList& Result) const;
		void search(const int* const Rows, const int NRows, const QString& Text, QModelIndexList& Result, const int column=-1,
					const bool completeCell = false) const;
		FitDataBaseData* convertToFitDataCoreBaseData(const QStringList& L) const;
		QString cellToString(const int r, const int c) const override;

		inline BaseData* convertToBaseData(const QStringList& L) const override
		{
			return convertToFitDataCoreBaseData(L);
		}

		inline void setRow(const QStringList& L, const int row) override
		{
			mData[row] = convertToFitDataCoreBaseData(L);
		}

		inline void setData(const int row, FitDataBaseData * const data)
		{
			mData[row] = data;
		}

		inline QModelIndex getIndex(const int row, const int column) const
		{
			return createIndex(row, column);
		}

		inline int getNSources() const
		{
			return NSources;
		}

		inline void setNSources(const int N)
		{
			NSources = N;
		}

		inline bool isRWErr()
		{
			return !RWError.isEmpty();
		}
		
	private:

		inline void setRow(BaseData* const data, const int row) override
		{
			mData[row] = data;
		}

		inline void setData(const int row, BaseData* const data)
		{
			mData[row] = reinterpret_cast<FitDataBaseData* const>(data);
		}

		inline void addRow(BaseData* const data) override
		{
			TableViewWindowCore::addRow(data);
		}

		int NSources = 0;
		QString RWError;
};

#endif
