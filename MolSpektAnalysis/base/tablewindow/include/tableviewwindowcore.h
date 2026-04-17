//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_VIEW_WINDOW_CORE_H
#define TABLE_VIEW_WINDOW_CORE_H


#include <QAbstractTableModel>
<<<<<<< HEAD
#include <QRegularExpression>
=======
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
#include <cmath>

#include "basedata.h"

class TermEnergy;
class Spektrum;
class Molecule;

class QTextStream;

class TableViewWindowCore : public QAbstractTableModel
{
	public:
<<<<<<< HEAD
		TableViewWindowCore(Molecule* mol = nullptr, QObject *parent = 0, QRegularExpression readSpecialPartRegex = QRegularExpression());
=======
		TableViewWindowCore(Molecule* mol = nullptr, QObject *parent = 0, QRegExp readSpecialPartRegex = QRegExp());
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
		~TableViewWindowCore();
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		void writeData(QTextStream& S);
		int rowCount(const QModelIndex &parent = QModelIndex()) const override;
		void setRowCount(const int count);
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;
		int getMaxJ();
		int getMaxv();
		void addRow(BaseData* const data);
		void addRow(const QStringList& L);
		void setRow(BaseData* const data, const int row);
		void setRow(const QStringList& L, const int row);
		BaseData* getRow(const int row) const;
		void insertRows(const int startRow, const std::vector<BaseData*>& newRows);
		void deleteRow(const int index);
		void deleteRows(const int *indices, const int numRows);
		int get_v(const int row) const;
		int getJ(const int row) const;
		int getJs(const int row) const;
		virtual void setJs(const int, const int){};
		int getIso(const int row) const;
		virtual void setIso(const int, const int){};
		const QString& getSourceFile(const int row) const;
		void setSourceFile(const int, const QString&){};
		int getProgression(const int row) const;
		void setProgression(const int row, const int progression);
		double getEnergy(const int row) const;
		double getUncertainty(const int row) const;
		void setUncertainty(const int row, const double uncertainty);
		double getObsCalc(const int row) const;
		void setObsCalc(const int, const double){};
		void setIsoIcon(const int row, const QPixmap* const Icon);
		void setMolecule(Molecule* const mol);
		virtual void shrinkAllSpectRefs();
        void EmitDataChanged(const int row, const int column);
		void EmitDataChanged(const int upperRow, const int lowerRow, const int leftColumn, const int rightColumn);
		virtual QString cellToString(const int, const int) const;
		void RemoveEmptyRows();
		void setData(const int row, BaseData * const data);
		void setHorizontalHeader(const QStringList &Labels);
		void setVerticalHeader(const QStringList &Labels);
		void sortTab(int* S2);

		inline virtual QString readData(QTextStream&)
		{
			return QString();
		}

		inline virtual BaseData* convertToBaseData(const QStringList& L) const
		{
			return new BaseData;
		}

		inline int getNSources() const
		{
			return NSources;
		}

		inline BaseData* getData(const int row) const
		{
			return mData[row];
		}

		inline std::vector<BaseData*> getData()
		{
			return mData;
		}

		inline QModelIndex getIndex(const int row, const int column) const
		{
			return createIndex(row, column);
		}

		inline Molecule* getMolecule() const
		{
			return molecule;
		}

		inline static int getNumDecimalPlaces(const double uncertainty)
		{
			return 2 - static_cast<int>(log10(uncertainty));
		}

	protected:
		void startSearch(int& N, int *& Rows) const;
		void finishSearch(int *const Rows, const QModelIndexList& Result) const;

		QStringList horizontalHeders, verticalHeaders;
		int NSources = 0;
		std::vector<BaseData*> mData;
<<<<<<< HEAD
		const QRegularExpression mStartSpecialPart;
=======
		const QRegExp mStartSpecialPart;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
		QPixmap *NewPix;
		Molecule* molecule;
};

#endif
