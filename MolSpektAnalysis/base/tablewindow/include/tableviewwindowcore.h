//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_VIEW_WINDOW_CORE_H
#define TABLE_VIEW_WINDOW_CORE_H


#include <QAbstractTableModel>


struct BaseData;

class TermEnergy;
class Spektrum;
class Molecule;

class QTextStream;

class TableViewWindowCore : public QAbstractTableModel
{
	public:
		TableViewWindowCore(Molecule* mol = nullptr, QObject *parent = 0, QRegExp readSpecialPartRegex = QRegExp());
		~TableViewWindowCore();
		QString readData(QTextStream& S);
		void writeData(QTextStream& S);
		int rowCount(const QModelIndex &parent = QModelIndex()) const override;
		void setRowCount(const int count);
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
		int getMaxJ();
		int getMaxv();
		virtual void addRow(BaseData* const data);
		virtual void addRow(const QStringList& L);
		virtual void setRow(BaseData* const data, const int row);
		virtual void setRow(const QStringList& L, const int row);
		virtual BaseData* getRow(const int row) const;
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
		void shrinkAllSpectRefs();
		virtual BaseData* convertToBaseData(const QStringList& L) const;
        void EmitDataChanged(const int row, const int column);

		inline BaseData* getData(const int row) const
		{
			return mData[row];
		}

		inline std::vector<BaseData*> getData()
		{
			return mData;
		}

		inline void setData(const int row, BaseData * const data)
		{
			mData[row] = data;
		}

		inline QModelIndex getIndex(const int row, const int column) const
		{
			return createIndex(row, column);
		}

		inline Molecule* getMolecule() const
		{
			return molecule;
		}

	protected:
		std::vector<BaseData*> mData;
		const QRegExp mStartSpecialPart;
		QPixmap *NewPix;
		Molecule* molecule;
};

#endif
