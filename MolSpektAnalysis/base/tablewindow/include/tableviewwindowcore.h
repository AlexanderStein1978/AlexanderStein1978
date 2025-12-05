//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef TABLE_VIEW_WINDOW_CORE_H
#define TABLE_VIEW_WINDOW_CORE_H


struct BaseData;

class TermEnergy;
class Spektrum;
class Molecule;

class QTextStream;

class TableViewWindowCore : public QAbstractTableModel
{
	public:
		TableViewWindowCore(Molecule* mol = nullptr, QObject *parent = 0, QString readSpecialPartRegex = "");
		~TableViewWindowCore();
		QString readData(QTextStream& S);
		void writeData(QTextStream& S);
		int rowCount(const QModelIndex &parent = QModelIndex()) const override;
		void setRowCount(const int count);
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		std::vector<BaseData*> getData();
		bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
		int getMaxJ();
		int getMaxv();
		int addRow(const int cr);
		virtual void addRow(BaseData* const data);
		virtual void addRow(const QStringList& L);
		virtual void setRow(BaseData* const data, const int row);
		virtual void setRow(const QStringList& L, const int row);
		virtual BaseData* getRow(const int row) const;
		void deleteRow(const int index);
		void deleteRows(const int *indices, const int numRows);
		const std::string& get_vs(const int row) const;
		void set_vs(const int row, const std::string& vs);
		int getJs(const int row) const;
		void setJs(const int row, const int Js);
		int getIso(const int row) const;
		void setIso(const int row, const int iso);
		const std::string& getSourceFile(const int row) const;
		void setSourceFile(const int row, const std::string& filename);
		int getProgression(const int row) const;
		void setProgression(const int row, const int progression);
		double getUncertainty(const int row) const;
		void setUncertainty(const int row, const double uncertainty);
		double getObsCalc(const int row) const;
		void setObsCalc(const int row, const double obsCalc);
		void setIsoIcon(const int row, const QPixmap* const Icon);
		void setMolecule(Molecule* const mol);
		void shrinkAllSpectRefs();
		virtual BaseData* convertToBaseData(const QStringList& L) const;
        void EmitDataChanged(QModelIndex& index);

		inline BaseData* getData(const int row) const
		{
			return mData[row];
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
