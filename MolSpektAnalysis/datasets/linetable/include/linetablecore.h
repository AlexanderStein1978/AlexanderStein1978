//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef FITDATACORE_H
#define FITDATACORE_H


#include "tableviewwindowcore.h"
#include "linetablebasedata.h"

#include <vector>
#include <cmath>


class TermEnergy;
class Spektrum;
class Molecule;
struct vsOListElement;
struct Marker;

class QTextStream;


class LineTableCore : public TableViewWindowCore
{
	public:

		enum TableCols {CPN, Cvs, CJs, Cvss, CJss, CF, CWN, Cerr, CIso, CFile, CSNR, CDev, CC, CFCF, CEUp, CEav, CEUma, CEdJ, CCalc, COmC};
		static const int TableNormCols = 12;

		static std::vector<TableCols> convertHeaderStringsToColumnVector(const QStringList& headerStrings);
		static QStringList convertColumnVectorToHeaderStrings(const std::vector<TableCols>& headerVector);
		static QString convertHeaderEnumToHeaderString(const TableCols enumValue);

		LineTableCore(LineTable* lt, Molecule* mol = nullptr, QObject *parent = 0);
		~LineTableCore();
		QString readData(QTextStream& S) {};
		bool readData(const int numLines, QString& Buffer, const bool FCA, int& MaxPN, const int d);
		void writeData(QTextStream& S);
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
		void setData(const std::vector<LineTableCore::TableCols>& columns, const std::vector<BaseData*> data);
		void setColumnCount(const int count);
		BaseData* convertToBaseData(const QStringList& L) const override;
		LineTableBaseData* convertStringArrayToLineTableBaseData(const QString* const array, const int length) const;
		void WriteTFGS(QTextStream& S, vsOListElement *vsOList, const int vs0);
		void Assign_vs(double ****const UD, const double AT, const int NumC, const int mvs, const int* const SO);
		void AssignFC(const int *const LO, const int *const XIT, const int* const EIT, const int ENv, const int ENJ, double ****EData, const int XNv, const int XNJ,
					  double ****XData, const int XNC, const int cTX[], const int cTE[], const int cTL[]);
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;
		int getMarkedLine(const Marker& marker, const QString SFN) const;
		void AddMarked(Marker** Lines, const Marker* const LaserLine, const int AnzahlMarker, int& Offset, int &NpProg, const int MaxPN, const QString& SpektFile);
		void ShowUpTerm(const int* const SO, const int MCT, const int* const CT, const int* const IsoT, double ****ELU, const int Mv, const int MJ, double ****UT, const int mvs,
							   const int Jeo, const int* const UIsoT, const float S, const int UNC, const int UMCT, const int* const UCT);
		void ShowUpTermtable(const int * const SO, int &n, QString** Data);
		void RemoveDoubled(const int *const sortArray);
		void setUncertainty(const int row, const double Error, const double OvError);
		void set_vs(const int row, const int vs);
		void set_vss(const int row, const int vss);
		void setJss(const int row, const int Jss);
		void setFineStructureQN(const int row, const int FQN);

		void setUncertainty(const int row, const double uncertainty)
		{
			TableViewWindowCore::setUncertainty(row, uncertainty);
		}

		inline double getWaveNumber(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->waveNumber;
		}

		inline double getSNR(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->SNR;
		}

		inline int get_vs(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->vs;
		}

		inline int get_vss(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->vss;
		}

		inline int getJss(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->Jss;
		}

		inline void setJss(const int row, const int Jss)
		{
			reinterpret_cast<LineTableBaseData*>(mData[row])->Jss = Jss;
		}

		inline int getFineStructureQN(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->F;
		}

		inline void setFineStructureQN(const int row, const int FQN)
		{
			reinterpret_cast<LineTableBaseData*>(mData[row])->F = FQN;
		}

		inline double getUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->upperEnergy;
		}

		inline double getAverageUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->averageUpperEnergy;
		}

		inline double getDiffToAverageUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->diffToAverageEnergy;
		}

		inline QString getComment(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->Comment;
		}

		inline void setComment(const int row, const QString& comment)
		{
			reinterpret_cast<LineTableBaseData*>(mData[row])->Comment = comment;
		}

		inline double getFCF(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->FCF;
		}

		inline void setFCF(const int row, const double FCF)
		{
			reinterpret_cast<LineTableBaseData*>(mData[row])->FCF = FCF;
		}

		inline double getCalculatedUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->calculatedEnergy;
		}

		inline double getDeviationToCalculatedUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->diffToCalculatedEnergy;
		}

		inline double getDiffToNextJ(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->energyDiffToNextJ;
		}

		inline const std::vector<TableCols>& getColumnVector() const
		{
			return mColumns;
		}

	private:
		std::vector<TableCols> mColumns;
		const QRegExp mStartSpecialPart = QRegExp("");
		LineTable* const lTab;
};

#endif
