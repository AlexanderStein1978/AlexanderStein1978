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

class QTextStream;


class LineTableCore : public TableViewWindowCore
{
	public:

		enum TableCols {CPN, Cvs, CJs, Cvss, CJss, CF, CWN, Cerr, CIso, CFile, CSNR, CDev, CC,
				CFCF, CEUp, CEav, CEUma, CEdJ, CCalc, COmC};
		const int TableNormCols = 12;

		static std::vector<TableCols> convertHeaderStringsToColumnVector(const QStringList& headerStrings);
		static QStringList convertColumnVectorToHeaderStrings(const std::vector<TableCols>& headerVector);

		LineTableCore(LineTable* lt, Molecule* mol = nullptr, QObject *parent = 0);
		~LineTableCore();
		QString readData(QTextStream& S) {};
		bool readData(const int numLines, QString& Buffer, const bool FCA, int& MaxPN, const int d);
		void writeData(QTextStream& S);
		QVariant data (const QModelIndex &index, int role = Qt::DisplayRole) const override;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;
		bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
		void setColumnCount(const int count);
		BaseData* convertToBaseData(const QStringList& L) const override;
		LineTableBaseData* convertStringArrayToLineTableBaseData(const QString* const array, const int length) const;
		void WriteTFGS(QTextStream& S, vsOListElement *vsOList, const int vs0);
		void Assign_vs(double ****const UD, const double AT, const int NumC, const int mvs, const int* const SO);
		void AssignFC(const int *const LO, const int *const XIT, const int* const EIT, const int ENv, const int ENJ, double ****EData, const int XNv, const int XNJ,
					  double ****XData, const int XNC, const int cTX[], const int cTE[], const int cTL[]);
		int columnCount(const QModelIndex &parent = QModelIndex()) const override;

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

		inline int getFineStructureQN(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->F;
		}

		inline double getUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->upperEnergy;
		}

		inline double getAverageUpperEnergy(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->averageUpperEnergy;
		}

		inline QString getComment(const int row) const
		{
			return reinterpret_cast<LineTableBaseData*>(mData[row])->Comment;
		}

		inline void setFCF(const int row, const double FCF)
		{
			reinterpret_cast<LineTableBaseData*>(mData[row])->FCF = FCF;
		}

	private:
		std::vector<TableCols> mColumns;
		const QRegExp mStartSpecialPart = QRegExp("");
		LineTable* const lTab;
};

#endif
