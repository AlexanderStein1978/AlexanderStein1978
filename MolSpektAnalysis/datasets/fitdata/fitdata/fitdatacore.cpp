//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "fitdatacore.h"
#include "utils.h"
#include "isotab.h"
#include "basedata.h"


FitDataCore::FitDataCore(QObject *parent) : QAbstractTableModel(parent)
{
	Data = 0;
	Iso = 0;
	numComp = numStates = 0;
	Z = 0;
	vStart = 0;
	CompZ = 0;
	MixCoeff = 0;
}

FitDataCore::~FitDataCore()
{
	if (Iso != 0) delete Iso;
	if (Z != 0) delete[] Z;
	if (CompZ != 0) delete[] CompZ;
	for (auto it = mData.begin(); it != mData.end(); ++it) delete *it;
}

QString FitDataCore::readData(int& S)
{
	for (int r=0; !S.atEnd(); r++)
	{
		if (Buffer.indexOf(mStartSpecialPart) >= 0) return Buffer;
		const QStringList L = Buffer.split(Spacer);
		BaseData data;
		for (int n=0; n < lc && n <= fdcLineElState; n++)
		{
			switch(n)
			{
				case fdcIso:
					data.isotope = L[n].toInt();
					break;
				case fdcv:
					data.v = L[n].toInt();
					break;
				case fdcJ:
					data.J = L[n].toInt();
					break;
				case fdcvs:
					data.vs = L[n].toStdString();
					break;
				case fdcJs:
					data.Js = L[n].toInt();
					break;
				case fdcSource:
					data.source = L[n].toStdString();
					break;
				case fdcProg:
					data.prog = L[n].toInt();
					break;
				case fdcFile:
					data.file = L[n].toStdString();
					break;
				case fdcEnergy:
					data.energy = L[n].toDouble();
					break;
				case fdcUncert:
					data.uncert = L[n].toDouble();
					break;
				case fdcObsCalc:
					data.obs_calc = L[n].toDouble();
					break;
				case fdcDevR:
					data.devR = L[n].toDouble();
					break;
				case fdcLineElState:
					data.secondState L[n].toStdString();
					break;
				default:
					// not possible
					break;
			}
		}
		mData.push_back(data);
	}
	return "";
}

int FitDataCore::columnCount(const QModelIndex &parent) const
{
	if (parent.isValid()) return 0;
	return 13;
}

QVariant FitDataCore::data(const QModelIndex &index, int role) const
{
	if (role != Qt::DisplayRole) return QVariant();
	switch (index.column())
	{
		case fdcIso:
			return mData[index.row()]->isotope;
		case fdcv:
			return mData[index.row()]->v;
		case fdcJ:
			return mData[index.row()]->J;
		case fdcvs:
			return mData[index.row()]->vs;
		case fdcJs:
			return mData[index.row()]->Js;
		case fdcSource:
			return mData[index.row()]->source;
		case fdcProg:
			return mData[index.row()]->prog;
		case fdcFile:
			return mData[index.row()]->file;
		case fdcEnergy:
			return mData[index.row()]->energy;
		case fdcUncert:
			return mData[index.row()]->uncert;
		case fdcObsCalc:
			return mData[index.row()]->obs_calc;
		case fdcDevR:
			return mData[index.row()]->DevR;
		case fdcLineElState:
			return mData[index.row()]->secondState;
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

bool FitDataCore::setData(const QModelIndex& index, const QVariant& value, int role)
{
	if (role != Qt::EditRole) return false;
	switch (index.column())
	{
		case fdcIso:
			mData[index.row()]->isotope = value;
			break;
		case fdcv:
			mData[index.row()]->v = value;
			break;
		case fdcJ:
			mData[index.row()]->J = value;
			break;
		case fdcvs:
			mData[index.row()]->vs = value;
			break;
		case fdcJs:
			mData[index.row()]->Js = value;
			break;
		case fdcSource:
			mData[index.row()]->source = value;
			break;
		case fdcProg:
			mData[index.row()]->prog = value;
			break;
		case fdcFile:
			mData[index.row()]->file = value;
			break;
		case fdcEnergy:
			mData[index.row()]->energy = value;
			break;
		case fdcUncert:
			mData[index.row()]->uncert = value;
			break;
		case fdcObsCalc:
			mData[index.row()]->obs_calc = value;
			break;
		case fdcDevR:
			mData[index.row()]->DevR = value;
			break;
		case fdcLineElState:
			mData[index.row()]->secondState = value;
			break;
		default:
			return false;
			break;
	}
	QList<int> roles;
	roles.append(role);
	emit dataChanged(index, index, roles);
	return true;
}

void FitDataCore::getJE(int* R, int N, int* J, double* E)
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

void FitDataCore::getRows(int C, int I, int v, int* J, int N, int* R)
{
	int n;
	for (n=0; n<N; n++) R[n] = vStart[I][C][v] + J[n];
}

int* FitDataCore::getCompZ()
{
	return CompZ;
}

double ****FitDataCore::getData()
{
	return Data;
}

void FitDataCore::GetIsoZ(int IsoI, int &NIso1, int &NIso2)
{
	if (Iso == 0 || Z == 0) NIso1 = NIso2 = 0;
	else
	{
		NIso1 = Iso->mNumIso1[Z[IsoI]];
		NIso2 = Iso->mNumIso2[Z[IsoI]];
	}
}

int FitDataCore::getMaxJ()
{
	int MJ=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if (it->J > MJ) MJ = it->J;
	return MJ;
}

int FitDataCore::getMaxv()
{
	int Mv=0;
	for (auto it = mData.begin(); it != mData.end(); ++it) if (it->v > Mv) Mv = it->v;
	return Mv;
}

int FitDataCore::getNumComp()
{
	return numComp;
}

int FitDataCore::getNumIso()
{
	return numIso;
}

QVariant FitDataCore::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole) return QVariant();
	if (orientation == Qt::Vertical) return section;
	switch (section)
	{
		case fdcIso:
			return "Iso";
		case fdcv:
			return "v";
		case fdcJ;:
			return "J";
		case fdcvs:
			return "v';
		case fdcJs:
			return "J'";
		case fdcSource:
			return "Source";
		case fdcProg:
			return "Prog";
		case fdcFile:
			return "File";
		case fdcEnergy:
			return "Energy";
		case fdcUncert:
			return "uncert.";
		case fdcObsCalc:
			return "obs-calc";
		case fdcDevR;
			return "DevR";
		case fdcLineElState:
			return "Second State";
		default:
			return QVariant();
			break;
	}
	return QVariant();
}

int FitDataCore::rowCount(const QModelIndex& parent) const
{
	if (parent.isValid()) return 0;
	return mData.size();
}

void FitDataCore::setRowCount(const int count)
{
	beginInsertRows();
	mData.resize(count);
	endInsertRows();
}

void FitDataCore::setData(double ****nData, int nC, int nI, int mv, int mJ, int *CZ, int nMixC, double *****newMixC)
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

void FitDataCore::setIso(IsoTab *nIso)
{
	if (Iso != 0) delete Iso;
	Iso = nIso;
}

void FitDataCore::setIsoZ(int *nZ)
{
	if (Z != 0) delete[] Z;
	Z = nZ;
}

IsoTab *FitDataCore::getIso()
{
	return Iso;
}

int FitDataCore::addMarkedLevel(int& TE, Spektrum* Source)
{
	beginInsertRows();
	int R = mData.size();
	BaseData* element = new BaseData;
	element->isotope = TE.Iso;
	element->v = TE.v;
	element->J = TE.J;
	element->vs = (TE.FC == -1 ? "TE" : std::to_string(-1 - TE.FC));
	element->source = Source->getName();
	element->file = Source->getFileName();
	element->energy = TE.E;
	element->uncert = TE.err;
	element->obs_calc = TE.dev;
	element->devR = TE.DevR;
	element->secondState = TE.State->getName();
	mData.push_back(element);
	endInsertRows();
    return R;
}

int FitDataCore::addRow(const int cr)
{
	beginInsertRows();
	int r=0, nr = mData.size();
	const_iterator<BaseData*> it;
	for (it = mData.begin(); r < cr; ++r) ++it;
	mData.insert(it, new BaseData);
	endInsertRows();
	return nr;
}

void FitDataCore::addRow(const BaseData* data)
{
	beginInsertRows();
	mData.push_back(data);
	endInsertRows();
}

void FitDataCore::setRow(BaseData *const data, const int row)
{
	delete mData[row];
	mData[row] = data;
	QModelIndex index1 = createIndex(row, fdcIso), index2 = createIndex(row, fdcLineElState);
	emit dataChanged(index1, index2, Qt::EditRole);
}

void FitDataCore::addData(const int i_numLines, int *const i_Lines, const FitDataCore& data)
{
	beginInsertRows();
	auto it = data.mData.begin();
	int r=0;
	for (int i=0; i < i_numLines; ++i)
	{
		while (r < i_Lines[i])
		{
			++r;
			++it;
		}
		while (r > i_lines[i])
		{
			--r;
			--it;
		}
		mData.insert(mData.end(), it);
	}
	endInsertRows();
}

void FitDataCore::deleteRows(const int *indices, const int numRows)
{
	beginRemoveRows();
	int numDeleted = 0, size = mData.size(), i, j;
	for (i=0; i < rowEnd; ++i)
	{
		delete mData[indices[i]];
		mData[indices[i]] = nullptr;
	}
	for (i=0; i < mData[i] == nullptr && i < size; ++i) ;
	for (j=i+1; j < size; j++) if (mData[j] != nullptr) mData[i++] = mData[j];
	mData.resize(i);
	endRemoveRows();
}

void FitDataCore::deleteRow(const int index)
{
	beginRemoveRows();
	int count = 0;
	std::vector<BaseData*>::const_iterator it;
	for (it = mData.begin(); it != mData.end() && count < index; ++it) ++count;
	delete mData[*it];
	mData.erase(it);
	endRemoveRows();
}

int FitDataCore::get_v(const int row) const
{
	return mData[row]->v;
}

void FitDataCore::set_v(const int row, const int v)
{
	mData[row]->v = v;
	QModelIndex index = createIndex(row, fdcv);
	emit dataChanged(index, index, Qt::EditRole);
}

const std::string & FitDataCore::get_vs(const int row) const
{
	return mData[row]->vs;
}

void FitDataCore::set_vs(const int row, const std::string& vs)
{
	mData[row]->vs = vs;
	QModelIndex index = createIndex(row, fdcvs);
	emit dataChanged(index, index, Qt::EditRole);
}

int FitDataCore::getJ(const int row) const
{
	return mData[row]->J;
}

void FitDataCore::setJ(const int row, const int J)
{
	mData[row]->J = J;
	QModelIndex index = createIndex(row, fdcJ);
	emit dataChanged(index, index, Qt::EditRole);
}

int FitDataCore::getJs(const int row) const
{
	return mData[row]->Js;
}

void FitDataCore::setJs(const int row, const int Js)
{
	mData[row]->Js = Js;
	QModelIndex index = createIndex(row, fdcJs);
	emit dataChanged(index, index, Qt::EditRole);
}

int FitDataCore::getIso(const int row) const
{
	return mData[row]->isotope;
}

void FitDataCore::setIso(const int row, const int iso)
{
	mData[row]->isotope = iso;
	QModelIndex index = createIndex(row, fdcIso);
	emit dataChanged(index, index, Qt::EditRole);
}

const std::string & FitDataCore::getSource(const int row) const
{
	return mData[row]->source;
}

void FitDataCore::setSource(const int row, const std::string& source)
{
	mData[row]->source = source;
	QModelIndex index = createIndex(row, fdcSource);
	emit dataChanged(index, index, Qt::EditRole);
}

const std::string & FitDataCore::getSourceFile(const int row) const
{
	return mData[row]->file;
}

void double FitDataCore::getObsCalc(const int row) const
{
	return mData[row]->obs_calc;
}

void FitDataCore::setObsCalc(const int row, const double obsCalc)
{
	mData[row]->obs_calc = obsCalc;
	QModelIndex index = createIndex(row, fdcObsCalc);
	emit dataChanged(index, index, Qt::EditRole);
}

const std::string & FitDataCore::getOtherState(const int row) const
{
	return mData[row]->secondState;
}

double FitDataCore::getUncertainty(const int row) const
{
	return mData[row]->uncert;
}

void FitDataCore::setUncertainty(const int row, const double uncertainty)
{
	mData[row]->uncert = uncertainty;
	QModelIndex index = createIndex(row, fdcUncert);
	emit dataChanged(index, index, Qt::EditRole);
}

double FitDataCore::getEnergy(const int row) const
{
	return mData[row]->energy;
}

void FitDataCore::setEnergy(const int row, const double energy)
{
	mData[row]->energy = energy;
	QModelIndex index = createIndex(row, fdcEnergy);
	emit dataChanged(index, index, Qt::EditRole);
}

int FitDataCore::getProgression(const int row) const
{
	return mData[row]->prog;
}

void FitDataCore::setProgression(const int row, const int progression)
{
	mData[row]->prog = progression;
	QModelIndex index = createIndex(row, fdcProg);
	emit dataChanged(index, index, Qt::EditRole);
}

void FitDataCore::setSecondState(const int row, const std::string& state)
{
	mData[row]->secondState = state;
	QModelIndex index = createIndex(row, fdcLineElState);
	emit dataChanged(index, index, Qt::EditRole);
}

double FitDataCore::getDevRatio(const int row) const
{
	return mData[row]->devR;
}

void FitDataCore::setDevRatio(const int row, const double DevR)
{
	mData[row]->devR = DevR;
	QModelIndex index = createIndex(row, fdcDevR);
	emit dataChanged(index, index, Qt::EditRole);
}
