//
// C++ Implementation: termtable
//
// Description: 
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "termtable.h"
#include "molecule.h"
#include "utils.h"
#include "constants.h"

#include <math.h>

#include <QMessageBox>
#include <QTextStream>
#include <QLabel>
#include <QFile>
#include <QFileDialog>
#include <QHeaderView>
#include <QGridLayout>
#include <QListWidget>
#include <QRadioButton>

TermTable::TermTable(MainWindow *MW) : TableWindow(TermEnergyTable, MW)
{
	State = 0;
	NPerturbations = 0;
	Perturbations = 0;
	tDat = new TermData();
	table->setModel(tDat);
	table->horizontalHeader()->setVisible(true);
	table->verticalHeader()->setVisible(true);
	setWindowTitle("New term energy table");
	connect(table, SIGNAL(SelChanged()), this, SIGNAL(SelChanged()));
}

TermTable::~TermTable()
{
	delete tDat;
	if (Perturbations != 0) delete[] Perturbations;
}

double ****TermTable::getData()
{
	return tDat->getData();
}

int TermTable::getMaxJ()
{
	return tDat->getMaxJ();
}

int TermTable::getMaxv()
{
	return tDat->getMaxv();
}

int TermTable::getNPerturbations(int c, int Iso, int v)
{
	if (c==-1 || Iso == -1 || v==-1 || NPerturbations == 0) return NPerturbations;
	int n, m;
	for (n=m=0; n < NPerturbations; n++) 
		if (Perturbations[n].Comp == c && Perturbations[n].Iso == Iso 
			&& Perturbations[n].v == v) m++;
	return m;
}

int TermTable::getNumComp()
{
	return tDat->getNumComp();
}

int TermTable::getNumIso()
{
	return tDat->getNumIso();
}

Perturbation* TermTable::getPerturbation(int n)
{
	if (Perturbations[n].Perturber == 0)
	{
		int s, NS, t, NT, m;
		QString N;
		ElState *S;
		for (s=0, NS = molecule->getNumStates(); s < NS && N != Perturbations[n].PName; s++)
		{
			NT = (S = molecule->getStateP(s))->getNumTermTables();
			for (t=0; (t < NT ? 
				(N = S->getTermTableName(t)) != Perturbations[n].PName : false); t++) ;
			if (N == Perturbations[n].PName) Perturbations[n].Perturber = S->getTermTable(t);
		}
		for (m=0; m < NPerturbations; m++) 
			if (Perturbations[m].PName == Perturbations[n].PName)
				Perturbations[m].Perturber = Perturbations[n].Perturber;
	}
	return Perturbations + n;
}

Perturbation* TermTable::getPerturbation(int c, int Iso, int v, int Number)
{
	int n, m;
	for (n=0, m=-1; (n < NPerturbations ? m < Number : false); n++)
		if (Perturbations[n].Comp == c && Perturbations[n].Iso == Iso
			&& Perturbations[n].v == v) m++;
	if (m < Number) return 0;
	return getPerturbation(--n);
}

bool TermTable::readData(QString FileName)
{
	QFile Datei(FileName);
	if (!read(&Datei))  return false;
	QTextStream stream( &Datei );
	QString Buffer = stream.readLine();
	if (Buffer.left(3) == "AZ=") 
	{
		Datei.close();
		if (!OpenXnfitTO(FileName)) 
		{
			QString Fehlermeldung = "Die Datei " + FileName + " konnte nicht geöffnet werden.";
			QMessageBox::information( this, "Application name", 
									  "Fehler beim Öffnen der Datei.", Fehlermeldung);
			return false;
		}
		Filename = FileName;
		setName(Filename);
		setImported();
		return true;
	}
	int *Iso = new int[MaxTermRows], *Ju = new int[MaxTermRows], *CompZ;
	int *vu = new int[MaxTermRows], i=0, j=1, k, l, m=0, n, iBuff, IsoM=0, JuM=0, vuM=0;
	int RC=0, *Iso2 = new int[MaxTermRows], *Comp = new int[MaxTermRows], cc=0, iB, nMixC=0;
	double *TE = new double[MaxTermRows], **MixC = 0;
	bool AI;
	bool Success = true;
	QStringList SList;
	if (Buffer.indexOf("Source", 0, Qt::CaseInsensitive) != -1)
	{
		//printf("TermTable::readData, new type\n");
		QStringList SL;
		int nT = 0;
		int *Piso1 = 0, *Piso2 = 0;
		for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
		setSource(Buffer.right(Buffer.length() - n));
		while (!stream.atEnd() && RC < MaxTermRows)
		{
			Buffer = stream.readLine();
			switch (nT)
			{
				case 0:
					if (Buffer.indexOf("Name", 0, Qt::CaseInsensitive) != -1)
					{
						for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
						setName(Buffer.right(Buffer.length() - n));
					}
					else if (Buffer.indexOf("Iso 1:", 0, Qt::CaseInsensitive) != -1) 
					{
						for (nMixC = n = 0; (n = Buffer.indexOf("Mix C S", n+1)) != -1; nMixC++) ;
						if (nMixC > 0) MixC = Create(MaxTermRows, nMixC);
						else MixC = 0;
						nT++;
					}
					break;
				case 1:
					if ((n = Buffer.indexOf("|")) > 0)
					{
						for (m=Buffer.indexOf("|", n+1); 
										Buffer[m].isSpace() || Buffer[m] == '|'; m--) ;
						while (Buffer[n].isSpace() || Buffer[n] == '|') n--;
						Iso[RC] = Buffer.left(n + 1).toInt();
						for (n++; Buffer[n].isSpace() || Buffer[n] == '|'; n++) ;
						Iso2[RC] = Buffer.mid(n, m - n + 1).toInt();
						for (m++; Buffer[m].isSpace() || Buffer[m] == '|'; m++) ;
						for (n=Buffer.indexOf("|", m); 
										Buffer[n].isSpace() || Buffer[n] == '|'; n--) ;
						Comp[RC] = iB = Buffer.mid(m, n - m + 1).toInt();
						if (iB >= cc) cc = iB + 1;
						for (n++; Buffer[n].isSpace() || Buffer[n] == '|'; n++) ;
						for (m=Buffer.indexOf("|", n); 
										Buffer[m].isSpace() || Buffer[m] == '|'; m--) ;
						vu[RC] = iB = Buffer.mid(n, m - n + 1).toInt();
						if (iB >= vuM) vuM = iB+1;
						for (m++; Buffer[m].isSpace() || Buffer[m] == '|'; m++) ;
						for (n=Buffer.indexOf("|", m); Buffer[n].isSpace() 
										|| Buffer[n] == '|'; n--) ;
						Ju[RC] = iB = Buffer.mid(m, n - m + 1).toInt();
						if (iB >= JuM) JuM = iB+1;
						for (n++; Buffer[n].isSpace() || Buffer[n] == '|'; n++) ;
						SList =  Buffer.right(Buffer.length() - n).split('|');
						TE[RC] = SList[0].toDouble();
						for (n=1; n <= nMixC && n < SList.count(); n++) MixC[RC][n-1] = SList[n].toDouble();
						RC++;
						//printf("Buffer = %s\n", Buffer.ascii());
					}
					else if (Buffer.indexOf("Perturbations:", 0, Qt::CaseInsensitive) != -1) 
					{
						if (Perturbations != 0) delete[] Perturbations;
						NPerturbations = 
							Buffer.right(Buffer.length() - Buffer.indexOf('=') -1).toInt();
						Perturbations = new Perturbation[NPerturbations];
						stream.readLine();
						nT++;
						n=0;
						Piso1 = new int[NPerturbations];
						Piso2 = new int[NPerturbations];
					}
					break;
				case 2: 
					if (Buffer.indexOf(" | ") == 0 || n == NPerturbations) continue;
					SL = Buffer.split(" | ");
					if (SL.count() < 8) continue;
					Piso1[n] = SL[0].toInt();
					Piso2[n] = SL[1].toInt();
					Perturbations[n].Iso = m;
					Perturbations[n].Comp = SL[2].toInt();
					Perturbations[n].v = SL[3].toInt();
					Perturbations[n].J = SL[4].toInt();
					Perturbations[n].PName = SL[5];
					Perturbations[n].PComp = SL[6].toInt();
					Perturbations[n].Pv = SL[7].toInt();
					Perturbations[n++].Perturber = 0;
					break;
			}
		}
		if (RC >= MaxTermRows) 
		{
			QString Fehlermeldung = "The table of term energies in the File " 
									+ FileName + 
							" has more rows than the maximum amount of allowed rows of "; 
			Fehlermeldung += QString::number(MaxTermRows);
			Fehlermeldung += '!';
			QMessageBox::information( this, "MolSpektAnalysis", 
											  "Error opening file.", Fehlermeldung);
		}
		//printf("Vor Iso\n");
		if (molecule == 0)
		{
			//printf("no molecule\n");
			int iso1[100], iso2[100], ic=0;
			for (i=0, iso1[0] = 0; i<RC; i++)
			{
				j = 0;
				//printf("i=%d, RC=%d\n", i, RC);
				while (iso1[j] != Iso[i] || iso2[j] != Iso2[i]) if (++j >= ic) 
				{
					if (ic == 100)
					{
						QString Fehlermeldung = "The table of term energies in the File " +
						  FileName + 
						  " contains an unrealistic high amount of different isotopologues!";
						QMessageBox::information( this, "MolSpektAnalysis", 
								"Error opening file:" + Fehlermeldung);
						delete[] Iso;
						delete[] Iso2;
						delete[] vu;
						delete[] Ju;
						delete[] Comp;
						delete[] TE;
						return false;
					}
					iso1[ic] = Iso[i];
					iso2[ic] = Iso2[i];
					if (j > ic) j = ic;
					ic++;
				}
				Iso[i] = j;
			}
			for (i=0; i < NPerturbations; i++)
			{
				for (j=0; (j < ic ? Piso1[i] != iso1[j] || Piso2[i] != iso2[j] : false); j++) ;
				if (j < ic) Perturbations[i].Iso = j;
				else
				{
					QMessageBox::information( this, "MolSpektAnalysis", 
								"Error reading perturbation data: unrealistic isotopologue!");
					delete[] Perturbations;
					Perturbations = 0;
					NPerturbations = 0;
					delete[] Piso1;
					delete[] Piso2;
					break;
				}
			}
			IsoM = ic;
			IsoTab *IsoT = new IsoTab(ic);
			int *Z = new int[ic];
			for (i=0; i < ic; i++)
			{
				Z[i] = i;
				IsoT->mNumIso1[i] = iso1[i];
				IsoT->mNumIso2[i] = iso2[i];
			}
			tDat->setIso(IsoT);
			tDat->setIsoZ(Z);
		}
		else
		{
			//printf("molecule=%d\n", molecule);
			IsoTab *IsoT = tDat->getIso();
			IsoM = IsoT->numIso;
			int *Z = new int[IsoM];
			for (i=0; i<IsoM; i++) Z[i] = i;
			for (i=0; i<RC; i++)
			{
				for (j=0; (j < IsoM ? IsoT->mNumIso1[j] != Iso[i] 
								 || IsoT->mNumIso2[j] != Iso2[i] : false); j++) ;
				if (j == IsoM)
				{
					QString Fehlermeldung = 
							"The table of term energies in the File " + FileName + 
						" contains an isotopomer which doesn't belong to the molecule " 
							+ molecule->getName() + '!';
					QMessageBox::information( this, "MolSpektAnalysis", 
											  "Error opening file:" + Fehlermeldung);
					delete[] Iso;
					delete[] Iso2;
					delete[] vu;
					delete[] Ju;
					delete[] Comp;
					delete[] TE;
					return false;
				}
				Iso[i] = j;
			}
			for (i=0; i < NPerturbations; i++)
			{
				for (j=0; (j < IsoM ? IsoT->mNumIso1[j] != Piso1[i] || IsoT->mNumIso2[j] != Piso2[i] : false); j++) ;
				if (j == IsoM)
				{
					QMessageBox::information( this, "MolSpektAnalysis", 
								"Error reading perturbation data: unrealistic isotopologue!");
					delete[] Perturbations;
					Perturbations = 0;
					NPerturbations = 0;
					delete[] Piso1;
					delete[] Piso2;
					break;
				}
				Perturbations[i].Iso = j;
			}
			tDat->setIsoZ(Z);
		}
		AI = false;
		if (NPerturbations > 0)
		{
			delete[] Piso1;
			delete[] Piso2;
		}
	}
	else
	{
		if (Buffer.left(4) == "Iso=") 
		{
			//printf("Vor read\n");
			int I = Buffer.right(Buffer.length() - 4).toInt(), J=0;
			while (!stream.atEnd())
			{
				Buffer = stream.readLine();
				if (Buffer.left(2) == "v=")
				{
					n = Buffer.indexOf("\t");
					if ((vu[i] = Buffer.mid(2, n-2).toInt()) > vuM) vuM = vu[i];
					Ju[i] = J;
					Iso[i] = I;
					TE[i++] = Buffer.right(Buffer.length() - n - 3).toDouble();
					//if (TE[i-1] == 20.22246015) 
					//printf("i=%d, I=%d, J=%d, v=%d, T=%f\n", i-1, I, J, vu[i-1], TE[i-1]);
					if (i==MaxTermRows)
					{
						printf("TermTable::ReadData: Error, number of rows > MaxTermRows=%d\n", MaxTermRows);
						break;
					}
				}
				else if (Buffer.left(2) == "J=") 
				{
					if ((J = Buffer.right(Buffer.length() - 2).toInt()) > JuM) JuM = J;
				}
				else if (Buffer.left(4) == "Iso=") 
					if ((I = Buffer.right(Buffer.length() - 4).toInt()) > IsoM) IsoM = I;
			}
			JuM++;
			vuM++;
			IsoM++;
			//printf("JuM=%d, vuM=%d, IsoM=%d\n", JuM, vuM, IsoM);
		}
		else
		{
			//printf("Import TermTable\n");
			Buffer = stream.readLine();
			setSource(Buffer);
			setName("Imported TermE");
			iBuff = Buffer.left(5).toInt();
			while (!stream.atEnd())
			{	
				if (iBuff == 1 && Buffer.indexOf("!level") == -1)
				{
					while (iBuff == j++)
					{
						//printf("Buffer = %s, i=%d\n", Buffer.ascii(), i);
						for (m=0; Buffer[m].isSpace() && !Buffer[m].isNull(); m++) ;
						while (!Buffer[m].isSpace() && !Buffer[m].isNull()) m++;
						while (Buffer[m].isSpace() && !Buffer[m].isNull()) m++;
						for (n=m; !Buffer[n].isSpace() && !Buffer[n].isNull(); n++) ;
						vu[i] = Buffer.mid(m, n - m).toInt();
						if (vu[i] >= vuM) vuM = vu[i] + 1; 
						while (Buffer[n].isSpace() && !Buffer[n].isNull()) n++;
						for (m=n; !Buffer[m].isSpace() && !Buffer[m].isNull(); m++) ;
						Ju[i] = Buffer.mid(n, m - n).toInt();
						if (Ju[i] >= JuM) JuM = Ju[i] + 1;
						while (Buffer[m].isSpace() && !Buffer[m].isNull()) m++;
						for (n=m; !Buffer[n].isSpace() && !Buffer[n].isNull(); n++) ;
						TE[i] = Buffer.mid(m, n - m).toDouble();
						while (Buffer[n].isSpace() && !Buffer[n].isNull()) n++;
						for (m=n; !Buffer[m].isSpace() && !Buffer[m].isNull(); m++) ;
						Comp[i] = Buffer.mid(n, m - n).toInt();
						if (Comp[i] >= cc) cc = Comp[i] + 1;
						while (Buffer[m].isSpace() && !Buffer[m].isNull()) m++;
						while (!Buffer[m].isSpace() && !Buffer[m].isNull()) m++;
						while (Buffer[m].isSpace() && !Buffer[m].isNull()) m++;
						for (n=m; !Buffer[n].isSpace() && !Buffer[n].isNull(); n++) ;
						Iso[i] = (Buffer.mid(n, m - n).toInt() - 1) / 10;
						if (Iso[i] >= IsoM) IsoM = Iso[i] + 1;
						//printf("vu=%d, Ju=%d, Iso=%d, TE[%d]=%f\n", vu[i], Ju[i], 
						//Iso[i], i, TE[i]);
						i++;
						if (stream.atEnd()) break;
						Buffer = stream.readLine();
						iBuff = Buffer.left(5).toInt();
					}
					j = 1;
				}
				else
				{
					Buffer = stream.readLine();
					iBuff = Buffer.left(5).toInt();
				}
			}
			RC = i;
		}
		if (cc == 0) cc = 1;
		AI = true;
		setImported();
	}
	if (cc > 1)
	{
		int CompC[cc];
		for (n=0; n < cc; n++) CompC[n] = 0;
		for (n=0; n < RC; n++) 
		{
			if (Comp[n] < 0 || Comp[n] >= cc) 
			{
				printf("Comp[%d]=%d, cc=%d\n", n, Comp[n], cc);
				Comp[n] = 0;
				Success = false;
			}
			CompC[Comp[n]]++; 
		}
		for (n=m=0; n < cc; n++)
		{
			if (CompC[n] > 0) CompC[n] = m++;
			else CompC[n]--;
		}
		for (n=0; n < RC; n++) Comp[n] = CompC[Comp[n]];
		CompZ = new int[m];
		for (n=0; n < cc; n++) if (CompC[n] >= 0) CompZ[CompC[n]] = n;
		cc = m;
	}
	else
	{
		CompZ = new int[1];
		CompZ[0] = 0;
	}
	//printf("Vor Ende\n");
	bool IsoNotAv[IsoM];
	int *IsoZ = tDat->getIsoZ();
	for (n=0; n < IsoM; ++n) IsoNotAv[n] = true;
	for (n=0; n < RC; ++n) if (Iso[n] >= 0 && Iso[n] < IsoM && IsoNotAv[Iso[n]]) IsoNotAv[Iso[n]] = false;
	for (n=0; n < IsoM && !IsoNotAv[Iso[n]]; ++n) ;
	if (n < IsoM)
	{
		int *NewZ = new int[n];
		for (m=n=0; n < IsoM; ++n) if (!IsoNotAv[n]) NewZ[m++] = IsoZ[n];
		IsoM = m;
		tDat->setIsoZ(NewZ);
	}
	double ****Data = Create(cc, IsoM, vuM, JuM);
	double *****MixCoeff = (nMixC > 0 ? Create(cc, IsoM, vuM, JuM, nMixC) : 0);
	//printf("cc=%d, IsoM=%d, vuM=%d, JuM=%d, RC=%d\n", cc, IsoM, vuM, JuM, RC);
	for (m=0; m<cc; m++) for (j=0; j<IsoM; j++) for (k=0; k<vuM; k++) for (l=0; l<JuM; l++)
		Data[m][j][k][l] = 0.0;
	//printf("Vor Data, MaxTermRows=%d\n", MaxTermRows);
	for (j=0; j<RC; j++) 
	{
		if (Comp[j] < 0 || Comp[j] >= cc) 
		{
			printf("Comp[%d]=%d, cc=%d\n", j, Comp[j], cc);
			Comp[j] = 0;
			Success = false;
		}
		if (Iso[j] < 0 || Iso[j] >= IsoM) 
		{
			printf("Iso[%d]=%d, IsoM=%d\n", j, Iso[j], IsoM);
			Iso[j] = 0;
			Success = false;
		}
		if (vu[j] < 0 || vu[j] >= vuM) 
		{
			printf("vu[%d]=%d, vuM=%d\n", j, vu[j], vuM);
			vu[j] = 0;
			Success = false;
		}
		if (Ju[j] < 0 || Ju[j] >= JuM) 
		{
			printf("Ju[%d]=%d, JuM=%d\n", j, Ju[j], JuM);
			Ju[j] = 0;
			Success = false;
		}
		Data[Comp[j]][Iso[j]][vu[j]][Ju[j]] = TE[j];
		for (n=0; n < nMixC; n++) MixCoeff[Comp[j]][Iso[j]][vu[j]][Ju[j]][n] = MixC[j][n];
	}
    /*for (k=0; k<10; k++) for (l=0; l<JuM; l++) 
	printf("ELU[0][%d][%d]=%f\n", k, l, ELU[0][k][l]);*/
	//printf("Anzahl Zeilen=%d\n", RC);
	tDat->setData(Data, cc, IsoM, vuM - 1, JuM - 1, CompZ, nMixC, MixCoeff);
	//printf("Nach setData\n");
	Datei.close();
	if (getName().isEmpty()) setName(Buffer);
	if (molecule != 0 && AI) AssignIso();
	Saved();
	//printf("Vor  delete\n");
	delete[] vu;
	delete[] Ju;
	delete[] Iso;
	delete[] Iso2;
	delete[] Comp;
	delete[] TE;
	if (nMixC > 0) Destroy(MixC, MaxTermRows);
	//printf("Nach delete\n");
	return Success;
}

bool TermTable::OpenXnfitTO(QString Dateiname)
{
	QFile Datei(Dateiname);
	if (!Datei.open(QIODevice::ReadOnly)) return false;
	QTextStream S(&Datei);
	QString Buffer, B2;
	int AZ, i, j, J=0, v, MaxJ, Maxv, NumComp, NumIso;
	Buffer = S.readLine();
	B2 = Buffer.left((i = Buffer.indexOf(' ')));
	AZ = B2.right(B2.length() - 3).toInt();
	i++;
	B2 = Buffer.mid(i, (j = Buffer.indexOf(' ', i)) - i);
	MaxJ = B2.right(B2.length() - 3).toInt();
	B2 = Buffer.right(Buffer.length() - j - 1);
	Maxv = B2.right(B2.length() - 3).toInt();
	NumComp = NumIso = 1;
	double ****Data = Create(NumComp, NumIso, Maxv + 1, MaxJ + 1);
	for (v=0; v<=Maxv; v++) for (J=0; J<=MaxJ; J++) Data[0][0][v][J] = 0.0;
	for (i=1; i<AZ; i++)
	{
		Buffer = S.readLine();
		if (Buffer.left(2) == "J=") J = Buffer.right(Buffer.length() - 2).toInt();
		else 
		{
			v = Buffer.left(5).toInt();
			Data[0][0][v][J] = Buffer.mid(5, 18).toDouble();
		}
	}
	tDat->setData(Data, NumComp, NumIso, Maxv, MaxJ);
	if (molecule != 0) AssignIso();
	else
	{
		IsoTab *IsoT = new IsoTab(1);
		tDat->setIso(IsoT);
		int *Z = new int;
		*Z = 1;
		tDat->setIsoZ(Z);
	}
	return true;
}

void TermTable::AssignIso()
{
	//printf("TermTable: Beginn AssignIso()\n");
	if (molecule == NULL)
	{
		printf("Error: TermTable::AssignIso(): molecule == NULL!\n");
		return;
	}
	IsoTab *Iso = tDat->getIso();
	int NumIso = tDat->getNumIso();
	if (Iso->numIso < NumIso)
	{
		printf("TermTable::AssignIso: error: the table has more isotopomers than the molecule!\n");
		return;
	}
	int i, j, k, SIso[Iso->numIso], STI[NumIso], MO[NumIso];	
	tDat->GetIsoZ(0, i, j);
	//printf("i=%d, j=%d\n", i, j);
	for (k=0; k < Iso->numIso; k++) if (i == Iso->mNumIso1[k] && j == Iso->mNumIso2[k]) return;
	double IsoZ[NumIso], MIso[Iso->numIso], D, B=-1.0, ****Data = tDat->getData();
	//bool S = true, init = true;
	//printf("Nach init\n");
	int *Z = new int[NumIso];
	//printf("Nach Z\n");
	if (NumIso > 1)
	{
		for (i = 0; i < Iso->numIso; i++) 
		{
			MIso[i] = Iso->redMass[i];
			SIso[i] = i;
		}
		for (i = 0; i < NumIso; i++) 
		{
			IsoZ[i] = Data[0][i][0][2] - Data[0][i][0][1];
			STI[i] = MO[i] = i;
		}
		for (j=0; j!=-1;) 
		{
			j=-1;
			for (i=1; i < NumIso; i++) if (IsoZ[MO[i]] < IsoZ[MO[i-1]])
			{
				j = MO[i];
				MO[i] = MO[i-1];
				MO[i-1] = j;
			}
		}
		for (j=0; j!=-1;)
		{
			j=-1;
			for (i=1; i < Iso->numIso; i++) if (MIso[SIso[i]] > MIso[SIso[i-1]])
			{
				j = SIso[i];
				SIso[i] = SIso[i-1];
				SIso[i-1] = j;
			}
		}
		for (j = Iso->numIso - 1; j>=0; j--) MIso[SIso[j]] /= MIso[SIso[0]];
		for (B=0.0, j = NumIso - 1; j >= 0; j--) 
		{
			IsoZ[MO[j]] = IsoZ[MO[0]] / IsoZ[MO[j]];
			B += fabs(IsoZ[MO[j]] - MIso[SIso[j]]);
			Z[MO[j]] = SIso[j];
			//printf("j=%d, MO[j]=%d, SIso[j]=%d, Z[%d]=%d, IsoZ=%f, MIso=%f\n", 
				//   j, MO[j], SIso[j], MO[j], Z[MO[j]], IsoZ[MO[j]], MIso[SIso[j]]);
		}
		//printf("B=%f\n", B);
		while (STI[0] + NumIso < Iso->numIso)
		{
			if (STI[NumIso - 1] < Iso->numIso - 1) STI[NumIso - 1]++;
			else
			{
				for (i = NumIso - 2; (i>=0 ? STI[i] + NumIso - i == Iso->numIso : false); i--) ;
				for (STI[i++]++; i < NumIso; i++) STI[i] = STI[i-1] + 1;
			}
			for (j=STI[0] + 1; j < Iso->numIso; j++) MIso[SIso[j]] /= MIso[SIso[STI[0]]];
			for (D=0.0, j=0; j < NumIso; j++) D += fabs(IsoZ[MO[j]] - MIso[SIso[STI[j]]]);
			//printf("D=%f\n", D);
			if (D < B || B == -1)
			{
				B = D;
				for (j=0; j < NumIso; j++) Z[MO[j]] = SIso[STI[j]];
			}
		}
		//printf("B=%f, N=%d, B/(N-1)=%f\n", B, NumIso, B / (NumIso - 1));
	}
	else Z[0] = Iso->refIso;
	tDat->setIsoZ(Z);
	//printf("TermTable: Ende AssignIso()\n");
	Changed();
}

void TermTable::getCompT(int &Max, int *&Trans)
{
	int n, NC = tDat->getNumComp(), *CompZ = tDat->getCompZ();
	for (n = Max = 0; n < NC; n++) if (CompZ[n] > Max) Max = CompZ[n];
	Trans = new int[Max + 1];
	for (n=0; n <= Max; n++) Trans[n] = -1;
	for (n=0; n < NC; n++) Trans[CompZ[n]] = n;
}

int* TermTable::getIsoT()
{
	IsoTab *Iso = tDat->getIso();
	if (Iso == 0)
	{
		printf("Error: TermTable not assigned to a molecule!");
		return 0;
	}
	int n, *R = new int[Iso->numIso], NI = tDat->getNumIso(), *Z = tDat->getIsoZ();
	for (n=0; n < Iso->numIso; n++) R[n] = -1;
	for (n=0; n < NI; n++) R[Z[n]] = n;
	return R;
}

void TermTable::GetIsoZ(int IsoI, int &I1, int &I2)
{
	tDat->GetIsoZ(IsoI, I1, I2);
	/*if (IsoI < 0 || IsoI >= NumIso)
	{
		printf("GetIsoZ Error: IsoI=%d, NumIso=%d\n", IsoI, NumIso);
		I1 = I2 = 0;
		return;
	}	
	if (molecule == 0 || Z == 0) 
	{
		I1 = Tab->item(0, 0)->text().toInt();
		I2 = Tab->item(0, 1)->text().toInt();
		if (IsoI == 0) return;
		int i=0, r, M=1, I[IsoI][2];
		I[0][0] = I1;
		I[0][1] = I2;
		for (r=1; r < numVRows; r++)
		{
			I1 = Tab->item(r, 0)->text().toInt();
			I2 = Tab->item(r, 1)->text().toInt();
			if (I[i][0] != I1 || I[i][1] != I2)
			{
				for (i=0; (i<M ? I[i][0] != I1 || I[i][1] != I2 : false); i++) ;
				if (i==M) 
				{
					M++;
					I[i][0] = I1;
					I[i][1] = I2;
				}
				if (i==IsoI) return;
			}
		}
		return;
	}
	printf("IsoI=%d,Z=%d\n", IsoI, Z[IsoI]);
	IsoTab *Iso = molecule->getIso();
	I1 = Iso->mNumIso1[Z[IsoI]];
	I2 = Iso->mNumIso2[Z[IsoI]];
	delete Iso;*/
}

ElState* TermTable::getElState()
{
	return State;
}

void TermTable::setData(double ****nData, int numComp, int numIso, int maxv, int maxJ, int *CompZ, int nMixC, double *****newMixC)
{
	int *Z, i;
	Z = new int[numIso];
	for (i=0; i < numIso; i++) Z[i] = i;
	tDat->setIsoZ(Z);
	tDat->setData(nData, numComp, numIso, maxv, maxJ, CompZ, nMixC, newMixC);
	setNewCreated();
	Changed();
}

void TermTable::setMolecule(Molecule *mol, ElState *nState)
{
	Molecule *aMol = molecule;
	TableWindow::setMolecule(mol);
	State = nState;
	if (mol != 0)
	{
        if (aMol != mol)
        {
            IsoTab *IsoT = tDat->getIso(), *NIsoT = mol->getIso();
            int *NIsoZ = 0, NI = tDat->getNumIso(), *AIsoZ = tDat->getIsoZ();
            if (0 != IsoT && 0 != NIsoT && 0 != AIsoZ)
            {
                int n, a;
                for (NIsoZ = new int[NI], a=0; a < NI; ++a)
                {
                    for (NIsoZ[a] = -1, n=0; n < NIsoT->numIso; ++n) if (IsoT->mNumIso1[AIsoZ[a]] == NIsoT->mNumIso1[n] && IsoT->mNumIso2[AIsoZ[a]] == NIsoT->mNumIso2[n])
                        NIsoZ[a] = n;
                    if (NIsoZ[a] == -1)
                    {
                        delete[] NIsoZ;
                        NIsoZ = 0;
                        break;
                    }
                }
            }
            tDat->setIso(NIsoT);
            if  (0 != NIsoZ) tDat->setIsoZ(NIsoZ);
            else if (tDat->getData() != 0) AssignIso();
        }
        else tDat->setIso(mol->getIso());
	}
}

void TermTable::setName(QString name)
{
    QString oldname = getName();
	TableWindow::setName(name);
	if (NPerturbations == 0) return;
	int n, m, np, NS, NT, s, t, p, NP;
	TermTable *P[MaxTermTables];
	QString PN[MaxTermTables], N;
	ElState *S;
	Perturbation *Pt;
	for (n = np = 0; n < NPerturbations; n++)
	{
		for (m=0; (m < np ? 
			(Perturbations[n].Perturber != 0 ? Perturbations[n].Perturber != P[m]
											 : Perturbations[n].PName != PN[m]) : false); m++) ;
		if (m == np && m < MaxTermTables)
		{
			PN[np++] = Perturbations[n].PName;
			if (Perturbations[n].Perturber == 0)
			{
				if (molecule != 0) for (s=0, NS = molecule->getNumStates(); s < NS && N != PN[m]; s++)
				{
					NT = (S = molecule->getStateP(s))->getNumTermTables();
					for (t=0; (t < NT ? (N = S->getTermTableName(t)) != PN[m] : false); t++) ;
					if (N == PN[m]) Perturbations[n].Perturber = S->getTermTable(t);
				}
			}
			P[m] = Perturbations[n].Perturber;
			if (P[m] != 0) for (p=0, NP = P[m]->NPerturbations; p < NP; p++) 
				if ((Pt = P[m]->Perturbations + p)->PName == oldname)
			{
				Pt->PName = name;
				Pt->Perturber = this;
			}
		}
	}
}

void TermTable::setPerturbations(int NPert, Perturbation* NewPerturbations)
{
	if (NPerturbations > 0) delete[] Perturbations;
	Perturbations = NewPerturbations;
	NPerturbations = NPert;
}

bool TermTable::writeData(QString NFilename)
{
	int NI = tDat->getNumIso(), NC = tDat->getNumComp(), Nv = tDat->getMaxv() + 1, nMixC = tDat->getnumStates();
	int NJ = tDat->getMaxJ() + 1, c, v, J, I, M1, M2, n, *CompZ = tDat->getCompZ();
	double ****Data = tDat->getData(), *****MixC = tDat->getMixCoeff();
	QFile Datei(NFilename);
	QString SM1, SM2;
	if (!write(&Datei)) return false; 
	QTextStream S(&Datei);
	S << "Source of data: " << getSource() << "\n";
	S << "Name: " << getName() << "\n";
	S << "Iso 1: | Iso 2: | component: | v: | J: | term energy:";
	for (n=0; n < nMixC; n++) S << " | Mix C S" + QString::number(n) + ":";
	S << '\n';
	for (I=0; I < NI; I++) 
	{
		tDat->GetIsoZ(I, M1, M2);
		SM1 = QString::number(M1);
		SM2 = QString::number(M2);
		for (c=0; c < NC; c++) for (v=0; v < Nv; v++) for (J=0; J < NJ; J++) 
					if (Data[c][I][v][J] != 0.0) 
		{
			S << SM1 << " | " << SM2 << " | " << QString::number(CompZ[c]) << " | ";
			S << QString::number(v) << " | " << QString::number(J) << " | ";
			S << QString::number(Data[c][I][v][J], 'f', 4);
			for (n=0; n < nMixC; n++) S << " | " << QString::number(MixC[c][I][v][J][n], 'f', 5);
			S << '\n';
		}
	}
	if (NPerturbations > 0)
	{
		S << "\nPerturbations: N=" + QString::number(NPerturbations) + "\n";
		S << "Iso 1: | Iso 2: | component: | v: | J: | Perturber: | Pcomp: | Pv:\n";
		for (n=0; n < NPerturbations; n++)
		{
			tDat->GetIsoZ(Perturbations[n].Iso, M1, M2);
			SM1 = QString::number(M1);
			SM2 = QString::number(M2);
			if (Perturbations[n].PName.isEmpty()) Perturbations[n].PName = Perturbations[n].Perturber->getName();
			S << SM1 + " | " + SM2 + " | " + QString::number(Perturbations[n].Comp) + " | "
			  << QString::number(Perturbations[n].v) + " | " + QString::number(Perturbations[n].J)
			  << " | " + Perturbations[n].PName + " | " + QString::number(Perturbations[n].PComp)
			  << " | " + QString::number(Perturbations[n].Pv) + "\n";
		}
	}
	return true;
}

void TermTable::getThermPopulation(double **&Pop, int &nv, int &nJ, double T, int Iso, int Comp)
{
	int v, J, EA, Maxv = tDat->getMaxv(), MaxJ = tDat->getMaxJ();
	double UF = -100 * C_h * C_c / (C_kB * T), Sum = 0.0, ****Data = tDat->getData();
	Pop = Create(nv = Maxv + 1, nJ = MaxJ + 1);
	for (J=0, EA = 1; J < nJ; J++, EA += 2) for (v=0; v < nv; v++) 
			Sum += (Pop[v][J] = (Data[Comp][Iso][v][J] != 0.0 ? EA * exp(Data[Comp][Iso][v][J] * UF) 
															  : 0.0));
	for (J=0, UF = 1.0 / Sum; J < nJ; J++) for (v=0; v < nv; v++) Pop[v][J] *= UF;
}

void TermTable::getSelE(int*& Js, double*& E, int& N)
{
	int *Rows = 0;
	table->getSelectedRows(Rows, N);
	Js = new int[N];
	E = new double[N];
	tDat->getJE(Rows, N, Js, E);
	delete[] Rows;
}

void TermTable::getViewnE(int*& Js, double*& E, int& N)
{
	int NR = tDat->rowCount(), n;
	bool *RV = new bool[NR];
	getViewnRows(RV);
	for (n=N=0; n < NR; n++) if (RV[n]) N++;
	int *Rows = new int[N];
	for (n=N=0; n < NR; n++) if (RV[n]) Rows[N++] = n;
	Js = new int[N];
	E = new double[N];
	tDat->getJE(Rows, N, Js, E);
	delete[] Rows;
}

void TermTable::setViewnLevels(MDIChild* Viewer, int C, int I, int v, int* J, int N)
{
	int *Rows = new int[N];
	tDat->getRows(C, I, v, J, N, Rows);
	setViewnRows(Viewer, N, Rows);
}

void TermTable::getMaxAndAverageDeviation(TermTable* Other, int comp, int iso, int vMax, 
										  int JMax, double& MaxDev, double& AvDev)
{
	int c, I, v, J, Nc = (comp >= 0 ? 1 : tDat->getNumComp()), NcO = Other->getNumComp();
	int NI = (iso >= 0 ? 1 : tDat->getNumIso()), NIO = Other->getNumIso();
	int cA[Nc], IA[NI], cAO[Nc], IAO[Nc], *compZ = tDat->getCompZ(), *compZO = Other->getCompZ();
	int *isoZ = tDat->getIsoZ(), *isoZO = Other->getIsoZ(), n, m; 
	double ****Data = tDat->getData(), ****DataO = Other->getData(), curDev;
	if (comp >= 0)
	{
		for (cA[0] = 0; compZ[cA[0]] != comp; cA[0]++) ;
		for (cAO[0] = 0; compZO[cAO[0]] != comp; cAO[0]++) ;
	}
	else 
	{
		for (c=m=0; c < Nc; c++)
		{
			for (n=0; (n < NcO ? compZ[c] != compZO[n] : false); n++) ;
			if (n < NcO)
			{
				cA[m] = c;
				cAO[m++] = n;
			}
		}
		Nc = m;
	}
	if (iso >= 0)
	{
		for (IA[0] = 0; isoZ[IA[0]] != iso; IA[0]++) ;
		for (IAO[0] = 0; isoZO[IAO[0]] != iso; IAO[0]++) ;
	}
	else
	{
		for (I=m=0; I < NI; I++)
		{
			for (n=0; (n < NIO ? isoZ[I] != isoZO[I] : false); n++) ;
			if (n < NIO)
			{
				IA[m] = I;
				IAO[m++] = n;
			}
		}
		NI = m;
	}
	for (c=m=0; c < Nc; c++) for (I=0; I < NI; I++) for (v=0; v <= vMax; v++)
		for (J=0; J <= JMax; J++) 
			if (Data[cA[c]][IA[I]][v][J] != 0.0 && DataO[cAO[c]][IAO[I]][v][J] != 0.0)
	{
		curDev = fabs(Data[cA[c]][IA[I]][v][J] - DataO[cAO[c]][IAO[I]][v][J]);
		if (curDev > MaxDev) MaxDev = curDev;
		AvDev += curDev;
		m++;
	}
	AvDev /= m;
}
