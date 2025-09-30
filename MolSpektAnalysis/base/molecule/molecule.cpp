//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "molecule.h"
#include "termtable.h"
#include "duntable.h"
#include "linetable.h"
#include "potential.h"
#include "utils.h"
#include "fcftab.h"
#include "fitdata.h"
#include "isotab.h"
#include "TableItem.h"

#include <math.h>

#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QFileDialog>
#include <QGridLayout>


Molecule::Molecule(MainWindow *mw) : MDIChild(MDIChild::MolData, mw, "Molecules (*.mol)", ".mol")
{
	//printf("Molecule::Molecule\n");
	MW = mw;
	UNI = true;
	numStates = numTransitions = 0;
	atom1 = atom2 = 0;
	setMinimumSize(400, 490);
	QGridLayout *Layout = new QGridLayout(this);
	Layout->addWidget(new QLabel("Atom1:", this), 0, 0);
	Layout->addWidget(BAtom1 = new QComboBox(this), 0, 1);
	BAtom1->setEditable(false);
	connect(BAtom1, SIGNAL(currentIndexChanged(int)), this, SLOT(updateNI()));
	Layout->addWidget(new QLabel("Atom2:", this), 0, 2);
	Layout->addWidget(BAtom2 = new QComboBox(this), 0, 3);
	BAtom2->setEditable(false);
	connect(BAtom2, SIGNAL(currentIndexChanged(int)), this, SLOT(updateNI()));
	Layout->addWidget(new QLabel("Ref. iso.:", this), 1, 0);
	Layout->addWidget(RefIso = new QComboBox(this), 1, 1);
	RefIso->setEditable(false);
	connect(RefIso, SIGNAL(currentIndexChanged(int)), this, SLOT(Changed()));
	Layout->addWidget(Name = new QLabel("name of molecule:", this), 1, 2, 1, 2);
	Layout->setRowMinimumHeight(2, 20);
	Layout->addWidget(new QLabel("Known electronic states:", this), 3, 0, 1, 4);
	Layout->addWidget(TStates = new QTableWidget(MaxStates, 5, this), 4, 0, 1, 4);
	TStates->setHorizontalHeaderLabels(QStringList() << "name" 
			<< "term energies" << "dunham coeff." << "potential" << "fit data");
	connect(TStates, SIGNAL(currentCellChanged(int, int, int, int)), 
			this, SLOT(TStatesClicked(int, int)));
	connect(TStates, SIGNAL(itemChanged(QTableWidgetItem*)), 
			this, SLOT(TStatesChanged(QTableWidgetItem*)));
	Layout->setRowMinimumHeight(5, 20);
	Layout->addWidget(new QLabel("Used transitions:", this), 6, 0, 1, 4);
	Layout->addWidget(TTransitions = new QTableWidget(MaxTransitions, 4, this), 7, 0, 1, 4);
	TTransitions->setHorizontalHeaderLabels(QStringList() << "lower state" << "upper state"
			<< "lines" << "FCF");
	connect(TTransitions, SIGNAL(currentCellChanged(int, int, int, int)), 
			this, SLOT(TTransitionsClicked(int, int)));
	setWindowTitle("New molecule");
	//printf("Vor update Atoms\n");
	updateAtoms();
	Saved();
}


Molecule::~Molecule()
{

}

bool Molecule::readData(QString FileName)
{
	if (MW == 0) return false;
	int n, m, k, L=0, SI = 0;
	int ST = 0;
	for (n=0; n < numStates; n++) if (States[n] != 0) States[n]->close();
	for (n=0; n < numTransitions; n++) if (Transitions[n] != 0) Transitions[n]->close();
	numStates = numTransitions = 0;
	QString Buffer, name, TName, TFile, TSource, TState1, TState2;
	QFile Datei(FileName);
	TermTable *term;
	DunTable *Dun;
	Potential *Pot;
	FitData *FD;
	LineTable *Line;
	FCFTab *FCF;
	if (!read(&Datei)) return false;
	QTextStream S(&Datei);
	QStringList SList;
	ComboBox *Box;
	ElState *St;
	Atom *atom;
	//TableItem *DebT = 0;
	while (!S.atEnd())
	{
		//printf("read line %d, i=%d\n", line++, numStates);
		//if (DebT != 0) 
			//printf("StateText = %s\n", DebT->text().ascii());
		Buffer = S.readLine();
		if (Buffer.indexOf("Table of states", 0, Qt::CaseInsensitive) != -1) 
		{
			ST = 1;
			continue;
		}
		if (Buffer.indexOf("Table of transitions", 0, Qt::CaseInsensitive) != -1) 
		{
			ST = 2;
			continue;
		}
		if (Buffer.indexOf("List of term energy tables", 0, Qt::CaseInsensitive) != -1) 
		{
			ST = 3;
			continue;
		}
		if (Buffer.indexOf("List of dunham coefficient sets", 0, Qt::CaseInsensitive) != -1) 
		{
			ST = 4;
			continue;
		}
		if (Buffer.indexOf("List of potentials", 0, Qt::CaseInsensitive) != -1) 
		{
			ST = 5;
			continue;
		}
		if (Buffer.indexOf("List of line tables", 0, Qt::CaseInsensitive) != -1) 
		{
			ST = 6;
			continue;
		}
		if (Buffer.indexOf("List of FCF tables", 0, Qt::CaseInsensitive) != -1)
		{
			ST = 7;
			continue;
		}
		if (Buffer.indexOf("List of fit datasets", 0, Qt::CaseInsensitive) != -1)
		{
			ST = 8;
			continue;
		}
		switch (ST)
		{
			case 0:
				if (Buffer.indexOf("atom 1", 0, Qt::CaseInsensitive) != -1)
				{
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
                    Buffer = Buffer.right(Buffer.length() - n);
                    atom = MW->getAtom(getAbsolutePath(Buffer, FileName));
					if (atom != 0)
					{
                        //printf("atom1:\n");
                        BAtom1->addItem(atom->getName());
						BAtom1->setCurrentIndex(BAtom1->count() - 1);
					}
					//printf("atom1=%d\n", atom1);
				}
				else if (Buffer.indexOf("atom 2", 0, Qt::CaseInsensitive) != -1)
				{
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
                    Buffer = Buffer.right(Buffer.length() - n);
                    atom = MW->getAtom(getAbsolutePath(Buffer, FileName));
					if (atom != 0)
					{
						//printf("atom2:\n");
						BAtom2->addItem(atom->getName());
						BAtom2->setCurrentIndex(BAtom2->count() - 1);
					}
					//printf("atom2=%d\n", atom2);
				}
                else if (Buffer.indexOf("reference isotopo", 0, Qt::CaseInsensitive) != -1)
				{
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					RIN = Buffer.right(Buffer.length() - n);
					//printf("RIN=%s\n", RIN.ascii());
					for (n=0; (n < RefIso->count() ? RefIso->itemText(n) != RIN : false); n++) ;
					RefIso->setCurrentIndex(n);
				}
				else if (Buffer.indexOf("name of molecule", 0, Qt::CaseInsensitive) != -1) 
				{
					for (n = Buffer.indexOf(":") + 1; Buffer[n].isSpace(); n++) ;
					name = Buffer.right(Buffer.length() - n);
					//printf("name=%s\n", name.ascii());
					Name->setText("name of molecule: " + name);
					setName(name);
				}
				break;
			case 1:
				if (Buffer.isEmpty()) break;
                if (Buffer.indexOf("Name:") != -1)
				{
					if (Buffer.indexOf("lambda", 0, Qt::CaseInsensitive) != -1) L=9;
					else L=4;
					break;
				}
				SList = Buffer.split('|');
				if (SList.count() < L) return false;
                else if (SList.count() >= 11) L = 11;
				for (k=0; k<L; k++)
				{
					Buffer = SList[k];
					for (n=0; n < Buffer.length() && Buffer[n].isSpace(); n++) ;
					for (m = Buffer.length() - 1; m >= n && Buffer[m].isSpace(); m--) ;
					//printf("k=%d\n", k);
					if (n <= m)
					{
						Buffer = Buffer.mid(n, m - n + 1);
						//printf("numStates=%d, k=%d, Buffer=%s!\n", numStates, k, Buffer.ascii());
						switch (k)
						{
							case 0:
								for (SI=0; (SI < numStates ? States[SI]->getName() != Buffer : false);
															 SI++) ;
                                if (SI == numStates) addState(Buffer, true);
								break;
							case 1:
								if (L>=9)
								{
									States[SI]->setLambda(Buffer.toInt());
									break;
								}
								/* Falls through. */
							case 6:
                                if (!Buffer.isEmpty()) if (!States[SI]->setMainTermTable(Buffer = getAbsolutePath(Buffer, FileName)))
                                    if ((term = MW->getTermTable(Buffer, this, States[SI])) != 0)
								{
									//printf("States[%d].termTable=%d\n", SI, term);
									States[SI]->setMainTermTable(term);
									if ((Box = States[SI]->getTermBox()) == 0)
									{
										States[SI]->setTermBox(Box = new ComboBox);
										Box->setEditable(false);
										TStates->setCellWidget(SI, 1, Box);
									}
								}
								break;
							case 2:
								if (L>=9)
								{
									if (Buffer.indexOf('.') >= 0) States[SI]->setS(Buffer.toFloat());
									else States[SI]->setS(Buffer.toInt());
									break;
								}
								/* Falls through. */
							case 7:
                                if (!Buffer.isEmpty()) if (!States[SI]->setMainDunTable(Buffer = getAbsolutePath(Buffer, FileName)))
									if ((Dun = MW->getDunTable(Buffer, this)) != 0)
								{
									States[SI]->setMainDunTable(Dun);
									if ((Box = States[SI]->getDunBox()) == 0)
									{
										States[SI]->setDunBox(Box = new ComboBox);
										Box->setEditable(false);
										TStates->setCellWidget(SI, 2, Box);
									}
								}
								break;
							case 3:
								if (L>=9)
								{
									if (Buffer == "g") States[SI]->setParity(1);
									else if (Buffer == "u") States[SI]->setParity(2);
									break;
								}
								/* Falls through. */
							case 8:
                                if (!Buffer.isEmpty()) if (!States[SI]->setMainPotential(Buffer = getAbsolutePath(Buffer, FileName)))
									if ((Pot = MW->getPotential(Buffer, this)) != 0)
								{
									States[SI]->setMainPotential(Pot);
									if ((Box = States[SI]->getPotBox()) == 0)
									{
										States[SI]->setPotBox(Box = new ComboBox(0));
										Box->setEditable(false);
										TStates->setCellWidget(SI, 3, Box);
									}
								}
								break;
							case 4:
								if (Buffer == "+") States[SI]->setSymmetry(1);
								else if (Buffer == "-") States[SI]->setSymmetry(-1);
								break;
							case 5:
								States[SI]->setOmega(Buffer.toInt());
								break;
							case 9:
                                if (!Buffer.isEmpty()) if (!States[SI]->setMainFitData(Buffer = getAbsolutePath(Buffer, FileName)))
									if ((FD = MW->getFitData(Buffer, this)) != 0)
								{
									States[SI]->setMainFitData(FD);
									if ((Box = States[SI]->getFitDataBox()) == 0)
									{
										States[SI]->setFitDataBox(Box = new ComboBox(0));
										Box->setEditable(false);
										TStates->setCellWidget(SI, 4, Box);
									}
								}
								break;
                            case 10:
                                States[SI]->setBe(Buffer.toDouble());
                                break;
						}
						if (SI == MaxStates) break;
					}
				}
				break;
			case 2:
				//printf("Transitions[0]=%d\n", Transitions[0]);
                if (Buffer.indexOf("lines:") != -1 || Buffer.isEmpty()) break;
				SList = Buffer.split('|');
				if (SList.count() < 4) break;
				for (k=0; k < SList.count(); k++)
				{
					Buffer = SList[k];
					for (n=0; n < Buffer.length() && Buffer[n].isSpace(); n++) ;
					if (n == Buffer.length()) SList[k] = "";
					else 
					{
						for (m = Buffer.length() - 1; Buffer[m].isSpace(); m--) ;
						SList[k] = Buffer.mid(n, m - n + 1);
					}
				}
				//printf("Vor Zustandssuche, j=%d\n", numTransitions);
				for (k=0; (k < numTransitions ? ((St = Transitions[k]->getUpperState()) != 0 ? 
								 St->getName() : "") != SList[1] 
						|| ((St = Transitions[k]->getLowerState()) != 0 ? St->getName() : "") 
										 != SList[0] : false); k++) ;
				//printf("Nach Zustandssuche, k=%d\n", k);
				//printf("UState=%d, LState=%d\n", 
					//	   Transitions[0]->getUpperState(), Transitions[0]->getLowerState());
				//printf("Transitions[0]=%d\n", Transitions[0]);
				if (k == numTransitions) 
				{
					//printf("Erzeuge Transition\n");
					if (numTransitions == MaxTransitions) break;
					Transitions[numTransitions] = MW->CreateTransition();
					Transitions[numTransitions]->setMW(MW, this);
					Transitions[numTransitions]->setTransNum(numTransitions);
					if ((Box = Transitions[numTransitions]->getLSB()) == 0)
					{
						Transitions[numTransitions]->setLSB(Box = new ComboBox(0));
						Box->setEditable(false);
						TTransitions->setCellWidget(numTransitions, 0, Box);
					}
					if ((Box = Transitions[numTransitions]->getUSB()) == 0)
					{
						Transitions[numTransitions]->setUSB(Box = new ComboBox(0));
						Box->setEditable(false);
						TTransitions->setCellWidget(numTransitions, 1, Box);
					}
					//printf("numStates=%d\n", numStates);
					for (n=0; n < numStates; n++)
					{
						if (States[n]->getName() == SList[1]) 
						{
							Transitions[k]->setUpperState(States[n]);
							//printf("UpperState=%d\n", n);
						}
						else if (States[n]->getName() == SList[0]) 
						{
							Transitions[k]->setLowerState(States[n]);
							//printf("LowerState=%d\n", n);
						}
					}
					if (Transitions[numTransitions++]->getLowerState() == 0 && numStates < MaxStates
					    		&& !SList[0].isEmpty())
                        Transitions[k]->setLowerState(((n = addState(SList[0], true)) != -1 ?
								States[n] : 0));
					if (Transitions[k]->getUpperState() == 0 && numStates < MaxStates
					   			&& !SList[1].isEmpty())
                        Transitions[k]->setUpperState(((n = addState(SList[1], true)) != -1 ?
								States[n] : 0));
                    ElState *lowerState = Transitions[k]->getLowerState(), *upperState = Transitions[k]->getUpperState();
                    updateTransitionName(k, (lowerState != 0 ? lowerState->getStateNum() : -1),
                                            (upperState != 0 ? upperState->getStateNum() : -1));
				}
                if (!SList[2].isEmpty()) if (!Transitions[k]->setMainLineTable(Buffer = getAbsolutePath(SList[2], FileName)))
                    if ((Line = MW->getLineTable(Buffer, this)) != 0)
				{
					//printf("m=%d, k=%d, N=%d, numStates=%d, numTransitions=%d\n", 
						//   m, k, N, numStates, numTransitions);
					Transitions[k]->setMainLineTable(Line);
					if ((Box = Transitions[k]->getLineBox()) == 0)
					{
						Transitions[k]->setLineBox(Box = new ComboBox(0));
						Box->setEditable(false);
						TTransitions->setCellWidget(k, 2, Box);
					}
				}
                if (!SList[3].isEmpty()) if (!Transitions[k]->setMainFCFTable(Buffer = getAbsolutePath(SList[3], FileName)))
                    if ((FCF = MW->getFCFTable(Buffer, this)) != 0)
				{
					Transitions[k]->setMainFCFTable(FCF);
					if ((Box = Transitions[k]->getFCFBox()) == 0)
					{
						Transitions[k]->setFCFBox(Box = new ComboBox(0));
						Box->setEditable(false);
						TTransitions->setCellWidget(k, 3, Box);
					}
				}
				if (SList.count() >= 5) Transitions[k]->setTransitionStrengthText(SList[4]);
				break;
			case 3:
			case 4:
			case 5:
			case 8:
				if (Buffer.indexOf("State:", 0, Qt::CaseInsensitive) != -1 || Buffer.isEmpty()) break;
				SList = Buffer.split('|');
				for (k=0; k<4; k++)
				{
					if (k < SList.count()) 
					{
						Buffer = SList[k];
						for (n=0; n < Buffer.length() && Buffer[n].isSpace(); n++) ;
						for (m = Buffer.length() - 1; m >= n && Buffer[m].isSpace(); m--) ;
						if (n <= m) Buffer = Buffer.mid(n, m - n + 1);
					}
					else Buffer = "";
					switch (k)
					{
						case 0:
							TState1 = Buffer;
							break;
						case 1:
							TName = Buffer;
							break;
						case 2:
							TFile = Buffer;
							break;
						case 3:
							TSource = Buffer;
							break;
					}
				}
				for (k=0; (k < numStates ? States[k]->getName() != TState1 : false); k++) ;
                if ((k == numStates ? addState(TState1, true) == -1 : false) || TFile.trimmed().isEmpty()) break;
				//printf("Anfang3: StateText = %s\n", DebT->text().ascii());
				switch (ST)
				{
					case 3:
						//printf("ST=3\n");
                        States[k]->addTermTable(TName, getAbsolutePath(TFile, FileName), TSource);
						TStatesClicked(k, 1);
						break;
					case 4:
						//printf("ST=4\n");
                        States[k]->addDunTable(TName, getAbsolutePath(TFile, FileName), TSource);
						TStatesClicked(k, 2);
						break;
					case 5:
						//printf("ST=5\n");
                        States[k]->addPotential(TName, getAbsolutePath(TFile, FileName), TSource);
						TStatesClicked(k, 3);
						break;
					case 8:
                        States[k]->addFitData(TName, getAbsolutePath(TFile, FileName), TSource);
						TStatesClicked(k, 4);
						break;
				}
				//printf("Ende: StateText = %s\n", DebT->text().ascii());
				break;
			case 6:
			case 7:
                if (Buffer.indexOf("state:") != -1 || Buffer.isEmpty()) break;
				SList = Buffer.split('|');
				for (k=0; k<5; k++)
				{
					if (k < SList.count()) 
					{
						Buffer = SList[k];
						for (n=0; n < Buffer.length() && Buffer[n].isSpace(); n++) ;
						if (n == Buffer.length()) Buffer = "";
						else
						{
							for (m = Buffer.length() - 1; m >= n && Buffer[m].isSpace(); m--) ;
							if (n <= m) Buffer = Buffer.mid(n, m - n + 1);
						}
					}
					else Buffer = "";
					switch (k)
					{
						case 0:
							TState1 = Buffer;
							break;
						case 1:
							TState2 = Buffer;
							break;
						case 2:
							TName = Buffer;
							break;
						case 3:
							TFile = Buffer;
							break;
						case 4:
							TSource = Buffer;
							break;
					}
				}
				if (TFile.trimmed().isEmpty()) break;
				//printf("TState1=%s, TState2=%s\n", TState1.ascii(), TState2.ascii());
				for (k=0; (k < numTransitions ? (Transitions[k]->getUpperState() != 0 ?
								 Transitions[k]->getUpperState()->getName() : "") != TState1 
							|| (Transitions[k]->getLowerState() != 0 ? 
						Transitions[k]->getLowerState()->getName() : "") != TState2 : false); k++) ;
				if (k == numTransitions) 
				{
					//printf("Erzeuge Transition: k=%d, numTransitions=%d\n", k, numTransitions);
					if (numTransitions == MaxTransitions) break;
					Transitions[numTransitions] = MW->CreateTransition();
					Transitions[numTransitions]->setMW(MW, this);
					Transitions[numTransitions]->setTransNum(numTransitions);
					if ((Box = Transitions[numTransitions]->getLSB()) == 0)
					{
						Transitions[numTransitions]->setLSB(Box = new ComboBox);
						Box->setEditable(false);
						TTransitions->setCellWidget(numTransitions, 0, Box);
					}
					if ((Box = Transitions[numTransitions]->getUSB()) == 0)
					{
						Transitions[numTransitions]->setUSB(Box = new ComboBox);
						Box->setEditable(false);
						TTransitions->setCellWidget(numTransitions, 1, Box);
					}
					for (n=0; n < numStates; n++)
					{
						if (States[n]->getName() == TState1) Transitions[k]->setUpperState(States[n]);
						else if (States[n]->getName() == TState2) 
							Transitions[k]->setLowerState(States[n]);
					}
					if (Transitions[numTransitions++]->getLowerState() == 0 && !TState2.isEmpty())
                        Transitions[k]->setLowerState(((n = addState(TState2, true)) != -1 ? States[n] : 0));
					if (Transitions[k]->getUpperState() == 0 && !TState1.isEmpty())
                        Transitions[k]->setUpperState(((n = addState(TState1, true)) != -1 ? States[n] : 0));
                    ElState *lowerState = Transitions[k]->getLowerState(), *upperState = Transitions[k]->getUpperState();
                    updateTransitionName(k, (lowerState != 0 ? lowerState->getStateNum() : -1),
                                            (upperState != 0 ? upperState->getStateNum() : -1));
					//printf("NumStates=%d, k=%d\n", numStates, k);
					//printf("UState=%d, LState=%d\n", 
						//   Transitions[0]->getUpperState(), Transitions[0]->getLowerState());
					//printf("Transitions[0]=%d\n", Transitions[0]);
				}
				if (ST == 6)
				{
                    Transitions[k]->addLineTable(TName, getAbsolutePath(TFile, FileName), TSource);
					TTransitionsClicked(k, 2);
				}
				else
				{
                    Transitions[k]->addFCFTable(TName, getAbsolutePath(TFile, FileName), TSource);
					TTransitionsClicked(k, 3);
				}
				//printf("Transitions[0]=%d\n", Transitions[0]);
				break;
		}
	}
	//printf("numStates=%d, numTransitions=%d\n", numStates, numTransitions);
    for (n=0; n < numStates; ++n) States[n]->setReading(false);
    Saved();
	return true;
}

Atom *Molecule::getAtom1()
{
	return atom1;
}

Atom *Molecule::getAtom2()
{
	return atom2;
}

void Molecule::setAtoms(Atom* A1, Atom* A2)
{
	int n, N = BAtom1->count();
	QString Name = A1->getName();
	for (n=0; (n<N  ? BAtom1->itemText(n) != Name : false); n++) ;
	if (n<N)
	{
		BAtom1->blockSignals(true);
		BAtom1->setCurrentIndex(n);
		BAtom1->blockSignals(false);
	}
	for (n=0, Name = A2->getName(); (n<N ? BAtom2->itemText(n) != Name : false); n++) ;
	if (n<N) BAtom2->setCurrentIndex(n);
}

QString Molecule::getState(int Index)
{
	if (Index < 0 || Index >= numStates)
	{
		printf("Error: the state with Index = %d is not existing!\n", Index);
		return 0;
	}
	return States[Index]->getName();
}

ElState *Molecule::getState(QString Name)
{
	int i;
	for (i=0; i < numStates; i++) if (States[i]->getName() == Name) return States[i];
	printf("Molecule::getState: a state with the name %s does not exist!\n", Name.toLatin1().data());
	return 0;
}

ElState *Molecule::getStateP(int Index)
{
	if (Index < 0 || Index >= numStates)
	{
		printf("Error: the state with Index = %d is not existing!\n", Index);
		return 0;
	}
	return States[Index];
}

int Molecule::getNumTransitions()
{
	return numTransitions;
}

Transition *Molecule::getTransitionP(int i)
{
	if (i < 0 || i >= numTransitions)
	{
		printf("Error: the transition with Index = %d is not existing!\n", i);
		return 0;
	}
	return Transitions[i];
}

double Molecule::getTransitionStrength(ElState* State1, ElState* State2)
{
	if (State1 == 0 || State2 == 0) return 0.0;
	int n;
	for (n=0; n < numTransitions; n++) 
		if ((Transitions[n]->getLowerState() == State1 && Transitions[n]->getUpperState() == State2)
			|| (Transitions[n]->getUpperState() == State1 && Transitions[n]->getLowerState() == State2)) 
		return Transitions[n]->getTransitionStrength();
	return 0.0;
}

int Molecule::addState(QString Name, const bool i_reading)
{
	if (numStates == MaxStates)
	{
		QMessageBox::warning(this, "QT4MolSpektAn", 
						"Error: the maximum amount of states has been reached for this molecule!");
		return -1;
	}
	int n = numStates;
	TableItem *TI = new TableItem(Name);
	States[n] = MW->CreateElState();
	States[n]->setMW(MW, this);
	States[n]->setStateNum(n);
	numStates++;
	TStates->setCellWidget(n, 0, TI);
	//printf("TI=%s\n", TI->text().ascii());
	connect(TI, SIGNAL(RightClicked(QPoint)), States[n], SLOT(showShowMenu(QPoint)));
	connect(TI, SIGNAL(textEdited(QString)), States[n], SLOT(setName(QString)));
	connect(States[n], SIGNAL(nameChanged(QString)), TI, SLOT(setText(QString)));
	States[n]->setName(Name);
    States[n]->setReading(i_reading);
	return n;
}

int Molecule::addTransition(ElState *LS, ElState *US, LineTable *L, FCFTab *F)
{
	if (numTransitions == MaxTransitions)
	{
		QMessageBox::warning(this, "MolSpektAnalysis",
					"Error: the maximum amount of transitions has been reached for this molecule!");
		return -1;
	}
	int c, n;
	for (c=0, n = numTransitions; c<4; c++) TTransitionsClicked(n, c);
	Transitions[n]->setLowerState(LS);
	Transitions[n]->setUpperState(US);
	Transitions[n]->addLineTable(L);
	Transitions[n]->addFCFTable(F);
	return n;
}

void Molecule::addLineTable(LineTable *L, ElState *LS, ElState *US)
{
	int t;
	for (t=0; (t < numTransitions ? Transitions[t]->getLowerState() != LS 
			|| Transitions[t]->getUpperState() != US : false); t++) ;
	if (t < numTransitions) Transitions[t]->addLineTable(L);
	else addTransition(LS, US, L);
}

void Molecule::addFCFTable(FCFTab* Tab, ElState* LState, ElState* UState)
{
	int t;
	for (t=0; (t < numTransitions ? Transitions[t]->getLowerState() != LState
		|| Transitions[t]->getUpperState() != UState : false); t++) ;
	if (t < numTransitions) Transitions[t]->addFCFTable(Tab);
	else addTransition(LState, UState, 0, Tab);
}

void Molecule::addTerm(int i, TermTable *TermTab)
{
	printf("Molecule::addTerm\n");
	if (i < 0 || i >= numStates)
	{
		printf("Molecule::addTerm: the state with SIndex = %d is not existing!\n", i);
		return;
	}
	States[i]->addTermTable(TermTab);
	TStatesClicked(i, 1);
	Changed();
}

TermTable *Molecule::getTerm(int SIndex)
{
	if (SIndex >= 0 && SIndex < numStates) return States[SIndex]->getTermTable();
	printf("Molecule::getTerm: the state with SIndex = %d is not existing!\n", SIndex);
	return 0;
}

TermTable *Molecule::getTermTable(QString Name)
{
    int n, m, NT;
    for (n=0; n < numStates; ++n) for (m=0, NT = States[n]->getNumTermTables(); m < NT; ++m) if (States[n]->getTermTableName(m) == Name)
        return States[n]->getTermTable(m);
    return 0;
}

ElState *Molecule::getStateOfTermTable(QString Name)
{
    int n, m, NT;
    for (n=0; n < numStates; ++n) for (m=0, NT = States[n]->getNumTermTables(); m < NT; ++m) if (States[n]->getTermTableName(m) == Name)
        return States[n];
    return 0;
}

void Molecule::getTermData(int StateNum, int &numComp, int &numIso, int &numv, int &numJ, 
						   double ****&Data)
{
	printf("Molecule::getTermData: StateNum=%d\n", StateNum);
	TermTable *T = getTerm(StateNum);
	if (T == 0)
	{
		Data = 0;
		numComp = numIso = numv = numJ = 0;
	}
	else
	{
		Data = T->getData();
		numComp = T->getNumComp();
		numIso = T->getNumIso();
		numv = T->getMaxv() + 1;
		numJ = T->getMaxJ() + 1;
	}
}

void Molecule::addDunK(int SIndex, DunTable *DunhamK)
{
	if (SIndex >= 0 && SIndex < numStates) 
	{
		States[SIndex]->addDunTable(DunhamK);
		Changed();
	}
	else 
		printf("Molecule::addDunK: the state with SIndex = %d is not existing!\n", SIndex);
}

DunTable *Molecule::getDunK(int SIndex)
{
	if (SIndex >= 0 && SIndex < numStates) return States[SIndex]->getDunTable();
	printf("Molecule::getDunK: the state with SIndex = %d is not existing!\n", SIndex);
	return 0;
}

void Molecule::addPot(int SIndex, Potential *Pot)
{
	//printf("Molecule::addPot(%d, Potential\n", SIndex);
	if (SIndex >= 0 && SIndex < numStates) 
	{
		States[SIndex]->addPotential(Pot);
		Changed();
	}
	else printf("Molecule::addPot: the state with SIndex = %d is not existing!\n", SIndex);
}

Potential *Molecule::getPot(int SIndex)
{
	if (SIndex >= 0 && SIndex < numStates) return States[SIndex]->getPotential();
	printf("Molecule::getPot: the state with SIndex = %d is not existing!\n", SIndex);
	return 0;
}

QString Molecule::getPotName(int SIndex)
{
	if (SIndex >= 0 && SIndex < numStates) return States[SIndex]->getPotentialName();
	printf("Molecule::getPot: the state with SIndex = %d is not existing!\n", SIndex);
	return 0;
}

Potential *Molecule::getPot(QString SName)
{
	ElState *S = getState(SName);
	return (S!=0 ? S->getPotential() : 0);
}

void Molecule::addFitData(int SIndex, FitData* FD)
{
	if (SIndex >= 0 && SIndex < numStates)
	{
		States[SIndex]->addFitData(FD);
		Changed();
	}
	else printf("Molecule::addFitData: the state with SIndex = %d is not existing!\n", SIndex);
}

FitData* Molecule::getFitData(int SIndex)
{
	if (SIndex >= 0 && SIndex < numStates) return States[SIndex]->getFitData();
	printf("Molecule::getFitData: the state with SIndex = %d is not existing!\n", SIndex);
	return 0;
}

FitData* Molecule::getFitData(QString Name)
{
    int s, f, F;
    for (s=0; s < numStates; ++s) for (f=0, F = States[s]->getNumFitDataSets(); f<F; ++f) if (States[s]->getFitDataName(f) == Name)
        return States[s]->getFitData(f);
    return 0;
}

void Molecule::addLines(int SIndex1, int SIndex2, LineTable *Lines)
{
	int i;
	if ((i = getTransitionIndex(SIndex1, SIndex2)) != -1) 
	{
		Transitions[i]->addLineTable(Lines);
		Changed();
	}
}

void Molecule::addFCF(int SIndex1, int SIndex2, FCFTab* FCF)
{
	int i;
	if ((i = getTransitionIndex(SIndex1, SIndex2)) != -1)
	{
		Transitions[i]->addFCFTable(FCF);
		Changed();
	}
}

int Molecule::getTransitionIndex(int SIndex1, int SIndex2)
{
	if (SIndex1 < 0 || SIndex1 >= numStates)
	{
		printf("Molecule::getTransitionIndex: the state with SIndex1 = %d is not existing!\n",
			   SIndex1);
		return -1;
	}
	if (SIndex2 < 0 || SIndex2 >= numStates)
	{
		printf("Molecule::getTransitionIndex: the state with SIndex2 = %d is not existing!\n",
			   SIndex2);
		return -1;
	}
	int i = 0;
	while (Transitions[i]->getLowerState() != States[SIndex1] 
			  || Transitions[i]->getUpperState() != States[SIndex2])
		if (++i == numTransitions)
	{
		if (numTransitions == MaxTransitions)
		{
			QMessageBox::warning(this, "QT4MolSpektAn", 
				"Error: the maximum amount of transitions has been reached for this molecule!");
			return -1;
		}
		numTransitions++;
		Transitions[i]->setLowerState(States[SIndex1]);
		TStates->setItem(i, 0, new QTableWidgetItem(States[SIndex1]->getName()));
		Transitions[i]->setUpperState(States[SIndex2]);
		TStates->setItem(i, 1, new QTableWidgetItem(States[SIndex2]->getName()));
		break;
	}
	return i;
}

LineTable *Molecule::getLines(int SIndex1, int SIndex2)
{
	if (SIndex1 < 0 || SIndex1 >= numStates)
	{
		printf("Molecule::getLines: the state with SIndex1 = %d is not existing!\n", SIndex1);
		return 0;
	}
	if (SIndex2 < 0 || SIndex2 >= numStates)
	{
		printf("Molecule::getLines: the state with SIndex2 = %d is not existing!\n", SIndex2);
		return 0;
	}
	int i;
	for (i=0; i < numTransitions; i++) 
		if (Transitions[i]->getLowerState() == States[SIndex1] 
				  && Transitions[i]->getUpperState() == States[SIndex2]) 
			return Transitions[i]->getLineTable();
	printf("Molecule::getLines: the transition with SIndex1 = %d and SIndex2 = %d does not exist!\n",
		   SIndex1, SIndex2);
	return 0;
}

FCFTab* Molecule::getFCF(int SIndex1, int SIndex2)
{
	if (SIndex1 < 0 || SIndex2 >= numStates)
	{
		printf("Molecule::getFCF: the state with SIndex1 = %d is not existing!\n", SIndex1);
		return 0;
	}
	if (SIndex2 < 0 || SIndex2 >= numStates)
	{
		printf("Molecule::getLines: the state with SIndex2 = %d is not existing!\n", SIndex2);
		return 0;
	}
	int i;
	for (i=0; i < numTransitions; i++)
		if (Transitions[i]->getLowerState() == States[SIndex1]
				&& Transitions[i]->getUpperState() == States[SIndex2])
			return Transitions[i]->getFCFTable();
	printf("Molecule::getLines: the transition with SIndex1 = %d and SIndex2 = %d does not exist!\n",
		   SIndex1, SIndex2);
	return 0;
}

LineTable *Molecule::getLineTable(ElState *lState, ElState *uState)
{
	//printf("Molecule::getLineTable\n");
	int n, c;
	LineTable *L;
	for (n=0; n < numTransitions; n++) 
		if (Transitions[n]->getLowerState() == lState && Transitions[n]->getUpperState() == uState) 
	{
		//printf("n=%d, lState=%d, uState=%d\n", n, lState, uState);
		if ((L = Transitions[n]->getLineTable()) == 0) 
			if ((L = (MW != 0 ? MW->CreateLineTable() : 0)) != 0)
		{
			Transitions[n]->addLineTable(L);
			Changed();
			//printf("New LineTable\n");
		}
		if (!L->isVisible()) L->show();
		//printf("Ende2 getLineTable\n");
		return L;
	}
	if (numTransitions == MaxTransitions) return 0;
	if ((L = (MW != 0 ? MW->CreateLineTable() : 0)) == 0) return L;
	//printf("New Transition\n");
	for (n = numTransitions, c=0; c<3; c++) TTransitionsClicked(n, c);
	Transitions[n]->setLowerState(lState);
	Transitions[n]->setUpperState(uState);
	Transitions[n]->addLineTable(L);
	L->show();
	Changed();
	//printf("Ende getLineTable\n");
	return L;
}

LineTable* Molecule::getLineTable(QString Name)
{
	int t, l, N;
	for (t=0; t < numTransitions; t++) for (l=0, N = Transitions[t]->getNumLineTables(); l<N; l++)
		if (Transitions[t]->getLineTableName(l) == Name) return Transitions[t]->getLineTable(l);
	return 0;
}

FCFTab* Molecule::getFCFTable(ElState* lState, ElState* uState)
{
	int n;
	for (n=0; n < numTransitions; n++)
		if (Transitions[n]->getLowerState() == lState && Transitions[n]->getUpperState() == uState)
			return Transitions[n]->getFCFTable();
	return 0;
}

bool Molecule::checkAllConnections()
{
    int n, m;
    bool Ret = true;
    FitData *fitData = new FitData(0, MW, this);
    for (n=0; n < numStates; ++n)
    {
        for (m=0; m < States[n]->getNumTermTables(); ++m) if (States[n]->getTermTable(m) == 0) Ret = false;
        for (m=0; m < States[n]->getNumDunTables(); ++m) if (States[n]->getDunTable(m) == 0) Ret = false;
        for (m=0; m < States[n]->getNumPotentials(); ++m) if (States[n]->getPotential(m) == 0) Ret = false;
        for (m=0; m < States[n]->getNumFitDataSets(); ++m)
        {
            if (fitData->readData(States[n]->getFitDataFileName(m)))
            {
                if (!fitData->checkAllConnections()) Ret = false;
                if (!fitData->isSaved()) fitData->writeData();
            }
            else Ret = false;
        }
    }
    delete fitData;
    for (n=0; n < numTransitions; ++n)
    {
        for (m=0; m < Transitions[n]->getNumLineTables(); ++m)
        {
            if (Transitions[n]->getLineTable(m) == 0) Ret = false;
            else if (!Transitions[n]->getLineTable(m)->checkAllConnections()) Ret = false;
        }
        for (m=0; m < Transitions[n]->getNumFCFTables(); ++m) if (Transitions[n]->getFCFTable(m) == 0) Ret = false;
    }
    return Ret;
}

bool Molecule::FCFavailable()
{
	int n;
	for (n=0; n < numTransitions; n++) if (Transitions[n]->getNumFCFTables() > 0) return true;
	return false;
}

int Molecule::getNumStates()
{
	return numStates;
}

int Molecule::getNumIso()
{
	int n = atom1->getnIso(), m = atom2->getnIso();
	return (atom1 == atom2 ? (n * (n + 1)) / 2 : n * m);
}

IsoTab *Molecule::getIso()
{
	if (atom1 == 0)
	{
		printf("Molecule::getIso(): atom1 is 0, can't calculate the data!\n");
		return 0;
	}
	if (atom2 == 0)
	{
		printf("Molecule::getIso(): atom2 is 0, can't calculate the data!\n");
		return 0;
	}
	IsoTab *rIso;
	int i, j, n, m, k=0, numIso, ib;
	double NA, m1, m2;
	QString B;
	TermTable *xmTerm = (States[0]->getNumTermTables() > 0 ? States[0]->getTermTable() : 0);
	n = atom1->getnIso();
	if (atom1 == atom2)
	{
		numIso = (n * (n + 1)) / 2;
		rIso = new IsoTab(numIso);
		*rIso->chSymb1 = atom1->getChSymb();
		*rIso->chSymb2 = atom2->getChSymb();
		for (i=0; i<n; i++) for (j=i; j<n; j++)
		{
			rIso->mNumIso2[k] = atom1->getnNuc(j);
			rIso->mNumIso1[k] = rIso->mNumIso2[i];
			rIso->mIso1[k] = atom1->getIsoMass(i);
			rIso->mIso2[k] = atom2->getIsoMass(j);
			if (i == j)
			{
				rIso->texName[k] = "^{" + QString::number(rIso->mNumIso1[k]) + "}" 
						+ *rIso->chSymb1 + "_2";
				NA = atom1->getIsoNA(i);
				rIso->relNA[k] = NA * NA;
				if (atom1->getCoreSpin(i) == 0) rIso->JStep[k] = 2;
				else rIso->JStep[k] = 1; 
				rIso->redMass[k++] = atom1->getIsoMass(i) * .5;
			}
			else
			{
				rIso->texName[k] = "^{" + QString::number(rIso->mNumIso1[k]) + "}" + *rIso->chSymb1 
								 + "^{" + QString::number(rIso->mNumIso2[k]) + "}" + *rIso->chSymb1;
				rIso->relNA[k] = 2 * atom1->getIsoNA(i) * atom1->getIsoNA(j);
				rIso->JStep[k] = 1;
				m1 = atom1->getIsoMass(i);
				m2 = atom2->getIsoMass(j);
				rIso->redMass[k++] = m1 * m2 / (m1 + m2);
			}
		}
	}
	else
	{
		m = atom2->getnIso();
		numIso = n * m;
		rIso = new IsoTab(numIso);
		*rIso->chSymb1 = atom1->getChSymb();
		*rIso->chSymb2 = atom2->getChSymb();
		for (i=0; i<n; i++) for (j=0; j<m; j++)
		{
			rIso->mIso1[k] = atom1->getIsoMass(i);
			rIso->mIso2[k] = atom2->getIsoMass(j);
			rIso->mNumIso1[k] = atom1->getnNuc(i);
			rIso->mNumIso2[k] = atom2->getnNuc(j);
			rIso->texName[k] = "^{" + QString::number(rIso->mNumIso1[k]) + "}" + *rIso->chSymb1
					         + "^{" + QString::number(rIso->mNumIso2[k]) + "}" + *rIso->chSymb2;
			rIso->relNA[k] = atom1->getIsoNA(i) * atom2->getIsoNA(j);
			rIso->JStep[k] = 1;
			m1 = atom1->getIsoMass(i);
			m2 = atom2->getIsoMass(j);
			rIso->redMass[k++] = m1 * m2 / (m1 + m2);
		}
	}
	if ((k = RefIso->currentIndex()) == -1) for (i = 0, NA = 0.0; i < numIso; i++) 
			if (rIso->relNA[i] > NA) NA = rIso->relNA[k = i];
	rIso->refIso = k;
	if (xmTerm != 0)
	{
		n = xmTerm->getNumIso();
		for (m=0; m<n; m++)
		{
			xmTerm->GetIsoZ(m, i, j);
			if (i != rIso->mNumIso1[m] || j != rIso->mNumIso2[m])
			{
				for (k=m+1; (k < numIso ? i != rIso->mNumIso1[k] || j != rIso->mNumIso2[k] : false); k++) ;
				if (k < numIso)
				{
					m1 = rIso->mIso1[m];
					rIso->mIso1[m] = rIso->mIso1[k];
					rIso->mIso1[k] = m1;
					m1 = rIso->mIso2[m];
					rIso->mIso2[m] = rIso->mIso2[k];
					rIso->mIso2[k] = m1;
					ib = rIso->mNumIso1[m];
					rIso->mNumIso1[m] = rIso->mNumIso1[k];
					rIso->mNumIso1[k] = ib;
					ib = rIso->mNumIso2[m];
					rIso->mNumIso2[m] = rIso->mNumIso2[k];
					rIso->mNumIso2[k] = ib;
					B = rIso->texName[m];
					rIso->texName[m] = rIso->texName[k];
					rIso->texName[k] = B;
					m1 = rIso->relNA[m];
					rIso->relNA[m] = rIso->relNA[k];
					rIso->relNA[k] = m1;
					ib = rIso->JStep[m];
					rIso->JStep[m] = rIso->JStep[k];
					rIso->JStep[k] = ib;
					m1 = rIso->redMass[m];
					rIso->redMass[m] = rIso->redMass[k];
					rIso->redMass[k] = m1;
					if (rIso->refIso == m) rIso->refIso = k;
					else if (rIso->refIso == k) rIso->refIso = m;
				}
			}
		}
	}
	for (i = 0; i < numIso; i++) 
	{
		rIso->relRedMass[i] = rIso->redMass[rIso->refIso] / rIso->redMass[i];
		rIso->rootRRM[i] = sqrt(rIso->relRedMass[i]);
	}
	return rIso;
}

void Molecule::getRefIso(int &A1, int &A2)
{
	if (atom1 == 0 || atom2 == 0)
	{
		A1 = A2 = 0;
		return;
	}
	int N1 = atom2->getnIso(), rIso = RefIso->currentIndex();
	if (atom1 != atom2)
	{
		A1 = rIso / N1;
		A2 = rIso - A1 * N1;
	}
	else
	{
		int i, j;
		for (i=j=0; j <= rIso; i++) j += N1 - i;
		A1 = i - 1;
		A2 = N1 - j + rIso;
	}
}

void Molecule::setRefIso(int A1, int A2)
{
	if (atom1 == 0 || atom2 == 0 || A1 >= atom1->getnIso() || A2 >= atom2->getnIso())
	{
		printf("Error setRefIso!\n");
		return;
	}
	if (atom1 != atom2) RefIso->setCurrentIndex(A1 * atom2->getnIso() + A2);
	else
	{
		int i, j, N = atom1->getnIso();
		if (A1 > A2)
		{
			i = A1;
			A1 = A2;
			A2 = i;
		}
		for (i=j=0; i < A1; i++) j += N-i;
		RefIso->setCurrentIndex(j + A2 - i);
	}
}

int Molecule::getJStep(int Iso)
{
	if (atom1 == 0)
	{
		printf("Molecule::getJStep(int): atom1 is 0, can't calculate the data!\n");
		return 1;
	}
	if (atom2 == 0)
	{
		printf("Molecule::getJStep(int): atom2 is 0, can't calculate the data!\n");
		return 1;
	}
	if (atom1 != atom2) return 1;
	int i, j = atom1->getnIso(), k=0;
	for (i=0; i < Iso && j > 0; i += j--) k++;
	if (i == Iso && atom1->getCoreSpin(k) == 0) return 2;
	return 1;
}

double Molecule::getIsoMass(int Iso)
{
	if (atom1 == 0)
	{
		printf("Molecule::getIsoMass(int): atom1 is 0, can't calculate the data!\n");
		return 0.0;
	}
	if (atom2 == 0)
	{
		printf("Molecule::getIsoMass(int): atom2 is 0, can't calculate the data!\n");
		return 0.0;
	}	
	int n = atom1->getnIso(), m = atom2->getnIso();
	if (atom1 == atom2) 
	{
		while (Iso >= n) Iso -= n--;
		return atom1->getIsoMass(m-n) + atom2->getIsoMass(Iso);
	}
	m = Iso / n;
	return atom1->getIsoMass(m) + atom2->getIsoMass(Iso-m*n);
}

void Molecule::updateAtoms()
{
	//printf("Update Atoms\n");
	UNI = false;
	int i, i1 = -1, i2 = -1, N = MW->getNumAtoms();
	QString AN;
	Atom *A;
	BAtom1->clear();
	BAtom2->clear();
	for (i=0; i<N; i++)
	{
		if ((A = MW->getAtom(i)) == atom1) i1 = i;
		if (A == atom2) i2 = i;
		BAtom1->insertItem(i, AN = A->getName());
		BAtom2->insertItem(i, AN);
	}
	//printf("i1=%d, i2=%d, N=%d\n", i1, i2, N);
	if (i1 != -1) BAtom1->setCurrentIndex(i1);
	if (i2 != -1) BAtom2->setCurrentIndex(i2);
	UNI = true;
	updateNI();
	//printf("Vor Update TermBox\n");
	UpdateTermBox();
	//printf("Vor Update DunBox\n");
	UpdateDunBox();
	UpdatePotBox();
	UpdateFitDataBox();
	//printf("Vor Update LineBox\n");
	UpdateLineBox();
	UpdateFCFBox();
	//printf("Vor Update StateBox\n");
	UpdateStateBox();
	//printf("Ende von updateAtoms()\n");
}

bool Molecule::writeData(QString NFilename)
{
	int i, j;
	bool E;
    QString N = Name->text(), S[10], TN, FN;
	QFile Datei(NFilename);
	if (!write(&Datei)) return false;
	QTextStream St(&Datei);
    QString FileName = getFileName();
    St << "Atom 1: " << (atom1 != 0 ? atom1->getRelativePath(FileName) : "") << "\n";
    St << "Atom 2: " << (atom2 != 0 ? atom2->getRelativePath(FileName) : "") << "\n";
    St << "Reference isotopologue: " << RefIso->currentText() << "\n";
	St << Name->text() << "\n";
	St << "\nList of term energy tables:\n";
	St << "State: | name: | file name: | source:\n";
	for (i=0; i < numStates; i++)
	{
		TN = States[i]->getName();
		for (j=0; j < States[i]->getNumTermTables(); j++)
			if ((States[i]->isTermTableLoaded(j) ? States[i]->getTermTable(j)->isOnDisk() : true)
                    && !(FN = getRelativePath(FN = States[i]->getTermTableFileName(j), FileName)).isEmpty())
				St << TN << " | " << States[i]->getTermTableName(j) << " | " 
					<< FN << " | " << States[i]->getTermTableSource(j)
					<< "\n";
	}
	St << "\nList of Dunham coefficient sets:\n";
	St << "State: | name: | file name: | source:\n";
	for (i=0; i < numStates; i++) 
	{
		TN = States[i]->getName();
		for (j=0; j < States[i]->getNumDunTables(); j++)
			if ((States[i]->isDunTableLoaded(j) ? States[i]->getDunTable(j)->isOnDisk() : true)
                    && !(FN = getRelativePath(FN = States[i]->getDunTableFileName(j), FileName)).isEmpty())
				St << TN << " | " << States[i]->getDunTableName(j) << " | " 
					<< FN << " | " << States[i]->getDunTableSource(j) 
					<< "\n";
	}
	St << "\nList of potentials:\n";
	St << "State: | name: | file name: | source:\n";
	for (i=0; i < numStates; i++)
	{
		TN = States[i]->getName();
		for (j=0; j < States[i]->getNumPotentials(); j++)
			if ((States[i]->isPotentialLoaded(j) ? States[i]->getPotential(j)->isOnDisk() : true)
                    && !(FN = getRelativePath(FN = States[i]->getPotentialFileName(j), FileName)).isEmpty())
				St << TN << " | " << States[i]->getPotentialName(j) << " | " 
					<< FN << " | " << States[i]->getPotentialSource(j)
					<< "\n";
	}
	St << "\nList of fit datasets:\n";
	St << "State: | name: | file name: | source:\n";
	for (i=0; i < numStates; i++)
	{
		TN = States[i]->getName();
		for (j=0; j < States[i]->getNumFitDataSets(); j++)
			if ((States[i]->isFitDataLoaded(j) ? States[i]->getFitData(j)->isOnDisk() : true)
                    && !(FN = getRelativePath(FN = States[i]->getFitDataFileName(j), FileName)).isEmpty())
				St << TN << " | " << States[i]->getFitDataName(j) << " | "
				   << FN << " | " << States[i]->getFitDataSource(j) << "\n";
	}
	St << "\nList of line tables:\n";
	St << "Upper state: | lower state: | name: | file name: | source:\n";
	for (i=0; i < numTransitions; i++)
	{
		TN = (Transitions[i]->getUpperState() != 0 ? Transitions[i]->getUpperState()->getName() : "")
		   + " | " 
		   + (Transitions[i]->getLowerState() != 0 ? Transitions[i]->getLowerState()->getName() : "")
		   + " | ";
		for (j=0; j < Transitions[i]->getNumLineTables(); j++)
			if ((Transitions[i]->isLineTableLoaded(j) ? Transitions[i]->getLineTable(j)->isOnDisk() : true)
                    && !(FN = getRelativePath(FN = Transitions[i]->getLineTableFileName(j), FileName)).isEmpty())
				St << TN << Transitions[i]->getLineTableName(j) << " | " 
					<< FN << " | " << Transitions[i]->getLineTableSource(j) << "\n";
	}
	St << "\nList of FCF tables:\n";
	St << "Upper state: | lower state: | name: | file name: | source:\n";
	for (i=0; i < numTransitions; i++)
	{
		TN = (Transitions[i]->getUpperState() != 0 ? Transitions[i]->getUpperState()->getName() : "")
		   + " | "
		   + (Transitions[i]->getLowerState() != 0 ? Transitions[i]->getLowerState()->getName() : "")
		   + " | ";
		for (j=0; j < Transitions[i]->getNumFCFTables(); j++)
			if ((Transitions[i]->isFCFTableLoaded(j) ? Transitions[i]->getFCFTable(j)->isOnDisk() : true)
                    && !(FN = getRelativePath(FN = Transitions[i]->getFCFTableFileName(j), FileName)).isEmpty())
				St << TN << Transitions[i]->getFCFTableName(j) << " | "
					<< FN << " | " << Transitions[i]->getFCFTableSource(j) << "\n";
	}
	St << "\nTable of states:\n";
    St << "Name: | lambda: | S: | +/- - symmetry: | g/u - symmetry: | Omega: | term energies: | dunham coeff.: | ";
	St << "potential: | fit data:\n";
	for (i=0; i<numStates; i++)
	{
		S[0] = States[i]->getName();
		S[1] = QString::number(States[i]->getLambda());
		S[2] = QString::number(States[i]->getS(), 'f', 1);
		switch (States[i]->getParity())
		{
			case 0:
				S[3] = "";
				break;
			case 1:
				S[3] = "g";
				break;
			case 2:
				S[3] = "u";
				break;
		}
		switch (States[i]->getSymmetry())
		{
			case -1:
				S[4] = "-";
				break;
			case 0:
				S[4] = "";
				break;
			case 1:
				S[4] = "+";
				break;
		}
		S[5] = QString::number(States[i]->getOmega());
		S[6] = ((States[i]->isTermTableLoaded() ? States[i]->getTermTable()->isOnDisk() 
                    : true) ? getRelativePath(FN = States[i]->getTermTableFileName(), FileName) : "");
		S[7] = ((States[i]->isDunTableLoaded() ? States[i]->getDunTable()->isOnDisk()
                    : true) ? getRelativePath(FN = States[i]->getDunTableFileName(), FileName) : "");
		S[8] = ((States[i]->isPotentialLoaded() ? States[i]->getPotential()->isOnDisk()
                    : true) ? getRelativePath(FN = States[i]->getPotentialFileName(), FileName) : "");
		S[9] = ((States[i]->isFitDataLoaded() ? States[i]->getFitData()->isOnDisk()
                    : true) ? getRelativePath(FN = States[i]->getFitDataFileName(), FileName) : "");
		for (j=0, E = true; j<10; j++) 
		{
			if (!S[j].isEmpty()) E = false;
			else S[j] = " ";
		} 
		if (!E) 
			St << S[0] << " | " << S[1] << " | " << S[2] << " | " << S[3] << " | " << S[4] << " | " 
                    << S[5] << " | " << S[6] << " | " << S[7] << " | " << S[8] << " | " << S[9] << " | " << QString::number(States[i]->getBe(), 'f', 12) << "\n";
	}
	St << "\nTable of transitions:\n";
	St << "State 1: | State 2: | table of lines: | FCF: | Strength:\n";
	//printf("Vor Transitions\n");
	for (i=0; i<numTransitions; i++)
	{
		S[0] = (Transitions[i]->getLowerState() != 0 ? Transitions[i]->getLowerState()->getName() 
													 : "");
		S[1] = (Transitions[i]->getUpperState() != 0 ? Transitions[i]->getUpperState()->getName() 
													 : "");
		//printf("Vor Lines\n");
		S[2] = ((Transitions[i]->isLineTableLoaded() ? 
				 Transitions[i]->getLineTable()->isOnDisk() : true) ?
                 getRelativePath(FN = Transitions[i]->getLineTableFileName(), FileName) : "");
		S[3] = ((Transitions[i]->isFCFTableLoaded() ? Transitions[i]->getFCFTable()->isOnDisk()
                 : true) ? getRelativePath(FN = Transitions[i]->getFCFTableFileName(), FileName) : "");
		for (j=0, E = true; j<3; j++) 
		{
			if (!S[j].isEmpty()) E = false;
			else S[j] = " ";
		}
		if (!E) 
			St << S[0] << " | " << S[1] << " | " << S[2] << " | " << S[3] << " | " 
			   << Transitions[i]->getTransitionStrengthText() << " \n";
	}
	return true;
}

void Molecule::updateNI()
{
	if (!UNI) return;
	//printf("Molecule::updateNI\n");
	QString AName, N, name;
	bool c = false;
	if ((AName = BAtom1->currentText()).isEmpty()) return;
	//printf("Atom1=%s\n", AName.ascii());
	if ((atom1 != 0 ? atom1->getName() != AName : true)) 
	{
		atom1 = MW->getAtomN(AName);
		c = true;
	}
	if ((AName = BAtom2->currentText()).isEmpty()) return;
	//printf("Atom2=%s\n", AName.ascii());
	if ((atom2 != 0 ? atom2->getName() != AName : true)) 
	{
		atom2 = MW->getAtomN(AName);
		c = true;
	} 
	//printf("atom1=%d, atom2=%d\n", atom1, atom2);
	if (atom1 == 0 || atom2 == 0 || (!c && RefIso->count() > 0)) return;
	//printf("Nach Test: Name Atom1: %s, Name Atom2: %s\n", 
		//   atom1->getName().ascii(), atom2->getName().ascii());
	int ci = RefIso->currentIndex(), i, j, k=0, n = atom1->getnIso(), m = atom2->getnIso();
	//double MIso1[n], MIso2[m];
	//printf("ci=%d\n", ci);
	QString CS1 = atom1->getChSymb(), CS2 = atom2->getChSymb(), V;
	//printf("2\n");
	RefIso->clear();
	//printf("3\n");
	if (atom1 == atom2)
	{
		printf("atom1==atom2\n");
		if ((N = "name of molecule: " + (name = CS1 + "2")) != Name->text())
		{
			Name->setText(N);
			setName(name);
		}
		//printf("Nach Name\n");
		for (i=0; i<n; i++) for (j=i; j<n; j++) 
		{
			RefIso->insertItem(k++, V = (i==j ? 
					QString::number(atom1->getnNuc(i)) + CS1 + "2" :
					QString::number(atom1->getnNuc(i)) + CS1 + 
					QString::number(atom1->getnNuc(j)) + CS1));
			if (V == RIN && ci == -1) ci = (i * (2 * n - i + 1)) / 2 + j - i;
		}
	}
	else
	{
		printf("atom1!=atom2\n");
		if ((N = "name of molecule: " + (name = CS1 + CS2)) != Name->text())
		{
			Name->setText(N);
			setName(name);
		}
		for (i=0; i<n; i++) for (j=0; j<m; j++)
		{ 
			RefIso->insertItem(k++, V = QString::number(atom1->getnNuc(i)) + CS1 +
					QString::number(atom2->getnNuc(j)) + CS2);
			if (V == RIN && ci == -1) ci = (i * m + j);
		}
	}
	//printf("RIN = %s, V = %s, ci = %d\n", RIN.ascii(), V.ascii(), ci);
	if (ci != -1) RefIso->setCurrentIndex(ci);
	//for (i=0; i<n; i++) MIso1[i] = atom1->getIsoMass(i);
	//for (i=0; i<m; i++) MIso2[i] = atom2->getIsoMass(i);
	/*printf("numStates=%d\n", numStates);
	for (i=0; i < numStates; i++) if (States[i]->getPotential() != 0) 
			States[i]->getPotential()->setMolData(n, MIso1, m, MIso2);*/
	Changed();
	//printf("Ende von updateNI\n");
}

void Molecule::TStatesClicked(int r, int c)
{
	//printf("Molecule::TStatesClicked\n");
	if (r >= numStates) return;
	ComboBox *B;
	switch (c)
	{
		case 1:
			if ((B = States[r]->getTermBox()) == 0) 
			{
				States[r]->setTermBox(B = new ComboBox);
				B->setEditable(false);
				TStates->setCellWidget(r, c, B);
				Changed();
			}
			UpdateTermBox(B);
			break;
		case 2:
			if ((B = States[r]->getDunBox()) == 0)
			{
				States[r]->setDunBox(B = new ComboBox);
				B->setEditable(false);
				TStates->setCellWidget(r, c, B);
				Changed();
			}
			UpdateDunBox(B);
			break;
		case 3:
			if ((B = States[r]->getPotBox()) == 0)
			{
				States[r]->setPotBox(B = new ComboBox);
				B->setEditable(false);
				TStates->setCellWidget(r, c, B);
				Changed();
			}
			UpdatePotBox(B);
			break;
		case 4:
			if ((B = States[r]->getFitDataBox()) == 0)
			{
				States[r]->setFitDataBox(B = new ComboBox);
				B->setEditable(false);
				TStates->setCellWidget(r, c, B);
				Changed();
			}
			UpdateFitDataBox(B);
			break;
	}
}

void Molecule::setStateName(int n, QString Name)
{
    QTableWidgetItem *I = TStates->verticalHeaderItem(n), *J;
    QString N = Name.left(1);
    int i;
    for (i=0; i < numStates; ++i) if (i != n && ((J = TStates->verticalHeaderItem(i)) != 0) && J->text().left(1) == N) break;
    if (i < numStates)
    {
        QString Name2 = States[i]->getName();
        int j;
        for (j=0; j < Name.length() && j < Name2.length() && Name[j] == Name2[j]; ++j) ;
        if (j < Name.length() && j < Name2.length())
        {
            J->setText(N + Name2[j]);
            N += Name[j];
        }
    }
    if (I != 0) I->setText(N);
    else TStates->setVerticalHeaderItem(n, new QTableWidgetItem(N));
	//printf("State[%d] Name=%s\n", n, Name.toAscii().data());
}

void Molecule::updateTransitionName(int transIndex, int lowerStateIndex, int upperStateIndex)
{
    QTableWidgetItem *LI = TStates->verticalHeaderItem(lowerStateIndex);
    QTableWidgetItem *UI = TStates->verticalHeaderItem(upperStateIndex);
    QString Name = QString(LI != 0 ? LI->text() : "") + (LI != 0 || UI != 0 ? "<->" : "") + (UI != 0 ? UI->text() : "");
    QTableWidgetItem *I = TTransitions->verticalHeaderItem(transIndex);
	if (I != 0) I->setText(Name);
    else TTransitions->setVerticalHeaderItem(transIndex, new QTableWidgetItem(Name));
}

void Molecule::TTransitionsClicked(int r, int c)
{
	//printf("Molecule::TTransitionsClicked\n");
	if (r > numTransitions || MW == 0) return;
	if (r == numTransitions)
	{
		if (r == MaxTransitions)
		{
			printf("The maximum number of transitions has been reached!\n");
			return;
		}
		Transitions[numTransitions] = MW->CreateTransition();
        Transitions[numTransitions]->setTransNum(numTransitions);
		Transitions[numTransitions++]->setMW(MW, this);
	}
	ComboBox *B;
	switch (c)
	{
		case 0:
			if ((B = Transitions[r]->getLSB()) == 0) 
			{
				Transitions[r]->setLSB(B = new ComboBox);
				B->setEditable(false);
				TTransitions->setCellWidget(r, c, B);
			}
			UpdateStateBox();
			break;
		case 1:
			if ((B = Transitions[r]->getUSB()) == 0)
			{
				Transitions[r]->setUSB(B = new ComboBox);
				B->setEditable(false);
				TTransitions->setCellWidget(r, c, B);
			}
			UpdateStateBox();
			break;
		case 2:
			if ((B = Transitions[r]->getLineBox()) == 0)
			{
				Transitions[r]->setLineBox(B = new ComboBox);
				B->setEditable(false);
				TTransitions->setCellWidget(r, c, B);
			}
			UpdateLineBox();
			break;
		case 3:
			if ((B = Transitions[r]->getFCFBox()) == 0)
			{
				Transitions[r]->setFCFBox(B = new ComboBox);
				B->setEditable(false);
				TTransitions->setCellWidget(r, c, B);
			}
	}
	Changed();
}

void Molecule::ChangeObjName(MDIChild::Type type)
{
	switch (type)
	{
		case TermEnergyTable:
			UpdateTermBox();
			break;
		case DunhamTable:
			UpdateDunBox();
			break;
		case PotData:
			UpdatePotBox();
			break;
		case FitDataSet:
			UpdateFitDataBox();
			break;
		case LineTab:
			UpdateLineBox();
			break;
		case FranckCondonTable:
			UpdateFCFBox();
			break;
		default:
			printf("Molecule::ChangeObjName Error: The type %d gets not handeld by a molecule!", type);
			break;
	}
}

void Molecule::UpdateTermBox(QComboBox*)
{
	int n;
	for (n=0; n < numStates; n++) States[n]->refreshTermBox();
	/*if (MW == 0) return;
	//printf("T=%d\n", MW->getTermTable(0));
	int i, j, N, k, l;
	bool s;
	if (B == 0) s = true;
	else s = false;
	N = MW->getNumTermTables();
	TermTable *T;
	QString N1, N2;
	for (k=0; k < numStates || (s && k==0); k++)
	{
		//printf("Beginn Schleife\n");
		if (s) B = States[k].getTermBox();
		if (B != 0)
		{
			if (s) N1 = (States[k].termTable != 0 ? States[k].termTable->getName() : "");
			else N1 = B->currentText();
			//printf("States[%d].termTable=%d\n", k, States[k].termTable);
			B->clear();
			for (i=l=0, j=-1; i<N; i++) 
			{
				//printf("k=%d, i=%d, N=%d, N1=%s, N2=%s, j=%d\n", k, i, N, N1.ascii(), N2.ascii(), j); 
				N2 = (T = MW->getTermTable(i))->getName();
				if (!N1.isEmpty() && N1 == N2) j = l;
				//printf("T=%d\n", T);
				printf("k=%d, i=%d, N1=%s, N2=%s, j=%d\n", k, i, N1.ascii(), N2.ascii(), j); 
				if (!T->isAssigned() || j==l) 
				{
					B->addItem(N2);
					l++;
				}
			}
			//printf("Vor addItem\n");
			B->addItem("");
			printf("B->count=%d\n", B->count());
			if (!N1.isEmpty()) B->setCurrentIndex(j);
			else B->setCurrentIndex(B->count() - 1);
		}
	}
	//printf("Ende UpdatetermBox\n");*/
}

void Molecule::UpdateDunBox(QComboBox*)
{
	int n;
	for (n=0; n < numStates; n++) States[n]->refreshDunBox();
	/*if (MW == 0) return;
	int i, j, N, k, l;
	//printf("UpdateDunBox\n");
	bool s;
	if (B == 0) s = true;
	else s = false;
	N = MW->getNumDunTables();
	DunTable *T;
	QString N1, N2;
	for (k=0; k < numStates || (s && k==0); k++)
	{
		if (s) B = States[k].getDunBox();
		if (B != 0)
		{
			N1 = (s && States[k].DunK != 0 ? States[k].DunK->getName() : "");
			B->clear();
			for (i=l=0, j=-1; i<N; i++) 
			{
				N2 = (T = MW->getDunTable(i))->getName();
				if (!N1.isEmpty() && N1 == N2) j = l;
				if (!T->isAssigned() || j==l) 
				{
					B->addItem(N2);
					l++;
				}
			}
			B->addItem("");
			if (!N1.isEmpty()) B->setCurrentIndex(j);
			else B->setCurrentIndex(B->count() - 1);
		}
	}*/
}	

void Molecule::UpdatePotBox(QComboBox*)
{
	int n;
	for (n=0; n < numStates; n++) States[n]->refreshPotBox();
	/*if (MW == 0) return;
	int i, j, N, k, l;
	bool s;
	if (B == 0) s = true;
	else s = false;
	N = MW->getNumPotentials();
	Potential *P;
	QString N1, N2;
	for (k=0; k < numStates || (s && k==0); k++)
	{
		if (s) B = States[k].getPotBox();
		if (B != 0)
		{
			N1 = (s && States[k].Pot != 0 ? States[k].Pot->getName() : "");
			B->clear();
			for (i=l=0, j=-1; i<N; i++)
			{
				N2 = (P = MW->getPotential(i))->getName();
				if (!N1.isEmpty() && N1 == N2) j=l;
				if (!P->isAssigned() || j==l)
				{
					B->addItem(N2);
					l++;
				}
			}
			B->addItem("");
			if (!N1.isEmpty()) B->setCurrentIndex(j);
			else B->setCurrentIndex(B->count() - 1);
		}
	}*/
}

void Molecule::UpdateFitDataBox(QComboBox*)
{
	int n;
	for (n=0; n < numStates; n++) States[n]->refreshFitDataBox();
}

void Molecule::UpdateLineBox(QComboBox*)
{
	int n;
	for (n=0; n < numTransitions; n++) Transitions[n]->refreshLineBox();
	/*if (MW == 0) return;
	int i, j, N, k, l;
	bool s;
	if (B == 0) s = true;
	else s = false;
	N = MW->getNumLineTables();
	LineTable *T;
	QString N1, N2;
	for (k=0; k < numTransitions || (s && k==0); k++)
	{
		if (s) B = Transitions[k].getLineBox();
		if (B != 0)
		{
			//printf("Vor getName\n");
			N1 = (s && Transitions[k].Lines != 0 ? Transitions[k].Lines->getName() : "");
			B->clear();
			for (i=l=0, j=-1; i<N; i++) 
			{
				//printf("i=%d\n", i);
				N2 = (T = MW->getLineTable(i))->getName();
				//printf("Nach get Name\n");
				if (!N1.isEmpty() && N1 == N2) j = l;
				if (!T->isAssigned() || j==l) 
				{
					B->addItem(N2);
					l++;
				}
			}
			//printf("Vor Ende\n");
			B->addItem("");
			if (!N1.isEmpty()) B->setCurrentIndex(j);
			else B->setCurrentIndex(B->count() - 1);
		}
	}*/
}

void Molecule::UpdateFCFBox(QComboBox*)
{
	int n;
	for (n=0; n < numTransitions; n++) Transitions[n]->refreshFCFBox();
}

void Molecule::UpdateStateBox()
{
	printf("Molecule::UpdateStateBox(): numTransitions = %d\n", numTransitions);
	int n;
	for (n=0; n < numTransitions; n++) Transitions[n]->refreshStateBox();
	/*printf("Update StateBox()\n");
	int k;
	QString T[numStates + 1];
	for (k=1; k <= numStates; k++) T[k] = States[k-1].Name;
	T[0] = "unknown";
	for (k=0; k < numTransitions; k++) 
	{
		UpdateStateBox(Transitions[k].getUSB(), T, 
					   (Transitions[k].upperState != 0 ? Transitions[k].upperState->Name : ""));
		UpdateStateBox(Transitions[k].getLSB(), T, 
					   (Transitions[k].lowerState != 0 ? Transitions[k].lowerState->Name : ""));
	}
	printf("Ende von Update StateBox\n");*/
}

/*void Molecule::UpdateStateBox(QComboBox *B, QString *T, QString C)
{
	printf("Update StateBox(B, T)\n");
	if (B == 0) return;
	bool s;
	int i, j, n = B->count();
	printf("Vor erster Schleife\n");
	if (T==0)
	{
		T = new QString[numStates + 1];
		for (i=1; i <= numStates; i++) T[i] = States[i-1].Name;
		T[0] = "unknown";
		s = true;
	}
	else s = false;
	printf("Vor zweiter Schleife\n");
	B->clear();
	for (i=0; i <= numStates; i++) 
	{
		B->addItem(T[i]);
		if (T[i] == C) j = i;
	}
	B->setCurrentIndex(j);
	printf("Vor delete\n");
	if (s) delete[] T;
	printf("End of update StateBox\n");
}*/ 

void Molecule::shrinkAllSpectRefs()
{
    int n, m;
    for (n=0; n < numStates; ++n) for (m=0; m < States[n]->getNumFitDataSets(); ++m) if (States[n]->getFitData(m) != 0) States[n]->getFitData(m)->shrinkAllSpectRefs();
    for (n=0; n < numTransitions; ++n) for (m=0; m < Transitions[n]->getNumLineTables(); ++m) if (Transitions[n]->getLineTable(m) != 0)
        Transitions[n]->getLineTable(m)->shrinkAllSpectRefs();
}

void Molecule::TStatesChanged(QTableWidgetItem *I)
{
	int r = I->row(), c = I->column(), N, n;
	QString Name, T = I->text();
	QTableWidgetItem *HI;
	if ((r >= numStates && (T.isEmpty() || c > 0)) || MW == 0) return;
	//printf("TStateChanged: T=%s\n", T.ascii());
	switch (c) 
	{
		case 0:
			if (r > numStates) 
			{
				r = numStates;
				I->setText("");
			}
            if (r == numStates) addState(T, false);
			else States[r]->setName(T);
			if ((HI = TStates->verticalHeaderItem(r)) != 0) HI->setText(T.left(1));
			else TStates->setVerticalHeaderItem(r, new QTableWidgetItem(T.left(1)));
			break;
		case 2:
			N = MW->getNumDunTables();
			Name = T;
			for (n=0; n<N; n++) if (Name == MW->getDunTable(n)->getName())
					MW->getDunTable(n)->setElState(States[r]);
			break;
			//printf("TStatesChanges, c=%d\n", c);
		case 3:
			N = MW->getNumPotentials();
			Name = T;
			for (n=0; n<N; n++) if (Name == MW->getPotential(n)->getName())
					MW->getPotential(n)->setElState(States[r]);
			break;
	}
	Changed();
}

void Molecule::getKnownLevels(IsoTab *&Iso, int &NS, int &mv, int &mJ, bool ****&L)
{
	int mvs[numTransitions], mJs[numTransitions], mvss[numTransitions], mJss[numTransitions];
	bool ***uL[numTransitions], ***lL[numTransitions];
	Iso = getIso();
	ElState *St;
	int t, v, J, I, s;
	for (t=0, mv = 0, mJ = 0; t < numTransitions; t++) 
	{
		if (Transitions[t]->getLineTable() != 0)
		{
			Transitions[t]->getLineTable()->getKnownLevels(Iso->numIso, mvs[t], mvss[t], mJs[t],
										                  mJss[t], uL[t], lL[t]);
			if (mvs[t] > mv) mv = mvs[t];
			if (mvss[t] > mv) mv = mvss[t];
			if (mJs[t] > mJ) mJ = mJs[t];
			if (mJss[t] > mJ) mJ = mJss[t];
		}
		else 
		{
			mvs[t] = mJs[t] = mvss[t] = mJss[t] = -1;
			uL[t] = lL[t] = 0;
		}
	}
	//printf("1.\n");
	for (NS=0; (NS < numStates ? !States[NS]->getName().isEmpty() : false); NS++) ;
	L = CreateBool(Iso->numIso, NS, mv + 1, mJ + 1);
	for (I=0; I < Iso->numIso; I++) for (s=0; s < NS; s++) for (v=0; v <= mv; v++) 
				for (J=0; J <= mJ; J++) L[I][s][v][J] = false;
	//printf("2. numTransitions=%d\n", numTransitions);
	for (t=0; t < numTransitions; t++)
	{
		//printf("t=%d\n", t);
		for (s=0, St = Transitions[t]->getLowerState(); (s < NS ? States[s] != St : false); s++) ;
		//printf("t=%d, s1=%d, NS=%d\n", t, s, NS);
		if (s < NS) for (I=0; I < Iso->numIso; I++) for (v=0; v <= mvss[t]; v++) 
		{ //if (I==0) //printf("v=%d, mvss[t]=%d\n", v, mvss[t]);
					for (J=0; J <= mJss[t]; J++) 
					{ //if (v==50) //printf("J=%d, mJss=%d\n", J, mJss[t]);
						if (lL[t][I][v][J]) L[I][s][v][J] = true; }}
		for (s=0, St = Transitions[t]->getUpperState(); (s < NS ? States[s] != St : false); s++) ;
		//printf(" s2=%d\n", s);
		if (s < NS) for (I=0; I < Iso->numIso; I++) for (v=0; v <= mvs[t]; v++)
					for (J=0; J <= mJs[t]; J++) if (uL[t][I][v][J]) L[I][s][v][J] = true;
	}
	//printf("3.\n");
	//for (J=0; J < mJ; J++) if (uL[0][0][J]) printf("1.:mJ=%d,v'=0, J=%d vorhanden!\n", mJ, J);
	for (t=0; t < numTransitions; t++)
	{
		if (lL[t] != 0) Destroy(lL[t], Iso->numIso, mvss[t] + 1);
		if (uL[t] != 0) Destroy(uL[t], Iso->numIso, mvs[t] + 1);
	}
	//for (J=0; J <= mJ; J++) if (L[0][1][0][J]) printf("2.:v'=0, J=%d vorhanden!\n", J);
}
