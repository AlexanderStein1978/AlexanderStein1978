//
// C++ Implementation: aipsidiag
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2008 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "aipsidiag.h"
#include "MainWindow.h"
#include "potential.h"
#include "molecule.h"
#include "utils.h"
#include "constants.h"

#include <QLineEdit>
#include <QListWidget>
#include <QComboBox>
#include <QRadioButton>
#include <QCheckBox>
#include <QPushButton>
#include <QLabel>
#include <QFileDialog>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>

AIPSIDialog::AIPSIDialog(MainWindow *mw)
{
	int n, N;
	MW = mw;
	NFiles = 0;
	NCol = 0;
	NRow = 0;
	R = 0;
	E = 0;
	QSize Size(420, 452);
	setAttribute(Qt::WA_DeleteOnClose);
	setMinimumSize(Size);
	setMaximumSize(Size);
	QLabel *SLabel = new QLabel("Source description:", this);
	SLabel->setGeometry(10, 12, 150, 20);
	Source = new QLineEdit(this);
	Source->setGeometry(160, 5, 250, 30);
	QLabel *PathLabel = new QLabel("Input directory:", this);
	PathLabel->setGeometry(10, 50, 180, 20);
	Path = new QLineEdit(this);
	Path->setGeometry(5, 70, 180, 30);
	QLabel *FileLabel = new QLabel("Input files:", this);
	FileLabel->setGeometry(220, 50, 90, 20);
	QPushButton *FSButton = new QPushButton("Select files...", this);
	FSButton->setGeometry(215, 70, 90, 30);
	connect(FSButton, SIGNAL(clicked(bool)), this, SLOT(showFileDialog()));
	FileList = new QListWidget(this);
	FileList->setGeometry(310, 50, 100, 50);
	QLabel *EULabel = new QLabel("Energy unit:", this);
	EULabel->setGeometry(10, 118, 95, 20);
	EUBox = new QComboBox(this);
	EUBox->setGeometry(95, 112, 100, 30);
	EUBox->addItem("hartree");
	EUBox->addItem("cm^-1");
	QLabel *LULabel = new QLabel("Length unit:", this);
	LULabel->setGeometry(225, 118, 95, 20);
	LUBox = new QComboBox(this);
	LUBox->setGeometry(310, 112, 100, 30);
	LUBox->addItem("a0");
	LUBox->addItem("Angstrom");
	QLabel *RefLabel = new QLabel("Reference state:    File:", this);
	RefLabel->setGeometry(10, 156, 130, 20);
	PotFileBox = new QComboBox(this);
	PotFileBox->setGeometry(145, 150, 100, 30);
	connect(PotFileBox, SIGNAL(currentIndexChanged(int)), this, SLOT(PotFileBoxChanged(int)));
	QLabel *ColLabel = new QLabel("Column:", this);
	ColLabel->setGeometry(255, 156, 70, 20);
	PotNumBox = new QComboBox(this);
	PotNumBox->setGeometry(310, 150, 100, 30);
	QLabel *EOffset = new QLabel("Energy offset of the potential set (estimated by the reference state):",
								 this);
	EOffset->setGeometry(10, 194, 400, 20);
	PotOffsetMinSero = new QRadioButton("Set the potential minimum to sero", this);
	PotOffsetMinSero->setGeometry(10, 224, 400, 20);
	PotOffsetASero = new QRadioButton("Set the asymptote to sero", this);
	PotOffsetASero->setGeometry(10, 249, 400, 20);
	PotOffsetMatchA = new QRadioButton("Match the asymptote with the asymptote of:", this);
	PotOffsetMatchA->setGeometry(10, 274, 400, 20);
	PotentialBox = new QComboBox(this);
	PotentialBox->setGeometry(284, 268, 126, 30);
	for (n=0, N = MW->getNumPotentials(); n<N; n++) PotentialBox->addItem(MW->getPotential(n)->getName());
	PotentialBox->setEnabled(false);
	connect(PotOffsetMatchA, SIGNAL(toggled(bool)), PotentialBox, SLOT(setEnabled(bool)));
	PotOffsetConst = new QRadioButton("Do not change the energy offset of the potentials", this);
	PotOffsetConst->setGeometry(10, 299, 400, 20);
	PotOffsetMinSero->setChecked(true);
	AtMol = new QCheckBox("Add for each potential an electronic state to the molecule:", this);
	AtMol->setGeometry(10, 337, 400, 20);
	AtMol->setChecked(true);
	MoleculeBox = new QComboBox(this);
	MoleculeBox->setGeometry(284, 357, 126, 30);
	MoleculeBox->addItem("new molecule");
	for (n=0, N = MW->getNumMolecules(); n<N; n++) MoleculeBox->addItem(MW->getMolecule(n)->getName());
	connect(AtMol, SIGNAL(stateChanged(int)), this, SLOT(AtMolChecked(int)));
	QPushButton *OKButton = new QPushButton("OK", this);
	OKButton->setGeometry(10, 407, 150, 30);
	connect(OKButton, SIGNAL(clicked(bool)), this, SLOT(ImportData()));
	QPushButton *CancelButton = new QPushButton("Cancel", this);
	CancelButton->setGeometry(260, 407, 150, 30);
	connect(CancelButton, SIGNAL(clicked(bool)), this, SLOT(close()));
}

AIPSIDialog::~AIPSIDialog()
{
	if (E!=0) 
	{
		int i;
		Destroy(R, NFiles);
		for (i=0; i < NFiles; i++) Destroy(E[i], NCol[i]);
		delete[] E;
		delete[] NCol;
		delete[] NRow;
	}
}

void AIPSIDialog::AtMolChecked(int S)
{
	if (S == Qt::Unchecked) MoleculeBox->setEnabled(false);
	else MoleculeBox->setEnabled(true);
}

void AIPSIDialog::PotFileBoxChanged(int n)
{
	int N = NCol[n], i;
	for (i = PotNumBox->count() - 1; i>=N; i--) PotNumBox->removeItem(i);
	for (i = PotNumBox->count() + 1; i<=N; i++) PotNumBox->addItem(QString::number(i));
}

void AIPSIDialog::openFile(QString FileName, bool &GF, int &anRow, int &anCol, double *&aR, double **&aE)
{
	//printf("File=%s\n", FileName.ascii());
	int l, j, k, c, ac = -1, mc = 0, cc, lc;
	double lV, aV;
	bool ok;
	QStringList *BL = new QStringList[100000];
	QFile Datei(FileName);
	if (!Datei.open(QIODevice::ReadOnly)) return;
	QTextStream S(&Datei);
	//printf("Vor split\n");
	for (l=0; l < 100000 && !S.atEnd(); l++) //{ printf("l=%d\n", l);
		BL[l] = S.readLine().split(" ", QString::SkipEmptyParts);//}
	//printf("Nach split\n");
	if (l == 100000) printf("File %s is too long!", FileName.toLatin1().data());
	for (j=cc=c=0, aV = 0.0; j<l; j++)
	{
		lV = aV;
		lc = ac;
		if ((ac = BL[j].count()) == 0) 
		{
			ac = lc;
			continue;
		}
		aV = BL[j][0].toDouble(&ok);
		if (ac == lc && ok && aV > lV)
		{
			c++;
			cc++;
			if (ac > mc) mc = ac;
		}
		else
		{
			if (cc < 5 && c == cc) 
			{
				c=1;
				mc = ac;
			}
			cc = 1;
			if (j!=0) ac = lc;
			aV = lV;
		}
	}
	//printf("Nach erster Schleife\n");
	if (c < 5) 
	{
		delete[] BL;
		return;
	}
	GF = true;
	anRow = c;
	anCol = mc - 1;
	aR = new double[c];
	aE = Create(mc - 1, c);
	//printf("Vor zweiter Schleife, c=%d\n", c);
	for (j=cc=0, c=-1, aV = 0.0; j<l; j++)
	{
		lV = aV;
		lc = ac;
		if ((ac = BL[j].count()) == 0) 
		{
			ac = lc;
			continue;
		}
		aV = BL[j][0].toDouble(&ok);
		if (ac == lc && ok && aV > lV)
		{
			c++;
			cc++;
		}
		else
		{
			//printf("j=%d, ac=%d, lc=%d, ok=%d, aV=%f, lV=%f\n", j, ac, lc, ok, aV, lV);
			if (cc < 5 && c == cc) c=0;
			cc = 1;
			if (j!=0) ac = lc;
			aV = lV;
		}
		//printf("j=%d, c=%d\n", j, c);
		aR[c] = aV;
		for (k=1; k < BL[j].count() && k < mc; k++) aE[k-1][c] = BL[j][k].toDouble();
		for ( ; k < mc; k++) aE[k-1][c] = 0.0;
	}
	Datei.close();
	//printf("File closed\n");
	delete[] BL;
	//printf("Nach delete\n");
}

void AIPSIDialog::showFileDialog()
{
	int i, j, d, l, n, c;
	QString Dir = Path->text(), Text;
	QStringList rFiles;
	QStringList nFiles = 
			QFileDialog::getOpenFileNames(this, "Select one or more files with potential data to import",
										  Dir, "Data files (*)");
	if ((n = nFiles.count()) == 0) return;
	bool GF[n];
	int anRow[n], anCol[n], *nRow, *nCol;
	double *aR[n], **aE[n], **nR, ***nE;
	Files << nFiles;
	for (i=1, d = Files[0].lastIndexOf(DIRSEP); i < NFiles && d != -1; i++)
		while (d != -1 ? Files[i].left(d) != Files[0].left(d) : false) 
			d = Files[0].lastIndexOf(DIRSEP, d-1);
	Path->setText(Dir = Files[0].left(d+1));
	for (i=0; i<n; i++) GF[i] = false;
	for (i=0; i<n; i++) openFile(nFiles[i], GF[i], anRow[i], anCol[i], aR[i], aE[i]);
	for (i=n-1, l = Dir.length(); i>=0; i--) if (!GF[i]) 
	{
		Files.removeAt(NFiles + i);
		rFiles << nFiles[i].right(nFiles[i].length() - l);
	}
	if ((c = rFiles.count()) > 0)
	{
		if (c==1) Text = "file " + rFiles[0] + " does";
		else 
		{
			for (i=1, Text = "files " + rFiles[0]; i<c-1; i++) Text += ", " + rFiles[i];
			Text += ", and " + rFiles[i] + " do";
		}
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
								 "The " + Text + " not have the right file format!");
	}
	printf("Vor ende\n");
	if (c < n)
	{
		l = Files.count();
		nCol = new int[l];
		nRow = new int[l];
		nR = new double*[l];
		nE = new double**[l];
		for (i=0; i<NFiles; i++)
		{
			nCol[i] = NCol[i];
			nRow[i] = NRow[i];
			nR[i] = R[i];
			nE[i] = E[i];
		}
		delete[] NCol;
		delete[] NRow;
		delete[] E;
		delete[] R;
		NCol = nCol;
		NRow = nRow;
		R = nR;
		E = nE;
		for (i=0, j=NFiles, l = Dir.length(); i<n; i++) if (GF[i])
		{
			NCol[j] = anCol[i];
			NRow[j] = anRow[i];
			R[j] = aR[i];
			E[j++] = aE[i];
			Text = nFiles[i].right(nFiles[i].length() - l);
			FileList->addItem(Text);
			PotFileBox->addItem(Text);
		}
		NFiles = Files.count();
	}
}

void AIPSIDialog::ImportData()
{
	printf("Beginn importdata\n");
	if (NFiles == 0)
	{
		close();
		return;
	}
	printf("R[0][1]=%f\n", R[0][1]);
	int i, j, k, n, m, rf = PotFileBox->currentIndex(), rc = PotNumBox->currentIndex();
	double EOff;
	bool AM;
	Molecule *mol = 0;
	QString SText = Source->text(), File, FileN, Name;
	Potential RefPot(0, 0), *pot;
	printf("Vor LU\n");
	if (LUBox->currentIndex() == 0) 
		for (i=0; i < NFiles; i++) for (j=0; j < NRow[i]; j++) R[i][j] *= a0_Angstrom;
	printf("Vor eu\n");
	if (EUBox->currentIndex() == 0) for (i=0; i < NFiles; i++) for (j=0; j < NRow[i]; j++)
				for (k=0; k < NCol[i]; k++) E[i][k][j] *= hartree_cm;
	printf("Vor refpot, rf=%d, rc=%d\n", rf, rc);
	RefPot.setPointData(NRow[rf], R[rf], E[rf][rc]);
	printf("Vor potoffset\n");
	if (!PotOffsetConst->isChecked())
	{
		if (PotOffsetASero->isChecked()) EOff = RefPot.getAsymptote();
		else if (PotOffsetMinSero->isChecked()) EOff = RefPot.getMinimumE();
		else EOff = RefPot.getAsymptote() 
					- MW->getPotential(PotentialBox->currentIndex())->getAsymptote();
		for (i=0; i < NFiles; i++) for (j=0; j < NRow[i]; j++) for (k=0; k < NCol[i]; k++) 
					E[i][k][j] -= EOff;
		printf("EOff=%f\n", EOff);
	}
	printf("Vor atmol\n");
	if (AtMol->isChecked())
	{
		AM = true;
		if ((i = MoleculeBox->currentIndex()) == 0)
		{
			mol = MW->CreateMolecule();
			mol->show();
		}
		else mol = MW->getMolecule(i-1);
	}
	else AM = false;
	printf("Vor letzter Schleife\n");
	for (i=0; i < NFiles; i++) 
	{
		File = Files[i];
		if ((j = File.indexOf(".")) != -1) File = File.left(j);
		for (j=0; j < NCol[i]; j++)
		{
			FileN = File + "_" + QString::number(j+1);
			for (n=0; (m = File.indexOf(DIRSEP, n)) != -1; n=m+1) ;
			Name = File.right(File.length() - n);
			pot = MW->CreatePotential();
			pot->setPointData(NRow[i], R[i], E[i][j]);
			pot->setSource(SText);
			pot->setName(Name);
			pot->writeData(FileN + ".dat");
			pot->show();
			if (AM) 
			{
                k = mol->addState(QString::number(j+1) + Name, false);
				if (k != -1) mol->addPot(k, pot);
			}
		}
	}
	if (!MW->isVisible()) MW->show();
	printf("Vor close\n");
	close();
}
