//
// C++ Implementation: mdichild
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#include "mdichild.h"
#include "MainWindow.h"

#include <QMessageBox>
#include <QCloseEvent>
#include <QHideEvent>
#include <QShowEvent>
#include <QFile>
#include <QFileDialog>
#include <QPainter>

MDIChild::MDIChild(Type ntype, MainWindow *mw, QString filter, QString FE) : QWidget(nullptr), m_fileNameChanged(false)
{
	//printf("MDIChild::MDIChild, mw=%d\n", mw);
	changed = true;
	MW = mw;
	Filter = filter;
	FileExt = FE;
	type = ntype;
	imported = false;
	newCreated = false;
	blockChangeSignal = false;
	onDisk = false;
	setWindowTitle(getTypeString());
}

MDIChild::~MDIChild()
{
}

bool MDIChild::askForQuit()
{
	return true;
}

void MDIChild::setMainWindow(MainWindow *mw)
{
	MW = mw;
}

void MDIChild::closeEvent(QCloseEvent *i_event)
{
    if (0 != MW) MW->MdiChildClosed(this);
    i_event->accept();
}

void MDIChild::Changed()
{
	//printf("MDIChild::Changed() %s %s\n", Type.toAscii().data(), Name.toAscii().data());
	changed = true;
	if (!blockChangeSignal) emit propertiesChanged();
}

void MDIChild::Print(QPrinter &)
{
	QMessageBox::information(this, tr("MolSpektAnalysis"), 
						 "Printing is not implemented for the active window");
}					 

void MDIChild::setBlockChangeSignal(bool Block)
{
	blockChangeSignal = Block;
}

bool MDIChild::getBlockChangeSignal()
{
	return blockChangeSignal;
}

bool MDIChild::isSaved()
{
	if (!changed) return true;
	return false;
}

QString MDIChild::getFileName(bool save)
{
	//printf("getFilename: Filename=%s\n", FileName.ascii());
	if (imported && save) writeData();
	return FileName;
}

void MDIChild::setFileName(QString N)
{
	FileName = N;
	emit fileNameChanged();
	Changed();
}

QString MDIChild::getFName()
{
	return FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
}

QString MDIChild::getFilter()
{
	return Filter;
}

QString MDIChild::getName() const
{
	return Name;
}

MDIChild::Type MDIChild::getType()
{
	return type;
}

QString MDIChild::getTypeString()
{
	switch (type)
	{
		case AddSpect:
			return "AddSpectrum";
		case AtomData:
			return "Atom";
        case CouplingFunction:
            return "CouplingFunction";
		case DataSetPlot:
			return "DataPlot";
		case DunhamTable:
			return "DunhamTable";
		case FCFJDependency:
			return "FCFJDependency";
		case FitDataSet:
			return "FitDataSet";
		case FitSeriesResultTable:
			return "FitSeriesResultTable";
		case FranckCondonTable:
			return "FCFTable";
		case FranckCondonView:
			return "FCFView";
		case IntensityHistogramPlot:
			return "IntensityHistogram";
		case LineTab:
			return "LineTable";
		case MolData:
			return "Molecule";
		case Pict:
			return "Picture";
		case PotData:
			return "Potential";
		case ResidPlot:
			return "ResidualPlot";
		case SimpleDiagWindow:
			return "DiagWindow";
		case SimpleMDIChild:
			return "MDIChild";
		case Spect:
			return "Spectrum";
		case SpectrumSimulation:
			return "SpectrumSimulation";
		case TermEnergyPlot:
			return "TermEnergyPlot";
		case TermEnergyTable:
			return "TermEnergyTable";
		case TermEnergyView:
			return "TermEnergyView";
		case TextTable1:
		case TextTable2:
			return "TextTable";
		case WaveFunctionPlot:
			return "WaveFunctionPlot";
			break;
        default:
            printf("Error: MDIChild::getTypeString: for the current type %d is no name string defined!", type);
            break;
	}

	return "";
}

void MDIChild::hideEvent(QHideEvent *E)
{
	printf("MDIChild::HideEvent\n");
	emit hidden(true);
	E->accept();
}

bool MDIChild::readData(QString /*Filename*/)
{
	printf("Warning: For the object %s of type %s is readData not reimplemented and thus the data can't be read!\n", 
		   Name.toLatin1().data(), getTypeString().toLatin1().data());
	return true;
}

bool MDIChild::read(QFile *Datei)
{
	int n, m, l, nc;
	QString Dir, Filename = Datei->fileName(), fileName, EB, TypeS = getTypeString();
	if (Filename != FileName) 
	{
		if (Filename.isEmpty()) Datei->setFileName(FileName);
		else
		{
			FileName = Filename;
			emit fileNameChanged();
		}
	}
	if (!FileName.isEmpty() && !Datei->exists() && MW != 0)
	{
		for (n = MW->getNDirRep(), m=0; !Datei->exists() && m<n; m++)
		{
			MW->getDirRep(m, Dir, nc);
			if (nc < (l = FileName.length())) 
				Datei->setFileName(fileName = Dir + FileName.right(l - nc));
		}
		if (Datei->exists()) FileName = fileName;
	}
	if (FileName.isEmpty() || !Datei->exists())
	{
		Dir = (MW != 0 ? MW->getDir(type) : "");
		if (FileName.isEmpty())
			fileName = QFileDialog::getOpenFileName(this, "Open " + TypeS, Dir + FileName, Filter);
		else
		{
			EB = FileName.right(FileName.length() - FileName.lastIndexOf('.'));
			fileName = QFileDialog::getOpenFileName(this, "Please show me " + FileName, Dir,
							(Filter.indexOf(EB) != -1 ? Filter : (EB = TypeS + " (*" + EB + ");;" + TypeS + " (*.dat);;" + Filter)));
			//printf(EB.ascii());
		}
		if (fileName.isEmpty()) return false;
		if (MW != 0)
		{
			for (n = fileName.length() - 1, m = FileName.length() - 1; 
				 (n>=0 && m>=0 ? fileName[n] == FileName[m] : false) ; n--, m--) ;
			if (fileName.length() - n > 1) MW->addDirRep(fileName.left(n+1), m+1); 
			if (!MW->checkNotOpened(fileName)) return false;
		}
		Datei->setFileName(FileName = fileName);
        m_fileNameChanged = true;
		emit fileNameChanged();
	}
    else m_fileNameChanged = false;
	if (!Datei->open(QIODevice::ReadOnly))
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), tr("Error while opening file!"));
		return false;
	}
#ifdef Q_WS_WIN
	for (n = Filename.indexOf("\\"); -1 < (m = Filename.indexOf("\\", n + 1)); n = m) ;
#else
	for (n = Filename.indexOf('/'); -1 < (m = Filename.indexOf('/', n + 1)); n = m) ;
#endif
	setWindowTitle(Filename.right(Filename.length() - n - 1));
	QString nDir = Filename.left(n+1);
	if (MW != 0 && Dir != nDir) MW->setDir(nDir, type);
	changed = false;
	onDisk = true;
	imported = false;
    return true;
}

QString MDIChild::getRelativePath(const QString& CurrentPath, const QString &MolPath)
{
    QStringList MolList = MolPath.split(QRegExp("[\\/]"), Qt::SkipEmptyParts), LocalList = CurrentPath.split(QRegExp("[\\/]"), Qt::SkipEmptyParts);
    int n, m;
    for (n=0; n < MolList.length() && n < LocalList.length() && MolList[n] == LocalList[n]; ++n) ;
    if (n==0) return CurrentPath;
    QString Result;
    for (m=n; m < MolList.length() - 1; ++m) Result += QString("../");
    for (m=n; m < LocalList.length(); ++m) Result += LocalList[m] + (m < LocalList.length() - 1 ? "/" : "");
    return Result;
}

QString MDIChild::getAbsolutePath(QString &CurrentPath, QString &MolPath)
{
    if (CurrentPath.indexOf('(') > 0 || (CurrentPath.length() > 0 && CurrentPath[0] == '/') || (CurrentPath.length() > 1 && CurrentPath[1] == ':'))
        return CurrentPath;
    int m, n, s, l = MolPath.length();
    for (n=s=0; CurrentPath.mid(s, 3) == "../"; s+=3, ++n) ;
    for (m=0; m<=n; ++m) l = MolPath.lastIndexOf(QRegExp("[\\/]"), l-1);
	if (l == -1) return CurrentPath;
    QStringList L = CurrentPath.split(QRegExp("[\\/]"));
    QString Result = MolPath.left(l);
    while (n < L.length()) Result += DIRSEP + L[n++];
    return Result;
}

void MDIChild::Saved()
{
	//printf("MDIChild::Saved %s %s\n", Type.toAscii().data(), Name.toAscii().data());
	changed = false;
	onDisk = true;
}

bool MDIChild::isAssigned()
{
	return false;
}

void MDIChild::show()
{
	if (isMinimized()) setWindowState(windowState() & (~Qt::WindowMinimized | Qt::WindowActive));
	else if (isVisible() && MW != 0) MW->setActive(this);
	else QWidget::show();
}

void MDIChild::showEvent(QShowEvent *E)
{
	//printf("MDIChild::ShowEvent\n");
	emit hidden(false);
    if (MW != 0) MW->MdiChildShown(this);
	E->accept();
}

bool MDIChild::writeData(QString /*Filename*/)
{
	QMessageBox::information(this, "MolSpektAnalysis", "The active window cannot be saved!");
	return true;
}

bool MDIChild::write(QFile *Datei)
{
	int n, m, l = FileExt.length();
	//printf("MDIChild::write: FileName=%s\n", FileName.ascii());
	QString Dir, NFilename = Datei->fileName(), TypeS = getTypeString();
	if (!NFilename.isEmpty()) 
	{
		if (l>0 ? NFilename.right(l) != FileExt : false) 
		{
			FileName = (NFilename.right(4) == ".dat" ? NFilename.left(NFilename.length() - 4)
							: NFilename) + FileExt;
			Datei->setFileName(FileName);
		}
		else FileName = NFilename;
		imported = false;
		//printf("new Filename=%s\n", FileName.ascii());
		emit fileNameChanged();
	}
	else if (!FileName.isEmpty() && !imported && !newCreated) 
	{
		if (l>0 ? FileName.right(l) != FileExt : false)
		{
			FileName = (FileName.right(4) == ".dat" ? FileName.left(FileName.length() - 4)
							: FileName) + FileExt;
			emit fileNameChanged();
		}
		Datei->setFileName(FileName);
	}
	else 
	{
		Dir = (FileName.isEmpty() ? (MW != 0 ? MW->getDir(type) + Name : Name) + FileExt : FileName); 
		QString fileName = QFileDialog::getSaveFileName(this, "Save " + TypeS + ' ' + Name, Dir, Filter);
		if (fileName.isEmpty()) return false;
		if (l>0 ? fileName.right(l) != FileExt : false) 
			FileName = (fileName.right(4) == ".dat" ? fileName.left(fileName.length() - 4)
							: fileName) + FileExt;
		else FileName = fileName;
		emit fileNameChanged();
		Datei->setFileName(FileName);
	}
	if (!Datei->open(QIODevice::WriteOnly | QIODevice::Truncate))
	{
		QMessageBox::warning(this, "QT4MolSpektAn", "Error while opening/creating file.\n " + TypeS + " " + Name + " has not been saved!");
		return false;
	}
#ifdef Q_WS_WIN
	for (n = FileName.indexOf("\\"); -1 < (m = FileName.indexOf("\\", n + 1)); n = m) ;
#else
	for (n = FileName.indexOf('/'); -1 < (m = FileName.indexOf('/', n + 1)); n = m) ;
#endif
	setWindowTitle(FileName.right(FileName.length() - n - 1));
	QString nDir = FileName.left(n+1);
	if (nDir != Dir && MW != 0) MW->setDir(nDir, type);
	emit nameChanged(Name);
	changed = false;
	imported = false;
	newCreated = false;
	onDisk = true;
	return true;
}

void MDIChild::setName(QString name)
{
	bool C = (name != Name ? true : false);
	if (MW != 0) Name = MW->getChildName(name, this);
	else Name = name;
	if (C || name != Name)
	{
		Changed();
		emit nameChanged(Name);
	}
}

void MDIChild::setNewCreated()
{
	newCreated = true;
	onDisk = false;
}

void MDIChild::setFileExt(QString FE)
{
	FileExt = FE;
}

void MDIChild::setFilter(QString filter)
{
	Filter = filter;
}

void MDIChild::setType(Type ntype)
{
	type = ntype;
}

void MDIChild::setImported()
{
	if (!imported) 
	{
		Changed();
		imported = true;
	}
	onDisk = false;
}


int MDIChild::TextHeight(QFont F, QString T)
{
	//printf("Beginn TextHeight\n");
	int w, h;
	QPainter P;
	//F = QFont();
	Text(P, 0, 0, w, h, T, F, 0, false);
	//printf("Ende TextHeight\n");
	return h;
}

int MDIChild::TextWidth(QFont F, QString T)
{
	QPainter P;
	int w, h;
	Text(P, 0, 0, w, h, T, F, 0, false);
	return w;
}

void MDIChild::WriteText(QPainter &P, int x, int y, QString T, QFont Font, int Orientation)
{
	int w, h;
	Text(P, x, y, w, h, T, Font, Orientation, true);
}

void MDIChild::Text(QPainter &PD, int x, int y, int &w, int &h, QString T, const QFont F, int O,
							 bool write)
{
	//printf("Text: T=%s\n", T.ascii());
	w=h=0;
	if (T.length() == 0) return;
	int n=0, m=0, o=0, p=0, b, nw, nh, s=0, shift;
	double B;
	//F = QFont();
	QFontMetricsF FM(F);
	QFont SF = F;
	QString L;
	if ((n = T.indexOf("^")) != -1) s=-1;
	else if ((n = T.indexOf("_")) != -1) s=-2;
	else if ((n = T.indexOf("\\Alpha")) != -1) s=1;
	else if ((n = T.indexOf("\\Beta")) != -1) s=2;
	else if ((n = T.indexOf("\\Gamma")) != -1) s=3;
	else if ((n = T.indexOf("\\Delta")) != -1) s=4;
	else if ((n = T.indexOf("\\Epsilon")) != -1) s=5;
	else if ((n = T.indexOf("\\Zeta")) != -1) s=6;
	else if ((n = T.indexOf("\\Eta")) != -1) s=7;
	else if ((n = T.indexOf("\\Theta")) != -1) s=8;
	else if ((n = T.indexOf("\\Iota")) != -1) s=9;
	else if ((n = T.indexOf("\\Kappa")) != -1) s=10;
	else if ((n = T.indexOf("\\Lambda")) != -1) s=11;
	else if ((n = T.indexOf("\\Mu")) != -1) s=12;
	else if ((n = T.indexOf("\\Nu")) != -1) s=13;
	else if ((n = T.indexOf("\\Xi")) != -1) s=14;
	else if ((n = T.indexOf("\\Omicron")) != -1) s=15;
	else if ((n = T.indexOf("\\Pi")) != -1) s=16;
	else if ((n = T.indexOf("\\Rho")) != -1) s=17;
	else if ((n = T.indexOf("\\Sigma")) != -1) s=18;
	else if ((n = T.indexOf("\\Tau")) != -1) s=19;
	else if ((n = T.indexOf("\\Upsilon")) != -1) s=20;
	else if ((n = T.indexOf("\\Phi")) != -1) s=21;
	else if ((n = T.indexOf("\\Chi")) != -1) s=22;
	else if ((n = T.indexOf("\\Psi")) != -1) s=23;
	else if ((n = T.indexOf("\\Omega")) != -1) s=24;	
	else if ((n = T.indexOf("\\alpha")) != -1) s=25;
	else if ((n = T.indexOf("\\beta")) != -1) s=26;
	else if ((n = T.indexOf("\\gamma")) != -1) s=27;
	else if ((n = T.indexOf("\\delta")) != -1) s=28;
	else if ((n = T.indexOf("\\epsilon")) != -1) s=29;
	else if ((n = T.indexOf("\\zeta")) != -1) s=30;
	else if ((n = T.indexOf("\\eta")) != -1) s=31;
	else if ((n = T.indexOf("\\theta")) != -1) s=32;
	else if ((n = T.indexOf("\\iota")) != -1) s=33;
	else if ((n = T.indexOf("\\kappa")) != -1) s=34;
	else if ((n = T.indexOf("\\lambda")) != -1) s=35;
	else if ((n = T.indexOf("\\mu")) != -1) s=36;
	else if ((n = T.indexOf("\\nu")) != -1) s=37;
	else if ((n = T.indexOf("\\xi")) != -1) s=38;
	else if ((n = T.indexOf("\\omicron")) != -1) s=39;
	else if ((n = T.indexOf("\\pi")) != -1) s=40;
	else if ((n = T.indexOf("\\rho")) != -1) s=41;
	else if ((n = T.indexOf("\\sigma")) != -1) s=42;
	else if ((n = T.indexOf("\\tau")) != -1) s=43;
	else if ((n = T.indexOf("\\upsilon")) != -1) s=44;
	else if ((n = T.indexOf("\\phi")) != -1) s=45;
	else if ((n = T.indexOf("\\chi")) != -1) s=46;
	else if ((n = T.indexOf("\\psi")) != -1) s=47;
	else if ((n = T.indexOf("\\omega")) != -1) s=48;
	else if ((n = T.indexOf("\\rightarrow")) != -1) s=49;
	else if ((n = T.indexOf("\\leftarrow")) != -1) s=50;
	if (s!=0)
	{
		if (n>0) Text(PD, x, y, w, h, T.left(n), F, O, write);
		else w=h=0;
		if (s < 0 && (T.length() > n+2 ? T[n+1] == '{' : false))
		{
			m = T.indexOf('}', n+2);
			o = T.indexOf('{', m+1);
			//printf("n=%d, m=%d, o=%d\n", n, m, o);
			while (o<m && o!=-1 && m!=-1)
			{
				m = T.indexOf('}', m+1);
				o = T.indexOf('{', o+1);
			}
			//printf("m=%d, o=%d\n", m, o);
			if (m<0) m = T.length();
			o = n + 2;
			p = m + 1;
		}
		else
		{
            if (s<0)
			{
				o = n + 1;
				p = m = n + 2;
			}
			else switch (s)
			{
				case 1:
					p=m=(o=n)+6;
					break;
				case 2:
					p=m=(o=n)+5;
					break;
				case 3:
					p=m=(o=n)+6;
					break;
				case 4:
					p=m=(o=n)+6;
					break;
				case 5:
					p=m=(o=n)+8;
					break;
				case 6:
					p=m=(o=n)+5;
					break;
				case 7:
					p=m=(o=n)+4;
					break;
				case 8:
					p=m=(o=n)+6;
					break;
				case 9:
					p=m=(o=n)+5;
					break;
				case 10:
					p=m=(o=n)+6;
					break;
				case 11:
					p=m=(o=n)+7;
					break;
				case 12:
					p=m=(o=n)+3;
					break;
				case 13:
					p=m=(o=n)+3;
					break;
				case 14:
					p=m=(o=n)+3;
					break;
				case 15:
					p=m=(o=n)+8;
					break;
				case 16:
					p=m=(o=n)+3;
					break;
				case 17:
					p=m=(o=n)+4;
					break;
				case 18:
					p=m=(o=n)+6;
					break;
				case 19:
					p=m=(o=n)+4;
					break;
				case 20:
					p=m=(o=n)+8;
					break;
				case 21:
					p=m=(o=n)+4;
					break;
				case 22:
					p=m=(o=n)+4;
					break;
				case 23:
					p=m=(o=n)+4;
					break;
				case 24:
					p=m=(o=n)+6;
					break;
				case 25:
					p=m=(o=n)+6;
					break;
				case 26:
					p=m=(o=n)+5;
					break;
				case 27:
					p=m=(o=n)+6;
					break;
				case 28:
					p=m=(o=n)+6;
					break;
				case 29:
					p=m=(o=n)+8;
					break;
				case 30:
					p=m=(o=n)+5;
					break;
				case 31:
					p=m=(o=n)+4;
					break;
				case 32:
					p=m=(o=n)+6;
					break;
				case 33:
					p=m=(o=n)+5;
					break;
				case 34:
					p=m=(o=n)+6;
					break;
				case 35:
					p=m=(o=n)+7;
					break;
				case 36:
					p=m=(o=n)+3;
					break;
				case 37:
					p=m=(o=n)+3;
					break;
				case 38:
					p=m=(o=n)+3;
					break;
				case 39:
					p=m=(o=n)+8;
					break;
				case 40:
					p=m=(o=n)+3;
					break;
				case 41:
					p=m=(o=n)+4;
					break;
				case 42:
					p=m=(o=n)+6;
					break;
				case 43:
					p=m=(o=n)+4;
					break;
				case 44:
					p=m=(o=n)+8;
					break;
				case 45:
					p=m=(o=n)+4;
					break;
				case 46:
					p=m=(o=n)+4;
					break;
				case 47:
					p=m=(o=n)+4;
					break;
				case 48:
					p=m=(o=n)+6;
					break;
				case 49:
					p=m=(o=n)+11;
					break;
				case 50:
					p=m=(o=n)+10;
					break;
			}
		}
		if (n > 1000000) return;
		if (m>o)
		{
			if (s<0)
			{
				if ((b = F.pixelSize()) != -1) SF.setPixelSize(b = (3 * b) / 4);
				else SF.setPointSizeF(B = 0.75 * F.pointSizeF());
				if (s == -1) shift = - FM.ascent() / 2;
				else shift = FM.ascent() / 2;
				if (O==0) Text(PD, x+w, y + shift, nw, nh, T.mid(o, m-o), SF, O, write);
				else Text(PD, x + shift, y - w, nw, nh, T.mid(o, m-o), SF, O, write);
				w += nw;
				if ((b = FM.ascent() / 2 + nh) > h) h = b;
				//printf("T.length=%d, p=%d\n", T.length(), p);
			}
			else if (s == 49 || s == 50)
			{
				nw = 3 * (nh = TextHeight(F, T));
				nh /= 2;
				if (write)
				{
					if (O==0) 
					{
						PD.drawLine(x+w, y - nh, x + w + nw, y - nh);
						if (s==49)
						{
							PD.drawLine(x+w+5*nh, y-2*nh, x+w+nw, y-nh);
							PD.drawLine(x+w+5*nh, y, x+w+nw, y-nh);
						}
						else
						{
							PD.drawLine(x+w, y-nh, x+w+nh, y-2*nh);
							PD.drawLine(x+w, y-nh, x+w+nh, y);
						}
					}
					else 
					{
						PD.drawLine(x - nh, y - w, x - nh, y - w - nw);
						if (s==49)
						{
							PD.drawLine(x-2*nh, y-w-nw+nh, x-nh, y-w-nw);
							PD.drawLine(x, y-w-nw+nh, x-nh, y-w-nw);
						}
						else
						{
							PD.drawLine(x-nh, y-w, x-2*nh, y-w-nh);
							PD.drawLine(x-nh, y-w, x, y-w-nh);
						}
					}
				}
			}
			else 
			{
				switch (s)
				{
				case 1:
					L="A";
					break;
				case 2:
					L="B";
					break;
				case 3:
					L="G";
					break;
				case 4:
					L="D";
					break;
				case 5:
					L="E";
					break;
				case 6:
					L="Z";
					break;
				case 7:
					L="H";
					break;
				case 8:
					L=QChar(0x398);
					break;
				case 9:
					L="I";
					break;
				case 10:
					L="K";
					break;
				case 11:
					L="L";
					break;
				case 12:
					L="M";
					break;
				case 13:
					L="N";
					break;
				case 14:
					L="X";
					break;
				case 15:
					L="O";
					break;
				case 16:
					L="P";
					break;
				case 17:
					L="R";
					break;
				case 18:
					L="S";
					break;
				case 19:
					L="T";
					break;
				case 20:
					L="Y";
					break;
				case 21:
					L="F";
					break;
				case 22:
					L="C";
					break;
				case 23:
					L=QChar(0x3A8);
					break;
				case 24:
					L="W";
					break;
				case 25:
					L="a";
					break;
				case 26:
					L="b";
					break;
				case 27:
					L="g";
					break;
				case 28:
					L="d";
					break;
				case 29:
					L="e";
					break;
				case 30:
					L="z";
					break;
				case 31:
					L="h";
					break;
				case 32:
					L=QChar(0x3B8);
					break;
				case 33:
					L="i";
					break;
				case 34:
					L="k";
					break;
				case 35:
					L="l";
					break;
				case 36:
					L="m";
					break;
				case 37:
					L="n";
					break;
				case 38:
					L="x";
					break;
				case 39:
					L="o";
					break;
				case 40:
					L="p";
					break;
				case 41:
					L="r";
					break;
				case 42:
					L="s";
					break;
				case 43:
					L="t";
					break;
				case 44:
					L="y";
					break;
				case 45:
					L="f";
					break;
				case 46:
					L="c";
					break;
				case 47:
					L=QChar(0x3C8);
					break;
				case 48:
					L="w";
					break;
				}
				SF.setFamily("Symbol");
				if (O==0) Text(PD, x+w, y, nw, nh, L, SF, O, write);
				else Text(PD, x, y-w, nw, nh, L, SF, O, write);
				w += nw;
				if (nh > h) h = nh;
			}
		}
		if (p < T.length())
		{
			if (O==0) Text(PD, x+w, y, nw, nh, T.right(T.length() - p), F, O, write);
			else Text(PD, x, y-w, nw, nh, T.right(T.length() - p), F, O, write);
			w += nw;
			if (nh > h) h = nh;
		}
	}
	else
	{
		//printf("Vor write, h=%d\n", h);
		h = FM.height();
		//printf("Nach height\n");
		w = FM.width(T);
		//printf("Nach Metrics\n");
		if (write)
		{
			PD.setFont(F);
			if (O==1) 
			{
				PD.rotate(270);
				PD.drawText(-y, x, T);
				PD.rotate(90);
			}
			else PD.drawText(x, y, T);
		}
	}
	//printf("Text Ende\n");
}

