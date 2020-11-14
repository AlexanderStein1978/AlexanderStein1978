//
// C++ Interface: mdichild
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that comes with this source code
//
//

#ifndef MDICHILD_H
#define MDICHILD_H

#include <QWidget>
#include <QString>

class QFile;
class QPrinter;
class MainWindow;

/**
	@author Alexander Stein <AlexanderStein@t-online.de>
*/

class MDIChild : public QWidget
{
	Q_OBJECT	
public:
	typedef enum Typ {FitDataSet = -3, TextTable2, TextTable1, TermEnergyTable, DunhamTable, TermEnergyView, FranckCondonView, LineTab, 
		       FranckCondonTable, FitSeriesResultTable, SimpleMDIChild, SimpleDiagWindow, TermEnergyPlot, WaveFunctionPlot,
			   DataSetPlot, ResidPlot, AddSpect, SpectrumSimulation, AtomData, MolData, IntensityHistogramPlot, Pict, PotData, Spect,
               FCFJDependency, CouplingFunction, External} Type;
	
   	MDIChild(Type type = SimpleMDIChild, MainWindow *MW = 0,
			 QString filter = "files (*.*)", QString FileExt = ".dat");
   	virtual ~MDIChild();
	QString getFileName(bool save);
	void setFileName(QString Name);
	QString getFName();
	QString getName() const;
	void setType(Type type);
	Type getType();
	QString getFilter();
	bool isSaved();
	virtual bool isAssigned();
	void setMainWindow(MainWindow *MW);
	int TextHeight(QFont Font, QString Text);
	int TextWidth(QFont Font, QString Text);
	void WriteText(QPainter &P, int x, int y, QString Text, QFont Font, int Orientation);
	virtual void Print(QPrinter &P);
    QString getTypeString();
	virtual bool askForQuit();
    QString getAbsolutePath(QString &CurrentPath, QString &MolPath);
    QString getRelativePath(const QString &CurrentPath, const QString &MolPath);

    QString getRelativePath(const QString &MolPath)
    {
        return getRelativePath(FileName, MolPath);
    }
	
	inline bool isOnDisk() const
	{
		return onDisk;
	}
	
	inline QString getFileName() const
	{
		return FileName;
	}
	
	inline bool isImported() const
	{
		return imported;
	}

    inline bool didFileNameChange() const
    {
        return m_fileNameChanged;
    }
	
public slots:
	virtual bool writeData(QString Filename = "");
	virtual bool readData(QString Filename = "");
	virtual void setName(QString name);
	void show();
	void Changed();
    void Saved();

protected:
	void setImported();
	void setNewCreated();
	void setFilter(QString filter);
	void setFileExt(QString FileExt);
	bool read(QFile *File);
	bool write(QFile *File);
	void hideEvent(QHideEvent *E);
	void showEvent(QShowEvent *E);
	void setBlockChangeSignal(bool Block);
	bool getBlockChangeSignal();
    void closeEvent(QCloseEvent *i_event);
		
	MainWindow *MW;

signals:
	void hidden(bool);
	void nameChanged(QString);
	void fileNameChanged();
	void propertiesChanged();
private:
	void Text(QPainter &P, int x, int y, int &w, int &h, QString T, const QFont F, int O, 
					 bool write);
	
    bool m_fileNameChanged, blockChangeSignal;
	bool changed;
	bool imported;
	bool newCreated;
	bool onDisk;
	QString FileName;
	QString Name;
	QString Filter;
	Type type;
	QString FileExt;
};

#endif
