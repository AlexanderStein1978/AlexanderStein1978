//
// C++ Implementation: MainWindow
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2016
//
// Copyright: See README file that is delivered together with this source code
//
//

#include <QtGui>
#include <QTextStream>
#include <QFileDialog>
#include <QTextEdit>
#include <QMessageBox>
#include <QTableWidget>
#include <QStatusBar>
#include <QGridLayout>
#include <QComboBox>
#include <QMdiArea>
#include <QMenuBar>
#include <QMdiSubWindow>
#include <QApplication>
#include <QRadioButton>

#include <stdio.h>
#include <limits.h>

#include "MainWindow.h"
#include "molecule.h"
#include "termtable.h"
#include "duntable.h"
#include "fitdata.h"
#include "linetable.h"
#include "Spektrum.h"
#include "termview.h"
#include "TermPlot.h"
#include "potential.h"
#include "potentialplot.h"
#include "utils.h"
#include "aipsidiag.h"
#include "fcftable.h"
#include "wavefuncplot.h"
#include "spectsimulation.h"
#include "elstate.h"
#include "cutdialog.h"
#include "adddialog.h"
#include "dataplot.h"
#include "residualplot.h"
#include "intensityhistogram.h"
#include "fcfjdep.h"
#include "fcftab.h"
#include "montecarlosim.h"
#include "about.h"
#include "tableline.h"
#include "showaction.h"
#include "measuredtermenergies.h"
#include "spectlist.h"
#include "isotab.h"
#include "feldialog.h"
#include "CouplingFuncs.h"
#include "exportlineprofiledialog.h"
#include "SFQSCalcDialog.h"
#include "SFQSCalcControl.h"
#include "ImprovePotSeriesDialog.h"
#include "ImprovePotSeriesControl.h"
#include "CreateAnaPotSeriesFromMCSSplinePotSeriesDialog.h"
#include "MonteCarloSimDialog.h"
#include "MonteCarloSimControl.h"
#include "mcfsettingsdialog.h"
#include "sourceoffsetdialog.h"
#include "vassigndialog.h"
#include "assign_vs_comptrans.h"
#include "progression.h"
#include "fspdiag.h"
#include "coupledpotfitimportdialog.h"
#include "heapsort.h"
#include "tablelinesortfunctor.h"
#include "tlref.h"
#include "cttidialog.h"
#include "perturbation.h"
#include "CSWFImportDialog.h"
#include "CoupledSineWaveFunc.h"
#include "selectlinetabledialog.h"
#include "CreateAnaPotSeriesControl.h"


MainWindow::MainWindow()
{
	//printf("Begin MainWindow::MainWindow\n");
	numAtoms = 0;
	numMolecules = 0;
	numTermTables = 0;
	numMeasuredTermTables = 0;
	numDunTables = 0;
	numPotentials = 0;
	numFitDataSets = 0;
	numLineTables = 0;
	numFCFTables = 0;
	numSpectra = 0;
	numSpectLists = 0;
	numSpectSimulations = 0;
	numNAAtoms = numNAMolecules = numNATermTables = numNADunTables = numNAPotentials = 0;
	numNALineTables = numNAFCFTables = numNASpectra = numNAFitDataSets = 0;
	nCR = nCC = 0;
	cuttedRows = 0;
	NumWindowsToAskForQuit = 0;
	windowsToAskForQuit = 0;
	
	MaxDirRep = 100;
	numDirRep = 0;
	DirRepNC = new int[MaxDirRep];
	DirRep = new QString[MaxDirRep];
		
	setWindowTitle(QString("MolSpektAnalysis ") + MAVersion);
	
	Dir = QDir::currentPath();
	
	workspace = new QMdiArea;
	setCentralWidget(workspace);
	connect(workspace, SIGNAL(subWindowActivated(QMdiSubWindow*)), 
			this, SLOT(WindowActivated(QMdiSubWindow*)));
	statusBar();
	
	newAtomAct = new QAction(tr("&Atom"), this);
	//newAtomAct->setShortcut(tr("Ctrl+A"));
	newAtomAct->setStatusTip(tr("Create a window for the data of a new atom"));
	connect(newAtomAct, SIGNAL(triggered()), this, SLOT(newAtom()));
	
	newMoleculeAct = new QAction(tr("&Molecule"), this);
	//newMoleculeAct->setShortcut(tr("Ctrl+M"));
	newMoleculeAct->setStatusTip("Create a window for the data of a new molecule");
	connect(newMoleculeAct, SIGNAL(triggered()), this, SLOT(newMolecule()));
	
	newTermTableAct = new QAction("&Term energy table", this);
	//newTermTableAct->setShortcut(tr("Ctrl+T"));
	newTermTableAct->setStatusTip("Create a window for a table of term energies");
	connect(newTermTableAct, SIGNAL(triggered()), this, SLOT(newTermTable()));
	
	newMeasuredTermTableAct = new QAction("Measured term &energy table", this);
	newMeasuredTermTableAct->setStatusTip("Create a window for a table of measured term energies");
	connect(newMeasuredTermTableAct, SIGNAL(triggered()), this, SLOT(newMeasuredTermTable()));
	
	newDunTableAct = new QAction("&Dunham coefficients", this);
	//newDunTableAct->setShortcut(tr("Ctrl+D"));
	newDunTableAct->setStatusTip("Create a window for a table of Dunham coefficients");
	connect(newDunTableAct, SIGNAL(triggered()), this, SLOT(newDunTable()));
	
	newPotentialAct = new QAction("&Potential", this);
	//newPotentialAct->setShortcut(tr("Ctrl+P"));
	newPotentialAct->setStatusTip("Create a table window for a potenital");
	connect(newPotentialAct, SIGNAL(triggered()), this, SLOT(newPotential()));
	
	newFitDataAct = new QAction("F&it dataset", this);
	newFitDataAct->setStatusTip("Create a table window for a feit dataset");
	connect(newFitDataAct, SIGNAL(triggered()), this, SLOT(newFitData()));
	
	newLineTableAct = new QAction("&Line table", this);
	//newLineTableAct->setShortcut(tr("Ctrl+L"));
	newLineTableAct->setStatusTip("Create a window for a table of line tables");
	connect(newLineTableAct, SIGNAL(triggered()), this, SLOT(newLineTable()));
	
	newFCFTableAct = new QAction("&FCF table...", this);
	newFCFTableAct->setStatusTip("Calculate a new table of Franck-Condon factors");
	newFCFTableAct->setEnabled(false);
	connect(newFCFTableAct, SIGNAL(triggered()), this, SLOT(newFCFTable()));
	
	newSpectrumAct = new QAction("&Spectrum", this);
	//newSpectrumAct->setShortcut(tr("Shift+Ctrl+S"));
	newSpectrumAct->setStatusTip("Create a window for a spectrum");
	connect(newSpectrumAct, SIGNAL(triggered()), this, SLOT(newSpectrum()));
	
	newSpectListAct = new QAction("Spe&ctList", this);
	newSpectListAct->setStatusTip("Create a list for spectra");
	connect(newSpectListAct, SIGNAL(triggered()), this, SLOT(newSpectList()));
	
	openAct = new QAction("&Open...", this);
	openAct->setShortcut(tr("Ctrl+O"));
	openAct->setStatusTip("Open an atom file, a molecule file or a spectra...");
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));
	
	reloadAct = new QAction("&Reload file", this);
	reloadAct->setShortcut(tr("Ctrl+D"));
	reloadAct->setStatusTip("Reload the file displayed in the current window to discard all changes done since its last saving");
	connect(reloadAct, SIGNAL(triggered()), this, SLOT(reload()));
	
	saveAct = new QAction("&Save", this);
	saveAct->setShortcut(tr("Ctrl+S"));
	saveAct->setStatusTip("Save the contents of the active window");
	saveAct->setEnabled(false);
	connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));
	
	saveAsAct = new QAction("Save &As...", this);
	saveAsAct->setStatusTip("Save the contents of the active window under a new name");
	saveAsAct->setEnabled(false);
	connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));
	
	writePotDataAct = new QAction("&Write potential data...", this);
	writePotDataAct->setStatusTip(
					"Write an analytical potential file as it is used in the group of prof. Tiemann");
	connect(writePotDataAct, SIGNAL(triggered()), this, SLOT(writePotData()));
	
	writeTFGSAct = new QAction("Write &TFGS...", this);
	writeTFGSAct->setStatusTip("Write the line table as an input file for potential fitting as it is used in the group of prof. Tiemann");
	connect(writeTFGSAct, SIGNAL(triggered()), this, SLOT(writeTFGS()));
	
	writeExpPotFitInputAct = new QAction("Write &ES PotFit Input...", this);
	writeExpPotFitInputAct->setStatusTip("Write the line table as an input file for fitting a potential for the excited state as it is used in the group of prof. Tiemann");
	connect(writeExpPotFitInputAct, SIGNAL(triggered()), this, SLOT(writeExpPotFitInput()));
	
	write2AtInputAct = new QAction("Write &ZweiAt input.dat...", this);
	write2AtInputAct->setStatusTip("Write an input file for the program ZweiAt of the group of prof. Tiemann");
	connect(write2AtInputAct, SIGNAL(triggered()), this, SLOT(write2AtInput())); 
	
	exportTFDunAct = new QAction("E&xport to TF...", this);
	exportTFDunAct->setStatusTip("Export the table of Dunham coefficients to a format as it is used in the group of prof. Tiemann");
	connect(exportTFDunAct, SIGNAL(triggered()), this, SLOT(exportTFDun()));
	
	exportPotPointsAct = new QAction("Export potential p&oints...", this);
	exportPotPointsAct->setStatusTip("Export potential as a list of points");
	connect(exportPotPointsAct, SIGNAL(triggered()), this, SLOT(exportPotPoints()));
	
	exportAsymptoticLevelsAct = new QAction("Export &asymptotic levels...", this);
	exportAsymptoticLevelsAct->setStatusTip("Calculate asymptotic levels from a potential an export it to a text file");
	connect(exportAsymptoticLevelsAct, SIGNAL(triggered()), this, SLOT(exportAsymptoticLevels()));
	
	exportPictureAct = new QAction("Export &picture...", this);
	exportPictureAct->setStatusTip("Export the window content as a picture");
	connect(exportPictureAct, SIGNAL(triggered()), this, SLOT(exportPicture()));
	
	exportTableDataAct = new QAction("Expo&rt table data...", this);
	exportTableDataAct->setStatusTip("Export data from the table to a text file");
	connect(exportTableDataAct, SIGNAL(triggered()), this, SLOT(exportTableData()));
	
	exportWaveFunctionAct = new QAction("Export wa&vefunction...", this);
	exportWaveFunctionAct->setStatusTip("Export a wavefunction of a coupled system to a text file");
	connect(exportWaveFunctionAct, SIGNAL(triggered()), this, SLOT(exportWaveFunction()));

    exportLineProfileAct = new QAction("Export fitted line profile...", this);
    exportLineProfileAct->setStatusTip("Export a fitted line profile as point table in a text file.");
    connect(exportLineProfileAct, SIGNAL(triggered()), this, SLOT(exportLineProfile()));
	
	exportObservedLevelListsAct = new QAction("Export ob&served levels...", this);
	exportAsymptoticLevelsAct->setStatusTip(
		"Export lists of the quantum numbers of the observed levels in files seperated for electronic states, isotopologues and e/f components");
	connect(exportObservedLevelListsAct, SIGNAL(triggered()), 
			this, SLOT(exportObservedLevelLists()));
	
	importTermTableAct = new QAction("&Term energy table...", this);
	importTermTableAct->setStatusTip("Import term energy table");
	connect(importTermTableAct, SIGNAL(triggered()), this, SLOT(importTermTable()));
	
	importDunTableAct = new QAction("&Dunham coefficients...", this);
	importDunTableAct->setStatusTip("Import table with Dunham coefficients");
	connect(importDunTableAct, SIGNAL(triggered()), this, SLOT(importDunTable()));
	
	importPotentialAct = new QAction("&Potential...", this);
	importPotentialAct->setStatusTip("Import a potential");
	connect(importPotentialAct, SIGNAL(triggered()), this, SLOT(importPotential()));
	
	importFitDataAct = new QAction("&Fit data...", this);
	importFitDataAct->setStatusTip("Import a table of data for a potenital fit");
	connect(importFitDataAct, SIGNAL(triggered()), this, SLOT(importFitData()));
	
	importLineTableAct = new QAction("&Line table...", this);
	importLineTableAct->setStatusTip("Import a line table");
	connect(importLineTableAct, SIGNAL(triggered()), this, SLOT(importLineTable()));
	
	importSpectrumAct = new QAction("&Spectrum...", this);
	importSpectrumAct->setStatusTip("Import a spectrum");
	connect(importSpectrumAct, SIGNAL(triggered()), this, SLOT(importSpectrum()));
	
	importSpectListAct = new QAction("Sp&ectList...", this);
	importSpectListAct->setStatusTip("Import a list of spectra");
	connect(importSpectListAct, SIGNAL(triggered()), this, SLOT(importSpectList()));
	
	read2AtOutputAct = new QAction("&ZweiAt output...", this);
	read2AtOutputAct->setStatusTip(
							"Import a ZweiAt output file as it is used in the group of prof. Tiemann");
	connect(read2AtOutputAct, SIGNAL(triggered()), this, SLOT(read2AtOutput()));
	
	importAbInitioPotSetAct = new QAction("&Ab initio potentials...", this);
	importAbInitioPotSetAct->setStatusTip(
										"Import a whole set of ab initio potentials at the same time");
	connect(importAbInitioPotSetAct, SIGNAL(triggered()), this, SLOT(importAbInitioPotSet()));
	
	importMTTPotentialsAct = new QAction("&Import MTT potentials...", this);
	importMTTPotentialsAct->setStatusTip("Import modified Tang Tönnies potentials as used in arXiv:1203.4524v1");
	connect(importMTTPotentialsAct, SIGNAL(triggered()), this, SLOT(importMTTPotentials()));
	
	importDunhamAsenAct = new QAction("Dunham &coefficients from Asen...", this);
	importDunhamAsenAct->setStatusTip(
					"Import a Dunham type coefficient set in the format as received from Asen Pashov");
	connect(importDunhamAsenAct, SIGNAL(triggered()), this, SLOT(importDunhamAsen()));
	
	importCoupledPotfitOutputAct = new QAction("Coupled potfit &output...", this);
	importCoupledPotfitOutputAct->setStatusTip(
		"Import an output file of a potential fit program for coupled states as it is used in the group of prof. Tiemann");
	connect(importCoupledPotfitOutputAct, SIGNAL(triggered()), 
			this, SLOT(importCoupledPotfitOutput()));
	
	importCoupledTermTableAct = new QAction("Co&upled term energy table...", this);
	importCoupledTermTableAct->setStatusTip(
		"Import an output file with term energy data of a program for potential fits of coupled systems.");
	connect(importCoupledTermTableAct, SIGNAL(triggered()), 
			this, SLOT(importCoupledTermTable()));
	
	importCoupledWaveFunctionsAct = new QAction("Coupled &wave functions...", this);
	importCoupledWaveFunctionsAct->setStatusTip("Import file with data of wave functions from coupled channel calculations");
	connect(importCoupledWaveFunctionsAct, SIGNAL(triggered()), this, SLOT(importCoupledWaveFunctions()));
	
	printAct = new QAction("&Print...", this);
	printAct->setStatusTip("Prints the contents of the active window if possible");
	connect(printAct, SIGNAL(triggered()), this, SLOT(print()));
		
	exitAct = new QAction(tr("&Quit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("End the program"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(quit()));
	
	addRowAct = new QAction("&Add row", this);
	addRowAct->setShortcut(tr("Ctrl+N"));
	addRowAct->setStatusTip("Add a new row to the table");
	connect(addRowAct, SIGNAL(triggered()), this, SLOT(addRow()));
	
	deleteRowsAct = new QAction("&Delete rows", this);
	deleteRowsAct->setShortcut(tr("Ctrl+Del"));
	deleteRowsAct->setStatusTip("Delete complete rows from the table");
	connect(deleteRowsAct, SIGNAL(triggered()), this, SLOT(deleteRows()));
	
	copyRowsAct = new QAction("&Copy rows", this);
	copyRowsAct->setStatusTip("Copy complete rows from the table");
	connect(copyRowsAct, SIGNAL(triggered()), this, SLOT(copyRows()));
	
	cutRowsAct = new QAction("Cu&t rows", this);
	cutRowsAct->setShortcut(tr("Ctrl+X"));
	cutRowsAct->setStatusTip("Cut comlete rows from the table");
	connect(cutRowsAct, SIGNAL(triggered()), this, SLOT(cutRows()));
	
	insertRowsAct = new QAction("&Insert rows", this);
	insertRowsAct->setShortcut(tr("Ctrl+V"));
	insertRowsAct->setStatusTip(
		"Insert the content of rows previously cutted from a table as new lines to the selected table");
	connect(insertRowsAct, SIGNAL(triggered()), this, SLOT(insertRows()));	
	
	searchAct = new QAction("&Search...", this);
	searchAct->setStatusTip("Search for a text or a value larger, smaller or equal a given number in the current table");
	connect(searchAct, SIGNAL(triggered()), this, SLOT(search()));
	
	shiftValuesAct = new QAction("Sh&ift values...", this);
	shiftValuesAct->setStatusTip("Shift the values of the selected column by a value to enter");
	connect(shiftValuesAct, SIGNAL(triggered()), this, SLOT(shiftValues()));

    setTextInCellsAct = new QAction("Set text in cells...", this);
    setTextInCellsAct->setStatusTip("Type a text to be set into the selected columns");
    connect(setTextInCellsAct, SIGNAL(triggered()), this, SLOT(setTextInCells()));
	
	updateAct = new QAction("&Update", this);
	updateAct->setStatusTip(
		"Update the current table or its fitdata if a fit is necessary for a update of the table itself");
	connect(updateAct, SIGNAL(triggered()), this, SLOT(update()));
	
	findWrongDataAct = new QAction("&Find bad levels", this);
	findWrongDataAct->setStatusTip(
		"Searches the table for levels with obs-calc >= 3.5 or whose J values not exist or with J'-J''= 0 for sigma states");
	connect(findWrongDataAct, SIGNAL(triggered()), this, SLOT(findWrongData()));
	
	addCalculatedLevelsAct = new QAction("Add calculated &levels...", this);
	addCalculatedLevelsAct->setStatusTip("Add calculated energy levels to the table");
	connect(addCalculatedLevelsAct, SIGNAL(triggered()), this, SLOT(addCalculatedLevels()));
	
	showNumLevelsAct = new QAction("S&how level numbers", this);
	showNumLevelsAct->setStatusTip(
		"Show the number of levels in the current table subdivided into component and isotopologue");
	connect(showNumLevelsAct, SIGNAL(triggered()), this, SLOT(showNumLevels()));
	
	showUncertaintyStatsAct = new QAction("Show u&ncertainty stats", this);
	showUncertaintyStatsAct->setStatusTip("Shows how many levels have which uncertainties");
	connect(showUncertaintyStatsAct, SIGNAL(triggered()), this, SLOT(showUncertaintyStats()));
	
	changeSourceOffsetAct = new QAction("Change source &offset...", this);
	changeSourceOffsetAct->setStatusTip("View and change energy offsets given for the different data sources");
	connect(changeSourceOffsetAct, SIGNAL(triggered()), this, SLOT(changeSourceOffset()));
	
	compareLevelEnergiesAct = new QAction("Compare &level energies...", this);
	compareLevelEnergiesAct->setStatusTip(
		"Compare the level energies of the current table with the energies of another term energy table, Dunham coefficient set or potential and give the average energy offset and the remaining standard deviation");
	connect(compareLevelEnergiesAct, SIGNAL(triggered()), this, SLOT(compareLevelEnergies()));
	
	fixJOffsetsAct = new QAction("Fix &J offsets", this);
	fixJOffsetsAct->setStatusTip("Use the current fitdata set to create a new fitdata set where the energy offsets of the progressions are fixed to values calculated with the current term energy tabe of the electronic state.");
	connect(fixJOffsetsAct, SIGNAL(triggered()), this, SLOT(fixJOffsets()));
	
	removeSingleLinesAct = new QAction("Remove sing&gle lines", this);
	removeSingleLinesAct->setStatusTip("Remove all progressions from the table which consist only of a single line");
	connect(removeSingleLinesAct, SIGNAL(triggered()), this, SLOT(removeSingleLines()));
	
	removeDataFSourceAct = new QAction("&Remove data from source...", this);
	removeDataFSourceAct->setStatusTip("Deletes the lines from the table which stem from the selected source");
	connect(removeDataFSourceAct, SIGNAL(triggered()), this, SLOT(removeDataFSource()));
	
	selectDataFSourceAct = new QAction("S&elect data from source...", this);
	selectDataFSourceAct->setStatusTip("Selects the lines in the table which stem from the choosen source");
	connect(selectDataFSourceAct, SIGNAL(triggered()), this, SLOT(selectDataFSource()));
	
	showIsotopologuesAct = new QAction("&Show isotopologues", this);
	//showIsotopomersAct->setShortcut(tr("Ctrl+I"));
	showIsotopologuesAct->setStatusTip("Show a table of the isotopologues of this molecule");
	connect(showIsotopologuesAct, SIGNAL(triggered()), this, SLOT(showIsotopologues()));
	
	showLevelNumbersAct = new QAction("Show &level numbers", this);
	showLevelNumbersAct->setStatusTip(
							"Show the number of known levels subdivided into states and isotopomers");
	connect(showLevelNumbersAct, SIGNAL(triggered()), this, SLOT(showLevelNumbers()));
	
	extractDataWithComponentAct = new QAction("&Extract data from table...", this);
	extractDataWithComponentAct->setStatusTip(
		"Extract the lines/levels with e or f symmetry, of a special fine structure component or a special level range to a new fit data set");
	connect(extractDataWithComponentAct, SIGNAL(triggered()), this, SLOT(extractDataWithComponent()));
    
    extractChangedDataAct = new QAction("Extract data changed compared to fit data...", this);
    extractChangedDataAct->setStatusTip("Extract the lines/levels which are changed compared to a fit data set to select.");
    connect(extractChangedDataAct, SIGNAL(triggered()), this, SLOT(extractChangedData()));
    
    extractNewDataAct = new QAction("Extract data new compared to fit data...", this);
    extractNewDataAct->setStatusTip("Extract this lines/levels which are new compared to a fit data set to select.");
    connect(extractNewDataAct, SIGNAL(triggered()), this, SLOT(extractNewData()));
	
	calcFCFTableAct = new QAction("Calculate &FCF table...", this);
	calcFCFTableAct->setStatusTip("Calculate Franck-Condon factors to save in a FCF table");
	connect(calcFCFTableAct, SIGNAL(triggered()), this, SLOT(calcFCFTable()));

    checkAllConnectionsAct = new QAction("&Check all connections", this);
    checkAllConnectionsAct->setStatusTip(
        "Verify that all links of this molecule to term energy tables, tables of Dunham coefficients, potentials, fit data sets, line tables and FCF tables and the links of fit data sets and line tables to spectra are working.");
    connect(checkAllConnectionsAct, SIGNAL(triggered()), this, SLOT(checkAllConnections()));

    shrinkAllSpektRefsAct = new QAction("S&rink all spect refs", this);
    shrinkAllSpektRefsAct->setStatusTip(
        "Reduce the references of the lines in the fit data sets and line tables of this molecule, which contain usually the full path of the file containing the spectrum they stem from to the pure file name.");
    connect(shrinkAllSpektRefsAct, SIGNAL(triggered()), this, SLOT(shrinkAllSpectRefs()));

	addDunTableLineAct = new QAction("&Add line", this);
	addDunTableLineAct->setStatusTip("Add a line to the table of Dunham coefficients");
	connect(addDunTableLineAct, SIGNAL(triggered()), this, SLOT(addDunTableLine()));
	
	calcTermEnergiesAct = new QAction("&Calculate term energies", this);
	//calcTermEnergiesAct->setShortcut(tr("Ctrl+C"));
	calcTermEnergiesAct->setStatusTip(
								"Calculate term energies using the actual set of Dunham coefficients");
	connect(calcTermEnergiesAct, SIGNAL(triggered()), this, SLOT(calcTermEnergies()));
	
	updateDunhamAct = new QAction("&Update coefficients", this);
	updateDunhamAct->setShortcut(tr("Ctrl+U"));
	updateDunhamAct->setStatusTip(
		"Update the coefficients of the set by using the data from the available line tables assigned to the same electronic state");
	connect(updateDunhamAct, SIGNAL(triggered()), this, SLOT(UpdateDunham()));
	
	improveDunhamAct = new QAction("&Improve coefficient set", this);
	improveDunhamAct->setStatusTip(
							"Improve the Dunham coefficient set by adjusting the coefficient number");
	connect(improveDunhamAct, SIGNAL(triggered()), this, SLOT(ImproveDunham()));

	fitBorderLineAct = new QAction("&Fit border line...", this);
	fitBorderLineAct->setStatusTip(
		"See and/or change the points through which the border line gets drawn below which the available levels get used for the fit.");
	connect(fitBorderLineAct, SIGNAL(triggered()), this, SLOT(fitBorderLine()));
	
	removeFineStructureAct = new QAction("&Remove fine structure", this);
	removeFineStructureAct->setStatusTip(
		"Use the spin rotation coefficients to create a new fitdata set with removed fine structure");
	connect(removeFineStructureAct, SIGNAL(triggered()), 
			this, SLOT(removeFineStructure()));
	
	removeLambdaDoublingAct = new QAction("Remove &lambda doubling", this);
	removeLambdaDoublingAct->setStatusTip(
		"Use the lambda doubling coefficients to construct a new fitdata set with removed lambda doubling");
	connect(removeLambdaDoublingAct, SIGNAL(triggered()), this, SLOT(removeLambdaDoubling()));
	
	applyAdCorrAct = new QAction("A&pply adiabatic correction", this);
	applyAdCorrAct->setStatusTip(
		"Construct a new fitdata set with applied adiabatic correction to remove an energetic splitting between different isotopologues");
	connect(applyAdCorrAct, SIGNAL(triggered()), this, SLOT(applyAdCorr()));
	
	testKratzerAct = new QAction("&Test Kratzers relation", this);
	testKratzerAct->setStatusTip("Use the Dunham coefficients to test the Kratzer relation");
	connect(testKratzerAct, SIGNAL(triggered()), this, SLOT(testKratzer()));
	
	calcReAct = new QAction("Calculate R&e", this);
	calcReAct->setStatusTip("Use the rotational constant to calculate R_e");
	connect(calcReAct, SIGNAL(triggered()), this, SLOT(calcRe()));
	
	calcY00TeAct = new QAction("Calculate &Y00 and Te", this);
	calcY00TeAct->setStatusTip("Use the Dunham coefficients to calculate Y_00 and T_e");
	connect(calcY00TeAct, SIGNAL(triggered()), this, SLOT(calcY00Te()));
	
	showBandConstantsAct = new QAction("&Show band constants", this);
	showBandConstantsAct->setStatusTip(
		"Show a table of the band constants for the different vibrational bands calculated from the Dunham coefficients");
	connect(showBandConstantsAct, SIGNAL(triggered()), this, SLOT(showBandConstants()));
	
	calcTermEnergiesPotAct = new QAction("&Calculate term energies", this);
	calcTermEnergiesPotAct->setStatusTip(
			"Calculate all term energies by solving the radial Schrödinger equation for this potenital");
	connect(calcTermEnergiesPotAct, SIGNAL(triggered()), this, SLOT(calcTermEnergiesPot()));
	
	calcScatWaveFuncAct = new QAction("Calc. &Scattering WF", this);
	calcScatWaveFuncAct->setStatusTip("Calculate the scattering wave function and write it to a file");
	connect(calcScatWaveFuncAct, SIGNAL(triggered()), this, SLOT(calcScatWaveFunc()));
	
	showScatNodePosAct = new QAction("Show scatNodePositions", this);
	showScatNodePosAct->setStatusTip(
		"Show the node positions of the s-wave scattering wave functions for the different isotopomers");
	connect(showScatNodePosAct, SIGNAL(triggered()), this, SLOT(showScatNodePos()));
	
	showTeReAct = new QAction("Give Re and Te", this);
	showTeReAct->setStatusTip("Gives Re and Te of the current potential");
	connect(showTeReAct, SIGNAL(triggered()), this, SLOT(showTeRe()));

    calcRoughBeAct = new QAction("Calc roughly Be", this);
    calcRoughBeAct->setStatusTip("Calculates the rotational constant Be simply from the internuclear distance of the minimum of the potential function which is a fast but not very precise method.");
    connect(calcRoughBeAct, SIGNAL(triggered()), this, SLOT(calcRoughBe()));
	
	fitScatNodePosAct = new QAction("Fit scatNodePositions", this);
	fitScatNodePosAct->setStatusTip(
		"Fit a certain node position of the s-wave scattering wave function by varying a single potential parameter");
	connect(fitScatNodePosAct, SIGNAL(triggered()), this, SLOT(fitScatNodePos()));
	
	autoCalcScatLengthsPotentialSetAct = new QAction("Calc scattering lengths...", this);
	autoCalcScatLengthsPotentialSetAct->setStatusTip(
		"Use the program 'scat_all_asym9.exe' to calculate scattering lengths for a whole set of potentials and write it to a file");
	connect(autoCalcScatLengthsPotentialSetAct, SIGNAL(triggered()), 
			this, SLOT(autoCalcScatLengthsPotentialSet()));
	
	createAnaPotsFromMCSSplineSeriesAct = new QAction("Create AnaPot series...", this);
	createAnaPotsFromMCSSplineSeriesAct->setStatusTip(
		"Creates analytical potentials from a larger number of spline potentials collected in a single directory. Usefull e.g. for the application on results of Monte Carlo simulations.");
	connect(createAnaPotsFromMCSSplineSeriesAct, SIGNAL(triggered()), this, SLOT(createAnaPotsFromMCSSplineSeries()));
	
	fitControlAct = new QAction("Run potential fit series...", this);
	fitControlAct->setStatusTip(
		"Control an external fitting program to estimate the uncertainty of special fit parameters or fit a parameter which cannot be directly fitted by the external program");
	connect(fitControlAct, SIGNAL(triggered()), this, SLOT(fitControl()));
	
	monteCarloSimAct = new QAction("Monte Carlo simulation...", this);
	monteCarloSimAct->setStatusTip(
		"Run a series of Monte Carlo simulations with randomly varied measurement data to estimate the uncertainties of coefficients of the current potential");
	connect(monteCarloSimAct, SIGNAL(triggered()), this, SLOT(monteCarloSim()));
	
	summarizePotInfoAct = new QAction("Summarize potential info...", this);
	summarizePotInfoAct->setStatusTip("Combine the data of all potentials in a selectable directory in one table.");
	connect(summarizePotInfoAct, SIGNAL(triggered()), this, SLOT(summarizePotInfo()));
	
	calcSFQSMCSSeriesAct = new QAction("Calc. SFQS series...", this);
	calcSFQSMCSSeriesAct->setStatusTip(
		"Calculate chiSq of all levels and special chiSq of levels with outer turing points larger than a specified core radius for all potentials created by a MCS fit series in a directory to choose.");
	connect(calcSFQSMCSSeriesAct, SIGNAL(triggered()), this, SLOT(calcSFQSMCSSeries()));
	
	improvePotSeriesAct = new QAction("Improve pot series...", this);
	improvePotSeriesAct->setStatusTip("Fit a whole series of potentials collected in the same directory to the same fit data set");
	connect(improvePotSeriesAct, SIGNAL(triggered()), this, SLOT(improvePotSeries()));
	
	plotPotentialAct = new QAction("Plot potential", this);
	plotPotentialAct->setStatusTip("Show a plot of the current potential");
	connect(plotPotentialAct, SIGNAL(triggered()), this, SLOT(plotPotential()));
	
	showFQSAct = new QAction("Show chiSq", this);
	showFQSAct->setStatusTip(
		"Show the sum of the squared obs-calc values for all line tables connected via the assigned electronic state and a transition");
	connect(showFQSAct, SIGNAL(triggered()), this, SLOT(showFQS()));
	
	showBadListAct = new QAction("Show badlist", this);
	showBadListAct->setStatusTip(
		"Show a list of line frequencies and term energies, which do not fit sufficiently well to the current potential");
	connect(showBadListAct, SIGNAL(triggered()), this, SLOT(showBadList()));
	
	showFitDataAct = new QAction("Show fit data", this);
	showFitDataAct->setStatusTip("Show the data for the potential fit");
	connect(showFitDataAct, SIGNAL(triggered()), this, SLOT(showFitData()));
	
	updateFitDataAct = new QAction("Update fit data", this);
	updateFitDataAct->setStatusTip("Update the fit data for the potenital");
	connect(updateFitDataAct, SIGNAL(triggered()), this, SLOT(updateFitData()));
	
	fitSplinePotAct = new QAction("Fit spline potential", this);
	fitSplinePotAct->setStatusTip("Do a singlie fitting iteration to fit the spline potential to the data available for the connected electronic state");
	connect(fitSplinePotAct, SIGNAL(triggered()), this, SLOT(fitSplinePot()));
	
	fitAnaPotAct = new QAction("Fit analytical potential", this);
	fitAnaPotAct->setStatusTip("Use the current potential as initial curve to fit an analytical potential using the description of the group of prof. Tiemann or fit the current analytical potential to the available experimental data");
	connect(fitAnaPotAct, SIGNAL(triggered()), this, SLOT(fitAnaPot())); 
	
	fitTangToenniesPotAct = new QAction("Fit Tang-Toennies potential", this);
	fitTangToenniesPotAct->setStatusTip("Fit a Tang-Toennies potential");
	connect(fitTangToenniesPotAct, SIGNAL(triggered()), this, SLOT(fitTangToenniesPot()));
	
	fitWRobustWeightingAct = new QAction("Fit with robust weighting", this);
	fitWRobustWeightingAct->setStatusTip("Fit analytical potentials using a robust weighting algorithm");
	connect(fitWRobustWeightingAct, SIGNAL(triggered()), this, SLOT(fitWRobustWeighting()));
	
	fitIsoMassAct = new QAction("Fit iso mass", this);
	fitIsoMassAct->setStatusTip("Fit the mass of an isotopologue to the available levels/lines");
	connect(fitIsoMassAct, SIGNAL(triggered()), this, SLOT(fitIsoMass()));
	
	showSFuncsAct = new QAction("Show S functions", this);
	showSFuncsAct->setStatusTip("Plot the S funtions needed for fitting potentials using the IPA method");
	connect(showSFuncsAct, SIGNAL(triggered()), this, SLOT(showSFuncs()));
	
	setEnergyOffsetAct = new QAction("Set energy offset...", this);
	setEnergyOffsetAct->setStatusTip("Set the energy offset of the potential by specifying a certain offset value for the minimum or the asymptote");
	connect(setEnergyOffsetAct, SIGNAL(triggered()), this, SLOT(setEnergyOffset()));
	
	selPotCoeffAct = new QAction("Select coefficients...", this);
	selPotCoeffAct->setStatusTip(
							"Shows a dialog to change the coefficients used for the potential description");
	connect(selPotCoeffAct, SIGNAL(triggered()), this, SLOT(selPotCoeff()));
	
	scalePotentialAct = new QAction("Scale potential...", this);
	scalePotentialAct->setStatusTip(
			"Produce a new potential by scaling the current potentials to selectable Re and De values");
	connect(scalePotentialAct, SIGNAL(triggered()), this, SLOT(scalePotential()));
	
	cdConnectSRAct = new QAction("Connect S&R", this);
	cdConnectSRAct->setStatusTip(
			"Calculate a continuously differentiable short range extension to the current potential");
	connect(cdConnectSRAct, SIGNAL(triggered()), this, SLOT(cdConnectSR()));
	
	cdConnectLRAct = new QAction("Connect &LR...", this);
	cdConnectLRAct->setStatusTip(
			"Calculate a continuously differentiable long range extension to the current potential");
	connect(cdConnectLRAct, SIGNAL(triggered()), this, SLOT(cdConnectLR()));

	calcyssAct = new QAction("Recalculate y''", this);
	calcyssAct->setStatusTip(
		"Recalculate the second deviations of the points of the spline potential");
	connect(calcyssAct, SIGNAL(triggered()), this, SLOT(calcyss()));
	
	fixCoefficientsAct = new QAction("&fix coefficients", this);
	fixCoefficientsAct->setStatusTip(
		"Set the selected potential coefficients to be fixed in potential fits");
	connect(fixCoefficientsAct, SIGNAL(triggered()), this, SLOT(fixCoefficients()));
	
	varyCoefficientsAct = new QAction("&vary coefficients", this);
	varyCoefficientsAct->setStatusTip(
		"Set the selected potential coefficients to be variable in potential fits");
	connect(varyCoefficientsAct, SIGNAL(triggered()), this, SLOT(varyCoefficients()));
	
	createMLRPotAct = new QAction("Create MLR Potential", this);
	createMLRPotAct->setStatusTip("Calculate a MLR potential from the data of selected potential");
	connect(createMLRPotAct, SIGNAL(triggered()), this, SLOT(createMLRPot()));
	
	plotShowMouseCrossAct = new QAction("Show &mouse cross", this);
	plotShowMouseCrossAct->setStatusTip("Show a cross in the diagram connected to the mouse arrow");
	plotShowMouseCrossAct->setCheckable(true);
	connect(plotShowMouseCrossAct, SIGNAL(triggered(bool)), this, SLOT(plotShowMouseCross(bool)));
	
	plotAddPotentialAct = new QAction("&Add potential", this);
	plotAddPotentialAct->setStatusTip("Add a potential to the current plot");
	connect(plotAddPotentialAct, SIGNAL(triggered()), this, SLOT(plotAddPotential()));
	
	plotClearHistoryAct = new QAction("&Clear history", this);
	plotClearHistoryAct->setStatusTip("Delete all potentials and potential versions from the plot except the current version of the potential plotted first");
	connect(plotClearHistoryAct, SIGNAL(triggered()), this, SLOT(plotClearHistory()));
	
	plotShowDiagFuncsAct = new QAction("Show &diagnostic functions", this);
	plotShowDiagFuncsAct->setStatusTip("Addes to functions calculated from sums of the vibrational wave functions belonging to the observed levels first weigthed with the uncertainties and than weighted with the deviations of the levels in the fit");
	plotShowDiagFuncsAct->setCheckable(true);
	connect(plotShowDiagFuncsAct, SIGNAL(triggered(bool)), this, SLOT(plotShowDiagFuncs(bool)));
	
	plotShowHistoryAct = new QAction("Show &history", this);
	plotShowHistoryAct->setStatusTip(
		"Keep previous versions of the potential in the plot while the current potential gets updated");
	plotShowHistoryAct->setCheckable(true);
	plotShowHistoryAct->setChecked(true);
	connect(plotShowHistoryAct, SIGNAL(triggered(bool)), this, SLOT(plotShowHistory(bool)));
	
	plotShowPointsAct = new QAction("Show &points", this);
	plotShowPointsAct->setStatusTip(
						"Shows the point set of a spline potential and allows its modification");
	plotShowPointsAct->setCheckable(true);
	connect(plotShowPointsAct, SIGNAL(triggered(bool)), this, SLOT(plotShowPoints(bool)));
	
	plotPotSnapShotAct = new QAction("&Snap shot", this);
	plotPotSnapShotAct->setStatusTip("Alternative to the history function which keeps only the current version of the potential as an additional potential in the plot");
	connect(plotPotSnapShotAct, SIGNAL(triggered()), this, SLOT(plotPotSnapShot()));
	
	plotShowMarkerLabelsAct = new QAction("Sh&ow marker labels", this);
	plotShowMarkerLabelsAct->setStatusTip(
		"If this option is not checked, in a spectrum only the labels of already marked lines are shown");
	plotShowMarkerLabelsAct->setCheckable(true);
	plotShowMarkerLabelsAct->setChecked(true);
	connect(plotShowMarkerLabelsAct, SIGNAL(triggered(bool)), this, SLOT(plotShowMarkerLabels(bool)));
	
	SpectrumClearMarkedAct = new QAction("&Clear marked", this);
	SpectrumClearMarkedAct->setShortcut(tr("Ctrl+C"));
	SpectrumClearMarkedAct->setStatusTip("Remove all line selections from screen");
	connect(SpectrumClearMarkedAct, SIGNAL(triggered()), this, SLOT(SpectrumClear_Marked()));
	
	AssignvsAct = new QAction("&Assign v'", this);
	AssignvsAct->setStatusTip(
				"Uses term energies of the upper state to assign v' of the lines inside the table");
	connect(AssignvsAct, SIGNAL(triggered()), this, SLOT(Assignvs()));
	
	AssignFCAct = new QAction("A&ssign F", this);
	AssignFCAct->setStatusTip("Use term energies for both elelectronic states to assign the progressions in the table to fine structure components");
	connect(AssignFCAct, SIGNAL(triggered()), this, SLOT(AssignFC()));
	
	DeleteAct = new QAction("&Delete", this);
	DeleteAct->setShortcut(tr("Del"));
	DeleteAct->setStatusTip("Delete the contents from the selected cells of the table");
	connect(DeleteAct, SIGNAL(triggered()), this, SLOT(Delete()));
	
	FindBigDiffAct = new QAction("&Find big differences", this);
	FindBigDiffAct->setStatusTip(
			"Finds big differences between energies calculated from lines of the same progression");
	connect(FindBigDiffAct, SIGNAL(triggered()), this, SLOT(FindBigDiff()));
	
	ShowGSDeviationsAct = new QAction("Show GS Deviations", this);
	ShowGSDeviationsAct->setShortcut(tr("Shift+G"));
	ShowGSDeviationsAct->setStatusTip("Calculates and shows differences > 4 times the uncertainty of line frequency differences compared to the current term enrgy table of the lower electronic state");
	connect(ShowGSDeviationsAct, SIGNAL(triggered()), this, SLOT(ShowGSDeviations()));
	
	MarkSelectedAct = new QAction("&Mark selected", this);
	MarkSelectedAct->setShortcut(tr("Shift+M"));
	MarkSelectedAct->setStatusTip("Show and mark the selected lines from the table in the spectrum");
	connect(MarkSelectedAct, SIGNAL(triggered()), this, SLOT(MarkSelected()));
	
	TestProgressionsAct = new QAction("Test progressions", this);
	TestProgressionsAct->setStatusTip(
		"Run the assignment function 'Test Progression' on all progressions inside the current line table to update the table with the current ground state term energy set");
	connect(TestProgressionsAct, SIGNAL(triggered()), this, SLOT(TestProgressions()));
	
	RemoveDoubledAct = new QAction("&Remove doubled", this);
	RemoveDoubledAct->setStatusTip("Search table for doubled lines and delete them");
	RemoveDoubledAct->setEnabled(false);
	connect(RemoveDoubledAct, SIGNAL(triggered()), this, SLOT(RemoveDoubled()));
	
	SelectFoundAct = new QAction("Select marked", this);
	SelectFoundAct->setShortcut(tr("Shift+R"));
	SelectFoundAct->setStatusTip("Reselect lines marked in the spectra using Mark selected");
	connect(SelectFoundAct, SIGNAL(triggered()), this, SLOT(SelectFound()));
	
	FindSimilarProgressionAct = new QAction("Find similar progression...", this);
	FindSimilarProgressionAct->setStatusTip(
		"Search LineTables for progressions similar to the progression selected in the current LineTable");
	connect(FindSimilarProgressionAct, SIGNAL(triggered()), 
			this, SLOT(FindSimilarProgressions()));
	
	SetErrorAct = new QAction("Set error...", this);
	SetErrorAct->setShortcut(tr("Shift+E"));
	SetErrorAct->setStatusTip("Set the error for the selected lines...");
	connect(SetErrorAct, SIGNAL(triggered()), this, SLOT(SetError()));
	
	SetPNAct = new QAction("Set prog number", this);
	SetPNAct->setStatusTip("Set new progression numbers for all lines inside the table");
	connect(SetPNAct, SIGNAL(triggered()), this, SLOT(setPN()));
	
	setvsAct = new QAction("&Set v'...", this);
	setvsAct->setShortcut(tr("Shift+Z"));
	setvsAct->setStatusTip("Set the quantum number v' for the selected lines...");
	connect(setvsAct, SIGNAL(triggered()), this, SLOT(setvs()));
	
	setFCAct = new QAction("Set &F...", this);
	setFCAct->setStatusTip(
		"Set the fine structure quantum number F for the selected lines");
	connect(setFCAct, SIGNAL(triggered()), this, SLOT(setFC()));
	
	SetvssAscendingAct = new QAction("Set v'' &ascending", this);
	SetvssAscendingAct->setShortcut(tr("Shift+A"));
	SetvssAscendingAct->setStatusTip(
							"Set the quantum numbers v'' of the selected lines to ascending values");
	connect(SetvssAscendingAct, SIGNAL(triggered()), this, SLOT(SetvssAscending()));
	
	ShiftIsoAct = new QAction("Shift &isotope", this);
	ShiftIsoAct->setShortcut(tr("Shift+I"));
	ShiftIsoAct->setStatusTip("Shifts the isotope number of the selected lines");
	connect(ShiftIsoAct, SIGNAL(triggered()), this, SLOT(ShiftIso()));
	
	ShiftJdownAct = new QAction("Shift &J down", this);
	ShiftJdownAct->setShortcut(tr("Shift+J"));
	ShiftJdownAct->setStatusTip("Subtracts one from all J quantum numbers of the selected lines");
	connect(ShiftJdownAct, SIGNAL(triggered()), this, SLOT(ShiftJdown()));
	
	ShiftJupAct = new QAction("Shift J &up", this);
	ShiftJupAct->setShortcut(tr("Shift+U"));
	ShiftJupAct->setStatusTip("Adds one to all J quantum numbers of the selected lines");
	connect(ShiftJupAct, SIGNAL(triggered()), this, SLOT(ShiftJup()));
	
	ShiftvdownAct = new QAction("Shift v'' &down", this);
	ShiftvdownAct->setShortcut(tr("Shift+V"));
	ShiftvdownAct->setStatusTip("Subtracts one from the v'' quantum numbers of the selected lines");
	connect(ShiftvdownAct, SIGNAL(triggered()), this, SLOT(Shiftvdown()));
	
	ShiftvsdownAct = new QAction("Shi&ft v' down", this);
	ShiftvsdownAct->setShortcut(tr("Shift+S"));
	ShiftvsdownAct->setStatusTip("Subtracts one from the v' quantum numbers of the selected lines");
	connect(ShiftvsdownAct, SIGNAL(triggered()), this, SLOT(Shiftvsdown()));
	
	ShiftvsupAct = new QAction("Shif&t v' up", this);
	ShiftvsupAct->setShortcut(tr("Shift+W"));
	ShiftvsupAct->setStatusTip("Adds one to the v' quantum numbers of the selected lines");
	connect(ShiftvsupAct, SIGNAL(triggered()), this, SLOT(Shiftvsup()));
	
	ShiftvupAct = new QAction("S&hift v'' up", this);
	ShiftvupAct->setShortcut(tr("Shift+F"));
	ShiftvupAct->setStatusTip("Adds one to the v'' quantum numbers of the selected lines");
	connect(ShiftvupAct, SIGNAL(triggered()), this, SLOT(Shiftvup()));
	
	ShowUpTermAct = new QAction("Show upper term", this);
	ShowUpTermAct->setStatusTip(
							"Calculates the upper term energies and adds some columns to the table");
	connect(ShowUpTermAct, SIGNAL(triggered()), this, SLOT(ShowUpTerm()));
	
	ShowUpTermTableAct = new QAction("Show UpTermTable", this);
	ShowUpTermTableAct->setStatusTip(
						"Opens an additional table for the upper term energies of the progressions");
	connect(ShowUpTermTableAct, SIGNAL(triggered()), this, SLOT(ShowUpTermTable()));
	
	ShowCalcRelIntAct = new QAction("Calculate rel. intensities", this);
	ShowCalcRelIntAct->setStatusTip("Show an additional column with simulated relative intensities for the individual lines of a progression");
	connect(ShowCalcRelIntAct, SIGNAL(triggered()), this, SLOT(ShowCalcRelInt()));
	
	ShowWeakProgressionsAct = new QAction("Show weak progressions", this);
	ShowWeakProgressionsAct->setShortcut(tr("Shift+P"));
	ShowWeakProgressionsAct->setStatusTip("Searches for and selects weak progression in the table");
	connect(ShowWeakProgressionsAct, SIGNAL(triggered()), this, SLOT(ShowWeakProgressions()));
	
	sortbyvsAct = new QAction("Sort by &v'", this);
	sortbyvsAct->setStatusTip("Sorts the table by the quantum number v'");
	connect(sortbyvsAct, SIGNAL(triggered()), this, SLOT(sortbyvs()));
	
	SortProgAct = new QAction("Sort by &progressions", this);
	SortProgAct->setStatusTip("Sorts the table by progressions");
	connect(SortProgAct, SIGNAL(triggered()), this, SLOT(SortProg()));
	
	sortUpTermIvJAct = new QAction("Sort &UpTermIvJ", this);
	sortUpTermIvJAct->setStatusTip(
		"Sorts the UpTerm table by the isotope and the upper quantum numbers v' and J' and the upper term energy"); 
	connect(sortUpTermIvJAct, SIGNAL(triggered()), this, SLOT(sortUpTermIvJ()));
	
	SortIJvPAct = new QAction("Sort I&JvProg", this);
	SortIJvPAct->setStatusTip("Sorts the table first by the isotopomer, than by J' and v' and finally by the progression and v'' and J''");
	connect(SortIJvPAct, SIGNAL(triggered()), this, SLOT(SortIJvP()));
	
	SortIvPJAct = new QAction("Sort &IvJProg", this);
	SortIvPJAct->setStatusTip("Sorts the table first by the isotopomer, than by v' and J' and finally by the progression and v'' and J''");
	connect(SortIvPJAct, SIGNAL(triggered()), this, SLOT(SortIvPJ()));
	
	SortFPIntAct = new QAction("Sort &FPInt", this);
	SortFPIntAct->setStatusTip("Sorts the table first by the spectrum, than by the progression and finally by the signal-to-noise ratio");
	connect(SortFPIntAct, SIGNAL(triggered()), this, SLOT(SortFPInt()));
	
	sortIvJfFAct = new QAction("&Sort IvJf", this);
	sortIvJfFAct->setStatusTip("Sorts the data set first by the isotopologue, than by v'', J'', v', J' and finally by the level energy / line frequency");
	connect(sortIvJfFAct, SIGNAL(triggered()), this, SLOT(sortIvJfF()));
	
	sortvsIvJAct = new QAction("Sort &vsIvJ", this);
	sortvsIvJAct->setStatusTip("Sorts the data set firs by v', than by the isotopologue, v'', J' and finally by J''");
	connect(sortvsIvJAct, SIGNAL(triggered()), this, SLOT(sortvsIvJ()));
	
	SortSpectrumAct = new QAction("Sort Sp&ectrum", this);
	SortSpectrumAct->setStatusTip("Sorts the lines in the table by the file name of the source spectrum");
	connect(SortSpectrumAct, SIGNAL(triggered()), this, SLOT(SortSpectrum()));
	
	sortTabByDeviationAct = new QAction("Sort by &deviation", this);
	sortTabByDeviationAct->setStatusTip("Sorts the fit dataset by the obs-calc values");
	connect(sortTabByDeviationAct, SIGNAL(triggered()), this, SLOT(sortTabByDeviation()));
	
	sortTabByDevRatioAct = new QAction("Sort by Dev&Ratio", this);
	sortTabByDevRatioAct->setStatusTip(
		"Sorts the fit dataset by the ratios of the obs-calc values and the uncertainties");
	connect(sortTabByDevRatioAct, SIGNAL(triggered()), this, SLOT(sortTabByDevRatio()));

    sortTabByElStateAct = new QAction("Sort by electronic s&tate", this);
    sortTabByElStateAct->setStatusTip("Sorts a fit dataset containing data for different electronic states by the electronic state");
    connect(sortTabByElStateAct, SIGNAL(triggered()), this, SLOT(sortTabByElState()));
	
	SortfRemDoubledAct = new QAction("Sort f.rem. &doubled", this);
	SortfRemDoubledAct->setStatusTip("Sorts the lines in the table as it is temporaly done for removing doubled lines from the table");
	connect(SortfRemDoubledAct, SIGNAL(triggered()), this, SLOT(SortfRemDoubled()));

    sortTabByElStateAndProgNrAct = new QAction("Sort by state and prog num", this);
    sortTabByElStateAndProgNrAct->setStatusTip("Sorts a fit dataset containing data for different electronic states by electronic state and progression number.");
    connect(sortTabByElStateAndProgNrAct, SIGNAL(triggered()), this, SLOT(sortTabByElStateAndProgNr()));

    SortByProgNrAct = new QAction("Sort by progression num", this);
    SortByProgNrAct->setStatusTip("Sort the lines in the line table by the progression number.");
    connect(SortByProgNrAct, SIGNAL(triggered()), this, SLOT(SortByProgNr()));
	
	SplitLineTableAct = new QAction("Split table", this);
	SplitLineTableAct->setStatusTip("Shift all lines which are assigned to v'=-1 to a new line table created for an unknown upper state of the same molecule\n");
	connect(SplitLineTableAct, SIGNAL(triggered()), this, SLOT(SplitLineTable()));
	
	TakeOnChangesAct = new QAction("Take on changes", this);
	TakeOnChangesAct->setShortcut(tr("Shift+O"));
	TakeOnChangesAct->setStatusTip(
					"Takes on to the table changes in the assignment after calling 'mark selected'");
	connect(TakeOnChangesAct, SIGNAL(triggered()), this, SLOT(TakeOnChanges()));
	
	showTermViewAct = new QAction("Show term &view", this);
	showTermViewAct->setShortcut(tr("Shift+T"));
	showTermViewAct->setStatusTip("Shows a table which gives an overview of the term energies, doublet separations or vibrational differences");
	connect(showTermViewAct, SIGNAL(triggered()), this, SLOT(showTermView()));
	
	showTermPlotAct = new QAction("Show term plot", this);
	showTermPlotAct->setStatusTip("Show a plot of the term energies of the actual electronic state");
	connect(showTermPlotAct, SIGNAL(triggered()), this, SLOT(showTermPlot()));
	
	showTexTableAct = new QAction("Show &tex table...", this);
	showTexTableAct->setShortcut(tr("F8"));
	showTexTableAct->setStatusTip(
		"Shows a text field containing latex code for a table of the contents of the current window");
	connect(showTexTableAct, SIGNAL(triggered()), this, SLOT(showTexTable()));
	
	showFCFTableAct = new QAction("Show &FCF view", this);
	showFCFTableAct->setStatusTip(
		"Shows a table with Franck-Condon factors and transitions frequencies for a selectable electronic transition and J");
	connect(showFCFTableAct, SIGNAL(triggered()), this, SLOT(showFCFTable()));
	
	showFCFJDependencyAct = new QAction("Show FCF &JDep...", this);
	showFCFJDependencyAct->setStatusTip("Show a J dependent plot of Franck-Condon factors");
	connect(showFCFJDependencyAct, SIGNAL(triggered()), this, SLOT(showFCFJDependency()));
	
	showWaveFuncPlotAct = new QAction("Show &wave functions", this);
	showWaveFuncPlotAct->setStatusTip(
					"Show a plot of a selectable wave function together with the according potential");
	connect(showWaveFuncPlotAct, SIGNAL(triggered()), this, SLOT(showWaveFuncPlot()));
	
	showSpectSimulationAct = new QAction("Show spectrum simulation", this);
	showSpectSimulationAct->setStatusTip("Open a window for a simulation of an absorption, light induced fluorescence or filtered absorption spectrum");
	connect(showSpectSimulationAct, SIGNAL(triggered()), this, SLOT(showSpectSimulation()));
	
	showResidualPlotAct = new QAction("Show &residuals", this);
	showResidualPlotAct->setStatusTip("Show a plot of the energy differences of a selectable band for two different data sets as a function of J");
	connect(showResidualPlotAct, SIGNAL(triggered()), this, SLOT(showResidualPlot()));
	
	showDataPlotAct = new QAction("Show data f&ield", this);
	showDataPlotAct->setStatusTip(
	"Open a window for a plot of the data field included in a term energy or given by a line table");
	connect(showDataPlotAct, SIGNAL(triggered()), this, SLOT(showDataPlot()));
	
	SpectrumAutoSLPAct = new QAction("&Automatically analyse spectra", this);
	SpectrumAutoSLPAct->setShortcut(tr("F7"));
	SpectrumAutoSLPAct->setStatusTip(
				"Automatically assign all lines which can be assigned using the current data sets");
	connect(SpectrumAutoSLPAct, SIGNAL(triggered()), this, SLOT(SpectrumAutoSLP()));
	
	SpectrumChangeSettingsAct = new QAction("C&hange settings", this);
	SpectrumChangeSettingsAct->setShortcut(tr("F2"));
	SpectrumChangeSettingsAct->setStatusTip(
				"Changes settings as the minimum height of a peak and the search tolerancy for lines");
	connect(SpectrumChangeSettingsAct, SIGNAL(triggered()), this, SLOT(SpectrumChange_Settings()));
	
	SpectrumDisplayMarkedAct = new QAction("&Display marked", this);
	SpectrumDisplayMarkedAct->setStatusTip("Show the assigned lines");
	connect(SpectrumDisplayMarkedAct, SIGNAL(triggered()), this, SLOT(SpectrumDisplay_Marked()));
	
	SpectrumFindPeaksAct = new QAction("&Find peaks", this);
	SpectrumFindPeaksAct->setShortcut(tr("Ctrl+F"));
	SpectrumFindPeaksAct->setStatusTip("Find and mark the peaks inside the spectrum");
	connect(SpectrumFindPeaksAct, SIGNAL(triggered()), this, SLOT(SpectrumFindPeaks()));
	
	SpectrumFind_SatellitesAct = new QAction("F&ind satellites", this);
	SpectrumFind_SatellitesAct->setShortcut(tr("Ctrl+I"));
	SpectrumFind_SatellitesAct->setStatusTip(
									"Find satellites fitting to the currently selected progression");
	connect(SpectrumFind_SatellitesAct, SIGNAL(triggered()), this, SLOT(SpectrumFind_Satellites()));
	
	SpectrumFindEmissionLinesAct = new QAction("Find e&mmission lines...", this);
	SpectrumFindEmissionLinesAct->setStatusTip(
		"Find emission or absorption bands fitting to calculated term energies of ground and excited state");
	connect(SpectrumFindEmissionLinesAct, SIGNAL(triggered()), this, SLOT(SpectrumFindEmissionLines()));
	
	SpectrumFindProgressionsAct = new QAction("Find p&rogressions", this);
	connect(SpectrumFindProgressionsAct, SIGNAL(triggered()), this, SLOT(SpectrumFindProgressions()));
	
	SpectrumMultiSLPAct = new QAction("A&utomatically assign the lines", this);
	SpectrumMultiSLPAct->setShortcut(tr("F6"));
	SpectrumMultiSLPAct->setStatusTip(
		"Automatically assign all lines which can be assigned using the currently available data set");
	connect(SpectrumMultiSLPAct, SIGNAL(triggered()), this, SLOT(SpectrumMultiSLP()));
	
	SpectrumSingleSLPAct = new QAction("Show &longest progression", this);
	SpectrumSingleSLPAct->setShortcut(tr("F5"));
	SpectrumSingleSLPAct->setStatusTip(
			"Search for the longest unassigned progression in spectrum or to the selected line(s)");
	connect(SpectrumSingleSLPAct, SIGNAL(triggered()), this, SLOT(SpectrumSingleSLP()));
	
	SpectrumPrev_ProgressionAct = new QAction("Show &previous progression", this);
	SpectrumPrev_ProgressionAct->setShortcut(tr("P"));
	SpectrumPrev_ProgressionAct->setStatusTip("Show the previous assignment suggestion");
	connect(SpectrumPrev_ProgressionAct, SIGNAL(triggered()), this, SLOT(SpectrumPrev_Progression()));
	
	SpectrumNext_ProgressionAct = new QAction("Show &next progression", this);
	SpectrumNext_ProgressionAct->setShortcut(tr("N"));
	SpectrumNext_ProgressionAct->setStatusTip("Show the next assignment suggestion");
	connect(SpectrumNext_ProgressionAct, SIGNAL(triggered()), this, SLOT(SpectrumNext_Progression()));
	
	SpectrumAcceptAssignmentsAct = new QAction("Accept assi&gnments", this);
	SpectrumAcceptAssignmentsAct->setShortcut(tr("Ctrl+A"));
	SpectrumAcceptAssignmentsAct->setStatusTip("Accept the current suggestion for assignment");
	connect(SpectrumAcceptAssignmentsAct, SIGNAL(triggered()), 
			this, SLOT(SpectrumAcceptAssignments()));
	
	SpectrumSetLaserFrequencyAct = new QAction("&Set laser frequency", this);
	SpectrumSetLaserFrequencyAct->setShortcut(tr("Ctrl+L"));
	SpectrumSetLaserFrequencyAct->setStatusTip("Change the laser frequency set");
	connect(SpectrumSetLaserFrequencyAct, SIGNAL(triggered()), 
			this, SLOT(SpectrumSetLaserFrequency()));
	
	SpectrumShowAssignmentsOnTopAct = new QAction("Show assignments on &top", this);
	SpectrumShowAssignmentsOnTopAct->setStatusTip("Show all assignments on top of the spectrum");
	connect(SpectrumShowAssignmentsOnTopAct, SIGNAL(triggered()), 
			this, SLOT(SpectrumShowAssignmentsOnTop()));
	
	SpectrumTestProgressionAct = new QAction("T&est progression", this);
	SpectrumTestProgressionAct->setShortcut(tr("Ctrl+R"));
	SpectrumTestProgressionAct->setStatusTip(
			"Test if the currently selected lines form a progression and find more lines");
	connect(SpectrumTestProgressionAct, SIGNAL(triggered()), this, SLOT(SpectrumTest_Transition()));
	
	SpectrumAssignBandsAct = new QAction("Assign bands", this);
	SpectrumAssignBandsAct->setStatusTip(
						"Switch the mode for manual assignments to the assignmet of band structures");
	connect(SpectrumAssignBandsAct, SIGNAL(triggered()), this, SLOT(SpectrumAssignBands()));
	
	SpectrumAssignBandsByDPAct = new QAction("Assign bands by DP", this);
	SpectrumAssignBandsByDPAct->setStatusTip(
		"Assign gapless series of lines of a P or R band structure by the number of findable doublet partners using groundstate term energies.");
	connect(SpectrumAssignBandsByDPAct, SIGNAL(triggered()), this, SLOT(SpectrumAssignBandsByDP()));
	
	SpectrumFindLinesFromTableAct = new QAction("Find observed lines...", this);
	SpectrumFindLinesFromTableAct->setStatusTip(
		"Search for lines already observed (in other spectra) taken from a selectable linetable");
	connect(SpectrumFindLinesFromTableAct, SIGNAL(triggered()), this, SLOT(SpectrumFindLinesFromTable()));
	
	SpectrumContinueProgressionsAct = new QAction("Continue progressions", this);
	SpectrumContinueProgressionsAct->setStatusTip(
			"Searches for additional fitting lines to all partially selected progressions");
	connect(SpectrumContinueProgressionsAct, SIGNAL(triggered()), 
			this, SLOT(SpectrumContinueProgressions()));

    SpectrumFitGaussianLineProfileAct = new QAction("Fit Gaussian line profile", this);
    SpectrumFitGaussianLineProfileAct->setStatusTip("Fit a Gaussian line profile to the marked, red colored region of the spectrum.");
    connect(SpectrumFitGaussianLineProfileAct, SIGNAL(triggered()), this, SLOT(SpectrumFitGaussianLineProfile()));
	
	SpectrumCutAct = new QAction("Cu&t...", this);
	SpectrumCutAct->setStatusTip("Cut the spectrum");
	connect(SpectrumCutAct, SIGNAL(triggered()), this, SLOT(SpectrumCut()));
	
	SpectrumCutAssignedAct = new QAction("Cut assigned", this);
	SpectrumCutAssignedAct->setStatusTip(
									"Remove all assigned lines from the spectrum by cutting them out");
	connect(SpectrumCutAssignedAct, SIGNAL(triggered()), this, SLOT(SpectrumCutAssigned()));
	
	
	SpectrumCutStrongAct = new QAction("Cut strong lines...", this);
	SpectrumCutStrongAct->setStatusTip("Cut out all lines stronger than a selectable intensity");
	connect(SpectrumCutStrongAct, SIGNAL(triggered()), this, SLOT(SpectrumCutStrong()));
	
	SpectrumAddAct = new QAction("Add...", this);
	SpectrumAddAct->setStatusTip("Add shifted versions of the active spectrum to overlap the lines of a specific band or progression");
	connect(SpectrumAddAct, SIGNAL(triggered()), this, SLOT(SpectrumAdd()));
	
	SpectrumSetTypeAct = new QAction("Set t&ype...", this);
	SpectrumSetTypeAct->setStatusTip("The type of a spectrum can either be a laser induced fluorescence spectrum, an absorption spectrum or a normalized absorption spectrum");
	connect(SpectrumSetTypeAct, SIGNAL(triggered()), this, SLOT(SpectrumSetType()));
	
	SpectrumNormalizeAct = new QAction("N&ormalize...", this);
	SpectrumNormalizeAct->setStatusTip("Takes the difference with a selected reference spectrum and than divides it by the reference spectrum");
	connect(SpectrumNormalizeAct, SIGNAL(triggered()), this, SLOT(SpectrumNormalize()));
	
	SpectrumShowIntensityHistogramAct = new QAction("Show intensity histogram", this);
	SpectrumShowIntensityHistogramAct->setStatusTip(
	   "Show a histogram of the peak intensities calculated relative to the deepest of the neighboring minima");
	connect(SpectrumShowIntensityHistogramAct, SIGNAL(triggered()), 
			this, SLOT(SpectrumShowIntensityHistogram()));

    aboutAct = new QAction("&About", this);
    aboutAct->setStatusTip("About this program");
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(ShowAboutWindow()));
	
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileNewMenu = fileMenu->addMenu(tr("&New"));
	fileNewMenu->addAction(newAtomAct);
	fileNewMenu->addAction(newMoleculeAct);
	fileNewMenu->addSeparator();
	fileNewMenu->addAction(newTermTableAct);
	fileNewMenu->addAction(newMeasuredTermTableAct);
	fileNewMenu->addAction(newDunTableAct);
	fileNewMenu->addAction(newPotentialAct);
	fileNewMenu->addAction(newFitDataAct);
	fileNewMenu->addSeparator();
	fileNewMenu->addAction(newLineTableAct);
	fileNewMenu->addAction(newFCFTableAct);
	fileNewMenu->addSeparator();
	fileNewMenu->addAction(newSpectrumAct);
	fileNewMenu->addAction(newSpectListAct);
	
	fileMenu->addAction(openAct);
	fileMenu->addAction(reloadAct);
	
	fileMenu->addSeparator();
	
	fileMenu->addAction(saveAct);
	fileMenu->addAction(saveAsAct);
	
	fileMenu->addSeparator();
	fileMenu->addAction(printAct);
	fileMenu->addSeparator();
	
	fileMenu->addSeparator();
	importMenu = fileMenu->addMenu(tr("&Import"));
	importMenu->addAction(importTermTableAct);
	importMenu->addAction(importDunTableAct);
	importMenu->addAction(importPotentialAct);
	importMenu->addAction(importFitDataAct);
	importMenu->addAction(importLineTableAct);
	importMenu->addAction(importSpectrumAct);
	importMenu->addAction(importSpectListAct);
	importMenu->addSeparator();
	importMenu->addAction(importAbInitioPotSetAct);
	importMenu->addAction(importMTTPotentialsAct);
	importMenu->addAction(read2AtOutputAct);
	importMenu->addAction(importCoupledPotfitOutputAct);
	importMenu->addAction(importCoupledTermTableAct);
	importMenu->addAction(importCoupledWaveFunctionsAct);
	importMenu->addAction(importDunhamAsenAct);
	
	exportMenu = fileMenu->addMenu(tr("&Export"));
	exportMenu->addAction(exportPictureAct);
	exportMenu->addAction(exportTableDataAct);
	exportMenu->addSeparator();
	exportMenu->addAction(writePotDataAct);
	exportMenu->addAction(exportPotPointsAct);
	exportMenu->addAction(exportAsymptoticLevelsAct);
	exportMenu->addAction(exportWaveFunctionAct);
    exportMenu->addSeparator();
    exportMenu->addAction(exportLineProfileAct);
	exportMenu->addSeparator();
	exportMenu->addAction(writeTFGSAct);
	exportMenu->addAction(writeExpPotFitInputAct);
	exportMenu->addAction(write2AtInputAct);
	exportMenu->addAction(exportTFDunAct);
	exportMenu->addSeparator();
	exportMenu->addAction(exportObservedLevelListsAct);
	exportMenu->setEnabled(false);
	
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);
	
	tableMenu = menuBar()->addMenu("&Table");
	tableMenu->addAction(addRowAct);
	tableMenu->addAction(insertRowsAct);
	tableMenu->addAction(copyRowsAct);
	tableMenu->addAction(cutRowsAct);
	tableMenu->addAction(deleteRowsAct);
	tableMenu->addSeparator();
	tableMenu->addAction(searchAct);
	tableMenu->addSeparator();
	setMenu = tableMenu->addMenu("Set colu&mn value");
    setMenu->addAction(setTextInCellsAct);
    setMenu->addSeparator();
	setMenu->addAction(SetErrorAct);
	setMenu->addAction(setvsAct);
	setMenu->addAction(SetvssAscendingAct);
	setMenu->addAction(setFCAct);
	setMenu->addAction(SetPNAct);
	tableMenu->addSeparator();
	tableMenu->addAction(AssignvsAct);
	tableMenu->addAction(RemoveDoubledAct);
	tableMenu->addAction(findWrongDataAct);
	tableMenu->addSeparator();
    tableMenu->addAction(extractDataWithComponentAct);
    tableMenu->addAction(extractNewDataAct);
    tableMenu->addAction(extractChangedDataAct);
    tableMenu->addSeparator();
	tableMenu->addAction(selectDataFSourceAct);
	tableMenu->addAction(removeDataFSourceAct);
	tableMenu->addAction(removeSingleLinesAct);
	tableMenu->addSeparator();
	sortMenu = tableMenu->addMenu("Sort ta&ble");
	tableMenu->addAction(shiftValuesAct);
	tableMenu->addSeparator();
	tableMenu->addAction(updateAct);
	tableMenu->addAction(addCalculatedLevelsAct);
	tableMenu->addAction(fixJOffsetsAct);
	tableMenu->addSeparator();
	tableMenu->addAction(showNumLevelsAct);
	tableMenu->addAction(showUncertaintyStatsAct);
	tableMenu->addAction(compareLevelEnergiesAct);
	tableMenu->addSeparator();
	tableMenu->addAction(changeSourceOffsetAct);
	tableMenu->setEnabled(false);
	
	SpectrumMenu = menuBar()->addMenu("&Spectrum");
	SpectrumMenu->addAction(SpectrumSetTypeAct);
	SpectrumMenu->addAction(SpectrumChangeSettingsAct);
	SpectrumMenu->addAction(SpectrumSetLaserFrequencyAct);
	SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumFindPeaksAct);
	SpectrumMenu->addAction(SpectrumShowIntensityHistogramAct);
	SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumDisplayMarkedAct);
	SpectrumMenu->addAction(SpectrumShowAssignmentsOnTopAct);
	SpectrumMenu->addAction(SpectrumAcceptAssignmentsAct);
	SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumFindLinesFromTableAct);
	SpectrumMenu->addAction(SpectrumFindEmissionLinesAct);
	SpectrumMenu->addAction(SpectrumAssignBandsByDPAct);
	SpectrumMenu->addAction(SpectrumAssignBandsAct);
	SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumPrev_ProgressionAct);
	SpectrumMenu->addAction(SpectrumNext_ProgressionAct);
	SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumTestProgressionAct);
	SpectrumMenu->addAction(SpectrumFind_SatellitesAct);
	SpectrumMenu->addAction(SpectrumFindProgressionsAct);
	SpectrumMenu->addAction(SpectrumContinueProgressionsAct);
	SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumSingleSLPAct);
	SpectrumMenu->addAction(SpectrumMultiSLPAct);
	SpectrumMenu->addAction(SpectrumAutoSLPAct);
	SpectrumMenu->addSeparator();
    SpectrumMenu->addAction(SpectrumFitGaussianLineProfileAct);
    SpectrumMenu->addSeparator();
	SpectrumMenu->addAction(SpectrumNormalizeAct);
	SpectrumMenu->addAction(SpectrumCutAct);
	SpectrumMenu->addAction(SpectrumCutStrongAct);
	SpectrumMenu->addAction(SpectrumCutAssignedAct);
	SpectrumMenu->addAction(SpectrumAddAct);
	SpectrumMenu->setEnabled(false);
	
	MoleculeMenu = menuBar()->addMenu("&Molecule");
	MoleculeMenu->addAction(showIsotopologuesAct);
	MoleculeMenu->addAction(showLevelNumbersAct);
	MoleculeMenu->addSeparator();
	MoleculeMenu->addAction(calcFCFTableAct);
    MoleculeMenu->addSeparator();
    MoleculeMenu->addAction(checkAllConnectionsAct);
    MoleculeMenu->addAction(shrinkAllSpektRefsAct);
	MoleculeMenu->setEnabled(false);
	
	DunhamMenu = menuBar()->addMenu("&Dunham");
	DunhamMenu->addAction(calcTermEnergiesAct);
	DunhamMenu->addSeparator();
	DunhamMenu->addAction(updateDunhamAct);
	DunhamMenu->addAction(improveDunhamAct);
	DunhamMenu->addAction(fitBorderLineAct);
	DunhamMenu->addSeparator();
	DunhamMenu->addAction(removeLambdaDoublingAct);
	DunhamMenu->addAction(removeFineStructureAct);
	DunhamMenu->addAction(applyAdCorrAct);
	DunhamMenu->addSeparator();
	DunhamMenu->addAction(testKratzerAct);
	DunhamMenu->addAction(calcReAct);
	DunhamMenu->addAction(calcY00TeAct);
	DunhamMenu->addAction(showBandConstantsAct);
	DunhamMenu->setEnabled(false);
	
	PotentialMenu = menuBar()->addMenu("&Potential");
	PotentialMenu->addAction(calcTermEnergiesPotAct);
	PotentialMenu->addAction(showTeReAct);
    PotentialMenu->addAction(calcRoughBeAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(calcScatWaveFuncAct);
	PotentialMenu->addAction(showScatNodePosAct);
	PotentialMenu->addAction(fitScatNodePosAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(summarizePotInfoAct);
	PotentialMenu->addAction(calcSFQSMCSSeriesAct);
	PotentialMenu->addAction(improvePotSeriesAct);
	PotentialMenu->addAction(createAnaPotsFromMCSSplineSeriesAct);
	PotentialMenu->addAction(autoCalcScatLengthsPotentialSetAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(monteCarloSimAct);
	PotentialMenu->addAction(fitControlAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(plotPotentialAct);
	PotentialMenu->addAction(setEnergyOffsetAct);
	PotentialMenu->addAction(scalePotentialAct);
	PotentialMenu->addAction(selPotCoeffAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(cdConnectSRAct);
	PotentialMenu->addAction(cdConnectLRAct);
	PotentialMenu->addAction(calcyssAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(fixCoefficientsAct);
	PotentialMenu->addAction(varyCoefficientsAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(showFitDataAct);
	PotentialMenu->addAction(updateFitDataAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(showFQSAct);
	PotentialMenu->addAction(showBadListAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(fitSplinePotAct);
	PotentialMenu->addAction(fitAnaPotAct);
	PotentialMenu->addAction(fitTangToenniesPotAct);
	PotentialMenu->addAction(createMLRPotAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(fitIsoMassAct);
	PotentialMenu->addAction(fitWRobustWeightingAct);
	PotentialMenu->addSeparator();
	PotentialMenu->addAction(showSFuncsAct);
	enablePotMenu(false, NoPotential);
	improvePotSeriesAct->setEnabled(false);
	
	PlotMenu = menuBar()->addMenu("Pl&ot");
	PlotMenu->addAction(plotShowMouseCrossAct);
	PlotMenu->addSeparator();
	PlotMenu->addAction(SpectrumClearMarkedAct);
	PlotMenu->addAction(plotShowMarkerLabelsAct);
	PlotMenu->addSeparator();
	PlotMenu->addAction(plotShowPointsAct);
	PlotMenu->addAction(plotShowDiagFuncsAct);
	PlotMenu->addAction(plotShowHistoryAct);
	PlotMenu->addSeparator();
	PlotMenu->addAction(plotClearHistoryAct);
	PlotMenu->addAction(plotPotSnapShotAct);
	PlotMenu->addAction(plotAddPotentialAct);
	PlotMenu->setEnabled(false);
	
	lineTableMenu = menuBar()->addMenu("&Line table");
	shiftMenu = lineTableMenu->addMenu("Shift column values");
	shiftMenu->addAction(ShiftvupAct);
	shiftMenu->addAction(ShiftvdownAct);
	shiftMenu->addSeparator();
	shiftMenu->addAction(ShiftvsupAct);
	shiftMenu->addAction(ShiftvsdownAct);
	shiftMenu->addSeparator();
	shiftMenu->addAction(ShiftJupAct);
	shiftMenu->addAction(ShiftJdownAct);
	shiftMenu->addSeparator();
	shiftMenu->addAction(ShiftIsoAct);
	lineTableMenu->addSeparator();
	lineTableMenu->addAction(AssignFCAct);
	lineTableMenu->addAction(SplitLineTableAct);
	lineTableMenu->addSeparator();
	lineTableMenu->addAction(MarkSelectedAct);
	lineTableMenu->addAction(TakeOnChangesAct);
	lineTableMenu->addAction(SelectFoundAct);
	lineTableMenu->addAction(TestProgressionsAct);
	lineTableMenu->addSeparator();
	lineTableMenu->addAction(DeleteAct);
	lineTableMenu->addSeparator();
	lineTableMenu->addAction(ShowCalcRelIntAct);
	lineTableMenu->addAction(ShowUpTermAct);
	lineTableMenu->addAction(ShowUpTermTableAct);
	lineTableMenu->addSeparator();
	lineTableMenu->addAction(FindBigDiffAct);
	lineTableMenu->addAction(ShowGSDeviationsAct);
	lineTableMenu->addAction(ShowWeakProgressionsAct);
	lineTableMenu->addAction(FindSimilarProgressionAct);
	lineTableMenu->setEnabled(false);
	
	ShowMenu = menuBar()->addMenu("&Show");
	ShowMenu->setEnabled(false);
	ShowAtoms = ShowMenu->addMenu("&Atom");
	ShowAtoms->setEnabled(false);
	ShowMolecules = ShowMenu->addMenu("&Molecule");
	ShowMolecules->setEnabled(false);
	ShowMenu->addSeparator();
	ShowTermTables = ShowMenu->addMenu("&TermTable");
	ShowTermTables->setEnabled(false);
	ShowMeasuredTermTables = ShowMenu->addMenu("M&esuredTermTable");
	ShowMeasuredTermTables->setEnabled(false);
	ShowDunTables = ShowMenu->addMenu("&Dunham coefficients");
	ShowDunTables->setEnabled(false);
	ShowPotentials = ShowMenu->addMenu("&Potential");
	ShowPotentials->setEnabled(false);
	ShowFitDataSets = ShowMenu->addMenu("F&it data");
	ShowFitDataSets->setEnabled(false);
	ShowMenu->addSeparator();
	ShowLineTables = ShowMenu->addMenu("&Line table");
	ShowLineTables->setEnabled(false);
	ShowFCFTables = ShowMenu->addMenu("&FCF table");
	ShowFCFTables->setEnabled(false);
	ShowMenu->addSeparator();
	ShowSpectra = ShowMenu->addMenu("&Spectrum");
	ShowSpectLists = ShowMenu->addMenu("Spe&ctList");
	ShowMenu->addSeparator();
	ShowMenu->addAction(showFCFTableAct);
	ShowMenu->addAction(showFCFJDependencyAct);
	ShowMenu->addSeparator();
	ShowMenu->addAction(showWaveFuncPlotAct);
	ShowMenu->addAction(showSpectSimulationAct);
	ShowMenu->addSeparator();
	ShowMenu->addAction(showDataPlotAct);
	ShowMenu->addAction(showResidualPlotAct);
	ShowMenu->addSeparator();
	ShowMenu->addAction(showTermViewAct);
	showTermViewAct->setEnabled(false);
	ShowMenu->addAction(showTermPlotAct);
	showTermPlotAct->setEnabled(false);
	ShowMenu->addSeparator();
	ShowMenu->addAction(showTexTableAct);
	showTexTableAct->setEnabled(false);

    menuBar()->addAction(aboutAct);
	//printf("End MainWindow::MainWindow\n");
}

MainWindow::~MainWindow()
{
	delete[] DirRep;
	delete[] DirRepNC;
	if (cuttedRows != 0) Destroy(cuttedRows, nCR);
	if (windowsToAskForQuit != 0) delete[] windowsToAskForQuit;
}

void MainWindow::enablePotMenu(bool enable, PotentialType PotType)
{
	calcTermEnergiesPotAct->setEnabled(enable);
	showTeReAct->setEnabled(enable);
    calcRoughBeAct->setEnabled(enable);
	calcScatWaveFuncAct->setEnabled(enable);
	showScatNodePosAct->setEnabled(enable);
	fitScatNodePosAct->setEnabled(enable);
	autoCalcScatLengthsPotentialSetAct->setEnabled(enable);
	monteCarloSimAct->setEnabled(enable);
	plotPotentialAct->setEnabled(enable);
	setEnergyOffsetAct->setEnabled(enable);
	scalePotentialAct->setEnabled(enable);
	selPotCoeffAct->setEnabled(enable);
	fixCoefficientsAct->setEnabled(enable);
	varyCoefficientsAct->setEnabled(enable);
	showFitDataAct->setEnabled(enable);
	updateFitDataAct->setEnabled(enable);
	showFQSAct->setEnabled(enable);
	fitSplinePotAct->setEnabled(enable);
	fitAnaPotAct->setEnabled(enable);
	fitTangToenniesPotAct->setEnabled(enable);
	fitWRobustWeightingAct->setEnabled(enable);
	fitIsoMassAct->setEnabled(enable);
	createMLRPotAct->setEnabled(enable);
	showBadListAct->setEnabled(enable);
	showSFuncsAct->setEnabled(enable);
	cdConnectLRAct->setEnabled(enable);
	cdConnectSRAct->setEnabled(enable);
	calcyssAct->setEnabled(enable && PotType == SplinePotential);
}

int MainWindow::getNumAtoms()
{
	return numAtoms;
}

Atom *MainWindow::getAtom(int Index)
{
	if (Index >= 0 && Index < numAtoms) return atoms[Index];
	printf("MainWindow::getAtom: the Atom with Index = %d does not exist!", Index);
	return 0;
}

Atom *MainWindow::getAtom(QString FileName)
{
	int i;
	//printf("MainWindow::getAtom: Filename=%s\n", FileName.ascii());
	for (i=0; i<numAtoms; i++) if (atoms[i]->getFileName() == FileName) return atoms[i];
	QString FN = FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numAtoms; i++) if (FN == atoms[i]->getFName()) return atoms[i];
	for (i=0; i < numNAAtoms; i++) if (FileName == NAAtoms[i]) return 0;
	Atom *nAtom = CreateAtom();
	if (nAtom != 0) if (!nAtom->readData(FileName))
	{
		nAtom->close();
		numAtoms--;
		if (numNAAtoms < MaxAtoms) NAAtoms[numNAAtoms++] = FileName;
		return 0;
	}
	return nAtom;
}	

Atom *MainWindow::getAtomN(QString Name)
{
	//printf("MainWindow::getAtomN(%s)\n", Name.ascii());
	int i;
	for (i=0; i<numAtoms; i++) if (atoms[i]->getName() == Name) 
	{
		printf("return atoms[%d]\n", i);
		return atoms[i];
	}
	printf("MainWindow::getAtomN: the atom with Name = %s can't be found!\n", Name.toLatin1().data());
	return 0;
}

int MainWindow::getNumMolecules()
{
	return numMolecules;
}

Molecule *MainWindow::getMolecule(int Index)
{
	if (Index >= 0 && Index < numMolecules) return molecules[Index];
	printf("MainWindow::getMolecule: the molecule with Index = %d does not exist!\n", Index);
	return 0;
}

Molecule *MainWindow::getMolecule(QString Name)
{
	int i;
	for (i=0; i < numMolecules; i++) if (molecules[i]->getName() == Name) return molecules[i];
	printf("MainWindow::getMolecule: a molecule with the name %s does not exist!\n", Name.toLatin1().data());
	return 0;
}

Molecule* MainWindow::getMoleculewFM(QString FileName)
{
	int n;
	for (n=0; n < numMolecules; n++) if (molecules[n]->getFileName() == FileName) return molecules[n];
	Molecule *nMol = CreateMolecule();
	if (nMol == 0) return nMol;
	if (nMol->readData(FileName)) return nMol;
	return 0;
}

int MainWindow::getNumTermTables()
{
	return numTermTables;
}

TermTable *MainWindow::getTermTable(int Index)
{
	if (Index >= 0 && Index < numTermTables) return termTables[Index];
	printf("MainWindow::getTermTable: the table with Index = %d does not exist!", Index);
	return 0;
}

TermTable *MainWindow::getTermTable(QString FileName, Molecule *molecule, ElState *State)
{
	int i;
	QString S1;
	for (i=0; i<numTermTables; i++) /*{*/if ((S1 = termTables[i]->getFileName()) == FileName) 
		return termTables[i]; 
		//printf("TermTables[%d] = %s, Filename = %s!\n", i, S1.ascii(), FileName.ascii());}
	QString FN = FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numTermTables; i++) if (FN == termTables[i]->getFName()) return termTables[i];
	for (i=0; i < numNATermTables; i++) if (NATermTables[i] == FileName) return 0;
	TermTable *nTerm = CreateTermTable();
	if (nTerm == 0) return 0;
	if (State != 0 ? (State == molecule->getStateP(0) ? State->getTermTableFileName() != FileName : true) : false) 
		nTerm->setMolecule(molecule, State);
	if (!nTerm->readData(FileName))
	{
		nTerm->close();
		numTermTables--;
		if (numNATermTables < MaxTermTables) NATermTables[numNATermTables++] = FileName;
		return 0;
	}
	return nTerm;
}

int MainWindow::getNumMeasuredTermTables()
{
	return numMeasuredTermTables;
}

MeasuredTermEnergies *MainWindow::getMeasuredTermTable(int Index)
{
	if (Index >= 0 && Index < numMeasuredTermTables) return measuredTermTables[Index];
	printf("MainWindow::getMeasuredTermTable: the table with Index = %d does not exist!", Index);
	return 0;
}

MeasuredTermEnergies *MainWindow::getMeasuredTermTable(QString FileName, Molecule *molecule)
{
	int i;
	QString S1;
	for (i=0; i < numMeasuredTermTables; i++) 
		if ((S1 = measuredTermTables[i]->getFileName()) == FileName) return measuredTermTables[i]; 
	MeasuredTermEnergies *nTerm = CreateMeasuredTermTable();
	if (nTerm == 0) return 0;
	nTerm->setMolecule(molecule);
	if (!nTerm->readData(FileName))
	{
		nTerm->close();
		numMeasuredTermTables--;
		return 0;
	}
	return nTerm;
}

int MainWindow::getNumDunTables()
{
	//printf("NumDunTables=%d\n", numDunTables);
	return numDunTables;
}

DunTable *MainWindow::getDunTable(int Index)
{
	if (Index >= 0 && Index < numDunTables) return dunTables[Index];
	printf("MainMenu::getDunTable: the table with Index = %d does not exist!", Index);
	return 0;
}

DunTable *MainWindow::getDunTable(QString FileName, Molecule *molecule)
{
	int i;
	QString S1;
	for (i=0; i<numDunTables; i++) if ((S1 = dunTables[i]->getFileName()) == FileName)
			return dunTables[i];
	QString FN = FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numDunTables; i++) if (FN == dunTables[i]->getFName()) return dunTables[i];
	for (i=0; i < numNADunTables; i++) if (NADunTables[i] == FileName) return 0;
	DunTable *Dunham = CreateDunTable();
	if (Dunham == 0) return 0;
	Dunham->setMolecule(molecule);
	if (!Dunham->readData(FileName))
	{
		Dunham->close();
		numDunTables--;
		if (numNADunTables < MaxDunTables) NADunTables[numNADunTables++] = FileName;
		return 0;
	}
	return Dunham;
}

int MainWindow::getNumPotentials()
{
	return numPotentials;
}

Potential *MainWindow::getPotential(int Index)
{
	if (Index >= 0 && Index < numPotentials) return potentials[Index];
	printf("MainWindow::getPotential: the potential with Index = %d does not exist!", Index);
	return 0;
}

Potential *MainWindow::getPotential(QString Name)
{
	int i;
	for (i=0; i < numPotentials; i++) if (potentials[i]->getName() == Name) return potentials[i];
	int m, s, NS, p, NP;
	ElState *S;
	for (m=0; m < numMolecules; m++) for (s=0, NS = molecules[m]->getNumStates(); s < NS; s++)
			for (p=0, NP = (S = molecules[m]->getStateP(s))->getNumPotentials(); p < NP; p++)
				if (S->getPotentialName(p) == Name) return S->getPotential(p);
	printf("MainWindow::getPotential: the potential with the name %s does not exist!", Name.toLatin1().data());
	return 0;
}

Potential *MainWindow::getPotential(QString FileName, Molecule *molecule)
{
	int i;
	for (i=0; i<numPotentials; i++) if (potentials[i]->getFileName() == FileName)
			return potentials[i];
	QString FN = FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numPotentials; i++) if (FN == potentials[i]->getFName()) return potentials[i];
	for (i=0; i < numNAPotentials; i++) if (NAPotentials[i] == FileName) return 0;
	Potential *Pot = CreatePotential();
	if (Pot == 0) return 0;
	Pot->setMolecule(molecule);
	if (!Pot->readData(FileName))
	{
		Pot->close();
		numPotentials--;
		if (numNAPotentials < MaxPotentials) NAPotentials[numNAPotentials++] = FileName;
		return 0;
	}
	return Pot;
}

int MainWindow::getNumFitDataSets()
{
	return numFitDataSets;
}

FitData* MainWindow::getFitData(int Index)
{
	if (Index >= 0 && Index < numFitDataSets) return fitDataSets[Index];
	printf("MainWindow::getFitData: the fit dataset with Index = %d does not exist!", Index);
	return 0;
}

FitData* MainWindow::getFitData(QString Filename, Molecule* molecule)
{
	int i;
	for (i=0; i < numFitDataSets; i++) if (fitDataSets[i]->getFileName() == Filename)
		return fitDataSets[i];
	QString FN = Filename.right(Filename.length() - Filename.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numFitDataSets; i++) if (FN == fitDataSets[i]->getFName()) return fitDataSets[i];
	for (i=0; i < numNAFitDataSets; i++) if (NAFitDataSets[i] == Filename) return 0;
	FitData *Set = CreateFitData();
	if (Set == 0) return 0;
	Set->setMolecule(molecule);
	if (!Set->readData(Filename))
	{
		Set->close();
		numFitDataSets--;
		if (numFitDataSets < MaxFitDataSets) NAFitDataSets[numNAFitDataSets++] = Filename;
		return 0;
	}
	return Set;
}

int MainWindow::getNumLineTables()
{
	return numLineTables;
}

LineTable *MainWindow::getLineTable(int Index)
{
	//printf("lineTables[%d]=%d\n", Index, lineTables[Index]);
	if (Index >= 0 && Index < numLineTables) return lineTables[Index];
	printf("MainWindow::getLineTable: the table with Index = %d does not exist!\n", Index);
	return 0;
}

LineTable *MainWindow::getLineTable(QString FileName, Molecule *molecule)
{
	printf("MainWindow::getLineTable\n");
	int i;
	QString S1;
	for (i=0; i < numLineTables; i++) if ((S1 = lineTables[i]->getFileName()) == FileName) 
			return lineTables[i];
	QString FN = FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numLineTables; i++) if (FN == lineTables[i]->getFName()) return lineTables[i];
	for (i=0; i < numNALineTables; i++) if (NALineTables[i] == FileName) return 0;
	LineTable *Table = CreateLineTable();
	printf("Nach CreateLineTable\n");
	if (Table == 0) return 0;
	Table->setMolecule(molecule);
	printf("Nach setMolecule\n");
	if (!Table->readData(FileName))
	{
		Table->close();
		numLineTables--;
		if (numNALineTables < MaxLineTables) NALineTables[numNALineTables++] = FileName;
		return 0;
	}
	printf("Ende MainWindow::getLineTable\n");
	return Table;
}

int MainWindow::getNumFCFTables()
{
	return numFCFTables;
}

FCFTab* MainWindow::getFCFTable(int Index)
{
	if (Index >= 0 && Index < numFCFTables) return fcfTables[Index];
	printf("MainWindow::getFCFTable: the table with Index = %d does not exist!\n", Index);
	return 0;
}

FCFTab* MainWindow::getFCFTable(QString Filename, Molecule* molecule)
{
	int i;
	QString S1;
	for (i=0; i < numFCFTables; i++) if ((S1 = fcfTables[i]->getFileName()) == Filename)
		return fcfTables[i];
	QString FN = Filename.right(Filename.length() - Filename.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numFCFTables; i++) if (FN == fcfTables[i]->getFName()) return fcfTables[i];
	for (i=0; i < numNAFCFTables; i++) if (NAFCFTables[i] == Filename) return 0;
	FCFTab *Table = CreateFCFTable();
	if (Table == 0) return 0;
	Table->setMolecule(molecule);
	if (!Table->readData(Filename))
	{
		Table->close();
		numFCFTables--;
		if (numNAFCFTables < MaxFCFTables) NAFCFTables[numNAFCFTables++] = Filename;
		return 0;
	}
	return Table;
}

int MainWindow::getNumSpectra()
{
	return numSpectra;
}

Spektrum *MainWindow::getSpectrum(int Index)
{
	if (Index >= 0 && Index < numSpectra) return spectra[Index];
	printf("MainWindow::getSpectrum: the spectrum with the Index %d does not exist!", Index);
	return 0;
}

Spektrum *MainWindow::getSpectrum(QString FileName)
{
	int i;
	QString S1;
	for (i=0; i < numSpectra; i++) if ((S1 = spectra[i]->getFileName()) == FileName) return spectra[i];
	QString FN = FileName.right(FileName.length() - FileName.lastIndexOf(QRegExp("[\\/]")) - 1);
	for (i=0; i < numSpectra; i++) if (FN == spectra[i]->getFName()) return spectra[i];
	if (numSpectra > 0) printf("FN=%s, FName=%s\n", FN.toLatin1().data(), spectra[0]->getFName().toLatin1().data());
	for (i=0; i < numNASpectra; i++) if (NASpectra[i] == FileName) return 0;
	Spektrum *spectrum = CreateSpectrum();
	if (spectrum == 0) return 0;
	if (!spectrum->readData(FileName))
	{
		spectrum->close();
		numSpectra--;
		if (numNASpectra < MaxSpectra) NASpectra[numNASpectra++] = FileName;
		return 0;
	}
	return spectrum;
}

void MainWindow::quit()
{
	if (checkSaved()) emit quitApp();
}

void MainWindow::closeEvent(QCloseEvent *E)
{
	if (checkSaved()) E->accept();
	else E->ignore();
}

bool MainWindow::checkNotOpened(QString File, QString Filter)
{
	int i, n = getNumMDIChilds();
	MDIChild *C;
	for (i=0; i<n; i++) if ((C = getMDIChild(i))->getFileName() == File) 
		if (Filter.isEmpty() || Filter == C->getFilter())
	{
		if (QMessageBox::question(this, "MolSpektAnalysis",
					"The file '" + File + "' is already opened! Do you really want to open it again?",
					QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes)
			return true;
		return false;
	}
	return true;
}

int MainWindow::getNDirRep()
{
	return numDirRep;
}

void MainWindow::addDirRep(QString DR, int NC)
{
	if (numDirRep == MaxDirRep) return;
	DirRep[numDirRep] = DR;
	DirRepNC[numDirRep++] = NC;
}

void MainWindow::getDirRep(int &N, QString &DR, int &NC)
{
	DR = DirRep[N];
	NC = DirRepNC[N];
}

int MainWindow::getNumMDIChilds()
{
	return numAtoms + numTermTables + numMeasuredTermTables + numDunTables 
			+ numPotentials + numFitDataSets + numLineTables + numFCFTables 
			+ numMolecules + numSpectra + numSpectLists + numSpectSimulations;
}

MDIChild *MainWindow::getMDIChild(int i)
{
	if (i < 0) 
	{
		printf("Error: MainWindow::getMDIChild(i): i=%d is below 0!", i);
		return 0;
	}
	if (i < numAtoms) return atoms[i];
	int j = numAtoms + numTermTables, k;
	if (i < j) return termTables[i - numAtoms];
	if (i < (k = j + numDunTables)) return dunTables[i - j];
	if (i < (j = k + numPotentials)) return potentials[i - k];
	if (i < (k = j + numLineTables)) return lineTables[i - j];
	if (i < (j = k + numFCFTables)) return fcfTables[i-k];
	if (i < (k = j + numFitDataSets)) return fitDataSets[i-j];
	if (i < (j = k + numMolecules)) return molecules[i - k];
	if (i < (k = j + numSpectra)) return spectra[i-j];
	if (i < (j = k + numSpectLists)) return spectLists[i-k];
	if (i < (k = j + numSpectSimulations)) return spectSimulations[i-j];
	if (i < k + numMeasuredTermTables) return measuredTermTables[i-k];
	printf("Error: MainWindow::getMDIChild(i): i=%d is greater than the number of MDIChilds!", i);
	return 0;
}

QString MainWindow::getChildName(QString S, MDIChild *T)
{
	//printf("MainWindow::getChildname(S=%s, T=%s)\n", S.ascii(), T.ascii());
	int i, j, k, N = getNumMDIChilds(), I;
	QString R = S;
	MDIChild::Type type = T->getType();
	bool U = false;
	MDIChild *M;
	Transition * Tr;
	ElState * St;
	for (i = S.size() - 1; (i >= 0 ? S[i] >= '0' && S[i] <= '9' : false); i--) ; 
	i = S.size() - i - 1;
	if (i>0)
	{
		I = S.right(i).toInt();
		S = S.left(S.size() - i);
	}
	else I=0;
	while (!U)
	{
		U = true;
		for (i=0; i<N; i++) 
		{
			M = getMDIChild(i);
			if (M->getName() == R && M->getType() == type && M != T) U = false;
			//printf("i=%d, N=%s, T=%s\n", i, M->getName().ascii(), M->getType().ascii());
		}
		if (U)
		{
			if(type == MDIChild::PotData)
            {
                for (i=0; i < numMolecules; i++) for (j=0; j < molecules[i]->getNumStates(); j++)
					for (k=0, St = molecules[i]->getStateP(j); k < St->getNumPotentials(); k++) 
						if (R == St->getPotentialName(k) 
							&& T->getFileName() != St->getPotentialFileName(k)) U = false;
            }
			else if (type == MDIChild::DunhamTable)
            {
                for (i=0; i < numMolecules; i++) for (j=0; j < molecules[i]->getNumStates(); j++)
					for (k=0, St = molecules[i]->getStateP(j); k < St->getNumDunTables(); k++)
						if (R == St->getDunTableName(k)
							&& T->getFileName() != St->getDunTableFileName(k)) U = false;
            }
			else if (type == MDIChild::TermEnergyTable)
            {
                for (i=0; i < numMolecules; i++) for (j=0; j < molecules[i]->getNumStates(); j++)
					for (k=0, St = molecules[i]->getStateP(j); k < St->getNumTermTables(); k++)
						if (R == St->getDunTableName(k)
							&& T->getFileName() != St->getDunTableFileName(k)) U = false;
            }
			else if (type == MDIChild::FitDataSet)
            {
                for (i=0; i < numMolecules; i++) for (j=0; j < molecules[i]->getNumStates(); j++)
					for (k=0, St = molecules[i]->getStateP(j); k < St->getNumFitDataSets(); k++)
						if (R == St->getFitDataName(k)
							&& T->getFileName() != St->getFitDataFileName(k)) U = false;
            }
			else if (type == MDIChild::LineTab)
            {
                for (i=0; i < numMolecules; i++) for (j=0; j < molecules[i]->getNumTransitions(); j++)
					for (k=0, Tr = molecules[i]->getTransitionP(j); k < Tr->getNumLineTables(); k++)
						if (R == Tr->getLineTableName(k)
							&& T->getFileName() != Tr->getLineTableFileName(k)) U = false;
            }
			else if (type == MDIChild::FranckCondonTable) for (i=0; i < numFCFTables; i++)
				for (j=0; j < molecules[i]->getNumTransitions(); j++)
					for (k=0, Tr = molecules[i]->getTransitionP(j); k < Tr->getNumFCFTables(); k++)
						if (R == Tr->getFCFTableName(k)
							&& T->getFileName() != Tr->getFCFTableFileName(k)) U = false;
		}
		if (!U) R = S + QString::number(I++);
	}
	return R;
}

void MainWindow::setAskForQuit(MDIChild* Window)
{
	if (NumWindowsToAskForQuit > 0)
	{
		int n;
		MDIChild **NW = new MDIChild*[NumWindowsToAskForQuit + 1];
		for (n=0; n < NumWindowsToAskForQuit; n++) NW[n] = windowsToAskForQuit[n];
		delete[] windowsToAskForQuit;
		windowsToAskForQuit = NW;
	}
	else windowsToAskForQuit = new MDIChild*[1];
	windowsToAskForQuit[NumWindowsToAskForQuit++] = Window;
}

void MainWindow::setActive(MDIChild *W)
{
	QList<QMdiSubWindow*> subWindowList = workspace->subWindowList();
	for (QList<QMdiSubWindow*>::iterator it = subWindowList.begin(); it != subWindowList.end(); ++it)
		if ((*it)->widget() == W)
	{
		workspace->setActiveSubWindow(*it);
		return;
	}
}

bool MainWindow::checkSaved()
{
	int i, n = getNumMDIChilds();
	for (i=0; i < NumWindowsToAskForQuit; i++) if (!windowsToAskForQuit[i]->askForQuit()) return false;
	for (i=0; i<n; i++) if (!checkSaved(getMDIChild(i))) return false;
	return true;
}

bool MainWindow::checkSaved(MDIChild *C)
{
	if (!C->isSaved())
	{
		int R = QMessageBox::question(this, "MolSpektAnalysis", 
						"The data of the " + C->getTypeString() + " " + C->getName() + 
						"\nhas not been saved. Do you want\nto save it now?",
						QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel,
						QMessageBox::Yes);
		if (R == QMessageBox::Cancel) return false;
		if (R == QMessageBox::Yes) return C->writeData();
	}
	return true;
}

void MainWindow::newAtom()
{
	Atom *nAtom;
	if ((nAtom = CreateAtom()) == 0) return;
	printf("Vor Atom::show\n");
	nAtom->show();
}

void MainWindow::newMolecule()
{
	Molecule *nMol;
	if ((nMol = CreateMolecule()) == 0) return;
	nMol->show();
}

void MainWindow::newTermTable()
{
	TermTable *nTerm;
	if ((nTerm = CreateTermTable()) == 0) return;
	nTerm->show();
}

void MainWindow::newMeasuredTermTable()
{
	MeasuredTermEnergies *nTerm;
	if ((nTerm = CreateMeasuredTermTable()) == 0) return;
	nTerm->show();
}

void MainWindow::newDunTable()
{
	DunTable *Dunham;
	if ((Dunham = CreateDunTable()) == 0) return;
	Dunham->show();
}

void MainWindow::newPotential()
{
	Potential *Pot;
	if ((Pot = CreatePotential()) == 0) return;
	Pot->show();
}

void MainWindow::newFitData()
{
	FitData *Set;
	if ((Set = CreateFitData()) == 0) return;
	Set->show();
}

void MainWindow::newLineTable()
{
	LineTable *Table;
	if ((Table = CreateLineTable()) == 0) return;
	Table->show();
}

void MainWindow::newFCFTable()
{
    calcFCFDiag(NumFCF_WFPoints, false);
}

void MainWindow::calcFCFDiag(int NumWFPoints, bool update)
{
	FCFTab *Table = 0, *F = dynamic_cast<FCFTab*> (workspace->activeSubWindow()->widget());
	int i, i1 = 0, i2 = 1, N;
	if (update) Table = F;
	else if (numFCFTables == MaxFCFTables) return;
	Molecule *M = (F==0 ? dynamic_cast<Molecule*> (workspace->activeSubWindow()->widget()) : F->getMolecule());
	if (M==0)
	{
		printf("MainWindow::calcFCFDiag: error: There is no molecule selected to calculate FCF for!");
		return;
	}
	ElState *S, *St1 = (F!=0 ? F->getState1() : 0), *St2 = (F!=0 ? F->getState2() : 0);
	QDialog *D = new QDialog(this);
	QComboBox *S1 = new QComboBox(D), *S2 = new QComboBox(D);
	for (i=0, N = M->getNumStates(); i<N; i++) if ((S = M->getStateP(i))->getPotential() != 0)
	{
		S1->addItem(S->getName());
		S2->addItem(S->getName());
		if (F != 0)
		{
			if (S == St1) i1 = S1->count() - 1;
			if (S == St2) i2 = S2->count() - 1;
		}
	}
	S1->setCurrentIndex(i1);
	if (i2 < S2->count()) S2->setCurrentIndex(i2);
	if (update)
	{
		S1->setEnabled(false);
		S2->setEnabled(false);
	}
	else
	{
		S1->setEditable(false);
		S2->setEditable(false);
	}
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("Upper el. state:", D), 0, 0);
	L->addWidget(S1, 0, 1);
	L->addWidget(new QLabel("Lower el. state:", D), 0, 2);
	L->addWidget(S2, 0, 3);
	QGridLayout *L1 = new QGridLayout;
	L1->addWidget(new QLabel("Max v':", D), 0, 0);
	QLineEdit *Mvs = new QLineEdit("100", D);
	L1->addWidget(Mvs, 0, 1);
	L1->addWidget(new QLabel("Max v'':", D), 0, 2);
	QLineEdit *Mvss = new QLineEdit("100", D);
	L1->addWidget(Mvss, 0, 3);
	L1->addWidget(new QLabel("Max J:", D), 0, 4);
	QLineEdit *MJ = new QLineEdit("350", D);
	L1->addWidget(MJ, 0, 5);
	L->addLayout(L1, 1, 0, 1, 6);
	L->setRowMinimumHeight(2, 20);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	L->addWidget(OK, 3, 0);
	L->addWidget(Cancel, 3, 3);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		for (i=0; i<N; i++)
		{
			if (M->getState(i) == S1->currentText()) St1 = M->getStateP(i);
			if (M->getState(i) == S2->currentText()) St2 = M->getStateP(i);
		}
		if (Table == 0) Table = CreateFCFTable(M);
        Table->calcFCF(St1, St2, Mvs->text().toInt(), Mvss->text().toInt(), MJ->text().toInt(), NumWFPoints);
		Table->show();
	}
	delete D;
}

void MainWindow::newSpectrum()
{
	Spektrum *Spectrum;
	if ((Spectrum = CreateSpectrum()) == 0) return;
	Spectrum->show();
}

void MainWindow::newSpectList()
{
	SpectList *S = CreateSpectList();
	if (S != 0) S->show();
}

Atom *MainWindow::CreateAtom()
{
	if (numAtoms == MaxAtoms)
	{
		QMessageBox::information(this, tr("MolSpektAnalysis"), 
								 tr("The maximum amount of atoms has been reached!"));
		return 0;
	}
	Atom *nAtom = atoms[numAtoms] = new Atom;
	ShowAtom[numAtoms] = new ShowAction("New atom" + QString::number(numAtoms), this);
	ShowAtom[numAtoms]->setStatusTip("Show the window with the data of this atom");
	connect(ShowAtom[numAtoms], SIGNAL(triggered()), nAtom, SLOT(show()));
	//connect(nAtom, SIGNAL(hidden(bool)), ShowAtom[numAtoms], SLOT(setVisible(bool)));
	//connect(nAtom, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(nAtom, SIGNAL(nameChanged(QString)), 
			ShowAtom[numAtoms], SLOT(setName(QString)));
	QMdiSubWindow *subWindow = new QMdiSubWindow;
	subWindow->setWidget(nAtom);
	workspace->addSubWindow(subWindow);
    connect(nAtom, SIGNAL(closeThis()), subWindow, SLOT(close()));
	ShowAtoms->addAction(ShowAtom[numAtoms++]);
	ShowAtoms->setEnabled(true);
	ShowMenu->setEnabled(true);
	enableShowMenu(true);
	return nAtom;
}

Molecule *MainWindow::CreateMolecule()
{
	if (numMolecules == MaxMolecules)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "The maximum amount of molecules has been reached!");
		return 0;
	}
	Molecule *Mol = new Molecule(this);
	ShowMolecule[numMolecules] = new ShowAction("New molecule" + QString::number(numMolecules), this);
	ShowMolecule[numMolecules]->setStatusTip("Show the window with the data of this molecule");
	connect(ShowMolecule[numMolecules], SIGNAL(triggered()), Mol, SLOT(show()));
	//connect(Mol, SIGNAL(hidden(bool)), ShowMolecule[numMolecules], SLOT(setVisible(bool)));
	//connect(Mol, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(Mol, SIGNAL(nameChanged(QString)),
			ShowMolecule[numMolecules], SLOT(setName(QString)));
	connect(Mol, SIGNAL(nameChanged(QString)), this, SIGNAL(MoleculesChanged()));
	ShowMolecules->addAction(ShowMolecule[numMolecules]);
	ShowMolecules->setEnabled(true);
	ShowMenu->setEnabled(true);
	molecules[numMolecules++] = Mol;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(Mol);
	workspace->addSubWindow(subWindow);
	emit MoleculesChanged();
	enableShowMenu(true);
	printf("MainWindow::CreateMolecule: End\n");
	return Mol;
}

QTableWidget *MainWindow::createTable()
{
	QTableWidget *Tab = new QTableWidget;
	workspace->addSubWindow(Tab);
	return Tab;
}

TableWindow *MainWindow::createTableWindow()
{
	TableWindow *Tab = new TableWindow(MDIChild::TextTable1, this);
	workspace->addSubWindow(Tab);
	return Tab;
}

TermTable *MainWindow::CreateTermTable()
{
	if (numTermTables == MaxTermTables)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "The maximum amount of term energy tables has been reached!");
		return 0;
	}
	TermTable *nTerm = new TermTable(this);
	ShowTermTable[numTermTables] = new ShowAction("Show term energy table " 
			+ QString::number(numTermTables), this);
	ShowTermTable[numTermTables]->setStatusTip("Show this term energy table");
	connect(ShowTermTable[numTermTables], SIGNAL(triggered()), nTerm, SLOT(show()));
	//connect(nTerm, SIGNAL(hidden(bool)), ShowTermTable[numTermTables], SLOT(setVisible(bool)));
	//connect(nTerm, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(nTerm, SIGNAL(nameChanged(QString)), 
			ShowTermTable[numTermTables], SLOT(setName(QString)));
	connect(nTerm, SIGNAL(nameChanged(QString)), this, SIGNAL(TermTablesChanged()));
	ShowTermTables->addAction(ShowTermTable[numTermTables]);
	ShowTermTables->setEnabled(true);
	ShowMenu->setEnabled(true);
	termTables[numTermTables++] = nTerm;
	QMdiSubWindow *subWindow = new QMdiSubWindow;
	subWindow->setWidget(nTerm);
	workspace->addSubWindow(subWindow);
	//printf("TermTables[%d]=%d, nTerm=%d\n", numTermTables -1, termTables[numTermTables-1], nTerm);
	return nTerm;
}

MeasuredTermEnergies *MainWindow::CreateMeasuredTermTable()
{
	if (numMeasuredTermTables == MaxTermTables)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of measured term energy tables has been reached!");
		return 0;
	}
	MeasuredTermEnergies *nTerm = new MeasuredTermEnergies(this);
	ShowMeasuredTermTable[numMeasuredTermTables] = new ShowAction("Show measured term energy table "
			+ QString::number(numMeasuredTermTables), this);
	ShowMeasuredTermTable[numMeasuredTermTables]->setStatusTip(
			"Show this table of directly measured term energies");
	connect(ShowMeasuredTermTable[numMeasuredTermTables], SIGNAL(triggered()), nTerm, SLOT(show()));
	connect(nTerm, SIGNAL(nameChanged(QString)),
			ShowMeasuredTermTable[numMeasuredTermTables], SLOT(setName(QString)));
	ShowMeasuredTermTables->addAction(ShowMeasuredTermTable[numMeasuredTermTables]);
	ShowMeasuredTermTables->setEnabled(true);
	ShowMenu->setEnabled(true);
	measuredTermTables[numMeasuredTermTables++] = nTerm;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(nTerm);
	workspace->addSubWindow(subWindow);
	return nTerm;
}

DunTable *MainWindow::CreateDunTable()
{
	if (numDunTables == MaxDunTables)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of Dunham coefficient tables has been reached!");
		return 0;
	}
	DunTable *Dunham = new DunTable(this);
	ShowDunTable[numDunTables] =  new ShowAction("Show Dunham coefficient table "
			+ QString::number(numDunTables), this);
	ShowDunTable[numDunTables]->setStatusTip("Show this table of Dunham coefficients");
	connect(ShowDunTable[numDunTables], SIGNAL(triggered()), Dunham, SLOT(show()));
	//connect(Dunham, SIGNAL(hidden(bool)), ShowDunTable[numDunTables], SLOT(setVisible(bool)));
	//connect(Dunham, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(Dunham, SIGNAL(nameChanged(QString)),
			ShowDunTable[numDunTables], SLOT(setName(QString)));
	connect(Dunham, SIGNAL(nameChanged(QString)), this, SIGNAL(DunTablesChanged()));
	ShowDunTables->addAction(ShowDunTable[numDunTables]);
	ShowDunTables->setEnabled(true);
	ShowMenu->setEnabled(true);
	dunTables[numDunTables++] = Dunham;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(Dunham);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	return Dunham;
}

Potential *MainWindow::CreatePotential(int ThreadNum)
{
	if (numPotentials == MaxPotentials)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of potentials has been reached!");
		return 0;
	}
	Potential *Pot = new Potential(this, 0, ThreadNum);
	ShowPotential[numPotentials] = new ShowAction("Show potential table " 
			+ QString::number(numPotentials), this);
	ShowPotential[numPotentials]->setStatusTip("Show this table of potential coefficients");
	connect(ShowPotential[numPotentials], SIGNAL(triggered()), Pot, SLOT(show()));
	//connect(Pot, SIGNAL(hidden(bool)), ShowPotential[numPotentials], SLOT(setVisible(bool)));
	//connect(Pot, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(Pot, SIGNAL(nameChanged(QString)),
			ShowPotential[numPotentials], SLOT(setName(QString)));
	connect(Pot, SIGNAL(nameChanged(QString)), this, SIGNAL(PotentialsChanged()));
	ShowPotentials->addAction(ShowPotential[numPotentials]);
	ShowPotentials->setEnabled(true);
	ShowMenu->setEnabled(true);
	potentials[numPotentials++] = Pot;
	QMdiSubWindow *subWindow = new QMdiSubWindow;
	subWindow->setWidget(Pot);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	return Pot;
}

FitData* MainWindow::CreateFitData()
{
	if (numFitDataSets == MaxFitDataSets)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of fit datasets has been reached!");
		return 0;
	}
	FitData *Set = new FitData(0, this);
	ShowFitData[numFitDataSets] = new ShowAction("Show fit dataset "
			+ QString::number(numFitDataSets), this);
	ShowFitData[numFitDataSets]->setStatusTip("Show this table with fit data");
	connect(ShowFitData[numFitDataSets], SIGNAL(triggered()), Set, SLOT(show()));
	connect(Set, SIGNAL(nameChanged(QString)), ShowFitData[numFitDataSets], SLOT(setName(QString)));
	connect(Set, SIGNAL(nameChanged(QString)), this, SIGNAL(FitDataChanged()));
	ShowFitDataSets->addAction(ShowFitData[numFitDataSets]);
	ShowFitDataSets->setEnabled(true);
	ShowMenu->setEnabled(true);
	fitDataSets[numFitDataSets++] = Set;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(Set);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	return Set;
}

LineTable *MainWindow::CreateLineTable()
{
	//printf("MainWindow::CreateLineTable\n");
	if (numLineTables == MaxLineTables)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of line tables has already been reached!");
		return 0;
	}
	LineTable *Table = new LineTable(this);
	//printf("Table=%d\n", Table);
	ShowLineTable[numLineTables] = new ShowAction("Show line table " 
								 + QString::number(numLineTables), this);
	ShowLineTable[numLineTables]->setStatusTip("Shows this line table");
	//printf("Vor connect\n");
	connect(ShowLineTable[numLineTables], SIGNAL(triggered()), Table, SLOT(show()));
	//connect(Table, SIGNAL(hidden(bool)), ShowLineTable[numLineTables], SLOT(setVisible(bool)));
	//connect(Table, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(Table, SIGNAL(nameChanged(QString)), ShowLineTable[numLineTables], SLOT(setName(QString)));
	connect(Table, SIGNAL(nameChanged(QString)), this, SIGNAL(LineTablesChanged()));
	connect(this, SIGNAL(SpectrumChanged(Spektrum*)), Table, SLOT(updateMarker(Spektrum*)));
	connect(Table, SIGNAL(DataChanged()), this, SIGNAL(LineTableChanged()));
	//printf("Nach connect\n");
	ShowLineTables->addAction(ShowLineTable[numLineTables]);
	ShowLineTables->setEnabled(true);
	ShowMenu->setEnabled(true);
	lineTables[numLineTables++] = Table;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(Table);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	//printf("Ende CreateLineTable\n");
	return Table;
}

LineTable *MainWindow::getUnassignedLineTable()
{
	int n;
	for (n=0; (n < numLineTables ? lineTables[n]->isAssigned() : false); n++) ;
	if (n == numLineTables) CreateLineTable();
	if (n < numLineTables) return lineTables[n];
	return 0;
}

void MainWindow::loadLineTables()
{
	int l, nl, t, nt, m;
	Transition *T;
	for (m=0; m < numMolecules; m++) for (t=0, nt = molecules[m]->getNumTransitions(); t < nt; t++)
			for (l=0, nl = (T = molecules[m]->getTransitionP(t))->getNumLineTables(); l < nl; l++)
				T->getLineTable(l);
}

void MainWindow::LineTableSaved(ElState *UState, ElState *LState)
{
	int n;
	for (n=0; n < numSpectra; n++) spectra[n]->LineTableSaved(UState, LState);
}

FCFTab* MainWindow::CreateFCFTable(Molecule* Mol)
{
	if (numFCFTables == MaxFCFTables)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of FCF tables has already been reached!");
		return 0;
	}
	FCFTab *Table = new FCFTab(this, Mol);
	ShowFCFTable[numFCFTables] = new ShowAction("Show FCF table "
								+ QString::number(numFCFTables), this);
	ShowFCFTable[numFCFTables]->setStatusTip("Show this table of Franck-Condon factors");
	connect(ShowFCFTable[numFCFTables], SIGNAL(triggered()), Table, SLOT(show()));
	connect(Table, SIGNAL(nameChanged(QString)), 
			ShowFCFTable[numFCFTables], SLOT(setName(QString)));
	connect(Table, SIGNAL(propertiesChanged()), this, SIGNAL(FCFTablesChanged()));
	ShowFCFTables->addAction(ShowFCFTable[numFCFTables]);
	ShowFCFTables->setEnabled(true);
	ShowMenu->setEnabled(true);
	fcfTables[numFCFTables++] = Table;
	QMdiSubWindow *subWindow = new QMdiSubWindow;
	subWindow->setWidget(Table);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	return Table;
}

Spektrum *MainWindow::CreateSpectrum()
{
	if (numSpectra == MaxSpectra)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
								 "The maximum amount of spectra has already bean reached!");
		return 0;
	}
	Spektrum *Spectrum = new Spektrum(this);
	ShowSpectrum[numSpectra] = new ShowAction("Show spectrum " 
							 + QString::number(numSpectra), this);
	ShowSpectrum[numSpectra]->setStatusTip("Shows this spectra");
	connect(ShowSpectrum[numSpectra], SIGNAL(triggered()), Spectrum, SLOT(show()));
	//connect(Spectrum, SIGNAL(hidden(bool)), ShowSpectrum[numSpectra], SLOT(setVisible(bool)));
	//connect(Spectrum, SIGNAL(hidden(bool)), this, SLOT(enableShowMenu(bool)));
	connect(Spectrum, SIGNAL(nameChanged(QString)), ShowSpectrum[numSpectra], SLOT(setName(QString)));
	connect(Spectrum, SIGNAL(SpectrumChanged(Spektrum*)), this, SIGNAL(SpectrumChanged(Spektrum*)));
	Spectrum->setName("New spectrum");
	ShowSpectra->addAction(ShowSpectrum[numSpectra]);
	ShowSpectra->setEnabled(true);
	ShowMenu->setEnabled(true);
	spectra[numSpectra++] = Spectrum;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(Spectrum);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	emit newSpectrum(Spectrum);
	return Spectrum;
}

SpectList* MainWindow::CreateSpectList()
{
	if (numSpectLists == MaxSpectra)
	{
		QMessageBox::information(this, "MolSpectAnalysis",
						"The maximum amount of spectLists has already been reached!");
		return 0;
	}
	SpectList *SL = new SpectList(this);
	ShowSpectList[numSpectLists] = new ShowAction("Show spectList " 
												+ QString::number(numSpectLists), this);
	ShowSpectList[numSpectLists]->setStatusTip("Show this list of spectra");
	connect(ShowSpectList[numSpectLists], SIGNAL(triggered()), SL, SLOT(show()));
	SL->setName("New SpectList");
	ShowSpectLists->addAction(ShowSpectList[numSpectLists]);
	spectLists[numSpectLists++] = SL;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(SL);
	workspace->addSubWindow(subWindow);
	enableShowMenu(true);
	return SL;
}

SpectSimulation *MainWindow::CreateSpectSimulation()
{
	if (numSpectSimulations == MaxSpectSimulations)
	{
		QMessageBox::information(this, "MolSpektAnalysis",
					"The maximum amount of windows for simulated spectra has already been reached!");
		return 0;
	}
	spectSimulations[numSpectSimulations] = new SpectSimulation(this);
	//printf("spectSimulations[%d]=%d\n", numSpectSimulations, spectSimulations[numSpectSimulations]);
	connect(this, SIGNAL(MoleculesChanged()), 
			spectSimulations[numSpectSimulations], SLOT(RefreshMolecules()));
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(spectSimulations[numSpectSimulations]);
	workspace->addSubWindow(subWindow);
	return spectSimulations[numSpectSimulations++];
}

ElState *MainWindow::CreateElState()
{
	ElState *S = new ElState;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(S);
	workspace->addSubWindow(subWindow);
	return S;
}

Transition *MainWindow::CreateTransition()
{
	Transition *T = new Transition;
	QMdiSubWindow* subWindow = new QMdiSubWindow;
	subWindow->setWidget(T);
	workspace->addSubWindow(subWindow);
	return T;
}

DiagWindow *MainWindow::CreateDiagWindow(MDIChild::Type type)
{
	DiagWindow *W = new DiagWindow(type, this);
	workspace->addSubWindow(W);
	return W;
}

void MainWindow::enableShowMenu(bool enable)
{
	//printf("MainWindow::enableShowMenu\n");
	if (enable)
	{
		if (numAtoms > 0 && !ShowAtoms->isEnabled()) ShowAtoms->setEnabled(true);
		if (numMolecules > 0 && !ShowMolecules->isEnabled()) ShowMolecules->setEnabled(true);
		if (numTermTables > 0 && !ShowTermTables->isEnabled()) ShowTermTables->setEnabled(true);
		if (numMeasuredTermTables > 0 && !ShowMeasuredTermTables->isEnabled()) 
			ShowMeasuredTermTables->setEnabled(true);
		if (numDunTables > 0 && !ShowDunTables->isEnabled()) ShowDunTables->setEnabled(true);
		if (numPotentials > 0 && !ShowPotentials->isEnabled()) ShowPotentials->setEnabled(true);
		if (numFitDataSets > 0 && !ShowFitDataSets->isEnabled()) ShowFitDataSets->setEnabled(true);
		if (numLineTables > 0 && !ShowLineTables->isEnabled()) ShowLineTables->setEnabled(true);
		if (numFCFTables > 0 && !ShowFCFTables->isEnabled()) ShowFCFTables->setEnabled(true);
		if (numSpectra > 0 && !ShowSpectra->isEnabled()) ShowSpectra->setEnabled(true);
		if (numSpectLists > 0 && !ShowSpectLists->isEnabled())
			ShowSpectLists->setEnabled(true);
		if ((numAtoms > 0 || numMolecules > 0 || numTermTables > 0 || numDunTables > 0 
			   || numPotentials > 0 || numLineTables > 0 || numFCFTables > 0 
			   || numSpectra > 0 || numSpectLists > 0) 
			   && !ShowMenu->isEnabled()) ShowMenu->setEnabled(true);
	}
	else
	{
		if (numAtoms == 0 && ShowAtoms->isEnabled()) ShowAtoms->setEnabled(false);
		if (numMolecules == 0 && ShowMolecules->isEnabled()) ShowMolecules->setEnabled(false);
		if (numTermTables == 0 && ShowTermTables->isEnabled()) ShowTermTables->setEnabled(false);
		if (numMeasuredTermTables == 0 && ShowMeasuredTermTables->isEnabled()) 
			ShowMeasuredTermTables->setEnabled(false);
		if (numDunTables == 0 && ShowDunTables->isEnabled()) ShowDunTables->setEnabled(false);
		if (numPotentials == 0 && ShowPotentials->isEnabled()) ShowPotentials->setEnabled(false);
		if (numFitDataSets == 0 && ShowFitDataSets->isEnabled()) ShowFitDataSets->setEnabled(false);
		if (numLineTables == 0 && ShowLineTables->isEnabled()) ShowLineTables->setEnabled(false);
		if (numFCFTables == 0 && ShowFCFTables->isEnabled()) ShowFCFTables->setEnabled(false);
		if (numSpectra == 0 && ShowSpectra->isEnabled()) ShowSpectra->setEnabled(false);
		if (numSpectLists == 0 && ShowSpectLists->isEnabled()) 
			ShowSpectLists->setEnabled(false);
		if (numAtoms == 0 && numMolecules == 0 && numTermTables == 0 && numFitDataSets == 0
			      && numMeasuredTermTables == 0 && numDunTables == 0 && numPotentials == 0 
			      && numLineTables == 0 && numFCFTables == 0 && numSpectra == 0 
				  && ShowMenu->isEnabled()) ShowMenu->setEnabled(false);
	}
	//printf("Ende enableShowMenu\n");
}

void MainWindow::open()
{
	QString Filter;
	QString Filename = QFileDialog::getOpenFileName(this, "Choose a file to open", Dir, 
		"Molecules (*.mol);;Atoms (*.atom);;Term energy tables (*.term);;Measured term tables (*.mterm);;Dunham coefficients (*.dun);;Potentials (*.pot);;Fit datasets (*.fdat);;Line tables (*.lines);;FCF tables (*.fcf);;Spectrum (*.spect);;SpectList (*.slist);;Simulated spectrum (*.sspk);;Text table (*.dat *.txt);;Log file (*.log)",
			 &Filter);
	int n, m;
	if (Filename.isEmpty()) return;
	for (n=m=0; m!=-1; n=m+1) m = Filename.indexOf(DIRSEP, n);
	Dir = Filename.left(n);
	if (!checkNotOpened(Filename, Filter)) return;
	//printf("Filename = %s \n", Filename.ascii());
	if (Filter == "Atoms (*.atom)")
	{
		Atom *nAtom;
		if ((nAtom = CreateAtom()) == 0) return;
		if (nAtom->readData(Filename)) nAtom->show();
	}
	else if (Filter == "Molecules (*.mol)")
	{
		Molecule *nMol = CreateMolecule();
		if (nMol == 0) return;
		if (nMol->readData(Filename)) nMol->show();
	}
	else if (Filter == "Term energy tables (*.term)")
	{
		TermTable *nTerm = CreateTermTable();
		if (nTerm == 0) return;
		if (nTerm->readData(Filename)) nTerm->show();
		//printf("termTables[0]=%d\n", termTables[0]);
	}
	else if (Filter == "Measured term tables (*.mterm)")
	{
		MeasuredTermEnergies *nTerm = CreateMeasuredTermTable();
		if (nTerm == 0) return;
		if (nTerm->readData(Filename)) nTerm->show();
	}
	else if (Filter == "Dunham coefficients (*.dun)")
	{
		DunTable *Dunham = CreateDunTable();
		if (Dunham == 0) return;
		if (Dunham->readData(Filename)) Dunham->show();
	}
	else if (Filter == "Potentials (*.pot)")
	{
		Potential *Pot = CreatePotential();
		if (Pot == 0) return;
		//printf("Vor readData\n");
		if (Pot->readData(Filename)) 
		{
			//printf("Nach readData\n");
			Pot->show();
			//printf("Nach Show\n");
		}
	}
	else if (Filter == "Fit datasets (*.fdat)")
	{
		FitData *Set = CreateFitData();
		if (Set == 0) return;
		if (Set->readData(Filename)) Set->show();
	}
	else if (Filter == "Line tables (*.lines)")
	{
		LineTable *Table = CreateLineTable();
		if (Table == 0) return;
		if (Table->readData(Filename)) Table->show();
	}
	else if (Filter == "FCF tables (*.fcf)")
	{
		FCFTab *Table = CreateFCFTable();
		if (Table == 0) return;
		if (Table->readData(Filename)) Table->show();
	}
	else if (Filter == "Spectrum (*.spect)")
	{
		Spektrum *Spectrum = CreateSpectrum();
		if (Spectrum == 0) return;
		if (Spectrum->readData(Filename)) Spectrum->show();
	}
	else if (Filter == "SpectList (*.slist)")
	{
		SpectList *SL = CreateSpectList();
		if (SL != 0) if (SL->readData(Filename)) SL->show();
	}
	else if (Filter == "Text table (*.dat *.txt)")
	{
		TableWindow *Tab = createTableWindow();
		if (Tab->readData(Filename)) Tab->show();
	}
	else if (Filter == "Simulated spectrum (*.sspk)")
	{
		SpectSimulation *SM = CreateSpectSimulation();
		if (SM != 0 ? SM->readData(Filename) : false) SM->show();
	}
	else if (Filter == "Log file (*.log)") ResumeCalculation(Filename);
}

void MainWindow::importTermTable()
{
	QString Filter;
	QString Filename = QFileDialog::getOpenFileName(this, "Choose a term energy table to import",
													TermDir, "Term energy tables (*.dat *.term)");
	int n, m;
	if (Filename.isEmpty()) return;
	for (n=m=0; m!=-1; n=m+1) m = Filename.indexOf(DIRSEP, n);
	TermDir = Filename.left(n);
	if (!checkNotOpened(Filename, Filter)) return;
	TermTable *nTerm = CreateTermTable();
	if (nTerm == 0) return;
	if (nTerm->readData(Filename)) nTerm->show();
}

void MainWindow::importDunTable()
{
	QString Filter;
	QString Filename = QFileDialog::getOpenFileName(this, 
											"Choose a table of Dunham coefficients to import",
											DunDir, "Dunham coefficients (*.dat *.txt *.dun)");
	int n, m;
	if (Filename.isEmpty()) return;
	for (n=m=0; m!=-1; n=m+1) m = Filename.indexOf(DIRSEP, n);
	DunDir = Filename.left(n);
	if (!checkNotOpened(Filename, Filter)) return;
	DunTable *Dunham = CreateDunTable();
	if (Dunham == 0) return;
	if (Dunham->readData(Filename)) Dunham->show();
}

void MainWindow::importPotential()
{
	QString Filter;
	QString Filename = QFileDialog::getOpenFileName(this, 
									"Choose a potential to import", PotDir, "Potentials (*)");
	int n, m;
	if (Filename.isEmpty()) return;
	for (n=m=0; m!=-1; n=m+1) m = Filename.indexOf(DIRSEP, n);
	PotDir = Filename.left(n);
	if (!checkNotOpened(Filename, Filter)) return;
	Potential *Pot = CreatePotential();
	if (Pot == 0) return;
	if (Pot->readData(Filename)) Pot->show();
}

void MainWindow::importFitData()
{
	int n, N, I1, I2;
	QString Filename = QFileDialog::getOpenFileName(this, "Choose a fit dataset to import",
											FitDataDir, "Fit datasets (*)");
	if (Filename.isEmpty()) return;
	FitDataDir = Filename.left(n = Filename.lastIndexOf(DIRSEP) + 1);
	N = Filename.indexOf('.', n);
	if (N == -1) N = Filename.length();
	QString Name = "FitData_" + Filename.mid(n, N-n); 
	FitData *FD = CreateFitData();
	if (FD == 0) return;
	QFile F(Filename);
	F.open(QIODevice::ReadOnly);
	QTextStream S(&F);
	QString Buffer = S.readLine();
    QString Buffer2 = S.readLine();
    if (Buffer.left(7) == "Source:" || Buffer.left(5).toInt() > 0 || (Buffer2.left(5).toInt() > 0 && Buffer2.mid(15, 5).toInt() > 0 && Buffer2.mid(20, 5).toInt() > 0))
	{
		F.close();
		if (FD->readData(Filename)) FD->show();
		return;
	}
	while (Buffer != " ISOTOP   V*   J*   V   J  EMESS        ETOL    OBS-CAL ratio"
			&& !S.atEnd())
		Buffer = S.readLine();
	if (S.atEnd())
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "Reading the current file format is not implemented!");
		return;
	}
	IsoTab *Iso = selectIso();
	if (Iso == 0) return;
	int Error = 0;
	TableLine *L = new TableLine[10000];
	Buffer = S.readLine();
	QStringList SL = Buffer.split(' ', QString::SkipEmptyParts);
	double M = pow(10, SL[6].indexOf('.') - SL[6].length() + 1);
	for (N=0; !S.atEnd() && Buffer.left(11) != "  There are" && N < 10000; N++)
	{
		SL = Buffer.split(' ', QString::SkipEmptyParts);
		if (SL.count() < 9)
		{
			Error = 2;
			N--;
			continue;
		}
		I1 = SL[0].toInt();
		I2 = SL[1].toInt();
		for (n=0; (n < Iso->numIso ? I1 != Iso->mNumIso1[n] || I2 != Iso->mNumIso2[n]
								   : false); n++) ;
		if (n == Iso->numIso && Error == 0) Error = 1;
		L[N].Iso = n;
		L[N].vs = SL[2].toInt();
		L[N].vss = SL[4].toInt();
		L[N].Jss = SL[5].toInt();
		L[N].Js = L[N].Jss + SL[3].toInt();
		L[N].WN = -1.0 * SL[6].toDouble();
		L[N].err = double(SL[7].toInt()) * M;
		if (L[N].err < M) L[N].err = 9.0;
		L[N].dev = double(SL[8].toInt()) * M;
		L[N].DevR = L[N].dev / L[N].err;
		L[N].LTab = 0;
		L[N].PN = -1;
		L[N].isTE = true;
		Buffer = S.readLine();
	}
	if (S.atEnd()) QMessageBox::information(this, "MolSpektAnalysis", 
								 "The file is incomplete!");
	if (N == 10000) QMessageBox::information(this, "MolSpektAnalysis", 
								 "The file is too long!");
	if (Error == 2) QMessageBox::information(this, "MolSpektAnalysis", 
				"The file contains invalid lines!");
	else if (Error == 1) QMessageBox::information(this, "MolSpektAnalysis", 
				"At least one isotopologue is not fitting to the selected molecule!");
	FD->setData(L, N);
	FD->setName(Name);
	FD->setSource(Filename);
	FD->show();
}

void MainWindow::importLineTable()
{
	QString Filter;
	QString Filename = QFileDialog::getOpenFileName(this, 
						"Choose a line table to import", LineDir, "Line tables (*.dat *.txt *.lines)");
	int n, m;
	if (Filename.isEmpty()) return;
	for (n=m=0; m!=-1; n=m+1) m = Filename.indexOf(DIRSEP, n);
	LineDir = Filename.left(n);
	if (!checkNotOpened(Filename, Filter)) return;
	LineTable *Table = CreateLineTable();
	if (Table == 0) return;
	if (Table->readData(Filename)) Table->show();
}

void MainWindow::importSpectrum()
{
	QString Filter;
	QString Filename = QFileDialog::getOpenFileName(this, 
						"Choose a spectrum to import", SpectDir, "Spectra (*.dat *.spect)");
	int n, m;
	if (Filename.isEmpty()) return;
	for (n=m=0; m!=-1; n=m+1) m = Filename.indexOf(DIRSEP, n);
	SpectDir = Filename.left(n);
	if (!checkNotOpened(Filename, Filter)) return;
	Spektrum *Spectrum = CreateSpectrum();
	if (Spectrum == 0) return;
	if (Spectrum->readData(Filename)) Spectrum->show();
}

void MainWindow::importSpectList()
{
	QString Filter = "SpectList (*.slist *.dat)";
	QString Filename = QFileDialog::getOpenFileName(this,
						"Choose a spectList to import", Dir, Filter);
	if (Filename.isEmpty()) return;
	Dir = Filename.left(Filename.lastIndexOf(QRegExp("[\\/]")));
	if (!checkNotOpened(Filename, "")) return;
	SpectList *SL = CreateSpectList();
	if (SL != 0) if (SL->readData(Filename)) SL->show();
}

void MainWindow::reload()
{
	MDIChild *AW = activeMDIChild();
	if (AW == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "There is no active window which could be reloaded!");
		return;
	}
	QString T = AW->getTypeString(), N = AW->getName();
	QString WN = (T.isEmpty() ? (N.isEmpty() ? "active window" : "window " + N) 
							  : (N.isEmpty() ? "active " + T : T + " " + N));
	QFile Datei(AW->getFileName());
	if (!Datei.exists())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "The " + WN 
	+ " cannot be reloaded since it never has been saved or the file got deleted in the mean time!");
		return;
	}
	if (!AW->isSaved()) 
		if (QMessageBox::question(this, "MolSpektAnalysis", 
		   					  "Are you sure that you want to discard all changes done on the " + WN 
									  + " since its last saving?", 
							  QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::No)
			return;
	if (!AW->readData())
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "The reloading of the " + WN + " failed!");
}	

void MainWindow::save()
{
	MDIChild *C = activeMDIChild();
	if (C == 0)	QMessageBox::information(this, "MolSpektAnalysis", "The active window cannot be saved!");
	else C->writeData();
}

void MainWindow::saveAs()
{
	QString fileName;
	MDIChild *C = activeMDIChild();
	if (C == 0)	QMessageBox::information(this, "MolSpektAnalysis", "The active window cannot be saved!");
	else fileName = QFileDialog::getSaveFileName(this, "Save " + C->getType(), C->getFileName(), C->getFilter()); 
	if (!fileName.isEmpty()) C->writeData(fileName);
}

void MainWindow::print()
{
	MDIChild *W = activeMDIChild();
	if (W == 0) return;
	QPrinter Printer;
	Printer.setOutputFileName(PrintFile);
	W->Print(Printer);
	PrintFile = Printer.outputFileName();
}

void MainWindow::disableMenues()
{
    saveAct->setEnabled(false);
    saveAsAct->setEnabled(false);
    tableMenu->setEnabled(false);
    SpectrumMenu->setEnabled(false);
    MoleculeMenu->setEnabled(false);
    DunhamMenu->setEnabled(false);
    PlotMenu->setEnabled(false);
    lineTableMenu->setEnabled(false);
    writePotDataAct->setEnabled(false);
    exportMenu->setEnabled(false);
    exportTFDunAct->setEnabled(false);
    exportPotPointsAct->setEnabled(false);
    write2AtInputAct->setEnabled(false);
    writeExpPotFitInputAct->setEnabled(false);
    newFCFTableAct->setEnabled(false);
    exportAsymptoticLevelsAct->setEnabled(false);
    exportWaveFunctionAct->setEnabled(false);
    exportObservedLevelListsAct->setEnabled(false);
}

void MainWindow::WindowActivated(QMdiSubWindow *SubW)
{
	//printf("MainWindow::WindowActivated\n");
	int i;
	bool E = false, D = false, P = false, S = false, L = false, F = false, TW = false;
    bool M = false, exportLineProfile = false;
	QWidget *W = (SubW != 0 ? SubW->widget() : 0);
    if (W != 0 && !W->isVisible()) W = 0;
	PotentialPlot *Plot;
	DiagWindow *DW;
	Potential *Pot;
	ElState *State;
    FitData* fDat;
    CouplingFuncs* CF;
    disableMenues();
	if (dynamic_cast<MDIChild*>(W) != 0)
	{
		saveAct->setEnabled(true);
		saveAsAct->setEnabled(true);
	}
	for (i=0; i < numSpectra; i++) if (W == spectra[i])
	{
		S = true;
		SpectrumMenu->setEnabled(true);
		SpectrumShowAssignmentsOnTopAct->setChecked(spectra[i]->GetShowAssignmentsOnTop());
		SpectrumSetTypeAct->setEnabled(true);
		SpectrumChangeSettingsAct->setEnabled(true);
		SpectrumSetLaserFrequencyAct->setEnabled(true);
		SpectrumFindPeaksAct->setEnabled(true);
		SpectrumShowIntensityHistogramAct->setEnabled(true);
		SpectrumDisplayMarkedAct->setEnabled(true);
		SpectrumShowAssignmentsOnTopAct->setEnabled(true);
		SpectrumAcceptAssignmentsAct->setEnabled(true);
		SpectrumAssignBandsAct->setEnabled(true);
		SpectrumPrev_ProgressionAct->setEnabled(true);
		SpectrumNext_ProgressionAct->setEnabled(true);
		SpectrumTestProgressionAct->setEnabled(true);
		SpectrumFind_SatellitesAct->setEnabled(true);
		SpectrumFindProgressionsAct->setEnabled(true);
		SpectrumContinueProgressionsAct->setEnabled(true);
		SpectrumSingleSLPAct->setEnabled(true);
		SpectrumMultiSLPAct->setEnabled(true);
		SpectrumNormalizeAct->setEnabled(true);
		SpectrumCutAct->setEnabled(true);
		SpectrumCutStrongAct->setEnabled(true);
		SpectrumCutAssignedAct->setEnabled(true);
		SpectrumAddAct->setEnabled(true);
		if ((spectra[i]->getType() == NormalizedAbsorption || spectra[i]->getType() == ThermalEmission) && FELDialog::DataA(this))
			SpectrumFindEmissionLinesAct->setEnabled(true);
		else SpectrumFindEmissionLinesAct->setEnabled(false);
        if (spectra[i]->GetNumFittedLines() > 0) exportLineProfile = true;
	}
	if (dynamic_cast<SpectList*>(W) != 0)
	{
		SpectrumMenu->setEnabled(true);
		SpectrumSetTypeAct->setEnabled(false);
		SpectrumChangeSettingsAct->setEnabled(false);
		SpectrumSetLaserFrequencyAct->setEnabled(false);
		SpectrumFindPeaksAct->setEnabled(false);
		SpectrumShowIntensityHistogramAct->setEnabled(false);
		SpectrumDisplayMarkedAct->setEnabled(false);
		SpectrumShowAssignmentsOnTopAct->setEnabled(false);
		SpectrumAcceptAssignmentsAct->setEnabled(false);
		SpectrumAssignBandsAct->setEnabled(false);
		SpectrumPrev_ProgressionAct->setEnabled(false);
		SpectrumNext_ProgressionAct->setEnabled(false);
		SpectrumTestProgressionAct->setEnabled(false);
		SpectrumFind_SatellitesAct->setEnabled(false);
		SpectrumFindProgressionsAct->setEnabled(false);
		SpectrumContinueProgressionsAct->setEnabled(false);
		SpectrumSingleSLPAct->setEnabled(false);
		SpectrumMultiSLPAct->setEnabled(false);
		SpectrumNormalizeAct->setEnabled(false);
		SpectrumCutAct->setEnabled(false);
		SpectrumCutStrongAct->setEnabled(false);
		SpectrumCutAssignedAct->setEnabled(false);
		SpectrumAddAct->setEnabled(false);
		SpectrumFindEmissionLinesAct->setEnabled(false);
	}
	for (i = 0; i < numMolecules; i++) if (W == molecules[i]) 
	{
		molecules[i]->updateAtoms();
		MoleculeMenu->setEnabled(true);
		newFCFTableAct->setEnabled(true);
		calcFCFTableAct->setText("Calculate new &FCF table...");
		M = true;
		improvePotSeriesAct->setEnabled(true);
	}
	if (dynamic_cast<FCFTab*>(W) != 0)
	{
		MoleculeMenu->setEnabled(true);
		newFCFTableAct->setEnabled(true);
		calcFCFTableAct->setText("Update &FCF table...");
	}
    if ((fDat = dynamic_cast<FitData*>(W)) != 0)
	{
		writeExpPotFitInputAct->setEnabled(true);
		F = true;
	}
	for (i=0; i < numDunTables; i++) if (W == dunTables[i]) 
	{
		D = true;
		DunhamMenu->setEnabled(true);
		exportTFDunAct->setEnabled(true);
		if (dunTables[i]->isSpinRAv() 
				&& ((State = dunTables[i]->getElState()) != 0 ? 
					State->getNumFitDataSets() > 0 : false)) 
			removeFineStructureAct->setEnabled(true);
		else removeFineStructureAct->setEnabled(false);
		if (dunTables[i]->isLambdaDoublingAv()
				&& ((State = dunTables[i]->getElState()) != 0 ?
				    State->getNumFitDataSets() >  0 : false))
			removeLambdaDoublingAct->setEnabled(true);
		else removeLambdaDoublingAct->setEnabled(false);
	}
    if ((CF = dynamic_cast<CouplingFuncs*>(W)) != 0)
    {
        enablePotMenu(false, NoPotential);
        showTexTableAct->setEnabled(true);
    }
    else if ((Pot = dynamic_cast<Potential*>(W)) != 0)
	{
		enablePotMenu(true, Pot->getPotType());
		writePotDataAct->setEnabled(true);
		exportPotPointsAct->setEnabled(true);
		exportAsymptoticLevelsAct->setEnabled(true);
		if (Pot->getCoupledSineWaveFunc() != 0) exportWaveFunctionAct->setEnabled(true);
		P = true;
	}
	else enablePotMenu(false, NoPotential);
	if ((DW = dynamic_cast<DiagWindow*>(W)) != 0)
	{
		PlotMenu->setEnabled(true);
		plotShowMouseCrossAct->setChecked(DW->getShowMouseCross());
		exportPictureAct->setEnabled(true);
		if ((Plot = dynamic_cast<PotentialPlot*>(W)) != 0)
		{
			plotAddPotentialAct->setEnabled(true);
			plotClearHistoryAct->setEnabled(true);
			plotPotSnapShotAct->setEnabled(true);
			plotShowDiagFuncsAct->setEnabled(true);
			plotShowDiagFuncsAct->setChecked(Plot->getShowDiagFuncs());
			plotShowHistoryAct->setEnabled(true);
			plotShowHistoryAct->setChecked(Plot->getShowHistory());
			plotShowPointsAct->setEnabled(true);
			plotShowPointsAct->setChecked(Plot->getShowPoints());
		}
		else
		{
			plotAddPotentialAct->setEnabled(false);
			plotClearHistoryAct->setEnabled(false);
			plotPotSnapShotAct->setEnabled(false);
			plotShowDiagFuncsAct->setEnabled(false);
			plotShowHistoryAct->setEnabled(false);
			plotShowPointsAct->setEnabled(false);
		}
		if (S)
		{
			plotShowMarkerLabelsAct->setEnabled(true);
			plotShowMarkerLabelsAct->setChecked(DW->getShowMarkerLabels());
		}
		else plotShowMarkerLabelsAct->setEnabled(false);
	}
	else exportPictureAct->setEnabled(false);
	for (i=0; i < numLineTables; i++) if (W == lineTables[i]) 
	{
		lineTableMenu->setEnabled(true);
		writeExpPotFitInputAct->setEnabled(true);
		write2AtInputAct->setEnabled(true);
		L = true;
	}
	for (i=0; i < numTermTables; i++) if (W == termTables[i]) E = true;
	if (dynamic_cast<TableWindow*>(W) != 0 && !E)
	{
		TW = true;
		exportTableDataAct->setEnabled(true);
	}
	else exportTableDataAct->setEnabled(false);
	if (D || P || L || F || E) 
	{
		tableMenu->setEnabled(true);
		if (L) updateAct->setEnabled(false);
		else updateAct->setEnabled(true);
		if (L || F) 
		{
			findWrongDataAct->setEnabled(true);
			RemoveDoubledAct->setEnabled(true);
			AssignvsAct->setEnabled(true);
			sortMenu->clear();
			setMenu->setEnabled(true);
			writeTFGSAct->setEnabled(true);
			if (L)
			{
				sortMenu->addAction(sortUpTermIvJAct);
				sortMenu->addAction(SortProgAct);
				sortMenu->addAction(SortIJvPAct);
				sortMenu->addAction(SortIvPJAct);
				sortMenu->addAction(sortbyvsAct);
				sortMenu->addAction(SortFPIntAct);
				sortMenu->addAction(SortSpectrumAct);
				sortMenu->addAction(SortfRemDoubledAct);
                sortMenu->addAction(SortByProgNrAct);
				SetvssAscendingAct->setEnabled(true);
				setvsAct->setEnabled(true);
				SetPNAct->setEnabled(true);
			}
			else 
			{
				sortMenu->addAction(sortTabByDevRatioAct);
				sortMenu->addAction(sortTabByDeviationAct);
				sortMenu->addAction(sortIvJfFAct);
				sortMenu->addAction(sortvsIvJAct);
				sortMenu->addAction(SortProgAct);
                if (fDat->containsDataForMoreThanOneState()) sortMenu->addAction(sortTabByElStateAct);
                sortMenu->addAction(sortTabByElStateAndProgNrAct);
				SetvssAscendingAct->setEnabled(false);
				setvsAct->setEnabled(false);
				SetPNAct->setEnabled(false);
			}
			sortMenu->setEnabled(true);
		}
		else 
		{
			findWrongDataAct->setEnabled(false);
			AssignvsAct->setEnabled(false);
			RemoveDoubledAct->setEnabled(false);
			sortMenu->setEnabled(false);
			setMenu->setEnabled(false);
			writeTFGSAct->setEnabled(false);
		}
		if (F) 
		{
			addCalculatedLevelsAct->setEnabled(true);
			removeDataFSourceAct->setEnabled(true);
			showNumLevelsAct->setEnabled(true);
			extractDataWithComponentAct->setEnabled(true);
            extractChangedDataAct->setEnabled(true);
            extractNewDataAct->setEnabled(true);
			showUncertaintyStatsAct->setEnabled(true);
			changeSourceOffsetAct->setEnabled(true);
			selectDataFSourceAct->setEnabled(true);
			removeSingleLinesAct->setEnabled(true);
			fixJOffsetsAct->setEnabled(true);
		}
		else 
		{
			fixJOffsetsAct->setEnabled(false);
			removeSingleLinesAct->setEnabled(false);
			addCalculatedLevelsAct->setEnabled(false);
			removeDataFSourceAct->setEnabled(false);
			showNumLevelsAct->setEnabled(false);
			extractDataWithComponentAct->setEnabled(false);
            extractChangedDataAct->setEnabled(false);
            extractNewDataAct->setEnabled(false);
			showUncertaintyStatsAct->setEnabled(false);
			changeSourceOffsetAct->setEnabled(false);
			selectDataFSourceAct->setEnabled(false);
		}
	}
	if (!E && showTermViewAct->isEnabled())
	{
		showTermViewAct->setEnabled(false);
		enableShowMenu(false);
	}
	else if (E && !showTermViewAct->isEnabled())
	{
		showTermViewAct->setEnabled(true);
		ShowMenu->setEnabled(true);
	}
	if ((L || E || F) && !showTermPlotAct->isEnabled()) showTermPlotAct->setEnabled(true);
	else if (!L && !E && !F && showTermPlotAct->isEnabled()) showTermPlotAct->setEnabled(false);
    if (!D && !P && CF == 0 && showTexTableAct->isEnabled())
	{
		showTexTableAct->setEnabled(false);
		enableShowMenu(false);
	}
	else if ((D || P) && !showTexTableAct->isEnabled())
	{
		showTexTableAct->setEnabled(true);
		ShowMenu->setEnabled(true);
	}
	if (L || M || F) exportObservedLevelListsAct->setEnabled(true);
    if ((L || P || D || F || DW != 0 || TW || M || exportLineProfile) && !exportMenu->isEnabled())
		exportMenu->setEnabled(true);
    else if (!L && !P && !D && !F && DW == 0 && !M && !TW && !exportLineProfile && exportMenu->isEnabled())
		exportMenu->setEnabled(false);
	if (D || P || E) compareLevelEnergiesAct->setEnabled(true);
	else compareLevelEnergiesAct->setEnabled(false);
    exportLineProfileAct->setEnabled(exportLineProfile);
	//printf("Ende WindowActivated\n");
}

void MainWindow::MdiChildClosed(QWidget* i_closingWindow)
{
    QMdiSubWindow* subWindow = workspace->activeSubWindow();
    if (0 != subWindow)
    {
        QWidget *W = subWindow->widget();
        if (W == i_closingWindow) disableMenues();
    }
}

void MainWindow::MdiChildShown(QWidget* i_shownWindow)
{
    QMdiSubWindow* subWindow = workspace->activeSubWindow();
    if (0 != subWindow)
    {
        QWidget *W = subWindow->widget();
        if (W == i_shownWindow) WindowActivated(subWindow);
    }
}

void MainWindow::deleteRows()
{
	TableWindow *Table = dynamic_cast<TableWindow*> (workspace->activeSubWindow()->widget());
	if (Table == 0) 
	{
		printf("MainWindow::deleteRows() error: the active window is no TableWindow!");
		return;
	}
	Table->DeleteRows();
}

void MainWindow::addRow()
{
	TableWindow *Table = dynamic_cast<TableWindow*> (workspace->activeSubWindow()->widget());
	if (Table == 0)
	{
		printf("MainWindow::addRow() error: the active window is no TableWindow!");
		return;
	}
	Table->AddRow();
}

void MainWindow::cutRows()
{
	TableWindow *Table = dynamic_cast<TableWindow*> (workspace->activeSubWindow()->widget());
	if (Table == 0)
	{
		printf("MainWindow::cutRows() error: the active window is no TableWindow!");
		return;
	}
	QClipboard *CB = QApplication::clipboard();
	int r, c;
	QStringList Row, Text;
	Table->cutRows(nCR, nCC, cuttedRows);
	for (r=0; r < nCR; r++)
	{
		Row.clear();
		for (c=0; c < nCC; c++) Row << cuttedRows[r][c];
		Text << Row.join("\t");
	}
	CB->setText(Text.join("\n"));
}

void MainWindow::copyRows()
{
	TableWindow *Table = dynamic_cast<TableWindow*> (workspace->activeSubWindow()->widget());
	if (Table == 0) 
	{
		printf("MainWindow::copyRows() error: the active Window is no TableWindow!");
		return;
	}
	QClipboard *CB = QApplication::clipboard();
	int r, c;
	QStringList Row, Text;
	Table->copyRows(nCR, nCC, cuttedRows);
	for (r=0; r < nCR; r++)
	{
		Row.clear();
		for (c=0; c < nCC; c++) Row << cuttedRows[r][c];
		Text << Row.join("\t");
	}
	CB->setText(Text.join("\n"));
}

void MainWindow::insertRows()
{
	TableWindow *Table = dynamic_cast<TableWindow*> (workspace->activeSubWindow()->widget());
	if (Table == 0)
	{
		printf("MainWindow::insertRows() error: the active window is no TableWindow!");
		return;
	}
	QClipboard *CB = QApplication::clipboard();
	QString text = CB->text();
	if (!text.isEmpty())
	{
		if (cuttedRows != 0) Destroy(cuttedRows, nCR);
		QStringList Text = text.replace(',', '.').split('\n'), Rows[Text.count()];
		int r, c;
		for (r = nCC = 0, nCR = Text.count(); r < nCR; r++)
		{
			Rows[r] = Text[r].split('\t');
			if (Rows[r].count() > nCC) nCC = Rows[r].count();
		}
		cuttedRows = CreateQString(nCR, nCC);
		for (r=0; r < nCR; r++) for (c=0; c < Rows[r].count(); c++) cuttedRows[r][c] = Rows[r][c];
	}
	Table->insertRows(nCR, nCC, cuttedRows);
}

void MainWindow::search()
{
	int c, tsle;
	TableWindow *Tab = dynamic_cast<TableWindow*> (workspace->activeSubWindow()->widget());
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QStringList ColH = Tab->getHorizontalHeaderLabels();
	QLineEdit *Value = new QLineEdit("", D);
	QComboBox *Columns = new QComboBox(D), *Tsle = new QComboBox(D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	D->setWindowTitle("Search");
	L->addWidget(new QLabel("Search for", D), 0, 0);
	L->addWidget(Tsle, 0, 1);
	L->addWidget(Value, 0, 2);
	L->addWidget(new QLabel("in column", D), 1, 0);
	L->addWidget(Columns, 1, 1, 1, 2);
	L->setRowMinimumHeight(2, 20);
	L->addWidget(OK, 3, 0);
	L->addWidget(Cancel, 3, 2);
	Columns->addItems(ColH);
	Columns->addItem("All columns");
	Columns->setEditable(false);
	Tsle->addItems(QStringList() << "cell text" << "text in cell" << "int =" << "int <" << "int >" << "fabs(int) <"
								<< "fabs(int) >" << "double =" << "double <" << "double >" << "fabs(double) <" << "fabs(double) >");
	Tsle->setEditable(false);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		if ((c = Columns->currentIndex()) == ColH.count()) c=-1;
		tsle = Tsle->currentIndex();
		switch(tsle)
		{
			case 0:
				Tab->search(Value->text(), c, true);
				break;
			case 1:
				Tab->search(Value->text(), c, false);
				break;
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
				Tab->search(c, Value->text().toInt(), tsle - 2);
				break;
			default:
				Tab->search(c, Value->text().toDouble(), tsle - 7);
				break;
		}
	}
	delete D;
}

void MainWindow::shiftValues()
{
	TableWindow *T = dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget());
	if (T == 0) return;
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("Value to shift by:", D), 0, 0);
	QLineEdit *E = new QLineEdit(D);
	L->addWidget(E, 0, 1);
	L->setRowMinimumHeight(1, 20);
	QPushButton *O = new QPushButton("OK", D), *C = new QPushButton("Cancel", D);
	L->addWidget(O, 2, 0);
	L->addWidget(C, 2, 1);
	connect(O, SIGNAL(clicked()), D, SLOT(accept()));
	connect(C, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) T->shiftCellValue(E->text().toInt());
	delete D;
}

void MainWindow::setTextInCells()
{
    QDialog D(this);
    QGridLayout *L = new QGridLayout(&D);
    L->addWidget(new QLabel("Text to set:", &D), 0, 0);
    QLineEdit *Text = new QLineEdit(&D);
    L->addWidget(Text, 0, 1);
    L->setRowMinimumHeight(2, 20);
    QPushButton *OK = new QPushButton("OK", &D), *Cancel = new QPushButton("Cancel", &D);
    L->addWidget(OK, 2, 0);
    L->addWidget(Cancel, 2, 1);
    connect(OK, SIGNAL(clicked()), &D, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), &D, SLOT(reject()));
    if (D.exec() == QDialog::Accepted) dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget())->setCellText(Text->text());
}

void MainWindow::findWrongData()
{
	FitData *F = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	if (F!=0) F->findWrongData();
	LineTable *L = dynamic_cast<LineTable*>(workspace->activeSubWindow()->widget());
	if (L!=0) L->findErrors();
}

void MainWindow::extractDataWithComponent()
{
	FitData *F = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget()), *nFD;
	if (F==0) return;
	ElState *S = F->getElState();
	int l = (S != 0 ? S->getLambda() : 1), nv, *av, v1, v2, J1, J2, v, vmin, *vmax, Jmin, Jmax, vM;
	QList<int> FList = F->getaFC();
	int NF = FList.count(), n, m, NLines, f, J;
	double *Offsets, vI = 0.0, dv = 0.0;
	QStringList Names;
	QString B;
	TableLine *Lines;
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QComboBox *ef = new QComboBox(D), *FC = new QComboBox(D), *v1B = new QComboBox(D), *v2B = new QComboBox(D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	QRadioButton *MaxMin = new QRadioButton("Max/min", D), *BelowLine = new QRadioButton("Below line (J)", D);
	QLineEdit *J1E = new QLineEdit("0", D), *J2E = new QLineEdit(QString::number(F->getJMax()), D);
	L->addWidget(new QLabel("Component to extract:", D), 0, 0, 1, 4);
	L->addWidget(new QLabel("Parity:", D), 1, 0);
	L->addWidget(ef, 1, 1);
	L->addWidget(new QLabel("F:", D), 1, 2);
	L->addWidget(FC, 1, 3);
	L->setRowMinimumHeight(2, 20);
	L->addWidget(new QLabel("Level range:", D), 3, 0, 1, 3);
	L->addWidget(MaxMin, 4, 0, 1, 2);
	L->addWidget(BelowLine, 4, 2, 1, 2);
	L->addWidget(new QLabel("v1:", D), 5, 0);
	L->addWidget(v1B, 5, 1);
	L->addWidget(new QLabel("J1:", D), 5, 2);
	L->addWidget(J1E, 5, 3);
	L->addWidget(new QLabel("v2:", D), 6, 0);
	L->addWidget(v2B, 6, 1);
	L->addWidget(new QLabel("J2:", D), 6, 2);
	L->addWidget(J2E, 6, 3);
	L->setRowMinimumHeight(7, 20);
	L->addWidget(OK, 8, 0, 1, 2);
	L->addWidget(Cancel, 8, 2, 1, 2);
	ef->addItem("e");
	if (l >= 1)
	{
		ef->addItem("f");
		ef->addItem("both");
	}
	for (n=0; n < NF; n++) FC->addItem(QString::number(FList[n]));
	if (NF > 1) FC->addItem("all");
	MaxMin->setChecked(true);
	F->getav(nv, av);
	for (n=0; n < nv; n++) 
	{
		v1B->addItem(B = QString::number(av[n]));
		v2B->addItem(B);
	}
	v1B->setEditable(false);
	v2B->setEditable(false);
	v2B->setCurrentIndex(nv - 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		v1 = av[v1B->currentIndex()];
		v2 = av[v2B->currentIndex()];
		J1 = J1E->text().toInt();
		J2 = J2E->text().toInt();
		if (MaxMin->isChecked())
		{
			if (v1 > v2)
			{
				v = v1;
				v1 = v2;
				v2 = v;
			}
			if (J1 > J2)
			{
				J = J1;
				J1 = J2;
				J2 = J;
			}
			vM = (v2 >= 0 ? v2 : F->getvMax());
			Jmax = (J2 >= 0 ? J2 : F->getJMax()); 
			vmin = (v1 >= 0 ? v1 : 0);
			Jmin = (J1 >= 0 ? J1 : 0);
		}
		else 
		{
			if (J1 > J2)
			{
				v = v1;
				v1 = v2;
				v2 = v;
				J = J1;
				J1 = J2;
				J2 = J;
			}
			vI = double(v2 - v1) / (J2 - J1);
			if (vI < 0)
			{
				if ((vM = v1 - int(double(J1) * vI)) > (v = F->getvMax())) vM = v;
				dv = double(v1) - double(J1) * vI;
				if ((Jmax = J2 - int(double(v2) / vI)) > (J = F->getJMax())) Jmax = J;
				Jmin = 0;
				vmin = 0;
			}
			else
			{
				vM = F->getvMax();
				Jmax = F->getJMax();
				vmin = 0;
				if ((Jmin = J1 + int(double(v1) / vI)) < 0) Jmin = 0;
			}
		}
		vmax = new int[Jmax + 1];
		if (MaxMin->isChecked()) for (J = Jmin; J <= Jmax; J++) vmax[J] = vM;
		else for (J = Jmin; J <= Jmax; J++, dv += vI) vmax[J] = (dv <= vM ? int(dv) : vM);
		F->getData(Lines, NLines, (l == 0 ? -1 : 1 - ef->currentIndex()), 
				   (NF == 1 || (f = FC->currentIndex()) == NF ? -2 : FList[f]),
				   (MaxMin->isChecked() && v1 == v2 ? v1 : -2), Jmax);
		if (BelowLine->isChecked() || v1 != v2)
		{
			for (n=0; (n < NLines ? Lines[n].vss >= vmin && Lines[n].vss <= vmax[Lines[n].Jss] : false); n++) ;
			for (m=n+1; m < NLines; m++) if (Lines[m].vss >= vmin && Lines[m].vss <= vmax[Lines[m].Jss]) Lines[n++] = Lines[m];
			NLines = n;
		}
		nFD = CreateFitData();
		if (S != 0) S->addFitData(nFD);
		nFD->setData(Lines, NLines);
		nFD->setName("part" + F->getName());
		nFD->setFileName("Part" + F->getFName());
		nFD->setSource("Part of " + F->getName() + "with source " + F->getSource());
		nFD->setvMax(vM);
		nFD->setJMax(Jmax);
		F->getSourceOffset(Names, Offsets);
		nFD->setSourceOffset(Names, Offsets);
		nFD->show();
	}
	delete D;
	if (nv > 0) delete[] av;
}

FitData* MainWindow::showSelectFitDataDialog(ElState * const i_elStateToSelectFitDataFrom, const QString &i_windowTitle, const QString &i_text)
{
    if (0 == i_elStateToSelectFitDataFrom) return 0;
    QDialog* D = new QDialog(this);
    QGridLayout *L = new QGridLayout(D);
    D->setWindowTitle(i_windowTitle);
    L->addWidget(new QLabel(i_text, D), 0, 0, 1, 2);
    L->addWidget(new QLabel("Selected fit data:", D), 1, 0);
    QComboBox *FDBox = new QComboBox(D);
    for (int n=0; n < i_elStateToSelectFitDataFrom->getNumFitDataSets(); ++n) FDBox->addItem(i_elStateToSelectFitDataFrom->getFitDataName(n));
    FDBox->setEditable(false);
    L->addWidget(FDBox, 1, 1);
    L->setRowMinimumHeight(2, 20);
    QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
    L->addWidget(OK, 3, 0);
    L->addWidget(Cancel, 3, 1);
    connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
    connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
    FitData* rValue = (D->exec() == QDialog::Accepted ? i_elStateToSelectFitDataFrom->getFitData(FDBox->currentIndex()) : 0);
    delete D;
    return rValue;
}

void MainWindow::extractNewData()
{
    FitData* ToCall = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
    FitData* ToCompare = showSelectFitDataDialog(ToCall->getElState(), "Extract new data", "Please select the fit data to compare width.");
    if (0 != ToCompare) ToCall->extractNewData(ToCompare, CreateFitData());
}

void MainWindow::extractChangedData()
{
    FitData* ToCall = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
    FitData* ToCompare = showSelectFitDataDialog(ToCall->getElState(), "Extract changed data", "Please select the fit data to compare width.");
    if (0 != ToCompare) ToCall->extractChangedData(ToCompare, CreateFitData());
}

void MainWindow::update()
{
    int NumWFPoints = NumPoints;
	FitData *F = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	TermTable *T = dynamic_cast<TermTable*>(workspace->activeSubWindow()->widget());
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	Potential *P = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	ElState *S;
	if (F != 0) F->updateData();
	else if (D != 0)
	{
		if ((S = D->getElState()) != 0) F = S->getFitData();
		else return;
		if (F == 0) 
		{
			F = CreateFitData();
			S->addFitData(F);
		}
		F->updateData();
		if (!F->isVisible()) F->show();
	}
	else if (P != 0) P->updateFitData();
	else if (T != 0)
	{
		if ((S = T->getElState()) == 0) return;
		QString Source = T->getSource();
		int n, N;
		for (n=0, N = S->getNumPotentials(); 
			 (n<N ? Source.indexOf(S->getPotentialName(n)) == -1 : false); n++) ;
		if (n<N) 
		{
			P = S->getPotential(n);
            if (P!= 0) P->calcTermEnergies(T, NumWFPoints);
            else if ((P = S->getPotential()) != 0) P->calcTermEnergies(T, NumWFPoints);
            else if ((D = S->getDunTable()) != 0) D->calcTermEnergies(T, NumWFPoints);
		}
		else
		{
			for (n=0, N = S->getNumDunTables(); 
				 (n<N ? Source.indexOf(S->getDunTableName(n)) == -1 : false); n++) ;
            if (n<N ? (D = S->getDunTable(n)) != 0 : false) D->calcTermEnergies(T, NumWFPoints);
            else if ((P = S->getPotential()) != 0) P->calcTermEnergies(T, NumWFPoints);
            else if ((D = S->getDunTable()) != 0) D->calcTermEnergies(T, NumWFPoints);
		}
	}
}

QString MainWindow::getDir(MDIChild::Type type)
{
	switch (type)
	{
		case MDIChild::AtomData:
			if (!AtomDir.isEmpty()) return AtomDir;
			break;
		case MDIChild::MolData:
			if (!MolDir.isEmpty()) return MolDir;
			break;
		case MDIChild::TermEnergyTable:
			if (!TermDir.isEmpty()) return TermDir;
			break;
		case MDIChild::LineTab:
			if (!LineDir.isEmpty()) return LineDir;
			break;
		case MDIChild::FranckCondonTable:
			if (!FCFDir.isEmpty()) return FCFDir;
			break;
		case MDIChild::DunhamTable:
			if (!DunDir.isEmpty()) return DunDir;
			break;
		case MDIChild::PotData:
			if (!PotDir.isEmpty()) return PotDir;
			break;
		case MDIChild::FitDataSet:
			if (!FitDataDir.isEmpty()) return FitDataDir;
			break;
		case MDIChild::Spect:
			if (!SpectDir.isEmpty()) return SpectDir;
			break;
		case MDIChild::Pict:
			if (!PictureDir.isEmpty()) return PictureDir;
			break;
		default:
			return Dir;
			break;
	}
	return Dir;
}

void MainWindow::setDir(QString dir, MDIChild::Type type)
{
	Dir = dir;
	switch (type)
	{
		case MDIChild::AtomData:
			AtomDir = Dir;
			break;
		case MDIChild::MolData:
			MolDir = Dir;
			break;
		case MDIChild::TermEnergyTable:
			TermDir = Dir;
			break;
		case MDIChild::LineTab:
			LineDir = Dir;
			break;
		case MDIChild::FranckCondonTable:
			FCFDir = Dir;
			break;
		case MDIChild::DunhamTable:
			DunDir = Dir;
			break;
		case MDIChild::PotData:
			PotDir = Dir;
			break;
		case MDIChild::FitDataSet:
			FitDataDir = Dir;
			break;
		case MDIChild::Spect:
			SpectDir = Dir;
			break;
		case MDIChild::Pict:
			PictureDir = Dir;
			break;
		default:
			break;
	}
}

void MainWindow::importAbInitioPotSet()
{
	AIPSIDialog *D = new AIPSIDialog(this);
	D->show();
}

void MainWindow::showIsotopologues()
{
	int i;
    Molecule *M = getCurrentMolecule();
	IsoTab *Iso = M->getIso();
	QTableWidget *Tab = new QTableWidget(Iso->numIso, 5);
	QString S1 = M->getAtom1()->getChSymb(), S2 = M->getAtom2()->getChSymb();
	workspace->addSubWindow(Tab);
	Tab->setHorizontalHeaderLabels(QStringList() << "isotop of " + S1 << "isotop of " + S2 
			<< "rel. NA" << "red. mass [u]" << "rel. red. mass");
	for (i=0; i < Iso->numIso; i++)
	{
		Tab->setItem(i, 0, new QTableWidgetItem(QString::number(Iso->mNumIso1[i])));
		Tab->setItem(i, 1, new QTableWidgetItem(QString::number(Iso->mNumIso2[i])));
		Tab->setItem(i, 2, new QTableWidgetItem(QString::number(Iso->relNA[i], 'g', 10)));
		Tab->setItem(i, 3, new QTableWidgetItem(QString::number(Iso->redMass[i], 'g', 10)));
		Tab->setItem(i, 4, new QTableWidgetItem(QString::number(Iso->relRedMass[i], 'g', 10)));
	}
	Tab->setWindowTitle("Isotopomers of " + S1 + S2);
	Tab->show();
	delete Iso;
}

void MainWindow::showLevelNumbers()
{
	//printf("MainWindow::showLevelNumbers()\n");
	int mv, mJ, NS, S, gS, v, J, I, s;
    Molecule *M = getCurrentMolecule();
	IsoTab *Iso;
	bool ****Level;
	//printf("Vor getKnownLevels\n");
	M->getKnownLevels(Iso, NS, mv, mJ, Level);
	//printf("Nach getKnownLevels\n");
	int **Numbers = CreateInt(Iso->numIso, NS);
	QTableWidget *Tab = new QTableWidget(Iso->numIso + 1, NS + 1);
	QStringList hL, vL;
	for (I=0; I < Iso->numIso; I++) 
		vL << QString::number(Iso->mNumIso1[I]) + *Iso->chSymb1 
				+ QString::number(Iso->mNumIso2[I]) + *Iso->chSymb2;
	vL << "total";
	for (s=0; s < NS; s++) hL << M->getState(s);
	hL << "total";
	Tab->setHorizontalHeaderLabels(hL);
	Tab->setVerticalHeaderLabels(vL); 
	for (I=0; I < Iso->numIso; I++) 
	{
		for (S=s=0; s < NS; s++) 
		{
			for (v=0, Numbers[I][s] = 0; v <= mv; v++) for (J=0; J < mJ; J++) if (Level[I][s][v][J])
						Numbers[I][s]++;
			Tab->setItem(I, s, new QTableWidgetItem(QString::number(Numbers[I][s])));
			S += Numbers[I][s];
		}
		Tab->setItem(I, NS, new QTableWidgetItem(QString::number(S)));
	}
	for (gS=s=0; s < NS; s++)
	{
		for (S=I=0; I < Iso->numIso; I++) S += Numbers[I][s];
		Tab->setItem(Iso->numIso, s, new QTableWidgetItem(QString::number(S)));
		gS += S;
	}
	Tab->setItem(Iso->numIso, NS, new QTableWidgetItem(QString::number(gS)));
	Tab->setWindowTitle("Level numbers of " + *Iso->chSymb1 + *Iso->chSymb2);
	workspace->addSubWindow(Tab);
	Tab->show();
	//printf("Vor Destroy\n");
	Destroy(Level, Iso->numIso, NS, mv + 1);
	Destroy(Numbers, Iso->numIso);
	delete Iso;
	//printf("Ende showLevelNumbers\n");
}

void MainWindow::compareLevelEnergies()
{
	int n, v1, v2, J1, J2, Jmax, *vmax, Jmin, vmin, v, I, J, NP, ND, NT, CI = 0, NI, vM, NC, N, c, Mv = 0, MI = 0, Mc = 0, MJ = 0;
	double ****Data1, ****Data2, vI = 0.0, dv = 0.0, Offset, StdDev, Mdev = 0;
	TableWindow *ActTable = dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget());
	TermTable *Tab1, *Tab2;
	ElState *State = ActTable->getElState();
	MDIChild::Type type = ActTable->getType();
	QString Name = ActTable->getName(), Buffer, TypeS = ActTable->getTypeString();
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QComboBox *Table2B = new QComboBox(D);
	QRadioButton *MaxMin = new QRadioButton("Max and min", D), *BelowLine = new QRadioButton("Below line (J)", D);
	QLineEdit *v1E, *v2E, *J1E, *J2E;
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	D->setWindowTitle("Select data to compare with");
	L->addWidget(new QLabel("Data source:"), 0, 0, 1, 2);
	L->addWidget(Table2B, 0, 2, 1, 2);
	L->addWidget(new QLabel("Level range:"), 1, 0, 1, 4);
	L->addWidget(MaxMin, 2, 0, 1, 2);
	L->addWidget(BelowLine, 2, 2, 1, 2);
	L->addWidget(new QLabel("v:", D), 3, 0);
	L->addWidget(v1E = new QLineEdit("-1", D), 3, 1);
	L->addWidget(new QLabel("J:", D), 3, 2);
	L->addWidget(J1E = new QLineEdit("-1", D), 3, 3);
	L->addWidget(new QLabel("v:", D), 4, 0);
	L->addWidget(v2E = new QLineEdit("-1", D), 4, 1);
	L->addWidget(new QLabel("J:", D), 4, 2);
	L->addWidget(J2E = new QLineEdit("-1", D), 4, 3);
	L->setRowMinimumHeight(5, 20);
	L->addWidget(OK, 6, 0, 1, 2);
	L->addWidget(Cancel, 6, 2, 1, 2);
	for (n=0, NT = State->getNumTermTables(); n < NT; n++)
	{
		Buffer = State->getTermTableName(n);
		if (Buffer == Name && type == MDIChild::TermEnergyTable) CI = n;
		else Table2B->addItem(Buffer);
	}
	for (n=0, ND = State->getNumDunTables(); n < ND; n++)
	{
		Buffer = State->getDunTableName(n);
		if (Buffer == Name && type == MDIChild::DunhamTable) CI = n + NT;
		else Table2B->addItem(Buffer);
	}
	for (n=0, NP = State->getNumPotentials(); n < NP; n++)
	{
		Buffer = State->getPotentialName(n);
		if (Buffer == Name && type == MDIChild::PotData) CI = n + NT + ND;
		else Table2B->addItem(Buffer);
	}
	Table2B->setEditable(false);
	MaxMin->setChecked(true);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		v1 = v1E->text().toInt();
		v2 = v2E->text().toInt();
		J1 = J1E->text().toInt();
		J2 = J2E->text().toInt();
		if (MaxMin->isChecked())
		{
			if (v1 > v2)
			{
				v = v1;
				v1 = v2;
				v2 = v;
			}
			if (J1 > J2)
			{
				J = J1;
				J1 = J2;
				J2 = J;
			}
			vM = (v2 >= 0 ? v2 : cMaxv);
			Jmax = (J2 >= 0 ? J2 : cMaxJ); 
			vmin = (v1 >= 0 ? v1 : 0);
			Jmin = (J1 >= 0 ? J1 : 0);
		}
		else 
		{
			if (J1 > J2)
			{
				v = v1;
				v1 = v2;
				v2 = v;
				J = J1;
				J1 = J2;
				J2 = J;
			}
			vI = double(v2 - v1) / (J2 - J1);
			if (vI < 0)
			{
				vM = v1 - int(double(J1) * vI);
				dv = double(v1) - double(J1) * vI;
				Jmax = J2 - int(double(v2) / vI);
				Jmin = 0;
				vmin = 0;
			}
			else
			{
				vM = cMaxv;
				Jmax = cMaxJ;
				vmin = 0;
				Jmin = J1 + int(double(v1) / vI);
			}
		}
		if (type == MDIChild::TermEnergyTable) Tab1 = dynamic_cast<TermTable*>(ActTable);
		else if (type == MDIChild::DunhamTable) dynamic_cast<DunTable*>(ActTable)->calcTermEnergies(Tab1 = 0, false, vM, Jmax);
		else dynamic_cast<Potential*>(ActTable)->calcTermEnergies(Tab1 = 0, vM, Jmax, false);
		if ((v = Tab1->getMaxv()) < vM) vM = v;
		if ((J = Tab1->getMaxJ()) < Jmax) Jmax = J;
		NI = Tab1->getNumIso();
		NC = Tab1->getNumComp();
		Data1 = Tab1->getData();
		if ((n = Table2B->currentIndex()) >= CI) n++;
		if (n < NT) 
		{
			Tab2 = State->getTermTable(n);
			Buffer = "TermTable";
		}
		else if (n < NT + ND) 
		{
			State->getDunTable(n - NT)->calcTermEnergies(Tab2 = 0, false, vM, Jmax);
			Buffer = "DunhamSet";
		}
		else 
		{
			State->getPotential(n - NT - ND)->calcTermEnergies(Tab2 = 0, vM, Jmax, false);
			Buffer = "Potential";
		}
		if ((v = Tab2->getMaxv()) < vM) vM = v;
		if ((J = Tab2->getMaxJ()) < Jmax) Jmax = J;
		if ((n = Tab2->getNumIso()) < NI) NI = n;
		if ((n = Tab2->getNumComp()) != NC)
		{
			NC = 1;
			QMessageBox::warning(this, "MolSpectAnalysis", 
				"Attention: the numbers of fine structure or e/f components of both tables are different, so only the first component of each table will be used!");
		}
		Data2 = Tab2->getData();
		vmax = new int[Jmax + 1];
		if (MaxMin->isChecked()) for (J = Jmin; J <= Jmax; J++) vmax[J] = vM;
		else for (J = Jmin; J <= Jmax; J++, dv += vI) vmax[J] = (dv <= vM ? int(dv) : vM);
		for (c=N=0, Offset = 0.0; c < NC; c++) for (I=0; I < NI; I++) for (J = Jmin; J <= Jmax; J++) 
			for (v = vmin; v <= vmax[J]; v++) if (Data1[c][I][v][J] > 0.0 && Data2[c][I][v][J] > 0.0)
		{
			Offset += (dv = Data1[c][I][v][J] - Data2[c][I][v][J]);
			N++;
			if (fabs(dv) > fabs(Mdev))
			{
				Mc = c;
				MI = I;
				Mv = v;
				MJ = J;
				Mdev = dv;
			}
		}
		Offset /= N;
		for (c=0, StdDev = 0.0; c < NC; c++) for (I=0; I < NI; I++) for (J = Jmin; J <= Jmax; J++) 
			for (v = vmin; v <= vmax[J]; v++) if (Data1[c][I][v][J] > 0.0 && Data2[c][I][v][J] > 0.0)
		{
			dv = Data1[c][I][v][J] - Data2[c][I][v][J] - Offset;
			if (dv*dv > 1.0) printf("Data1[%d][%d][%d][%d]=%f, Data2[%d][%d][%d][%d]=%f\n", c, I, v, J, Data1[c][I][v][J],
										c, I, v, J, Data2[c][I][v][J]);
			StdDev += dv * dv;
		}
		StdDev = sqrt(StdDev / (N-1));
		QMessageBox::information(this, "MolSpectAnalysis", 
			"The offset of " + TypeS + ' ' + Name + " compared to " + Buffer + ' ' + Table2B->currentText() + " is " 
			+ QString::number(Offset, 'g', 4) + " cm^-1 and the residual standard deviation " + QString::number(StdDev, 'g', 3) 
			+ " cm^-1. The maximum deviation is " + QString::number(Mdev, 'g', 4) + " for c=" + QString::number(Mc) + ", I=" 
			+ QString::number(MI) + ", v=" + QString::number(Mv) + " and J=" + QString::number(MJ) + '.');
		delete[] vmax;
	}
	delete D;
}

void MainWindow::calcFCFTable()
{
	bool update;
	if (dynamic_cast<FCFTab*>(workspace->activeSubWindow()->widget()) == 0) update = false;
	else update = true;
    calcFCFDiag(NumFCF_WFPoints, update);
}

void MainWindow::checkAllConnections()
{
    Molecule* Mol = getCurrentMolecule();
    bool Result = Mol->checkAllConnections();
    if (Result) QMessageBox::information(this, "MolSpektAnalysis", "Successfully verified, that all connections are working!");
    else QMessageBox::information(this, "MolSpektAnalysis", "Verification of all connections was not successfull, there are broken connections!");
}

void MainWindow::shrinkAllSpectRefs()
{
    Molecule* Mol = getCurrentMolecule();
    Mol->shrinkAllSpectRefs();
}

Molecule* MainWindow::getCurrentMolecule()
{
    Molecule *Mol = dynamic_cast<Molecule*>(workspace->activeSubWindow()->widget());
    if (Mol == 0)
    {
        TableWindow *T = dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget());
        Mol = T->getMolecule();
    }
    return Mol;
}

void MainWindow::showTermView()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numTermTables ? W != termTables[i] : false); i++) ;
	if (i == numTermTables)
	{
		printf("MainWindow::showTermView() error: The active window doesn't belong to a table of term energies!\n");
		return;
	}
	TermView *T = new TermView(termTables[i]);
	workspace->addSubWindow(T);
	T->show();
}

void MainWindow::showTermPlot()
{
	int i, j;
	TableWindow *W = dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget());
	ElState *S = W->getElState();
	Molecule *Mol = W->getMolecule();
	if (Mol == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", "The current table has to be assigned to a molecule first!");
		return;
	}
	TermPlot *T = new TermPlot(this, Mol, S);
	Transition *Tr;
	workspace->addSubWindow(T);
	connect(this, SIGNAL(SpectrumChanged(Spektrum*)), T, SLOT(UpdateMarked()));
	connect(this, SIGNAL(LineTableChanged()), T, SLOT(UpdateLines()));
	for (i=0; i < Mol->getNumTransitions(); i++) if ((Tr = Mol->getTransitionP(i))->getUpperState() == S)
		for (j=0; j < Tr->getNumLineTables(); j++) 
			connect(Tr->getLineTable(j), SIGNAL(SelChanged()), T, SLOT(Paint()));
	if (S != 0) for (i=0; i < S->getNumFitDataSets(); i++) 
	{
		connect(S->getFitData(i), SIGNAL(SelChanged()), T, SLOT(Paint()));
		connect(S->getFitData(i), SIGNAL(propertiesChanged()), T, SLOT(Paint()));
	}
	T->show();
}

void MainWindow::showTexTable()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numDunTables ? W != dunTables[i] : false); i++) ;
	if (i < numDunTables) 
	{
		QTextEdit *window = new QTextEdit(this);
		window->setPlainText(dunTables[i]->getTexTable());
		workspace->addSubWindow(window);
		window->show();
	}
	else 
	{
        CouplingFuncs* CF;
        for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
        if (i < numPotentials) potentials[i]->getTexTable(NumPoints);
        else if ((CF = dynamic_cast<CouplingFuncs*>(W)) != 0) CF->getTexTable();
        else
		{
			printf("MainWindow::showTexTable() error: For the active window has no function defined for creating a tex table!\n");
			return;
		}
	}
}

void MainWindow::showFCFJDependency()
{
	FCFJDep *D = new FCFJDep(this);
	workspace->addSubWindow(D);
	D->show();
}

void MainWindow::showFCFTable()
{
	FCFTable *FCF = new FCFTable(this);
	//printf("Vor addWindow\n");
	workspace->addSubWindow(FCF);
	//printf("Vor Show\n");
	FCF->show();
	//printf("Ende show FCFTable\n");
}

void MainWindow::showFCFTable(Molecule *Mol, QString LS, QString US)
{
	FCFTable *FCF = new FCFTable(this, Mol, LS, US);
	workspace->addSubWindow(FCF);
	FCF->show();
}

void MainWindow::showWaveFuncPlot()
{
    WaveFuncPlot *WFP = new WaveFuncPlot(NumPoints, this);
	workspace->addSubWindow(WFP);
	WFP->show();
}

void MainWindow::showWaveFuncPlot(Potential *Pot, int Iso, int J, int v)
{
    WaveFuncPlot *WFP = new WaveFuncPlot(NumPoints, this, Pot, Iso, J, v);
	workspace->addSubWindow(WFP);
	WFP->show();
}

void MainWindow::showDataPlot()
{
	DataPlot *DP = new DataPlot(this);
	workspace->addSubWindow(DP);
	connect(this, SIGNAL(MoleculesChanged()), DP, SLOT(moleculesChanged()));
	DP->show();
}

void MainWindow::showResidualPlot()
{
	ResidualPlot *RP = new ResidualPlot(this);
	workspace->addSubWindow(RP);
	connect(this, SIGNAL(MoleculesChanged()), RP, SLOT(moleculesChanged()));
	RP->show();
}

void MainWindow::showSpectSimulation()
{
	SpectSimulation *S = CreateSpectSimulation();
	//printf("Vor show, S=%d\n", S);
	if (S != 0) S->show();
}

void MainWindow::addDunTableLine()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numDunTables ? W != dunTables[i] : false); i++) ;
	if (i == numDunTables)
	{
		printf("MainWindow::addDunTableLine() error: The active window doesn't belong to a table of Dunham coefficients\n");
		return;
	}
	dunTables[i]->addTableLine();
}

void MainWindow::calcTermEnergies()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numDunTables ? W != dunTables[i] : false); i++) ;
	if (i == numDunTables)
	{
		printf("MainWindow::calcTermEnergies() error: The active window doesn't belong to a table of Dunham coefficients!\n");
		return;
	}
	TermTable *T = 0;
	dunTables[i]->calcTermEnergies(T, true);
}

void MainWindow::removeFineStructure()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	if (D == 0) return;
	TableLine *TL;
	int N;
	D->removeFinestructure(TL, N);
	if (N==0) return;
	FitData *CFD = D->getElState()->getFitData(), *NFD = CreateFitData();
	NFD->setData(TL, N);
	NFD->setName(CFD->getName() + "FSRem");
	QString FName = CFD->getFileName();
	NFD->setFileName(FName.left(FName.length() - 5) + "FSRem.fdat");
	NFD->setSource(CFD->getSource() + " with fine structure removed by " + D->getName());
	NFD->setvMax(CFD->getMaxv());
	NFD->setJMax(CFD->getJMax());
	D->getElState()->addFitData(NFD);
	NFD->show();
}

void MainWindow::removeLambdaDoubling()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	if (D==0) return;
	TableLine *TL;
	int N;
	D->removeLambdaDoubling(TL, N);
	if (N==0) return;
	ElState *S = D->getElState();
	FitData *CFD = S->getFitData(), *NFD = CreateFitData();
	NFD->setData(TL, N);
	NFD->setName(CFD->getName() + "LDRem");
	QString FName = CFD->getFileName();
	NFD->setFileName(FName.left(FName.length() - 5) + "LDRem.fdat");
	NFD->setSource(CFD->getSource() + " with lambda doubling removed by " + D->getName());
	NFD->setvMax(CFD->getMaxv());
	NFD->setJMax(CFD->getJMax());
	D->getElState()->addFitData(NFD);
	NFD->show();
}

void MainWindow::applyAdCorr()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	TableLine *TL;
	int N;
	D->applyAdiabaticCorrection(TL, N);
	if (N==0) return;
	ElState *S = D->getElState();
	FitData *CFD = S->getFitData(), *NFD = CreateFitData();
	NFD->setData(TL, N);
	NFD->setName(CFD->getName() + "adCorr");
	QString FName = CFD->getFileName();
	NFD->setFileName(FName.left(FName.length() - 5) + "adCorr.fdat");
	NFD->setSource(CFD->getSource() + " with adiabatic correction applied by " + D->getName());
	NFD->setvMax(CFD->getMaxv());
	NFD->setJMax(CFD->getJMax());
	D->getElState()->addFitData(NFD);
	NFD->show();
}

void MainWindow::calcRe()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	int n, e;
	double Re, Err;
	D->calcRe(Re, Err);
	if (Re == 0.0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
		 "Error: calculating Re is not possible with the current table, maybe it doesn't contain a rotational constant or it is not assigned to a molecule consisting of two sufficiently known atoms.");
		return;
	}
	if (Err == 0.0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", "The calculated value is Re=" + QString::number(Re) + '.');
		return;
	}
	n = int(ceil(-log10(Err)));
	e = int(round(Err * pow(10.0, n)));
	QMessageBox::information(this, "MolSpectAnalysis", 
							 "The calculated value is Re=" + QString::number(Re, 'f', n) + '(' + QString::number(e) + ").");
}

void MainWindow::testKratzer()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	double Y02c, Y02f;
	D->testKratzer(Y02f, Y02c);
	if (Y02f != 0.0 && Y02c != 0.0) 
		QMessageBox::information(this, "MolSpectAnalysis", 
			"The fitted value for Y02 is " + QString::number(Y02f) 
			+ ", while the value calculated from the rotational and vibrational constants is " + QString::number(Y02c) 
			+ ", thus the Kratzer relation holds by " + QString::number(1e2 * fabs((Y02f-Y02c) / Y02f), 'f', 2) + "%.");
	else 
		QMessageBox::information(this, "MolSpectAnalysis", 
			"Error: The Kratzer relation cannot be tested with the current table of Dunham coefficients, because it doesn't contain all the necessary coefficients.");
}

void MainWindow::showBandConstants()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	int v, c, ND, NC, NAD, NLD, NSR, nv;
	double **C, **err;
	QStringList HL;
	D->getBandConstants(C, err, NC, NLD, NSR, NAD, nv);
	TableWindow *Tab = new TableWindow(MDIChild::TextTable1, this);
	QString **Data = CreateQString(nv, 2 * NC);
	for (c=0, ND = NC - NAD - NLD - NSR; c < ND; c++)
	{
		switch (c)
		{
			case 0:
				HL << "G" << "err";
				break;
			case 1:
				HL << "B" << "err";
				break;
			case 2:
				HL << "D" << "err";
				break;
			case 3:
				HL << "H" << "err";
				break;
			default:
				HL << QString::number(c) << "err";
				break;
		}
	}
	for (c=0; c < NLD; c++) HL << "LD" + QString::number(c) << "err";
	for (c=0; c < NSR; c++) HL << "SR" + QString::number(c) << "err";
	for (c=0; c < NAD; c++) HL << "AD" + QString::number(c) << "err";
	Tab->setHorizontalHeader(HL);
	Tab->setWindowTitle("Band constants calculated from " + D->getName());
	for (v=0; v < nv; v++) for (c=0; c < NC; c++)
	{
		Data[v][2*c] = QString::number(C[v][c], 'g', 12);
		Data[v][2*c+1] = QString::number(err[v][c], 'g', 12);
	}
	Tab->setData(Data, nv, 2 * NC);
	Destroy(C, nv);
	Destroy(err, nv);
	Destroy(Data, nv);
	workspace->addSubWindow(Tab);
	Tab->show();
}

void MainWindow::calcY00Te()
{
	DunTable *D = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	double Y00, Y00Err, Te, TeErr;
	int n, m;
	D->calcTeY00(Te, TeErr, Y00, Y00Err);
	if (Y00 == 0.0) 
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
			"Error: The current table of Dunham coefficients doesn't contain all coefficients necessary for calculating Y_00.");
		return;
	}
	if (Y00Err == 0.0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
			"The calculated values are Y_00=" + QString::number(Y00) + " and T_e=" + QString::number(Te) + '.');
		return;
	}
	n = int(ceil(-log10(Y00Err)));
	m = int(ceil(-log10(TeErr)));
	QMessageBox::information(this, "MolSpectAnalysis",
		"The calculated values are Y_00=" + QString::number(Y00, 'f', n) + '(' 
		+ QString::number(int(round(Y00Err * pow(10.0, n)))) + ") and T_e=" + QString::number(Te, 'f', m) + '('
		+ QString::number(int(round(TeErr * pow(10.0, m)))) + ").");
}

void MainWindow::exportTFDun()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numDunTables ? W != dunTables[i] : false); i++) ;
	if (i == numDunTables)
	{
		printf("MainWindow::exportTFDun() error: The active window doesn't belong to a table of Dunham coefficients!\n");
		return;
	}
	dunTables[i]->exportTF();
}

void MainWindow::exportPicture()
{
	DiagWindow *D = dynamic_cast<DiagWindow*>(workspace->activeSubWindow()->widget());
	if (D!=0) D->exportPicture();
}

void MainWindow::exportTableData()
{
	TableWindow *T = dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget());
	if (T == 0) return;
	QString Di = (TableDataDir.isEmpty() ? Dir : TableDataDir), F;
	if ((F = QFileDialog::getSaveFileName(this, "Select file name", Di, 
		                           "Text files (*.dat)")).isEmpty()) return;
	TableDataDir = F.left(F.lastIndexOf(QRegExp("[\\/]")));
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	D->setWindowTitle("Data to export:");
	QRadioButton *SC = new QRadioButton("Selected cells", D);
	L->addWidget(SC, 0, 0);
	QRadioButton *CT = new QRadioButton("Complete table", D);
	L->addWidget(CT, 0, 1);
	QCheckBox *EC = new QCheckBox("Exchange rows and columns", D);
	L->addWidget(EC, 1, 0, 1, 2);
	SC->setChecked(true);
	EC->setChecked(false);
	L->setRowMinimumHeight(2, 20);
	QPushButton *OK = new QPushButton("OK", D), *Ca = new QPushButton("Cancel", D);
	L->addWidget(OK, 3, 0);
	L->addWidget(Ca, 3, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Ca, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) 
		T->exportTableData(F, SC->isChecked(), EC->isChecked());
	delete D;
}

void MainWindow::exportWaveFunction()
{
	Potential *P = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
    if (P != 0) P->exportWaveFunction(NumPoints);
}

void MainWindow::exportLineProfile()
{
    ExportLineProfileDialog* Dialog = new ExportLineProfileDialog(activeSpectrum(), this);
    QMdiSubWindow* subWindow = new QMdiSubWindow;
    subWindow->setWidget(Dialog);
    workspace->addSubWindow(subWindow);
    connect(Dialog, SIGNAL(closeThis()), subWindow, SLOT(close()));
    Dialog->show();
}

void MainWindow::UpdateDunham()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numDunTables ? W != dunTables[i] : false); i++) ;
	if (i == numDunTables)
	{
		printf("MainWindow::UpdateDunham() error: The active window doesn't belong to a table of Dunham coefficients!\n");
		return;
	}
	dunTables[i]->Fit();
}

void MainWindow::ImproveDunham()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numDunTables ? W != dunTables[i] : false); i++) ;
	if (i == numDunTables)
	{
		printf("MainWindow::ImproveDunham() error: The active window doesn't belong to a table of Dunham coefficients!\n");
		return;
	}
	dunTables[i]->Fit(3, 2, false, 47, 79, 17, 219);
}

void MainWindow::fitBorderLine()
{
	DunTable *dunT = dynamic_cast<DunTable*>(workspace->activeSubWindow()->widget());
	QDialog *D = new QDialog(this);
	D->setWindowTitle("Set the fit border line points");
	int vp1, vp2, Jp1, Jp2;
	dunT->getLinePoints(vp1, Jp1, vp2, Jp2);
	QLineEdit *vp1E = new QLineEdit(QString::number(vp1), D), *Jp1E = new QLineEdit(QString::number(Jp1), D);
	QLineEdit *vp2E = new QLineEdit(QString::number(vp2), D), *Jp2E = new QLineEdit(QString::number(Jp2), D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("v_1:", D), 0, 0);
	L->addWidget(vp1E, 0, 1);
	L->addWidget(new QLabel("J_1:", D), 0, 2);
	L->addWidget(Jp1E, 0, 3);
	L->addWidget(new QLabel("v_2:", D), 1, 0);
	L->addWidget(vp2E, 1, 1);
	L->addWidget(new QLabel("J_2:", D), 1, 2);
	L->addWidget(Jp2E, 1, 3);
	L->setRowMinimumHeight(2, 20);
	L->addWidget(OK, 3, 0, 1, 2);
	L->addWidget(Cancel, 3, 2, 1, 2);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) 
		dunT->setLinePoints(vp1E->text().toInt(), Jp1E->text().toInt(), vp2E->text().toInt(), Jp2E->text().toInt());
	delete L;
	delete D;
}

void MainWindow::calcTermEnergiesPot()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
	  printf(
		"MainWindow::calcTermEnergiesPot() error: The active window doesn't belong to a potential!");
	  return;
	}
	TermTable *T = 0;
    potentials[i]->calcTermEnergies(T, NumPoints);
	if (T == 0) QMessageBox::information(this, tr("MolSpektAnalysis"), tr("No bound levels found!"));
}

void MainWindow::showTeRe()
{
	Potential *T = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	double Te, Re;
	T->getMinimum(Re, Te);
	QMessageBox::information(this, "MolSpectAnalysis", 
		"The potential minimum is at R_e=" + QString::number(Re, 'f', 8) + " and T_e=" + QString::number(Te, 'f', 6) + '.');
}

void MainWindow::calcRoughBe()
{
    Potential *P = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
    double Be = P->guessBe();
    QMessageBox::information(this, "MolSpektAnalysis", "The rotational constant Be is approximately " + QString::number(Be, 'f', 4) + ".");
}

void MainWindow::calcScatWaveFunc()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
	  printf("MainWindow::calcScatWaveFunc() error: The active window doesn't belong to a potential!");
	  return;
	}
	potentials[i]->calcScatWaveFunc();
}

void MainWindow::autoCalcScatLengthsPotentialSet()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
	  printf("MainWindow::calcScatWaveFunc() error: The active window doesn't belong to a potential!");
	  return;
	}
    potentials[i]->autoCalcScatLengthsPotentialSet(NumPoints);
}

void MainWindow::summarizePotInfo()
{
	int n, N, row = 0, nr, nc, r;
	QString dir = QFileDialog::getExistingDirectory(this, "Select the directory with the potentials to summaraize", Dir, 0);
	if (dir.isEmpty()) return;
	QDir D(dir);
	QStringList PotList = D.entryList(QStringList() << "*.pot" << "*FTT.dat", QDir::Files);
	if ((N = PotList.count()) == 0) 
	{
		QMessageBox::information(this, "MolSpectAnalysis", "The choosen directory contains no potentials (*.pot)!");
		return;
	}
	Potential *Pot = new Potential(this);
	TableWindow *ResultTab = createTableWindow();
	ResultTab->setWindowTitle(Dir);
	for (n=0; n<N; n++)
	{
		//printf("%s\n", D.absoluteFilePath(PotList[n]).toAscii().data());
		if (PotList[n].right(4) == ".pot") 
		{
			Pot->readData(D.absoluteFilePath(PotList[n]));
			fillSummarizePotInfoRow(row, ResultTab, Pot);
		}
		else if (PotList[n].right(7) == "FTT.dat")
		{
			QFile file(D.absoluteFilePath(PotList[n]));
			file.open(QIODevice::ReadOnly);
			QTextStream S(&file);
			QString Line = S.readLine();
			QStringList Header = Line.split('\t'), Row;
			QString **PotData = CreateQString(nr = Header.size(), 2);
			for (r=0; r < nr; r++) PotData[r][0] = Header[r];
			while (!S.atEnd())
			{
				Line = S.readLine();
				Row = Line.split('\t');
				if (Row.size() < nr) continue;
				for (r=0; r < nr; r++) PotData[r][1] = Row[r];
				Pot->setData(PotData, nr, 2);
				fillSummarizePotInfoRow(row, ResultTab, Pot, Row[nr - 5], Row[nr - 4], Row[nr - 3], Row[nr - 2], Row[nr - 1]);
			}
			Destroy(PotData, nr);
		}
	}
	ResultTab->getTabDimensions(nr, nc);
	ResultTab->setTabDimensions(row, nc);
	ResultTab->show();
	delete Pot;
}

void MainWindow::fillSummarizePotInfoRow(int& row, TableWindow* ResultTab, Potential* CPot, QString FQS, QString FQS_Bad, QString N_Bad, 
										 QString N_Bad_PAL, QString sigma)
{
	QStringList HL;
	QString **Data, *RD;
	double Re, De;
	int m, NR = 0, NC = 0, ADC = 2;
	if (!FQS.isEmpty()) ADC++;
	if (!FQS_Bad.isEmpty()) ADC++;
	if (!N_Bad.isEmpty()) ADC++;
	if (!N_Bad_PAL.isEmpty()) ADC++;
	if (!sigma.isEmpty()) ADC++;
	Data = CPot->getData(NR, NC);
	RD = new QString[NR + ADC];
	CPot->getReDe(Re, De);
	if (row == 0)
	{
		ResultTab->setTabDimensions(100, NR + ADC);
		for (m=0; m < NR; m++) HL << Data[m][0];
		HL << "R_e" << "D_e";
		if (!FQS.isEmpty()) HL << "FQS";
		if (!FQS_Bad.isEmpty()) HL << "FQS_Bad";
		if (!N_Bad.isEmpty()) HL << "N_Bad";
		if (!N_Bad_PAL.isEmpty()) HL << "N_Bad_PAL";
		if (!sigma.isEmpty()) HL << "sigma";
		ResultTab->setHorizontalHeader(HL);
	}
	else
	{
		int NTR, NTC;
		ResultTab->getTabDimensions(NTR, NTC);
		if (row >= NTR || NR + ADC > NTC) 
		{
			if (row >= NTR) NTR += 100;
			if (NR + ADC > NTC) NTC = NR + ADC;
			ResultTab->setTabDimensions(NTR, NTC);
		}
	}
	for (m=0; m < NR; m++) RD[m] = Data[m][1];
	RD[m++] = QString::number(Re, 'f', 8);
	RD[m++] = QString::number(De, 'f', 6);
	if (!FQS.isEmpty()) RD[m++] = FQS;
	if (!FQS_Bad.isEmpty()) RD[m++] = FQS_Bad;
	if (!N_Bad.isEmpty()) RD[m++] = N_Bad;
	if (!N_Bad_PAL.isEmpty()) RD[m++] = N_Bad_PAL;
	if (!sigma.isEmpty()) RD[m] = sigma;
	ResultTab->setRowData(row++, RD);
	Destroy(Data, NR);
	delete[] RD;
}

void MainWindow::calcSFQSMCSSeries()
{
	if (numMolecules == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", "Error: Molecule data has to be loaded first!");
		return;
	}
	ElState *St;
	Potential *Pot;
	QString PotDir, FDatDir;
	int NumParFits;
	double SFQSRad;
	SFQSCalcDialog *D = new SFQSCalcDialog(this, getDir(MDIChild::PotData));
	D->exec(St, Pot, PotDir, FDatDir, SFQSRad, NumParFits);
	delete D;
	if (St == 0) return;
	SFQSCalcControl *CalcC = new SFQSCalcControl(this, St, Pot, PotDir, FDatDir, SFQSRad,  NumParFits);
	workspace->addSubWindow(CalcC);
	CalcC->show();
}

void MainWindow::improvePotSeries()
{
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	ImprovePotSeriesDialog *D = new ImprovePotSeriesDialog(this, PotDir, Pot);
	ElState *St;
	FitData *Fd;
	QString Dir;
	double Threshold;
	int NParFits;
	if (D->exec(St, Fd, Dir, Threshold, NParFits))
	{
		ImprovePotSeriesControl *C = new ImprovePotSeriesControl(this, St, Pot, Fd, Dir, Threshold, NParFits);
		workspace->addSubWindow(C);
		C->show();
	}
	delete D;
}

void MainWindow::createAnaPotsFromMCSSplineSeries()
{
	CreateAnaPotSeriesFromMCSSplinePotSeriesDialog *D = new CreateAnaPotSeriesFromMCSSplinePotSeriesDialog(this);
	bool Success = false, improveAnaPots, UseSvd, UseLeveMarq;
	QString MolFN, StateN, FDDir, SPDir, APDir; 
	int MaxIt; 
	double Prec;
	while(!Success)
	{
		if (D->exec() == QDialog::Rejected) break;
		D->getData(MolFN, StateN, FDDir, SPDir, APDir, improveAnaPots, UseSvd, UseLeveMarq, MaxIt, Prec);
		Success = CreateAnaPotSeriesFromMCSSplinePotSeries(MolFN, StateN, FDDir, SPDir, APDir, improveAnaPots,
														   UseSvd, UseLeveMarq, MaxIt, Prec, 0);
	}
	delete D;
	if (Success) QMessageBox::information(this, "MolSpektAnalysis", "The current series of fits has finished.");
}

void MainWindow::monteCarloSim()
{
	MonteCarloSimDialog *D = new MonteCarloSimDialog(this);
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	ElState *S;
	FitData *FD;
	if ((S = Pot->getElState()) != 0 ? ((FD = S->getFitData()) != 0 ? FD->getNumLines() == 0 : true) : true)
	{
		QMessageBox::information(this, "MolSpectAnalysis", "Error: the potential has to be assigned to an electronic state with fitdata first!");
		delete D;
		return;
	}
	int N, P;
	double UncFact;
	QString Dir;
	if (D->exec(Dir, N, P, UncFact))
	{
        MonteCarloSimControl *C = new MonteCarloSimControl(this, Pot, Dir, N, P, UncFact, NumPoints);
		workspace->addSubWindow(C);
		C->show();
	}
	delete D;
}

void MainWindow::fitControl()
{
	MCFSettingsDialog *D = new MCFSettingsDialog(this);
	if (D->exec() == QDialog::Rejected) return;
	MonteCarloSim *Window = new MonteCarloSim(this, D);
	//ProcessView *Window = new ProcessView(this, Dir, 0);
	workspace->addSubWindow(Window);
	Window->show();
	//Window->start(0);
}

void MainWindow::showScatNodePos()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
	  printf("MainWindow::showScatNodePos() error: The active window doesn't belong to a potential!");
	  return;
	}
	Molecule *M = potentials[i]->getMolecule();
	if (M==0) 
	{
		QMessageBox::information(this, tr("QT4MolSpektAn"), 
								 tr("showScatNodePos: fatal error, no isotopic data available!"));
		return;
	}
	IsoTab *IsoT = M->getIso();
	int NI = IsoT->numIso, Mv, n;
	double **NodePos = Create(NI, cMaxv+1);
	potentials[i]->calcScatWaveFunc(false, NodePos);
	for (Mv = 0, i=0; i < NI; i++) 
	{
		for (n=0; (n <= cMaxv ? NodePos[i][n] != 0.0 : false); n++) ;
		if (n > Mv) Mv = n;
	}
	QTableWidget *Tab = new QTableWidget(Mv, NI);
	for (i=0; i<NI; i++) 
		Tab->setHorizontalHeaderItem(i, new QTableWidgetItem(QString::number(IsoT->mNumIso1[i]) 
				+ *IsoT->chSymb1 + QString::number(IsoT->mNumIso2[i]) + *IsoT->chSymb2));
	for (n=0; n < Mv; n++) for (i=0; i < NI; i++) 
			Tab->setItem(n, i, new QTableWidgetItem(QString::number(NodePos[i][n], 'f', 3)));
	workspace->addSubWindow(Tab);
	Tab->show();
	Destroy(NodePos, NI);
}

void MainWindow::fitScatNodePos()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
	  printf("MainWindow::fitScatNodePos() error: The active window doesn't belong to a potential!");
	  return;
	}
	double **NodePos = potentials[i]->fitNodePos();
	if (NodePos == 0) return;
	Molecule *M = potentials[i]->getMolecule();
	IsoTab *IsoT = M->getIso();
	int NI = IsoT->numIso, Mv, n;
	for (Mv = 0, i=0; i < NI; i++) 
	{
		for (n=0; (n <= cMaxv ? NodePos[i][n] != 0.0 : false); n++) ;
		if (n > Mv) Mv = n;
	}
	QTableWidget *Tab = new QTableWidget(Mv, NI);
	for (i=0; i<NI; i++) 
		Tab->setHorizontalHeaderItem(i, new QTableWidgetItem(QString::number(IsoT->mNumIso1[i]) 
				+ *IsoT->chSymb1 + QString::number(IsoT->mNumIso2[i]) + *IsoT->chSymb2));
	for (n=0; n < Mv; n++) for (i=0; i < NI; i++) 
			Tab->setItem(n, i, new QTableWidgetItem(QString::number(NodePos[i][n], 'f', 3)));
	workspace->addSubWindow(Tab);
	Tab->show();
	delete IsoT;
	//printf("Vor Destroy\n");
	Destroy(NodePos, NI);
}

void MainWindow::writePotData()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i== numPotentials)
	{
	  printf("MainWindow::writePotData() error: The active window doesn't belong to a potential!");
	  return;
	}
	potentials[i]->writePotData();
}

void MainWindow::showFitData()
{
	Potential *Pot = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (Pot != 0) Pot->showFitData();
}

void MainWindow::updateFitData()
{
	Potential *Pot = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (Pot != 0) Pot->updateFitData();
}

void MainWindow::fixJOffsets()
{
	FitData *NFD = CreateFitData(), *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	ElState *St = FD->getElState();
	TermTable *TT = (St != 0 ? St->getTermTable() : 0);
	if (TT == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
			"The fitdata set has to be assigned to an electronic state with an available term energy table first.");
		return;
	}
	if (NFD == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
						 "Error: The maximum number of fitdata sets for the current electronic state has already been reached.");
		return;
	}
	double ****TData = TT->getData(), UT, WC, W, *Offs;
	int NI = TT->getNumIso(), Mv = TT->getMaxv(), MJ = TT->getMaxJ(), NL, *RN, n, m, r, c;
	int *CompT, MCT;
	TT->getCompT(MCT, CompT);
	TableLine *LData;
	FD->getData(LData, NL, RN, sortByProg);
	delete[] RN;
	for (n=1, m=0; n<=NL; n++) if (n<NL ? LData[n].LTab != LData[n-1].LTab || LData[n].PN != LData[n-1].PN : true)
	{
		if (!LData[m].isTE && LData[m].Iso < NI)
		{
			c = (LData[m].FC >= 0 && LData[m].FC <= MCT ? (CompT[LData[m].FC] >= 0 ? CompT[LData[m].FC] : 0) : 0);
			for (r=m, UT = WC = 0.0; r<n; r++) 
				if (LData[r].vss >= 0 && LData[r].vss <= Mv && LData[r].Jss >= 0 && LData[r].Jss <= MJ)
			{
				WC += (W = 1.0 / (LData[r].err * LData[r].err));
				UT += ((TData[c][LData[m].Iso][LData[r].vss][LData[r].Jss] + LData[r].WN) * W);
			}
			if (WC > 0.0)
			{
				UT /= WC;
				for (r=m; r<n; r++) 
				{
					LData[r].WN = UT - LData[r].WN;
					LData[r].isTE = true;
					LData[r].vs = -1;
				}
			}
		}
		m=n;
	}
	delete[] CompT;
	NFD->setData(LData, NL);
	NFD->setName(FD->getName() + "FixedOffsets");
	QString Buffer = FD->getFileName();
	NFD->setFileName(Buffer.left(Buffer.lastIndexOf('.')) + "FixedOffsets.fdat");
	NFD->setSource("FitData " + FD->getSource());
	QStringList Names;
	FD->getSourceOffset(Names, Offs);
	NFD->setSourceOffset(Names, Offs);
	NFD->setJMax(FD->getJMax());
	NFD->setvMax(FD->getvMax());
	NFD->show();
}

void MainWindow::removeSingleLines()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	FD->removeSingleLines();
}

void MainWindow::removeDataFSource()
{
	FitData *FD = dynamic_cast<FitData*> (workspace->activeSubWindow()->widget());
	if (FD != 0) FD->removeDataFSource();
}

void MainWindow::selectDataFSource()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	QDialog *D = new QDialog(this);
	D->setWindowTitle("Select data from source...");
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("Source:", D), 0, 0);
	QComboBox *B = new QComboBox(D);
	B->addItems(FD->getSources());
	B->setEditable(false);
	L->addWidget(B, 0, 1);
	L->setRowMinimumHeight(1, 20);
	QPushButton *O = new QPushButton("OK", D), *C = new QPushButton("Cancel", D);
	L->addWidget(O, 2, 0);
	L->addWidget(C, 2, 1);
	connect(O, SIGNAL(clicked()), D, SLOT(accept()));
	connect(C, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) FD->selectDataFSource(B->currentText());
	delete D;
}

void MainWindow::addCalculatedLevels()
{
	FitData *FD = dynamic_cast<FitData*> (workspace->activeSubWindow()->widget());
	if (FD == 0) return;
	ElState *St = FD->getElState();
	Molecule *Mol = FD->getMolecule();
    int n, N = (St != 0 ? St->getNumTermTables() : 0), NumC;
	if (N == 0 || Mol == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
		 "The fidata set has to be assignend to an electronic state for which at least one table with term energies has to be available!"); 
		return;
	}
	IsoTab *Iso = Mol->getIso();
	QDialog *D = new QDialog(this);
	QGridLayout *La = new QGridLayout(D), *La2 = new QGridLayout;
	QLineEdit *vE, *JMaxE, *ErrE;
    QComboBox *TTB, *IsoB, *CompB;
	QCheckBox *MaxB = new QCheckBox("Maximum", D);
	QPushButton *OK, *Cancel;
	D->setWindowTitle("Please select levels to add");
	La->addWidget(new QLabel("Term energy table:", D), 0, 0);
	La->addWidget(TTB = new QComboBox(D), 0, 1);
	La->addWidget(new QLabel("Isotopologue:", D), 1, 0);
	La->addWidget(IsoB = new QComboBox(D), 1, 1);
    La->addWidget(new QLabel("Component:", D), 2, 0);
    La->addWidget(CompB = new QComboBox(D), 2, 1);
    La->addLayout(La2, 3, 0);
	La2->addWidget(MaxB, 0, 0);
	La2->addWidget(new QLabel("v:", D), 0, 1);
    La->addWidget(vE = new QLineEdit("0", D), 3, 1);
    La->addWidget(new QLabel("Maximum J:", D), 4, 0);
    La->addWidget(JMaxE = new QLineEdit("335", D), 4, 1);
    La->addWidget(new QLabel("Uncertainty:", D), 5, 0);
    La->addWidget(ErrE = new QLineEdit("9.99", D), 5, 1);
    La->setRowMinimumHeight(6, 20);
    La->addWidget(OK = new QPushButton("OK", D), 7, 0);
    La->addWidget(Cancel = new QPushButton("Cancel", D), 7, 1);
	for (n=0; n<N; n++) TTB->addItem(St->getTermTableName(n));
	for (n=0; n < Iso->numIso; n++) IsoB->addItem(Iso->getIsoName(n));
	IsoB->addItem("All isotopologues");
    CompB->addItem("e levels");
    if (St->getLambda() == 1)
    {
        CompB->addItem("f levels");
        CompB->addItem("All");
        NumC = 2;
    }
    else NumC = 1;
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Rejected)
	{
		delete D;
		delete Iso;
		return;
	}
    int v, vM = vE->text().toInt(), I = IsoB->currentIndex(), JM = JMaxE->text().toInt(), comp = CompB->currentIndex();
    if (NumC == 2 && comp < 2) comp = 1 - comp;
	TermTable *TT = St->getTermTable(TTB->currentIndex());
	bool Max = MaxB->isChecked(), AllIso;
	double Err = ErrE->text().toDouble();
	delete D;
	if (TT == 0) 
	{
		delete Iso;
		return;
	}
	if (I == Iso->numIso)
	{
		I = Iso->numIso - 1;
		AllIso = true;
	}
	else AllIso = false;
	if (I >= TT->getNumIso())
	{
		if (AllIso) I = TT->getNumIso() - 1;
		else
		{
			QMessageBox::information(this, "MolSpectAnalysis", "Error: Levels for the isotopologue " + Iso->getIsoName(I)
		                         + " are not available in the term energy table \"" + TT->getName() + "\"!");
			delete Iso;
			return;
		}
	}
	delete Iso;
	if (vM > TT->getMaxv()) 
	{
		if (!Max)
		{
			QMessageBox::information(this, "MolSpectAnalysis", "Error: Levels with v=" + QString::number(vM) 
								 + " are not available in the term energy table \"" + TT->getName() + "\"!");
			return;
		}
		else vM = TT->getMaxv();
	}
	if ((n = TT->getMaxJ()) < JM) JM = n;
	int JStep, c, NC = TT->getNumComp(), JStart, l, J;
    int i, CompI;
	double ****Data = TT->getData(), S = St->getS();
	TableLine *TLines = new TableLine[(Max ? NC * (vM + 1) : NC) * (JM + 1) * (AllIso ? I + 1 : 1)];
	QString FN = TT->getFileName();
    for (c=l=0; c < NumC; c++)
    {
        if (comp < 2 && comp != c) continue;
        for (i = (AllIso ? 0 : I); i<=I; i++) for (v = (Max ? 0 : vM); v <= vM; v++)
        {
            JStart = St->getJStart(i, c);
            JStep = Mol->getJStep(i);
            CompI = (NC == 2 ? c : 0);
            for (J = JStart; J <= JM; J += JStep) if (Data[CompI][i][v][J] != 0.0)
            {
                TLines[l].dev = TLines[l].DevR = 0.0;
                TLines[l].err = Err;
                TLines[l].FC = (S != 0.0 ? c : -1);
                TLines[l].File = FN;
                TLines[l].Iso = i;
                TLines[l].isTE = true;
                TLines[l].Js = J + (S == 0 ? (1 - c) : 1);
                TLines[l].Jss = J;
                TLines[l].LTab = 0;
                TLines[l].PN = -1;
                TLines[l].Row = 0;
                TLines[l].SourceN = l;
                TLines[l].vs = -1;
                TLines[l].vss = v;
                TLines[l++].WN = Data[CompI][i][v][J];
            }
        }
    }
	FD->addData(TLines, l);
	FD->setJMax(JM);
	FD->setvMax(vM);
}

void MainWindow::showNumLevels()
{
	FitData *F = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	Molecule *Mol = F->getMolecule();
    bool ef = (F->getElState() != 0 ? (F->getElState()->getLambda() > 0 ? true : false) : false);
	if (Mol == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
								 "The fit dataset has to be assinged to a molecule first!");
		return;
	}
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	D->setWindowTitle("Select type of levels");
	L->addWidget(new QLabel("Level type:", D), 0, 0);
	QComboBox *B = new QComboBox(D);
    QLineEdit* LE = new QLineEdit("-1");
    QCheckBox* efB = new QCheckBox("Distinguish between e and f levels", D);
    QCheckBox* disComp = new QCheckBox("Distinguis components", D);
    QCheckBox* disElStates = new QCheckBox("Distinguis electronic states", D);
    LE->setValidator(new QIntValidator(LE));
	B->setEditable(false);
	B->addItem("All levels");
	B->addItem("Fluorescence levels");
	B->addItem("Absorption levels");
    B->addItem("Fluorescence lines");
	L->addWidget(B, 0, 1);
	L->setRowMinimumHeight(1, 20);
    L->addWidget(new QLabel("Max v:"), 2, 0);
    L->addWidget(LE, 2, 1);
    L->addWidget(efB, 3, 0, 1, 2);
    L->addWidget(disComp, 4, 0, 1, 2);
    L->addWidget(disElStates, 5, 0, 1, 2);
    efB->setChecked(ef);
    L->setRowMinimumHeight(6, 20);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
    L->addWidget(OK, 7, 0);
    L->addWidget(Cancel, 7, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Rejected)
	{
		delete D;
		return;
	}
	IsoTab *Iso = Mol->getIso();
    int **LN, NC, nc = 0, NI, I, IC, c, m, GS, type = B->currentIndex();
    QString Title;
    switch (type)
    {
    case 0:
        Title = "Levels";
        break;
    case 1:
        Title = "Fluorescence levels";
        break;
    case 2:
        Title = "Absorption levels";
        break;
    case 3:
        Title = "Fluorescence lines";
        break;
    }
    QList<ElState*> statesList;
    F->getNumLevels(LN, NI, NC, type, LE->text().toInt(), efB->isChecked(), disComp->isChecked(), disElStates->isChecked(), &statesList);
	if (NI > Iso->numIso)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
			"Error: The fit dataset contains too much different isotopologues, thus the resulting sums over all isotopologues include not all levels!");
		NI = Iso->numIso;
	}
	QTableWidget *W = new QTableWidget(Iso->numIso + 1, NC + 1);
    W->setWindowTitle(Title + " FitData " + F->getName());
    QStringList ColumnHeads;
    if (disElStates->isChecked())
    {
        nc = NC / statesList.count();
        for (m=0; m < statesList.count(); ++m) for (c=0; c < nc; ++c) ColumnHeads << statesList[m]->getName();
    }
    if (disComp->isChecked())
    {
        if (ColumnHeads.isEmpty()) for (c=0; c < NC; c++) ColumnHeads << QString::number(c);
        else for (m=0; m < statesList.count(); ++m) for (c=0; c < nc; ++c) ColumnHeads[m * nc + c] += ' ' + QString::number(c);
    }
    if (efB->isChecked())
    {
        if (ColumnHeads.isEmpty()) ColumnHeads << "e" << "f";
        else for (c=0; c < NC; c+=2)
        {
            ColumnHeads[c] += " e";
            ColumnHeads[c+1] += " f";
        }
    }
    ColumnHeads << "Sum";
    W->setHorizontalHeaderLabels(ColumnHeads);
	for (I=0; I < Iso->numIso; I++)
		W->setVerticalHeaderItem(I, new QTableWidgetItem(Iso->getIsoName(I)));
	W->setVerticalHeaderItem(Iso->numIso, new QTableWidgetItem("Sum"));
	for (I = GS = 0; I < NI; I++)
	{
		for (c = IC = 0; c < NC; c++)
		{
			W->setItem(I, c, new QTableWidgetItem(QString::number(LN[I][c])));
			IC += LN[I][c];
		}
		W->setItem(I, NC, new QTableWidgetItem(QString::number(IC)));
		GS += IC;
	}
	for (c=0; c < NC; c++)
	{
		for (I = IC = 0; I < NI; I++) IC += LN[I][c];
		W->setItem(Iso->numIso, c, new QTableWidgetItem(QString::number(IC)));
	}
	for (I = NI; I < Iso->numIso; I++) for (c=0; c <= NC; c++) W->setItem(I, c, new QTableWidgetItem("0"));
	W->setItem(Iso->numIso, NC, new QTableWidgetItem(QString::number(GS)));
	workspace->addSubWindow(W);
	W->show();
	delete Iso;
	Destroy(LN, NI);
}

void MainWindow::showUncertaintyStats()
{
	FitData *D = dynamic_cast<FitData*> (workspace->activeSubWindow()->widget());
	if (D==0) return;
	QList<double> Unc;
	QList<int> Count;
	D->getUncertaintyStats(Unc, Count);
	QTableWidget *W = new QTableWidget(Unc.count(), 2);
	W->setWindowTitle("Uncertainties FitData " + D->getName());
	W->setHorizontalHeaderLabels(QStringList() << "Uncertainty" << "Count");
	int n, N = Unc.count();
	for (n=0; n<N; n++)
	{
		W->setItem(n, 0, new QTableWidgetItem(QString::number(Unc[n])));
		W->setItem(n, 1, new QTableWidgetItem(QString::number(Count[n])));
	}
	workspace->addSubWindow(W);
	W->show();
}

void MainWindow::changeSourceOffset()
{
	FitData *D = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	QStringList SourceNames;
	double *Offsets;
	D->getSourceOffset(SourceNames, Offsets);
	if (SourceNames.count() == 0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", "Error: the fit dataset does not contain any source information");
		return;
	}
	SourceOffsetDialog *SD = new SourceOffsetDialog(SourceNames, Offsets, this);
	int accepted = SD->exec();
	delete SD;
	if (accepted == QDialog::Rejected) delete[] Offsets;
	else D->setSourceOffset(SourceNames, Offsets);
}

void MainWindow::exportPotPoints()
{
	Potential *Pot = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (Pot == 0)
	{
	  printf("MainWindow::exportPotPoints() error: The active window doesn't belong to a potential!");
	  return;
	}
	int n;
	IsoTab *Iso = Pot->getIsoT();
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QLabel *L1 = new QLabel("Min R [A]:", D);
	L->addWidget(L1, 0, 0);
	QLineEdit *LRMin = new QLineEdit(QString::number(rmin, 'f', 1), D);
	L->addWidget(LRMin, 0, 1);
	QLabel *L2 = new QLabel("Max R [A]:", D);
	L->addWidget(L2, 1, 0);
	QLineEdit *LRMax = new QLineEdit(QString::number(rmax, 'f', 1), D);
	L->addWidget(LRMax, 1, 1);
	QLabel *L3 = new QLabel("Number of points:", D);
	L->addWidget(L3, 2, 0);
    QLineEdit *LNumPoints = new QLineEdit(QString::number(NumPoints), D);
	L->addWidget(LNumPoints, 2, 1);
	L->addWidget(new QLabel("Number of digits:", D), 3, 0);
	QLineEdit *LNumDigits = new QLineEdit("4", D);
	L->addWidget(LNumDigits, 3, 1);
	QCheckBox *FirstRow = new QCheckBox("First row with N and R_min", D);
	L->addWidget(FirstRow, 4, 0, 1, 2);
	FirstRow->setChecked(false);
	QCheckBox *ADCorr = new QCheckBox("Aditional column with AdCorr E", D);
	L->addWidget(ADCorr, 5, 0, 1, 2);
	ADCorr->setChecked(false);
	QLabel *AdCorrIsoL = new QLabel("for isotopologue:", D);
	L->addWidget(AdCorrIsoL, 6, 0);
	QComboBox *AdCorrIsoB = new QComboBox(D);
	L->addWidget(AdCorrIsoB, 6, 1);
	AdCorrIsoB->setEditable(false);
	if (Iso != 0 && Pot->isAdCorrA()) for (n=0; n < Iso->numIso; n++) AdCorrIsoB->addItem(Iso->getIsoName(n));
	else
	{
		ADCorr->setEnabled(false);
		AdCorrIsoL->setEnabled(false);
		AdCorrIsoB->setEnabled(false);
	}
	QCheckBox *UncertFromMCS = new QCheckBox("Estimate uncertainties using MCS potentials", D);
	L->addWidget(UncertFromMCS, 7, 0, 1, 2);
	UncertFromMCS->setChecked(false);
	L->setRowMinimumHeight(8, 20);
	QPushButton *BOK = new QPushButton("OK", D);
	L->addWidget(BOK, 9, 0);
	connect(BOK, SIGNAL(clicked()), D, SLOT(accept()));
	QPushButton *BCa = new QPushButton("Cancel", D);
	L->addWidget(BCa, 9, 1);
 	connect(BCa, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
		QString MCSDir, FileName = QFileDialog::getSaveFileName(this, tr("Select file for the points"), 
												Pot->getName() + "points.dat", tr("Text files (*.dat)"));
		if (UncertFromMCS->isChecked()) 
			MCSDir = QFileDialog::getExistingDirectory(this, "Please select the directory with MCS potentials", FileName);
		if (!FileName.isEmpty()) 
			Pot->writePoints(FileName, LRMin->text().toDouble(), LRMax->text().toDouble(), LNumPoints->text().toInt(), 
							 LNumDigits->text().toInt(), FirstRow->isChecked(), ADCorr->isChecked(), AdCorrIsoB->currentIndex(), MCSDir);
	}
	delete D;
}

void MainWindow::exportAsymptoticLevels()
{
	Potential *Pot = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (Pot == 0) return;
	QString FileName = QFileDialog::getSaveFileName(this, "Please select a file name for the levels", 
													Pot->getName() + "levels.dat", tr("Text files (*.dat)"));
	if (FileName.isEmpty()) return;
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("Number of v:", D), 0, 0);
	QLineEdit *Nv = new QLineEdit("10", D);
	L->addWidget(Nv, 0, 1);
	L->addWidget(new QLabel("Maximum J:", D), 1, 0);
	QLineEdit *MJ = new QLineEdit("2", D);
	L->addWidget(MJ, 1, 1);
	L->setRowMinimumHeight(2, 20);
	QPushButton *OK = new QPushButton("OK", D);
	L->addWidget(OK, 3, 0);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	QPushButton *Cancel = new QPushButton("Cancel", D);
	L->addWidget(Cancel, 3, 1);
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
    if (D->exec() == QDialog::Accepted) Pot->exportAsymptoticLevels(FileName, Nv->text().toInt(), MJ->text().toInt(), NumPoints);
	delete D;
}

void MainWindow::exportObservedLevelLists()
{
	int s, n, m, N, NStates, I, NI = 0, J, NJ, v, nv, l, NLT, NFD, NT, *FDNL = 0, *LTNL = 0;
	TableLine **FDTL = 0, **LTTL = 0;
	QWidget *W = workspace->activeSubWindow()->widget();
	Molecule *Mol = dynamic_cast<Molecule*>(W);
	FitData **FDL = 0;
	LineTable **LTL = 0;
	ElState **StL, *St;
	QString WDir = QFileDialog::getExistingDirectory(this, "MolSpectAnalysis", Dir);
	IsoTab *Iso = 0;
	bool *****Data;
	QFile File;
	QTextStream Str;
	if (WDir.isEmpty()) return;
	if (Mol != 0)
	{
		NT = Mol->getNumTransitions();
		StL = new ElState*[N = Mol->getNumStates()];
		for (s = NLT = NFD = 0; s<N; s++)
		{
			if ((n = (StL[s] = Mol->getStateP(s))->getNumFitDataSets()) == 0)
			{
				for (l=0; (l < NT ? (Mol->getTransitionP(l)->getLowerState() != StL[s] 
					 && Mol->getTransitionP(l)->getUpperState() != StL[s])
					  || Mol->getTransitionP(l)->getNumLineTables() == 0 : false); l++) ;
				if (l == NT) StL[s] = 0;
			}
			NFD += n;
		}
		for (s=0; (s<N ? StL[s] != 0 : false); s++) ;
		for (l=s+1; l<N; l++) if (StL[l] != 0) StL[s++] = StL[l];
		NStates = s;
		for (l=0; l < NT; l++) NLT += Mol->getTransitionP(l)->getNumLineTables();
		if (NFD > 0)
		{
			FDL = new FitData*[NFD];
			FDNL = new int[NFD];
			FDTL = new TableLine*[NFD];
			for (s=l=0; s < NStates; s++) 
				for (n=0, N = StL[s]->getNumFitDataSets(); n<N; n++)
					if ((FDL[l] = StL[s]->getFitData(n)) != 0)
			{
				FDL[l]->getData(FDTL[l], FDNL[l]);
				if (FDNL[l] > 0) l++;
			}
			if ((NFD = l) == 0)
			{
				delete[] FDL;
				delete[] FDNL;
				delete[] FDTL;
			}
		}
		if (NLT > 0)
		{
			LTL = new LineTable*[NLT];
			LTNL = new int[NLT];
			LTTL = new TableLine*[NLT];
			for (s=l=0; s < Mol->getNumTransitions(); s++)
				for (n=0, N = Mol->getTransitionP(s)->getNumLineTables(); n<N; n++)
					if ((LTL[l] = Mol->getTransitionP(s)->getLineTable(n)) != 0)
			{
				LTL[l]->getLines(LTTL[l], LTNL[l]);
				if (LTNL[l] > 0) l++;
			}
			if ((NLT = l) == 0)
			{
				delete[] LTL;
				delete[] LTNL;
				delete[] LTTL;
			}
		}
	}
	else
	{
		FDL = new FitData*[1];
		if ((FDL[0] = dynamic_cast<FitData*>(W)) != 0)
		{
			NLT = 0;
			NFD = 1;
			FDNL = new int[1];
			FDTL = new TableLine*[1];
			FDL[0]->getData(FDTL[0], FDNL[0]);
			StL = new ElState*[NStates = 1];
			StL[0] = FDL[0]->getElState();
			Mol = FDL[0]->getMolecule();
		}
		else
		{
			delete[] FDL;
			FDL = 0;
			LTL = new LineTable*[NLT = 1];
			LTL[0] = dynamic_cast<LineTable*>(W);
			NFD = 0;
			LTNL = new int[NLT];
			LTTL = new TableLine*[NLT];
			LTL[0]->getLines(LTTL[0], LTNL[0]);
			StL = new ElState*[NStates = 2];
			StL[0] = LTL[0]->getTransition()->getLowerState();
			StL[1] = LTL[0]->getTransition()->getUpperState();
			if (StL[0] == 0)
			{
				StL[0] = StL[1];
				NStates = 1;
			}
			else if (StL[1] == 0) NStates = 1;
			Mol = LTL[0]->getMolecule();
		}
		if (StL[0] == 0) NStates = 0;
	}
	if (Mol != 0)
	{
		Iso = Mol->getIso();
		NI = Iso->numIso;
	}
	for (l = nv = NJ = 0; l < NFD; l++) for (n=0; n < FDNL[l]; n++)
	{
		if (FDTL[l][n].vss >= nv) nv = FDTL[l][n].vss + 1;
		if (FDTL[l][n].Jss >= NJ) NJ = FDTL[l][n].Jss + 1;
	}
	for (l=0; l < NLT; l++) for (n=0; n < LTNL[l]; n++)
	{
		if (LTTL[l][n].vss >= nv) nv = LTTL[l][n].vss + 1;
		if (LTTL[l][n].Jss >= NJ) NJ = LTTL[l][n].Jss + 1;
		if (LTTL[l][n].vs >= nv) nv = LTTL[l][n].vs + 1;
		if (LTTL[l][n].Js >= NJ) NJ = LTTL[l][n].Js + 1;
	}
	if (NStates > 0 && NI > 0 && nv > 0 && NJ > 0)
	{
		Data = CreateBool(NStates, 4, NI, nv, NJ);
		for (s=0; s < NStates; s++) for (n=0; n<4; n++) for (v=0; v < nv; v++)
			for (I=0; I < NI; I++) for (J=0; J < NJ; J++) Data[s][n][I][v][J] = false;
		for (l=0; l < NFD; l++) if (FDNL[l] > 0) 
		{
			if ((St = FDL[l]->getElState()) != 0)
			{
				for (s=0; StL[s] != St; s++) ;
				for (n=0; n < FDNL[l]; n++) 
					if ((v = FDTL[l][n].vss) >= 0 && (J = FDTL[l][n].Jss) >= 0 
						&& (I = FDTL[l][n].Iso) >= 0)
				{
					m = (FDTL[l][n].Jss == FDTL[l][n].Js ? 1 : 0);
					if (FDTL[l][n].LTab == 0) m+=2;
					Data[s][m][I][v][J] = true;
				}
			}
			delete[] FDTL[l];
		}
		for (l=0; l < NLT; l++) if (LTNL[l] > 0)
		{
			if ((St = LTL[l]->getTransition()->getLowerState()) != 0)
			{
				for (s=0; StL[s] != St; s++) ;
				for (n=0; n < LTNL[l]; n++)
					if ((v = LTTL[l][n].vss) >= 0 && (J = LTTL[l][n].Jss) >= 0
						&& (I = LTTL[l][n].Iso) >= 0) Data[s][0][I][v][J] = true;
			}
			if ((St = LTL[l]->getTransition()->getUpperState()) != 0)
			{
				for (s=0; StL[s] != St; s++) ;
				for (n=0; n < LTNL[l]; n++)
					if ((v = LTTL[l][n].vs) >= 0 && (J = LTTL[l][n].Js) >= 0
						&& (I = LTTL[l][n].Iso) >= 0)
							Data[s][(LTTL[l][n].Jss == LTTL[l][n].Js ? 1 : 0)][I][v][J] = true;
			}
			delete[] LTTL[l];
		}
		if (WDir.right(1) != DIRSEP) WDir += DIRSEP;
		for (s=0; s < NStates; s++) for (m=0; m<4; m++) for (I=0; I < NI; I++)
		{
			for (v=0, J = NJ; v < nv && J == NJ; v++) 
				for (J=0; (J < NJ ? !Data[s][m][I][v][J] : false); J++) ;
			if (v < nv || J < NJ)
			{
				File.setFileName(WDir + "Level" + StL[s]->getName() + (m<2 ? "LIF" : "Abs") + (Iso != 0 ? Iso->getIsoName(I) : "") 
									+ (m==0 || m==2 ? "e.dat" : "f.dat"));
				File.open(QIODevice::WriteOnly);
				Str.setDevice(&File);
				Str << "v\tJ\n";
				for (v=0; v < nv; v++) for (J=0; J < NJ; J++) if (Data[s][m][I][v][J]) 
					Str << QString::number(v) << '\t' << QString::number(J) << '\n';
				File.close();
			}
		}
		Destroy(Data, NStates, 4, NI, nv);
	}
	else if (NStates == 0) 
		QMessageBox::information(this, "MolSpectAnalysis", "Error: The table has to be assinged to a molecule!");
	else QMessageBox::information(this, "MolSpectAnalysis", "Error: No data available to export!");
	if (NFD > 0)
	{
		delete[] FDL;
		delete[] FDNL;
		delete[] FDTL;
	}
	if (NLT > 0)
	{
		delete[] LTL;
		delete[] LTTL;
		delete[] LTNL;
	}
	delete[] StL;
	if (Iso != 0) delete Iso;
}

void MainWindow::plotPotential()
{
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	PotentialPlot *plot = new PotentialPlot(Pot, this);
	workspace->addSubWindow(plot);
	plot->show();
} 

void MainWindow::selPotCoeff()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; i < numPotentials; i++) 
	{
		if (W == potentials[i]) potentials[i]->setCoefficients();
		return;
	}
	printf("MainWindow::plotPotential() error: The active window doesn't belong to a potential!");
}

void MainWindow::setEnergyOffset()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::setEnergyOffset() error: The active window doesn't belong to a potential!");
		return;
	}
	QDialog *D = new QDialog(this);
	D->setWindowTitle("Change energy offset");
	QGridLayout *L = new QGridLayout(D);
	QRadioButton *RA = new QRadioButton("Set asymptote to...", D);
	RA->setChecked(true);
	L->addWidget(RA, 0, 0, 1, 2);
	QRadioButton *RM = new QRadioButton("Set potential minimum to...", D);
	L->addWidget(RM, 1, 0, 1, 2);
	QLabel *La = new QLabel("new energy:", D);
	L->addWidget(La, 2, 0);
	QLineEdit *E = new QLineEdit(D);
	E->setText("0");
	L->addWidget(E, 2, 1);
	L->setRowMinimumHeight(3, 20);
	QPushButton *OKB = new QPushButton("OK", D);
	L->addWidget(OKB, 4, 0);
	connect(OKB, SIGNAL(pressed()), D, SLOT(accept()));
	QPushButton *CancelB = new QPushButton("Cancel", D);
	L->addWidget(CancelB, 4, 1);
	connect(CancelB, SIGNAL(pressed()), D, SLOT(reject()));
	if (D->exec() == 1) 
	{
		if (RA->isChecked()) potentials[i]->setAsymptote(E->text().toDouble());
		else potentials[i]->setMinimum(E->text().toDouble());
	}
	delete D;
}

void MainWindow::showFQS()
{
	int i;
	double FQS, StdDev, sigma, FQS_PAL;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::showFQS() error: The active window doesn't belong to a potential!");
		return;
	}
    FQS = potentials[i]->getFQS(NumPoints, false, &StdDev, &sigma, &FQS_PAL);
	if (FQS == -1) QMessageBox::information(this, "Potential " + potentials[i]->getName(), 
							 "Error: the potential is not assigned to a molecule!");
	else QMessageBox::information(this, "Potential " + potentials[i]->getName(), 
							 "The chiSq of the potential is " + QString::number(FQS, 'g', 6) + ", sigma=" 
							 + QString::number(sigma, 'f', 3) + (FQS_PAL > 0.0 ? ", FQS_PAL=" + QString::number(FQS_PAL, 'f', 2) : "") 
							 + " and the standard deviation " + QString::number(StdDev, 'g', 3) + " cm^-1");
}

void MainWindow::showBadList()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::showBadList() error: The active window doesn't belong to a potential!");
		return;
	}
    potentials[i]->showBadList(NumPoints);
}


void MainWindow::createMLRPot()
{
	Potential *P1 = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (P1 != 0) P1->createMLRPot();
}

void MainWindow::fixCoefficients()
{
	Potential *P = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (P != 0) P->VaryCoefficients(false);
	else 
		printf(
		 "MainWindow::fixCoefficients() error: The active window doesn't belong to a potential!");
}

void MainWindow::varyCoefficients()
{
	Potential *P = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	if (P != 0) P->VaryCoefficients(true);
	else 
		printf(
		  "MainWindow::varyCoefficients() error: The active window doesn't belong to a potential!");
}

void MainWindow::fitSplinePot()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::fitSplinePot() error: The active window doesn't belong to a potential!");
		return;
	}
	potentials[i]->FitSplinePot();
}

void MainWindow::showSFuncs()
{
	Potential *pot = dynamic_cast<Potential*> (workspace->activeSubWindow()->widget());
	int n, N, r, NPoints;
	double **S, Rmin, Rmax;
	pot->getS(Rmin = rmin, Rmax = 20.0, NPoints = 1e4, N, S);
	if (S==0)
	{
		QMessageBox::information(this, "MolSpectAnalysis", 
					"The calculation of S functions is not implemented for the current type of potential");
		return;
	}
	double **Data = Create(NPoints, 2), h = (Rmax - Rmin) / (NPoints - 1);
	DiagWindow *W = new DiagWindow(MDIChild::SimpleDiagWindow, this);
	W->setWindowTitle("S functions of potential " + pot->getName());
	workspace->addSubWindow(W);
	for (r=1, Data[0][0] = Rmin; r < NPoints; r++) Data[r][0] = Data[r-1][0] + h;
	for (n=0; n<N; n++)
	{
		for (r=0; r < NPoints; r++) 
		{
			Data[r][1] = S[n][r];
			if (isnan(S[n][r])) printf("isnan n=%d, r=%d\n", n, r);
		}
		W->addData(Data, NPoints);
	}
	Destroy(Data, NPoints);
	Destroy(S, N);
	W->show();
}

void MainWindow::fitAnaPot()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::fitAnaPot() error: The active window doesn't belong to a potential!");
		return;
	}
	potentials[i]->FitAnaPot();
}

void MainWindow::fitWRobustWeighting()
{
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	Pot->FitAnaPot(true, true);
}

void MainWindow::fitIsoMass()
{
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	Molecule *Mol = Pot->getMolecule();
	ElState *State = Pot->getElState();
	FitData *FDat = State->getFitData();
	if (FDat == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: no fitdata available!");
		return;
	}
	int i, *IsoL, NI;
	double FQS, IsoMass;
	FDat->getavIso(IsoL, NI);
	if (NI == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: no data available for the fit!");
		return;
	}
	IsoTab *Iso = Mol->getIso();
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QComboBox *IsoB = new QComboBox(D);
	QLineEdit *SE = new QLineEdit("1.0E-7", D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	D->setWindowTitle("Select isotopologue and initial step size");
	L->addWidget(new QLabel("Isotopologue:", D), 0, 0);
	L->addWidget(IsoB, 0, 1);
	L->addWidget(new QLabel("Init step size:", D), 1, 0);
	L->addWidget(SE, 1, 1);
	L->setRowMinimumHeight(2, 20);
	L->addWidget(OK, 3, 0);
	L->addWidget(Cancel, 3, 1);
	for (i=0; i < NI; i++) IsoB->addItem(Iso->getIsoName(i));
	IsoB->setEditable(false);
	SE->setValidator(new QDoubleValidator(1e-10, 1e-2, 0, SE));
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted)
	{
        FQS = Pot->FitIsoMass(IsoL[IsoB->currentIndex()], SE->text().toDouble(), IsoMass, NumPoints);
		QMessageBox::information(this, "MolSpektAnalysis", "The resulting mass of the isotopologue " 
								+ Iso->getIsoName(IsoL[IsoB->currentIndex()]) + " is " 
								+ QString::number(IsoMass, 'f', 12) + " u, the obtained chisq is "
								+ QString::number(FQS, 'g', 3) + '.');
	}
	delete[] IsoL;
	delete D;
	delete Iso;
}

void MainWindow::fitTangToenniesPot()
{
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	Pot->FitTangToenniesPot();
}

void MainWindow::scalePotential()
{
	int i;
	double Re, De;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::scalePotential() error: The active window doesn't belong to a potential!");
		return;
	}
	potentials[i]->getReDe(Re, De);
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QLabel *L1 = new QLabel("new R_e:", D);
	L->addWidget(L1, 0, 0);
	QLineEdit *NRe = new QLineEdit(D);
	NRe->setText(QString::number(Re, 'f', 4));
	L->addWidget(NRe, 0, 1);
	QLabel *L2 = new QLabel("new D_e:", D);
	L->addWidget(L2, 1, 0);
	QLineEdit *NDe = new QLineEdit(D);
	NDe->setText(QString::number(De, 'f', 4));
	L->addWidget(NDe, 1, 1);
	L->setRowMinimumHeight(2, 10);
	QPushButton *OK = new QPushButton("OK", D);
	L->addWidget(OK, 3, 0);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	QPushButton *Cancel = new QPushButton("Cancel", D);
	L->addWidget(Cancel, 3, 1);
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) 
		potentials[i]->scalePotential(NRe->text().toDouble(), NDe->text().toDouble())->show();
	delete D;
}

void MainWindow::cdConnectSR()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::cdConnectSR() error: The active window does not belong to a potential!");
		return;
	}
	potentials[i]->cdConnectSR();
}

void MainWindow::cdConnectLR()
{
	int i;
	QWidget *W = workspace->activeSubWindow()->widget();
	for (i=0; (i < numPotentials ? W != potentials[i] : false); i++) ;
	if (i == numPotentials)
	{
		printf("MainWindow::cdConnectLR() error: The active  window does not belong to a potential!");
		return;
	}
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QComboBox *B = new QComboBox(D);
	QPushButton *O = new QPushButton("OK", D), *C = new QPushButton("Cancel", D);
	B->setEditable(false);
	B->addItems(QStringList() << "C6" << "C8" << "C10" << "C12");
	D->setWindowTitle("Select first coefficient for connection");
	L->addWidget(new QLabel("Coefficient:", D), 0, 0);
	L->addWidget(B, 0, 1);
	L->setRowMinimumHeight(1, 20);
	L->addWidget(O, 2, 0);
	L->addWidget(C, 2, 1);
	connect(O, SIGNAL(clicked()), D, SLOT(accept()));
	connect(C, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) potentials[i]->cdConnectLR(2 * B->currentIndex() + 6);
	delete D;
}

void MainWindow::calcyss()
{
	Potential *P = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	P->calcyss();
}

void MainWindow::plotShowMouseCross(bool show)
{
	DiagWindow *W = dynamic_cast<DiagWindow*>(workspace->activeSubWindow()->widget());
	if (W == 0)
	{
		printf("MainWindow::plotShowMouseCross() error: The active window is not a DiagWindow!");
		return;
	}
	W->setShowMouseCross(show);
}

void MainWindow::plotAddPotential()
{
	PotentialPlot *Plot = dynamic_cast<PotentialPlot*>(workspace->activeSubWindow()->widget());
	if (Plot == 0)
	{
		printf("MainWindow::plotAddPotential() error: The active window is not a potential plot!\n");
		return;
	}
	int m, s, p, i=0;
	QStringList L;
	QString MP, B;
	int RefL[MaxPotentials][3];
	ElState *SP;
	for (p=0; p < numPotentials; p++) 
	{
		L << potentials[p]->getName();
		RefL[i][0] = RefL[i][1] = -1;
		RefL[i++][2] = p;
	}
	for (m=0; m < numMolecules && i < MaxPotentials; m++) 
		for (s=0; s < molecules[m]->getNumStates() && i < MaxPotentials; s++)
	{ 
		L << molecules[m]->getPotName(s);
		RefL[i][0] = m;
		RefL[i][1] = s;
		RefL[i++][2] = -1;
	}
	for (m=0; m < numMolecules && i < MaxPotentials; m++) 
		for (s=0; s < molecules[m]->getNumStates() && i < MaxPotentials; s++)
			for (MP = (SP = molecules[m]->getStateP(s))->getPotentialName(), p=0; 
				   p < SP->getNumPotentials() && i < MaxPotentials; p++)
			if (!SP->isPotentialLoaded(p) && (B = SP->getPotentialName(p)) != MP) 
	{
		L << B;
		RefL[i][0] = m;
		RefL[i][1] = s;
		RefL[i++][2] = p;
	}
	QDialog *D = new QDialog(this);
	QGridLayout *Lay = new QGridLayout(D);
	Lay->addWidget(new QLabel("Potential to add:", D), 0, 0);
	QComboBox *CB = new QComboBox(D);
	Lay->addWidget(CB, 0, 1);
	CB->addItems(L);
	Lay->setRowMinimumHeight(1, 20);
	QPushButton *OK = new QPushButton("OK", D);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	Lay->addWidget(OK, 2, 0);
	QPushButton *Cancel = new QPushButton("Cancel", D);
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	Lay->addWidget(Cancel, 2, 1);
	if (D->exec() == QDialog::Accepted) 
	{
		i = CB->currentIndex();
		Plot->addPotential(RefL[i][0] == -1 ? potentials[RefL[i][2]] : (RefL[i][2] == -1 ?
				molecules[RefL[i][0]]->getPot(RefL[i][1])
				: molecules[RefL[i][0]]->getStateP(RefL[i][1])->getPotential(RefL[i][2])));
	}
	delete D;
}

void MainWindow::plotClearHistory()
{
	PotentialPlot *Plot = dynamic_cast<PotentialPlot*>(workspace->activeSubWindow()->widget());
	if (Plot == 0)
	{
		printf("MainWindow::plotClearHistory() error: The active window is not a potential plot!\n");
		return;
	}
	Plot->clearHistory();
}

void MainWindow::plotPotSnapShot()
{
	PotentialPlot *Plot = dynamic_cast<PotentialPlot*>(workspace->activeSubWindow()->widget());
	if (Plot == 0)
	{
		printf("MainWindow::plotPotSnapShot() error: The active window is not a potential plot!\n");
		return;
	}
	Plot->PotSnapShot();
	plotShowHistoryAct->setChecked(false);
}

void MainWindow::plotShowDiagFuncs(bool C)
{
	PotentialPlot *Plot = dynamic_cast<PotentialPlot*>(workspace->activeSubWindow()->widget());
	if (Plot == 0)
	{
		printf(
			 "MainWindow::plotShowDiagFuncs(bool) error: The active window is not a potential plot!\n");
		return;
	}
	Plot->setShowDiagFuncs(C);
}

void MainWindow::plotShowHistory(bool C)
{
	PotentialPlot *Plot = dynamic_cast<PotentialPlot*>(workspace->activeSubWindow()->widget());
	if (Plot == 0)
	{
		printf("MainWindow::plotShowHistory(bool) error: The activeWindow is not a potential plot!\n");
		return;
	}
	Plot->setShowHistory(C);
}

void MainWindow::plotShowPoints(bool C)
{
	PotentialPlot *Plot = dynamic_cast<PotentialPlot*>(workspace->activeSubWindow()->widget());
	if (Plot == 0)
	{
		printf("MainWindow::plotShowPoints(bool) error: The activeWindow is not a potential plot!\n");
		return;
	}
	Plot->setShowPoints(C);
}

void MainWindow::plotShowMarkerLabels(bool checked)
{
	DiagWindow *W = dynamic_cast<DiagWindow*>(workspace->activeSubWindow()->widget());
	if (W != 0) W->setShowMarkerLabels(checked);
}

LineTable *MainWindow::activeLineTable()
{
	int i;
	QWidget *AW = workspace->activeSubWindow()->widget();
	for (i=0; i<numLineTables; i++) if (AW == lineTables[i]) return lineTables[i];
	printf("MainWindow::activeLineTable Error: The active Window is no LineTable!");
	return 0;
}

MDIChild *MainWindow::activeMDIChild()
{
	return dynamic_cast<MDIChild*> (workspace->activeSubWindow()->widget());
}

DiagWindow *MainWindow::getActiveDiagWindow()
{
	return dynamic_cast<DiagWindow*> (workspace->activeSubWindow()->widget());
}

void MainWindow::AssignFC()
{
	LineTable *Table = activeLineTable();
	if (Table != 0) Table->AssignFC();
}

void MainWindow::Assignvs()
{
	LineTable *Table;
	FitData *FDat;
    bool Success = false, IsoNAv, JSm;
    int i, j, n, N, NI, NJ, NC, m, **IsoT = 0, **CompTA, *M, *Nv, nv, ni, MC, CC, KC, Nef, CMJ, LAvIso = 0, FAvIso;
    Assign_vs_CompTrans *CompT;
	double ****TE, ****TE1;
	double Tol, DATol;
    TermTable **TT, *IsoLimTT = 0, *JLimTT;
	if ((Table = activeLineTable()) != 0) 
	{
		Table->Assignvs();
		return;
	}
	if ((FDat = dynamic_cast<FitData*> (workspace->activeSubWindow()->widget())) == 0) return;
	ElState *State;
	vAssignDialog *D = new vAssignDialog(this, State = FDat->getElState());
	NI = FDat->getMolecule()->getNumIso();
    int nc, *vMin, *vMax;
    bool IsoAv[NI];
    int *UsedFor = 0;
    while (!Success)
	{
		if (D->exec() == QDialog::Rejected) break;
        if (!D->getSelectedTables(TT, vMin, vMax, N, Tol, DATol)) continue;
		if (N==0) continue;
        for (n=0; n < NI; ++n) IsoAv[n] = true;
        IsoT = new int*[N];
        CompTA = new int*[N];
        M = new int[N];
        UsedFor = new int[N];
        for (n = NC = 0; n<N; NC += M[n++] + 1) TT[n]->getCompT(M[n], CompTA[n]);
        for (n = 0, NJ = TT[0]->getMaxJ(), IsoNAv = JSm = false, JLimTT = 0; n<N; n++)
		{
            if (vMin[n] < 0) vMin[n] = 0;
            if (vMax[n] < vMin[n]) vMax[n] = TT[n]->getMaxv();
            else if (vMax[n] > TT[n]->getMaxv())
            {
                QMessageBox::information(this, "MolSpektAnalysis", "The table <" + TT[n]->getName() + "> does not have that much vibrational levels available as you selected!");
                break;
            }
            for (m=0, IsoT[n] = TT[n]->getIsoT(); m < NI; m++) if (IsoT[n][m] == -1)
            {
                IsoNAv = true;
                if (IsoAv[n]) IsoLimTT = TT[n];
                IsoAv[n] = false;
            }
            if ((CMJ = TT[n]->getMaxJ()) != NJ)
            {
                JSm = true;
                if (CMJ < NJ)
                {
                    NJ = CMJ;
                    JLimTT = TT[n];
                }
                else if (JLimTT == 0) JLimTT = TT[0];
            }
		}
        if (n==N)
        {
            for (n=0; n<N; ++n) UsedFor[n] = -1;
            for (n=0; n<N; ++n) if (vMin[n] > 0)
            {
                for (m=0; m<N; ++m) if (UsedFor[m] == -1 && n!=m && TT[n]->getElState() == TT[m]->getElState() && vMin[n] > vMax[m])
                {
                    for (i=0; i < M[n]; ++i)
                    {
                        for (j=0; j < M[m] && CompTA[n][i] != CompTA[m][j]; ++j) ;
                        if (j == M[m]) break;
                    }
                    if (i < M[n]) continue;
                    UsedFor[m] = n;
                    if (vMin[m] != 0 || vMax[m] + 1 != vMin[n])
                    {
                        bool Covered[vMin[n]];
                        for (i=0; i < vMin[n]; ++i) Covered[i] = false;
                        for (j=0; j<N; ++j) if (UsedFor[j] == n || UsedFor[j] == m) for (i = vMin[j]; i <= vMax[j]; ++i) Covered[i] = true;
                        for (j=m+1; j<N; ++j) if (UsedFor[j] == -1 && TT[j]->getElState() == TT[n]->getElState() && vMin[n] > vMax[j])
                        {
                            for (i = vMin[j]; i <= vMax[j] && !Covered[i]; ++i) ;
                            if (i == vMax[j] + 1)
                            {
                                for (i=0; i <= vMax[j]; ++j) Covered[i] = true;
                                UsedFor[j] = n;
                            }
                        }
                        for (i=0; i < vMin[n] && Covered[i]; ++i) ;
                        if (i == vMin[n])
                        {
                            for (i=0; i<N; ++i) if (UsedFor[i] == m) UsedFor[i] = n;
                            break;
                        }
                        for (i=0; i<N; ++i) if (UsedFor[i] == n) UsedFor[i] = -1;
                    }
                    else
                    {
                        for (i=0; i<N; ++i) if (UsedFor[i] == m) UsedFor[i] = n;
                        break;
                    }
                }
                if (m==N)
                {
                    QMessageBox::information(this, "MolSpektAnalysis", "The choosen minimum v for table <"
                                             + TT[n]->getName() + "> is larger than 0 and the gap is not exactly covered by the validity ranges of other tables.");
                    break;
                }
            }
            if (n==N)
            {
                if (IsoNAv || JSm)
                {
                    QString WarningMSG;
                    if (JSm) WarningMSG = "The table <" + JLimTT->getName() + "> limits the maximum used J to " + QString::number(NJ) + "!\n";
                    if (IsoNAv)
                    {
                        for (FAvIso = -1, n=0; n < NI; ++n) if (IsoAv[n])
                        {
                            LAvIso = n;
                            if (FAvIso == -1) FAvIso = n;
                        }
                        WarningMSG += "The table <" + IsoLimTT->getName() + "> limits the isotopes worked on to ";
                        IsoTab *isoTab = FDat->getMolecule()->getIso();
                        for (n=0; n < NI; ++n) if (IsoAv[n])
                        {
                            if (n != FAvIso)
                            {
                                if (n == LAvIso) WarningMSG += " and ";
                                else WarningMSG += ", ";
                            }
                            WarningMSG += isoTab->getIsoName(n);
                        }
                        delete isoTab;
                        WarningMSG += "!\n";
                    }
                    if (QMessageBox::question(this, "MolSpektAnalysis", WarningMSG + "Continue?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No)
						== QMessageBox::Yes) Success = true;
                }
                else Success = true;
            }
        }
        if (!Success)
        {
            for (n=0; n<N; ++n)
            {
                delete[] IsoT[n];
                delete[] CompTA[n];
            }
            delete[] IsoT;
            IsoT = 0;
            delete[] M;
            delete[] CompTA;
            M=0;
            CompTA = 0;
            delete[] vMin;
            delete[] vMax;
            vMin = 0;
            vMax = 0;
            delete[] UsedFor;
            UsedFor = 0;
        }
	}
	delete D;
	if (!Success) 
	{
		if (IsoT != 0) delete[] IsoT;
		return;
	}
    CompT = new Assign_vs_CompTrans[NC];
    for (n = NC = 0, MC = -1; n<N; n++) if (UsedFor[n] == -1)
	{
        Nef = (TT[n]->getElState() != 0 && TT[n]->getElState()->getLambda() > 0 ? 2 : 1);
        for (m=0, CC = MC + 1; m <= M[n]; m++)
        {
            if ((KC = CompT[NC].FC = CompTA[n][m] / Nef + CC) > MC) MC = KC;
            CompT[NC].State = TT[n]->getElState();
            CompT[NC++].ef = CompTA[n][m] % Nef;
        }
	}
	for (n = nc = 0; n<N; n++) nc += (M[n] = TT[n]->getNumComp());
	TE = new double***[nc];
    int IsoTrans[NI], l, k;
    for (ni = n = 0; n < NI; ++n) IsoTrans[n] = IsoAv[n] ? ni++ : -1;
	Nv = new int[nc];
    for (n = nc = 0; n<N; n++) if (UsedFor[n] == -1)
	{
		nv = TT[n]->getMaxv() + 1;
        if (vMax[n] < nv - 1) nv = vMax[n] + 1;
		TE1 = TT[n]->getData();
        for (m=0; m < M[n]; m++)
		{
            TE[nc] = new double**[ni];
            for (i=0; i < NI; ++i) if (IsoAv[i])
            {
                TE[nc][IsoTrans[i]] = new double*[nv];
                for (j = vMin[n]; j <= vMax[n]; ++j) TE[nc][IsoTrans[i]][j] = TE1[m][IsoT[n][i]][j];
            }
			Nv[nc++] = nv;
		}
        for (l=0; l<N; ++l) if (UsedFor[l] == n)
        {
            nc -= M[n];
            TE1 = TT[l]->getData();
            for (m=0; m < M[l]; ++m, ++nc)
            {
                for (k=0; CompTA[l][k] != CompTA[n][m]; ++k) ;
                for (i=0; i < NI; ++i) if (IsoAv[i]) for (j = vMin[l]; j <= vMax[l]; ++j) TE[nc][IsoTrans[i]][j] = TE1[k][IsoT[l][i]][j];
            }
        }
	}
    FDat->Assign_v(TE, NC, NI, NJ + 1, Nv, IsoTrans, CompT, Tol, DATol);
    for (n=0; n < nc; ++n)
    {
        for (i=0; i < ni; ++i) delete[] TE[n][i];
        delete[] TE[n];
    }
    for (n=0; n<N; ++n)
    {
        delete[] IsoT[n];
        delete[] CompTA[n];
    }
	delete[] M;
	delete[] CompTA;
	delete[] TT;
	delete[] IsoT;
	delete[] TE;
	delete[] CompT;
	delete[] Nv;
    delete[] UsedFor;
    delete[] vMin;
    delete[] vMax;
}

void MainWindow::SplitLineTable()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->splitTable();
}

void MainWindow::Delete()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->Delete();
}

void MainWindow::FindBigDiff()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->FindBigDiff();
}

void MainWindow::MarkSelected()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->MarkSelected();
}

void MainWindow::TestProgressions()
{
	LineTable *Table;
    if ((Table = activeLineTable()) != 0) Table->TestProgressions(NumPoints);
}

void MainWindow::RemoveDoubled()
{
	TableWindow *Table = dynamic_cast<TableWindow*>(workspace->activeSubWindow()->widget());
	if (Table != 0) Table->RemoveDoubled();
}

void MainWindow::SelectFound()
{
	//LineTable *Table;
	//if ((Table = activeLineTable()) != 0) Table->SelectFound();
}

void MainWindow::SetError()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) 
	{
		Table->SetError();
		return;
	}
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	QDialog *D = new QDialog(this);
	D->setWindowTitle("Select uncertainty to set");
	QGridLayout *L = new QGridLayout(D);
	L->addWidget(new QLabel("Uncertainty:", D), 0, 0);
	QLineEdit *U = new QLineEdit("0.0201", D);
	L->addWidget(U, 0, 1);
	QCheckBox *M = new QCheckBox("Set as minimum", D);
	M->setChecked(true);
	L->addWidget(M, 1, 0, 1, 2);
	L->setRowMinimumHeight(2, 20);
	QPushButton *OK = new QPushButton("OK", D), *C = new QPushButton("Cancel", D);
	L->addWidget(OK, 3, 0);
	L->addWidget(C, 3, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(C, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) FD->setUncertainty(U->text().toDouble(), M->isChecked());
	delete D;
}

void MainWindow::setvs()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->setvs();
}

void MainWindow::setFC()
{
	LineTable *Table = activeLineTable();
	QDialog *D = new QDialog(this);
	QGridLayout *L = new QGridLayout(D);
	QLineEdit *FE = new QLineEdit("-1", D);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	D->setWindowTitle("Select FC");
	L->addWidget(new QLabel("New F:", D), 0, 0);
	L->addWidget(FE, 0, 1);
	L->setRowMinimumHeight(1, 20);
	L->addWidget(OK, 2, 0);
	L->addWidget(Cancel, 2, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Accepted) 
	{
		if (Table != 0) Table->setFC(FE->text());
		else
		{
			FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
			FD->setFC(FE->text().toInt());
		}
	}
	delete D;
}

void MainWindow::SetvssAscending()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->SetvssAscending();
}

void MainWindow::setPN()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->SetPN();
}

void MainWindow::ShiftIso()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->ShiftIso();
}

void MainWindow::ShiftJdown()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->ShiftJdown();
}

void MainWindow::ShiftJup()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->ShiftJup();
}

void MainWindow::Shiftvdown()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->Shiftvdown();
}

void MainWindow::Shiftvsdown()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->Shiftvsdown();
}

void MainWindow::Shiftvsup()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->Shiftvsup();
}

void MainWindow::Shiftvup()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->Shiftvup();
}

void MainWindow::ShowCalcRelInt()
{
	LineTable *Table;
    if ((Table = activeLineTable()) != 0) Table->ShowCalcRelInt(NumPoints);
}

void MainWindow::ShowGSDeviations()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->ShowGSDeviations();
}

void MainWindow::ShowUpTerm()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->ShowUpTerm();
}

void MainWindow::ShowUpTermTable()
{
	LineTable *Table;
	if ((Table = activeLineTable()) == 0) return;
	QWidget *TT = Table->ShowUpTermTable();
	if (TT == 0) return;
	workspace->addSubWindow(TT);
	TT->show();
}

void MainWindow::ShowWeakProgressions()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->ShowWeakProgressions();
}

void MainWindow::FindSimilarProgressions()
{
	LineTable *Table;
	if ((Table = activeLineTable()) == 0) return;
	Progression P = Table->getSelectedProgression();
	if (P.N == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "There is no progression selected in the current LineTable!", 
								 QMessageBox::Ok);
		return;
	}
	FSPDiag *D = new FSPDiag(this, P);
	workspace->addSubWindow(D);
	D->show();
}

void MainWindow::sortbyvs()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->sortbyvs();
}

void MainWindow::SortProg()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	if (FD != 0) FD->sortProg();
	else
	{
		LineTable *Table;
		if ((Table = activeLineTable()) != 0) Table->SortProg();
	}
}

void MainWindow::SortFPInt()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->SortFPInt();
}

void MainWindow::SortIJvP()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->SortIJvP();
}

void MainWindow::SortIvPJ()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->SortIvPJ();
}

void MainWindow::sortUpTermIvJ()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->sortUpTermIvJ();
}

void MainWindow::SortSpectrum()
{
	LineTable *Table = activeLineTable();
	if (Table != 0) Table->SortSpectrum();
}

void MainWindow::SortfRemDoubled()
{
	LineTable *Tab = dynamic_cast<LineTable*>(workspace->activeSubWindow()->widget());
	if (Tab != 0) Tab->SortfRemDoubled();
}

void MainWindow::sortIvJfF()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	if (FD != 0) FD->sortIvJF();
}

void MainWindow::sortTabByDeviation()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	FD->sortByDev();
}

void MainWindow::sortTabByDevRatio()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	FD->sortByDevRatio();
}

void MainWindow::sortTabByElState()
{
    FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
    FD->sortByLineElState();
}

void MainWindow::sortvsIvJ()
{
	FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
    FD->sortvsIvJ();
}

void MainWindow::sortTabByElStateAndProgNr()
{
    FitData *FD = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
    FD->sortByLTabAndProg();
}

void MainWindow::SortByProgNr()
{
    LineTable* LTab = dynamic_cast<LineTable*>(workspace->activeSubWindow()->widget());
    LTab->sortByProgNumber();
}

void MainWindow::TakeOnChanges()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->TakeOnChanges();
}

void MainWindow::writeTFGS()
{
	FitData *FDat = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	if (FDat == 0)
	{
		LineTable *Table;
		if ((Table = activeLineTable()) != 0) Table->WriteTFGS();
	}
	else
	{
		QString FileName = QFileDialog::getSaveFileName(this, "Please select the file name for the exported fit data", Dir, 
														"All files (*.*)");
		if (FileName.isEmpty()) return;
		Dir = FileName.left(FileName.lastIndexOf(DIRSEP));
		if (!FDat->writeTFGS(FileName)) 
			QMessageBox::information(this, "MolSpektAnalysis", 
				"Error: Was not able to write the file, please check the necessary data and choose a valid path.");
	}
}

void MainWindow::writeExpPotFitInput()
{
	FitData *Set = dynamic_cast<FitData*>(workspace->activeSubWindow()->widget());
	LineTable *Table = activeLineTable();
	QString Name;
	if (Table == 0)
	{
		if (Set == 0) return;
		Name = Set->getName();
	}
	else Name = Table->getName();
	QString fileName = QFileDialog::getSaveFileName(this, 
													"Select file name for fit input file",
													"Term" + Name + ".dat", 
													"Text files (*.dat)", 0, QFileDialog::DontConfirmOverwrite);
	if (fileName.isEmpty()) return;
	if (Set != 0) {if (!Set->writeExPotFitInput(fileName)) QMessageBox::information(this, "MolSpektAnalysis", 
								 "An error occured while writing the file!");}
	else if (Table->writeExcPotFitInput(fileName) == false)
		QMessageBox::information(this, "MolSpektAnalysis", 
								 "An error occured while writing the file!"); 
}

void MainWindow::write2AtInput()
{
	LineTable *Table;
	if ((Table = activeLineTable()) != 0) Table->W2AI();
}

IsoTab* MainWindow::selectIso()
{
	int n;
	QDialog *D = new QDialog(this);
	D->setWindowTitle("Please select the molecule");
	QGridLayout *L = new QGridLayout(D);
	QComboBox *MB = new QComboBox(D);
	MB->setEditable(false);
	for (n=0; n < numMolecules; n++) MB->addItem(molecules[n]->getName());
	L->addWidget(new QLabel("Molecule:", this), 0, 0);
	L->addWidget(MB, 0, 1);
	L->setRowMinimumHeight(1, 20);
	QPushButton *OK = new QPushButton("OK", D), *Cancel = new QPushButton("Cancel", D);
	L->addWidget(OK, 2, 0);
	L->addWidget(Cancel, 2, 1);
	connect(OK, SIGNAL(clicked()), D, SLOT(accept()));
	connect(Cancel, SIGNAL(clicked()), D, SLOT(reject()));
	if (D->exec() == QDialog::Rejected)
	{
		delete D;
		return 0;
	}
	IsoTab *R = molecules[MB->currentIndex()]->getIso();
	delete D;
	return R;
}

void MainWindow::importMTTPotentials()
{
	QString Filename = QFileDialog::getOpenFileName(this, 
					"Choose a file with MTT potential data to import", Dir, "Text files (*.dat)");
	if (Filename.isEmpty()) return;
	QFile File(Filename);
	File.open(QIODevice::ReadOnly);
	QTextStream S(&File);
	QString Source = S.readLine();
	QStringList L, Names = S.readLine().split("\t");
	int N = Names.count() - 1, n, m;
	if (N==0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: invalid file format");
		return;
	}
	int NA[N], NLRC[N], *pLRC[N], p;
	double *A[N], *LRC[N], alpha[N], alpha1[N], B[N], beta[N], gamma[N], Uinf[N], UF = hartree_cm, da0 = 1.0 / a0_Angstrom;
	bool E = false;
	for (n=0; n<N; n++)
	{
		NLRC[n] = NA[n] = 0;
		pLRC[n] = new int[20];
		A[n] = new double[20];
		LRC[n] = new double[20];
	}
	L = S.readLine().split("\t");
	if (L[0].indexOf("U_inf") >= 0)
	{
		for (n=0; n<N && n+1 < L.count(); n++) Uinf[n] = L[n+1].toDouble() * UF;
		L = S.readLine().split("\t");
	}
	if (L[0].indexOf("B") >= 0)
	{
		for (n=0; n<N && n+1 < L.count(); n++) B[n] = L[n+1].toDouble() * UF;
		L = S.readLine().split("\t");
	}
	else for (n=0; n<N; n++) B[n] = 0.0;
	if (L[0].indexOf("alpha1", 0, Qt::CaseInsensitive) >= 0)
	{
		for (n=0; n<N && n+1 < L.count(); n++) alpha1[n] = L[n+1].toDouble() * da0;
		L = S.readLine().split("\t");
	}
	else for (n=0; n<N; n++) alpha1[n] = 0.0;
	for (m=0; ((L.count() > 0 ? (p = L[0].indexOf("A")) >= 0 : false) ? 
			L[0].mid(p+1, 1).toInt() == m : false); m++, L = S.readLine().split("\t"))
	{
		if (m==20)
		{
			QMessageBox::information(this, "MolSpektAnalysis", "Error: too much A coefficients!");
			E = true;
		}
		if (m < 20) for (n=0; n<N && n+1 < L.count(); n++) A[n][m] = L[n+1].toDouble() * UF;
		UF *= da0;
	}
	if (m==0) 
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: no A coefficients found!");
		E = true;
	}
	for (n=0; n<N; n++) for (NA[n] = m; A[n][NA[n]-1] == 0.0; NA[n]--) ;
	if (L.count() > 0 ? L[0].indexOf("alpha", 0, Qt::CaseInsensitive) >= 0 : false) 
		for (n=0; n<N && n+1 < L.count(); n++) alpha[n] = L[n+1].toDouble() * da0;
	else if (!E)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: no alpha coefficient found!");
		E = true;
	}
	L = S.readLine().split("\t");
	if (L.count() > 0 ? L[0].indexOf("beta", 0, Qt::CaseInsensitive) >= 0 : false) 
		for (n=0; n<N && n+1 < L.count(); n++) beta[n] = L[n+1].toDouble() * da0;
	else if (!E)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: no beta coefficient found!");
		E = true;
	}
	L = S.readLine().split("\t");
	if (L.count() > 0 ? L[0].indexOf("gamma", 0, Qt::CaseInsensitive) >= 0 : false)
	{
		for (n=0, UF = da0 * da0; n<N && n+1 < L.count(); n++) gamma[n] = L[n+1].toDouble() * UF;
		L = S.readLine().split("\t");
	}
	else for (n=0; n<N; n++) gamma[n] = 0.0;
	for (m=0; (L.count() > 0 ? (p = L[0].indexOf("C")) >= 0 : false); 
		 m++, L = S.readLine().split("\t"))
	{
		if (m==20 && !E)
		{
			QMessageBox::information(this, "MolSpektAnalysis", "Error: too much long range coefficients!");
			E = true;
		}
		p = L[0].right(L[0].length() - p - 1).toInt();
		UF = hartree_cm * pow(a0_Angstrom, p);
		for (n=0; n<N && n+1 < L.count(); n++) if (NLRC[n] < 20)
		{
			LRC[n][NLRC[n]] = L[n+1].toDouble() * UF;
			if (LRC[n][NLRC[n]] != 0.0) 
			{
				pLRC[n][NLRC[n]] = p;
				NLRC[n]++;
			}
		}
	}
	if (m==0 && !E) 
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: no long range coefficients found!");
		E = true;
	}
	if (L.count() > 0 ? L[0].indexOf("Uinf", 0, Qt::CaseInsensitive) >= 0 : false) for (n=0; n<N && n+1 < L.count(); n++) 
		Uinf[n] = L[n+1].toDouble();
	Potential *Pot;
	for (n=0; n<N; n++)
	{
		if ((Pot = CreatePotential()) != 0)
		{
			Pot->setName(Names[n+1]);
			Pot->setSource(Source);
			Pot->setMTTCoefficients(B[n], alpha1[n], NA[n], A[n], alpha[n], beta[n], gamma[n], NLRC[n], pLRC[n], LRC[n], Uinf[n]);
			Pot->show();
		}
		else
		{
			delete[] A[n];
			delete[] pLRC[n];
			delete[] LRC[n];
		}
	}
}

void MainWindow::importCoupledPotfitOutput()
{
	QString Filename = QFileDialog::getOpenFileName(this, 
					"Choose an output file to import", Dir, "Output files (out_*)");
	if (Filename.isEmpty()) return;
    CoupledPotFitOutputImportDialog Dialog(this);
	Dir = Filename.left(Filename.lastIndexOf(QRegExp("[\\/]")) + 1);
    if (Dialog.exec() == QDialog::Rejected) return;
    int n, N, NUpdateFitData;
    Molecule* Mol;
    double Tolerancy;
    FitData** UpdateFitDataArray;
    Dialog.getResult(Mol, UpdateFitDataArray, NUpdateFitData, Tolerancy);
    IsoTab *Iso = Mol->getIso();
	if (Iso == 0) return;
	QFile File(Filename);
	File.open(QIODevice::ReadOnly);
	QTextStream S(&File);
	QString Buffer;
	while (Buffer.left(18) != " iso   vo- vu type" && Buffer.left(18) != " iso   vo st  type" && !S.atEnd()) 
		Buffer = S.readLine();
	if (S.atEnd())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: invalid file format");
		return;
	}
	int m, i1, i2, s, NL[10];
	Buffer = S.readAll();
	QStringList SL = Buffer.split('\n');
	for (n=0; n<10; n++) NL[n] = 0;
	for (N=0; (N < SL.count() ? SL[N].left(11) != " bad lines:" : false); N++)
		NL[SL[N].mid(11, 1).toInt()]++;
	TableLine *TL[10];
	int Mv[10], MJ[10];
	for (n=0; n < 10; n++) if (NL[n] > 0) 
	{
		TL[n] = new TableLine[NL[n]];
		NL[n] = Mv[n] = MJ[n] = 0;
	}
	for (n=0; n<N; n++)
	{
		s = SL[n].mid(11, 1).toInt();
		i1 = SL[n].left(3).toInt();
		i2 = SL[n].mid(3, 3).toInt();
		for (m=0; (m < Iso->numIso ? (Iso->mNumIso1[m] != i1 || Iso->mNumIso2[m] != i2)
					&& (Iso->mNumIso1[m] != i2 || Iso->mNumIso2[m] != i1) : false); m++) ;
		TL[s][NL[s]].Iso = m;
		if ((TL[s][NL[s]].vss = SL[n].mid(6, 3).toInt()) > Mv[s]) Mv[s] = TL[s][NL[s]].vss;
		TL[s][NL[s]].vs = SL[n].mid(9, 3).toInt();
        if ((TL[s][NL[s]].Jss = SL[n].mid(14, 4).toInt()) > MJ[s]) MJ[s] = TL[s][NL[s]].Jss;
		TL[s][NL[s]].Js = (SL[n][13] == 'f' ? TL[s][NL[s]].Jss : TL[s][NL[s]].Jss + 1);
        TL[s][NL[s]].WN = SL[n].mid(18, 19).toDouble();
        TL[s][NL[s]].dev = TL[s][NL[s]].WN - SL[n].mid(37, 19).toDouble();
        TL[s][NL[s]].err = SL[n].mid(56, 10).toDouble() * 0.0001;
		TL[s][NL[s]].DevR = TL[s][NL[s]].dev / TL[s][NL[s]].err;
		TL[s][NL[s]].PN = 0;
		TL[s][NL[s]].LTab = 0;
		TL[s][NL[s]].FC = 0;
        TL[s][NL[s]].isTE = true;
        TL[s][NL[s]++].isSelected = false;
	}
    if (NUpdateFitData > 0)
    {
        int *isort = utils::heapSort(TableLineSortFunctor(TL, NL), N);
        TLRef sort[N];
        for (sort[isort[0]].state = 0; NL[sort[isort[0]].state] == 0; ++sort[isort[0]].state) ;
        for (n = 1, sort[isort[0]].line = 0; n<N; ++n)
        {
            sort[isort[n]].line = sort[isort[n-1]].line + 1;
            if (sort[isort[n]].line == NL[sort[isort[n-1]].state])
            {
                sort[isort[n]].state = sort[isort[n-1]].state + 1;
                sort[isort[n]].line = 0;
            }
            else sort[isort[n]].state = sort[isort[n-1]].state;
        }
        delete[] isort;
        for (n=0; n < NUpdateFitData; ++n) UpdateFitDataArray[n]->setDev(TL, sort, N, Tolerancy);
        for (s=0; s < 10; ++s)
        {
            for (n=m=0; n < NL[s]; ++n) if (!TL[s][n].isSelected) ++m;
            if (m < NL[s])
            {
                TableLine* LB = 0;
                if (m>0)
                {
                    LB = new TableLine[m];
                    for (n = NL[s] = 0; NL[s] < m; n++) if (!TL[s][n].isSelected) LB[NL[s]++] = TL[s][n];
                }
                delete[] TL[s];
                TL[s] = LB;
				NL[s] = m;
            }
        }
    }
	FitData *FD;
	QString SName = Filename.right(Filename.length() 
									- Filename.lastIndexOf(QRegExp("[\\/]")) - 1);
	QString Name = SName.left(SName.indexOf('.'));
	for (n=0; n < 10; n++) if (NL[n] > 0)
	{
		FD = CreateFitData();
		FD->setName(Name + QString::number(n));
		FD->setSource(SName);
		FD->setJMax(MJ[n]);
		FD->setvMax(Mv[n]);
		FD->setData(TL[n], NL[n]);
		FD->show();
	}
    if (UpdateFitDataArray != 0) delete[] UpdateFitDataArray;
}

void MainWindow::importCoupledTermTable()
{
	QString Filename = QFileDialog::getOpenFileName(this, 
					"Choose an term energy file to import", Dir, "Term energy files (T5_E_Bv)");
	if (Filename.isEmpty()) return;
	Dir = Filename.left(Filename.lastIndexOf(QRegExp("[\\/]")) + 1);
	int Nv, NJ, s, c=0, i, n, I=0, v=0, J=0, NStates, NI, *Nc, R, *S, *L;
	bool **DA, singleTab;
	double Diff, D, Mu;
	Molecule *Mol;
	ElState **States;
	CTTIDialog *Dialog = new CTTIDialog(this);
	if ((R = Dialog->exec()) == QDialog::Accepted) Dialog->getData(Mol, States, NStates, Nv, NJ, singleTab);
	delete Dialog;
	Nv++;
	NJ++;
	if (R == QDialog::Rejected || NStates == 0 || Nv <= 0 || NJ <= 0) return;
	IsoTab *Iso = Mol->getIso();
	QFile F(Filename);
	F.open(QIODevice::ReadOnly);
	QTextStream St(&F);
	QStringList IData;
	St.readLine();
	NI = Iso->numIso;
	if (singleTab)
	{
		double ****Data = Create(1, NI, Nv, NJ), *****MixC = Create(1, NI, Nv, NJ, NStates);
		for (I=0; I < NI; I++) for (v=0; v < Nv; v++) for (J=0; J < NJ; J++) for (Data[0][I][v][J] = 0.0, n=0; n < NStates; n++)
			MixC[0][I][v][J][n] = 0.0;
		while (!St.atEnd())
		{
			IData = St.readLine().split(' ', QString::SkipEmptyParts);
			if (IData.count() < 5 + NStates)
			{
				if (IData[0] == "mass:") 
				{
					for (i=0, Diff = 1.0, Mu = IData[1].toDouble(); i < NI; i++)
						if ((D = fabs(Mu - Iso->redMass[i])) < Diff)
					{
						Diff = D;
						I=i;
					}
					v=0;
				}
				continue;
			}
			if (v==0) J = IData[0].toInt();
			if (v < Nv && J < NJ)
			{
				Data[0][I][v][J] = IData[2].toDouble();
				for (n=0; n < NStates; n++) MixC[0][I][v][J][n] = IData[5+n].toDouble() * 0.00001;
			}
			v++;
		}
		TermTable *TT = CreateTermTable();
		TT->setName("Mixed TermTable");
		States[0]->addTermTable(TT);
		TT->setData(Data, 1, NI, Nv - 1, NJ - 1, 0, NStates, MixC);
		TT->show();
		delete[] States;
		delete Iso;
		return;
	}
	Nc = new int[NStates];
	S = new int[NStates];
	L = new int[NStates];
	DA = new bool*[NStates];
	double db[NStates];
	double *****TED = new double****[NStates], ******MC = new double*****[NStates];
	for (s=0; s < NStates; s++)
	{
		Nc[s] = ((L[s] = States[s]->getLambda()) == 0 ? 1 : 2) 
		      * (2 * (S[s] = States[s]->getS()) + 1);
		TED[s] = Create(Nc[s], NI, Nv, NJ);
		MC[s] = Create(Nc[s], NI, Nv, NJ, NStates);
		DA[s] = new bool[Nc[s]];
		for (c=0; c < Nc[s]; c++) for (I=0, DA[s][c] = false; I < NI; I++) for (v=0; v < Nv; v++)
			for (J=0; J < NJ; J++) for(i=0, TED[s][c][I][v][J] = 0.0; i < NStates; i++)
				MC[s][c][I][v][J][i] = 0.0;
	}
	while (!St.atEnd())
	{
		IData = St.readLine().split(' ', QString::SkipEmptyParts);
		if (IData.count() < 5 + NStates)
		{
			if (IData[0] == "mass:")
			{
				for (i=0, Diff = 1.0, Mu = IData[1].toDouble(); i < NI; i++) 
					if ((D = fabs(Mu - Iso->redMass[i])) < Diff)
				{
					Diff = D;
					I = i;
				}
				continue;
			}
			else if (IData.count() == 6 && Nc[0] == 1 && Nc[1] == 2) s=c=1;
		}
		else if (IData.count() == 5 + NStates)
		{
			for (s=c=0, i=1, db[0] = D = IData[5].toDouble(); i < NStates; i++)
				if ((db[i] = IData[5+i].toDouble()) > D)
			{
				D = db[i];
				s=i;
			}
		}
		else continue;
		J = IData[0].toInt();
		v = IData[1].toInt();
		if (v < Nv && J < NJ)
		{
			if (!DA[s][c]) DA[s][c] = true;
			TED[s][c][I][v][J] = IData[2].toDouble();
			for (i=0; i < NStates; i++) MC[s][c][I][v][J][i] = db[i];
		}
	}
	for (s=0; s < NStates; s++) if (Nv > 0 && NJ > 0)
	{
		for (c=0; (c < Nc[s] ? DA[s][c] : false); c++) ;
		for (i=c; i < Nc[s]; i++) if (!DA[s][i])
		{
			Destroy(TED[s][i], NI, Nv);
			Destroy(MC[s][i], NI, Nv, NJ);
		}
		for (i=c+1; i < Nc[s]; i++) if (DA[s][i])
		{
			TED[s][c] = TED[s][i];
			MC[s][c++] = MC[s][i];
		}
		Nc[s] = c;
		if (c==0) 
			QMessageBox::information(this, "MolSpektAnalysis", 
							"Error: no data available for state " + States[s]->getName() + "!");
	}
	TermTable *TT[NStates];
	for (s=0, DA[0][0] = true; s < NStates; s++) if (Nc[s] == 0) DA[0][0] = false;
	for (i=0; i < NStates; i++) TT[i] = CreateTermTable();
	if (TT[NStates-1] != 0 && DA[0][0])
	{
		int NP[NStates];
		int ib[3];
		Perturbation *Pert[NStates];
		for (s=0; s < NStates; s++) 
		{
			TT[s]->setData(TED[s], Nc[s], NI, Nv - 1, NJ - 1);
			NP[s] = 0;
			for (c=0; c < Nc[s]; c++) for (I=0; I < NI; I++) for (v=0; v < Nv; v++)
				for (n=0; n < NStates; n++) if (s!=n)
			{
				for (i=0; i<3; i++) ib[i]=0;
				for (J=0; J < NJ; J++) if (MC[s][c][I][v][J][n] != 0.0)
				{
					ib[0] = ib[1];
					ib[1] = ib[2];
					ib[2] = J;
					if (MC[s][c][I][v][ib[0]][n] < MC[s][c][I][v][ib[1]][n]
							&& MC[s][c][I][v][ib[2]][n] < MC[s][c][I][v][ib[1]][n])
					{
						for (i=0; (i < Nv ? TED[n][0][I][i][ib[1]] < TED[s][c][I][v][ib[1]]
										  : false); i++) ;
						if ((i < Nv ? TED[n][0][I][i][ib[0]] < TED[s][c][I][v][ib[0]] 
								|| TED[n][0][I][i][ib[2]] < TED[s][c][I][v][ib[2]] : false)
							 || (i>0 ? TED[n][0][I][i-1][ib[0]] > TED[s][c][I][v][ib[0]]
								|| TED[n][0][I][i-1][ib[2]] > TED[s][c][I][v][ib[2]] : false)) NP[s]++;
					}
				}
			}
			Pert[s] = new Perturbation[NP[s]];
			NP[s] = 0;
			for (c=0; c < Nc[s]; c++) for (I=0; I < NI; I++) for (v=0; v < Nv; v++)
				for (n=0; n < NStates; n++) if (s!=n)
			{
				for (i=0; i<3; i++) ib[i]=0;
				for (J=0; J < NJ; J++) if (MC[s][c][I][v][J][n] != 0.0)
				{
					ib[0] = ib[1];
					ib[1] = ib[2];
					ib[2] = J;
					if (MC[s][c][I][v][ib[0]][n] < MC[s][c][I][v][ib[1]][n]
							&& MC[s][c][I][v][ib[2]][n] < MC[s][c][I][v][ib[1]][n])
					{
						for (i=0; (i < Nv ? TED[n][0][I][i][ib[1]] < TED[s][c][I][v][ib[1]]
										  : false); i++) ;
						if (i < Nv ? TED[n][0][I][i][ib[0]] < TED[s][c][I][v][ib[0]] : false)
						{
							Pert[s][NP[s]].Comp = c;
							Pert[s][NP[s]].Iso = I;
							Pert[s][NP[s]].v = v;
							Pert[s][NP[s]].J = ib[0];
							Pert[s][NP[s]].PComp = 0;
							Pert[s][NP[s]].Perturber = TT[n];
							Pert[s][NP[s]].PName = TT[n]->getName();
							Pert[s][NP[s]++].Pv = i;
						}
						else if (i < Nv ? TED[n][0][I][i][ib[2]] < TED[s][c][I][v][ib[2]] : false)
						{
							Pert[s][NP[s]].Comp = c;
							Pert[s][NP[s]].Iso = I;
							Pert[s][NP[s]].v = v;
							Pert[s][NP[s]].J = ib[1];
							Pert[s][NP[s]].PComp = 0;
							Pert[s][NP[s]].Perturber = TT[n];
							Pert[s][NP[s]].PName = TT[n]->getName();
							Pert[s][NP[s]++].Pv = i;
						}
						else if (i>0 ? TED[n][0][I][i-1][ib[0]] > TED[s][c][I][v][ib[0]] : false)
						{
							Pert[s][NP[s]].Comp = c;
							Pert[s][NP[s]].Iso = I;
							Pert[s][NP[s]].v = v;
							Pert[s][NP[s]].J = ib[0];
							Pert[s][NP[s]].PComp = 0;
							Pert[s][NP[s]].Perturber = TT[n];
							Pert[s][NP[s]].PName = TT[n]->getName();
							Pert[s][NP[s]++].Pv = i-1;
						}
						else if (i>0 ? TED[n][0][I][i-1][ib[2]] > TED[s][c][I][v][ib[2]] : false)
						{
							Pert[s][NP[s]].Comp = c;
							Pert[s][NP[s]].Iso = I;
							Pert[s][NP[s]].v = v;
							Pert[s][NP[s]].J = ib[1];
							Pert[s][NP[s]].PComp = 0;
							Pert[s][NP[s]].Perturber = TT[n];
							Pert[s][NP[s]].PName = TT[n]->getName();
							Pert[s][NP[s]++].Pv = i-1;
						}
					}
				}
			}
			TT[s]->setPerturbations(NP[s], Pert[s]);
		}
		for (s=0; s < NStates; s++) 
		{
			TT[s]->setName("Coupled" + States[s]->getName());
			TT[s]->setSource(Filename);
			States[s]->addTermTable(TT[s]);
			TT[s]->show();
		}
	}
	else
	{
		for (s=0; s < NStates; s++) 
		{
			if (TT[s] != 0)
			{
				delete TT[s];
				numTermTables--;
			}
			if (Nv > 0 && NJ > 0) Destroy(TED[s], Nc[s], NI, Nv);
			else delete[] TED[s];
		}
	}
	delete[] TED;
	for (s=0; s < NStates; s++) 
	{
		if (Nv > 0 && NJ > 0) Destroy(MC[s], Nc[s], NI, Nv, NJ);
		else delete[] MC[s];
		delete[] DA[s];
	}
	delete[] MC;
	delete[] DA;
	delete Iso;
	delete[] Nc;
	delete[] S;
	delete[] L;
	delete[] States;
}

void MainWindow::importCoupledWaveFunctions()
{
	Potential *Pot = CreatePotential();
	if (Pot == 0) return;
	Molecule *Mol;
	ElState **States;
	int n, m, r, c, NStates, *Iso, NIso, C, Nv, NJ, NCoeff, NChan = 0, J=0, v=0, *Components;
	QString *IsoDirs, PotFile, T5File, ErrFile;
	QFile PFile, TFile, WFile;
	QDir WFDir;
	QTextStream S;
	QStringList L;
	double *R = new double[1000];
	IsoTab *IsoT = 0;
	CSWFImportDialog Diag(this, Dir);
	while (true)
	{
		if (Diag.exec() == QDialog::Rejected) 
		{
			delete[] R;
			return;
		}
		Diag.getData(Mol, States, Components, NStates, Iso, IsoDirs, NIso, PotFile, T5File);
		if (NStates == 0)
		{
			QMessageBox::information(this, "MolSpectAnalysis", "Error: There are no electronic states selected!");
			continue;
		}
		if (NIso == 0)
		{
			QMessageBox::information(this, "MolSpectAnalysis", "Error: There are no isotopologues selected!");
			continue;
		}
		PFile.setFileName(PotFile);
		if (!PFile.exists())
		{
			QMessageBox::information(this, "MolSpectAnalysis",
			 "Error: Potential file does not exist!");
			continue;
		}
		TFile.setFileName(T5File);
		if (!TFile.exists())
		{
			QMessageBox::information(this, "MolSpectAnalysis",
					"Error: The file T5_E_Bv does not exist!");
			continue;
		}
		TFile.open(QIODevice::ReadOnly);
		S.setDevice(&TFile);
		S.readLine();
		S.readLine();
		while (L.count() != 5 + NStates && !S.atEnd()) L = S.readLine().split(' ', QString::SkipEmptyParts);
		TFile.close();
		if (L.count() != 5 + NStates)
		{
			QMessageBox::information(this, "MolSpectAnalysis", "Error: The selected number of channels is " 
				+ QString::number(NStates) + ", while the file \"T5_E_Bv\" contains data from a calculation with " 
				+ QString::number(L.count() - 5) + " channels!");
			continue;
		}
		for (n=0; n < NIso; n++)
		{
			if (IsoDirs[n].right(1) != DIRSEP) IsoDirs[n] += DIRSEP;
			WFDir.setPath(IsoDirs[n]);
			for (m=0, C = WFDir.count(); (m<C ? WFDir[m].left(3) != "wf " : false); m++) ;
			if (m==C)
			{
				QMessageBox::information(this, "MolSpectAnalysis",
					"Error: The directory " + IsoDirs[n] 
					 + "does not contain any files with wave functions of the correct type!");
				break;
			}
			for (m=n+1; m < NIso; m++) if (IsoDirs[n] == IsoDirs[m])
			{
				QMessageBox::information(this, "MolSpectAnalysis",
					"Error: The directory " + QString::number(n) 
					 + " is the same as the directory " + QString::number(m) + '!');
				break;
			}
			if (m < NIso) break;
		}
		if (n < NIso) continue;
		if (!Pot->readData(PotFile))
		{
			QMessageBox::information(this, "MolSpectAnalysis",
				"Error reading potential file \"" + PotFile + "\"!");
			continue;
		}
		for (n=m=1; n < NStates; n++) if (States[n] != States[n-1]) m++;
		if ((n = Pot->getNCoupledPotentials() + 1) < m)
		{
			QMessageBox::information(this, "MolSpectAnalysis",
				"Error: There are more electronic states selected (" + QString::number(m) 
				+ ") than the potential file contains potentials (" + QString::number(n) + ")!");
			continue;
		}
		WFDir.setPath(IsoDirs[0]);
		for (n=0; WFDir[n].left(3) != "wf "; n++) ;
		WFile.setFileName(IsoDirs[0] + WFDir[n]);
		WFile.open(QIODevice::ReadOnly);
		S.setDevice(&WFile);
		if (IsoT != 0) delete IsoT;
		IsoT = Mol->getIso();
		S.readLine();
		for (NCoeff = 0; !S.atEnd() && NCoeff < 1000; NCoeff++)
		{
			L = S.readLine().split(' ', QString::SkipEmptyParts);
			if (NCoeff == 0) NChan = L.count() - 2;
			else if (L.count() <= NChan) break;
			R[NCoeff] = L[0].toDouble() * a0_Angstrom;
		}
		if (NChan != NStates)
		{
			QMessageBox::information(this, "MolSpectAnalysis", 
				"Error: The number of selected states (" + QString::number(NStates) + 
				") is different compared to the number of channels in the wave function files (" + QString::number(NChan) + ")!");
			continue;
		}
		if (!S.atEnd()) 
			QMessageBox::information(this, "MolSpectAnalysis", 
							"The wave functions have more coefficients than expected!");
		WFile.close();
		break;
	}
	Pot->show();
	for (Nv = NJ = n = 0; n < NIso; n++)
	{
		WFDir.setPath(IsoDirs[n]);
		for (m=0, C = WFDir.count(); m<C; m++) if (WFDir[m].left(3) == "wf ")
		{
			r = WFDir[m].indexOf('_');
			if ((v = WFDir[m].mid(2, r-2).toInt()) > Nv) Nv = v;
			if ((J = WFDir[m].right(WFDir[m].length() - r - 1).toInt()) >= NJ) NJ = J+1; 
		}
	}
	double *****WF = new double****[NIso], ****MC = new double***[NIso]; 
	double ***E = new double**[NIso], MD, IM, ID;
	int ITR[IsoT->numIso];
	for (n=0; n < IsoT->numIso; n++) ITR[n] = -1;
	for (n=0; n < NIso; n++) 
		for (J=0, WF[n] = new double***[NJ], MC[n] = new double**[NJ], E[n] = new double*[NJ]; J < NJ; J++)
			for (v=0, WF[n][J] = new double**[Nv], MC[n][J] = new double*[Nv],
					E[n][J] = new double[Nv]; v < Nv; v++)
	{
		WF[n][J][v] = 0;
		MC[n][J][v] = 0;
		E[n][J][v] = 0.0;
	}
	for (n=0; n < NIso; n++)
	{
		ITR[Iso[n]] = n; 
		WFDir.setPath(IsoDirs[n]);
		for (m=0, C = WFDir.count(); m<C; m++) if (WFDir[m].left(3) == "wf ")
		{
			WFile.setFileName(IsoDirs[n] + WFDir[m]);
			r = WFDir[m].indexOf('_');
			v = WFDir[m].mid(2, r-2).toInt() - 1;
			J = WFDir[m].right(WFDir[m].length() - r - 1).toInt();
			for (c=0, WF[n][J][v] = new double*[NChan]; c < NChan; c++)
				WF[n][J][v][c] = new double[NCoeff];
			WFile.open(QIODevice::ReadOnly);
			S.setDevice(&WFile);
			S.readLine();
			for (r=0; r < NCoeff; r++)
			{
				L = S.readLine().split(' ', QString::SkipEmptyParts);
				if (L.count() <= NChan)
				{
					ErrFile = IsoDirs[n] + WFDir[m];
					break;
				}
				for (c=0; c < NChan; c++) WF[n][J][v][c][r] = L[c+1].toDouble();
			}
			WFile.close();
		}
	}
	if (!ErrFile.isEmpty()) 
		QMessageBox::information(this, "MolSpectAnalysis", "Error reading file \""
			+ ErrFile + "\"!");
	TFile.open(QIODevice::ReadOnly);
	S.setDevice(&TFile);
	S.readLine();
	for ( ; !S.atEnd(); v++)
	{
		L = S.readLine().split(' ', QString::SkipEmptyParts);
		if (L.count() < 5 + NChan)
		{
			if (L[0] == "mass:")
			{
				for (m=1, MD = ID = fabs((IM = L[1].toDouble()) - IsoT->redMass[0]); 
					 m < IsoT->numIso; m++) if ((ID = fabs(IM - IsoT->redMass[m])) < MD) 
				{
					MD = ID;
					n = ITR[m];
				}
				L = S.readLine().split(' ', QString::SkipEmptyParts);
				//for (m=0; m < L.count(); m++) printf("L[%d]=%s\n", m, L[m].toAscii().data());
				if (L.count() < 5 + NChan) continue;
				J = L[0].toInt();
				v=0;
			}
			else continue;
		}
		if (n>=0 && J < NJ && v < Nv ? WF[n][J][v] != 0 : false)
		{
			E[n][J][v] = L[2].toDouble();
			MC[n][J][v] = new double[NChan];
			for (m=0; m < NChan; m++) MC[n][J][v][m] = L[5+m].toDouble();
		}
	}
	bool Err = false;
	for (n=0; n < NIso; n++) for (J=0; J < NJ; J++) for (v=0; v < Nv; v++) if (WF[n][J][v] != 0 && E[n][J][v] == 0.0) 
		Err = true;
	if (Err) 
		QMessageBox::information(this, "MolSpectAnalysis", 
						"Error: The file T5_E_Bv does not contain energies for levels for which wave functions are available!"); 
	CoupledSineWaveFunc *SWF = new CoupledSineWaveFunc(this, Mol);
	SWF->setData(NIso, NJ, Nv, NChan, NCoeff, States, Components, Iso, R, WF, MC, E);
	SWF->setSource(PotFile);
	Pot->setCouplingData(SWF, NChan, States, Components);
	delete[] IsoDirs;
	delete IsoT;
}

void MainWindow::importDunhamAsen()
{
	QString InpFile = QFileDialog::getOpenFileName(this, "Choose a coefficient file to import", Dir,
			"Dunham coefficient files (*.dat; *.DAT)");
	if (InpFile.isEmpty()) return;
	QFile Datei(InpFile);
	if (!Datei.open(QIODevice::ReadOnly))
	{
		
		return;
	}
	int i, N=0, Type[MaxDunCoefficients], k[MaxDunCoefficients], l[MaxDunCoefficients];
	double C[MaxDunCoefficients];
	QString Buffer;
	QTextStream S(&Datei);
	DunTable *D = CreateDunTable();
	if (D==0) return;
	for (i=0; i<3; i++) S.readLine();
	Buffer = S.readLine();
	if (Buffer.size() > 22) D->setSource(Buffer.right(Buffer.size() - 22));
	while (!S.atEnd() && N < MaxDunCoefficients)
	{
		Buffer = S.readLine();
		if (Buffer.size() < 42) continue;
		if (Buffer.mid(30, 1) != "Y" || Buffer.mid(41, 1) != "X") continue;
		Type[N] = 1;
		l[N] = Buffer.mid(36, 2).toInt();
		k[N] = Buffer.mid(33, 2).toInt();
		C[N++] = Buffer.mid(1, 25).replace("D", "E").toDouble();
	}
	if (N == MaxDunCoefficients && !S.atEnd()) 
		QMessageBox::information(this, "MolSpektAnalysis", 
						"The number of Dunham coefficients is too high to get imported completely!",
	   					QMessageBox::Ok); 
	D->setData(N, Type, k, l, C);
	D->show();
}

void MainWindow::read2AtOutput()
{
	QString InpFile = QFileDialog::getOpenFileName(this, "Choose a ZweiAt-Output.dat", Dir,
			 									"Zweiat output files (*.dat)");
	if (InpFile.isEmpty()) return;
	printf("Vor Datei\n");
	QFile Datei(InpFile);
	if (!Datei.open(QIODevice::ReadOnly)) 
    {
		QMessageBox::information(this, "MolSpektAnalysis", "The file " + InpFile + 
				"can not be opened!", QMessageBox::Ok);
		return;
    }
	int i, nu = 0, nl = 0;
    QTextStream S(&Datei);
	for (i=0; i<3; i++) S.readLine();
    QString Source = S.readLine().mid(3, 80);
	for (i=4; i<24; i++) S.readLine();
    QString Buffer;
	bool upper;
    int kl[60], ll[60], tl[60], ku[60], lu[60], tu[60];
    double Kl[60], Ku[60];
    for (i=0; i<60; i++)
    {
		upper = true;
		Buffer = S.readLine();
		if (Buffer.mid(2, 1) == "Y") tu[nu] = 1;
		else if(Buffer.mid(2, 3) == "LAY") tu[nu] = 6;
		else if(Buffer.right(3).toInt() == i + 1) upper = false;
		else break;
		if (upper)
		{
			ku[nu] = Buffer.mid(5, 2).toInt();
			lu[nu] = Buffer.mid(7, 3).toInt();
		}
		else
		{
			kl[nl] = ku[nu - 1];
			ll[nl] = lu[nu - 1];
			tl[nl] = tu[nu - 1];
		}
		if (Buffer.length() == 13)
		{
			Buffer = S.readLine();
			upper = false;
			tl[nl] = tu[nu];
			kl[nl] = ku[nu];
			ll[nl] = lu[nu];
		}
		if (upper) Ku[nu++] = Buffer.mid(13, 22).toDouble();
		else Kl[nl++] = Buffer.mid(60, 22).toDouble();
    }
    printf("Nach Lesen und Löschen\n");
	if (nu > 0)
	{
		DunTable *Du = CreateDunTable();
		Du->setSource(Source);
		Du->setName("upper_" + InpFile);
		Du->setData(nu, tu, ku, lu, Ku);
		Du->show();
	}
	if (nl > 0)
	{
		DunTable *Dl = CreateDunTable();
		Dl->setSource(Source);
		Dl->setName("lower_" + InpFile);
		Dl->setData(nl, tl, kl, ll, Kl);
		Dl->show();
	}
}

Spektrum *MainWindow::activeSpectrum()
{
	int i;
	QWidget *AW = workspace->activeSubWindow()->widget();
	for (i=0; i<numSpectra; i++) if (AW == spectra[i]) return spectra[i];
	printf("MainWindow::activeSpectra Error: The active window is no spectrum!");
	return 0;
}

void MainWindow::SpectrumNormalize()
{
	int i, j=-1;
	QWidget *W = workspace->activeSubWindow()->widget();
	QDialog D(this);
	D.setWindowTitle("Select reference spectrum for the normalization");
	QGridLayout *L = new QGridLayout(&D);
	QLabel *La = new QLabel("Reference spectrum:", &D);
	L->addWidget(La, 0, 0);
	QComboBox *B = new QComboBox(&D);
	B->setEditable(false);
	for (i=0; i < numSpectra; i++)
	{
		if (W == spectra[i]) j=i;
		else B->addItem(spectra[i]->getName());
	}
	if (j==-1)
	{
		printf("MainWindow::SpectrumNormalize error: The active window is no spectrum!");
		return;
	}
	L->addWidget(B, 0, 1);
	L->setRowMinimumHeight(1, 20);
	QPushButton *O = new QPushButton("OK", &D);
	connect(O, SIGNAL(clicked()), &D, SLOT(accept()));
	L->addWidget(O, 2, 0);
	QPushButton *C = new QPushButton("Cancel", &D);
	connect(C, SIGNAL(clicked()), &D, SLOT(reject()));
	L->addWidget(C, 2, 1);
	if (D.exec() == QDialog::Rejected) return;
	if ((i = B->currentIndex()) >= j) i++;
	spectra[j]->normalize(spectra[i]);
}

void MainWindow::SpectrumAcceptAssignments()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->AcceptAssignments();
}

void MainWindow::SpectrumAutoSLP()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->AutoSLP();
}

void MainWindow::SpectrumChange_Settings()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->CSLSettings();
}

void MainWindow::SpectrumClear_Marked()
{
	DiagWindow *W = dynamic_cast<DiagWindow*> (workspace->activeSubWindow()->widget());
	if (W!=0) W->ClearMarked();
}

void MainWindow::SpectrumDisplay_Marked()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->DisplayMarkers();
}

void MainWindow::SpectrumFindPeaks()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->editFind();
}

void MainWindow::SpectrumFind_Satellites()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->FindSat();
}

void MainWindow::SpectrumFindEmissionLines()
{
	Spektrum *S = dynamic_cast<Spektrum*>(workspace->activeSubWindow()->widget());
	FELDialog *D = new FELDialog(this, S->getMinH());
	if (D->exec() != QDialog::Rejected) 
		S->FindEmissionLines(D->getMol(), D->getIso(), D->getState(), D->getMinH(), D->getTol(), D->getMLPB());
	delete D;
}

void MainWindow::SpectrumFindProgressions()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->FindProgressions();
}

void MainWindow::SpectrumMultiSLP()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->MultiSLP();
}

void MainWindow::SpectrumNext_Progression()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->FoundF();
}

void MainWindow::SpectrumPrev_Progression()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->FoundB();
}

void MainWindow::SpectrumSetLaserFrequency()
{
	Spektrum *S = activeSpectrum();
	if (S==0) return;
	QDialog SDialog(this);
	SDialog.setWindowTitle("Set laser frequency");
	SDialog.setMaximumSize(200, 130);
	SDialog.setMinimumSize(200, 130);
	QLabel L1("Please specify the new laser frequency to be set.", &SDialog);
	L1.setGeometry(10, 10, 180, 40);
	QLabel L2("New frequency:", &SDialog);
	L2.setGeometry(10, 60, 65, 20);
	QLineEdit NE(QString::number(S->GetLaserFrequency(), 'f', 4), &SDialog);
	NE.setGeometry(100, 60, 90, 20);
	NE.setAlignment(Qt::AlignRight);
	QPushButton OK("OK", &SDialog), Cancel("Cancel", &SDialog);
	OK.setGeometry(10, 100, 70, 20);
	Cancel.setGeometry(120, 100, 70, 20);
	connect(&OK, SIGNAL(clicked()), &SDialog, SLOT(accept()));
	connect(&Cancel, SIGNAL(clicked()), &SDialog, SLOT(reject()));
	if (SDialog.exec() != QDialog::Rejected) S->SetLaserFrequency(NE.text().toDouble());
}

void MainWindow::SpectrumSetType()
{
	Spektrum *S = activeSpectrum();
	if (S == 0) return;
	QDialog Dialog(this);
	Dialog.setWindowTitle("Select the type of the spectrum");
	QGridLayout *L = new QGridLayout(&Dialog);
	QLabel *L1 = new QLabel("Type of spectrum:", &Dialog);
	L->addWidget(L1, 0, 0);
	QComboBox *B = new QComboBox(&Dialog);
	B->addItem("Laser induced fluorescence");
	B->addItem("Absorption");
	B->addItem("Normalized absorption");
	B->addItem("Thermal emission");
	B->setEditable(false);
	B->setCurrentIndex(S->getType());
	L->addWidget(B, 0, 1);
	L->setRowMinimumHeight(1, 20);
	QPushButton *O = new QPushButton("OK", &Dialog);
	connect(O, SIGNAL(clicked()), &Dialog, SLOT(accept()));
	L->addWidget(O, 2, 0);
	QPushButton *C = new QPushButton("Cancel", &Dialog);
	connect(C, SIGNAL(clicked()), &Dialog, SLOT(reject()));
	L->addWidget(C, 2, 1);
	if (Dialog.exec() != QDialog::Rejected) S->setType(B->currentIndex());
}

void MainWindow::SpectrumShowAssignmentsOnTop()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->SetShowAssignmentsOnTop(SpectrumShowAssignmentsOnTopAct->isChecked());
}

void MainWindow::SpectrumShowIntensityHistogram()
{
	Spektrum *S = activeSpectrum();
	if (S==0) return;
	IntensityHistogram *I = new IntensityHistogram(S, this);
	workspace->addSubWindow(I);
	I->show();
}

void MainWindow::SpectrumShow_Found()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->ShowFound();
}

void MainWindow::SpectrumSingleSLP()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->SingleSLP();
}

void MainWindow::SpectrumTest_Transition()
{
	Spektrum *S = activeSpectrum();
    if (S!=0) S->TTransition(NumPoints);
}

void MainWindow::SpectrumFitGaussianLineProfile()
{
    Spektrum *S = activeSpectrum();
    double I, E, Width, Offset, Sigma;
    S->FitGaussianLineProfile(E, I, Width, Offset, Sigma);
    if (Sigma > 0.0)
    {
        QTableWidget* Tab = new QTableWidget(5, 2, this);
        Tab->setHorizontalHeaderLabels(QStringList() << "Name" << "Value");
        Tab->setItem(0, 0, new QTableWidgetItem("Energy [cm^-1]:"));
        Tab->setItem(0, 1, new QTableWidgetItem(QString::number(E, 'f', 4)));
        Tab->setItem(1, 0, new QTableWidgetItem("Intensity [a.u.]:"));
        Tab->setItem(1, 1, new QTableWidgetItem(QString::number(I, 'g', 3)));
        Tab->setItem(2, 0, new QTableWidgetItem("Width [cm^-1]:"));
        Tab->setItem(2, 1, new QTableWidgetItem(QString::number(Width, 'f', 4)));
        Tab->setItem(3, 0, new QTableWidgetItem("Intensity offset:"));
        Tab->setItem(3, 1, new QTableWidgetItem(QString::number(Offset, 'g', 3)));
        Tab->setItem(4, 0, new QTableWidgetItem("Sigma"));
        Tab->setItem(4, 1, new QTableWidgetItem(QString::number(Sigma, 'f', 2)));
        workspace->addSubWindow(Tab);
        Tab->show();
    }
}

void MainWindow::SpectrumCut()
{
	Spektrum *S = activeSpectrum();
	if (S==0) return;
	CutDialog *D = new CutDialog(this, S);
	workspace->addSubWindow(D);
	D->show();
}

void MainWindow::SpectrumCutAssigned()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->cutAssignedLines();
}

void MainWindow::SpectrumCutStrong()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->cutStrongLines();
}

void MainWindow::SpectrumAdd()
{
	Spektrum *S = activeSpectrum();
	if (S==0) return;
	AddDialog *D = new AddDialog(this, S);
	workspace->addSubWindow(D);
	D->show();
}

void MainWindow::SpectrumAssignBands()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->AssignBands();
}

void MainWindow::SpectrumAssignBandsByDP()
{
	Spektrum *S = dynamic_cast<Spektrum*>(workspace->activeSubWindow()->widget());
	S->assignBandByDoubletPartners();
}

void MainWindow::SpectrumFindLinesFromTable()
{
	Spektrum *S = dynamic_cast<Spektrum*>(workspace->activeSubWindow()->widget());
	SelectLineTableDialog *D = new SelectLineTableDialog(this);
	LineTable *LT = (D->exec() == QDialog::Accepted ? D->getLineTable() : 0);
	if (LT != 0) S->FindLinesFromTable(LT);
}

void MainWindow::SpectrumContinueProgressions()
{
	Spektrum *S = activeSpectrum();
	if (S!=0) S->ContinueProgressions();
}

void MainWindow::showMDIChild(QWidget* C)
{
	workspace->addSubWindow(C);
	C->show();
}

void MainWindow::showStatusText(QString Text)
{
	//printf("showStatusText: %s\n", Text.ascii());
	QStatusBar *B = statusBar();
	if (Text.isEmpty()) B->clearMessage();
	else B->showMessage(Text);
}

bool MainWindow::CreateAnaPotSeriesFromMCSSplinePotSeries(QString MolFN, QString StateN, QString FDDir, 
					QString SPDir, QString APDir, bool improveAnaPots, bool UseSvd, bool UseLeveMarq, int MaxIt, 
					double Prec, int nStart)
{
	QFileInfo FI(FDDir);
	QFile File, Logfile;
	FitData *FDat;
	if (!FDDir.isEmpty() && (!FI.isDir() || !FI.isReadable()))
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: \"" + FDDir + "\" is not an existing directory!");
		return false;
	}
	if (!improveAnaPots)
	{
		FI.setFile(SPDir);
		if (!FI.isDir() || !FI.isReadable())
		{
			QMessageBox::information(this, "MolSpektAnalysis", "Error: \"" + SPDir + "\" is not an existing directory!");
			return false;
		}
	}
	FI.setFile(APDir);
	if (!FI.isDir() || !FI.isReadable())
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: \"" + APDir + "\" is not an existing directory!");
		return false;
	}
	QDir SplinePotDir(SPDir);
	QStringList EntryList = SplinePotDir.entryList(QStringList() << "*.pot", QDir::Files);
	if (EntryList.count() == 0)
	{
		QMessageBox::information(this, "MolSpektAnalysis", "Error: The directory \"" + SPDir + "\" does not contain any potentials!");
		return false;
	}
	Molecule *Mol = getMoleculewFM(MolFN);
	if (Mol == 0) return false;
	ElState *State = Mol->getState(StateN);
	if (State == 0) return false;
	Potential *Pot = dynamic_cast<Potential*>(workspace->activeSubWindow()->widget());
	if (Pot == 0) Pot = State->getPotential();
	FDat = (FDDir.isEmpty() ? State->getFitData() : CreateFitData());
	CreateAnaPotSeriesControl *C = new CreateAnaPotSeriesControl(this, State, Pot, FDat, SPDir, FDDir, APDir, Prec, 8, MaxIt, nStart, improveAnaPots, UseSvd, UseLeveMarq);
	workspace->addSubWindow(C);
	C->show();
	return true;
}

void MainWindow::ResumeCalculation(QString Logfile)
{
	QFile LogF(Logfile);
	if (!LogF.open(QIODevice::ReadOnly)) return;
	QTextStream S(&LogF);
	QString Line = S.readLine();
	QStringList L;
	int n=0;
	if (Line.left(41) == "CreateAnaPotSeriesFromMCSSplinePotSeries(")
	{
		Line = Line.mid(41, Line.indexOf(')') - 41);
		L = Line.split(", ");
		while (!S.atEnd())
		{
			Line = S.readLine();
			if (Line.length() > 3) n = Line.mid(3, Line.indexOf(' ') - 3).toInt() + 1;
		}
		LogF.close();
		CreateAnaPotSeriesFromMCSSplinePotSeries(L[0].mid(1, L[0].length() - 2), L[1].mid(1, L[1].length() - 2), 
					L[2].mid(1, L[2].length() - 2), L[3].mid(1, L[3].length() - 2), L[4].mid(1, L[4].length() - 2),
					(L[5] == "true" ? true : false), (L[6] == "true" ? true : false), 
					(L[7] == "true" ? true : false), L[8].toInt(), L[9].toDouble(), n);
	}
}

void MainWindow::setResumeComputationLogfile(QString File)
{
	ResumeComputationLogfile = File;
	startTimer(0);
}

void MainWindow::timerEvent(QTimerEvent* event)
{
    killTimer(event->timerId());
	ResumeCalculation(ResumeComputationLogfile);
}

void MainWindow::ShowAboutWindow()
{
    About* W = new About;
    QMdiSubWindow* subW = workspace->addSubWindow(W);
    W->show();
    connect(W, SIGNAL(closeThis()), subW, SLOT(close()));
}

void MainWindow::OpenProject(QString /*FileName*/)
{
	/*
	//Achtung: Es muss noch an eT und fT gedacht werden, da sonst möglicherweise Speicherzugriffsfehler beim Öffnen eines neuen Projektes
	QString Fehlermeldung;
	QFile Datei(FileName);
	if (!Datei.open(QIODevice::ReadOnly))
	{
		Fehlermeldung = "Die Datei "+FileName+" konnte nicht geöffnet werden.";
		QMessageBox::information( this, "Application name", "Fehler beim Öffnen der Datei.",
								  Fehlermeldung);
	}
	Project = FileName;
	Q3TextStream S(&Datei);
	QString Buffer = S.readLine();
	int i, j;
	if (Buffer.left(6) == "Number")
	{
		if ((i = Buffer.find("=")) != -1) NIsoA = Buffer.right(Buffer.length() - i - 2).toInt();
		if (NIsoA > 0) MIsoA = new double[NIsoA];
		Buffer = S.readLine();
		for (j=0; j<NIsoA; j++)
		{
			if ((i = Buffer.find("=")) != -1) 
				MIsoA[j] = Buffer.right(Buffer.length() - i - 2).toDouble();
			Buffer = S.readLine();
			printf("Buffer=%s\n", Buffer.ascii());
		}
		
	if ((i = Buffer.find("=")) != -1) NIsoB = Buffer.right(Buffer.length() - i - 2).toInt();
		printf("NIsoA=%d, NIsoB=%d\n", NIsoA, NIsoB);
		if (NIsoB > 0) MIsoB = new double[NIsoB];
		Buffer = S.readLine();
		for (j=0; j<NIsoB; j++)
		{
			if ((i = Buffer.find("=")) != -1) 
				MIsoB[j] = Buffer.right(Buffer.length() - i - 2).toDouble();
			Buffer = S.readLine();
		}
	}
	if ((i = Buffer.find("="))  != -1) Dun = Buffer.right(Buffer.length() - i - 2);
	Buffer = S.readLine();
	if ((i = Buffer.find("=")) != -1) Term = Buffer.right(Buffer.length() - i - 2);
	Buffer = S.readLine();
	if ((i = Buffer.find("=")) != -1) upDun = Buffer.right(Buffer.length() - i - 2);
	Buffer = S.readLine();
	if ((i = Buffer.find("=")) != -1) upTerm = Buffer.right(Buffer.length() - i - 2);
	Buffer = S.readLine();
	if ((i = Buffer.find("=")) != -1) Lines = Buffer.right(Buffer.length() - i - 2);
	Buffer = S.readLine();
	if ((i = Buffer.find("=")) != -1) OpenSpekt = Buffer.right(Buffer.length() - i - 2);
	if (!OpenSpekt.isEmpty()) Open(OpenSpekt);
	S.readLine();
	SpektList = S.read();
	printf(SpektList.ascii());
	if (!Term.isEmpty()) OpenTermenergies(Term);
	if (!Dun.isEmpty()) 
		if (!RDunham(Dun, NI, LMAXA, KMAXA, LMAXO, KMAXO, veu, veo, Jeu, YA, YO, ELU))
	{
		Fehlermeldung = "Die Datei "+FileName+" konnte nicht geöffnet werden.";
		QMessageBox::information( this, "Application name", "Fehler beim Öffnen der Datei.",
								  Fehlermeldung);
	}
	if (!upTerm.isEmpty())
	{
		OpenTermenergies(upTerm);
		ShowTermPlot();
	}
	printf("Nach ShowTermPlot\n");
	if (!Lines.isEmpty())
	{
		ShowLineTable();
		tabelle->Open(Lines);
	}
	if (upDun == "output.dat") 
	{
		ShowLineTable();
		tabelle->Read2AtOutput();
		//if (TP != NULL) tabelle->PlotDunham();
		//printf("Nach Read2AtOutput\n");
	}
	printf("Ende von OpenProject\n");*/
}


