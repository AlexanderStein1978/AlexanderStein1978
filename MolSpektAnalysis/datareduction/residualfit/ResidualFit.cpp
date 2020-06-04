//
// C++ Implementation: ResidualFit
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2011 - 2019
//
// Copyright: See README file that comes with this source code
//
//


#include "ResidualFit.h"


ResidualFit::ResidualFit(const ElState *i_state) : cycleDetected(false), NumSplines(0), cycleIndex(0), iso(-1), v(-1), comp(-1), splines(0),
    cycleDetectList(new QList<SplinePoint*>), stateName(0), m_state(i_state), m_localPerturbations(), m_connectSplines()
{
}

ResidualFit::~ResidualFit()
{
	clear();
	clearCycleDetectList();
	delete cycleDetectList;
    if (stateName != 0) delete stateName;
}

void ResidualFit::clear()
{
	int n;
	if (NumSplines > 0)
	{
		for (n=0; n < NumSplines; n++) delete splines[n];
		delete[] splines;
		NumSplines = 0;
		splines = 0;
        ConnectSplines();
	}
    if (m_localPerturbations.size() > 0)
    {
        for (std::vector<LocalPerturbation*>::const_iterator it = m_localPerturbations.begin(); it != m_localPerturbations.end(); ++it) delete *it;
        m_localPerturbations.clear();
    }
}

void ResidualFit::FitTo(int N, Point* points)
{
	if (NumSplines > 0) clear();
	if (N==0 || points == 0) return;
	Point B;
	double step, cp;// devm2 = 0.0, devm1 = 0.0, cdev = 0.0;
	//double devp1 = 0.0;
	QList<int> jumpList, orderList;
	Point *DataPoints = 0; //new Point[N];
	int n, m, j, k=0, /*l=0, p,*/ Nsps, i;// *pos;
	SplinePoint *sps;
	//bool toImprove = true, moved, pointDeleted;
	B.x = 1.0;
	while (B.x != -1.0)
	{
		B.x = -1.0;
		for (n=0; n<N-1; n++) 
			if (points[n].x > points[n+1].x 
				|| (points[n].x == points[n+1].x && n>0 
				 && fabs(points[n].y - points[n-1].y) 
				     > fabs(points[n-1].y - points[n+1].y))) 
		{
			B = points[n];
			points[n] = points[n+1];
			points[n+1] = B;
		}
	}
	for (n=0; n<N; n++) if (points[n].unc > 9.0)
	{
		while (points[n].unc > 9.0) points[n].unc -= 9.0;
		points[n].sig = 1.0 / points[n].unc;
	}
	splines = new Spline*[NumSplines = 1];
	sps = new SplinePoint[2];
	sps[0].x = points[0].x;
	sps[1].x = points[N-1].x;
	splines[0] = new Spline(2, sps, N<10);
	printf("Step 1: Fit spline 1 of 1\n");
	splines[0]->FitTo(N, points);
    /*for (n=0; n<N; n++)
	{
		devm2 = devm1;
		devm1 = cdev;
		cdev = devp1;
		devp1 = points[n].y - splines[0]->gety(points[n].x);
		if (cdev * devm1 < 0.0 && devm1 * devm2 > 0.0 && cdev * devp1 > 0.0 
				&& (fabs(devm1) > fabs(devm2) || fabs(cdev) > fabs(devp1)) 
				&& (fabs(cdev) > 2.0 * points[n-1].unc 
					|| fabs(devm1) > 2.0 * points[n-2].unc)) 
			jumpList.append(n-1);
	}
	if (jumpList.size() > 0)
	{
		delete splines[0];
		delete[] splines;
		splines = new Spline*[NumSplines = jumpList.size() + 1];
		for (n=0; n < NumSplines; n++)
		{
			sps = new SplinePoint[2];
			sps[0].x = points[m = (n>0 ? jumpList[n-1] : 0)].x;
			sps[1].x = points[(k = (n < jumpList.size() ? jumpList[n] : N)) - 1].x;
			splines[n] = new Spline(2, sps, k-m<10);
			if ((j=k-m) == 1) 
			{
				sps[0].y = sps[1].y = points[m].y;
				sps[0].yss = sps[1].yss = 0.0;
			}
			else 
			{
				printf("Step 2: Fit spline %d of %d\n", n+1, NumSplines);
				splines[n]->FitTo(j, points + m);
			}
		}
	}
	for (n=0; n < NumSplines; n++) 
	{
		i = (n>0 ? jumpList[n-1] : 0);
		j = (n < jumpList.size() ? jumpList[n] : N);
		pointDeleted = false;
		while (!pointDeleted && 4 * (splines[n]->getNPoints() - 1) < j-i)
		{
			for (m=i, cdev = 0.0, toImprove = false; m<j; m++)
			{
				devm1 = cdev;
				cdev = points[m].y - splines[n]->gety(points[m].x);
				if (cdev * devm1 > 0.0) k++;
				else k=l=0;
				if (fabs(cdev) > points[m].unc) l++;
				//printf("k=%d, NPoints=%d\n", k, splines[n]->getNPoints());
				if (k>=4 && l>=3) 
				{
					toImprove = true;
					break;
				} 
			}
			if (!toImprove) break;
			if (splines[n]->isNatural()) splines[n]->setNatural(false);
			else
			{
				Nsps = splines[n]->getNPoints() + 1;
				*for (m=i, devm1 = 0.0; m<j; m++) 
					devm1 += fabs(points[m].y - splines[n]->gety(points[m].x));
                devm1 /= Nsps;*
				sps = new SplinePoint[Nsps];
				pos = new int[Nsps];
				sps[0].x = points[pos[0] = i].x;
				sps[Nsps - 1].x = points[pos[Nsps - 1] = j-1].x;
				*for (m=i, k=1, devp1 = devm1, cdev = 0.0; k < Nsps - 1; m++)
				{
					cdev += (devm2 = fabs(points[m].y - splines[n]->gety(points[m].x)));
					if (cdev >= devp1)
					{
						sps[k].x = points[pos[k] = (sps[k-1].x < points[m-1].x 
										&& (points[m].x == sps[Nsps - 1].x
										 || devp1 - cdev + devm2 < cdev - devp1) ?
										      m-1 : m)].x;
						devp1 += devm1;
						k++;
					}
				}
				for (m=1; m < Nsps - 1 && pos[m-1] < pos[Nsps - 1] - 5 * (Nsps - m - 1); 
																					m++)
					if (pos[m] < pos[m-1] + 5) 
						sps[m].x = points[pos[m] = pos[m-1] + 5].x;
				for (m = Nsps - 2; m>0 && pos[m+1] > pos[0] + 5*m; m--)
					if (pos[m] > pos[m+1] - 5) 
                        sps[m].x = points[pos[m] = pos[m+1] - 5].x;*
				for (m=i, devm2 = devm1 = cdev = 0.0; m<j; m++)
				{
					devp1 = cdev;
					devm1 += fabs(cdev = points[m].y - splines[n]->gety(points[m].x));
					if (cdev * devp1 > 0.0) k++;
					else 
					{
						if (k>=4 && l>=3) devm2 += devm1;
						devm1 = 0.0;
						k=l=0;
					}
					if (fabs(cdev) > points[m].unc) l++;
				}
				cp = devm2 / (Nsps - 1);
				for (m=i, devm2 = devm1 = cdev = 0.0, p=1; m<j && p < Nsps - 1; m++)
				{
					devp1 = cdev;
					devm1 += fabs(cdev = points[m].y - splines[n]->gety(points[m].x));
					if (cdev * devp1 > 0.0) k++;
					else 
					{
						if (k>=4 && l>=3) devm2 += devm1;
						devm1 = 0.0;
						k=l=0;
					}
					if (devm2 + (k>=4 && l>=3 ? devm1 : 0.0) >= cp * double(p))
					{
						sps[p].x = points[pos[p]=m].x;
						p++;
					}
					if (fabs(cdev) > points[m].unc) l++;
				}
				if (p < Nsps - 1 || pos[p] >= j - 2)
				{
					for (k=p; k < Nsps - 1; k++) sps[k].x = points[pos[k] = j + 2 * (k - Nsps)].x;
					for (k=p-1; pos[k] >= pos[k+1] - 1; k--) sps[k].x = points[pos[k] = j + 2 * (k - Nsps)].x;
				}
				clearCycleDetectList();
				for (moved = true; moved && !inCycle(Nsps, sps); ) for (m=1, moved = false; m < Nsps - 1; m++)
				{
					k = pos[m];
					if (sps[m].x - sps[m-1].x > 2 * (sps[m+1].x - sps[m].x))
					{
						//printf("sps[%d].x = %f\n", m, points[pos[m]-1].x);
						for (sps[m].x = points[--pos[m]].x;
							sps[m].x - sps[m-1].x > 2 * (sps[m+1].x - sps[m].x);
							sps[m].x = points[--pos[m]].x) ;// printf("sps[%d].x = %f\n", m, sps[m].x);
						if (2 * (sps[m].x - sps[m-1].x) < sps[m+1].x - sps[m].x
							&& sps[m+1].x - points[pos[m]+1].x > sps[m].x - sps[m-1].x)
						{
							sps[m].x = points[++pos[m]].x;
							//printf("sps[%d].x = %f\n", m, sps[m].x);
						}
					}
					else if (2 * (sps[m].x - sps[m-1].x) < sps[m+1].x - sps[m].x)
					{
						//printf("sps[%d].x = %f\n", m, points[pos[m]+1].x);
						for (sps[m].x = points[++pos[m]].x;
							2 * (sps[m].x - sps[m-1].x) < sps[m+1].x - sps[m].x;
							sps[m].x = points[++pos[m]].x) ; //printf("sps[%d].x = %f\n", m, sps[m].x);
						if (sps[m].x - sps[m-1].x > 2 * ( sps[m+1].x - sps[m].x)
							&& sps[m+1].x - sps[m].x < points[pos[m]-1].x - sps[m-1].x)
						{
							sps[m].x = points[--pos[m]].x;
							//printf("sps[%d].x = %f\n", m, sps[m].x);
						}
					}
					if (k != pos[m]) moved = true;
				}
				if (cycleDetected) setAverage(Nsps, sps);
				for (m=1; m < Nsps; m++) 
				{
					for (k=m-1; sps[k].x == -1.0; k--) ;
					for (p=k+1; p<m && sps[p].x == sps[k].x; p++) ;
					if (sps[p].x == sps[m].x) 
					{
						sps[m].x = -1.0;
						pointDeleted = true;
					}
				}
				if (sps[Nsps - 1].x == -1.0) sps[Nsps - 2].x = splines[n]->getxN(); 
				delete[] pos;
				if (pointDeleted)
				{
					for (m=1; sps[m].x != -1.0; m++) ;
					for (p=m+1; p < Nsps; p++) if (sps[p].x != -1.0) sps[m++] = sps[p];
					Nsps = m;
				}
				splines[n]->setSplinePoints(Nsps, sps);
			}
			printf("Step3: Fit spline %d of %d\n", n+1, NumSplines);
			splines[n]->FitTo(j-i, points + i);
		}
	}
	for (n=m=0; n<N; n++)
	{
		if (m < jumpList.size() && n == jumpList[m])
		{
			k = (m>0 ? jumpList[m-1] : 0);
			printf("Step4: Fit spline %d of %d\n", m, jumpList.size() + 1);
			splines[m]->FitTo(jumpList[m] - k, points + k);
			m++;
		}
		points[n].SplineNum = m;
		if ((cdev = fabs(points[n].y - splines[m]->gety(points[n].x))) 
				> 4.0 * points[n].unc)
		{
			while (cdev > points[n].unc) points[n].unc += 9.0;
			points[n].sig = 1.0 / points[n].unc;
		}
	}
	k = (m>0 ? jumpList[m-1] : 0);
	printf("Step4: Fit spline %d of %d\n", m+1, jumpList.size() + 1);
    splines[m]->FitTo(N-k, points + k);*/
	return;
	
	
	
	
	for (n=0; n <= jumpList.size(); n++) orderList.append(n);
	for (n=0; n < jumpList.size(); n++)
	{
		for (m=0; orderList[m] != n; m++) ;
		for (k=n+1; k < jumpList.size() 
			&& points[jumpList[k]].x - points[jumpList[n]-1].x <= 6.0
			&& fabs(points[jumpList[k]].y - points[jumpList[n]-1].y) 
			  > 4.0 * MAX(points[jumpList[k]].sig, points[jumpList[n]-1].sig); k++) ;
		if (k < jumpList.size() 
			&& points[jumpList[k]].x - points[jumpList[n]-1].x <= 6.0)
		{
			for (j=m+1; j < orderList.size() && orderList[j] != k+1; j++) ;
			if (j < orderList.size() && j > m+1) orderList.move(j, m+1);
		}
	}
	for (n=m=0; n < orderList.size(); n++) 
		for (j=(orderList[n] > 0 ? jumpList[orderList[n]-1] : 0); 
				j < (orderList[n] < jumpList.size() ? jumpList[orderList[n]] : N); j++) 
			DataPoints[m++] = points[j];
	jumpList.clear();
	for (n=1; n<N; n++) 
		if (DataPoints[n].x - DataPoints[n-1].x > 6.0 
		  || fabs(DataPoints[n].y - DataPoints[n-1].y) 
				> 4.0 * MAX(DataPoints[n].sig, DataPoints[n-1].sig))
			jumpList.append(n);
	splines = new Spline*[NumSplines = jumpList.size() + 1];
	for (n=0; n < NumSplines; n++) 
	{
		m = (k = (n == NumSplines - 1 ? N : jumpList[n])) 
		  - (j = (n==0 ? 0 : jumpList[n-1]));
		if (m==1)
		{
			sps = new SplinePoint[Nsps = 1];
			sps[0].x = DataPoints[j].x;
		}
		else if (m < 15)
		{
			sps = new SplinePoint[Nsps = 2];
			sps[0].x = DataPoints[j].x;
			sps[1].x = DataPoints[k-1].x;
		}
		else
		{
			sps = new SplinePoint[Nsps = (m + 15) / 10];
			step = double(m) / (Nsps - 1);
			for (i=0, cp = j; i < Nsps; cp += step, i++) 
				sps[i].x = DataPoints[(i < Nsps - 1 ? int(cp + 1e-5) : k-1)].x;
		}
		splines[n] = new Spline(Nsps, sps, false);
		splines[n]->FitTo(m, DataPoints + j);
	}
	delete[] DataPoints;
}

void ResidualFit::JoinSplines(int i_indexFirstSpline, const std::vector<Js> &i_data)
{
    if (-1 == i_indexFirstSpline)
    {
        for (i_indexFirstSpline = 0; i_indexFirstSpline < NumSplines && splines[i_indexFirstSpline]->getxN() < i_data[0].J; ++i_indexFirstSpline) ;
        if (splines[i_indexFirstSpline]->getxN() < i_data[0].J || i_indexFirstSpline + 1 >= NumSplines || splines[i_indexFirstSpline + 1]->getx0() > i_data.rbegin()->J) return;
    }
    splines[i_indexFirstSpline]->JoinForDeperturbation(splines[i_indexFirstSpline + 1], i_data[0].J, i_data.rbegin()->J);
    for (int j = i_indexFirstSpline + 2; j < NumSplines; ++j) splines[j-1] = splines[j];
    --NumSplines;
    ConnectSplines();
}

int ResidualFit::addLocalPerturbation(const double i_IsoF, const double i_Omega, const std::vector<Js> &i_data, const int i_JStep,
                    const double i_BeCurState)
{
    if (i_data.size() == 0) return -1;
    int i, j, NData;
    for (i=0; i < NumSplines && splines[i]->getxN() < i_data[0].J; ++i) ;
    if (i == NumSplines) return -1;
    Point *backupData/*, *splineFitData*/;
    if (splines[i]->getxN() >= i_data[0].J && i+1 < NumSplines && splines[i+1]->getx0() <= i_data.rbegin()->J)
    {
        int NData1, NData2;
        Point *tempFitData1, *tempFitData2;
        splines[i]->getFitData(NData1, tempFitData1);
        splines[i+1]->getFitData(NData2, tempFitData2);
        backupData = new Point[NData = NData1 + NData2];
        for (j=0; j < NData1; ++j) backupData[j] = tempFitData1[j];
        for (j=0; j < NData2; ++j) backupData[NData1 + j] = tempFitData2[j];
        JoinSplines(i, i_data);
    }
    else
    {
        Point *tempFitData;
        splines[i]->getFitData(NData, tempFitData);
        for (j=0, backupData = new Point[NData]; j < NData; ++j) backupData[j] = tempFitData[j];
    }
    std::vector<Js> pertFitData(i_data);
    for (std::vector<Js>::iterator it = pertFitData.begin(); it != pertFitData.end(); ++it) it->E_calc -= splines[i]->gety(it->J);
    LocalPerturbation* perturbation = new LocalPerturbation(i_IsoF, i_Omega, pertFitData, i_JStep, i_BeCurState);

    //perturbation->InitDebugLogging("LocalPerturbatioonFitDebugLogFile.txt");

    AutoFitLocalPerturbation(*perturbation, 1);
    /*for (j=0, splineFitData = new Point[NData]; j < NData; ++j)
    {
        splineFitData[j] = backupData[j];
        splineFitData[j].y -= perturbation->GetEnergyDiff(splineFitData[j].x);
    }
    splines[i]->FitTo(NData, splineFitData);
    for (j=0; j < static_cast<int>(i_data.size()); ++j) pertFitData[j].E_calc = i_data[j].E_calc - splines[i]->gety(i_data[j].J);
    perturbation->SetData(pertFitData);
    AutoFitLocalPerturbation(*perturbation, 1);
    for (j=0; j < NData; ++j) splineFitData[j].y = backupData[j].y - perturbation->GetEnergyDiff(splineFitData[j].x);
    splines[i]->FitTo(NData, splineFitData);
    for (j=0; j < static_cast<int>(i_data.size()); ++j) pertFitData[j].E_calc = i_data[j].E_calc - splines[i]->gety(i_data[j].J);
    perturbation->SetData(pertFitData);
    AutoFitLocalPerturbation(*perturbation, 10);
    DetermineUncertaintiesByBeVariation(*perturbation, *(splines[i]), NData, backupData, splineFitData, i_data,
                                        pertFitData);*/
    double centerJ = perturbation->GetCenter();
    std::vector<LocalPerturbation*>::iterator it = m_localPerturbations.begin();
    for (j=0; it != m_localPerturbations.end() && (*it)->GetCenter() < centerJ; ++j, ++it) ;
    m_localPerturbations.insert(it, perturbation);
    splines[i]->addLocalPerturbation(perturbation);

    //perturbation->EndDebugLogging();

    connect(perturbation, SIGNAL(Changed()), this, SIGNAL(Changed()));
    emit Changed();
    delete[] backupData;
    if (m_connectSplines.size() > 0) ConnectSplines();
    return j;
}

void ResidualFit::DetermineUncertaintiesByBeVariation(LocalPerturbation &io_perturbation)
{
    int i=0;
    while (!splines[i]->containsLocalPerturbation(&io_perturbation)) ++i;
    int NData;
    Point* splineFitData;
    splines[i]->getFitData(NData, splineFitData);
    Point*  backupSplineFitData = new Point[NData];
    for (int n=0; n < NData; ++n) backupSplineFitData[n] = splineFitData[n];
    std::vector<Js> perturbationFitData = io_perturbation.GetFitData();
    std::vector<Js> data = perturbationFitData;
    for (std::vector<Js>::iterator it = data.begin(); it != data.end(); ++it) it->E_calc += splines[i]->gety(it->J);
    DetermineUncertaintiesByBeVariation(io_perturbation, *(splines[i]), NData, backupSplineFitData, splineFitData, data, perturbationFitData);
    delete[] backupSplineFitData;
}

void ResidualFit::DetermineUncertaintiesByBeVariation(LocalPerturbation& io_perturbation,
     const Spline& i_spline, const int i_nData, Point *i_backupData, Point *i_splineFitData,
     const std::vector<Js> &i_data, std::vector<Js> i_pertFitData) const
{
    int NBandC = io_perturbation.GetNBandC();
    double H12min, H12max, *minBandC = new double[NBandC], *maxBandC = new double[NBandC];
    io_perturbation.SetBeFixed(true);
    bool result = (DetermineExtremalBeOneDirection(H12min, minBandC, io_perturbation, i_spline, -1.0, i_nData, i_backupData,
                                    i_splineFitData, i_data, i_pertFitData)
                && DetermineExtremalBeOneDirection(H12max, maxBandC, io_perturbation, i_spline, 1.0, i_nData, i_backupData,
                                    i_splineFitData, i_data, i_pertFitData));
    io_perturbation.SetBeFixed(false);
    if (result)
    {
        io_perturbation.SetH12Uncertainty(0.5 * abs(H12max - H12min));
        for (int n=0; n < NBandC; ++n) io_perturbation.SetBandCUncertainty(n, 0.5 * abs(maxBandC[n] - minBandC[n]));
    }
    else QMessageBox::warning(0, "MolSpektAnalysis", "The calculation of the uncertainties by Be variation failed!", QMessageBox::Ok);
    delete[] minBandC;
    delete[] maxBandC;
}

bool ResidualFit::DetermineExtremalBeOneDirection(double& o_H12Val, double* o_BandCVal,
    LocalPerturbation i_perturbation, Spline i_spline, const double i_dir, const int i_nData,
    Point *i_backupData, Point *i_splineFitData, const std::vector<Js> &i_data,
    std::vector<Js> i_pertFitData) const
{
    double currentSigma = i_perturbation.GetFitSigma(), expectedSigma = 1.35 * currentSigma;
    double startBe, B;
    i_perturbation.GetBandC(1, startBe, B);
    double initStep = i_dir * 0.01 * startBe, lastStartSigma = currentSigma, penultimateStartSigma = currentSigma;
    double currentBe = startBe, lastBe = startBe - initStep, penultimateBe = 0.0;
    const int maxIt = 100;
    int count, n, m;
    bool maxReached = false;
    for (count = 0; !maxReached && count < maxIt && (i_dir * (currentBe - lastBe) > 0.0 || i_dir * (lastBe - penultimateBe) > 0.0) &&
         (lastStartSigma < penultimateStartSigma || lastStartSigma < expectedSigma || penultimateStartSigma < expectedSigma); ++count)
    {
        for (int j=0; j < i_nData; ++j) i_splineFitData[j].y =
            i_backupData[j].y - i_perturbation.GetEnergyDiff(i_splineFitData[j].x);
        i_spline.FitTo(i_nData, i_splineFitData);
        for (int j=0; j < static_cast<int>(i_data.size()); ++j)
        {
            i_pertFitData[j].E_calc = i_data[j].E_calc - i_spline.gety(i_data[j].J);
            i_pertFitData[j].isUp = i_perturbation.IsDataUp(j);
        }
        i_perturbation.SetData(i_pertFitData);
        double minDUp = 1e99, maxDDown = -1e99;
        DetermineMinDUpMaxDDown(minDUp, maxDDown, n, m, i_perturbation);
        if (minDUp != 1e99 && maxDDown != -1e99 && minDUp - maxDDown <= 2e-3)
        {
            if (count > 0) break;
            else return false;
        }
        penultimateStartSigma = lastStartSigma;
        if (!DetermineExtremalBeOneStep(o_H12Val, o_BandCVal, maxReached, lastStartSigma, i_perturbation, initStep, startBe)) return false;
        currentSigma = i_perturbation.GetFitSigma();
        penultimateBe = lastBe;
        lastBe = currentBe;
        i_perturbation.GetBandC(1, currentBe, B);
    }
    return (count < maxIt);
}

void ResidualFit::clearCycleDetectList()
{
	int n, N = cycleDetectList->size();
	for (n=0; n<N; n++) delete[] (*cycleDetectList)[n];
	cycleDetectList->clear();
	cycleDetected = false;
}

bool ResidualFit::inCycle(int N, SplinePoint* spoints)
{
	int n, l, L = cycleDetectList->size();
	for (l=0; !cycleDetected && l<L; l++) for (cycleDetected = true, n=0; cycleDetected && n<N; n++) if (spoints[n] != (*cycleDetectList)[l][n]) 
		cycleDetected = false;
	if (!cycleDetected)
	{
		SplinePoint *deepCopy = new SplinePoint[N];
		for (n=0; n<N; n++) deepCopy[n] = spoints[n];
		cycleDetectList->append(deepCopy);
	}
	else cycleIndex = l;
	return cycleDetected;
}

bool ResidualFit::readData(QTextStream *i_stream)
{
    bool result = true, endFound = false;
    int numReadSplines = 0, SplineN = -1;
    if (splines != 0) clear();
    if (stateName != 0)
    {
        delete stateName;
        stateName = 0;
    }
    while (!i_stream->atEnd())
    {
        QString Buffer = i_stream->readLine();
        if (Buffer.indexOf("Electronic state:", 0, Qt::CaseInsensitive) >= 0)
        {
            if (stateName != 0)
            {
                result = false;
                break;
            }
            stateName = new QString(Buffer.right(Buffer.length() - 17).trimmed());
        }
        else if (Buffer.indexOf("Iso:", 0, Qt::CaseInsensitive) >= 0) iso = Buffer.right(Buffer.length() - 4).toInt();
        else if (Buffer.indexOf("v:") >= 0) v = Buffer.right(Buffer.length() - 2).toInt();
        else if (Buffer.indexOf("Comp:", 0, Qt::CaseInsensitive) >= 0) comp = Buffer.right(Buffer.length() - 5).toInt();
        else if (Buffer.indexOf("Number of splines:", 0, Qt::CaseInsensitive) >= 0)
        {
            NumSplines = Buffer.right(Buffer.length() - 18).toInt();
            if (splines != 0 || NumSplines <= 0)
            {
                result = false;
                break;
            }
            splines = new Spline*[NumSplines];
        }
        else if (Buffer.indexOf("Begin Spline", 0, Qt::CaseInsensitive) >= 0)
        {
            if (numReadSplines == NumSplines)
            {
                result = false;
                break;
            }
            splines[numReadSplines] = new Spline;
            if (!splines[numReadSplines]->readData(i_stream))
            {
                if (splines[numReadSplines]->GetError() == Spline::SEFileCorrupted) result = false;
                delete splines[numReadSplines];
                --NumSplines;
            }
            else ++numReadSplines;
            if (!result) break;
        }
        else if (Buffer.indexOf("Number of LocalPerturbations:", 0, Qt::CaseInsensitive) >= 0)
        {
            int N = Buffer.right(Buffer.length() - 29).toInt();
            m_localPerturbations.reserve(N);
        }
        else if (Buffer.contains("Spline:", Qt::CaseInsensitive)) SplineN = Buffer.right(Buffer.length() - 7).toInt();
        else if (Buffer.contains("Begin LocalPerturbation", Qt::CaseInsensitive))
        {
            if (m_localPerturbations.size() == m_localPerturbations.capacity())
            {
                result = false;
                break;
            }
            LocalPerturbation* pert = new LocalPerturbation;
            if (!pert->readData(i_stream))
            {
                result = false;
                break;
            }
            m_localPerturbations.push_back(pert);
            if (numReadSplines > SplineN) splines[SplineN]->addLocalPerturbation(pert);
        }
        else if (Buffer.indexOf("End ResidualFit", 0, Qt::CaseInsensitive) >= 0)
        {
            endFound = true;
            break;
        }
    }
    if (!endFound) result = false;
    if (!result)
    {
        for (int n=0; n < numReadSplines; ++n) delete splines[n];
        delete[] splines;
        splines = 0;
        NumSplines = 0;
    }
    return result;
}

void ResidualFit::RemoveLocalPerturbation(const int i_index)
{
    for (int i=0; i < NumSplines; ++i) splines[i]->RemoveLocalPerturbation(m_localPerturbations[i_index]);
    delete m_localPerturbations[i_index];
    m_localPerturbations.erase(m_localPerturbations.begin() + i_index);
}

void ResidualFit::setAssignment(const QString* i_stateName, const int i_iso, const int i_v, const int i_comp)
{
    iso = i_iso;
    v = i_v;
    comp = i_comp;
    if (stateName != 0) delete stateName;
    stateName = new QString(*i_stateName);
}

void ResidualFit::setAverage(int N, SplinePoint* spoints)
{
	int l, L = cycleDetectList->size(), cycleLength = L - cycleIndex, n;
	if (!cycleDetected) return;
	for (n=0; n<N; n++)
	{
		for (l = cycleIndex + 1; l<L; l++) spoints[n] += (*cycleDetectList)[l][n];
		spoints[n] /= cycleLength;
	}
}

void ResidualFit::setFitData(const int N, const Point * const points, const double i_IsoF, const double i_Omega, const int i_JStep, const double i_BeCurState)
{
    int n, s, l;
    for (n=s=l=0; s < NumSplines && n<N; ++n) if (points[n].x > splines[s]->getxN())
    {
        splines[s++]->copyData(n-l, points + l);
        l=n;
    }
    if (s < NumSplines) splines[s]->copyData(N-l, points + l);
    n=0;
    for (std::vector<LocalPerturbation*>::iterator it = m_localPerturbations.begin(); it != m_localPerturbations.end(); ++it)
    {
        while (n>0 && points[n].x > (*it)->GetMinJ() + 0.1) --n;
        while (n<N-1 && points[n].x < (*it)->GetMinJ() - 0.1) ++n;
        std::vector<Js> pertFitData;
        while (n<N && points[n].x < (*it)->GetMaxJ() + 0.1)
        {
            pertFitData.push_back(Js(points[n].x, points[n].obs - points[n].y, points[n].obs, points[n].unc));
            ++n;
        }
        (*it)->SetData(pertFitData);
        (*it)->SetElStateSpecificValues(i_IsoF, i_Omega, i_JStep, i_BeCurState);
    }
}

void ResidualFit::writeData(QTextStream *i_stream) const
{
    *i_stream << "Begin ResidualFit\n";
    *i_stream << "Electronic state: " << *stateName << '\n';
    *i_stream << "Iso: " << QString::number(iso) << '\n';
    *i_stream << "v: " << QString::number(v) << '\n';
    *i_stream << "Comp: " << QString::number(comp) << '\n';
    *i_stream << "Number of splines: " << QString::number(NumSplines) << '\n';
    for (int n=0; n < NumSplines; ++n) splines[n]->writeData(i_stream);
    *i_stream << "Number of LocalPerturbations: " << QString::number(m_localPerturbations.size()) << '\n';
    for (std::vector<LocalPerturbation*>::const_iterator it = m_localPerturbations.begin(); it != m_localPerturbations.end(); ++it)
    {
        int n;
        for (n=0; n < NumSplines && !splines[n]->containsLocalPerturbation(*it); ++n) ;
        *i_stream << "Spline: " << QString::number(n) << '\n';
        (*it)->writeData(i_stream);
    }
    *i_stream << "End ResidualFit\n";
}

bool ResidualFit::DetermineExtremalBeOneStep(double &o_H12Val, double * const o_BandCValues, bool& o_maxReached, double& o_startSigma,
     LocalPerturbation& i_localPert, double i_step, const double i_startBe) const
{
    LocalPerturbation fitPert = i_localPert;
    double currentSigma = i_localPert.GetFitSigma(), plateauSigma = 0.0, expectedSigma = 1.35 * currentSigma;
    double lastSigma = currentSigma, acceptedSigmaDev = 0.008 * expectedSigma, Be, lastBe = 0.0, B, T, iscDiff = 0.01 * acceptedSigmaDev;
    i_localPert.GetBandC(1, Be, B);
    const int maxIt = 100, maxIsc = 5;
    int count, nmUp, nmDown, iscCount = 0;
    bool stepWasTooLarge = false, notReachedBarrier = true;
    bool ilocalPertChanged = false, wasReverseDirection = false;
    o_maxReached = false;
    o_startSigma = currentSigma;
    for (count = 0; count < maxIt && abs(currentSigma - expectedSigma) >= acceptedSigmaDev; ++count)
    {
        double minDUp = 1e99, maxDDown = -1e99;
        if (!stepWasTooLarge && fabs(lastSigma - currentSigma) < iscDiff && (Be - lastBe) * i_step < 0.0) ++iscCount;
        else iscCount = 0;
        if (iscCount == maxIsc)
        {
            if (!ilocalPertChanged) expectedSigma *= currentSigma / o_startSigma;
            plateauSigma = currentSigma;
        }
        lastBe = Be;
        if (lastSigma == currentSigma && !stepWasTooLarge) Be += i_step;
        else
        {
            if (stepWasTooLarge || currentSigma > expectedSigma)
            {
                if (fabs(i_step) < 0.05 * fabs(Be - i_startBe)) break;
                Be -= (i_step *= 0.5);
                fitPert = i_localPert;
                if (plateauSigma > 0.0) notReachedBarrier = false;
                wasReverseDirection = true;
            }
            else
            {
                if (fabs(currentSigma - plateauSigma) < iscDiff)
                {
                    if (notReachedBarrier) Be += (i_step *= 2.0);
                    else if (fabs(i_step) < 0.05 * fabs(Be - i_startBe)) break;
                    else Be += (i_step *= 0.5);
                }
                else if (wasReverseDirection) Be += (i_step *= 0.5);
                else if (expectedSigma - currentSigma >= 2.0 * (currentSigma - lastSigma)) Be += (i_step *= 2.0);
                else if (expectedSigma - currentSigma >= currentSigma - lastSigma) Be += i_step;
                else Be += (i_step *= 0.5);
                i_localPert = fitPert;
                ilocalPertChanged = true;
                wasReverseDirection = false;
            }
        }
        fitPert.SetBandC(1, Be, 0.0);
        DetermineMinDUpMaxDDown(minDUp, maxDDown, nmUp, nmDown, fitPert);
        if (minDUp != 1e99 && maxDDown != -1e99)
        {
            if (minDUp - maxDDown > 2e-3)
            {
                fitPert.GetBandC(0, T, B);
                T -= 0.5 * (minDUp + maxDDown);
                fitPert.SetBandC(0, T, B);
                stepWasTooLarge = false;
            }
            else
            {
                o_maxReached = true;
                double R1 = 0.0, R2 = 0.0, IsoF = fitPert.GetIsoF(), O = fitPert.GetOmega(), J1 = fitPert.GetDataJValue(nmUp);
                double J2 = fitPert.GetDataJValue(nmDown);
                double JF1 = IsoF * (J1 * (J1 + 1.0) - O * O), JF2 = IsoF * (J2 * (J2 + 1.0) - O * O);
                for (int n = fitPert.GetNBandC() - 1; n >= 2; --n)
                {
                    fitPert.GetBandC(n, T, B);
                    R1 += T;
                    R2 += T;
                    R1 *= JF1;
                    R2 *= JF2;
                }
                double newBe = (minDUp - maxDDown + 2e-3 - JF1 * R1 + JF2 * R2) / (IsoF * (J1 * (J1 + 1.0) - J2 * (J2 + 1.0)));
                double avBe = 0.5 * (lastBe + Be);
                if (abs(avBe - newBe) < abs(Be - avBe))
                {
                    stepWasTooLarge = false;
                    Be = newBe;
                }
                else
                {
                    stepWasTooLarge = true;
                    continue;
                }
                T = minDUp + 1e-3 - JF1 * (Be + R1);
                fitPert.SetBandC(0, T, 0.0);
                fitPert.SetBandC(1, Be, 0.0);
            }
        }
        else if (minDUp < 0.0 || maxDDown > 0.0)
        {
            fitPert.GetBandC(0, T, B);
            T += (minDUp < 0.0 ? minDUp - 1e-3 : maxDDown + 1e-3);
            fitPert.SetBandC(0, T, B);
        }
        //FitLocalPerturbation(i_localPert, Be);
        fitPert.LevenbergMarquardt(100, 1e-6);
        lastSigma = currentSigma;
        currentSigma = fitPert.GetFitSigma();
        if (o_maxReached)
        {
            if (currentSigma < expectedSigma + acceptedSigmaDev) break;
            else
            {
                o_maxReached = false;
                i_step = Be - lastBe;
            }
        }
    }
    fitPert.GetH12(o_H12Val, B);
    for (int n=0; n < fitPert.GetNBandC(); ++n) fitPert.GetBandC(n, o_BandCValues[n], B);
    return (count < maxIt);
}

void ResidualFit::DetermineMinDUpMaxDDown(double &o_minDUp, double &o_maxDDown, int &o_nmUp, int& o_nmDown, const LocalPerturbation &i_localPert) const
{
    for (int n=0; n < i_localPert.GetNData(); ++n)
    {
        bool isUp = i_localPert.IsDataUp(n);
        double EDiff = i_localPert.GetDataEobs(n) - i_localPert.GetPointCalcBand(i_localPert.GetDataJValue(n));
        if (isUp && EDiff < o_minDUp)
        {
            o_minDUp = EDiff;
            o_nmUp = n;
        }
        else if (!isUp && EDiff > o_maxDDown)
        {
            o_maxDDown = EDiff;
            o_nmDown = n;
        }
    }
}

double ResidualFit::GetPointPerturber(const double JValue, const int pertIndex)
{
    return m_localPerturbations[pertIndex]->GetPointCalcBand(JValue) - m_localPerturbations[pertIndex]->GetPointUpOBCa(JValue)
            + GetPointFromSplines(JValue);
}

double ResidualFit::GetPoint(const double JValue)
{
    if (m_localPerturbations.size() == 0) return GetPointFromSplines(JValue);
    Molecule* mol = (m_state != 0 ? m_state->getMolecule() : 0);
    int n, JStep = (mol != 0 ? mol->getJStep(iso) : 1);
    for (n=0; n < static_cast<int>(m_localPerturbations.size()) && m_localPerturbations[n]->GetCenter() < JValue - JStep; ++n) ;
    if (n < static_cast<int>(m_localPerturbations.size()) && abs(m_localPerturbations[n]->GetCenter() - JValue) <= JStep)
    {
        double R1 = GetPoint(JValue, n), R2 = GetPoint(JValue, n+1);
        return (fabs(R1) < fabs(R2) ? R1 : R2);
    }
    return GetPoint(JValue, n);
}

double ResidualFit::GetPoint(const double i_JValue, const int i_pertIndex)
{
    LocalPerturbation* perturbation, *nextPerturbation;
    double Jval, ref;
    bool nextPertUp, pertUp;
    if (i_pertIndex == 0)
    {
        perturbation = m_localPerturbations[0];
        Jval = perturbation->GetDataJValue(0);
        ref = perturbation->GetPointUpOBCa(Jval);
        pertUp = (abs(ref - perturbation->GetPointUpPert(Jval)) < abs(ref - perturbation->GetPointLowPert(Jval)));
        return GetPointFromSplines(i_JValue) + perturbation->GetPointPert(pertUp, i_JValue) - perturbation->GetPointUpOBCa(i_JValue);
    }
    if (i_pertIndex < static_cast<int>(m_localPerturbations.size()))
    {
        perturbation = m_localPerturbations[i_pertIndex - 1];
        Jval = perturbation->GetDataJValue(0);
        ref = perturbation->GetPointUpOBCa(Jval);
        pertUp = (abs(ref - perturbation->GetPointUpPert(Jval)) < abs(ref - perturbation->GetPointLowPert(Jval)));
        nextPerturbation = m_localPerturbations[i_pertIndex];
        Jval = nextPerturbation->GetDataJValue(0);
        ref = nextPerturbation->GetPointUpOBCa(Jval);
        nextPertUp = (abs(ref - nextPerturbation->GetPointUpPert(Jval)) < abs(ref - nextPerturbation->GetPointLowPert(Jval)));
        return GetPointFromSplines(i_JValue) + perturbation->GetPointPert(!pertUp, i_JValue) - perturbation->GetPointUpOBCa(i_JValue) +
                             nextPerturbation->GetPointPert(nextPertUp, i_JValue) - nextPerturbation->GetPointUpOBCa(i_JValue);
    }
    nextPerturbation = m_localPerturbations[i_pertIndex - 1];
    Jval = nextPerturbation->GetDataJValue(0);
    ref = nextPerturbation->GetPointUpOBCa(Jval);
    nextPertUp = (abs(ref - nextPerturbation->GetPointUpPert(Jval)) < abs(ref - nextPerturbation->GetPointLowPert(Jval)));
    return GetPointFromSplines(i_JValue) + nextPerturbation->GetPointPert(!nextPertUp, i_JValue) - nextPerturbation->GetPointUpOBCa(i_JValue);
}

void ResidualFit::ConnectSplines()
{
    if (m_connectSplines.size() > 0)
    {
        for (std::vector<ConnectSpline*>::const_iterator it = m_connectSplines.begin(); it != m_connectSplines.end(); ++it) delete *it;
        m_connectSplines.clear();
    }
    if (NumSplines == 0) return;
    Molecule* mol = (m_state != 0 ? m_state->getMolecule() : 0);
    double JStep = (mol != 0 ? mol->getJStep(iso) : 1.0) + 0.1;
    if (splines[0]->getx0() >= (m_state != 0 ? m_state->getJStart(iso, comp) + 0.1 : 0.5))
        m_connectSplines.push_back(new ConnectSpline(0, splines[0]));
    for (int i=1; i < NumSplines; ++i) if (splines[i]->getx0() - splines[i-1]->getxN() > JStep) m_connectSplines.push_back(new ConnectSpline(splines[i-1], splines[i]));
    m_connectSplines.push_back(new ConnectSpline(splines[NumSplines - 1], 0));
}

double ResidualFit::GetPointFromSplines(const double JValue)
{
    for (int i=0; i < NumSplines && splines[i]->getx0() <= JValue; ++i) if (splines[i]->getxN() >= JValue) return splines[i]->gety(JValue);
    if (m_connectSplines.size() == 0) ConnectSplines();
    for (std::vector<ConnectSpline*>::const_iterator it = m_connectSplines.begin(); it != m_connectSplines.end() && (*it)->GetJMin() <= JValue; ++it)
        if ((*it)->GetJMax() >= JValue) return (*it)->GetValue(JValue);
    return 1e99;
}
