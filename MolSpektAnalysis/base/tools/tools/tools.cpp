//
// C++ Implementation: tools
//
//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2006 - 2019
//
// Copyright: See README file that comes with this source code
//
//

#include "tools.h"
#include "utils.h"
#include "marker.h"
#include "progressions.h"

#include <math.h>
#include <stdio.h>

#include <qmessagebox.h>
#include <QLineEdit>


bool SolvLinEqS(double **G, int n)
{
	int i, j, k;
	double *B, F, T;
	for (i=0; i<n-1; i++)
	{
		for (B=G[0]; B != NULL;) 
		{
			B = NULL;
			for (j=i+1; j<n; j++) if (fabs(G[j-1][i]) < fabs(G[j][i]))
			{
				B = G[j-1];
				G[j-1] = G[j];
				G[j] = B;
			}
		}
		for (j=n-1; G[j][i]==0.0 && j>i; j--);
		if (j==i) return false;
		for (j--; j>=i; j--)
		{
			F = G[j+1][i] / G[j][i];
			for (k=i+1; k<=n; k++)
			{
				T = F * G[j][k];
				G[j+1][k] -= T;
				//if (fabs(G[j+1][k]) < fabs(Err * T)) G[j+1][k] = 0.0;
			}
		}
	}
	if (G[n-1][n-1] == 0.0) return false;
	for (i=n-1; i>=0; i--) 
	{
		for (j=n-1; j>i; j--) G[i][n] -= G[j][n] * G[i][j];
		G[i][n] /= G[i][i];
	}
	return true;
}	

bool SolvLinEqSwItImp(double** EqnS, int N, double *mFQS)
{
	double **tEqnS = Create(N, N+1), FQS = 1e98, aFQS = 1e99, tRes[N], bRes[N];
	int n, m;
	for (n=0; n<N; n++)
	{
		tEqnS[n][N] = EqnS[n][N];
		tRes[n] = 0.0;
	}
	while (FQS < aFQS)
	{
		aFQS = FQS;
		for (n=0; n<N; n++) for (m=0; m<N; m++) tEqnS[n][m] = EqnS[n][m];
		if (!SolvLinEqS(tEqnS, N)) return false;
		for (n=0; n<N; n++) tRes[n] += tEqnS[n][N];
		for (n=0, FQS = 0.0; n<N; n++)
		{
			for (tEqnS[n][N] = EqnS[n][N], m=0; m<N; m++) tEqnS[n][N] -= EqnS[n][m] * tRes[m];
			FQS += tEqnS[n][N] * tEqnS[n][N];
		}
		if (FQS < aFQS) for (n=0; n<N; n++) bRes[n] = tRes[n];
		printf("SolvLinEqSwItImp: FQS=%g\n", FQS);
	}
	for (n=0; n<N; n++) EqnS[n][N] = bRes[n];
	Destroy(tEqnS, N);
	if (mFQS != 0) *mFQS = aFQS;
	return true;
}

int Runden(const double &d)
{
    int i = (int)d;
    if (d - (double)i > 0.5) return i+1;
    return i;
}

void Test(double &d1, double &d2, QLineEdit *L1, QLineEdit *L2)
{
    double buffer;
    d1 = L1->text().toDouble();
    d2 = L2->text().toDouble();
    if (d1 >= d2)
    {
	if (d1 > d2)
	{
	    buffer = d1;
	    d1 = d2;
	    d2 = buffer;
	}
	else 
	{
	    d1 -= 0.5;
	    d2 += 0.5;
	}
	L1->setText(QString::number(d1, 'g', 11));
	L2->setText(QString::number(d2, 'g', 11));
    }
}

void GSearchLines(Marker *marker, Marker *LaserLine, const int &AnzahlMarker, const double &ST,
		  const double &GMH, Progressions &Prog, const int &veu, const int &Jeu, const int &NI,
		  double ***&ELU)
{
    Marker *FoundLines[veu+1];
    Marker **MB, *LL;
    int FoundCount = 0, i, j=0, FG=0, AGMarker = 0, Ii, viu, Jiu;
    for (i=0; i<AnzahlMarker; i++) 
    {
	if (marker[i].Marked) 
	{
	    j++;
	    FG = i;
	}
	if (marker[i].HFLM >= GMH) AGMarker++; 
    }
    if (j==1) LL = marker + FG;
    else LL = LaserLine;
    Marker *GMarker[AGMarker];
    j = 0;
    for (i=0; i<AnzahlMarker; i++) if (marker[i].HFLM >= GMH) GMarker[j++] = marker + i;
    double LES, ZS;
    for (Ii=0; Ii < NI; Ii++) 
    {
	for (Jiu=0; Jiu < Jeu; Jiu++) 
	{
	    for (viu=0; viu < veu; viu++)
	    {
		//printf("Ii=%d, Jiu=%d, viu=%d\n", Ii, Jiu, viu);
		if(ELU[Ii][viu][Jiu] == 0.0) break;
		LES = LL->Line[0] + ELU[Ii][viu][Jiu];
		j = AGMarker - 2;
		for (i=0; i<veu; i++)
		{
		    if(ELU[Ii][i][Jiu] == 0.0) break;
		    ZS = LES - ELU[Ii][i][Jiu];
		    while (ZS < GMarker[j]->Line[0] && j>0) j--;
		    if (GMarker[j+1]->Line[0] - ZS < ZS - GMarker[j]->Line[0]) j++;
		    if (fabs(GMarker[j]->Line[0] - ZS) < ST) FG++;
		    if (j==AGMarker-1) j--;
		} 
		if (FG > 2)
		{
		    printf("Progression gefunden.\n");
		    j = AnzahlMarker - 2;
		    for (i=0; i<veu; i++)
		    {
			if(ELU[Ii][i][Jiu] == 0.0) break;
			ZS = LES - ELU[Ii][i][Jiu];
			while (ZS < marker[j].Line[0] && j>0) j--;
			if (marker[j+1].Line[0] - ZS < ZS - marker[j].Line[0]) j++;
			if (fabs(marker[j].Line[0] - ZS) < ST)
			{
			    //printf("Gefunden: v=%d, J=%d, WN=%f\n", i, Jiu, marker[j].Line[0]);
			    FoundLines[FoundCount] = marker + j;
			    FoundLines[FoundCount]->DD = ZS - FoundLines[FoundCount]->Line[0];
			    FoundLines[FoundCount]->Iso = Ii;
			    FoundLines[FoundCount]->vss = i;
			    FoundLines[FoundCount++]->Jss = Jiu;
			}
			if (j==AnzahlMarker-1) j--;
		    } 
		    MB = new Marker*[FoundCount];
		    for (i=0; i<FoundCount; i++) MB[i] = FoundLines[i];
		    Prog.Insert(MB, FoundCount, FG);
		    printf("FoundCount=%d\n", FoundCount);
		}
		FoundCount = FG = 0;
	    }
	}
    }
}

void ParabInterpol(double x1, double y1, double x2, double y2, double x3, double y3, double &c, double &a)
{
	double b, xd21, xd32, x1q, x2q, x3q, dxq12, dxq32, dq, yd12, yd23, dxc;
	xd21 = x2 - x1;
	xd32 = x3 - x2;
	x1q = x1 * x1;
	x2q = x2 * x2;
	x3q = x3 * x3;
	dxq12 = x1q - x2q;
	dxq32 = x3q - x2q;
	dq = xd21 / xd32;
	yd12 = y1 - y2;
	yd23 = y2 - y3;
	b = (yd12 - yd23 * dq) / (dxq12 + dxq32 * dq);
	c = (yd23 + b * dxq32) / (2 * b * xd32);
	dxc = x1 - c;
	a = y1 - b * dxc * dxc;
}

