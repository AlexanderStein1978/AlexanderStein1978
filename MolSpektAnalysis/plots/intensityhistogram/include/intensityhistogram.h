//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#ifndef INTENSITYHISTOGRAM_H
#define INTENSITYHISTOGRAM_H

#include <DiagWindow.h>

class Spektrum;

class QPainter;

class IntensityHistogram : public DiagWindow
{
  public:
	IntensityHistogram(Spektrum *S, MainWindow* MW = 0);
	
  private:
    void PSpektrum(QPainter &P, const QRect& A, bool PrintFN);
	
	Spektrum *Spectrum;
};

#endif // INTENSITYHISTOGRAM_H
