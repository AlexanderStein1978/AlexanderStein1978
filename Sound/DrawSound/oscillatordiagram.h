//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once

#include "DiagWindow.h"
#include "oscillatorarray.h"

class SoundMainWindow;


class OscillatorDiagram : public DiagWindow
{
public:
    OscillatorDiagram(SoundMainWindow* MW, const QString& filename);
    ~OscillatorDiagram();

    void setData(OscillatorArray::Results& data);

    inline const OscillatorArray::Results& getData() const
    {
        return mData;
    }

protected:
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;

private:
    OscillatorArray::Results mData;
    int *mTimePixelAssignments, *mFrequencyPixelAssignments;
};
