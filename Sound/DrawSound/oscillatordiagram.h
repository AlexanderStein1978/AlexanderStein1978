#pragma once

#include "DiagWindow.h"
#include "oscillatorarray.h"

class SoundMainWindow;


class OscillatorDiagram : public DiagWindow
{
public:
    OscillatorDiagram(SoundMainWindow* MW);
    ~OscillatorDiagram();

    void setData(OscillatorArray::Results& data);

protected:
    void PSpektrum(QPainter &P, const QRect &A, bool PrintFN ) override;

private:
    OscillatorArray::Results mData;
    int *mTimePixelAssignments, *mFrequencyPixelAssignments;
};
