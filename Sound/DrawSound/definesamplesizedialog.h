//
// Author: Alexander Stein <AlexanderStein@t-online.de>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include <QWidget>

<<<<<<< HEAD
=======

>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
class SoundRecordAndDrawControl;
class QComboBox;


class DefineSampleSizeDialog : public QWidget
{
Q_OBJECT

public:
    DefineSampleSizeDialog(SoundRecordAndDrawControl* parent, const char *const inputData, const int nBytes);
    ~DefineSampleSizeDialog();

private slots:
    void Draw();
    void Save();
    void DeviceChanged(int index);

private:
<<<<<<< HEAD
    QComboBox *mDeviceBox, *mSampleSizeBox, *mSampleRateBox;
=======
    QComboBox *mDeviceBox, *mSampleSizeBox, *mSampleRateBox, *mSampleTypeBox;
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
    SoundRecordAndDrawControl* mControl;
    const char *const mInputData;
    const int mNBytes;
};
