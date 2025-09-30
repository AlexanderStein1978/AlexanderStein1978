//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include <QWidget>


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
    QComboBox *mDeviceBox, *mSampleSizeBox, *mSampleRateBox, *mSampleTypeBox;
    SoundRecordAndDrawControl* mControl;
    const char *const mInputData;
    const int mNBytes;
};
