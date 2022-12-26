#ifndef LOGWINDOW_H
#define LOGWINDOW_H


#include "tablewindow.h"
#include "logmodel.h"

class QFile;
class QTextStream;
class QLineEdit;
class QPushButton;
class QCheckBox;
class QMutex;

class MainWindow;


class LogWindow : public TableWindow
{
    Q_OBJECT
public:
    explicit LogWindow(MainWindow *parent = 0);
    ~LogWindow();

    void LogMessage(QStringList& message);
    void SetMessageBuffer(QMutex& loggerMutex, LogList &buffer);

private slots:
    void ShowFileDialog();
    void WriteLogFileChanged(int state);
    void MaxRowsCanged();

private:
    QLineEdit *mFilenameEdit, *mMaxTableSize;
    QPushButton* mFileDialogButton;
    QCheckBox* mWriteLogFileCheckBox;

    LogModel mModel;

    QMutex* mLoggerMutex;
    QTextStream *mLogStream;
    QFile *mLogFile;
    bool mWriteToFile;
};

#endif // LOGWINDOW_H
