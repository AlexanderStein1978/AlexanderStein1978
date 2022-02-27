#include "logwindow.h"

#include <QGridLayout>
#include <QCheckBox>


LogWindow::LogWindow(MainWindow *parent) : TableWindow(External, parent, nullptr), mFilenameEdit(new QLineEdit(this)), mFileDialogButton(new QPushButton("...", this)),
    mWriteLogFileCheckBox(new QCheckBox("Write messages to file", this)), mModel(), mLoggerMutex(nullptr), mLogStream(nullptr), mLogFile(nullptr), mWriteToFile(false)
{
    QGridLayout L(this);
    L.addWidget(mWriteLogFileCheckBox, 0, 0);
    L.addWidget(new QLabel("Filename:", this), 0, 1);
    L.addWidget(mFilenameEdit, 0, 2);
    L.addWidget(mFileDialogButton, 0, 3);
    L.setColumnMinimumWidth(3, 30);
    L.setColumnStretch(0, 10);
    L.setColumnStretch(1, 10);
    L.setColumnStretch(2, 10);
    L.addWidget(table = new MTable(this), 1, 0, 1, 4);
    L.setRowMinimumHeight(0, 20);
    L.setRowStretch(1, 10);
    table->setModel(&mModel);
    connect(mWriteLogFileCheckBox, SIGNAL(Qt::CheckState), this, SLOT(WriteLogFileChanged(Qt::CheckState)));
    connect(mFileDialogButton, SIGNAL(clicked()), this, SLOT(ShowFileDialog()));
}
