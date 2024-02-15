#include "windowselectdialog.h"
#include "DiagWindow.h"

#include <QGridLayout>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>


WindowSelectDialog::WindowSelectDialog(const std::vector<DiagWindow*>& windows, QWidget *parent) : QDialog(parent), mSelectionBox(new QComboBox(this)), mOkButton(new QPushButton("OK", this)),
    mCancelButton(new QPushButton("Cancel", this))
{
    setWindowTitle("Please select window to add data");
    QGridLayout* layout = new QGridLayout(this);
    layout->addWidget(new QLabel("Window:", this), 0, 0);
    layout->addWidget(mSelectionBox, 0, 1);
    layout->setRowMinimumHeight(1, 20);
    layout->addWidget(mOkButton, 2, 0);
    layout->addWidget(mCancelButton, 2, 1);
    for (DiagWindow* w : windows) mSelectionBox->addItem(w->getName());
    mSelectionBox->addItem("New...");
    mSelectionBox->setEditable(false);
    connect(mOkButton, SIGNAL(clicked()), this, SLOT(accept()));
    connect(mCancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

int WindowSelectDialog::GetSelection()
{
    return mSelectionBox->currentIndex();
}


NameSelectionDialog::NameSelectionDialog(QWidget* parent) : QDialog(parent), mNameEdit(new QLineEdit("", this)), mOkButton(new QPushButton("OK", this)), mCancelButton(new QPushButton("Cancel", this))
{
    setWindowTitle("Please enter name");
    QGridLayout* layout = new QGridLayout(this);
    layout->addWidget(new QLabel("Name:", this), 0, 0);
    layout->addWidget(mNameEdit, 0, 1);
    layout->setRowMinimumHeight(1, 20);
    layout->addWidget(mOkButton, 2, 0);
    layout->addWidget(mCancelButton, 2, 1);
    connect(mOkButton, SIGNAL(clicked()), this, SLOT(accept()));
    connect(mCancelButton, SIGNAL(clicked()), this, SLOT(reject()));
}

QString NameSelectionDialog::GetName()
{
    return mNameEdit->text();
}
