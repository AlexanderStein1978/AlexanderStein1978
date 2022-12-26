#include "logmodel.h"
#include "loglist.h"


LogModel::LogModel(QObject *parent) : mMessageBuffer(nullptr)
{
}

LogModel::~LogModel()
{
    if (nullptr != mMessageBuffer) mMessageBuffer->clear();
}


int LogModel::columnCount(const QModelIndex &parent) const
{
    return (parent.isValid() ? 0 : 5);
}

int LogModel::rowCount(const QModelIndex &parent) const
{
    return (parent.isValid() || nullptr == mMessageBuffer ? 0 : mMessageBuffer->size());
}

QVariant LogModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || nullptr == mMessageBuffer || mMessageBuffer->size() == 0) return QVariant();
    if (role == Qt::TextAlignmentRole) return Qt::AlignLeft;
    if (role != Qt::DisplayRole) return QVariant();
    return mMessageBuffer->getElement(index.row())[index.column()];
}

QVariant LogModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role != Qt::DisplayRole) return QVariant();
    if (orientation == Qt::Vertical) return section;
    switch (section)
    {
        case 0:
            return "Date and time";
        case 1:
            return "Severity";
        case 2:
            return "Function";
        case 3:
            return "File and line";
        default:
            return "Message";
            break;
    }
}

void LogModel::AddLogMessage(const QStringList& message)
{
    bool maxReached = mMessageBuffer->isMaxReached();
    beginInsertRows(QModelIndex(), 0, 0);
    if (maxReached) beginRemoveRows(QModelIndex(), mMessageBuffer->size() - 1, mMessageBuffer->size() - 1);
    mMessageBuffer->addElement(message);
    if (maxReached) endRemoveRows();
    endInsertRows();
}

void LogModel::SetMessageBuffer(LogList& messageBuffer)
{
    beginInsertRows(QModelIndex(), 0, messageBuffer.size() - 1);
    mMessageBuffer = &messageBuffer;
    endInsertRows();
}

void LogModel::SetMaxRows(const int maxValue)
{
    if (mMessageBuffer->size() > maxValue)
    {
        beginRemoveRows(QModelIndex(), maxValue, mMessageBuffer->size());
        mMessageBuffer->SetMaxSize(maxValue);
        endRemoveRows();
    }
    else mMessageBuffer->SetMaxSize(maxValue);
}

