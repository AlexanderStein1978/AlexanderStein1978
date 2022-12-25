#include "logmodel.h"


LogModel::LogModel(QObject *parent) : mMaxRowCount(2147483647)
{

}

int LogModel::columnCount(const QModelIndex &parent) const
{
    return (parent.isValid() ? 0 : 5);
}

int LogModel::rowCount(const QModelIndex &parent) const
{
    return (parent.isValid() ? 0 : mMessageBuffer.size());
}

QVariant LogModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid() || mMessageBuffer.size() == 0) return QVariant();
    if (role == Qt::TextAlignmentRole) return Qt::AlignLeft;
    if (role != Qt::DisplayRole) return QVariant();
    return mMessageBuffer[index.row()][index.column()];
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

void LogModel::AddLogMessage(QStringList &message)
{
    if (mMessageBuffer.size() == mMaxRowCount)
    {
        beginRemoveRows(QModelIndex(), 0, 0);
        mMessageBuffer.pop_front();
        endInsertRows();
    }
    beginInsertRows(QModelIndex(), mMessageBuffer.size(), mMessageBuffer.size());
    mMessageBuffer.append(message);
    endInsertRows();
}

void LogModel::SetMessageBuffer(QList<QStringList> &messageBuffer)
{
    beginInsertRows(QModelIndex(), 0, messageBuffer.size() - 1);
    mMessageBuffer = messageBuffer;
    endInsertRows();
}

void LogModel::SetMaxRows(const int maxValue)
{
    if (mMessageBuffer.size() > maxValue)
    {
        mMaxRowCount = maxValue;
        int mR = mMessageBuffer.size() - maxValue - 1;
        beginRemoveRows(QModelIndex(), 0, mR);
        for (int n = mR; n>=0; n--) mMessageBuffer.removeAt(n);
        endRemoveRows();
    }
}

