#pragma once


#include <QStringList>


struct LogListItem
{
    LogListItem() : next(nullptr), prev(nullptr) {}

    LogListItem* next, *prev;
    QStringList message;
};
