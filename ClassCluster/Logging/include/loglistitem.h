//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#pragma once


#include <QStringList>


struct LogListItem
{
    LogListItem() : next(nullptr), prev(nullptr) {}

    LogListItem* next, *prev;
    QStringList message;
};
