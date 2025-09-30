//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "gtest/gtest.h"
#include <qapplication.h>

int main(int argc, char **argv)
{
	QApplication a( argc, argv );
	
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
