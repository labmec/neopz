//
//  TestMain.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

//#include "TestMain.h"

#ifdef USING_BOOST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN my application tests

#include <boost/test/unit_test.hpp>

#else

int main()
{
    std::cout << "Boost needs to be configured for testing to work\n";
    return 0;
}

#endif