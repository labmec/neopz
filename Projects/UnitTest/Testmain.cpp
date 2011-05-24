//
//  TestMain.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//



#ifdef USING_UNIT_TEST

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN my application tests

#include <boost/test/unit_test.hpp>

#else

#include <iostream>

int main()
{
    std::cout << "Boost needs to be configured for testing to work\n";
    return 0;
}

#endif