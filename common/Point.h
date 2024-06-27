#pragma once
#include<iostream>
#include<NTL/RR.h>

using namespace std;
using namespace NTL;

class Point {
public:
	RR x;
	RR y;

	long locmm;

	Point();
	Point(RR _x, RR _y);
};
