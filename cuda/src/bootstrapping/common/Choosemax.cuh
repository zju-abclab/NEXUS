#pragma once

#include <NTL/RR.h>

#include <iostream>

using namespace std;
using namespace NTL;

int MinIndex(RR* c, int num);
void MaxSubsetSum(RR* a, int m, int n, int* cur_index);
