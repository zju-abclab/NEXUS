#pragma once

#include <iostream>
#include <list>
#include <NTL/RR.h>
#include <cmath>
#include <vector>
#include "PolyUpdate.h"
#include "MinicompFunc.h"
#include "RemezApp.h"
#include "seal/seal.h"

using namespace std;
using namespace NTL;
using namespace minicomp;

void upgrade_oddbaby(long n, Tree& tree);
void upgrade_baby(long n, Tree& tree);