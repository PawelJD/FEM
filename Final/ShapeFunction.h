#pragma once
#include "Structures.h"

array<double, 4> Nvals(double ksi, double eta);

void ShapeFunction(int npc, ElemUniv& univ);