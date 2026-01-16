#pragma once
#include "Structures.h"

void computeHForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc, double k);

void computeHForAll(Grid& grid, const ElemUniv& univ, int npc, double k);

void addLocalToHG(Equation& eq, const int nodeIds[4], const double Hloc[4][4], const double Hbcloc[4][4]);