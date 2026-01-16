#pragma once
#include "Structures.h"

void computeCForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc, double ro, double Cp);

void computeCForAll(Grid& grid, const ElemUniv& univ, int npc, double ro, double Cp);

void addLocalToCG(Equation& eq, const int nodeIds[4], const double Cloc[4][4]);