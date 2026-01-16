#pragma once
#include "Structures.h"

void computePForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc, double alpha, double Tot);

void computePForAll(Grid& grid, const Grid& g, const ElemUniv& univ, int npc, double alpha, double Tot);

void addLocalToPG(Equation& eq, const int nodeIds[4], const double Ploc[4]);