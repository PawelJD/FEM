#pragma once
#include "Structures.h"

void computeHbcForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc, double alpha);

void computeHbcForAll(Grid& grid, const Grid& g, const ElemUniv& univ, int npc, double alpha);