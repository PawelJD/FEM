#pragma once
#include <iostream>

#include "Structures.h"

Jakobian computeJacobian(const Element& elem, const Grid& grid, const array<double, 4>& dN_dKsi, const array<double, 4>& dN_dEta);

void computeJacobiansForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc);

void computeJacobiansForAll(Grid& grid, const ElemUniv& univ, int npc);