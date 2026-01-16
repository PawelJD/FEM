#include <array>

#include "Jacobian.h"

using namespace std;

Jakobian computeJacobian(const Element& elem, const Grid& grid, const array<double, 4>& dN_dKsi, const array<double, 4>& dN_dEta) {

	Jakobian J;
	const Node& n1 = grid.nodes[elem.nodeIds[0] - 1];
	const Node& n2 = grid.nodes[elem.nodeIds[1] - 1];
	const Node& n3 = grid.nodes[elem.nodeIds[2] - 1];
	const Node& n4 = grid.nodes[elem.nodeIds[3] - 1];

	//cout << "Wezly elementu: "
 //   << "(" << n1.x << "," << n1.y << ") "
 //   << "(" << n2.x << "," << n2.y << ") "
 //   << "(" << n3.x << "," << n3.y << ") "
 //   << "(" << n4.x << "," << n4.y << ")\n";

	// Macierz J
	J.J[0][0] = dN_dKsi[0] * n1.x + dN_dKsi[1] * n2.x + dN_dKsi[2] * n3.x + dN_dKsi[3] * n4.x;
	J.J[0][1] = dN_dEta[0] * n1.x + dN_dEta[1] * n2.x + dN_dEta[2] * n3.x + dN_dEta[3] * n4.x;
	J.J[1][0] = dN_dKsi[0] * n1.y + dN_dKsi[1] * n2.y + dN_dKsi[2] * n3.y + dN_dKsi[3] * n4.y;
	J.J[1][1] = dN_dEta[0] * n1.y + dN_dEta[1] * n2.y + dN_dEta[2] * n3.y + dN_dEta[3] * n4.y;

	// Wyznacznik J
	J.detJ = J.J[0][0] * J.J[1][1] - J.J[0][1] * J.J[1][0];

	// Macierz odwrotna
	double invDet = 1.0 / J.detJ;
	J.invJ[0][0] = J.J[1][1] * invDet;
	J.invJ[0][1] = -J.J[0][1] * invDet;
	J.invJ[1][0] = -J.J[1][0] * invDet;
	J.invJ[1][1] = J.J[0][0] * invDet;

	return J;
}

void computeJacobiansForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc) {
	int nPoints = npc * npc;
	elem.jacobians.clear();
	elem.jacobians.reserve(nPoints);

	for (int p = 0; p < nPoints; ++p) {
		Jakobian J = computeJacobian(elem, grid, univ.dN_dKsi[p], univ.dN_dEta[p]);
		elem.jacobians.push_back(J);
	}
}

void computeJacobiansForAll(Grid& grid, const ElemUniv& univ, int npc)
{
	for (Element& elem : grid.elements) {
		computeJacobiansForElement(elem, grid, univ, npc);
	}
}