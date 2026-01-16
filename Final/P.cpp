#include "P.h"
#include "GaussIntegration.h"

void computePForElement(Element& elem, const Grid& grid, const ElemUniv& univ, int npc, double alpha, double Tot)
{
    for (int a = 0; a < 4; ++a)
        elem.P[a] = 0.0;

    GaussLegendreData data = Data(npc);

    int edgeNodes[4][2] = {
        {0, 1}, // krawêdŸ 0
        {1, 2}, // krawêdŸ 1
        {2, 3}, // krawêdŸ 2
        {3, 0}  // krawêdŸ 3
    };

    for (int e = 0; e < 4; ++e) {
        int loc1 = edgeNodes[e][0];
        int loc2 = edgeNodes[e][1];

        const Node& n1 = grid.nodes[elem.nodeIds[loc1] - 1];
        const Node& n2 = grid.nodes[elem.nodeIds[loc2] - 1];

        if (!(n1.BC && n2.BC))
            continue;

        double dx = n2.x - n1.x;
        double dy = n2.y - n1.y;
        double L = sqrt(dx * dx + dy * dy);
        double J1D = L / 2.0;

        for (int k = 0; k < npc; ++k) {
            const array<double, 4>& N = univ.surfaces[e].N[k];
            double w = data.w[k];
            double dS = w * J1D;

            for (int a = 0; a < 4; ++a)
                elem.P[a] += alpha * N[a] * Tot * dS;
        }
    }
}

void computePForAll(Grid& grid, const Grid& g, const ElemUniv& univ, int npc, double alpha, double Tot) {
    for (Element& elem : grid.elements)
        computePForElement(elem, g, univ, npc, alpha, Tot);
}

void addLocalToPG(Equation& eq, const int nodeIds[4], const double Ploc[4]) {
    int I[4];
    for (int a = 0; a < 4; ++a) {
        I[a] = nodeIds[a] - 1;
    }
    for (int a = 0; a < 4; ++a)
        eq.PG[I[a]] += Ploc[a];
}
