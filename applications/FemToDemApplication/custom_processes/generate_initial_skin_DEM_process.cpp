//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/generate_initial_skin_DEM_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

GenerateInitialSkinDEMProcess::GenerateInitialSkinDEMProcess(
    ModelPart& rModelPart,
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void GenerateInitialSkinDEMProcess::Execute() 
{
    auto nodal_neigh_process = FindNodalNeighboursProcess(mrModelPart, 5, 5);
    nodal_neigh_process.Execute();
    auto p_DEM_properties = mrDEMModelPart.pGetProperties(1);

    auto &r_submodel_part = mrModelPart.GetSubModelPart("SkinDEMModelPart");
    const auto it_node_begin = r_submodel_part.NodesBegin();
    double num_DEM = 0;
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_submodel_part.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;

        if (!it_node->GetValue(IS_DEM)) { // we have to generate its DEM
            auto& r_neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            Vector potential_radii(r_neigh_nodes.size());
            Vector distances(r_neigh_nodes.size());
    
            // Loop over the neighbours of that node
            bool has_dem_neigh = false;
            for (int neigh = 0; neigh < r_neigh_nodes.size(); neigh++) {
                auto& r_neighbour_node = r_neigh_nodes[neigh];
                const double dist_between_nodes = CalculateDistanceBetweenNodes(it_node, r_neighbour_node);
                distances(neigh) = dist_between_nodes;
                if (r_neighbour_node.GetValue(IS_DEM)) { // Has DEM
                    has_dem_neigh = true;
                    potential_radii(neigh) = dist_between_nodes - r_neighbour_node.GetValue(RADIUS);
                    if (potential_radii(neigh) < 0.0 || potential_radii(neigh) / dist_between_nodes < 0.2) {
                        const double new_radius = dist_between_nodes*0.5;
                        auto& pDEM_particle = r_neighbour_node.GetValue(DEM_PARTICLE_POINTER);
                        auto& r_radius_neigh_old = pDEM_particle->GetGeometry()[0].GetSolutionStepValue(RADIUS);
                        r_radius_neigh_old = new_radius;
                        r_neighbour_node.SetValue(RADIUS, new_radius);
                        potential_radii(neigh) = new_radius;
                    }
                } else {
                    potential_radii(neigh) = 0.0;
                }
            }
            // Let's compute the Radius of the new DEM
            double radius;
            if (has_dem_neigh) {
                radius = this->GetMinimumValue(potential_radii);
            } else {
                radius = this->GetMinimumValue(distances)*0.5;
            }
            const array_1d<double,3>& r_coordinates = it_node->Coordinates();
            const int id = this->GetMaximumDEMId() + 1;
            this->CreateDEMParticle(id, r_coordinates, p_DEM_properties, radius, it_node);
            num_DEM++;
        }
    }
}


/***********************************************************************************/
/***********************************************************************************/

void GenerateInitialSkinDEMProcess::CreateDEMParticle(
    const int Id,
    const array_1d<double, 3> Coordinates,
    const Properties::Pointer pProperties,
    const double Radius,
	NodeIteratorType rNode
)
{
    auto spheric_particle = mParticleCreator.CreateSphericParticleRaw(mrDEMModelPart, Id, Coordinates, pProperties, Radius, "SphericParticle3D");
    rNode->SetValue(IS_DEM, true);
    rNode->SetValue(RADIUS, Radius);
    rNode->SetValue(DEM_PARTICLE_POINTER, spheric_particle);
}


/***********************************************************************************/
/***********************************************************************************/

double GenerateInitialSkinDEMProcess::CalculateDistanceBetweenNodes(
    NodeIteratorType rNode1, 
    const NodeType& rNode2
    )
{
    const double X1 = rNode1->X();
    const double X2 = rNode2.X();
    const double Y1 = rNode1->Y();
    const double Y2 = rNode2.Y();
    const double Z1 = rNode1->Z();
    const double Z2 = rNode2.Z();
    return std::sqrt(std::pow(X1-X2, 2) + std::pow(Y1-Y2, 2) + std::pow(Z1-Z2, 2));
}

/***********************************************************************************/
/***********************************************************************************/

double GenerateInitialSkinDEMProcess::GetMinimumValue(
    const Vector& rValues
    )
{ // this method assumes that the Vector is NOT full of 0.0's
    double aux = 1.0e10;
    for (int i = 0; i < rValues.size(); i++) 
        if (aux > rValues[i] && rValues[i] != 0.0) 
            aux = rValues[i];
    return aux;
}

/***********************************************************************************/
/***********************************************************************************/

int GenerateInitialSkinDEMProcess::GetMaximumDEMId()
{
    int max_id = 0;
    const auto it_DEM_begin = mrDEMModelPart.ElementsBegin();
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrDEMModelPart.Elements().size()); i++) {
        auto it_DEM = it_DEM_begin + i;
        auto& r_geometry = it_DEM->GetGeometry();
        const int DEM_id = r_geometry[0].Id();
        max_id = (max_id < DEM_id) ? DEM_id : max_id;
    }
    return max_id;
}

}  // namespace Kratos
