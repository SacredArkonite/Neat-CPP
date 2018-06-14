#pragma once

#include "defines.h"
#include "genome.h"


namespace Population
{

	GEN_PTR CreateGenome(const N_SIZE nIns, const N_SIZE nOuts);
	POP_PTR CreatePop(const int popSize, const N_SIZE inputs, const N_SIZE outputs, C_SIZE& hist);
	float Compatibility(const float c1, const float c2, const float c3, const Genome& genome1, const Genome& genome2);
	bool CheckIfConnectionExists(const std::vector<GEN_PTR>::iterator it, const std::pair<N_SIZE, N_SIZE> connection);
	POP_PTR MutateWeights(POP_PTR pop, float genomeWeightMutation, float weightMutation, float minWeightMutation, float maxWeightMutation);
	POP_PTR MutateStructure(POP_PTR pop, float new_node_percent, float new_link_percent, C_SIZE& hist);
	POP_PTR ClassifyGenomes(POP_PTR encyclopedia, POP_PTR pop, std::vector<uint16_t>& census, float c1, const float c2, const float c3, const float dt);
	POP_PTR UpdateEncyclopedia(const POP_PTR& pop, const std::vector<uint16_t>& census);
	GEN_PTR Mate(const GEN_PTR& genome1, const GEN_PTR& genome2, const float enableChance);
	POP_PTR CreateOffsprings(POP_PTR pop, const std::vector<uint16_t>& census, const float noCrossover, const float enableGeneChance);

}