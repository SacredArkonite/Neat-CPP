#pragma once

#include "defines.h"
#include "genome.h"


class Population
{
public:
	Population(const unsigned int popSize, const N_SIZE inputs, const N_SIZE outputs);
	void ClassifyGenomes(float c1, const float c2, const float c3, const float dt);
	void RefreshEncyclopedia();
	void CalculateFitness();
	void ExplicitFitnessSharing(const float c1, const float c2, const float c3, const float dt);
	void CreateOffsprings(const float noCrossover, const float enableGeneChance);
	void MutateStructure(float new_node_percent, float new_link_percent);
	void MutateWeights(float GenomeMutationProb, float weightMutationProb, float minWeightMutation, float maxWeightMutation);
	void PrintFitness();
	void PrintHighestFitness();

private:
	void SortBySpeciesFitness();

	POP_PTR pop;
	POP_PTR encyclopedia;
	GEN_PTR best;
	std::vector<unsigned int> census;
	C_SIZE innovationNumber;
	unsigned short speciesNumber;
	float highestFitness;
	int generation;
	int VIPCount;
};