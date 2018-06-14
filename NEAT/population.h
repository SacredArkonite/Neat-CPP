#pragma once

#include "defines.h"
#include "genome.h"


namespace Population
{

	GEN_PTR CreateGenome(const N_SIZE nIns, const N_SIZE nOuts);
	POP_PTR CreatePop(const int popSize, const N_SIZE inputs, const N_SIZE outputs, C_SIZE& hist);
	float Compatibility(const float c1, const float c2, const float c3, const Genome& genome1, const Genome& genome2);
	bool CheckIfConnectionExists(const std::vector<GEN_PTR>::iterator it, const std::pair<N_SIZE, N_SIZE> connection);
	GEN_PTR Mate(const GEN_PTR& genome1, const GEN_PTR& genome2, const float enableChance);
	
	class Population
	{
	public:
		Population(const unsigned int popSize, const N_SIZE inputs, const N_SIZE outputs);
		void ClassifyGenomes(float c1, const float c2, const float c3, const float dt);
		void UpdateEncyclopedia();
		void CalculateFitness();
		void ExplicitFitnessSharing(const float c1, const float c2, const float c3, const float dt);
		void CreateOffsprings(const float noCrossover, const float enableGeneChance);
		void MutateStructure(float new_node_percent, float new_link_percent);
		void MutateWeights(float GenomeMutationProb, float weightMutationProb, float minWeightMutation, float maxWeightMutation);
		void PrintFitness();

	private:
		POP_PTR pop;
		POP_PTR encyclopedia;
		std::vector<unsigned int> census;
		C_SIZE innovationNumber;
	};
}