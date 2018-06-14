#pragma once

#include <vector>

#include<memory>

#include "defines.h"

struct Genome
{
	std::vector<N_SIZE> sourceNode;
	std::vector<N_SIZE> destNode;
	std::vector<N_SIZE> inputNode;
	std::vector<N_SIZE> outputNode;
	std::vector<C_SIZE> history;
	std::vector<C_SIZE> disabledIndex;
	std::vector<float> weights;
	N_SIZE nodes = 0;
	size_t evolutionHash = 0;
	uint16_t species = 0;
	float fitness = 0.0f;

	bool operator < (const Genome& str) const
	{
		return (species < str.species);
	}
};

typedef std::unique_ptr<Genome> GEN_PTR;
typedef std::unique_ptr<std::vector<GEN_PTR>> POP_PTR;

namespace GenomeUtil
{

	GEN_PTR CreateGenome(const N_SIZE nIns, const N_SIZE nOuts);
	float Compatibility(const float c1, const float c2, const float c3, const Genome& genome1, const Genome& genome2);
	bool CheckIfConnectionExists(const GEN_PTR& gen, const std::pair<N_SIZE, N_SIZE> connection);
	GEN_PTR Mate(const GEN_PTR& genome1, const GEN_PTR& genome2, const float enableChance);
}