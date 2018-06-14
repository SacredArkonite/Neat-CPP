#pragma once

#include <vector>

#include<memory>

#include "defines.h"

struct Genome
{
	//Genome(const Genome & src) = default;
	/*sourceNode(src.sourceNode),
	destNode(src.destNode),
	history(src.history),
	disabledIndex(src.disabledIndex),
	weights(src.weights),
	nodes(src.nodes),
	evolutionHash(src.evolutionHash),
	species(src.species),
	fitness(src.fitness)
	{}*/

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