#pragma once

#include "genome.h"

namespace Fitness
{
	float Simulate(const GEN_PTR& gen);
	POP_PTR CalculateFitness(POP_PTR pop);
	uint16_t SumSharing(const GEN_PTR& gen, const POP_PTR& pop, const float c1, const float c2, const float c3, const float dt);
	POP_PTR ExplicitFitnessSharing(POP_PTR pop, const float c1, const float c2, const float c3, const float dt);
}
