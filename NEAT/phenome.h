#pragma once

#include "genome.h"

#include <vector>

namespace Phenome
{
	inline float ActivationFunction(const float value);
	std::vector<float> Propagate(const std::vector<float>& input, const GEN_PTR& gen, uint32_t maxSteps = 5);
}