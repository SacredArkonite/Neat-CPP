#pragma once

#include<vector>

namespace Simulation
{

	std::vector<float> GetSimData(int frame);
	std::vector<float> GetSimFitness(int frame, std::vector<float> action);
	std::vector<float> GetSimFitness2(int frame, std::vector<float> action);

}