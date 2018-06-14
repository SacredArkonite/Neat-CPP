#include "simulation.h"


namespace Simulation
{

	//XOR
	std::vector<float> GetSimData(int frame)
	{
		switch (frame % 1)
		{
		case 0:
			return {1, -1,-1 };
		case 1:
			return {1, -1,1 };
		case 2:
			return {1, 1,-1 };
		case 3:
			return {1, 1,1 };
		default:
			return {};
		}
	}

	//XOR
	std::vector<float> GetSimFitness(int frame, std::vector<float> action)
	{
		switch (frame % 1)
		{
		case 0:
			return { action[0] > 0.0f ? action[0] : 1 };
		case 1:
			return { action[0] < 0.0f ? -action[0] : 1 };
		case 2:
			return { action[0] < 0.0f ? -action[0] : 1 };
		case 3:
			return { action[0] > 0.0f ? action[0] : 1 };
		default:
			return { 0 };
		}
	}

	//XOR
	std::vector<float> GetSimFitness2(int frame, std::vector<float> action)
	{
		switch (frame % 4)
		{
		case 0:
			return { 1.0f + action[0] };
		case 1:
			return { 1.0f - action[0] };
		case 2:
			return { 1.0f - action[0] };
		case 3:
			return { 1.0f + action[0] };
		default:
			return { 0 };
		}
	}

}