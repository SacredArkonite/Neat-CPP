#include "simulation.h"


namespace Simulation
{

	//XOR
	std::vector<float> GetSimData(int frame)
	{
		switch (frame % 4)
		{
		case 0:
			return { 0,0,0 };
		case 1:
			return { 0,0,1 };
		case 2:
			return { 0,1,0 };
		case 3:
			return { 0,1,1 };
		default:
			break;
		}
	}

	//XOR
	std::vector<float> GetSimFitness(int frame, std::vector<float> action)
	{
		switch (frame % 4)
		{
		case 0:
			return { action[0] > 0.5f ? 1.0f : 0.0f };
		case 1:
			return { action[0] < 0.5f ? 1.0f : 0.0f };
		case 2:
			return { action[0] < 0.5f ? 1.0f : 0.0f };
		case 3:
			return { action[0] > 0.5f ? 1.0f : 0.0f };
		default:
			break;
		}
	}

}