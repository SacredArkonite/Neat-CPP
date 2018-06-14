#include "fitness.h"

#include "population.h"
#include "simulation.h"
#include "phenome.h"
#include "genome.h"

#include <iostream>

namespace Fitness
{

	float Simulate(const GEN_PTR& gen)
	{
		float fitness = 0;

		for (int frame = 0; frame < 4; frame++) {
			std::vector<float> inputData = Simulation::GetSimData(frame);
			std::vector<float> actions = Phenome::Propagate(inputData, gen, 2);
			std::vector<float> delta = Simulation::GetSimFitness(frame, actions);

			for each (float d in delta)
				fitness += d*d;
		}
		return 1-fitness/8;
	}

	GEN_PTR CalculateFitness(const POP_PTR& pop)
	{
		float highestFitness = 0;
		GEN_PTR best;
		auto it_end = pop->end();
		for (auto it = pop->begin(); it < it_end; it++)
		{
			float fitness = Simulate(*it);
			(*it)->fitness = fitness +1 / (2 * (*it)->nodes + (*it)->history.size()); //add the size as a second parameter 
			if (highestFitness < fitness) {
				highestFitness = fitness;
				best = std::make_unique<Genome>(**it);
			}
		}
		return best;
	}

	uint16_t SumSharing(const GEN_PTR& gen, const POP_PTR& pop, const float c1, const float c2, const float c3, const float dt)
	{
		auto pop_it_begin = pop->begin();
		auto pop_it = pop_it_begin;
		auto pop_it_end = pop->end();
		uint16_t sum = 0;
		for (; pop_it < pop_it_end; pop_it++)
		{
			sum += (GenomeUtil::Compatibility(c1, c2, c3, *gen, **pop_it) < dt ? 1 : 0);
		}
		return sum;
	}

	void ExplicitFitnessSharing(const POP_PTR& pop, const float c1, const float c2, const float c3, const float dt)
	{
		//Generate a map to all the different species and calculate species global fitness
		auto pop_it_begin = pop->begin();
		auto pop_it = pop_it_begin;
		auto pop_it_end = pop->end();
		std::vector<float> newFitness = std::vector<float>();

		//Calculate new fitness
		for (; pop_it < pop_it_end; pop_it++)
		{
			newFitness.push_back((*pop_it)->fitness / SumSharing(*pop_it, pop, c1, c2, c3, dt));
		}

		//Update fitness
		int i = 0;
		for (pop_it = pop_it_begin; pop_it < pop_it_end; pop_it++, i++)
		{
			(*pop_it)->fitness = newFitness[i];
		}
	}
}
