#include<iostream>
#include<vector>
#include<memory>
#include<algorithm>
#include<random>
#include<map>

#include "RNG.h"
#include "defines.h"
#include "hasher.h"
#include "genome.h"
#include "mutate.h"
#include "population.h"
#include "simulation.h"
#include "phenome.h"

RNG rng;


float Simulate(const GEN_PTR& gen)
{
	float fitness = 0;

	for (int frame = 0; frame < 4; frame++) {
		std::vector<float> inputData = Simulation::GetSimData(frame);
		std::vector<float> actions = Phenome::Propagate(inputData, gen, 2);
		std::vector<float> delta = Simulation::GetSimFitness(frame, actions);

		for each (float d in delta)
			fitness += d;
	}

	return fitness;
}

POP_PTR CalculateFitness(POP_PTR pop)
{
	auto it_end = pop->end();
	for (auto it = pop->begin(); it < it_end; it++)
	{
		(*it)->fitness = Simulate(*it);
	}
	return pop;
}

uint16_t SumSharing(const GEN_PTR& gen, const POP_PTR& pop, const float c1, const float c2, const float c3, const float dt)
{
	auto pop_it_begin = pop->begin();
	auto pop_it = pop_it_begin;
	auto pop_it_end = pop->end();
	uint16_t sum = 0;
	for (; pop_it < pop_it_end; pop_it++)
	{
		sum += (Population::Compatibility(c1, c2, c3, *gen, **pop_it) < dt ? 1.0f : 0.0f);
	}
	return sum;
}

POP_PTR ExplicitFitnessSharing(POP_PTR pop, const float c1, const float c2, const float c3, const float dt)
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
	for (; pop_it < pop_it_end; pop_it++, i++)
	{
		(*pop_it)->fitness = newFitness[i];
	}

	return pop;
}

int main()
{
	std::cout << "hello bitches!" << std::endl;
/*
	C_SIZE hist;
	POP_PTR pop = GenerateExample(hist);
	GEN_PTR offspring = Mate(pop->at(0), pop->at(1),0.25);
	*/


	//float delta = Compatibility(1, 1, 1, *(*pop)[0], *(*pop)[1]);

	//Create initial population
	C_SIZE hist;
	POP_PTR pop = Population::CreatePop(150, 3, 1, hist);

	//Generate species dictionnary
	POP_PTR speciesEncyclopedia = std::make_unique<std::vector<GEN_PTR>>();;
	speciesEncyclopedia->push_back(std::make_unique<Genome>(*(pop->at(rng.LessThan(10)))));
	std::vector<uint16_t> census;

	for (int i = 0; i < 30; i++)
	{
		//Speciate
		pop = Population::ClassifyGenomes(std::move(speciesEncyclopedia), std::move(pop), census, 1, 1, 0, 4);

		//Update encyclopedia for next gen classification
		speciesEncyclopedia = Population::UpdateEncyclopedia(pop, census);

		//Calculate Fitness / Simulate
		pop = CalculateFitness(std::move(pop));

		//Adjust Fitness (Explicit fitness sharing)
		pop = ExplicitFitnessSharing(std::move(pop), 1.0, 1.0, 0.4, 3.0);

		//Reproduce
		pop = Population::CreateOffsprings(std::move(pop), census, 0.25, 0.25);

		//Mutate Structure
		pop = Population::MutateStructure(std::move(pop), 0.03, 0.05, hist);
		
		//Mutate weights
		pop = Population::MutateWeights(std::move(pop), 0.8, 0.9, 0.5, 1.5);


		if (i%10 == 0) std::cout << "GEN # " << i << std::endl;
	}

	//Calculate Fitness / Simulate
	pop = CalculateFitness(std::move(pop));

	for (auto it = pop->begin(); it < pop->end(); it++)
	{
		std::cout << (*it)->fitness << std::endl;
	}
	char w;
	std::cin>>w;
	return 0;
}