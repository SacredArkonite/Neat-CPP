
#include "defines.h"
#include "population.h"
#include "fitness.h"

#include<iostream>

int main()
{
	std::cout << "hello bitches!" << std::endl;

	//Create initial population
	C_SIZE hist;
	POP_PTR pop = Population::CreatePop(150, 3, 1, hist);

	//Generate species dictionnary
	POP_PTR speciesEncyclopedia = Population::CreateEncyclopedia((*pop)[0]);
	std::vector<uint16_t> census;

	for (int i = 0; i < 30; i++)
	{
		//Speciate
		pop = Population::ClassifyGenomes(std::move(speciesEncyclopedia), std::move(pop), census, 1, 1, 0, 4);

		//Update encyclopedia for next gen classification
		speciesEncyclopedia = Population::UpdateEncyclopedia(pop, census);

		//Calculate Fitness / Simulate
		pop = Fitness::CalculateFitness(std::move(pop));

		//Adjust Fitness (Explicit fitness sharing)
		pop = Fitness::ExplicitFitnessSharing(std::move(pop), 1.0, 1.0, 0.4, 3.0);

		//Reproduce
		pop = Population::CreateOffsprings(std::move(pop), census, 0.25, 0.25);

		//Mutate Structure
		pop = Population::MutateStructure(std::move(pop), 0.03, 0.05, hist);
		
		//Mutate weights
		pop = Population::MutateWeights(std::move(pop), 0.8, 0.9, 0.5, 1.5);


		if (i%10 == 0) std::cout << "GEN # " << i << std::endl;
	}

	//Calculate Fitness / Simulate
	pop = Fitness::CalculateFitness(std::move(pop));

	for (auto it = pop->begin(); it < pop->end(); it++)
	{
		std::cout << (*it)->fitness << std::endl;
	}
	char w;
	std::cin>>w;
	return 0;
}