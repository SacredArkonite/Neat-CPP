
#include "defines.h"
#include "population.h"
#include "fitness.h"

#include<iostream>

int main()
{
	std::cout << "hello bitches!" << std::endl;

	//Create initial population
	Population pop(150, 3, 1);

	//Speciate
	pop.ClassifyGenomes(1.0f, 1.0f, 0.4f, 3.0f);

	//Calculate Fitness / Simulate
	pop.CalculateFitness();

	//Adjust Fitness (Explicit fitness sharing)
	pop.ExplicitFitnessSharing(1.0f, 1.0f, 0.4f, 3.0f);
	
	for (int i = 0; i < 25; i++)
	{

		//Reproduce
		pop.CreateOffsprings(0.25, 0.25);
		
		//Mutate Structure
		pop.MutateStructure(0.03f, 0.05f);
		
		//Mutate weights
		pop.MutateWeights(0.8f, 0.9f, 0.5f, 1.5f);

		//Speciate
		pop.ClassifyGenomes(1.0f, 1.0f, 0.4f, 3.0f);

		//Calculate Fitness / Simulate
		pop.CalculateFitness();

		//Adjust Fitness (Explicit fitness sharing)
		pop.ExplicitFitnessSharing(1.0f, 1.0f, 0.4f, 3.0f);

		pop.PrintHighestFitness();
	}
	
	//Calculate Fitness / Simulate
	pop.PrintFitness();

	char w;
	std::cin>>w;
	return 0;
}