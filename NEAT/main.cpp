
#include "defines.h"
#include "population.h"
#include "fitness.h"

#include<iostream>

int main()
{
	std::cout << "hello bitches!" << std::endl;

	//Create initial population
	Population pop(150, 3, 1);

	//Generate species dictionnary
	
	for (int i = 0; i < 30; i++)
	{
		//Speciate
		pop.ClassifyGenomes(1.0, 1.0, 0.4, 3.0);

		//Update encyclopedia for next gen classification
		pop.UpdateEncyclopedia();
		
		//Calculate Fitness / Simulate
		pop.CalculateFitness();

		//Adjust Fitness (Explicit fitness sharing)
		pop.ExplicitFitnessSharing(1.0, 1.0, 0.4, 3.0);

		//Reproduce
		pop.CreateOffsprings(0.25, 0.25);
		
		//Mutate Structure
		pop.MutateStructure(0.003, 0.005);
		
		//Mutate weights
		pop.MutateWeights(0.8, 0.9, 0.5, 1.5);
		

		if (i%10 == 0) std::cout << "GEN # " << i << std::endl;
	}
	
	//Calculate Fitness / Simulate
	
	pop.PrintFitness();

	char w;
	std::cin>>w;
	return 0;
}