#include "population.h"

#include "mutate.h"
#include "fitness.h"
#include "RNG.h"
#include "phenome.h"

#include<map>

#include<iostream>
#include <algorithm>

	
Population::Population(const unsigned int popSize, const N_SIZE inputs, const N_SIZE outputs)
{
	//Generate a vector of genome
	pop = std::make_unique<std::vector<GEN_PTR>>();
	speciesNumber = 0;
	innovationNumber = 0;
	generation = 0;
	highestFitness = 0;
	VIPCount = 0;

	for (unsigned int i = 0; i < popSize; i++)
	{
		pop->push_back(GenomeUtil::CreateGenome(inputs, outputs));
	}

	//Actualize the noevlty number based on the initial genome
	innovationNumber = (*pop)[0]->history.size() + 1;

	//Create the initial encyclopedia
	encyclopedia = std::make_unique<std::vector<GEN_PTR>>();
	encyclopedia->push_back(std::make_unique<Genome>(*((*pop)[0])));

	//Update the census
	census.push_back(popSize);

}

void Population::SortBySpeciesFitness()
{
	std::sort(pop->begin(), pop->end(), [](const GEN_PTR& gen1, const GEN_PTR& gen2)
	{ return (gen1->species == gen2->species) ? (gen1->fitness > gen2->fitness) : (gen1->species < gen2->species); });
}
	
void Population::ClassifyGenomes(float c1, const float c2, const float c3, const float dt)
{
	//Update all the genome's species, update the encyclopedia and build a census
	census = std::vector<unsigned int>(speciesNumber + 1);

	auto it_end = pop->end();
	for (auto it = pop->begin(); it != it_end; it++)
	{
		auto dic_begin = encyclopedia->begin();
		auto dic_end = encyclopedia->end();
		bool found = false;
		for (auto dic = dic_begin; dic != dic_end; dic++)
		{
			if (GenomeUtil::Compatibility(c1, c2, c3, **it, **dic) < dt)
			{
				uint16_t speciesNb = (*dic)->species;
				(*it)->species = speciesNb;
				census[speciesNb] += 1;
				found = true;
				break;
			}
		}
		if (!found)
		{
			(*it)->species = ++speciesNumber;
			encyclopedia->push_back(std::make_unique<Genome>(**it));
			census.push_back(1);
		}
	}

	//Refresh to randomly migrate the the species references in the search space
	RefreshEncyclopedia();
}
	
void Population::RefreshEncyclopedia()
{
	//Refresh the encyclopedia by selecting a new random genome of each available species
	
	uint16_t specieLocation = 0;

	auto it_end = census.end();
	for (auto it = census.begin(); it != it_end; it++)
	{
		if (*it != 0) {
			GEN_PTR specimen = std::make_unique<Genome>(*(pop->at(specieLocation + RNG::LessThan(*it))));
			encyclopedia->at(specimen->species) = std::move(specimen);
			specieLocation += *it;
		}
	}
}

void Population::CalculateFitness()
{
	best = Fitness::CalculateFitness(pop);
	Fitness::CalculateFitness(encyclopedia);
	highestFitness = best->fitness -1 / (2 * best->nodes + best->history.size());
}

void Population::ExplicitFitnessSharing(const float c1, const float c2, const float c3, const float dt)
{
	Fitness::ExplicitFitnessSharing(pop, c1, c2, c3, dt);

	//Sort the population by species/fitness
	SortBySpeciesFitness();
}

void Population::CreateOffsprings(const float noCrossover, const float enableGeneChance)
{
	generation++;
	VIPCount = 0;

	//Create an entire new gen
	POP_PTR nexxgen = std::make_unique<std::vector<GEN_PTR>>();

	//Generate a map to all the different species and calculate species global fitness
	auto pop_it_begin = pop->begin();
	auto pop_it = pop_it_begin;
	auto pop_it_end = pop->end();
	std::vector<std::pair<uint16_t, uint16_t>> speciesIndex;
	std::vector<float> speciesFitness;
	speciesIndex.push_back({ (*pop_it)->species, 0 });
	speciesFitness.push_back(0);
	float populationFitness = 0;

	while (pop_it < pop_it_end)
	{
		if ((*pop_it)->species != speciesIndex.back().first)
		{
			speciesIndex.push_back({ (*pop_it)->species, std::distance(pop_it_begin,pop_it) });
			speciesFitness.push_back(0);
		}
		speciesFitness.back() += (*pop_it)->fitness;
		populationFitness += (*pop_it)->fitness;

		pop_it++;
	}

	//Select MVPS
	unsigned int nb_species = census.size();
	auto speciesIndex_it_begin = speciesIndex.begin();
	auto speciesIndex_it = speciesIndex_it_begin;
	auto speciesIndex_it_end = speciesIndex.end();
	for (uint16_t species = 0; species < nb_species; species++)
	{
		//Select the best from each species with at least 5
		if (census[species] >= 5)
		{	//Should be safe, we just created the index
			while (speciesIndex_it->first != species) speciesIndex_it++;
			auto left = pop_it_begin + speciesIndex_it->second;
			auto right = (speciesIndex_it + 1 == speciesIndex_it_end) ? pop_it_end : (pop_it_begin + (speciesIndex_it + 1)->second);

			auto mvp_index = std::max_element(left, right,
				[]
			(const GEN_PTR& gen1, const GEN_PTR& gen2)
			{return gen1->fitness < gen2->fitness; });

			auto mvp_clone = std::make_unique<Genome>(**mvp_index);
			nexxgen->push_back(std::move(mvp_clone));
			VIPCount++;
		}
	}

	//Calculate how many offsprings per species we want
	unsigned int pop_size = pop->size() - nexxgen->size();
	std::vector<uint16_t> requiredChilds = std::vector<uint16_t>();
	int16_t spotsLeft = pop_size;
	for each (float speciesFit in speciesFitness)
	{
		uint16_t quantity = (int)(speciesFit / populationFitness * pop_size) + 1;
		requiredChilds.push_back(quantity);
		spotsLeft -= quantity;
	}
	unsigned int speciesCount = requiredChilds.size();
	//Fill leftover or remove excess at random
	if (spotsLeft < 0)
	{
		for (; spotsLeft != 0; spotsLeft++)
		{
			requiredChilds[RNG::LessThan(speciesCount)] -= 1;
			//Species should DIE here!!
			requiredChilds.erase(std::remove(requiredChilds.begin(), requiredChilds.end(), 0), requiredChilds.end());
			speciesCount = requiredChilds.size();

		}
	}
	else if (spotsLeft > 0)
	{
		for (; spotsLeft != 0; spotsLeft--)
		{
			requiredChilds[RNG::LessThan(speciesCount)] += 1;
		}
	}

	//Create offsprings for each species
	for (uint16_t species = 0; species < speciesCount; species++)
	{
		for (int childCount = 0; childCount < requiredChilds[species]; childCount++)
		{
			//Mutation without crossover
			if (RNG::RngProb() < noCrossover)
			{
				//Select 1 genome at random based on fitness wheel
				auto left = pop->begin() + (speciesIndex[species].second);
				float fitnessArrow = RNG::RngProb() * speciesFitness[species] / 2; //Only chose from the top 50% of assets
				while (fitnessArrow >(*left)->fitness)
				{
					fitnessArrow -= (*left)->fitness;
					left++;
				}
				nexxgen->push_back(std::make_unique<Genome>(**left));
			}
			else
			{
				//Select 2 genomes at random based on fitness wheel
				auto left = pop->begin() + speciesIndex[species].second;
				float fitnessArrow = RNG::RngProb() * speciesFitness[species] / 2; //Only chose from the top 50% of assets
				while (fitnessArrow > (*left)->fitness)
				{
					fitnessArrow -= (*left)->fitness;
					left++;
				}
				GEN_PTR parent1 = std::make_unique<Genome>(**left);

				left = pop->begin() + speciesIndex[species].second;
				fitnessArrow = RNG::RngProb() * speciesFitness[species] / 2; //Only chose from the top 50% of assets
				while (fitnessArrow > (*left)->fitness)
				{
					fitnessArrow -= (*left)->fitness;
					left++;
				}
				GEN_PTR parent2 = std::make_unique<Genome>(**left);

				//make a child
				nexxgen->push_back(GenomeUtil::Mate(parent1, parent2, enableGeneChance));
			}
		}
	}

	pop = std::move(nexxgen);
}
	
void Population::MutateStructure(float new_node_percent, float new_link_percent)
{
	//making a new node and a new connection should be mutually exclusive
	new_link_percent = new_link_percent + new_node_percent;


	//Prepping next gen
	POP_PTR nexxgen = std::make_unique<std::vector<GEN_PTR>>();

	//We use a multimap to order the genome by their {type,operation} 
	//in order to identify common types after the transform
	//so we can give the same history number to the same transform
	std::multimap<std::pair<size_t, C_SIZE>, GEN_PTR> mm_node;
	std::multimap<std::pair<size_t, std::pair<N_SIZE, N_SIZE>>, GEN_PTR> mm_conn;

	auto end = (*pop).end();

	//Ignore MVPS
	for (auto it = (*pop).begin(); it < (*pop).begin() + VIPCount; it++)
		nexxgen->push_back(std::move(*it));

	//Distribute the genome with their transform into the maps
	for (auto it = (*pop).begin() + VIPCount; it < end; it++)
	{
		float rn = RNG::RngProb();
		if (rn < new_node_percent)
		{	// Mutate Add Node
			//Do not mutate if we reached the max number of nodes
			if ((*it)->nodes >= MAX_NODES)
			{
				nexxgen->push_back(std::move(*it));
			}
			else
			{
				size_t hash = (*it)->evolutionHash;
				C_SIZE connection_id = RNG::LessThan((*it)->history.size());
				auto key = std::pair<size_t, C_SIZE>(hash, connection_id);
				mm_node.insert({ key, std::move(*it) });
			}
		}
		else if (rn < new_link_percent)
		{	// Mutate Add Connection
			//Do not mutate if we reached the max number of connections
			if ((*it)->sourceNode.size() >= MAX_CONNECTIONS)
			{
				nexxgen->push_back(std::move(*it));
			}
			else
			{
				size_t hash = (*it)->evolutionHash;
				N_SIZE nb_nodes = (*it)->nodes;
				auto connection = std::pair<N_SIZE, N_SIZE>(RNG::LessThan(nb_nodes), RNG::LessThan(nb_nodes));

				//Check if this connection already exists
				if (GenomeUtil::CheckIfConnectionExists(*it, connection))
				{
					nexxgen->push_back(std::move(*it));
				}
				else
				{
					auto key = std::pair<size_t, std::pair<N_SIZE, N_SIZE>>(hash, connection);
					mm_conn.insert({ key, std::move(*it) });
				}
			}
		}
		else
		{
			// Keep as is
			nexxgen->push_back(std::move(*it));
		}

	}

	//Iterate over the new nodes multimap in key order
	if (!mm_node.empty())
	{
		std::multimap<std::pair<size_t, C_SIZE>, GEN_PTR>::iterator it = mm_node.begin();
		std::multimap<std::pair<size_t, C_SIZE>, GEN_PTR>::iterator end = mm_node.end();
		std::pair<size_t, C_SIZE> last_key;
		for (last_key = it->first; it != end; it++) {
			//Increment the mutation ID if when needed
			if (it->first != last_key)
			{
				last_key = it->first;
				innovationNumber += 2;
			}
			Mutate::AddNode(it->second, last_key.second, innovationNumber);
			nexxgen->push_back(std::move(it->second));
		}

		innovationNumber += 2;
	}

	//Iterate over the new conx multimap in key order
	if (!mm_conn.empty())
	{
		std::multimap<std::pair<size_t, std::pair<N_SIZE, N_SIZE>>, GEN_PTR>::iterator it = mm_conn.begin();
		std::multimap<std::pair<size_t, std::pair<N_SIZE, N_SIZE>>, GEN_PTR>::iterator end = mm_conn.end();
		std::pair<size_t, std::pair<N_SIZE, N_SIZE>> last_key = it->first;
		for (; it != end; it++) {
			//Increment the mutation ID if when needed
			if (it->first != last_key)
			{
				last_key = it->first;
				innovationNumber++;
			}
			Mutate::AddConnection(it->second, last_key.second.first, last_key.second.second, innovationNumber);
			nexxgen->push_back(std::move(it->second));
		}
		innovationNumber++;
	}

	pop = std::move(nexxgen);
}
	
void Population::MutateWeights(float genomeMutationProb, float weightMutationProb, float minWeightMutation, float maxWeightMutation)
{
	//We mutate by scaling randomly
	int c = 0;

	auto pop_end = pop->end();
	for (auto it = pop->begin() + VIPCount; it != pop_end; it++) {
		//Mutate the genome?
		if (RNG::RngProb() < genomeMutationProb)
		{
			auto i_w_end = (*it)->weights.end();
			for (auto i_w = (*it)->weights.begin(); i_w != i_w_end; i_w++)
			{
				//Mutate the weight by scaling or +/-
				if (RNG::RngProb() < weightMutationProb)
				{
					(*i_w) *= RNG::RngRange();
				}
				//Mutate random (same as initial)
				else
				{
					// change sing?
					(*i_w) = RNG::RngWeight();
				}
			}
		}
		c++;
	}
}
	
void Population::PrintFitness()
{
	CalculateFitness();
	for (auto it = pop->begin(); it < pop->end(); it++)
	{
		std::cout << (*it)->fitness * (2 * (*it)->nodes + (*it)->history.size()) << std::endl;
	}
}

void Population::PrintHighestFitness()
{
	std::cout << "GEN " << generation << " Highest fitness " << highestFitness << std::endl;
}
