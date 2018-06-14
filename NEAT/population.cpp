#include "population.h"

#include "mutate.h"
#include "hasher.h"
#include "RNG.h"

#include<map>

namespace Population
{
	RNG rng;

	GEN_PTR CreateGenome(const N_SIZE nIns, const N_SIZE nOuts)
	{
		auto genome = std::make_unique<Genome>();
		C_SIZE hist = 1;
		genome->nodes = 0;
		genome->species = 0;
		genome->fitness = 0.0;

		for (N_SIZE i = 0; i < nIns; i++)
			genome = Mutate::AddInput(std::move(genome), hist);

		for (N_SIZE j = 0; j < nOuts; j++)
			genome = Mutate::AddOutput(std::move(genome), hist);

		return genome;
	}

	POP_PTR CreatePop(const int popSize, const N_SIZE inputs, const N_SIZE outputs, C_SIZE& hist)
	{
		POP_PTR pop = std::make_unique<std::vector<GEN_PTR>>();

		for (int i = 0; i < popSize; i++)
		{
			pop->push_back(CreateGenome(inputs, outputs));
		}

		//Actualize the noevlty ID based on the initial size
		hist = (*pop)[0]->history.size() + 1;
		return pop;
	}

	float Compatibility(const float c1, const float c2, const float c3, const Genome& genome1, const Genome& genome2)
	{
		//Count the number of disjoint and excess between the two genes
		auto it1 = genome1.history.begin();
		auto it2 = genome2.history.begin();
		auto end1 = genome1.history.end();
		auto end2 = genome2.history.end();

		auto w1 = genome1.weights.begin();
		auto w2 = genome2.weights.begin();

		int disjoint = 0;
		int excess = 0;
		int matching = 0;
		float d_sum = 0;

		while (it1 != end1 && it2 != end2)
		{
			if (*it1 == *it2)
			{	//same
				//Calculate average weight difference of matching genes
				matching++;
				d_sum += std::abs(*w1 - *w2);

				it1++;
				it2++;
				w1++;
				w2++;
			}
			else if (*it1 < *it2)
			{	//gene1 has disjoint
				disjoint++;
				it1++;
				w1++;
			}
			else
			{	//gene2 has disjoint
				disjoint++;
				it2++;
				w2++;
			}
		}
		excess = (end1 - it1) + (end2 - it2);

		return (c1*excess + c2*disjoint + c3*d_sum / matching);
	}

	POP_PTR GenerateExample(C_SIZE& hist)
	{
		POP_PTR pop = CreatePop(2, 3, 1, hist);
		(*pop)[0] = Mutate::AddNode(std::move((*pop)[0]), 1, 4);
		(*pop)[1] = Mutate::AddNode(std::move((*pop)[1]), 1, 4);
		(*pop)[1] = Mutate::AddNode(std::move((*pop)[1]), 4, 6);
		(*pop)[0] = Mutate::AddConnection(std::move((*pop)[0]), 0, 4, 8);
		(*pop)[1] = Mutate::AddConnection(std::move((*pop)[1]), 2, 4, 9);
		(*pop)[1] = Mutate::AddConnection(std::move((*pop)[1]), 0, 5, 10);

		return pop;
	}

	bool CheckIfConnectionExists(const std::vector<GEN_PTR>::iterator it, const std::pair<N_SIZE, N_SIZE> connection)
	{
		auto sr_begin = (*it)->sourceNode.begin();
		auto sr_ptr = sr_begin;
		auto sr_end = (*it)->sourceNode.end();
		C_SIZE index = -1; // We crash if not found
		for (; sr_ptr != sr_end; sr_ptr++) {
			if ((*sr_ptr) == connection.first)
			{
				index = sr_ptr - (*it)->sourceNode.begin();
				if ((*it)->destNode[index] == connection.second)
				{	//FOUND IT!
					//Remove from the disabled list if its there
					auto pos_in_disabled = std::find((*it)->disabledIndex.begin(), (*it)->disabledIndex.end(), index);
					if (pos_in_disabled != (*it)->disabledIndex.end())
						(*it)->disabledIndex.erase(pos_in_disabled);
					return true;
				}
			}
		}
		return false;
	}

	POP_PTR MutateWeights(POP_PTR pop, float genomeWeightMutation, float weightMutation, float minWeightMutation, float maxWeightMutation)
	{
		//We mutate by scaling randomly
		int c = 0;

		auto pop_end = pop->end();
		for (auto it = pop->begin(); it != pop_end; it++) {
			//Mutate the genome?
			if (rng.RngProb() < genomeWeightMutation)
			{
				auto i_w_end = (*it)->weights.end();
				for (auto i_w = (*it)->weights.begin(); i_w != i_w_end; i_w++)
				{
					//Mutate the weight by scaling?
					if (rng.RngProb() < weightMutation)
					{
						(*i_w) *= rng.RngRange();
					}
					//Mutate random (same as initial)
					else
					{
						// change sing?
						(*i_w) = rng.RngWeight();
					}
				}
			}
			c++;
		}

		return pop;
	}

	POP_PTR MutateStructure(POP_PTR pop, float new_node_percent, float new_link_percent, C_SIZE& hist)
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

		//Distribute the genome with their transform into the maps
		for (auto it = (*pop).begin(); it < end; it++)
		{
			float rn = rng.RngProb();
			if (rn < new_node_percent)
			{	// Mutate Add Node
				//Do not mutate if we reached the max number of nodes
				if ((*it)->nodes == MAX_NODES)
				{
					nexxgen->push_back(std::move(*it));
				}
				else
				{
					size_t hash = (*it)->evolutionHash;
					C_SIZE connection_id = rng.LessThan((*it)->history.size());
					auto key = std::pair<size_t, C_SIZE>(hash, connection_id);
					mm_node.insert({ key, std::move(*it) });
				}
			}
			else if (rn < new_link_percent)
			{	// Mutate Add Connection
				//Do not mutate if we reached the max number of connections
				if ((*it)->sourceNode.size() == MAX_NODES)
				{
					nexxgen->push_back(std::move(*it));
				}
				else
				{
					size_t hash = (*it)->evolutionHash;
					N_SIZE nb_nodes = (*it)->nodes;
					auto connection = std::pair<N_SIZE, N_SIZE>(rng.LessThan(nb_nodes), rng.LessThan(nb_nodes));

					//Check if this connection already exists
					if (CheckIfConnectionExists(it, connection))
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
					hist += 2;
				}
				nexxgen->push_back(Mutate::AddNode(std::move(it->second), last_key.second, hist));
			}

			hist += 2;
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
					hist++;
				}
				nexxgen->push_back(Mutate::AddConnection(std::move(it->second), last_key.second.first, last_key.second.second, hist));
			}
			hist++;
		}

		return nexxgen;
	}

	int debug_counter = 0;
	POP_PTR ClassifyGenomes(POP_PTR encyclopedia, POP_PTR pop, std::vector<uint16_t>& census, float c1, const float c2, const float c3, const float dt)
	{
		//Update all the genome's species, build a census and order the population by species
		census = std::vector<uint16_t>(encyclopedia->size());
		auto it_end = pop->end();
		for (auto it = pop->begin(); it != it_end; it++)
		{
			auto dic_begin = encyclopedia->begin();
			auto dic_end = encyclopedia->end();
			bool found = false;
			for (auto dic = dic_begin; dic != dic_end; dic++)
			{
				if (Compatibility(c1, c2, c3, **it, **dic) < dt)
				{
					uint16_t specieNb = dic - dic_begin;
					(*it)->species = specieNb;
					census[specieNb] += 1;
					found = true;
					break;
				}
			}
			if (!found)
			{
				encyclopedia->push_back(std::make_unique<Genome>(**it));
				(*it)->species = encyclopedia->size();
				census.push_back(1);
			}
		}

		std::sort(pop->begin(), pop->end(), [](const GEN_PTR& gen1, const GEN_PTR& gen2)
		{ return (gen1->species == gen2->species) ? (gen1->fitness > gen2->fitness) : (gen1->species < gen2->species); });

		return pop;
	}

	POP_PTR UpdateEncyclopedia(const POP_PTR& pop, const std::vector<uint16_t>& census)
	{
		//Select a random individual for each species
		POP_PTR encyclopedia = std::make_unique<std::vector<GEN_PTR>>();
		uint16_t specieLocation = 0;

		auto it_end = census.end();
		for (auto it = census.begin(); it != it_end; it++)
		{
			if (*it != 0) {
				GEN_PTR specimen = std::make_unique<Genome>(*(pop->at(specieLocation + rng.LessThan(*it))));
				encyclopedia->push_back(std::move(specimen));
				specieLocation += *it;
			}
		}

		return encyclopedia;
	}

	GEN_PTR Mate(const GEN_PTR& genome1, const GEN_PTR& genome2, const float enableChance)
	{
		//Create offspring. genome1 has priority on weights
		auto it1 = 0;
		auto it2 = 0;
		auto end1 = genome1->history.size();
		auto end2 = genome2->history.size();

		auto dis_it1 = genome1->disabledIndex.begin();
		auto dis_it2 = genome2->disabledIndex.begin();
		auto dis_it1_end = genome1->disabledIndex.end();
		auto dis_it2_end = genome2->disabledIndex.end();
		uint16_t dis_offset_1 = 0;
		uint16_t dis_offset_2 = 0;

		GEN_PTR offspring = std::make_unique<Genome>();
		while (it1 < end1 && it2 < end2)
		{
			if (genome1->history[it1] == genome2->history[it2])
			{	//same

				//Inherits weights from random parent
				offspring->history.push_back(genome1->history[it1]);
				offspring->sourceNode.push_back(genome1->sourceNode[it1]);
				offspring->destNode.push_back(genome1->destNode[it1]);
				offspring->weights.push_back(rng.RngBool() ? genome1->weights[it1] : genome2->weights[it2]);

				//Inherits disabled links
				if (dis_it1 != dis_it1_end && *dis_it1 == it1)
				{
					offspring->disabledIndex.push_back(it1 + dis_offset_1);
					dis_it1++;
				}
				if (dis_it2 != dis_it2_end && *dis_it2 == it2)
				{
					//Do not add if already there from the other parent
					if (offspring->disabledIndex.size()>0 && offspring->disabledIndex.back() != it2 + dis_offset_2)
						offspring->disabledIndex.push_back(it2 + dis_offset_2);
					dis_it2++;
				}

				it1++;
				it2++;
			}
			else if (genome1->history[it1] < genome2->history[it2])
			{	//gene1 has disjoint

				//Inherits from genome 1
				offspring->history.push_back(genome1->history[it1]);
				offspring->sourceNode.push_back(genome1->sourceNode[it1]);
				offspring->destNode.push_back(genome1->destNode[it1]);
				offspring->weights.push_back(genome1->weights[it1]);

				//Inherits disabled links
				if (dis_it1 != dis_it1_end && *dis_it1 == it1)
				{
					offspring->disabledIndex.push_back(it1 + dis_offset_1);
					dis_it1++;
				}

				it1++;
				dis_offset_2++;
			}
			else
			{	//gene2 has disjoint

				//Inherits from genome 2
				offspring->history.push_back(genome2->history[it2]);
				offspring->sourceNode.push_back(genome2->sourceNode[it2]);
				offspring->destNode.push_back(genome2->destNode[it2]);
				offspring->weights.push_back(genome2->weights[it2]);

				//Inherits disabled links
				if (dis_it2 != dis_it2_end && *dis_it2 == it2)
				{
					offspring->disabledIndex.push_back(it2 + dis_offset_2);
					dis_it2++;
				}

				it2++;
				dis_offset_1++;
			}
		}

		//Then, add the excess genes
		//Try for genome1
		while (it1 < end1)
		{
			//Inherits from genome 1
			offspring->history.push_back(genome1->history[it1]);
			offspring->sourceNode.push_back(genome1->sourceNode[it1]);
			offspring->destNode.push_back(genome1->destNode[it1]);
			offspring->weights.push_back(genome1->weights[it1]);

			//Inherits disabled links
			if (dis_it1 != dis_it1_end && *dis_it1 == it1)
			{
				offspring->disabledIndex.push_back(it1 + dis_offset_1);
				dis_it1++;
			}

			it1++;
		}
		//Try for genome2
		while (it2 < end2)
		{
			//Inherits from genome 2
			offspring->history.push_back(genome2->history[it2]);
			offspring->sourceNode.push_back(genome2->sourceNode[it2]);
			offspring->destNode.push_back(genome2->destNode[it2]);
			offspring->weights.push_back(genome2->weights[it2]);

			//Inherits disabled links
			if (dis_it2 != dis_it2_end && *dis_it2 == it2)
			{
				offspring->disabledIndex.push_back(it2 + dis_offset_2);
				dis_it2++;
			}

			it2++;
		}

		//Finally, re-enable disabled links at x% rate
		offspring->disabledIndex.erase(
			std::remove_if(
				offspring->disabledIndex.begin(),
				offspring->disabledIndex.end(),
				[](uint16_t enableChance) {return (rng.RngProb()<enableChance); }),
			offspring->disabledIndex.end());


		//Collect the nodes because we forgot somehow...
		std::vector<N_SIZE> nodeCollector = offspring->sourceNode;
		nodeCollector.insert(nodeCollector.end(), offspring->destNode.begin(), offspring->destNode.end());

		std::vector<N_SIZE> inputNodeCollector = genome1->inputNode;
		inputNodeCollector.insert(inputNodeCollector.end(), genome2->inputNode.begin(), genome2->inputNode.end());

		std::vector<N_SIZE> outputNodeCollector = genome1->outputNode;
		outputNodeCollector.insert(outputNodeCollector.end(), genome2->outputNode.begin(), genome2->outputNode.end());

		//remove duplicates
		sort(nodeCollector.begin(), nodeCollector.end());
		nodeCollector.erase(unique(nodeCollector.begin(), nodeCollector.end()), nodeCollector.end());

		sort(inputNodeCollector.begin(), inputNodeCollector.end());
		inputNodeCollector.erase(unique(inputNodeCollector.begin(), inputNodeCollector.end()), inputNodeCollector.end());

		sort(outputNodeCollector.begin(), outputNodeCollector.end());
		outputNodeCollector.erase(unique(outputNodeCollector.begin(), outputNodeCollector.end()), outputNodeCollector.end());

		//Add to offspring
		offspring->inputNode = inputNodeCollector;
		offspring->outputNode = outputNodeCollector;
		offspring->nodes = nodeCollector.size();

		//Create the new Hash
		offspring->evolutionHash = Hash::HashGenetics(offspring->history);

		return offspring;
	}


	POP_PTR CreateOffsprings(POP_PTR pop, const std::vector<uint16_t>& census, const float noCrossover, const float enableGeneChance)
	{
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
		uint16_t nb_species = census.size();
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
			}
		}

		//Calculate how many offsprings per species we want
		uint16_t pop_size = pop->size() - nexxgen->size();
		std::vector<uint16_t> requiredChilds = std::vector<uint16_t>();
		int16_t spotsLeft = pop_size;
		for each (float speciesFit in speciesFitness)
		{
			uint16_t quantity = (int)(speciesFit / populationFitness * pop_size) + 1;
			requiredChilds.push_back(quantity);
			spotsLeft -= quantity;
		}
		uint16_t speciesCount = requiredChilds.size();
		//Fill leftover or remove excess at random
		if (spotsLeft < 0)
		{
			for (; spotsLeft != 0; spotsLeft++)
			{
				requiredChilds[rng.LessThan(speciesCount)] -= 1;
				//Species should DIE here!!
				requiredChilds.erase(std::remove(requiredChilds.begin(), requiredChilds.end(), 0), requiredChilds.end());
				speciesCount = requiredChilds.size();

			}
		}
		else if (spotsLeft > 0)
		{
			for (; spotsLeft != 0; spotsLeft--)
			{
				requiredChilds[rng.LessThan(speciesCount)] += 1;
			}
		}

		//Create offsprings for each species
		for (uint16_t species = 0; species < speciesCount; species++)
		{
			for (int childCount = 0; childCount < requiredChilds[species]; childCount++)
			{
				//Mutation without crossover
				if (rng.RngProb() < noCrossover)
				{
					//Select 1 genome at random based on fitness wheel
					auto left = pop->begin() + speciesIndex[species].first;
					float fitnessArrow = rng.RngProb() * speciesFitness[species];
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
					auto left = pop->begin() + speciesIndex[species].first;
					float fitnessArrow = rng.RngProb() * speciesFitness[species];
					while (fitnessArrow > (*left)->fitness)
					{
						fitnessArrow -= (*left)->fitness;
						left++;
					}
					GEN_PTR parent1 = std::make_unique<Genome>(**left);

					left = pop->begin() + speciesIndex[species].first;
					fitnessArrow = rng.RngProb() * speciesFitness[species];
					while (fitnessArrow > (*left)->fitness)
					{
						fitnessArrow -= (*left)->fitness;
						left++;
					}
					GEN_PTR parent2 = std::make_unique<Genome>(**left);

					//make a child
					nexxgen->push_back(Mate(parent1, parent2, enableGeneChance));
				}
			}
		}

		return nexxgen;
	}

}