#include<iostream>
#include<vector>
#include<memory>
#include<algorithm>
#include<random>
#include<map>

#include"RNG.h"

RNG rng;

typedef uint8_t N_SIZE;
typedef uint16_t N_SIZEx2;
typedef uint32_t C_SIZE;
#define MAX_NODES 0xFF
#define MAX_CONNECTIONS 1500

struct Genome
{
	//Genome(const Genome & src) = default;
		/*sourceNode(src.sourceNode),
		destNode(src.destNode),
		history(src.history),
		disabledIndex(src.disabledIndex),
		weights(src.weights),
		nodes(src.nodes),
		evolutionHash(src.evolutionHash),
		species(src.species),
		fitness(src.fitness)
	{}*/

	std::vector<N_SIZE> sourceNode;
	std::vector<N_SIZE> destNode;
	std::vector<C_SIZE> history;
	std::vector<C_SIZE> disabledIndex;
	std::vector<float> weights;
	N_SIZE nodes = 0;
	size_t evolutionHash = 0;
	uint16_t species = 0;
	float fitness = 0.0f;

	bool operator < (const Genome& str) const
	{
		return (species < str.species);
	}
};
typedef std::unique_ptr<Genome> GEN_PTR;
typedef std::unique_ptr<std::vector<GEN_PTR>> POP_PTR;

size_t hashGenetics(const std::vector<C_SIZE>& history) {
	std::hash<int> hasher;
	size_t newhash = 0;
	for (auto it = history.begin(); it != history.end(); it++) {
		newhash ^= hasher(*it) + 0x9e3779b9 + (newhash << 6) + (newhash >> 2);
	}
	return newhash;
}

size_t hashGenetics(const size_t oldHash, const std::vector<C_SIZE>& newHist) {
	std::hash<int> hasher;
	size_t newHash = oldHash;
	for (auto it = newHist.begin(); it != newHist.end(); it++)
	{
		newHash ^= hasher(*it) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
	}
	return newHash;
}

size_t hashGenetics(const size_t oldHash, const C_SIZE newHist1, const C_SIZE newHist2) {
	std::hash<int> hasher;
	size_t newHash = oldHash;
	newHash ^= hasher(newHist1) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
	newHash ^= hasher(newHist2) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
	return newHash;
}

GEN_PTR CreateGenome(const N_SIZE nIns, const N_SIZE nOuts)
{
	auto genome = std::make_unique<Genome>();
	C_SIZE hist = 1;
	genome->nodes = nIns + nOuts;
	genome->species = 0;
	genome->fitness = 0.0;
	for (N_SIZE i = 0; i < nIns; i++) {
		for (N_SIZE j = nIns; j < genome->nodes; j++) {
			genome->sourceNode.push_back(i);
			genome->destNode.push_back(j);
			genome->weights.push_back(rng.RngWeight());
			genome->history.push_back(hist++);
		}
	}

	genome->evolutionHash = hashGenetics(genome->history);
	return genome;
}

GEN_PTR MutateAddConnection(GEN_PTR genome, const N_SIZE from, const N_SIZE to, const C_SIZE histNb)
{
	genome->history.push_back(histNb);
	genome->sourceNode.push_back(from);
	genome->destNode.push_back(to);
	genome->weights.push_back(rng.RngWeight());

	return genome;
}

GEN_PTR MutateAddNode(GEN_PTR genome, const C_SIZE index, const C_SIZE histNb)
{
	//Don't forget to disable the old connection!!
	genome->disabledIndex.push_back(index);
	genome->nodes++;
	genome->history.push_back(histNb);
	genome->history.push_back(histNb+1);

	genome->sourceNode.push_back(genome->sourceNode[index]);
	genome->destNode.push_back(genome->nodes);
	genome->weights.push_back(1.0);

	genome->sourceNode.push_back(genome->nodes);
	genome->destNode.push_back(genome->destNode[index]);
	genome->weights.push_back(genome->weights[index]);

	genome->evolutionHash = hashGenetics(genome->evolutionHash, histNb, histNb+1);

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

	return (c1*excess + c2*disjoint + c3*d_sum/matching);
}

POP_PTR GenerateExample(C_SIZE& hist)
{
	POP_PTR pop = CreatePop(2, 3, 1, hist);
	(*pop)[0] = MutateAddNode(std::move((*pop)[0]), 1, 4);
	(*pop)[1] = MutateAddNode(std::move((*pop)[1]), 1, 4);
	(*pop)[1] = MutateAddNode(std::move((*pop)[1]), 4, 6);
	(*pop)[0] = MutateAddConnection(std::move((*pop)[0]), 1, 5, 8);
	(*pop)[1] = MutateAddConnection(std::move((*pop)[1]), 3, 5, 9);
	(*pop)[1] = MutateAddConnection(std::move((*pop)[1]), 1, 6, 10);
	
	return pop;
}

bool CheckIfConnectionExists(const std::vector<GEN_PTR>::iterator it, const std::pair<N_SIZE, N_SIZE> connection)
{
	auto sr_begin = (*it)->sourceNode.begin();
	auto sr_ptr = sr_begin;
	auto sr_end = (*it)->sourceNode.end();
	C_SIZE index = -1 ; // We crash if not found
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
		for (last_key = it->first ; it != end; it++) {
			//Increment the mutation ID if when needed
			if (it->first != last_key)
			{
				last_key = it->first;
				hist += 2;
			}
			nexxgen->push_back(MutateAddNode(std::move(it->second),last_key.second,hist));
		}

		hist += 2;
	}

	//Iterate over the new conx multimap in key order
	if (!mm_conn.empty())
	{
		std::multimap<std::pair<size_t, std::pair<N_SIZE,N_SIZE>>, GEN_PTR>::iterator it = mm_conn.begin();
		std::multimap<std::pair<size_t, std::pair<N_SIZE, N_SIZE>>, GEN_PTR>::iterator end = mm_conn.end();
		std::pair<size_t, std::pair<N_SIZE, N_SIZE>> last_key = it->first;
		for (; it != end; it++) {
			//Increment the mutation ID if when needed
			if (it->first != last_key)
			{
				last_key = it->first;
				hist ++;
			}
			nexxgen->push_back(MutateAddConnection(std::move(it->second), last_key.second.first, last_key.second.second, hist));
		}
		hist++;
	}

	return nexxgen;
}

int debug_counter = 0;
POP_PTR ClassifyGenomes(POP_PTR encyclopedia, POP_PTR pop, std::vector<uint16_t>& census,  float c1, const float c2, const float c3, const float dt)
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
			debug_counter++;
			if (debug_counter == 141) {
				std::cout << "WTF \n";
			}
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
	
	std::sort(pop->begin(), pop->end());

	return pop;
}

POP_PTR UpdateEncyclopedia(const POP_PTR& pop, const std::vector<uint16_t>& census)
{
	//Select a random individual for each species
	POP_PTR encyclopedia = std::make_unique<std::vector<GEN_PTR>>();
	uint16_t specieLocation = 0;

	auto it_end = census.end();
	for (auto it = census.begin();it != it_end; it++)
	{
		if (*it != 0) {
			GEN_PTR specimen = std::make_unique<Genome>(*(pop->at(specieLocation + rng.LessThan(*it))));
			encyclopedia->push_back(std::move(specimen));
			specieLocation += *it;
		}
	}

	return encyclopedia;
}

int main()
{
	std::cout << "hello bitches!" << std::endl;

	//POP_PTR pop = GenerateExample();

	//float delta = Compatibility(1, 1, 1, *(*pop)[0], *(*pop)[1]);

	//Create initial population
	C_SIZE hist;
	POP_PTR pop = CreatePop(1000, 2, 1, hist);

	//Generate species dictionnary
	POP_PTR speciesEncyclopedia = std::make_unique<std::vector<GEN_PTR>>();;
	speciesEncyclopedia->push_back(std::make_unique<Genome>(*(pop->at(rng.LessThan(10)))));
	std::vector<uint16_t> census;

	for (int i = 0; i < 30; i++)
	{
		//Calculate Fitness / Simulate


		//Adjust Fitness


		//Speciate
		pop = ClassifyGenomes(std::move(speciesEncyclopedia), std::move(pop), census, 1, 1, 0, 4);

		//Update encyclopedia for next gen classification
		speciesEncyclopedia = UpdateEncyclopedia(pop, census);

		//Reproduce


		//Mutate Structure
		pop = MutateStructure(std::move(pop), 0.03, 0.05, hist);
		
		//Mutate weights
		pop = MutateWeights(std::move(pop), 0.8, 0.9, 0.5, 1.5);


		if (i%10 == 0) std::cout << "GEN # " << i << std::endl;
	}

	char w;
	std::cin>>w;
	return 0;
}