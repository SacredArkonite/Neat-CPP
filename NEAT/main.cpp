#include<iostream>
#include<vector>
#include<memory>
#include<algorithm>
#include<random>
#include<map>

typedef uint8_t N_SIZE;
typedef uint16_t N_SIZEx2;
typedef uint32_t C_SIZE;
#define MAX_NODES 0xFF
#define MAX_CONNECTIONS 1500

struct Genome
{
	std::vector<N_SIZE> sourceNode;
	std::vector<N_SIZE> destNode;
	std::vector<C_SIZE> history;
	std::vector<C_SIZE> disabledIndex;
	std::vector<float> weights;
	N_SIZE nodes = 0;
	size_t evolutionHash = 0;
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
	for (N_SIZE i = 0; i < nIns; i++) {
		for (N_SIZE j = nIns; j < genome->nodes; j++) {
			genome->sourceNode.push_back(i);
			genome->destNode.push_back(j);
			genome->weights.push_back(0.0);
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
	genome->weights.push_back(0.0);

	return genome;
}

GEN_PTR MutateAddNode(GEN_PTR genome, const C_SIZE index, const C_SIZE histNb)
{
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
	excess = end1 - it1 + end2 - it2;

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
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<float> rngProb(0, 1);
	std::uniform_real_distribution<float> rngRange(minWeightMutation, maxWeightMutation);
	std::uniform_int_distribution<unsigned short> rngBool(0, 1);

	auto pop_end = pop->end();
	for (auto it = pop->begin(); it != pop_end; it++) {
		//Mutate the genome?
		if (rngProb(rng) < genomeWeightMutation)
		{
			auto i_w_end = (*it)->weights.end();
			for (auto i_w = (*it)->weights.begin(); i_w != i_w_end; i_w++)
			{
				//Mutate the weight by scaling?
				if (rngProb(rng) < weightMutation)
				{
					(*i_w) *= rngRange(rng);
				}
				//Mutate random
				else
				{
					// change sing?
					if (rngBool(rng))
					{
						(*i_w) = rngRange(rng);
					}
					else
					{
						(*i_w) = -rngRange(rng);
					}
				}
			}
		}
	}

	return pop;
}

POP_PTR MutatePop(POP_PTR pop, float new_node_percent, float new_link_percent, float genomeWeightMutation, float weightMutation, float minWeightMutation, float maxWeightMutation, C_SIZE& hist)
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

	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<unsigned int> rngDist(0, 0xFFFFFFFF);
	auto end = (*pop).end();

	//Distribute the genome with their transform into the maps
	for (auto it = (*pop).begin(); it < end; it++)
	{
		float rn = (rngDist(rng) % 256) / 256.0f;
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
				C_SIZE connection_id = rngDist(rng) % ((*it)->history.size());
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
				auto connection = std::pair<N_SIZE, N_SIZE>(rngDist(rng) % nb_nodes, rngDist(rng) % nb_nodes);

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

	//Mutate weights
	nexxgen = MutateWeights(std::move(nexxgen), genomeWeightMutation,  weightMutation,  minWeightMutation,maxWeightMutation);

	return nexxgen;
}

int main()
{
	std::cout << "hello bitches!" << std::endl;

	//POP_PTR pop = GenerateExample();

	//float delta = Compatibility(1, 1, 1, *(*pop)[0], *(*pop)[1]);

	C_SIZE hist;
	POP_PTR pop = CreatePop(100, 2, 1, hist);
	for (int i = 0; i < 20; i++)
	{
		pop = MutatePop(std::move(pop), 0.25, 0.25, 0.8, 0.9, 0.5, 1.5, hist);
		if (i%10 == 0) std::cout << "GEN # " << i << std::endl;
	}
	//size_t type0 = pop->at(0)->evolutionHash;
	//size_t type1 = pop->at(1)->evolutionHash;
	//size_t type2 = pop->at(2)->evolutionHash;
	//size_t type3 = pop->at(3)->evolutionHash;
	//size_t type4 = pop->at(4)->evolutionHash;

//	std::sort(pop.begin(), pop.end(), [](std::unique_ptr<Genome> a, std::unique_ptr<Genome> b) { return a->historic < b->historic; });


	char w;
	std::cin>>w;
	return 0;
}