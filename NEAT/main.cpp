#include<iostream>
#include<vector>
#include<memory>
#include<algorithm>

struct Genome
{
	std::vector<int> sourceGene;
	std::vector<int> destGene;
	std::vector<int> history;
	std::vector<int> disabledIndex;
	std::vector<float> weights;
	int nodes = 0;
	size_t evolutionHash = 0;
};
typedef std::unique_ptr<std::vector<std::unique_ptr<Genome>>> POP_PTR;

size_t hashGenetics(const std::vector<int>& history) {
	std::hash<int> hasher;
	size_t newhash = 0;
	for (auto it = history.begin(); it != history.end(); it++) {
		newhash ^= hasher(*it) + 0x9e3779b9 + (newhash << 6) + (newhash >> 2);
	}
	return newhash;
}

size_t hashGenetics(const size_t oldHash, const std::vector<int>& newHist) {
	std::hash<int> hasher;
	size_t newHash = oldHash;
	for (auto it = newHist.begin(); it != newHist.end(); it++)
	{
		newHash ^= hasher(*it) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
	}
	return newHash;
}

size_t hashGenetics(const size_t oldHash, const int newHist1, const int newHist2) {
	std::hash<int> hasher;
	size_t newHash = oldHash;
	newHash ^= hasher(newHist1) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
	newHash ^= hasher(newHist2) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
	return newHash;
}

std::unique_ptr<Genome> CreateGenome(const int nIns, const int nOuts)
{
	auto genome = std::make_unique<Genome>();
	int hist = 1;
	genome->nodes = nIns + nOuts;
	for (int i = 1; i <= nIns; i++) {
		for (int j = nIns + 1; j <= genome->nodes; j++) {
			genome->sourceGene.push_back(i);
			genome->destGene.push_back(j);
			genome->weights.push_back(0.0);
			genome->history.push_back(hist++);
		}
	}

	genome->evolutionHash = hashGenetics(genome->history);
	return genome;
}

std::unique_ptr<Genome> MutateAddConnection(std::unique_ptr<Genome> genome, const int from, const int to, const int histNb)
{
	genome->history.push_back(histNb);
	genome->sourceGene.push_back(from);
	genome->destGene.push_back(to);
	genome->weights.push_back(0.0);

	return genome;
}

std::unique_ptr<Genome> MutateAddNode(std::unique_ptr<Genome> genome, const int index, const int histNb1, const int histNb2)
{
	genome->disabledIndex.push_back(index);
	genome->nodes++;
	genome->history.push_back(histNb1);
	genome->history.push_back(histNb2);

	genome->sourceGene.push_back(genome->sourceGene[index]);
	genome->destGene.push_back(genome->nodes);
	genome->weights.push_back(1.0);

	genome->sourceGene.push_back(genome->nodes);
	genome->destGene.push_back(genome->destGene[index]);
	genome->weights.push_back(genome->weights[index]);

	genome->evolutionHash = hashGenetics(genome->evolutionHash, histNb1, histNb2);

	return genome;
}

POP_PTR CreatePop(const int popSize, const int inputs, const int outputs)
{
	POP_PTR pop = std::make_unique<std::vector<std::unique_ptr<Genome>>>();

	for (int i = 0; i < popSize; i++)
	{
		pop->push_back(CreateGenome(inputs, outputs));
	}

	return pop;
}

void MutatePop(std::vector<std::unique_ptr<Genome>>& pop)
{
	pop[2] = MutateAddNode(std::move(pop[2]), 2, 4, 5);
	pop[3] = MutateAddNode(std::move(pop[3]), 2, 4, 5);
	pop[4] = MutateAddNode(std::move(pop[4]), 0, 6, 7);
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

POP_PTR GenerateExample()
{
	POP_PTR pop = CreatePop(2, 3, 1);
	(*pop)[0] = MutateAddNode(std::move((*pop)[0]), 1, 4, 5);
	(*pop)[1] = MutateAddNode(std::move((*pop)[1]), 1, 4, 5);
	(*pop)[1] = MutateAddNode(std::move((*pop)[1]), 4, 6, 7);
	(*pop)[0] = MutateAddConnection(std::move((*pop)[0]), 1, 5, 8);
	(*pop)[1] = MutateAddConnection(std::move((*pop)[1]), 3, 5, 9);
	(*pop)[1] = MutateAddConnection(std::move((*pop)[1]), 1, 6, 10);
	
	return pop;
}

int main()
{
	std::cout << "hello bitches!" << std::endl;

	POP_PTR pop = GenerateExample();

	float delta = Compatibility(1, 1, 1, *(*pop)[0], *(*pop)[1]);
	
	//POP_PTR pop = CreatePop(100, 3, 1);
	//MutatePop(*pop);
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