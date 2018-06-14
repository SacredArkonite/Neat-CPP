#include "genome.h"

#include "mutate.h"
#include "hasher.h"
#include "RNG.h"

namespace GenomeUtil
{
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

	bool CheckIfConnectionExists(const GEN_PTR& gen, const std::pair<N_SIZE, N_SIZE> connection)
	{
		auto sr_begin = gen->sourceNode.begin();
		auto sr_ptr = sr_begin;
		auto sr_end = gen->sourceNode.end();
		C_SIZE index = -1; // We crash if not found
		for (; sr_ptr != sr_end; sr_ptr++) {
			if ((*sr_ptr) == connection.first)
			{
				index = sr_ptr - gen->sourceNode.begin();
				if (gen->destNode[index] == connection.second)
				{	//FOUND IT!
					//Remove from the disabled list if its there
					auto pos_in_disabled = std::find(gen->disabledIndex.begin(), gen->disabledIndex.end(), index);
					if (pos_in_disabled != gen->disabledIndex.end())
						gen->disabledIndex.erase(pos_in_disabled);
					return true;
				}
			}
		}
		return false;
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
				offspring->weights.push_back(RNG::RngBool() ? genome1->weights[it1] : genome2->weights[it2]);

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
				[](uint16_t enableChance) {return (RNG::RngProb()<enableChance); }),
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

}