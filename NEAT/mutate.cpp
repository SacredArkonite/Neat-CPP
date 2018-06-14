#include "mutate.h"
#include "genome.h"
#include "defines.h"
#include "hasher.h"
#include "RNG.h"

namespace Mutate
{
	void AddConnection(const GEN_PTR& genome, const N_SIZE from, const N_SIZE to, const C_SIZE histNb)
	{
		//add only if the destination is not an input or the source an output
		if (std::find(genome->inputNode.begin(), genome->inputNode.end(), to) == genome->inputNode.end()
			&& std::find(genome->outputNode.begin(), genome->outputNode.end(), from) == genome->outputNode.end())
		{
			genome->history.push_back(histNb);
			genome->sourceNode.push_back(from);
			genome->destNode.push_back(to);
			genome->weights.push_back(RNG::RngWeight());
			genome->evolutionHash = Hash::HashGenetics(genome->evolutionHash, histNb);
		}
	}

	void AddNode(const GEN_PTR& genome, const C_SIZE index, const C_SIZE histNb)
	{
		//Don't forget to disable the old connection!!
		genome->disabledIndex.push_back(index);
		genome->history.push_back(histNb);
		genome->history.push_back(histNb + 1);

		genome->sourceNode.push_back(genome->sourceNode[index]);
		genome->destNode.push_back(genome->nodes);
		genome->weights.push_back(1.0);

		genome->sourceNode.push_back(genome->nodes);
		genome->destNode.push_back(genome->destNode[index]);
		genome->weights.push_back(genome->weights[index]);

		genome->evolutionHash = Hash::HashGenetics(genome->evolutionHash, histNb, histNb + 1);
		genome->nodes++;
	}

	void AddInput(const GEN_PTR& genome, C_SIZE& histNb)
	{
		genome->inputNode.push_back(genome->nodes);
		for (int i = genome->outputNode.size() - 1; i >= 0; i--) {
			AddConnection(genome, genome->nodes, genome->outputNode[i], histNb++);
		}
		genome->nodes++;
	}

	void AddOutput(const GEN_PTR& genome, C_SIZE& histNb)
	{
		genome->outputNode.push_back(genome->nodes);
		for (int i = genome->inputNode.size() - 1; i >= 0; i--) {
			AddConnection(genome, genome->inputNode[i], genome->nodes, histNb++);
		}
		genome->nodes++;
	}

}
