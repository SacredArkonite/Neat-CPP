#include "phenome.h"


namespace Phenome
{

	inline float ActivationFunction(const float value)
	{
		return value / (1 + abs(value));
	}

	std::vector<float> Propagate(const std::vector<float>& input, const GEN_PTR& gen, uint32_t maxSteps)
	{
		//Fill inputs
		std::vector<float> nodes(gen->nodes, 0);
		std::vector<float> temp_nodes(gen->nodes, 0);
		N_SIZE input_it = 0;
		N_SIZE input_it_end = gen->inputNode.size();
		

		//Let propagate a few steps
		auto disabled_it = gen->disabledIndex.begin();
		auto disabled_it_end = gen->disabledIndex.end();
		for (unsigned int i = 0; i < maxSteps; i++) {
			//Feed inputs
			input_it = 0;
			for (; input_it < input_it_end; input_it++)
			{
				N_SIZE inputName = gen->inputNode[input_it];
				nodes[inputName] = input[input_it];
			}
			//1 step
			for (unsigned int c = 0; c < gen->history.size(); c++) {
				if (disabled_it != disabled_it_end && *disabled_it == c)
					disabled_it++;
				else
					temp_nodes[gen->destNode[c]] += (nodes[gen->sourceNode[c]] * gen->weights[c]);
			}
			//Activation function
			for (unsigned int c = 0; c < nodes.size(); c++) {
				nodes[c] = ActivationFunction(temp_nodes[c]);
				temp_nodes[c] = nodes[c];
			}
		}

		//Collect the outputs
		N_SIZE output_it = 0;
		N_SIZE output_it_end = gen->outputNode.size();

		std::vector<float> output(gen->outputNode.size(), 0);
		for (; output_it < output_it_end; output_it++)
		{
			N_SIZE outputName = gen->outputNode[output_it];
			output[output_it] = nodes[outputName];
		}

		return output;
	}
}