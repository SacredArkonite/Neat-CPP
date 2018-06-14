#pragma once

#include "genome.h"
#include "defines.h"

namespace Mutate
{
	void AddConnection(const GEN_PTR& genome, const N_SIZE from, const N_SIZE to, const C_SIZE histNb);
	void AddNode(const GEN_PTR& genome, const C_SIZE index, const C_SIZE histNb);
	void AddInput(const GEN_PTR& genome, C_SIZE& histNb);
	void AddOutput(const GEN_PTR& genome, C_SIZE& histNb);

}