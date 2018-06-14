#pragma once

#include "genome.h"
#include "defines.h"

namespace Mutate
{
	GEN_PTR AddConnection(GEN_PTR genome, const N_SIZE from, const N_SIZE to, const C_SIZE histNb);
	GEN_PTR AddNode(GEN_PTR genome, const C_SIZE index, const C_SIZE histNb);
	GEN_PTR AddInput(GEN_PTR genome, C_SIZE& histNb);
	GEN_PTR AddOutput(GEN_PTR genome, C_SIZE& histNb);

}