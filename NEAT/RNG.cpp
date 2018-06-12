#include "RNG.h"


RNG::RNG() :rng(FIXED_SEED), rngBool(0, 1), rngMax(0, 0xFFFFFFFF), rngProb(0, 1), rngRange(0.5,1.5)
{
}


RNG::~RNG()
{
}

bool RNG::RngBool()
{
	return rngBool(rng);
}

uint32_t RNG::LessThan(uint32_t value)
{
	return rngMax(rng) % value;
}

float RNG::RngProb()
{
	return rngProb(rng);
}

float RNG::RngRange()
{
	return rngRange(rng);
}

float RNG::RngWeight()
{
	return rngBool(rng) ? rngRange(rng) : -rngRange(rng);
}

