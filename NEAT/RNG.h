#pragma once
#include <random>

#define FIXED_SEED 0x42069420

namespace RNG
{

	class RNG
	{
	public:
		RNG();
		~RNG();
		bool RngBool();
		uint32_t LessThan(uint32_t value);
		float RngProb();
		float RngRange();
		float RngWeight();

	private:
		std::mt19937 rng;
		std::uniform_int_distribution<uint16_t> rngBool;
		std::uniform_int_distribution<uint32_t> rngMax;
		std::uniform_real_distribution<float> rngProb;
		std::uniform_real_distribution<float> rngRange;
		std::random_device rd;
	};

	bool RngBool();
	uint32_t LessThan(uint32_t value);
	float RngProb();
	float RngRange();
	float RngWeight();
}