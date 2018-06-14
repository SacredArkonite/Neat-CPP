#pragma once

#include"defines.h"

#include <vector>

namespace Hash
{
	size_t HashGenetics(const std::vector<C_SIZE>& history);
	size_t HashGenetics(const size_t oldHash, const std::vector<C_SIZE>& newHist);
	size_t HashGenetics(const size_t oldHash, const C_SIZE newHist);
	size_t HashGenetics(const size_t oldHash, const C_SIZE newHist1, const C_SIZE newHist2);
}