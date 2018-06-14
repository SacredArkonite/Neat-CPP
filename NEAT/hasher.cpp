#include "hasher.h"

#include "hasher.h"
#include "defines.h"

#include <vector>

namespace Hash
{
	size_t HashGenetics(const std::vector<C_SIZE>& history) {
		std::hash<int> hasher;
		size_t newhash = 0;
		for (auto it = history.begin(); it != history.end(); it++) {
			newhash ^= hasher(*it) + 0x9e3779b9 + (newhash << 6) + (newhash >> 2);
		}
		return newhash;
	}

	size_t HashGenetics(const size_t oldHash, const std::vector<C_SIZE>& newHist) {
		std::hash<int> hasher;
		size_t newHash = oldHash;
		for (auto it = newHist.begin(); it != newHist.end(); it++)
		{
			newHash ^= hasher(*it) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
		}
		return newHash;
	}

	size_t HashGenetics(const size_t oldHash, const C_SIZE newHist) {
		std::hash<int> hasher;
		size_t newHash = oldHash;
		newHash ^= hasher(newHist) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
		return newHash;
	}

	size_t HashGenetics(const size_t oldHash, const C_SIZE newHist1, const C_SIZE newHist2) {
		std::hash<int> hasher;
		size_t newHash = oldHash;
		newHash ^= hasher(newHist1) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
		newHash ^= hasher(newHist2) + 0x9e3779b9 + (newHash << 6) + (newHash >> 2);
		return newHash;
	}
}