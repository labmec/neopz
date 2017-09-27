
#include "TPZHash.h"

uint32_t Hash(std::string str) {

    std::string::size_type size = str.size();

    uint8_t *arr = new uint8_t[size];

    for (unsigned int i = 0; i < size; ++i) {
		arr[i] = str.at(i);
    }

    uint32_t out;

    MurmurHash3_x86_32(arr, size, 0xB0F57EE3, &out);

	delete[] arr;

    return out;
}