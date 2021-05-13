
#ifndef LIBNORMALIZ_NMZ_HASHLIBRARY_HPP
#define LIBNORMALIZ_NMZ_HASHLIBRARY_HPP

// We are using "hash-library" by Stephan Brumme
// to compute SHA-256 hashes.
// See
// https://create.stephan-brumme.com/hash-library/
// https://github.com/stbrumme/hash-library

#include <string>
#include <vector>

namespace libnormaliz {

using namespace std;

// Compute SHA-256 hash of a string, excluding final zero
// Return hex-values as string over '0',...,'9','a',...,'f'

string sha256str(const string& text, bool verbose = false);

// Compute SHA-256 hash of a string, excluding final zero
// Return as vector<char> of size sha256.HashBytes (== 32)
vector<unsigned char> sha256hexvec(const string& text, bool verbose = false);

}  // namespace libnormaliz
#endif

