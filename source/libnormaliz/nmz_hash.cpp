#include <hash-library/sha256.h>
#include "libnormaliz/nmz_hash.h"
// hash-library by Stephan Brumme
// https://create.stephan-brumme.com/hash-library/
// https://github.com/stbrumme/hash-library

namespace libnormaliz {

using namespace std;

// Compute SHA-256 hash of a string, excluding final zero
// Return hex-values as string over '0',...,'9','a',...,'f'
string sha256str(const string& text) {
    SHA256 sha256;
    return sha256(text);
}

// Compute SHA-256 hash of a string, excluding final zero
// Return as vector<char> of size sha256.HashBytes (== 32)
vector<unsigned char> sha256hexvec(const string& text) {
    SHA256 sha256;
    sha256.add(text.c_str(), text.size());
    unsigned char rawHash[sha256.HashBytes];
    sha256.getHash(rawHash);
    int s = sizeof(rawHash)/sizeof(rawHash[0]);
    // assuming sizeof(unsigned char) == 1, then s == sha256.HashBytes == 32
    vector<unsigned char> v(rawHash, rawHash + s);
    return v;
}
}  // namespace libnormaliz
