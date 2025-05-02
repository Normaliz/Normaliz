#ifndef LIBNORMALIZ_NMZ_VERSION_H
#define LIBNORMALIZ_NMZ_VERSION_H

#define NMZ_VERSION_MAJOR  3
#define NMZ_VERSION_MINOR  10
#define NMZ_VERSION_PATCH  5
#define NMZ_VERSION        3.10.5
#define NMZ_RELEASE (NMZ_VERSION_MAJOR * 10000 + NMZ_VERSION_MINOR * 100 + NMZ_VERSION_PATCH)

namespace libnormaliz {
unsigned int getVersion();
} //end namespace libnormaliz

#endif // NMZ_VERSION_H
