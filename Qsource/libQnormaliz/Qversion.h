#ifndef QNMZ_VERSION_H
#define QNMZ_VERSION_H

#define QNMZ_VERSION_MAJOR  3
#define QNMZ_VERSION_MINOR  6
#define QNMZ_VERSION_PATCH  1
#define QNMZ_VERSION        3.6.1
#define QNMZ_RELEASE (QNMZ_VERSION_MAJOR * 10000 + QNMZ_VERSION_MINOR * 100 + QNMZ_VERSION_PATCH)

namespace libQnormaliz {
inline unsigned int getVersion()
{
    return QNMZ_RELEASE;
}

} //end namespace libQnormaliz

#endif // QNMZ_VERSION_H
