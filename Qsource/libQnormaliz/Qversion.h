#ifndef QNMZ_VERSION_H
#define QNMZ_VERSION_H

#define QNMZ_VERSION_MAJOR  3
#define QNMZ_VERSION_MINOR  5
#define QNMZ_VERSION_PATCH  4
#define QNMZ_VERSION        3.5.4
#define QNMZ_RELEASE (QNMZ_VERSION_MAJOR * 10000 + QNMZ_VERSION_MINOR * 100 + QNMZ_VERSION_PATCH)

namespace libQnormaliz {
inline unsigned int getVersion()
{
    return QNMZ_RELEASE;
}

} //end namespace libQnormaliz

#endif // QNMZ_VERSION_H
