#include "libsais/include/libsais64.h"


int64_t libsais64(const uint8_t * T, int64_t * SA, int64_t n, int64_t fs, int64_t * freq);

int64_t libsais64_plcp(const uint8_t * T, const int64_t * SA, int64_t * PLCP, int64_t n);

int64_t libsais64_lcp(const int64_t * PLCP, const int64_t * SA, int64_t * LCP, int64_t n);