#ifndef DEDUP_CORE_H__
#define DEDUP_CORE_H__

namespace dashing2 {
extern int exhaustive_dedup;
extern int maxcand_global;

#ifndef EARLYSTOP
#define EARLYSTOP 1
#endif
static constexpr bool earlystop = EARLYSTOP;
size_t default_candidates(const size_t nitems);

}

#endif
