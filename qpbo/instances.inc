#ifndef _INSTANCES_INC_
#define _INSTANCES_INC_

#include "QPBO.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations

template class QPBO<int>;
template class QPBO<float>;
template class QPBO<double>;
template class QPBO<int64_t>;

#endif
