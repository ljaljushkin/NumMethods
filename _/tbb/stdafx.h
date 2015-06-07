#include <math.h>
#include <time.h>
#include "timer.h"
#include <stdlib.h>
#include <stdio.h>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/compat/thread"
using namespace tbb;

#include "WawePropagation.h"
#include "CalcError.h"