//
// Created by natalia on 10/05/19.
//

#ifndef INTPOINTSFEM_TIMER_H
#define INTPOINTSFEM_TIMER_H

#include <chrono>
using namespace std::chrono;

#ifdef USING_CUDA
#include <cuda.h>
#endif

#ifdef USING_TBB
#include "tbb/parallel_for.h"
#include "tbb/tick_count.h"
#endif

class Timer {

protected:
#ifdef USING_CUDA
    cudaEvent_t                         start_cuda;
    cudaEvent_t                         stop_cuda;
    float			                    time_span_cuda;
#endif

#ifdef USING_TBB
    tbb::tick_count                     start_tbb;
    tbb::tick_count                     stop_tbb;
    REAL                                time_span_tbb;
#endif

    high_resolution_clock::time_point 	start_chrono;
    high_resolution_clock::time_point 	stop_chrono;
    REAL			 					time_span_chrono;

    int time_unit;
    int timer_option;

public:
    enum WhichUnit{ENanoseconds, EMicroseconds, EMilliseconds, ESeconds};
    enum WhichTimer{EChrono, ECudaEvent, ETBBTimer};

    Timer() {
#ifdef USING_CUDA
        cudaEventCreate(&start_cuda);
        cudaEventCreate(&stop_cuda);
        time_span_cuda = 0;
#endif
        //default
        time_unit = ESeconds;
        timer_option = EChrono;
    }

    void TimeUnit(WhichUnit unit) {
    	time_unit = unit;
    }

    void TimerOption(WhichTimer option) {
        timer_option = option;
    }

    ~Timer() {
#ifdef USING_CUDA
        if (timer_option == ECudaEvent) {
            cudaEventDestroy(start_cuda);
            cudaEventDestroy(stop_cuda);
        }
#endif
    }

    void Start() {
#ifdef USING_CUDA
    	if (timer_option == ECudaEvent) cudaEventRecord(start_cuda);
#endif

#ifdef USING_TBB
        if (timer_option == ETBBTimer) start_tbb = tbb::tick_count::now();
#endif

    	if (timer_option == EChrono) start_chrono = high_resolution_clock::now();
    }

    void Stop() {
#ifdef USING_CUDA
        if (timer_option == ECudaEvent) {
            cudaEventRecord(stop_cuda);
            cudaEventSynchronize(stop_cuda);
        }
#endif

#ifdef USING_TBB
        if (timer_option == ETBBTimer) stop_tbb = tbb::tick_count::now();
#endif

    	if (timer_option == EChrono) stop_chrono = high_resolution_clock::now();

    }

    REAL ElapsedTime() {
#ifdef USING_CUDA
    if (timer_option == ECudaEvent) {
        cudaEventElapsedTime(&time_span_cuda, start_cuda, stop_cuda);
        if(time_unit == ENanoseconds) {
            return time_span_cuda*1000000;
        }
        else if(time_unit == EMicroseconds) {
            return time_span_cuda*1000;
        }
        else if(time_unit == EMilliseconds) {
            return time_span_cuda*1;
        }
        else if(time_unit == ESeconds) {
            return time_span_cuda/1000;
        }
    }
#endif

#ifdef USING_TBB
    if (timer_option == ETBBTimer) {
        time_span_tbb = (stop_tbb - start_tbb).seconds();
        if(time_unit == ENanoseconds) {
            return time_span_tbb*1000000000;
        }
        else if(time_unit == EMicroseconds) {
            return time_span_tbb*1000000;
        }
        else if(time_unit == EMilliseconds) {
            return time_span_tbb*1000;
        }
        else if(time_unit == ESeconds) {
            return time_span_tbb;
        }
    }
#endif

    if (timer_option == EChrono) {
        time_span_chrono = duration_cast<nanoseconds>(stop_chrono - start_chrono).count();
    	if(time_unit == ENanoseconds) {
    		return time_span_chrono;
    	}
    	else if(time_unit == EMicroseconds) {
    		return time_span_chrono/1000;
    	}
    	else if(time_unit == EMilliseconds) {
    		return time_span_chrono/1000000;
    	}
    	else if(time_unit == ESeconds) {
    		return time_span_chrono/1000000000;
    	}
    }
}

    std::string Unit() {
        if(time_unit == ENanoseconds) {
            return "\tns";
        }
        else if(time_unit == EMicroseconds) {
            return "\tus";
        }
        else if(time_unit == EMilliseconds) {
            return "\tms";
        }
        else if(time_unit == ESeconds) {
            return "\ts";
        }
    }
};


#endif //INTPOINTSFEM_TIMER_H
