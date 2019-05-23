//
// Created by natalia on 10/05/19.
//

#ifndef INTPOINTSFEM_TIMER_H
#define INTPOINTSFEM_TIMER_H

#include <chrono>
using namespace std::chrono;

#ifdef __CUDACC__
#include <cuda.h>
#endif

class Timer {
protected:
#ifdef __CUDACC__
    cudaEvent_t     start;
    cudaEvent_t     stop;
    float			time_span;
#else
    high_resolution_clock::time_point 	start;
    high_resolution_clock::time_point 	stop;
    REAL			 					time_span;
#endif
    int time_unit;
public:
    enum WhichUnit{ENanoseconds, EMicroseconds, EMilliseconds, ESeconds};

    Timer() {
#ifdef __CUDACC__
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
#endif
    }

    void TimerConfig(WhichUnit unit) {
    	time_unit = unit;
    }

    ~Timer() {
#ifdef __CUDACC__
    	cudaEventDestroy(start);
		cudaEventDestroy(stop);
#endif
    }

    void Start() {
#ifdef __CUDACC__
    	cudaEventRecord(start);
#else
    	start = high_resolution_clock::now();
#endif
    }

    void Stop() {
#ifdef __CUDACC__
    	cudaEventRecord(stop);
    	cudaEventSynchronize(stop);
#else
    	stop = high_resolution_clock::now();
#endif
    }

    REAL ElapsedTime() {
#ifdef __CUDACC__
    	cudaEventElapsedTime(&time_span, start, stop);
    	if(time_unit == ENanoseconds) {
    		return time_span*1000000;
    	}
    	else if(time_unit == EMicroseconds) {
    		return time_span*1000;
    	}
    	else if(time_unit == EMilliseconds) {
    		return time_span*1;
    	}
    	else if(time_unit == ESeconds) {
    		return time_span/1000;
    	}
    	else {
    		return time_span; //default milliseconds
    	}

#else
        time_span = duration_cast<nanoseconds>(stop - start).count();
    	if(time_unit == ENanoseconds) {
    		return time_span;
    	}
    	else if(time_unit == EMicroseconds) {
    		return time_span/1000;
    	}
    	else if(time_unit == EMilliseconds) {
    		return time_span/1000000;
    	}
    	else if(time_unit == ESeconds) {
    		return time_span/1000000000;
    	}
    	else {
    		return time_span/1000000; //default milliseconds
    	}
#endif
    }

    std::string Unit() {
        if(time_unit == ENanoseconds) {
            return "    ns";
        }
        else if(time_unit == EMicroseconds) {
            return "    us";
        }
        else if(time_unit == EMilliseconds) {
            return "    ms";
        }
        else if(time_unit == ESeconds) {
            return "    s";
        }
        else {
            return "    ms"; //default milliseconds
        }
    }
};


#endif //INTPOINTSFEM_TIMER_H
