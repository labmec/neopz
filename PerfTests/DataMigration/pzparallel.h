#include <unistd.h>
#include <thread>
#include <vector>
namespace pz {
    // util - return the number of cores on linux and macosx
    int get_number_of_cores() {
        return sysconf(_SC_NPROCESSORS_ONLN);
    }
    
    template <typename body_t>
    struct thread_arg {
        int i;
        body_t *obj;
    };
    
    template <typename body_t>
    static void *thread_work(void *parm) {
        thread_arg<body_t> *arg = (thread_arg<body_t>*) parm;
        (*arg->obj)(arg->i);
        return NULL;
    }
    
    template <typename body_t>
    void parallel_for(int n, body_t &obj) {
        // get the number of cores
        int nthreads = get_number_of_cores();
        // case the number of threads is lower than the cores
        if(nthreads > n)
            nthreads = n;
        // allocate threads
        std::vector<std::thread> threads;
        thread_arg<body_t> *args = new thread_arg<body_t>[nthreads];
        // create args and threads
        for (int t = 0; t < nthreads; t++) {
            args[t].i = t;
            args[t].obj = &obj;
            threads.push_back(std::thread(thread_work<body_t>,&args[t]));
        }
        // syncronize all threads
        for (int t = 0; t < nthreads; t++)
          {threads[t].join();}
        // free memory
        delete [] args;
    }
    
} // namespace
