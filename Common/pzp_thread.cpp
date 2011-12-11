#include "pzp_thread.h"

#include<map>

#ifdef THREAD_NOTES

using namespace std;

struct pz_pthread_info_t
{
  pz_pthread_info_t() : tid(0), last_op(0), last_ct(0)
  {}

  int tid;
  int last_op;
  int last_ct;
  // Create time
  // Join time
};

map<int,pz_pthread_info_t> threads_info;
static bool pz_pthread_logging = false;
static unsigned max_sim_threads = 0;
static unsigned sim_threads = 0;

pthread_mutex_t pz_pthread_mutex;

void pzp_thread_log_start()
{
  pz_pthread_logging = true;
  max_sim_threads = 0;
  sim_threads = 0;
}

void pzp_thread_log_stop()
{
  pz_pthread_logging = false;
}

void pzp_thread_log_report(ostream& o)
{
  cout << "pzp_thread_log_report:" << endl;
  cout << "Max number of simultaneous threads = " << max_sim_threads << endl;
}


int pz_pthread_join(pthread_t thread, void **value_ptr, const char* fn,
		    const char* file, unsigned line)
{
  int ret = pthread_join(thread, value_ptr);

  if (pz_pthread_logging && (ret == 0)) {

    pthread_mutex_lock(&pz_pthread_mutex);

    sim_threads--;
    std::cout << (int*) thread << " -> "
	      << (int*) pthread_self() << "; Join thread at "
	      << fn << " (" << file << ":" << line << ")."
	      << " sim_threads = " << sim_threads << std::endl;

    pthread_mutex_unlock(&pz_pthread_mutex);

  }

  return ret;
}

int pz_pthread_create(pthread_t *thread, const pthread_attr_t *attr, void
		      *(*start_routine)(void *), void *arg, const char* fn,
		      const char* file, unsigned line)
{
  int ret = pthread_create(thread, attr, start_routine, arg);

  if (pz_pthread_logging && (ret == 0)) {

    pthread_mutex_lock(&pz_pthread_mutex);

    sim_threads++;
    if (sim_threads > max_sim_threads) max_sim_threads = sim_threads;

    std::cout << (int*) pthread_self() << " - "		 
	      << (int*) *thread   << "; Create thread at " 
	      << fn << " (" << file << ":" << line << ")."
	      << " sim_threads = " << sim_threads << std::endl;

    pthread_mutex_unlock(&pz_pthread_mutex);

  }

  return ret;
}
  
#endif
