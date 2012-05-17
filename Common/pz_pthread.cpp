#include "pz_pthread.h"

#include<map>

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

class pz_pthread_mutex_t
{
public: 
  pz_pthread_mutex_t() 
  {
    pthread_mutex_init(&mutex, NULL);
  }
  ~pz_pthread_mutex_t()
  {
    pthread_mutex_destroy(&mutex);
  }

  int lock() {
    int ret = pthread_mutex_lock(&mutex);
    if (ret) {
      cerr << "ERROR: (" << __FILE__ << ":" << __LINE__ << ") pthread_mutex_lock(&mutex) returned " << ret << endl;
    }
    return ret;
  }

  int unlock() {
    return pthread_mutex_unlock(&mutex);
  }

private:
  
  pthread_mutex_t mutex;

};

#if (defined THREAD_NOTES || defined THREAD_MUTEX_NOTES || defined THREAD_COND_NOTES)

map<int,pz_pthread_info_t> threads_info;
static bool pz_pthread_logging = false;
static unsigned max_sim_threads = 0;
static unsigned sim_threads = 0;
pz_pthread_mutex_t pz_pmutex;

void pz_pthread_log_start()
{
  pz_pthread_logging = true;

#ifdef THREAD_NOTES
  max_sim_threads = 0;
  sim_threads = 0;
#endif //THREAD_NOTES
}

void pz_pthread_log_stop()
{
  pz_pthread_logging = false;
}

void pz_pthread_log_report(ostream& o)
{
  cout << "pz_pthread_log_report:" << endl;
#ifdef THREAD_NOTES
  cout << "Max number of simultaneous threads = " << max_sim_threads << endl;
#endif //THREAD_NOTES
}

#else // THREAD_NOTES || THREAD_MUTEX_NOTES || THREAD_COND_NOTES

void pz_pthread_log_stop() {}
void pz_pthread_log_start() {}
void pz_pthread_log_report(std::ostream& o) {}

#endif // THREAD_NOTES || THREAD_MUTEX_NOTES || THREAD_COND_NOTES



#ifdef THREAD_NOTES

int pz_pthread_join(pthread_t thread, void **value_ptr, const char* fn,
		    const char* file, unsigned line)
{
  int ret = pthread_join(thread, value_ptr);

  if (pz_pthread_logging && (ret == 0)) {

    pz_pmutex.lock();

    sim_threads--;
    std::cout << (int*) thread << " -> "
	      << (int*) pthread_self() << "; Join thread at "
	      << fn << " (" << file << ":" << line << ")."
	      << " sim_threads = " << sim_threads << std::endl;

    pz_pmutex.unlock();

  }

  return ret;
}

int pz_pthread_create(pthread_t *thread, const pthread_attr_t *attr, void
		      *(*start_routine)(void *), void *arg, const char* fn,
		      const char* file, unsigned line)
{
  int ret = pthread_create(thread, attr, start_routine, arg);

  if (pz_pthread_logging && (ret == 0)) {

    pz_pmutex.lock();

    sim_threads++;
    if (sim_threads > max_sim_threads) max_sim_threads = sim_threads;

    std::cout << (int*) pthread_self() << " - "		 
	      << (int*) *thread   << "; Create thread at " 
	      << fn << " (" << file << ":" << line << ")."
	      << " sim_threads = " << sim_threads << std::endl;

    pz_pmutex.unlock();

  }

  return ret;
}

#endif //THREAD_NOTES


#ifdef THREAD_MUTEX_NOTES

int pz_pthread_mutex_init(pthread_mutex_t* mutex, const pthread_mutexattr_t* attr,
			  const char* fn, const char* file, unsigned line)
{
  int ret = pthread_mutex_init(mutex, attr);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_mutex_init() ret : " << ret 
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }

  return ret;
}

int pz_pthread_mutex_destroy(pthread_mutex_t *mutex,
			     const char* fn, const char* file, unsigned line)
{
  int ret = pthread_mutex_destroy(mutex);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_mutex_destroy() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

int pz_pthread_mutex_lock(pthread_mutex_t *mutex,
			  const char* fn, const char* file, unsigned line)
{
  int ret = pthread_mutex_lock(mutex);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_mutex_lock() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

int pz_pthread_mutex_unlock(pthread_mutex_t *mutex,
			    const char* fn, const char* file, unsigned line)
{
  int ret = pthread_mutex_unlock(mutex);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_mutex_unlock() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

#endif


#ifdef THREAD_COND_NOTES

int pz_pthread_cond_init(pthread_cond_t* cond, const pthread_condattr_t* attr,
			 const char* fn, const char* file, unsigned line)
{
  int ret = pthread_cond_init(cond,attr);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_cond_init() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

int pz_pthread_cond_destroy(pthread_cond_t *cond,
			    const char* fn, const char* file, unsigned line)
{
  int ret = pthread_cond_destroy(cond);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_cond_destroy() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

int pz_pthread_cond_wait(pthread_cond_t* cond, pthread_mutex_t* mutex,
			 const char* fn, const char* file, unsigned line)
{
  int ret = pthread_cond_wait(cond, mutex);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_cond_wait() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

int pz_pthread_cond_broadcast(pthread_cond_t *cond,
			      const char* fn, const char* file, unsigned line)
{
  int ret = pthread_cond_broadcast(cond);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_cond_broadcast() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

int pz_pthread_cond_signal(pthread_cond_t *cond,
			   const char* fn, const char* file, unsigned line)
{
  int ret = pthread_cond_signal(cond);

  if (pz_pthread_logging) {

    pz_pmutex.lock();

    std::cout << "pthread_cond_signal() ret : " << ret
	      << " : self thread " << (int*) pthread_self() << " : at : "		 
	      << fn << " : " << file << " : " << line << std::endl;
    
    pz_pmutex.unlock();

  }
  return ret;
}

#endif //THREAD_COND_NOTES
