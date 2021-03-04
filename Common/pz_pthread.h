#ifndef PZ_PTHREAD_H
#define PZ_PTHREAD_H

//#include<pthread.h>
#include<iostream> // cout, ostream
#include <pz_config.h>

//#define THREAD_NOTES
//#define THREAD_MUTEX_NOTES
//#define THREAD_COND_NOTES

#ifdef THREAD_NOTES

#pragma message ( "warning: Thread notes are on... (see pz_pthread.h)" )

int pz_pthread_create(pthread_t *thread, const pthread_attr_t *attr, void
		      *(*start_routine)(void *), void *arg, const char* fn,
		      const char* file, unsigned line);

int pz_pthread_join(pthread_t thread, void **value_ptr, 
		    const char* fn, const char* file, unsigned line);

#define PZ_PTHREAD_CREATE(thread,attr,routine,args,fn)	\
  pz_pthread_create(thread,attr,routine,args,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_JOIN(thread,val,fn)				 \
  pz_pthread_join(thread,val,fn,__FILE__,__LINE__)

#else

#define PZ_PTHREAD_CREATE(thread,attr,routine,args,fn)	\
  pthread_create(thread,attr,routine,args)

#define PZ_PTHREAD_JOIN(thread,val,fn)		\
  pthread_join(thread,val)

#endif

#ifdef THREAD_MUTEX_NOTES

#pragma message ( "warning: Thread mutex notes are on... (see pz_pthread.h)" )

int pz_pthread_mutex_init(pthread_mutex_t* mutex, const pthread_mutexattr_t* attr,
			  const char* fn, const char* file, unsigned line);
int pz_pthread_mutex_destroy(pthread_mutex_t *mutex,
			     const char* fn, const char* file, unsigned line);
int pz_pthread_mutex_lock(pthread_mutex_t *mutex,
			  const char* fn, const char* file, unsigned line);
int pz_pthread_mutex_unlock(pthread_mutex_t *mutex,
			    const char* fn, const char* file, unsigned line);
// 11 thread_mutex_destroy
// 13 thread_mutex_init
// 49 thread_mutex_lock
// 51 thread_mutex_unlock

int pz_pthread_mutex_init(pthread_mutex_t* mutex, const pthread_mutexattr_t* attr,
const char* fn, const char* file, unsigned line);


#define PZ_PTHREAD_MUTEX_INIT(mutex,attr,fn)    \
  pz_pthread_mutex_init(mutex,attr,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_MUTEX_DESTROY(mutex,fn)      \
  pz_pthread_mutex_destroy(mutex,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_MUTEX_LOCK(mutex,fn)	        \
  pz_pthread_mutex_lock(mutex,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_MUTEX_UNLOCK(mutex,fn)       \
  pz_pthread_mutex_unlock(mutex,fn,__FILE__,__LINE__)

#else

#define PZ_PTHREAD_MUTEX_INIT(mutex,attr,fn)    \
  pthread_mutex_init(mutex,attr)

#define PZ_PTHREAD_MUTEX_DESTROY(mutex,fn)      \
  pthread_mutex_destroy(mutex)

#define PZ_PTHREAD_MUTEX_LOCK(mutex,fn)	        \
  pthread_mutex_lock(mutex)

#define PZ_PTHREAD_MUTEX_UNLOCK(mutex,fn)       \
  pthread_mutex_unlock(mutex)

#endif


#ifdef THREAD_COND_NOTES

#pragma message ( "warning: Thread cond notes are on... (see pz_pthread.h)" )

int pz_pthread_cond_init(pthread_cond_t*  cond, const pthread_condattr_t*  attr,
			 const char* fn, const char* file, unsigned line);

int pz_pthread_cond_destroy(pthread_cond_t *cond,
			    const char* fn, const char* file, unsigned line);

int pz_pthread_cond_wait(pthread_cond_t*  cond, pthread_mutex_t* mutex,
			 const char* fn, const char* file, unsigned line);

int pz_pthread_cond_broadcast(pthread_cond_t *cond,
			      const char* fn, const char* file, unsigned line);

int pz_pthread_cond_signal(pthread_cond_t *cond,
			   const char* fn, const char* file, unsigned line);

#define PZ_PTHREAD_COND_INIT(cond,attr,fn)		\
  pz_pthread_cond_init(cond,attr,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_COND_DESTROY(cond,fn)		\
  pz_pthread_cond_destroy(cond,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_COND_WAIT(cond,mutex,fn)			\
  pz_pthread_cond_wait(cond,mutex,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_COND_SIGNAL(cond,fn)			\
  pz_pthread_cond_signal(cond,fn,__FILE__,__LINE__)

#define PZ_PTHREAD_COND_BROADCAST(cond,fn)			\
  pz_pthread_cond_broadcast(cond,fn,__FILE__,__LINE__)

#else

#define PZ_PTHREAD_COND_INIT(cond,attr,fn)	\
  pthread_cond_init(cond,attr)

#define PZ_PTHREAD_COND_DESTROY(cond,fn)	\
  pthread_cond_destroy(cond)

#define PZ_PTHREAD_COND_WAIT(cond,mutex,fn)	\
  pthread_cond_wait(cond,mutex)

#define PZ_PTHREAD_COND_SIGNAL(cond,fn)		\
  pthread_cond_signal(cond)

#define PZ_PTHREAD_COND_BROADCAST(cond,fn)	\
  pthread_cond_broadcast(cond)

#endif


void pz_pthread_log_start();
void pz_pthread_log_stop();
void pz_pthread_log_report(std::ostream& o);

#endif // PZ_PTHREAD_H
