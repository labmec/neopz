#ifndef PZP_THREAD_H

#include<pthread.h>
#include<iostream> // cout, ostream

//#define THREAD_NOTES

#ifdef THREAD_NOTES

#warning "Thread notes are on... (see pzp_thread.h)"

int pz_pthread_create(pthread_t *thread, const pthread_attr_t *attr, void
		      *(*start_routine)(void *), void *arg, const char* fn,
		      const char* file, unsigned line);

int pz_pthread_join(pthread_t thread, void **value_ptr, 
		    const char* fn, const char* file, unsigned line);

#define PZP_THREAD_CREATE(thread,attr,routine,args,fn)	\
  pz_pthread_create(thread,attr,routine,args,fn,__FILE__,__LINE__)

#define PZP_THREAD_JOIN(thread,val,fn)				 \
  pz_pthread_join(thread,val,fn,__FILE__,__LINE__)

#else

#define PZP_THREAD_CREATE(thread,attr,routine,args,fn)	\
  pthread_create(thread,attr,routine,args)

#define PZP_THREAD_JOIN(thread,val,fn)		\
  pthread_join(thread,val)

#endif

void pzp_thread_log_start();
void pzp_thread_log_stop();
void pzp_thread_log_report(std::ostream& o);

#endif // PZP_THREAD_H
