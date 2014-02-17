/***************************************************************************
 *   Copyright (C) 2014 by Edson Borin                                     *
 *   edson@ic.unicamp.br                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "pznuma.h"

// TODO: replace TPZ_HWLOC_Manager by TPZ_PageMigration_Manager. 
//       allow mbing and move_pages to work even when hwloc is not linked (test for linux?)

#ifdef USING_HWLOC

#include <hwloc.h>
#include <string.h>
#include <errno.h>
#include <numaif.h>
//#include <getcpu.h> // getcpu

//#include <sys/syscall.h>

const long PAGE_EXP  = 12;
const long PAGE_SZ   = (1<<12);
const long PAGE_MSK  = ((1<<PAGE_EXP)-1);

#include "arglib.h"

clarg::argBool mig_mp("-mig_mp", "Use move_pages when migrating pages.");
clarg::argBool mig_mbind("-mig_mbind", "Use mbind when migrating pages.");
clarg::argBool mig_hwloc("-mig_hwloc", "Use hwloc when migrating pages.");

class TPZ_HWLOC_Manager
{
public:

  TPZ_HWLOC_Manager()
  {
    hwloc_topology_init(&hw_topo);
    hwloc_topology_load(hw_topo);

    hwloc_obj_t obj;
    hw_cache_sz = 0;
    for (obj = hwloc_get_obj_by_type(hw_topo, HWLOC_OBJ_PU, 0); obj; obj = obj->parent)
      if (obj->type == HWLOC_OBJ_CACHE)
	hw_cache_sz += obj->attr->cache.size;
  }

  ~TPZ_HWLOC_Manager()
  {
    hwloc_topology_destroy(hw_topo);
  }

  void migrate_to_local(char* start, unsigned long long sz_in_bytes) {
    if (mig_mp.was_set()) {
      migrate_to_local_move_pages(start, sz_in_bytes);
    }
    else if (mig_mbind.was_set()) {
      migrate_to_local_mbind(start, sz_in_bytes);
    }
    else if (mig_hwloc.was_set()) {
      migrate_to_local_hwloc(start, sz_in_bytes);
    }
  }

  void migrate_to_local_move_pages(char* start, unsigned long long sz_in_bytes) {
    int ret;
    unsigned cpu, node;
    if ( (ret=sched_getcpu()) < 0) {
      printf("sched_getcpu(...) returned %d\n", ret);
      return;
    }
    else {
      cpu = ret;
      node = cpu >> 3;
    }

    unsigned long nodemask = 1 << node;
    unsigned long maxnode = 9;
    
    long page_start = (long) start;
    page_start = page_start & ~PAGE_MSK;
    long last_page = (long) start + (long) sz_in_bytes;
    last_page = last_page & ~PAGE_MSK;
    long count = (last_page - page_start + PAGE_SZ) >> PAGE_EXP;

    void** pages  = (void**) malloc(count*sizeof(void*));
    int*   nodes  = (int*)   malloc(count*sizeof(int));
    int*   status = (int*)   malloc(count*sizeof(int));
    
    for (int i=0; i<count; i++) {
      pages[i] = (void*) (page_start + PAGE_SZ*i);
      nodes[i] = node;
      status[i] = 0;
    }

    if ( (ret = move_pages(0 /*pid*/, count , pages, nodes, status, MPOL_MF_MOVE)) != 0) {
      int err = errno;
      fprintf(stderr, "move_pages(0, npages(%d), pages, nodes (%d), status, MPOL_MF_MOVE) = %d (cpu=%d, node=%d)\n",
	      count, node, ret, cpu, node);
      fprintf(stderr, "errno (%d) = %s\n", err, strerror(err));
    }
  }

  void migrate_to_local_mbind(char* start, unsigned long sz_in_bytes) {

    int ret;

    unsigned cpu, node;
    if ( (ret=sched_getcpu()) < 0) {
      printf("sched_getcpu(...) returned %d\n", ret);
      return;
    }
    else {
      cpu = ret;
      node = cpu >> 3;
    }

    unsigned long nodemask = 1 << node;
    //unsigned long maxnode = 7;
    //unsigned long nodemask = 1;
    unsigned long maxnode = 9;

    long page_start = (long) start;
    page_start = page_start & ~PAGE_MSK;
    long last_page = (long) start + (long) sz_in_bytes;
    last_page = last_page & ~PAGE_MSK;
    long nbytes = last_page - page_start + PAGE_SZ;

#ifndef MPOL_F_STATIC_NODES
#define MPOL_F_STATIC_NODES  (1 << 15)
#endif

    if ( (ret = mbind((void*) page_start, nbytes, MPOL_BIND  , &nodemask, maxnode, /* flags */ MPOL_MF_MOVE) ) != 0 ) {
      // Error.
      int err = errno;
      fprintf(stderr,"ERROR: mbind(start=%x, sz=%ld, MPOL_BIND, nodemask = %08x (node = %d), maxnode = %d, MPOL_MF_MOVE) = %d\n",
	      page_start, nbytes, nodemask, node, maxnode, ret);
      fprintf(stderr, "errno (%d) = %s\n", err, strerror(err));
      
      return;
    }
    else {
      fprintf(stderr,"mbind(start=%x, sz=%ld, MPOL_BIND, nodemask = %08x (node = %d), maxnode = %d, MPOL_MF_MOVE) = %d\n",
	      page_start, nbytes, nodemask, node, maxnode, ret);

    }
  }  

  void migrate_to_local_hwloc(char* start, unsigned long long sz_in_bytes) {

    int ret;

    // hwloc_bitmap_t set = hwloc_bitmap_alloc();
    // if ( (ret=hwloc_get_cpubind(hw_topo, set, HWLOC_CPUBIND_THREAD)) != 0) {
    //   fprintf(stderr,"ERROR: hwloc_get_cpubind(...) returned %d\n", ret);
    // }


    hwloc_cpuset_t last = hwloc_bitmap_alloc();
    if ( (ret=hwloc_get_last_cpu_location(hw_topo, last, HWLOC_CPUBIND_THREAD)) != 0) {
      fprintf(stderr,"ERROR: hwloc_get_last_cpu_location(...) returned %d\n", ret);
    }
    if (hwloc_bitmap_iszero(last)) {
      fprintf(stderr,"ERROR: hwloc_get_last_cpu_location returned a zeroed bit vector.\n");
    }

    
    //hwloc_bitmap_singlify(set);

    long page_start = (long) start;
    long last_page = (long) start + (long) sz_in_bytes;
    page_start = page_start & ~PAGE_MSK;
    last_page = last_page & ~PAGE_MSK;
    long nbytes = last_page - page_start + PAGE_SZ;
    
    unsigned cpu, node;
    if ( (ret=sched_getcpu()) < 0) {
      printf("getcpu(...) returned %d\n",ret);
    }
    else {
      cpu = ret;
      node = -1;
    }

    char buffer[256];
    hwloc_bitmap_snprintf(buffer,256,last);
    printf("migrate: start = 0x%X, page_start = 0x%X, #bytes = %ld (nbytes = %ld), last = %s, getcpu(cpu=%d,node=%d)\n", 
	   (long) start, (long) page_start, (long) sz_in_bytes, nbytes, buffer,cpu,node);

    if ( (ret = hwloc_set_area_membind(hw_topo, (const void*) page_start, nbytes,
				       (hwloc_const_cpuset_t)last, HWLOC_MEMBIND_BIND,
				       HWLOC_MEMBIND_MIGRATE )) != 0 ) {
      int err = errno;
      fprintf(stderr, "hwloc_set_area_membind(hw_topo, page_start=0x%X, sz=%lld, set, HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_MIGRAGE) = %d\n",
	      page_start, nbytes, ret);
      
      fprintf(stderr, "errno (%d) = %s\n", err, strerror(err));

    }
    //hwloc_bitmap_free(set);
    hwloc_bitmap_free(last);
  }  

  unsigned long hw_cache_sz;

private:

  hwloc_topology_t hw_topo;

};

TPZ_HWLOC_Manager HWLOC_Manager;



#endif // USING_HWLOC


void migrate_to_local(char* start, unsigned long long sz_in_bytes)
{
#ifdef USING_HWLOC
  HWLOC_Manager.migrate_to_local(start, sz_in_bytes);
#endif
}

