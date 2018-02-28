/*******************************************************************************
 *   Copyright (C) 2014 by Edson Borin                                         *
 *   edson@ic.unicamp.br                                                       *
 *                                                                             *
 *   This program is free software; you can redistribute it and/or modify      *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation; either version 2 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program; if not, write to the                             *
 *   Free Software Foundation, Inc.,                                           *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 *
 ******************************************************************************/

#include "tpzpagemigrationmanager.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "arglib.h"

#ifdef USING_LIBNUMA
/**
 * @brief Command line arguments habilited by libnuma.
 */
clarg::argBool mig_mp("-mig_mp", "Use move_pages when migrating pages.");
clarg::argBool mig_mbind("-mig_mbind", "Use mbind when migrating pages.");
#endif

#ifdef USING_HWLOC
/**
 * @brief Command line arguments habilited by hwloc.
 */
clarg::argBool mig_hwloc("-mig_hwloc", "Use hwloc when migrating pages.");
#endif

/**
 * @brief Initialize hwloc tolopology when the lib is enabled.
 */
TPZPageMigrationManager::TPZPageMigrationManager() {
#ifdef USING_HWLOC
    hwloc_topology_init(&hw_topo);
    hwloc_topology_load(hw_topo);

    hwloc_obj_t obj;
    HwCacheSize = 0;
    for (obj = hwloc_get_obj_by_type(hw_topo, HWLOC_OBJ_PU, 0); obj; obj = obj->parent)
        if (obj->type == HWLOC_OBJ_CACHE)
            HwCacheSize += obj->attr->cache.size;
#endif
}

/**
 * @brief Destroy hwloc topology when the lib is enabled.
 */
TPZPageMigrationManager::~TPZPageMigrationManager() {
#ifdef USING_HWLOC
    hwloc_topology_destroy(hw_topo);
#endif
}

#if defined (USING_HWLOC) || defined (USING_LIBNUMA)
/**
 * @brief Interface to call a specific technique to migrate data to the
 * local node
 * @param Pointer to adress that will be migrate.
 * @param Size in bytes of the data to be migrate.
 */
void TPZPageMigrationManager::MigrateToLocal(char* start, uint64_t size_in_bytes) {
#ifdef USING_LIBNUMA
    if (mig_mp.was_set()) {
        MigrateToLocalMovePages(start, size_in_bytes);
		return;
    }
    else if (mig_mbind.was_set()) {
        MigrateToLocalMbind(start, size_in_bytes);
		return;
    }
#endif

#ifdef USING_HWLOC
    if (mig_hwloc.was_set()) {
        MigrateToLocalHwloc(start, size_in_bytes);
    }
#endif
}
#endif

#ifdef USING_LIBNUMA

const int64_t PAGE_EXP  = 12;
const int64_t PAGE_SZ   = (1<<12);
const int64_t PAGE_MSK  = ((1<<PAGE_EXP)-1);

/**
 * @brief Migrate data using the system call move_pages.
 * @param Pointer to adress that will be migrate.
 * @param Size in bytes of the data to be migrate.
 */
void TPZPageMigrationManager::MigrateToLocalMovePages(char* start, uint64_t size_in_bytes) {
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
    uint64_t nodemask = 1 << node;
    uint64_t maxnode = 9;
    int64_t page_start = (int64_t) start;
    page_start = page_start & ~PAGE_MSK;
    int64_t last_page = (int64_t) start + (int64_t) size_in_bytes;
    last_page = last_page & ~PAGE_MSK;
    int64_t count = (last_page - page_start + PAGE_SZ) >> PAGE_EXP;
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
        fprintf(stderr, "move_pages(0, npages(%ld), pages, nodes (%u), status, MPOL_MF_MOVE) = %d (cpu=%d, node=%d)\n",
                count, node, ret, cpu, node);
        fprintf(stderr, "errno (%d) = %s\n", err, strerror(err));
    }

}

/**
 * @brief Migrate data using the system call mbind.
 * @param Pointer to adress that will be migrate.
 * @param Size in bytes of the data to be migrate.
 */
void TPZPageMigrationManager::MigrateToLocalMbind(char* start, uint64_t size_in_bytes) {
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
    uint64_t nodemask = 1 << node;
    //uint64_t maxnode = 7;
    //uint64_t nodemask = 1;
    uint64_t maxnode = 9;
    int64_t page_start = (int64_t) start;
    page_start = page_start & ~PAGE_MSK;
    int64_t last_page = (int64_t) start + (int64_t) size_in_bytes;
    last_page = last_page & ~PAGE_MSK;
    int64_t nbytes = last_page - page_start + PAGE_SZ;
#ifndef MPOL_F_STATIC_NODES
#define MPOL_F_STATIC_NODES  (1 << 15)
#endif
    if ( (ret = mbind((void*) page_start, nbytes, MPOL_BIND  , &nodemask, maxnode, /* flags */ MPOL_MF_MOVE) ) != 0 ) {
        // Error.
        int err = errno;
        fprintf(stderr,"ERROR: mbind(start=%ld, sz=%ld, MPOL_BIND, nodemask = %08ld (node = %u), maxnode = %lu, MPOL_MF_MOVE) = %d\n",
                page_start, nbytes, nodemask, node, maxnode, ret);
        fprintf(stderr, "errno (%d) = %s\n", err, strerror(err));
        return;
    }
    else
	{
#ifdef PZDEBUG
        fprintf(stderr,"mbind(start=%ld, sz=%ld, MPOL_BIND, nodemask = %08ld (node = %d), maxnode = %lu, MPOL_MF_MOVE) = %d\n",
                page_start, nbytes, nodemask, node, maxnode, ret);
#endif
    }
}
#endif

#ifdef USING_HWLOC

/**
 * @brief Migrate data using the hwloc library.
 * @param Pointer to adress that will be migrate.
 * @param Size in bytes of the data to be migrate.
 */
void TPZPageMigrationManager::MigrateToLocalHwloc(char* start, uint64_t size_in_bytes) {
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
    int64_t page_start = (int64_t) start;
    int64_t last_page = (int64_t) start + (int64_t) size_in_bytes;
    page_start = page_start & ~PAGE_MSK;
    last_page = last_page & ~PAGE_MSK;
    int64_t nbytes = last_page - page_start + PAGE_SZ;
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
#ifdef PZDEBUG
    printf("migrate: start = 0x%ld, page_start = 0x%ld, #bytes = %ld (nbytes = %ld), last = %s, getcpu(cpu=%d,node=%d)\n", 
           (int64_t) start, (int64_t) page_start, (int64_t) size_in_bytes, nbytes, buffer,cpu,node);
#endif
    if ( (ret = hwloc_set_area_membind(hw_topo, (const void*) page_start, nbytes,
                                       (hwloc_const_cpuset_t)last, HWLOC_MEMBIND_BIND,
                                       HWLOC_MEMBIND_MIGRATE )) != 0 ) {
        int err = errno;
        fprintf(stderr, "hwloc_set_area_membind(hw_topo, page_start=0x%ld, sz=%ld, set, HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_MIGRAGE) = %d\n",
                page_start, nbytes, ret);
        fprintf(stderr, "errno (%d) = %s\n", err, strerror(err));
    }
    hwloc_bitmap_free(last);
}
#endif

TPZPageMigrationManager MigrationManager;

void migrate_to_local(char* start,uint64_t sz_in_bytes) {	
#if defined (USING_HWLOC) || defined (USING_LIBNUMA)
	MigrationManager.MigrateToLocal(start, sz_in_bytes);
#else
	fprintf(stderr, "PZ compiled without page migration support. Compile with USING_HWLOC or USING_LIBNUMA.");
#endif
}
