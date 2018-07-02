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

/**
 * @file
 * @brief Contains declaration of the TPZPageMigrationManager class which 
 * implements methods to migrate data in the NUMA architecture.
 */

#ifndef TPZ_PAGEMIGRATIONMANAGER_H
#define TPZ_PAGEMIGRATIONMANAGER_H

#include <stdint.h>

#ifdef USING_HWLOC
#include <hwloc.h>
#endif

#ifdef USING_LIBNUMA
#include <numaif.h>
#endif

#include <stdio.h>

class TPZPageMigrationManager {
	public:
		/**
		 * @brief Initialize hwloc tolopology when the lib is enabled.
		 */
		TPZPageMigrationManager();
		/**
		 * @brief Destroy hwloc topology when the lib is enabled.
		 */
		~TPZPageMigrationManager();
#if defined (USING_HWLOC) || defined (USING_LIBNUMA)
		/**
		 * @brief Interface to call a specific technique to migrate data to the
		 * local node
		 * @param start : Pointer to adress that will be migrate.
		 * @param size_in_bytes : Size in bytes of the data to be migrate.
		 */
		void MigrateToLocal(char* start, uint64_t size_in_bytes);
#endif

	private:
#ifdef USING_LIBNUMA
		/**
		 * @brief Migrate data using the system call move_pages.
		 * @param start : Pointer to adress that will be migrate.
		 * @param size_in_bytes : Size in bytes of the data to be migrate.
		 */
		void MigrateToLocalMovePages(char* start, uint64_t size_in_bytes);
		/**
		 * @brief Migrate data using the system call mbind.
		 * @param start : Pointer to adress that will be migrate.
		 * @param size_in_bytes : Size in bytes of the data to be migrate.
		 */
		void MigrateToLocalMbind(char* start, uint64_t size_in_bytes);
#endif
#ifdef USING_HWLOC
		/**
		 * @brief Migrate data using the hwloc library.
		 * @param start : Pointer to adress that will be migrate.
		 * @param size_in_bytes : Size in bytes of the data to be migrate.
		 */
		void MigrateToLocalHwloc(char* start, uint64_t size_in_bytes);
		/**
		 * @brief Hwloc Topology store data about the hardware.
		 */
		hwloc_topology_t hw_topo;
		/**
		 * @brief Store the sum of all the processing units caches sizes.
		 */
		uint64_t HwCacheSize;
#endif
}; // class TPZPageMigrationManager


void migrate_to_local(char* start, uint64_t sz_in_bytes);

#endif
