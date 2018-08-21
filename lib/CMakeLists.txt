# @file neopz/lib/CMakeList.txt  -- CMake file grouping pz modules and creating pz library. 

# creating pz lib and groups of the sources files
set(pzsources 
	${pzanalysis_STAT_SRCS} 
	${pzcommon_STAT_SRCS} 
	${pzexception_STAT_SRCS} 
	${pzexternal_STAT_SRCS} ${sloan_SRCS} ${FAD_SRCS} ${FAD_fad_SRCS} ${FAD_Tiny_SRCS} ${FAD_TinyET_SRCS}
	${pzfrontal_STAT_SRCS} 
	${pzgeom_STAT_SRCS} 
	${pzintegral_STAT_SRCS} 
    ${pzlinearsolvers_STAT_SRCS}
	${pzmaterial_STAT_SRCS} 
	${pzmatrix_STAT_SRCS} 
	${pzmesh_STAT_SRCS} 
	${pzmultigrid_STAT_SRCS} 
	${pzperfutil_STAT_SRCS} 
	${pzpre_STAT_SRCS} 
	${pzpost_STAT_SRCS}
	${pzrefine_STAT_SRCS} 
	${pzsave_STAT_SRCS} 
	${pzshape_STAT_SRCS} 
	${pzspecialmaps_STAT_SRCS} 
	${pzstrmatrix_STAT_SRCS} 
	${pzpublications_STAT_SRCS}
	${pzsubstruct_STAT_SRCS}
    ${pztopology_STAT_SRCS} 
	${pzutil_STAT_SRCS} 
	${pzrandom_STAT_SRCS} 
	${pzoptimization_STAT_SRCS}  
	)

target_sources( pz PRIVATE ${pzsources})
