add_definitions(-Dusing_pthread)
#add_definitions(-D__GNU__)

include_directories(
	"${PROJECT_SOURCE_DIR}/Projects/2017Arequipa/SimworxLib"
	"${PROJECT_SOURCE_DIR}/Projects/2017Arequipa/SimworxLib/AnothersToNeoPZ")


add_subdirectory(SimworxLib)

#add_subdirectory(ConsLawWithWavelets)
add_subdirectory(EulerConsLaw)
add_subdirectory(Torre)
add_subdirectory(Viga)


set(aqp2017 ${SWXTutorialTorre} ${SWXTutorialViga} ${EulerWithWavelets} PARENT_SCOPE )
