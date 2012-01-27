########### next target ###############

file(GLOB locheaders *.h)

file(GLOB locsources *.cpp)

configure_file(log4cxx.cfg log4cxx.cfg COPYONLY)



set(MarxLang_SRCS ConvertionMethods.cpp PropertiesTable.cpp
    ThermalMethodsTables.cpp TPBConservacao1DMarx.cpp 
    TSwxSteamInjectionData.cpp tpbrsolutionlist.cpp tpbrthermaldisc.cpp)

set(main_SRCS ${locheaders} ${MarxLang_SRCS})

add_executable(Steam_MarxLang ${main_SRCS} ${locheader})

target_link_libraries(Steam_MarxLang pz pthread)

set(steam_SRCS tpbrsolutionlist.cpp tpbrsteamflux.cpp tpbrsteammesh.cpp
    TPBrSteamSimulation.cpp tpbrthermaldisc.cpp tpbrthermalsolution.cpp ThermalMethodsTables.cpp ConvertionMethods.cpp tpbrcellconservation.cpp PropertiesTable.cpp
    ${locheaders})

add_executable(Steam ${steam_SRCS} ${locheader})

target_link_libraries(Steam pz pthread)

