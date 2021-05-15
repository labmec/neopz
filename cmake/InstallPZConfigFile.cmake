# We need this work around since CMake variables are not available inside install(SCRIPT)
# Reference: https://stackoverflow.com/questions/20792802/why-variables-are-not-accessed-inside-script-in-cmake
function(create_post_install_var var)
    install(CODE "set(${var} ${${var}})")
endfunction()

create_post_install_var(REAL_TYPE_DEF)
create_post_install_var(STATE_TYPE_DEF)

create_post_install_var(PZ_BRANCH)
create_post_install_var(PZ_REVISION)
create_post_install_var(PZ_REVISION_DATE)

create_post_install_var(VC)
create_post_install_var(PZ_LOG)
create_post_install_var(CMAKE_INSTALL_PREFIX)
create_post_install_var(CMAKE_INSTALL_INCLUDEDIR)
create_post_install_var(PZSOURCEDIR)

install(SCRIPT cmake/PostInstall.cmake)
