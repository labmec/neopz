#ifndef _PZ_SHARED_EXPORTS_H__
#define _PZ_SHARED_EXPORTS_H__

/*MSVC makes no symbols visible by default.
The CMake option

CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS

can help us with that. However, we still need to take care of
static and global variables.

The macro PZ_SHARED_EXPORT will expand to nothing
on UNIX systems, but will do the job on Windows.

CMake has a functionality to generate this file for us, but I think
that letting the macros easy to be found will be less confusing.

Usage for static member variables is as follows:
class TPZClass{
 . . . 
  static type PZ_EXPORT var;
}
*/


#ifdef PZ_STATIC_DEFINE
#  define PZ_EXPORT
#  define PZ_NO_EXPORT
#else
#  ifndef PZ_EXPORT
#    ifdef pz_EXPORTS
        /* We are building this library */
#      ifdef WIN32
#        define PZ_EXPORT __declspec(dllexport)
#      else
#        define PZ_EXPORT
#      endif
#    else
        /* We are using this library */
#      ifdef WIN32
#        define PZ_EXPORT __declspec(dllimport)
#      else
#        define PZ_EXPORT
#      endif
#    endif
#  endif

#  ifndef PZ_NO_EXPORT
#    define PZ_NO_EXPORT 
#  endif
#endif

#ifndef PZ_DEPRECATED
#  define PZ_DEPRECATED __declspec(deprecated)
#endif

#ifndef PZ_DEPRECATED_EXPORT
#  define PZ_DEPRECATED_EXPORT PZ_EXPORT PZ_DEPRECATED
#endif

#ifndef PZ_DEPRECATED_NO_EXPORT
#  define PZ_DEPRECATED_NO_EXPORT PZ_NO_EXPORT PZ_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef PZ_NO_DEPRECATED
#    define PZ_NO_DEPRECATED
#  endif
#endif


#endif /* _PZ_SHARED_EXPORTS_H__ */