# These programs use internal metis data structures.
include_directories(../libmetis)
link_directories(/home/karypis/local/lib)
# Build program.
add_executable(gpmetis_ve gpmetis_ve.cpp cmdline_gpmetis.c io.c stat.c)
add_executable(gpmetis_vlmv gpmetis_vlmv.cpp cmdline_gpmetis.c io.c stat.c)
add_executable(gpmetis gpmetis.c cmdline_gpmetis.c io.c stat.c)
foreach(prog gpmetis_ve gpmetis_vlmv gpmetis)
  target_link_libraries(${prog} metis)
#  target_link_libraries(${prog} metis profiler)
endforeach(prog)

if(METIS_INSTALL)
  install(TARGETS gpmetis_vlmv gpmetis_ve gpmetis
    RUNTIME DESTINATION bin)
endif()

# Try to find subversion revision.
set(SVNREV "")
file(TO_NATIVE_PATH ${PROJECT_SOURCE_DIR}/.svn svn_dir)
if(IS_DIRECTORY ${svn_dir})
  include(FindSubversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} metis)
    set(SVNREV ${metis_WC_REVISION})
  endif(Subversion_FOUND)
endif()
add_definitions(-DSVNINFO="${SVNREV}")
