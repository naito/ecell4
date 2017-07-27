file(GLOB SOURCES *.cpp)

foreach(SOURCE ${SOURCES})
    get_filename_component(PROGNAME ${SOURCE} NAME_WE)
    add_executable(${PROGNAME} ${SOURCE})
    target_link_libraries(${PROGNAME} ecell4-spatiocyte)
endforeach(SOURCE)
