add_executable(equilibrium equilibrium.cpp)
target_link_libraries(equilibrium ecell4-ode)

add_executable(equilibrium2 equilibrium2.cpp)
target_link_libraries(equilibrium2 ecell4-ode)

add_executable(odesimulator2 odesimulator2.cpp)
target_link_libraries(odesimulator2 ecell4-ode)

add_executable(dissociation dissociation.cpp)
target_link_libraries(dissociation ecell4-ode)


