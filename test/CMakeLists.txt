add_executable(test_frequencies.x test_frequencies.f90)
target_link_libraries(test_frequencies.x neo_rt)

add_executable(test_reslines.x test_reslines.f90)
target_link_libraries(test_reslines.x neo_rt)

add_executable(test_bounce.x test_bounce.f90)
target_link_libraries(test_bounce.x neo_rt)

add_executable(test_torque.x test_torque.f90)
target_link_libraries(test_torque.x neo_rt)

add_executable(test_misc.x test_misc.f90)
target_link_libraries(test_misc.x neo_rt)
