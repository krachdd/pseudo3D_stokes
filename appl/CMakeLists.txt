add_executable(pseudo3D_stokes EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(pseudo3D_stokes PUBLIC GRID_DIM=2)

dune_symlink_to_source_files(FILES "runStokesGeneric.py" "run_simulations.sh" "run_simulations_test.sh" "run_simulations_periodic_test.sh" "run_simulations_periodic_test_brinkman.sh" "genericInput.input" "genericInputPeriodic.input" "test_singlePrecipitate_2d_total" "test_periodic_BC" "test_singlePrecipitate_permCheck.py" "references") 

dumux_add_test(NAME test_stokesSinglePrecipitate_no
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_no-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_no-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_no.input")

dumux_add_test(NAME test_stokesSinglePrecipitate_quarter
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_quarter-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_quarter-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_quarter.input")

dumux_add_test(NAME test_stokesSinglePrecipitate_half
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_half-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_half-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_half.input")

dumux_add_test(NAME test_stokesSinglePrecipitate_threeQuarter
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_threeQuarter-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_threeQuarter-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_threeQuarter.input")

dumux_add_test(NAME test_stokesSinglePrecipitate_full
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_full-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_full-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_full.input")

dumux_add_test(NAME test_stokesSinglePrecipitate_fiveQuarter
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_fiveQuarter-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_fiveQuarter-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_fiveQuarter.input")

dumux_add_test(NAME test_stokesSinglePrecipitate_threeHalf
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/singlePrecipitate_threeHalf-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total/hx_singlePrecipitate_threeHalf-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes hx_singlePrecipitate_threeHalf.input")
                       
dumux_add_test(NAME test_stokesSinglePrecipitate_permeabilityTest
               TARGET pseudo3D_stokes
               COMMAND ./test_singlePrecipitate_permCheck.py)

dumux_add_test(NAME test_pipesPeriodicBrinkman
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/pipesPeriodicBrinkman-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_periodic_BC/pipesPeriodic-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes pipesPeriodicBrinkman.input")

dumux_add_test(NAME test_ellipsePeriodicBrinkman
              TARGET pseudo3D_stokes
              COMMAND ${CMAKE_SOURCE_DIR}/../dumux/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_CURRENT_BINARY_DIR}/references/ellipsePeriodicBrinkman-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_periodic_BC/ellipsePeriodic-00001.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/pseudo3D_stokes ellipsePeriodicBrinkman.input")

foreach(test_name_singlePrecipitate IN ITEMS test_stokesSinglePrecipitate_no 
                           test_stokesSinglePrecipitate_quarter 
                           test_stokesSinglePrecipitate_half 
                           test_stokesSinglePrecipitate_threeQuarter
                           test_stokesSinglePrecipitate_full
                           test_stokesSinglePrecipitate_fiveQuarter
                           test_stokesSinglePrecipitate_threeQuarter
                           test_stokesSinglePrecipitate_threeHalf)
set_tests_properties(${test_name_singlePrecipitate} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_singlePrecipitate_2d_total)
endforeach()

foreach(test_name_periodic IN ITEMS test_pipesPeriodicBrinkman
                           test_ellipsePeriodicBrinkman)
set_tests_properties(${test_name_periodic} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_periodic_BC)
endforeach()
