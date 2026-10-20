execute_process(COMMAND "${PROGRAM}" nonlinear
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0")
    message(FATAL_ERROR "Nonlinear electron initialization was accepted")
endif()
if(NOT "${output}${error}" MATCHES "nonlinear electron collisions are unsupported")
    message(FATAL_ERROR "Expected explicit nonlinear electron rejection: ${output}${error}")
endif()
