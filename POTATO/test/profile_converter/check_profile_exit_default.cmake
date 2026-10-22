execute_process(COMMAND "${BINARY}" off RESULT_VARIABLE result
  OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(result EQUAL 0)
  message(FATAL_ERROR "Default-off profile exit did not stop")
endif()
if(NOT error MATCHES "ERROR STOP POTATO orbit is outside tabulated profile support")
  message(FATAL_ERROR "Unexpected default-off failure: ${output}${error}")
endif()
