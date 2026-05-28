include(FetchContent)

function(find_or_fetch DEPENDENCY)

    if(DEFINED ENV{CODE} AND EXISTS $ENV{CODE}/${DEPENDENCY})
        set(SOURCE_DIR $ENV{CODE}/${DEPENDENCY})
        message(STATUS "Using ${DEPENDENCY} in $ENV{CODE}/${DEPENDENCY}")
    else()
        fetch(${DEPENDENCY} SOURCE_DIR)
    endif()

    add_subdirectory(${SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${DEPENDENCY}
        EXCLUDE_FROM_ALL
    )
endfunction()

function(fetch DEPENDENCY SOURCE_DIR)
    set(REPO_URL https://github.com/itpplasma/${DEPENDENCY}.git)
    # <DEP>_BRANCH (cache or env) overrides the ref, so an upstream release can
    # build this code against a candidate. Unset keeps the matching-branch default.
    string(TOUPPER ${DEPENDENCY} _DEP_UPPER)
    if(DEFINED ${_DEP_UPPER}_BRANCH AND NOT "${${_DEP_UPPER}_BRANCH}" STREQUAL "")
        set(REMOTE_BRANCH "${${_DEP_UPPER}_BRANCH}")
    elseif(DEFINED ENV{${_DEP_UPPER}_BRANCH} AND NOT "$ENV{${_DEP_UPPER}_BRANCH}" STREQUAL "")
        set(REMOTE_BRANCH "$ENV{${_DEP_UPPER}_BRANCH}")
    else()
        get_branch_or_main(${REPO_URL} REMOTE_BRANCH)
    endif()
    message(STATUS "Fetch ${DEPENDENCY} branch ${REMOTE_BRANCH} from ${REPO_URL}")


    FetchContent_Declare(
        ${DEPENDENCY}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        GIT_REPOSITORY ${REPO_URL}
        GIT_TAG ${REMOTE_BRANCH}
    )
    FetchContent_Populate(${DEPENDENCY})

    string(TOLOWER ${DEPENDENCY} DEPENDENCY_LOWER)
    message(STATUS "Fetched ${DEPENDENCY} to ${${DEPENDENCY_LOWER}_SOURCE_DIR}")

    set(SOURCE_DIR ${${DEPENDENCY_LOWER}_SOURCE_DIR} PARENT_SCOPE)
endfunction()

function(get_branch_or_main REPO_URL REMOTE_BRANCH)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(
        COMMAND git ls-remote --heads ${REPO_URL} ${BRANCH}
        OUTPUT_VARIABLE BRANCH_EXISTS
        ERROR_VARIABLE GIT_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    message(STATUS "Branch ${BRANCH} exists: ${BRANCH_EXISTS}")

    if(BRANCH_EXISTS)
        set(${REMOTE_BRANCH} ${BRANCH} PARENT_SCOPE)
    else()
        set(${REMOTE_BRANCH} "main" PARENT_SCOPE)
    endif()
endfunction()
