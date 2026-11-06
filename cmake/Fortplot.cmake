# fortplot revision shared by NEO-RT and POTATO. Include before libneo:
# FetchContent keeps the first declaration, so this pin wins over transitive ones.
include(FetchContent)
FetchContent_Declare(
    fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG 35a9cd26ab20cf50ba6f7347d7daf13566e0dc91
)
FetchContent_MakeAvailable(fortplot)
