# fortplot revision shared by NEO-RT and POTATO. Include before libneo, so this
# pin is the one built and differing requests are reported.
itp_add_dependency(fortplot REF 35a9cd26ab20cf50ba6f7347d7daf13566e0dc91
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    TARGET fortplot)
