# fortplot revision shared by NEO-RT and POTATO. Include before libneo, so this
# pin is the one built and differing requests are reported.
itp_add_dependency(fortplot REF 1645fa41c11b5786c71bc27e5b41dae7c7bbba23
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    TARGET fortplot)
