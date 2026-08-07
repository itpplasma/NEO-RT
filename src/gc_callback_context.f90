module neort_gc_callback_context
    !! Typed root for callback-owned state passed through the cylindrical
    !! transport adapters.  Concrete backends extend this type, which keeps
    !! the adapter boundary type-safe and avoids compiler-generated unlimited-
    !! polymorphic vtables in every consumer of an adapter-containing type.
    implicit none
    private

    type, public :: gc_callback_context_t
    end type gc_callback_context_t

end module neort_gc_callback_context

