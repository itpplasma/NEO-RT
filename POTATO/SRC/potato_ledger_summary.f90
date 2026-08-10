module potato_ledger_summary_mod
  implicit none

contains

  subroutine summarize_processed_ledger(first_energy,last_energy,bounds,certified, &
      maximum_bound,all_certified)
    integer, intent(in) :: first_energy,last_energy
    double precision, intent(in) :: bounds(:)
    logical, intent(in) :: certified(:)
    double precision, intent(out) :: maximum_bound
    logical, intent(out) :: all_certified

    maximum_bound=huge(1.d0)
    all_certified=.false.
    if(first_energy.lt.1 .or. last_energy.gt.size(bounds) .or. &
        last_energy.gt.size(certified) .or. first_energy.gt.last_energy) return
    maximum_bound=maxval(bounds(first_energy:last_energy))
    all_certified=all(certified(first_energy:last_energy))
  end subroutine summarize_processed_ledger

end module potato_ledger_summary_mod
