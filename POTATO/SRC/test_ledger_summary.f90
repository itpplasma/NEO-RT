program test_ledger_summary
  use potato_ledger_summary_mod, only : summarize_processed_ledger
  implicit none

  double precision :: bounds(3),maximum_bound
  logical :: certified(3),all_certified

  bounds=(/huge(1.d0),2.d-3,4.d-3/)
  certified=(/.false.,.true.,.true./)
  call summarize_processed_ledger(2,3,bounds,certified,maximum_bound,all_certified)
  if(abs(maximum_bound-4.d-3).gt.1.d-15) error stop 'processed ledger maximum is wrong'
  if(.not.all_certified) error stop 'processed ledger certification is wrong'

  call summarize_processed_ledger(1,3,bounds,certified,maximum_bound,all_certified)
  if(maximum_bound.ne.huge(1.d0)) error stop 'unprocessed ledger slot was ignored incorrectly'
  if(all_certified) error stop 'unprocessed ledger certification was ignored incorrectly'
end program test_ledger_summary
