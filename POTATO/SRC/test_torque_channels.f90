program test_torque_channels
  use resint_mod, only : NCHAN, ICH_TRAPPED, ICH_PASS_VPAR_POS, &
                         ICH_PASS_VPAR_NEG, channel_of
  implicit none

  real(8), parameter :: co_going = 1d0, counter_going = -1d0

  if (NCHAN /= 3) error stop 'expected trapped and two passing channels'
  if (channel_of(.true., co_going) /= ICH_TRAPPED) &
    error stop 'trapped co-going half must go to the trapped channel'
  if (channel_of(.true., counter_going) /= ICH_TRAPPED) &
    error stop 'trapped counter-going half must go to the trapped channel'
  if (channel_of(.false., co_going) /= ICH_PASS_VPAR_POS) &
    error stop 'passing orbit with v_par > 0 has the wrong channel'
  if (channel_of(.false., counter_going) /= ICH_PASS_VPAR_NEG) &
    error stop 'passing orbit with v_par < 0 has the wrong channel'
  print *, 'test_torque_channels: OK'
end program test_torque_channels
