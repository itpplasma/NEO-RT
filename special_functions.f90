subroutine comelp ( hk, ck, ce )

!*****************************************************************************80
!
!! COMELP computes complete elliptic integrals K(k) and E(k).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) HK, the modulus.  0 <= HK <= 1.
!
!    Output, real ( kind = 8 ) CK, CE, the values of K(HK) and E(HK).
!
  implicit none

  real ( kind = 8 ) ae
  real ( kind = 8 ) ak
  real ( kind = 8 ) be
  real ( kind = 8 ) bk
  real ( kind = 8 ) ce
  real ( kind = 8 ) ck
  real ( kind = 8 ) hk
  real ( kind = 8 ) pk

  pk = 1.0D+00 - hk * hk

  if ( hk == 1.0D+00 ) then

    ck = 1.0D+300
    ce = 1.0D+00

  else

    ak = ((( &
        0.01451196212D+00   * pk &
      + 0.03742563713D+00 ) * pk &
      + 0.03590092383D+00 ) * pk &
      + 0.09666344259D+00 ) * pk &
      + 1.38629436112D+00

    bk = ((( &
        0.00441787012D+00   * pk &
      + 0.03328355346D+00 ) * pk &
      + 0.06880248576D+00 ) * pk &
      + 0.12498593597D+00 ) * pk &
      + 0.5D+00

    ck = ak - bk * log ( pk )

    ae = ((( &
        0.01736506451D+00   * pk &
      + 0.04757383546D+00 ) * pk &
      + 0.0626060122D+00  ) * pk &
      + 0.44325141463D+00 ) * pk &
      + 1.0D+00

    be = ((( &
        0.00526449639D+00   * pk &
      + 0.04069697526D+00 ) * pk &
      + 0.09200180037D+00 ) * pk &
      + 0.2499836831D+00  ) * pk

    ce = ae - be * log ( pk )

  end if

  return
end subroutine comelp
