! charge_exchange_module.f90
module charge_exchange_module

  private  ! everything is private by default
  public :: load_ct_rates, charge_transfer_recombination, charge_transfer_ionization

  real(KIND=8):: CTRecomb(6,4,31)
  real(KIND=8):: CTIon(7,3,31)

  ! TODO(code) check indices in this routine

CONTAINS

SUBROUTINE load_ct_rates()
  !Load the charge transfer ionization
  !and recombination rates from file
  implicit none

  integer:: i, j, k, unit_num, ios

  write(*,*) 'Initializing charge exchange rates'

  ! Load Ionization
  open(newunit=unit_num, file='./data/charge_transfer/ct_ionization.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open CT Ionization file'
      return
  end if

  !Reading data from the file into the 3D array
  do i = 3, 31
     do j = 1, 3
        do k = 1, 7
           read(unit_num, *, iostat=ios) CTIon(k,j,i)
           if (ios /= 0) exit
        end do
        if (ios /= 0) exit
     end do
     if (ios /= 0) exit
  end do
  close(unit_num)

  ! Zero out helium
  do j = 1, 3
     do k = 1, 7
        CTIon(k,j,2) = 0.d0
     end do
  end do

  ! Load Recombination
  open(newunit=unit_num, file='./data/charge_transfer/ct_recombination.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open CT Recombination file'
      return
  end if

  !Reading data from the file into the 3D array
  do i = 3, 31
     do j = 1, 4
        do k = 1, 6
           read(unit_num, *, iostat=ios) CTRecomb(k,j,i)
           if (ios /= 0) exit
        end do
        if (ios /= 0) exit
     end do
     if (ios /= 0) exit
  end do
  close(unit_num)

END SUBROUTINE load_ct_rates

FUNCTION charge_transfer_recombination(ion, nelem, T) result(rate)
  ! ion is stage of ionization, 2 for the ion going to the atom
  ! nelem is atomic number of element, 2 up to 30
  ! Example:  O+ + H => O + H+ is HCTRecom(2,8,1e4)
  ! Note that temperature is in linear scale
  implicit none

  integer, intent(in):: ion, nelem
  real(KIND=8), intent(in):: T
  real(KIND=8):: rate
  real(KIND=8):: a_op, b_op, c_op, d_op, e_op, f_op
  real(KIND=8):: logT, tused
  integer:: ipIon

  rate = 0.0

  ! No charge transfer above 10^5 K
  if (T .gt. 1.d5) then
     return
  end if

  ! No recombination on the ground state
  if (ion .eq. 1) then
     return
  end if

  ! Deal with helium separately
  ! He+ + H --> He + H+
  ! if (nelem.eq.2 .and. ion.eq.1) then
  !    rate = 1.20d-15 * ((T/300.d0)**0.25d0)
  !    return
  ! end if

  ! deal with oxygen separately 
  if (nelem .eq. 8) then         
     if (ion .eq. 2) then

        if (T .lt. 10.d0) then
           rate = 3.744d-10
           return
        end if

        a_op = 2.3344302d-10
        b_op = 2.3651505d-10
        c_op = -1.3146803d-10
        d_op = 2.9979994d-11
        e_op = -2.8577012d-12
        f_op = 1.1963502d-13
        logT = log(T)
        rate = ((((f_op*logT + e_op)*logT + d_op)*logT + c_op)*logT + b_op)*logT + a_op

        return
     else if (ion .eq. 3) then
        if (T .le. 1500.0) then
           rate = 0.5337d-9 * ((T/100.d0)**(-0.076d0))
        else 
          rate = 0.4344d-9 + (0.6340d-9 * (log10(T/1500.0)**2.06d0));
        end if
        
        return
      end if
  end if

  ! Deal with nitrogen separately
  if (nelem.eq.7 .and. ion.eq.2) then
     ! N+2 + H -> N+ + H+ 
	 if (T .le. 1500d0) then
	    rate = 0.8692d-9 * ((T/1500d0)**0.17d0)
	 else if (T.le.20000.d0) then
	    rate = 0.9703d-9 * ((T/10000.d0)**0.058d0)
	 else
	    rate = 1.0101d-9 + 1.4589d-9 * ( log10(T/20000.d0)**2.06d0)
     endif

     return
  end if

  ipIon = ion - 1
  !use statistical charge transfer for ion > 4
  if (ion .ge. 4) then
     rate = 1.92d-9 * real(ipIon+1, kind=8)
     return
  end if
    
  !Make sure te is between temp. boundaries; set constant outside of range
  tused = 0.d0
  tused = min(max(T,CTRecomb(5,ipIon,nelem)),CTRecomb(6,ipIon,nelem))
  tused = tused * 1d-4

  ! The interpolation equation
  rate = CTRecomb(1,ipIon,nelem) * 1d-9 * (tused**CTRecomb(2,ipIon,nelem)) * (1.d0 + CTRecomb(3,ipIon,nelem) * exp(CTRecomb(4,ipIon,nelem)*tused) )

END FUNCTION charge_transfer_recombination

FUNCTION charge_transfer_ionization(ion, nelem, T) result(rate)
  ! ion is stage of ionization, 1 for atom
  ! nelem is atomic number of element, 2 up to 30
  ! Example:  O + H+ => O+ + H is HCTIon(1,8,1e4)
  ! Note that temperature is in linear scale
  implicit none
  
  integer, intent(in):: ion, nelem
  real(KIND=8), intent(in):: T
  real(KIND=8):: rate
  real(KIND=8):: a, b, c
  real(KIND=8):: a_o, b_o, c_o, d_o, e_o, f_o, g_o
  real(KIND=8):: logT, tused
  integer:: ipIon

  rate = 0.d0

  ! No charge transfer above 10^5 K
  if (T.ge.1d5) then
     return
  end if

  ! Deal with helium separately
  ! He + H+ --> He+ + H
  if (nelem.eq.2 .and. ion.eq.1) then 
     ! This particular rate seems to get very small at low te
     if (T.lt.3E3) then 
        return 
     end if

     if (T.lt.1E4) then 
        rate = 1.26d-9 * (T**(-0.75d0)) * exp(-1.275d5/T)
     else
        rate = 4d-37 * (T**4.74d0)
     end if

     return
  end if

  ! Deal with oxygen separately
  if (nelem.eq.8 .and. ion.eq.0) then
     if (T .le. 10.d0) then
        rate = 4.749d-20
     else if (T.gt.10.d0 .and. T.le.190.d0) then
        a = -21.134531d0
        b = -242.06831d0
        c = 84.761441d0
        rate = exp(a + (b/T) + (c/(T*T)))
     else if (T.gt.190.d0 .and. T.le.200.d0) then
        rate = 2.18733d-12*(T-190.d0) + 1.85823d-10
     else
        a_o = -7.6767404d-14
        b_o = -3.7282001d-13
        c_o = -1.488594d-12
        d_o = -3.6606214d-12 
        e_o = 2.0699463d-12
        f_o = -2.6139493d-13
        g_o = 1.1580844d-14
        logT = log(T)
        rate = (((((g_o*logT + f_o)*logT + e_o)*logT + d_o)*logT + c_o)*logT + b_o)*logT + a_o
     end if

     return
  end if

  ! Deal with iron separately
  ! Fe + H+ -> Fe+ + H
  if (nelem.eq.26 .and. ion.eq.0) then
     rate = 1.d-14
     return
  end if

  ! Deal with sulfur separately
  ! S + H+ -> S+ + H
  if (nelem.eq.16 .and. ion.eq.0) then 
     rate = 5.4d-9
     return
  end if

  ! Deal with magnesium separately
  ! Mg + H+ -> Mg+ + H
  if (nelem.eq.12 .and. ion .eq. 0) then 
     rate = 9.76d-12*((T/1d4)**3.14d0)*(1.d0 + 55.54d0*exp(1.12d0*T/1d4))
     return
  end if

  ! Deal with Silicon separately
  ! Si + H+ -> Si+ + H
  if (nelem.eq.14 .and. ion.eq.0) then
     rate = 0.92d-12*((T/1e4)**1.15d0)*(1.d0 + 0.80d0*exp(0.24d0*T/1d4))
     return 
  end if

  ipIon = ion;
  if (ipIon.gt.1) then 
     rate = 0.d0
     return
  end if

  ! Make sure te is between temp. boundaries; set constant outside of range
  tused = 0.d0
  tused = min(max(T,CTIon(5,ipIon,nelem)),CTIon(6,ipIon,nelem))
  tused = tused * 1d-4
  tused = max(tused,1d-10) ! harley added to prevent zero temperature

  ! the interpolation equation
  rate = CTIon(1,ipIon,nelem) * 1d-9 * (tused**CTIon(2,ipIon,nelem)) * (1.d0 + CTIon(3,ipIon,nelem) * exp(CTIon(4,ipIon,nelem)*tused) ) * exp(-1.d0 * CTIon(7,ipIon,nelem)/tused)

END FUNCTION charge_transfer_ionization

end module charge_exchange_module