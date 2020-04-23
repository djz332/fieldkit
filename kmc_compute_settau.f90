  module kmc_compute_rates

  use omp_lib
  
  implicit none

  integer, parameter :: DP_INT = SELECTED_INT_KIND(16)
  integer, parameter :: DP_EIGHT = SELECTED_INT_KIND(8)
  integer, parameter :: DP = SELECTED_REAL_KIND(16)

  integer :: i, j, k, x
  integer(kind = 2) :: proxy
  integer :: status_neighbor(6)
  real :: dsite, rhosite, dd_dx, dd_dy, dd_dz, dphi_dx, dphi_dy, dphi_dz
  real :: dsite_neighbor(6), rhosite_neighbor(6)
  real(kind = 8) :: dx, dy, dz
  real(kind = 8) :: delta(3) !input from python

  logical :: settauflag

  private :: i, j, k, DP_INT, DP, DP_EIGHT 


  contains


  subroutine settau_calc(domain, dsite_input, density_input, Lx, Ly, Lz, settau_rates)

  integer, intent(in) :: Lx, Ly, Lz
  integer, dimension(Lx,Ly,Lz), intent(in) :: domain
  real(kind = 8), dimension(Lx,Ly,Lz), intent(in) :: density_input, dsite_input
  real, dimension(6,Lx,Ly,Lz), intent(out) :: settau_rates

  dx = delta(1)
  dy = delta(2)
  dz = delta(3)

  do k = 1,Lz
   
   do j = 1,Ly
   
    do i = 1,Lx 

       rhosite = density_input(i,j,k)
       dsite = dsite_input(i,j,k)

       settauflag = .True.

       if (Domain(i,j,k) .eq. 0) then

          do x = 1,6
  
             settau_rates(x,i,j,k) = 0

          end do
  
          settauflag = .False.

       end if

       if (settauflag) then
                          
          proxy = proxy_func(i, Lx, 1) 
          dsite_neighbor(1) = dsite_input(proxy, j, k) 
          rhosite_neighbor(1) = density_input(proxy, j, k)
          status_neighbor(1) = Domain(proxy, j, k)

          proxy = proxy_func(i, Lx, 0)
          dsite_neighbor(2) = dsite_input(proxy, j, k)
          rhosite_neighbor(2) = density_input(proxy, j, k)
          status_neighbor(2) = Domain(proxy, j, k)

          proxy = proxy_func(j, Ly,1)
          dsite_neighbor(3) = dsite_input(i, proxy, k)
          rhosite_neighbor(3) = density_input(i ,proxy, k)
          status_neighbor(3) = Domain(i,  proxy, k)

          proxy = proxy_func(j, Ly, 0)
          dsite_neighbor(4) = dsite_input(i, proxy, k)
          rhosite_neighbor(4) = density_input(i, proxy, k)
          status_neighbor(4) = Domain(i, proxy, k)

          proxy = proxy_func(k, Lz, 1)
          dsite_neighbor(5) = dsite_input(i, j, proxy)
          rhosite_neighbor(5) = density_input(i, j, proxy)
          status_neighbor(5) = Domain(i, j, proxy)

          proxy = proxy_func(k, Lz, 0)
          dsite_neighbor(6) = dsite_input(i, j, proxy)
          rhosite_neighbor(6) = density_input(i, j, proxy)
          status_neighbor(6) = Domain(i, j, proxy)

          dd_dx = diff_calc(dsite, dsite_neighbor(1), dsite_neighbor(2), status_neighbor(1), status_neighbor(2), dx, 1)
          dd_dy = diff_calc(dsite, dsite_neighbor(3), dsite_neighbor(4), status_neighbor(3), status_neighbor(4), dy, 1)
          dd_dz = diff_calc(dsite, dsite_neighbor(5), dsite_neighbor(6), status_neighbor(5), status_neighbor(6), dz, 1)
          dphi_dx = diff_calc(rhosite, rhosite_neighbor(1), rhosite_neighbor(2), status_neighbor(1), status_neighbor(2), dx, 2)
          dphi_dy = diff_calc(rhosite, rhosite_neighbor(3), rhosite_neighbor(4), status_neighbor(3), status_neighbor(4), dy, 2)
          dphi_dz = diff_calc(rhosite, rhosite_neighbor(5), rhosite_neighbor(6), status_neighbor(5), status_neighbor(6), dz, 2)


          settau_rates(1,i,j,k) = 0.5*(dsite + 0.5*(dd_dx - (dsite)*(dphi_dx)))

          settau_rates(2,i,j,k) = 0.5*(dsite - 0.5*(dd_dx - (dsite)*(dphi_dx)))

          settau_rates(3,i,j,k) = 0.5*(dsite + 0.5*(dd_dy - (dsite)*(dphi_dy)))

          settau_rates(4,i,j,k) = 0.5*(dsite - 0.5*(dd_dy - (dsite)*(dphi_dy)))

          settau_rates(5,i,j,k) = 0.5*(dsite + 0.5*(dd_dz - (dsite)*(dphi_dz)))

          settau_rates(6,i,j,k) = 0.5*(dsite - 0.5*(dd_dz - (dsite)*(dphi_dz)))
  

          do x = 1,6

             if ((settau_rates(x,i,j,k)) .lt. 0) then

                PRINT*, 'SETTAU AT SITE COORDINATES', i, j, k, ' IN DIRECTION ', x, ' IS NEGATIVE' 
                PRINT*, STATUS_NEIGHBOR(x)
                PRINT*, 'DSITE: ', DSITE, 'DD_DX ', DD_DX, 'DD_DY ', DD_DY, 'rhosite', rhosite, 'DPHI_DX ', DPHI_DX, 'DPHI_DY ', DPHI_DY 

                STOP 'ERROR'

             end if

          end do  

       end if
 
    end do
 
   end do

  end do  
  
  end subroutine settau_calc
     

  integer function proxy_func(a, L, b)

     implicit none

     integer, intent(in) :: a, L, b

     if (b .eq. 1) then

        if ((a + 1) .le. L) then

            proxy_func = a + 1

        else

            proxy_func = a + 1 - L

        end if

      end if

      if (b .eq. 0) then

          if ((a - 1) .gt. 0) then

              proxy_func = a - 1

          else

              proxy_func = a - 1 + L

          end if

       end if  

       return

   end function proxy_func    


   real function diff_calc(site_a,site_b,site_c, status_b, status_c, del, d)

     implicit none

     real, intent(in) :: site_a, site_b, site_c
     integer, intent(in) :: status_b, status_c, d
     real(kind = 8), intent(in) :: del

     real :: dsiteb, dsitec, rhositeb, rhositec
     
     if (d .eq. 1) then

        dsite = site_a
        dsiteb = site_b
        dsitec = site_c

        if((status_b .eq. 0) .or. (status_c .eq. 0)) then
        
           if (status_b .ne. 0) then
              
              diff_calc = dble(dble(dsiteb) - dble(dsite))/dble(del)

           else if (status_c .ne. 0) then

              diff_calc = dble(dble(dsite) - dble(dsitec))/dble(del)

           else if((status_b .eq. 0) .and. (status_c .eq. 0)) then 
             
              diff_calc = 0
 
           end if
        
        else

           diff_calc = dble(dble(dsiteb) - dble(dsitec))/dble(2)/dble(del)

        end if


     else if (d .eq. 2) then

        rhosite = site_a
        rhositeb = site_b
        rhositec = site_c

        if((status_b .eq. 0) .or. (status_c.eq. 0)) then

           if (status_b .ne. 0) then

              diff_calc = dble(-1)*dble(dble(rhositeb) - dble(rhosite))/dble(del)/dble(rhosite)

           else if (status_c .ne. 0) then

              diff_calc = dble(-1)*dble(dble(rhosite) - dble(rhositec))/dble(del)/dble(rhosite)

           else if((status_b.eq. 0) .and. (status_c .eq. 0)) then

              diff_calc = 0

           end if

        else

           diff_calc = dble(-1)*dble(dble(rhositeb) - dble(rhositec))/dble(2)/dble(del)/dble(rhosite)

        end if

     else

        print*, 'WRONG INPUT PARAMETERS FOR DIFF_CALC FUNCTION'

        stop 'ERROR'

     end if

   end function diff_calc

  
  end module kmc_compute_rates
