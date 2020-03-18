  module kmc_RW

  use omp_lib

  implicit none

  integer, parameter :: DP_REAL = SELECTED_REAL_KIND(8)
  integer,parameter :: DP_INT = SELECTED_INT_KIND(16)
  
  !!! inputs from python
  integer :: L, L_unscale, numtracer, num_status, prob_freq, print_freq, msd_freq, dfit_points,&
 msdfile_write_freq, num_runs, origin_freq, normalizeflag, case_parameter
  integer(kind = DP_INT) :: totalmoves
  integer, dimension(3) :: msdflag
  real, allocatable, dimension(:,:) :: settau_input
  !!! inputs from python

  integer(kind = DP_INT) :: i, j, k
  integer(8) :: io, dvalue_number
  integer :: neighbor(6)
  integer(kind = 2) :: temp1, temp2, temp3
  integer(kind = 2) :: proxy
  integer :: delta = 1, totalsites, atom_number, num_status_void, counter, binarylevel
  real :: r3
  real(kind = DP_REAL) :: thousand = 1000
   character(len = 30) :: exec, polymerfile, input_flag, input_file, parameter_value_str, numrun_value_str, configuration_string,&
       xyz_string, csv_string, probability_file_string, msd_file_string, dvalues_string, density_string, trajectory_string, &
configuration_file, probability_file, msd_file, dvalues_file, density_file, trajectory_file(30)
  logical :: data_file, insertflag, settauflag
  real,allocatable :: PI(:,:), QI(:,:)

  integer(kind = 1), allocatable, dimension(:) :: siteid_statusid
  real(kind = 4), allocatable, dimension(:,:) :: siteid_settau !6 columns
  integer(kind = 2), allocatable, dimension(:,:,:) :: tracerid_pos, dummyid_pos !3 columns
  integer(kind = 4), allocatable, dimension(:,:) :: tracerid_posn, dummyid_posn
  integer(kind = 1), allocatable, dimension(:,:) :: tracerid_statusid, dummyid_statusid
  real(kind = DP_REAL), allocatable, dimension(:,:,:) :: tracerid_settau, dummyid_settau !6 columns
  real(kind = DP_REAL), allocatable, dimension(:,:,:) :: tracerid_tau, dummyid_tau !6 columns
  real(kind = DP_REAL), allocatable, dimension(:,:) :: tracerid_tausitesum, dummyid_tausitesum
  real(kind = DP_REAL), allocatable, dimension(:,:) :: tracerid_tautracersum, dummyid_tautracersum

  integer :: tracersite, point, replica, direction, prob_counter(30), time_counter(30), msdaverage_counter, filenumber
  integer(kind = 2) :: swap1, swap2
  integer :: swap3
  real(kind = DP_REAL) :: tautotal(30), add(30), tauadd(30), tausubtract(30), tausearch, tausum(30), r1
  real(kind = DP_REAL) :: sumx, sumy, sumxy, sumx2, counter_dfit, dvalue
  real(kind = DP_REAL), allocatable :: taulevel(:,:,:)
  integer(kind = DP_INT), allocatable :: status_prob(:,:)
  real,allocatable :: average_status_prob(:,:)
  real :: pointvalue, msdsum
  double precision :: start, finish
  logical :: randomflag, file_exists

  integer,allocatable:: x_prob(:,:), y_prob(:,:), z_prob(:,:)
  real,allocatable :: final_prob(:,:,:)


  contains 

 
  subroutine create_files()

    num_status_void = 1

     binarylevel = NINT(log(dble(numtracer))/log(dble(2)))


     write(parameter_value_str, '(I7)') case_parameter

     configuration_string = 'configuration_'

     xyz_string = '.xyz'

     csv_string = '.csv'

     probability_file_string = 'prob_'

     msd_file_string = 'msd_'

     dvalues_string = 'dvalues_'

     density_string = 'density_'

     trajectory_string = 'trajectory_'

     configuration_file = trim(trim(configuration_string) // trim(adjustl(parameter_value_str)) // trim(xyz_string))

     probability_file = trim(trim(probability_file_string) // trim(adjustl(parameter_value_str)) // trim(csv_string))

     msd_file = trim(trim(msd_file_string) // trim(adjustl(parameter_value_str)) // trim(csv_string))

     dvalues_file = trim(trim(dvalues_string) // trim(adjustl(parameter_value_str)) // trim(csv_string))

     density_file = trim(trim(density_string) // trim(adjustl(parameter_value_str)) // trim(csv_string))

   do i = 1,num_runs

     write(numrun_value_str, '(I7)') i

     trajectory_file(i) = trim(trim(trajectory_string) // trim(adjustl(parameter_value_str)) // trim("_") // &
trim(adjustl(numrun_value_str)) //trim(xyz_string))

   end do


   do i = 1,num_runs

     OPEN(UNIT = 2, FILE = trajectory_file(i))

         WRITE(2,*) 'TRAJECTORY FILE FOR TRACERS, REPLICA #:', i

     CLOSE(2)

   end do

  end subroutine create_files


  subroutine grid_generator()

     allocate(siteid_statusid((L**(3)+1)))
     allocate(siteid_settau((L**(3)),6))

     totalsites = L**3

  end subroutine grid_generator


  subroutine status_allocation()

      do i = 1,totalsites+1

          siteid_statusid(i) = 1 ! random walk

       end do

  end subroutine status_allocation


  subroutine set_tau()

   do i = 1, totalsites

    do j = 1,6

       siteid_settau(i,j) = settau_input(i,j) ! all ones 

     end do

   end do

  end subroutine set_tau



  subroutine randomly_insert_tracers()
   
   !f2py threadsafe
   allocate(tracerid_pos(2*numtracer,num_runs,3))
   allocate(dummyid_pos(2*numtracer,num_runs,3))
   allocate(tracerid_posn(2*numtracer,num_runs))
   allocate(dummyid_posn(2*numtracer,num_runs))
   allocate(tracerid_statusid(2*numtracer,num_runs))
   allocate(dummyid_statusid(2*numtracer,num_runs))
   allocate(tracerid_settau(2*numtracer,num_runs,6))
   allocate(dummyid_settau(2*numtracer,num_runs,6))
   allocate(tracerid_tau(2*numtracer,num_runs,6))
   allocate(dummyid_tau(2*numtracer,num_runs,6))
   allocate(tracerid_tausitesum(2*numtracer,num_runs))
   allocate(dummyid_tausitesum(2*numtracer,num_runs))
   allocate(tracerid_tautracersum(2*numtracer,num_runs))
   allocate(dummyid_tautracersum(2*numtracer,num_runs))

   CALL OMP_SET_NUM_THREADS(15)
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( i, insertflag, r3, j, temp1, temp2, temp3, k)
   !$OMP DO

  do replica = 1,num_runs
    
    do i = 1,numtracer

       call random_number(r3)
       j = floor(r3*(totalsites-2)) + 1

       call sitenumber_inverse(j,temp1, temp2, temp3)

       tracerid_pos(i,replica,1) = temp1
       tracerid_pos(i,replica,2) = temp2
       tracerid_pos(i,replica,3) = temp3
       dummyid_pos(i,replica,1) = temp1
       dummyid_pos(i,replica,2) = temp2
       dummyid_pos(i,replica,3) = temp3
       tracerid_posn(i,replica) = j

    end do

  end do

  !$OMP END DO
  !$OMP END PARALLEL

  end subroutine randomly_insert_tracers

  

  integer function sitenumber(a,b,c)


     implicit none

     integer(kind = 2),intent(in) :: a,b,c

     sitenumber = (a-1)*L*L + (b-1)*L + (c-1) + 1

     return

  end function sitenumber


  Subroutine sitenumber_inverse(a,b,c,d)

     implicit none

     integer (kind = DP_INT), intent(in) :: a

     integer(kind = 2), intent(out) :: b,c,d

     d = mod(a-1,L) + 1

     c = NINT(DBLE(mod(a-d,L**2))/DBLE(L)) + 1

     b = NINT(DBLE(mod(a-(c-1)*L-d,L**3))/DBLE(L**2)) + 1

  end subroutine sitenumber_inverse

  
  subroutine algorithm_variable_setting()
 
 do replica = 1,num_runs

  do i = 1,numtracer

     do j = 1,6

         tracerid_tau(i,replica,j) = siteid_settau(tracerid_posn(i,replica),j)

     end do

  end do

 end do

 do replica = 1,num_runs

  do i = 1,numtracer

     tracerid_tautracersum(i,replica) = 0

     do j = 1,6

        tracerid_tautracersum(i,replica) = tracerid_tautracersum(i,replica) + tracerid_tau(i,replica,j)

     end do

  end do

 end do


  allocate(taulevel(numtracer, binarylevel, num_runs))

 do replica = 1,num_runs

  do i = 1,numtracer

     do j = 1,binarylevel

        taulevel(i,j,replica) = 0

     end do

  end do

 end do


 do replica = 1,num_runs

  do j = 1,numtracer/2

     taulevel(j,1,replica) = tracerid_tautracersum(2*j,replica) + tracerid_tautracersum((2*j-1),replica)

  end do

 end do


 do replica = 1,num_runs

  do i = 2,binarylevel

     do j = 1,numtracer/(2**(i))

         taulevel(j,i,replica) = taulevel(2*j,i-1,replica) + taulevel(2*j-1,i-1,replica)

     end do

  end do

 end do


 do replica = 1,num_runs

  tautotal(replica) = taulevel(1,binarylevel,replica)

  tauadd(replica) = 0

  tausubtract(replica) = siteid_settau(tracerid_posn(1,replica),1)
 
  do k = 2,6

      tausubtract(replica) = tausubtract(replica) + siteid_settau(tracerid_posn(1,replica),k)

  end do

  add(replica) = tautotal(replica)

 end do

  end subroutine algorithm_variable_setting 


  subroutine dynamics()

    !f2py threadsafe
    allocate(status_prob(num_status,num_runs))
    allocate(average_status_prob(num_status,num_runs))
    allocate(x_prob(2*L,num_runs))
    allocate(y_prob(2*L,num_runs))
    allocate(z_prob(2*L,num_runs))
    allocate(final_prob(2*L,3,num_runs))

   do replica = 1,num_runs

    do k = 1,L

       x_prob(k,replica) = 0
       y_prob(k,replica) = 0
       z_prob(k,replica) = 0

    end do

   end do


   do replica = 1,num_runs

    do k = 1,(num_status)

       status_prob(k,replica) = 0

    end do

   end do

   do replica = 1,num_runs

    prob_counter(replica) = 0

    time_counter(replica) = 1

   end do

    start = omp_get_wtime()

  CALL OMP_SET_NUM_THREADS(15)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, k, tracersite, pointvalue, randomflag, tausearch, j, point, direction, swap1, swap2, swap3, r1, filenumber)
  !$OMP DO
  
  do replica = 1,num_runs
    
    tracersite = 1
    
    do i = 2,totalmoves
                   
       do k = 1,6

          tracerid_tau(tracersite,replica,k) = siteid_settau(tracerid_posn(tracersite,replica),k)
          
       end do

       tauadd(replica) = 0
       
       do k =1,6
             
           tauadd(replica) = tauadd(replica) + siteid_settau(tracerid_posn(tracersite,replica),k)
             
       end do

       tautotal(replica) = add(replica) + tauadd(replica) - tausubtract(replica)

       if (tautotal(replica) .eq. 0) then

           print*, i, replica

           STOP 'TAUTOTAL IS ZERO'

       end if

       tracerid_tautracersum(tracersite,replica) = tracerid_tautracersum(tracersite,replica) + tauadd(replica)&
 - tausubtract(replica)

       do k = 1,binarylevel

          pointvalue = dble(tracersite)/(2**(k))
          taulevel(ceiling(pointvalue),k,replica) = taulevel(ceiling(pointvalue),k,replica) + tauadd(replica) - tausubtract(replica)
          
       end do

       randomflag = .True.

       do while (randomflag)
         
           call random_number(r1)

           k = 1

           tausearch = r1

           do j = (binarylevel - 1), 1, (-1)

               if (tausearch .le. (taulevel(2*k-1,j,replica)/tautotal(replica))) then

                   point = 2*k - 1

               else

                   point = 2*k

                   tausearch = tausearch - taulevel(2*k-1,j,replica)/tautotal(replica)

               end if

               k = point
           
           end do

           if (tausearch .le. (tracerid_tautracersum((2*k-1),replica)/tautotal(replica)))then

               point = 2*k - 1

           else

               point = 2*k

               tausearch = tausearch - (tracerid_tautracersum((2*k-1),replica)/tautotal(replica))
            
           end if

           
           tracersite = point

           tausum(replica) = 0

           direction = 0

           do k = 1,6

              tausum(replica) = tausum(replica) + (tracerid_tau(tracersite,replica,k))/tautotal(replica)

              if (tausum(replica) .ge. tausearch) then

                 direction = k
                 
                 exit

              end if
              
          end do

          if ((tracerid_tau(tracersite,replica,direction) .ne.0) .and. (direction .ne. 0)) then

               randomflag = .False.

          end if
       
       end do

       tausubtract(replica) = 0 

       do k = 1,6

           tausubtract(replica) = tausubtract(replica) + siteid_settau(tracerid_posn(tracersite,replica),k)

       end do        

       add(replica) = tautotal(replica)
      
       call move_posn(direction, tracersite, replica, swap1, swap2, swap3)
       
       tracerid_pos(tracersite,replica,(direction/2) + mod(direction,2)) = swap1
       dummyid_pos(tracersite,replica,(direction/2) + mod(direction,2)) = swap2
       tracerid_posn(tracersite,replica) = swap3
      

    if (mod(i,msd_freq) .eq. 0) then
       
       time_counter(replica) = time_counter(replica) + 1

       if (time_counter(replica) .eq. 2) then
         
          filenumber = replica + 10

          OPEN(UNIT = filenumber, FILE = trajectory_file(replica), STATUS = 'old', ACTION = 'write', POSITION = 'append')

          WRITE(filenumber,*) time_counter(replica), tautotal(replica)

          do j = 1,numtracer
          
             WRITE(filenumber,*) dummyid_pos(j,replica,1), dummyid_pos(j,replica,2), dummyid_pos(j,replica,3)

          end do

       else

          WRITE(filenumber,*) time_counter(replica), tautotal(replica)

          do j = 1,numtracer

              WRITE(filenumber,*) dummyid_pos(j,replica,1), dummyid_pos(j,replica,2), dummyid_pos(j,replica,3)

          end do

       end if

       if (i .eq. totalmoves) then

          CLOSE(filenumber)

       end if

    end if


    if (mod(i,print_freq) .eq. 0) then

       print*, i

    end if


   end do
  
 end do

 !$OMP END DO
 !$OMP END PARALLEL

  finish = omp_get_wtime()
  

  end subroutine dynamics

 

  subroutine move_posn(a, b, c, d, e, f)
    
     implicit none
     
     integer, intent(in) :: a,b,c

     integer(kind = 2), intent(out) :: d,e

     integer, intent(out) :: f
 
     d = (tracerid_pos(b,c,a/2 + mod(a,2))) + (2*(mod(a,2))-1)*delta

     e = (dummyid_pos(b,c,a/2 + mod(a,2))) + (2*(mod(a,2))-1)*delta


     if (mod(a,2) .eq. 1) then

         if (d .gt. L) then

              d = d - L

         end if

         if ((a .eq. 1) .or.  (a .eq. 2)) then

             f = sitenumber(d, tracerid_pos(b,c,2), tracerid_pos(b,c,3))

         else if ((a .eq. 3) .or.  (a .eq. 4)) then

             f = sitenumber(tracerid_pos(b,c,1), d, tracerid_pos(b,c,3))

         else

             f = sitenumber(tracerid_pos(b,c,1), tracerid_pos(b,c,2), d)

         end if

      else

          if (d .lt. 1) then

              d = d + L

          end if


         if ((a .eq. 1) .or.  (a .eq. 2)) then

             f = sitenumber(d, tracerid_pos(b,c,2), tracerid_pos(b,c,3))

         else if ((a .eq. 3) .or.  (a .eq. 4)) then

             f = sitenumber(tracerid_pos(b,c,1), d, tracerid_pos(b,c,3))

         else

             f = sitenumber(tracerid_pos(b,c,1), tracerid_pos(b,c,2), d)

         end if

      end if

  end subroutine move_posn 

  end module kmc_RW 
