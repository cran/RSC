   subroutine cormadvecdp(matrix,nrow,ncol,res,ressize,correcteven)

      use selectionalgo

      implicit none
      
      !note kind dp is defined in selectionalgo
      integer, intent(in) :: nrow, ncol,ressize, correcteven     !.Fortran R function will not deal with deferred size arrays
      real(kind=dp), dimension(nrow,ncol), intent(inout) :: matrix
      real(kind=dp), dimension(ressize), intent(out) :: res

      integer :: i, j, n, p
      real(kind=dp)  :: med,mad, fresh
      real(kind=dp), dimension(:), allocatable :: U, V, A, B


      p=ubound(matrix,2)
      n=ubound(matrix,1)
      allocate(U(n),V(n))

      if (correcteven==1) then
         do i=1,p-3,3
            U=matrix(:,i)
            med=iselect(U,evencorrection=.true.)
            matrix(:,i)=(U-med)/(sqrt2*cost*iselect(abs(U-med),evencorrection=.true.))
            U=matrix(:,i+1)
            med=iselect(U,evencorrection=.true.)
            matrix(:,i+1)=(U-med)/(sqrt2*cost*iselect(abs(U-med),evencorrection=.true.))
            U=matrix(:,i+2)
            med=iselect(U,evencorrection=.true.)
            matrix(:,i+2)=(U-med)/(sqrt2*cost*iselect(abs(U-med),evencorrection=.true.))
         end do
         do i=i,p
            U=matrix(:,i)
            med=iselect(U,evencorrection=.true.)
            matrix(:,i)=(U-med)/(sqrt2*cost*iselect(abs(U-med),evencorrection=.true.))
         end do
         !unrolled loops
         n=1
         do i=1,p-1
            do j=i+1,p
               U=matrix(:,i)+matrix(:,j)
               V=-matrix(:,i)+matrix(:,j)
               mad=(cost*iselect(abs(U-iselect(U,evencorrection=.true.)),evencorrection=.true.))**2
               med=(cost*iselect(abs(V-iselect(V,evencorrection=.true.)),evencorrection=.true.))**2
               res(n)=(mad-med)/(mad+med)
               n=n+1
            end do 
         end do
      else
         do i=1,p-3,3
            U=matrix(:,i)
            med=iselect(U)
            matrix(:,i)=(U-med)/(sqrt2*cost*iselect(abs(U-med)))
            U=matrix(:,i+1)
            med=iselect(U)
            matrix(:,i+1)=(U-med)/(sqrt2*cost*iselect(abs(U-med)))
            U=matrix(:,i+2)
            med=iselect(U)
            matrix(:,i+2)=(U-med)/(sqrt2*cost*iselect(abs(U-med)))
         end do
         do i=i,p
            U=matrix(:,i)
            med=iselect(U)
            matrix(:,i)=(U-med)/(sqrt2*cost*iselect(abs(U-med)))
         end do
         !unrolled loops
         n=1
         do i=1,p-1
            do j=i+1,p
               U=matrix(:,i)+matrix(:,j)
               V=-matrix(:,i)+matrix(:,j)
               mad=(cost*iselect(abs(U-iselect(U))))**2
               med=(cost*iselect(abs(V-iselect(V))))**2
               res(n)=(mad-med)/(mad+med)
               n=n+1
            end do 
         end do
      end if


   end subroutine cormadvecdp
