module selectionalgo 
   use, intrinsic :: iso_fortran_env

   implicit none

   real, parameter :: cost=1.4826, sqrt2=sqrt(2.) 
   integer, parameter :: sp = real32, dp=real64

   interface qselect
      procedure quickselectsp, quickselectscalarsp,quickselectdp, quickselectscalardp
   end interface qselect

   interface iselect
      procedure introselectsp, introselectscalarsp, introselectdp, introselectscalardp
   end interface iselect

   contains

      !!!!!!! SINGLE PRECISION
   !from here there are the selectionalgo and selectionalgoimpl condensed in one file
      recursive subroutine quickselectrecursivesp(invector,tovector,k,output)
         real(kind=sp), dimension(:), allocatable, target, intent(inout) :: invector
         real(kind=sp),  intent(inout), pointer, contiguous :: tovector(:)
         integer, intent(inout) :: k
         real(kind=sp), intent(out) :: output
         real(kind=sp) :: swapper, pvt
         integer :: i, r, subst
         real(kind=sp), pointer :: a,b,c

         r=size(tovector)

         select case(r)
         case(1)
            output=tovector(1)
            return
         case(2)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(2)
                  output=maxval(tovector)
               end select
            return
         case(3)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(3)
                  output=maxval(tovector)
               case(2)
               i=r/2+1
               output=sum(tovector)-minval(tovector)-maxval(tovector)
            end select
            return
         case default
            i=r/2+1

            !pivoting section (pivot of 3)
            a=>tovector(1)
            b=>tovector(i)
            c=>tovector(r)
            pvt=a+b+c
            swapper=max(a,b,c)
            b=swapper
            pvt=pvt-swapper
            swapper=min(a,b,c)
            a=swapper
            pvt=pvt-swapper
            c=pvt

            !one pass section
            subst=1
            b=>tovector(subst)
            do i=1,(r-1-3),3
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+1)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+2)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            do i=i,(r-1)
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            tovector(r)=tovector(subst)
            tovector(subst)=pvt

            !decide next step
            select case(subst-k)
            case(0)
               output=tovector(k)
            case(1:)
               tovector => tovector(1:subst-1)
               call quickselectrecursivesp(invector,tovector,k, output)
            case(:-1)
               k=k-subst
               tovector => tovector(subst+1:)
               call quickselectrecursivesp(invector,tovector,k, output)
            end select
         end select
      end subroutine quickselectrecursivesp

      recursive subroutine introselectrecursivesp(invector,tovector,k,output)
         real(kind=sp), dimension(:), allocatable, target, intent(inout) :: invector
         real(kind=sp),  intent(inout), pointer, contiguous :: tovector(:)
         integer, intent(inout) :: k
         real(kind=sp), intent(out) :: output
         real(kind=sp) :: swapper, pvt
         integer, save :: switch
         integer :: i, r, subst, r_min, switch_after
         real(kind=sp), pointer :: a,b,c

         r=ubound(tovector,1)
         r_min = 3000 ! at least this bigger to switch to median of medians
         switch_after = 5 ! at least fails this amount of time to reorder
                          ! one third of the vector to consider med of med

         select case(r)
         case(1)
            output=tovector(1)
            switch=0
            return
         case(2)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(2)
                  output=maxval(tovector)
            end select
            switch=0
            return
         case(3)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(3)
                  output=maxval(tovector)
               case(2)
               i=r/2+1
               output=sum(tovector)-minval(tovector)-maxval(tovector)
            end select
            switch=0
            return
         case default
            i=r/2+1

            !pivoting section (pivot of 3)
            a=>tovector(1)
            b=>tovector(i)
            c=>tovector(r)
            pvt=a+b+c
            swapper=max(a,b,c)
            b=swapper
            pvt=pvt-swapper
            swapper=min(a,b,c)
            a=swapper
            pvt=pvt-swapper
            c=pvt

            !one pass section
            subst=1
            b=>tovector(subst)
            do i=1,(r-1-3),3
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+1)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+2)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            do i=i,(r-1)
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            tovector(r)=tovector(subst)
            tovector(subst)=pvt

            !decide next step
            select case(subst-k)
            case(0)
               output=tovector(k)
               switch=0
            case(1:)
               tovector => tovector(1:subst-1)
               if ((subst-1)>r*2/3) then
                  switch=switch+1
                  if (switch==switch_after .and. r>r_min) then
                     !reset switch and call medofmedsp
                     switch=0 
                     call medofmedsp(invector,tovector,k, output)
                  else
                     call introselectrecursivesp(invector,tovector,k, output)
                  end if
               else
                  switch=0
                  call introselectrecursivesp(invector,tovector,k, output)
               end if
            case(:-1)
               k=k-subst
               tovector => tovector(subst+1:)
               if ((r-subst)>r*2/3) then
                  switch=switch+1
                  if (switch==switch_after .and. r>r_min) then
                     !reset switch and call medofmedsp
                     switch=0 
                     call medofmedsp(invector,tovector,k, output)
                  else
                     call introselectrecursivesp(invector,tovector,k, output)
                  end if
               else
                  switch=0
                  call introselectrecursivesp(invector,tovector,k, output)
               end if
            end select
         end select
      end subroutine introselectrecursivesp

      recursive subroutine medofmedsp(invector,tovector,k,output)
         real(kind=sp), dimension(:), allocatable, target, intent(inout) :: invector
         real(kind=sp),  intent(inout), pointer, contiguous :: tovector(:)
         real(kind=sp), intent(out) :: output
         integer, intent(inout) :: k

         real(kind=sp), pointer, contiguous :: subsec(:)
         real(kind=sp),  dimension(:),allocatable :: meds
         integer :: i,subst,r
         real(kind=sp) ::  pvt, swapper
         real(kind=sp), pointer :: a,b

         r=ubound(tovector,1)
         select case(mod(r,5))
            case(0)
               allocate(meds(r/5))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
            case(1)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               meds(r/5+1)=tovector(i)
            case(2)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               subsec=>tovector(i:)
               subst=minloc(subsec,1)
               pvt=subsec(subst)
               subsec(subst)=subsec(1)
               subsec(1)=pvt
               meds(r/5+1)=subsec(2)
            case(3)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               subsec=>tovector(i:)
               subst=minloc(subsec,1)
               pvt=subsec(subst)
               subsec(subst)=subsec(1)
               subsec(1)=pvt
               subst=minloc(subsec(2:),1)+1
               pvt=subsec(subst)
               subsec(subst)=subsec(2)
               subsec(2)=pvt
               meds(r/5+1)=pvt
            case(4)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               subsec=>tovector(i:)
               subst=minloc(subsec,1)
               pvt=subsec(subst)
               subsec(subst)=subsec(1)
               subsec(1)=pvt
               subst=minloc(subsec(2:),1)+1
               pvt=subsec(subst)
               subsec(subst)=subsec(2)
               subsec(2)=pvt
               subst=minloc(subsec(3:),1)+2
               pvt=subsec(subst)
               subsec(subst)=subsec(3)
               subsec(3)=pvt
               meds(r/5+1)=pvt
            end select

            !set pivot
            pvt=qselect(meds)
            do i=1,r-3,3
               if (tovector(i)==pvt) then
                  tovector(i)=tovector(r)
                  tovector(r)=pvt
                  exit
               else if (tovector(i+1)==pvt) then
                  tovector(i+1)=tovector(r)
                  tovector(r)=pvt
                  exit
               else if (tovector(i+2)==pvt) then
                  tovector(i+2)=tovector(r)
                  tovector(r)=pvt
                  exit
               else 
                  cycle
               end if
            end do
            do i=i,r
               if (tovector(i)==pvt) then
                  tovector(i)=tovector(r)
                  tovector(r)=pvt
                  exit
               else 
                  cycle
               end if
            end do

            call introselectrecursivesp(invector,tovector,k, output)
      end subroutine medofmedsp

      function introselectsp(invector,ord,evencorrection)
         real(kind=sp), intent(in), dimension(:) :: invector
         logical, optional, intent(in) :: evencorrection
         integer, intent(in), optional :: ord
         real(kind=sp) :: introselectsp

         real(kind=sp), dimension(:), allocatable, target :: vector
         real(kind=sp), pointer, contiguous :: tovector(:) => null()
         integer ::  k1, k

         k1=size(invector)
         allocate(vector(k1))
         vector=invector
         tovector => vector

         k1=ubound(invector,1)
         if(present(ord)) then
            k=ord
         else
            k=k1/2+1
         end if

         if(present(evencorrection)) then
            if (k/=k1/2+1 .or. mod(k1,2)/=0) then
               !write(*,'(a)') 'Warning: evencorrection argument ignored.'
               k1=0
            else 
               k1=0
               if (evencorrection) k1=k-1
            end if
         else
            k1=0
         end if

         call introselectrecursivesp(vector,tovector,k, introselectsp)

         if (k1/=0) then
            block
               real(kind=sp) :: temp
               temp=introselectsp
               tovector=>vector
               call introselectrecursivesp(vector,tovector,k1, introselectsp)
               introselectsp=(introselectsp+temp)/2.0
            end block
         end if
      end function introselectsp

      function quickselectsp(invector,ord,evencorrection)
         real(kind=sp), intent(in), dimension(:) :: invector
         logical, optional, intent(in) :: evencorrection
         integer, intent(in), optional :: ord
         real(kind=sp) :: quickselectsp

         real(kind=sp), dimension(:), allocatable, target :: vector
         real(kind=sp), pointer, contiguous :: tovector(:) => null()
         integer ::  k1, k

         k1=size(invector)
         allocate(vector(k1))
         vector=invector
         tovector => vector

         k1=ubound(invector,1)
         if(present(ord)) then
            k=ord
         else
            k=k1/2+1
         end if

         if(present(evencorrection)) then
            if (k/=k1/2+1 .or. mod(k1,2)/=0) then
               !write(*,'(a)') 'Warning: evencorrection argument ignored.'
               k1=0
            else 
               k1=0
               if (evencorrection) k1=k-1
            end if
         else
            k1=0
         end if

         call quickselectrecursivesp(vector,tovector,k, quickselectsp)

         if (k1/=0) then
            block
               real(kind=sp) :: temp
               temp=quickselectsp
               tovector=>vector(1:k1+1)
               call quickselectrecursivesp(vector,tovector,k1, quickselectsp)
               quickselectsp=(quickselectsp+temp)/2.0
            end block
         end if
      end function quickselectsp

      !deal with scalar input
      function introselectscalarsp(invector,ord, evencorrection) 
         real(kind=sp), intent(in) :: invector
         integer, intent(in), optional :: ord
         logical, optional, intent(in) :: evencorrection
         real(kind=sp) :: introselectscalarsp

         introselectscalarsp=invector
      end function introselectscalarsp

      function quickselectscalarsp(invector,ord, evencorrection) 
         real(kind=sp), intent(in) :: invector
         integer, intent(in), optional :: ord
         logical, optional, intent(in) :: evencorrection
         real(kind=sp) :: quickselectscalarsp

         quickselectscalarsp=invector
      end function quickselectscalarsp

      !!!!!!! DOUBLE PRECISION
   !from here there are the selectionalgo and selectionalgoimpl condensed in one file
      recursive subroutine quickselectrecursivedp(invector,tovector,k,output)
         real(kind=dp), dimension(:), allocatable, target, intent(inout) :: invector
         real(kind=dp),  intent(inout), pointer, contiguous :: tovector(:)
         integer, intent(inout) :: k
         real(kind=dp), intent(out) :: output
         real(kind=dp) :: swapper, pvt
         integer :: i, r, subst
         real(kind=dp), pointer :: a,b,c

         r=size(tovector)

         select case(r)
         case(1)
            output=tovector(1)
            return
         case(2)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(2)
                  output=maxval(tovector)
               end select
            return
         case(3)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(3)
                  output=maxval(tovector)
               case(2)
               i=r/2+1
               output=sum(tovector)-minval(tovector)-maxval(tovector)
            end select
            return
         case default
            i=r/2+1

            !pivoting section (pivot of 3)
            a=>tovector(1)
            b=>tovector(i)
            c=>tovector(r)
            pvt=a+b+c
            swapper=max(a,b,c)
            b=swapper
            pvt=pvt-swapper
            swapper=min(a,b,c)
            a=swapper
            pvt=pvt-swapper
            c=pvt

            !one pass section
            subst=1
            b=>tovector(subst)
            do i=1,(r-1-3),3
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+1)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+2)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            do i=i,(r-1)
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            tovector(r)=tovector(subst)
            tovector(subst)=pvt

            !decide next step
            select case(subst-k)
            case(0)
               output=tovector(k)
            case(1:)
               tovector => tovector(1:subst-1)
               call quickselectrecursivedp(invector,tovector,k, output)
            case(:-1)
               k=k-subst
               tovector => tovector(subst+1:)
               call quickselectrecursivedp(invector,tovector,k, output)
            end select
         end select
      end subroutine quickselectrecursivedp

      recursive subroutine introselectrecursivedp(invector,tovector,k,output)
         real(kind=dp), dimension(:), allocatable, target, intent(inout) :: invector
         real(kind=dp),  intent(inout), pointer, contiguous :: tovector(:)
         integer, intent(inout) :: k
         real(kind=dp), intent(out) :: output
         real(kind=dp) :: swapper, pvt
         integer, save :: switch
         integer :: i, r, subst, r_min, switch_after
         real(kind=dp), pointer :: a,b,c

         r=ubound(tovector,1)
         r_min = 3000 ! at least this bigger to switch to median of medians
         switch_after = 5 ! at least fails this amount of time to reorder
                          ! one third of the vector to consider med of med

         select case(r)
         case(1)
            output=tovector(1)
            switch=0
            return
         case(2)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(2)
                  output=maxval(tovector)
            end select
            switch=0
            return
         case(3)
            select case(k)
               case(1)
                  output=minval(tovector)
               case(3)
                  output=maxval(tovector)
               case(2)
               i=r/2+1
               output=sum(tovector)-minval(tovector)-maxval(tovector)
            end select
            switch=0
            return
         case default
            i=r/2+1

            !pivoting section (pivot of 3)
            a=>tovector(1)
            b=>tovector(i)
            c=>tovector(r)
            pvt=a+b+c
            swapper=max(a,b,c)
            b=swapper
            pvt=pvt-swapper
            swapper=min(a,b,c)
            a=swapper
            pvt=pvt-swapper
            c=pvt

            !one pass section
            subst=1
            b=>tovector(subst)
            do i=1,(r-1-3),3
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+1)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
               a => tovector(i+2)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            do i=i,(r-1)
               a => tovector(i)
               if (a<pvt) then
                     swapper=b
                     b=a
                     a=swapper
                  subst=subst+1
                  b=>tovector(subst)
               end if
            end do
            tovector(r)=tovector(subst)
            tovector(subst)=pvt

            !decide next step
            select case(subst-k)
            case(0)
               output=tovector(k)
               switch=0
            case(1:)
               tovector => tovector(1:subst-1)
               if ((subst-1)>r*2/3) then
                  switch=switch+1
                  if (switch==switch_after .and. r>r_min) then
                     !reset switch and call medofmeddp
                     switch=0 
                     call medofmeddp(invector,tovector,k, output)
                  else
                     call introselectrecursivedp(invector,tovector,k, output)
                  end if
               else
                  switch=0
                  call introselectrecursivedp(invector,tovector,k, output)
               end if
            case(:-1)
               k=k-subst
               tovector => tovector(subst+1:)
               if ((r-subst)>r*2/3) then
                  switch=switch+1
                  if (switch==switch_after .and. r>r_min) then
                     !reset switch and call medofmeddp
                     switch=0 
                     call medofmeddp(invector,tovector,k, output)
                  else
                     call introselectrecursivedp(invector,tovector,k, output)
                  end if
               else
                  switch=0
                  call introselectrecursivedp(invector,tovector,k, output)
               end if
            end select
         end select
      end subroutine introselectrecursivedp

      recursive subroutine medofmeddp(invector,tovector,k,output)
         real(kind=dp), dimension(:), allocatable, target, intent(inout) :: invector
         real(kind=dp),  intent(inout), pointer, contiguous :: tovector(:)
         real(kind=dp), intent(out) :: output
         integer, intent(inout) :: k

         real(kind=dp), pointer, contiguous :: subsec(:)
         real(kind=dp),  dimension(:),allocatable :: meds
         integer :: i,subst,r
         real(kind=dp) ::  pvt, swapper
         real(kind=dp), pointer :: a,b

         r=ubound(tovector,1)
         select case(mod(r,5))
            case(0)
               allocate(meds(r/5))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
            case(1)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               meds(r/5+1)=tovector(i)
            case(2)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               subsec=>tovector(i:)
               subst=minloc(subsec,1)
               pvt=subsec(subst)
               subsec(subst)=subsec(1)
               subsec(1)=pvt
               meds(r/5+1)=subsec(2)
            case(3)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               subsec=>tovector(i:)
               subst=minloc(subsec,1)
               pvt=subsec(subst)
               subsec(subst)=subsec(1)
               subsec(1)=pvt
               subst=minloc(subsec(2:),1)+1
               pvt=subsec(subst)
               subsec(subst)=subsec(2)
               subsec(2)=pvt
               meds(r/5+1)=pvt
            case(4)
               allocate(meds(r/5+1))
               do i=1,r-4,5
                  subsec=>tovector(i:i+4)
                  subst=minloc(subsec,1)
                  pvt=subsec(subst)
                  subsec(subst)=subsec(1)
                  subsec(1)=pvt
                  subst=minloc(subsec(2:),1)+1
                  pvt=subsec(subst)
                  subsec(subst)=subsec(2)
                  subsec(2)=pvt
                  subst=minloc(subsec(3:),1)+2
                  pvt=subsec(subst)
                  subsec(subst)=subsec(3)
                  subsec(3)=pvt
                  meds(i/5+1)=pvt
               end do
               subsec=>tovector(i:)
               subst=minloc(subsec,1)
               pvt=subsec(subst)
               subsec(subst)=subsec(1)
               subsec(1)=pvt
               subst=minloc(subsec(2:),1)+1
               pvt=subsec(subst)
               subsec(subst)=subsec(2)
               subsec(2)=pvt
               subst=minloc(subsec(3:),1)+2
               pvt=subsec(subst)
               subsec(subst)=subsec(3)
               subsec(3)=pvt
               meds(r/5+1)=pvt
            end select

            !set pivot
            pvt=qselect(meds)
            do i=1,r-3,3
               if (tovector(i)==pvt) then
                  tovector(i)=tovector(r)
                  tovector(r)=pvt
                  exit
               else if (tovector(i+1)==pvt) then
                  tovector(i+1)=tovector(r)
                  tovector(r)=pvt
                  exit
               else if (tovector(i+2)==pvt) then
                  tovector(i+2)=tovector(r)
                  tovector(r)=pvt
                  exit
               else 
                  cycle
               end if
            end do
            do i=i,r
               if (tovector(i)==pvt) then
                  tovector(i)=tovector(r)
                  tovector(r)=pvt
                  exit
               else 
                  cycle
               end if
            end do

            call introselectrecursivedp(invector,tovector,k, output)
      end subroutine medofmeddp

      function introselectdp(invector,ord,evencorrection)
         real(kind=dp), intent(in), dimension(:) :: invector
         logical, optional, intent(in) :: evencorrection
         integer, intent(in), optional :: ord
         real(kind=dp) :: introselectdp

         real(kind=dp), dimension(:), allocatable, target :: vector
         real(kind=dp), pointer, contiguous :: tovector(:) => null()
         integer ::  k1, k

         k1=size(invector)
         allocate(vector(k1))
         vector=invector
         tovector => vector

         k1=ubound(invector,1)
         if(present(ord)) then
            k=ord
         else
            k=k1/2+1
         end if

         if(present(evencorrection)) then
            if (k/=k1/2+1 .or. mod(k1,2)/=0) then
               !write(*,'(a)') 'Warning: evencorrection argument ignored.'
               k1=0
            else 
               k1=0
               if (evencorrection) k1=k-1
            end if
         else
            k1=0
         end if

         call introselectrecursivedp(vector,tovector,k, introselectdp)

         if (k1/=0) then
            block
               real(kind=dp) :: temp
               temp=introselectdp
               !tovector=>vector(1:k1+1)
               tovector=>vector
               call introselectrecursivedp(vector,tovector,k1, introselectdp)
               introselectdp=(introselectdp+temp)/2.0
            end block
         end if
      end function introselectdp

      function quickselectdp(invector,ord,evencorrection)
         real(kind=dp), intent(in), dimension(:) :: invector
         logical, optional, intent(in) :: evencorrection
         integer, intent(in), optional :: ord
         real(kind=dp) :: quickselectdp

         real(kind=dp), dimension(:), allocatable, target :: vector
         real(kind=dp), pointer, contiguous :: tovector(:) => null()
         integer ::  k1, k

         k1=size(invector)
         allocate(vector(k1))
         vector=invector
         tovector => vector

         k1=ubound(invector,1)
         if(present(ord)) then
            k=ord
         else
            k=k1/2+1
         end if

         if(present(evencorrection)) then
            if (k/=k1/2+1 .or. mod(k1,2)/=0) then
               !write(*,'(a)') 'Warning: evencorrection argument ignored.'
               k1=0
            else 
               k1=0
               if (evencorrection) k1=k-1
            end if
         else
            k1=0
         end if

         call quickselectrecursivedp(vector,tovector,k, quickselectdp)

         if (k1/=0) then
            block
               real(kind=dp) :: temp
               temp=quickselectdp
               tovector=>vector(1:k1+1)
               call quickselectrecursivedp(vector,tovector,k1, quickselectdp)
               quickselectdp=(quickselectdp+temp)/2.0
            end block
         end if
      end function quickselectdp

      !deal with scalar input
      function introselectscalardp(invector,ord, evencorrection) 
         real(kind=dp), intent(in) :: invector
         integer, intent(in), optional :: ord
         logical, optional, intent(in) :: evencorrection
         real(kind=dp) :: introselectscalardp

         introselectscalardp=invector
      end function introselectscalardp

      function quickselectscalardp(invector,ord, evencorrection) 
         real(kind=dp), intent(in) :: invector
         integer, intent(in), optional :: ord
         logical, optional, intent(in) :: evencorrection
         real(kind=dp) :: quickselectscalardp

         quickselectscalardp=invector
      end function quickselectscalardp


end module selectionalgo
