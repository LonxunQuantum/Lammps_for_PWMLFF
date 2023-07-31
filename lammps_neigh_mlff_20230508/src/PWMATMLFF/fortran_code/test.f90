program main
   use li_ff_mod
   implicit none

   integer(4) :: ff_idx, slen
   real(8) :: ocut
   character(300) :: name_ptr(1)

   ! Initialize variables
   ff_idx = 1
   slen = 20
   name_ptr(1) = "./myforcefield.ff"

   ! Call the dp_ff_load subroutine
   call li_ff12_load(name_ptr, ff_idx, slen, ocut)

   ! ... Your other code ...

   ! Call the dp_ff_deallocate subroutine
!    call li_ff_deallocate(ff_idx)

end program main
