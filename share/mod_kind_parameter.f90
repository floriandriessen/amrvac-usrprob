!> Numerical storage size parameters for real and integer values.
!> Adopted from https://fortran-lang.org/learn/best_practices/floating_point/
module mod_kind_parameter

   implicit none
   public

   !> Single precision real numbers for 6 digits (32 bits)
   !> 6 digits; range 10^(-37) to 10^(37)-1
   integer, parameter :: sp = selected_real_kind(6, 37)

   !> Double precision real numbers (64 bits)
   !> 15 digits: range 10^(-307) to 10^(307)-1
   integer, parameter :: dp = selected_real_kind(15, 307)

   !> Quadruple precision real numbers (128 bits)
   !> 33 digits: range 10^(-4931) to 10^(4931)-1
   integer, parameter :: qp = selected_real_kind(33, 4931)

   !> Char length for integers (8 bits) 
   !> range -2^(7) to 2^(7)-1
   integer, parameter :: i1 = selected_int_kind(2)

   !> Short length for integers (16 bits) 
   !> range -2^(15) to 2^(15)-1
   integer, parameter :: i2 = selected_int_kind(4)

   !> Length of default integers (32 bits)
   !> range -2^(31) to 2^(31)-1
   integer, parameter :: i4 = selected_int_kind(9)

   !> Long length for integers (64 bits)
   !> range -2^(63) to 2^(63)-1
   integer, parameter :: i8 = selected_int_kind(18)

end module mod_kind_parameter