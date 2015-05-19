#include "cppdefs.h"

module output_manager_core

   use field_manager
   use time, only: calendar_date

   implicit none

   public type_used_field, type_file, write_time_string, output_manager_fatal_error

   private

   integer,parameter,public :: max_path = 256

   integer,parameter,public :: time_method_none          = 0  ! time-independent variable
   integer,parameter,public :: time_method_instantaneous = 1
   integer,parameter,public :: time_method_mean          = 2
   integer,parameter,public :: time_method_integrated    = 3

   integer,parameter,public :: time_unit_none   = 0
   integer,parameter,public :: time_unit_second = 1
   integer,parameter,public :: time_unit_day    = 2
   integer,parameter,public :: time_unit_month  = 3
   integer,parameter,public :: time_unit_year   = 4

   integer,parameter,public :: rk = kind(_ONE_)

   type type_used_field
      character(len=string_length)      :: output_name = ''
      type (type_field),pointer         :: source      => null()
      integer                           :: time_method = time_method_instantaneous
      real(rk)                          :: work_0d
      real(rk),allocatable,dimension(:) :: work_1d
      integer                           :: ncid        = -1
      type (type_used_field),pointer    :: next        => null()
   end type type_used_field

   type type_file
      type (type_field_manager),pointer :: field_manager
      character(len=max_path)        :: path          = ''
      integer                        :: time_unit     = time_unit_none
      integer                        :: time_step     = 0
      integer                        :: n             = 0  ! Number of model time steps processed so far for next output
      integer                        :: next_julian   = -1
      integer                        :: next_seconds  = -1
      type (type_used_field),pointer :: first_field   => null()
      class (type_file),pointer      :: next          => null()
   contains
      procedure :: initialize
      procedure :: save
      procedure :: finalize
   end type type_file

contains

   subroutine initialize(self,julianday,secondsofday)
      class (type_file),intent(inout) :: self
      integer,          intent(in)    :: julianday,secondsofday
      stop 'output_manager_core:initialize not implemented'
   end subroutine

   subroutine save(self,julianday,secondsofday)
      class (type_file),intent(inout) :: self
      integer,          intent(in)    :: julianday,secondsofday
      stop 'output_manager_core:save not implemented'
   end subroutine

   subroutine finalize(file)
      class (type_file),intent(inout) :: file
   end subroutine
   
   subroutine write_time_string(jul,secs,timestr)
      integer,         intent(in)  :: jul,secs
      character(len=*),intent(out) :: timestr

      integer :: ss,min,hh,dd,mm,yy

      hh   = secs/3600
      min  = (secs-hh*3600)/60
      ss   = secs - 3600*hh - 60*min

      call calendar_date(jul,yy,mm,dd)

      write(timestr,'(i4.4,a1,i2.2,a1,i2.2,1x,i2.2,a1,i2.2,a1,i2.2)')  &
                           yy,'-',mm,'-',dd,hh,':',min,':',ss
   end subroutine write_time_string
   
   subroutine output_manager_fatal_error(location,error)
      character(len=*),intent(in) :: location,error
      
      FATAL trim(location)//': '//trim(error)
      stop 'output_manager::output_manager_fatal_error'
   end subroutine

end module output_manager_core
