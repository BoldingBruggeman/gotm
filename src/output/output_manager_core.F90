#include "cppdefs.h"

module output_manager_core

   use field_manager

   implicit none

   public type_output_item,type_output_category,type_output_field, type_file, write_time_string, output_manager_fatal_error, host, type_host

   private

   integer,parameter,public :: max_path = 256

   integer,parameter,public :: time_method_none          = 0  ! time-independent variable
   integer,parameter,public :: time_method_instantaneous = 1
   integer,parameter,public :: time_method_mean          = 2
   integer,parameter,public :: time_method_integrated    = 3

   integer,parameter,public :: time_unit_none   = 0
   integer,parameter,public :: time_unit_second = 1
   integer,parameter,public :: time_unit_hour   = 2
   integer,parameter,public :: time_unit_day    = 3
   integer,parameter,public :: time_unit_month  = 4
   integer,parameter,public :: time_unit_year   = 5

   integer,parameter,public :: rk = kind(_ONE_)

   type type_host
   contains
      procedure :: julian_day
      procedure :: calendar_date
   end type

   type type_output_item
      integer :: time_method = time_method_instantaneous
   end type

   type,extends(type_output_item) ::  type_output_category
      character(len=string_length)         :: name = ''
      character(len=string_length)         :: prefix = ''
      character(len=string_length)         :: postfix = ''
      integer                              :: output_level = output_level_default
      class (type_category_node),   pointer :: source => null()
      class (type_output_category),pointer :: next => null()
   end type
   
   type,extends(type_output_item) :: type_output_field
      character(len=string_length) :: output_name = ''
      type (type_field),pointer    :: source => null()

      ! Work arrays (only allocated/used if storing non-instantaneous data)
      real(rk)                          :: work_0d
      real(rk),allocatable              :: work_1d(:)
      real(rk),allocatable              :: work_2d(:,:)
      real(rk),allocatable              :: work_3d(:,:,:)

      ! Pointers to data to store (either pointing to instantaneous data, or to the above work arrays)
      real(rk),pointer                  :: data_0d        => null()
      real(rk),pointer                  :: data_1d(:)     => null()
      real(rk),pointer                  :: data_2d(:,:)   => null()
      real(rk),pointer                  :: data_3d(:,:,:) => null()

      class (type_output_field),pointer :: next => null()
   end type type_output_field

   type type_file
      type (type_field_manager),pointer :: field_manager => null()
      character(len=max_path)        :: path          = ''
      integer                        :: time_unit     = time_unit_none
      integer                        :: time_step     = 0
      integer                        :: n             = 0  ! Number of model time steps processed so far for next output
      integer                        :: next_julian   = -1
      integer                        :: next_seconds  = -1
      class (type_output_category),pointer :: first_category => null()
      class (type_output_field),pointer    :: first_field    => null()
      class (type_file),pointer      :: next          => null()
   contains
      procedure :: initialize
      procedure :: save
      procedure :: finalize
      procedure :: create_field
   end type type_file

   class (type_host),pointer,save :: host => null()

contains

   subroutine initialize(self,julianday,secondsofday)
      class (type_file),intent(inout) :: self
      integer,          intent(in)    :: julianday,secondsofday
      stop 'output_manager_core:initialize not implemented'
   end subroutine

   function create_field(self) result(field)
      class (type_file),intent(inout) :: self
      class (type_output_field), pointer :: field
      allocate(field)
   end function create_field

   subroutine save(self,julianday,secondsofday)
      class (type_file),intent(inout) :: self
      integer,          intent(in)    :: julianday,secondsofday
      stop 'output_manager_core:save not implemented'
   end subroutine

   subroutine finalize(file)
      class (type_file),intent(inout) :: file
   end subroutine

   subroutine julian_day(self,yyyy,mm,dd,julian)
      class (type_host), intent(in) :: self
      integer, intent(in)  :: yyyy,mm,dd
      integer, intent(out) :: julian
      call output_manager_fatal_error('julian_day','The host of the output manager must implement julian_day subroutine in its type_host-derived class.')
      julian = -1
   end subroutine julian_day

   subroutine calendar_date(self,julian,yyyy,mm,dd)
      class (type_host), intent(in) :: self
      integer, intent(in)  :: julian
      integer, intent(out) :: yyyy,mm,dd
      call output_manager_fatal_error('calendar_date','The host of the output manager must implement calendar_date subroutine in its type_host-derived class.')
      yyyy = -1
      mm = -1
      dd = -1
   end subroutine calendar_date

   subroutine write_time_string(jul,secs,timestr)
      integer,         intent(in)  :: jul,secs
      character(len=*),intent(out) :: timestr

      integer :: ss,min,hh,dd,mm,yy

      hh   = secs/3600
      min  = (secs-hh*3600)/60
      ss   = secs - 3600*hh - 60*min

      call host%calendar_date(jul,yy,mm,dd)

      write(timestr,'(i4.4,a1,i2.2,a1,i2.2,1x,i2.2,a1,i2.2,a1,i2.2)')  &
                           yy,'-',mm,'-',dd,hh,':',min,':',ss
   end subroutine write_time_string

   subroutine output_manager_fatal_error(location,error)
      character(len=*),intent(in) :: location,error
      
      FATAL trim(location)//': '//trim(error)
      stop 'output_manager::output_manager_fatal_error'
   end subroutine

end module output_manager_core
