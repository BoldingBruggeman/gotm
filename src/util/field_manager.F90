#include"cppdefs.h"

module field_manager

   implicit none

   ! Public subrotutine and functions
   public field_manager_init, field_manager_clean, field_manager_register, field_manager_send_data, field_manager_select_for_output

   ! Public data types and variables
   public type_field,first_field,nx,ny,nz

   ! Public parameters
   public string_length,default_fill_value,default_minimum,default_maximum
   public id_dim_lon,id_dim_lat,id_dim_z,id_dim_z1,id_dim_time

   private

   integer,parameter :: string_length = 256
   integer,parameter :: rk = kind(_ONE_)

   integer, parameter :: id_dim_lon  = 1
   integer, parameter :: id_dim_lat  = 2
   integer, parameter :: id_dim_z    = 3
   integer, parameter :: id_dim_z1   = 4
   integer, parameter :: id_dim_time = 5

   integer            :: nx
   integer            :: ny
   integer            :: nz

   real(rk),parameter :: default_fill_value = -huge(_ONE_)
   real(rk),parameter :: default_minimum = default_fill_value + spacing(default_fill_value)
   real(rk),parameter :: default_maximum = huge(_ONE_)

   type type_field
      character(len=string_length)  :: name          = ''
      character(len=string_length)  :: units         = ''
      character(len=string_length)  :: long_name     = ''
      character(len=string_length)  :: standard_name = ''
      real(rk)                      :: fill_value    = default_fill_value
      real(rk)                      :: minimum       = default_minimum
      real(rk)                      :: maximum       = default_maximum
      logical                       :: in_output     = .false.
      integer                       :: status        = 0 ! 1 = metadata provided, 2 = data provided
      integer,allocatable           :: dimensions(:)
      real(rk),pointer              :: data_0d       => null()
      real(rk),pointer,dimension(:) :: data_1d       => null()
      type (type_field),pointer :: next    => null()
   end type type_field

   type (type_field),pointer :: first_field

   interface field_manager_register
      module procedure register_field_0d
      module procedure register_field_1d
   end interface

   interface field_manager_send_data
      procedure send_data_0d
      procedure send_data_1d
      procedure send_data_by_name_0d
      procedure send_data_by_name_1d
   end interface

contains

   subroutine field_manager_init(nlev)
      integer,intent(in) :: nlev
      nx = 1
      ny = 1
      nz = nlev
      nullify(first_field)
   end subroutine

   subroutine field_manager_clean()
      type (type_field), pointer :: field, next_field
      field => first_field
      do while (associated(field))
         next_field => field%next
         deallocate(field)
         field => next_field
      end do
   end subroutine

   function add_field(name, units, long_name, standard_name, fill_value, minimum, maximum, dimensions, used) result(field)
      character(len=*), intent(in)           :: name, units, long_name
      character(len=*), intent(in), optional :: standard_name
      integer,optional, intent(in)           :: dimensions(:)
      real(rk),         intent(in), optional :: fill_value, minimum, maximum
      logical,          intent(out),optional :: used
      type (type_field), pointer :: field

      ! Find existing field (possible created by select_for_output) or create new one.
      field => find_field(name,create=.true.)
      if (field%status>0) call fatal_error('add_field','Field with name "'//trim(name)//'" has already been registered.')
      field%status = 1

      ! Copy field configuration
      field%name = name
      field%units = units
      field%long_name = long_name
      if (present(standard_name)) field%standard_name = standard_name
      if (present(fill_value)) field%fill_value = fill_value
      if (present(minimum)) field%minimum = minimum
      if (present(maximum)) field%maximum = maximum
      if (present(dimensions)) then
         allocate(field%dimensions(size(dimensions)+3))
         field%dimensions(3:2+size(dimensions)) = dimensions
      else
         allocate(field%dimensions(3))
      end if

      ! Add implicit dimensions
      field%dimensions(1) = id_dim_lon
      field%dimensions(2) = id_dim_lat
      field%dimensions(size(field%dimensions)) = id_dim_time

      ! Look over used fields and determine whether they include the current available field.
      if (present(used)) used = field%in_output
   end function add_field

   function field_manager_select_for_output(name) result(field)
      character(len=*), intent(in) :: name
      type (type_field), pointer :: field

      field => find_field(name,create=.true.)
      field%in_output = .true.
   end function field_manager_select_for_output

   function find_field(name,create) result(field)
      character(len=*),intent(in) :: name
      logical,optional,intent(in) :: create
      type (type_field), pointer :: field

      logical :: create_eff

      field => first_field
      do while (associated(field))
         if (field%name==name) return
         field => field%next
      end do

      create_eff = .false.
      if (present(create)) create_eff = create
      if (create) then
         allocate(field)
         field%name = name
         field%next => first_field
         first_field => field
      end if
   end function find_field

   subroutine register_field_0d(name, units, long_name, standard_name, fill_value, minimum, maximum, data, used)
      character(len=*),          intent(in)  :: name, units, long_name
      character(len=*),optional, intent(in)  :: standard_name
      real(rk),        optional, intent(in)  :: fill_value, minimum, maximum
      real(rk),        optional, target      :: data
      logical,         optional, intent(out) :: used

      type (type_field), pointer :: field

      field => add_field(name, units, long_name, standard_name, fill_value, minimum, maximum, used=used)
      if (present(data)) call send_data_0d(field,data)
   end subroutine register_field_0d
   
   subroutine register_field_1d(name, dimension, units, long_name, standard_name, fill_value, minimum, maximum, data, used)
      character(len=*),          intent(in)  :: name, units, long_name
      integer,                   intent(in)  :: dimension
      character(len=*),optional, intent(in)  :: standard_name
      real(rk),        optional, intent(in)  :: fill_value, minimum, maximum
      real(rk),        optional, target      :: data(:)
      logical,         optional, intent(out) :: used

      type (type_field), pointer :: field

      field => add_field(name, units, long_name, standard_name, fill_value, minimum, maximum, (/dimension/), used)
      if (present(data)) call send_data_1d(field,data)
   end subroutine register_field_1d

   subroutine send_data_by_name_0d(name, data)
      character(len=*),intent(in) :: name
      real(rk),        target     :: data

      type (type_field), pointer :: field

      field => find_field(name)
      if (.not.associated(field)) call fatal_error('send_data_by_name_0d','Field "'//trim(name)//'" has not been registered.')
      call send_data_0d(field,data)
   end subroutine send_data_by_name_0d

   subroutine send_data_by_name_1d(name, data)
      character(len=*),intent(in) :: name
      real(rk),        target     :: data(:)

      type (type_field), pointer :: field

      field => find_field(name)
      if (.not.associated(field)) call fatal_error('send_data_by_name_1d','Field "'//trim(name)//'" has not been registered.')
      call send_data_1d(field,data)
   end subroutine send_data_by_name_1d

   subroutine send_data_0d(field, data)
      type (type_field), intent(inout) :: field
      real(rk),target                            :: data
      if (field%status>1) call fatal_error('add_field','Data for field "'//trim(field%name)//'" has already been provided.')
      field%status = 2
      field%data_0d => data
   end subroutine send_data_0d

   subroutine send_data_1d(field, data)
      type (type_field), intent(inout) :: field
      real(rk),target                            :: data(:)
      if (field%status>1) call fatal_error('add_field','Data for field "'//trim(field%name)//'" has already been provided.')
      field%status = 2
      if (size(data)/=nz) call fatal_error('send_data_1d', &
         'length of data array provided for variable '//trim(field%name)//' does not match extents of its spatial domain.')
      field%data_1d => data
   end subroutine send_data_1d
   
   subroutine fatal_error(location,error)
      character(len=*),intent(in) :: location,error
      
      FATAL trim(location)//': '//trim(error)
      stop 'field_manager::fatal_error'
   end subroutine

end module field_manager