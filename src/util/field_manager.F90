#include"cppdefs.h"

module field_manager

   implicit none

   ! Public subroutine and functions
   public field_manager_init, field_manager_list, field_manager_clean, field_manager_register, field_manager_send_data, field_manager_select_for_output

   ! Public data types and variables
   public type_field,first_field,nx,ny,nz

   ! Public parameters
   public string_length,default_fill_value,default_minimum,default_maximum
   public id_dim_lon,id_dim_lat,id_dim_z,id_dim_z1,id_dim_time

   private

   integer,parameter :: string_length = 256
   integer,parameter :: nmaxdims = 10
   integer,parameter :: rk = kind(_ONE_)

   integer, parameter :: id_dim_lon  = 1
   integer, parameter :: id_dim_lat  = 2
   integer, parameter :: id_dim_z    = 3
   integer, parameter :: id_dim_z1   = 4
   integer, parameter :: id_dim_time = 5

   integer            :: nx
   integer            :: ny
   integer            :: nz

   logical :: singleton_dims(nmaxdims)

   real(rk),parameter :: default_fill_value = -huge(_ONE_)
   real(rk),parameter :: default_minimum = default_fill_value + spacing(default_fill_value)
   real(rk),parameter :: default_maximum = huge(_ONE_)

   integer            :: counter=0

   type type_field
      integer                      :: id             = 0
      character(len=string_length) :: name           = ''
      character(len=string_length) :: units          = ''
      character(len=string_length) :: long_name      = ''
      character(len=string_length) :: standard_name  = ''
      real(rk)                     :: fill_value     = default_fill_value
      real(rk)                     :: minimum        = default_minimum
      real(rk)                     :: maximum        = default_maximum
      logical                      :: in_output      = .false.
      integer                      :: status         = 0 ! 1 = metadata provided, 2 = data provided
      integer,allocatable          :: dimensions(:)
      real(rk),pointer             :: data_0d        => null()
      real(rk),pointer             :: data_1d(:)     => null()
      real(rk),pointer             :: data_2d(:,:)   => null()
      real(rk),pointer             :: data_3d(:,:,:) => null()
      type (type_field),pointer    :: next           => null()
   end type type_field

   type (type_field),pointer :: first_field

   interface field_manager_send_data
      module procedure send_data_0d
      module procedure send_data_1d
      module procedure send_data_2d
      module procedure send_data_3d
      module procedure send_data_by_name_0d
      module procedure send_data_by_name_1d
   end interface

contains

   subroutine field_manager_init(nlev)
      integer,intent(in) :: nlev
      counter = 0
      singleton_dims = .true.
!KB      singleton_dims(id_dim_lon) = .false.
!KB      singleton_dims(id_dim_lat) = .false.
      singleton_dims(id_dim_z)   = .false.
      singleton_dims(id_dim_z1)  = .false.
      nx = 1
      ny = 1
      nz = nlev
      nullify(first_field)
   end subroutine

   subroutine field_manager_list()
      type (type_field), pointer :: field, next_field
      character(256) :: line
      field => first_field
      write(line,'(A8,4x,A12,4x,A40)') 'name','unit',adjustl('long_name')
      write(*,*) trim(line)
      write(line,'(A68)') '----------------------------------------------------------------'
      write(*,*) trim(line)
      do while (associated(field))
         write(line,'(I2,2x,A15,2x,A15,2x,A45)') field%id,adjustl(field%name),adjustl(field%units),adjustl(field%long_name)
         write(*,*) trim(line)
!KB         write(*,*) field%dimensions
         next_field => field%next
         field => next_field
      end do
!      stop 'kurt'
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
      counter = counter + 1
      field%id   = counter
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
!KBwrite(*,*) 'aa ',dimensions,size(field%dimensions)
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

   subroutine field_manager_register(name, units, long_name, standard_name, fill_value, minimum, maximum, dimensions, data0d, data1d, data2d, data3d, used)
      character(len=*),          intent(in)  :: name, units, long_name
      character(len=*),optional, intent(in)  :: standard_name
      real(rk),        optional, intent(in)  :: fill_value, minimum, maximum
      integer,         optional, intent(in)  :: dimensions(:)
      real(rk),        optional, target      :: data0d,data1d(:),data2d(:,:),data3d(:,:,:)
      logical,         optional, intent(out) :: used

      type (type_field), pointer :: field

      field => add_field(name, units, long_name, standard_name, fill_value, minimum, maximum, dimensions, used=used)
      if (present(data0d)) call send_data_0d(field,data0d)
      if (present(data1d)) call send_data_1d(field,data1d)
      if (present(data2d)) call send_data_2d(field,data2d)
      if (present(data3d)) call send_data_3d(field,data3d)
   end subroutine field_manager_register

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

   subroutine get_expected_extents(field,extents)
      type (type_field),  intent(in)  :: field
      integer,allocatable,intent(out) :: extents(:)

      integer :: i,n

      ! First count data dimensions (excludes time)
      n = 0
      do i=1,size(field%dimensions)
         if (.not.singleton_dims(field%dimensions(i))) n = n + 1
      end do

      ! Allocate array to hold expected extents
      allocate(extents(n))

      ! Fill array with extents
      n = 0
      do i=1,size(field%dimensions)
         if (.not.singleton_dims(field%dimensions(i))) then
            n = n + 1
            select case (field%dimensions(i))
               case (id_dim_lon); extents(n) = nx
               case (id_dim_lat); extents(n) = ny
               case (id_dim_z,id_dim_z1); extents(n) = nz
               case default
                  call fatal_error('get_expected_extents','Unknown dimension id.')
            end select
         end if
      end do

!KBwrite(*,*) 'bb ',extents

   end subroutine get_expected_extents

   subroutine prepare_send_data(field,extents)
      type (type_field), intent(inout) :: field
      integer,           intent(in)    :: extents(:)

      integer,allocatable              :: reqextents(:)
      integer                          :: i
      character(len=2)                 :: str1,str2,str3

      ! Get expected array extents
      call get_expected_extents(field,reqextents)

      ! Check array rank
      if (size(extents)/=size(reqextents)) then
         write (str1,'(i0)') size(extents)
         write (str2,'(i0)') size(reqextents)
         call fatal_error('prepare_send_data',trim(str1)//'D data provided for '//trim(field%name)//', but this field should have '//trim(str2)//' dimensions.')
      end if

      ! Check array extents
      do i=1,size(reqextents)
         if (extents(i)/=reqextents(i)) then
            write (str1,'(i0)') i
            write (str2,'(i0)') extents(i)
            write (str3,'(i0)') reqextents(i)
            call fatal_error('prepare_send_data', 'Field '//trim(field%name)//', dimension  '//trim(str1)//': &
               &extents of provided data ('//trim(str2)//') does not match expected value '//trim(str3)//'.')
         end if
      end do

      if (field%status>1) call fatal_error('prepare_send_data','Data for field "'//trim(field%name)//'" have already been provided.')
      field%status = 2
   end subroutine prepare_send_data

   subroutine send_data_0d(field, data)
      type (type_field), intent(inout) :: field
      real(rk),target                  :: data
      call prepare_send_data(field,shape(data))
      field%data_0d => data
   end subroutine send_data_0d

   subroutine send_data_1d(field, data)
      type (type_field), intent(inout) :: field
      real(rk),target                  :: data(:)
      call prepare_send_data(field,shape(data))
      field%data_1d => data
   end subroutine send_data_1d

   subroutine send_data_2d(field, data)
      type (type_field), intent(inout) :: field
      real(rk),target                  :: data(:,:)
      call prepare_send_data(field,shape(data))
      field%data_2d => data
   end subroutine send_data_2d

   subroutine send_data_3d(field, data)
      type (type_field), intent(inout) :: field
      real(rk),target                  :: data(:,:,:)
      call prepare_send_data(field,shape(data))
      field%data_3d => data
   end subroutine send_data_3d

   subroutine fatal_error(location,error)
      character(len=*),intent(in) :: location,error

      FATAL trim(location)//': '//trim(error)
      stop 'field_manager::fatal_error'
   end subroutine

end module field_manager
