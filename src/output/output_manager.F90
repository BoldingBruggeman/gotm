#include"cppdefs.h"

module output_manager

   use time, only: calendar_date,julian_day
   use netcdf
   use fabm_config_types
   use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length

   implicit none

   public output_manager_init, output_manager_save, output_manager_clean, output_manager_register_field
   public id_dim_z,id_dim_z1

   private

   integer,parameter :: string_length = 256
   integer,parameter :: max_path = 256

   integer,parameter :: rk = kind(_ONE_)
   REAL_4B,parameter :: dummy = 0.
   integer,parameter :: ncrk = kind(dummy)

   integer,parameter :: time_method_none          = 0  ! time-independent variable
   integer,parameter :: time_method_instantaneous = 1
   integer,parameter :: time_method_mean          = 2
   integer,parameter :: time_method_integrated    = 3

   integer,parameter :: time_unit_none   = 0
   integer,parameter :: time_unit_second = 1
   integer,parameter :: time_unit_day    = 2
   integer,parameter :: time_unit_month  = 3
   integer,parameter :: time_unit_year   = 4

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

   type type_available_field
      character(len=string_length)  :: name          = ''
      character(len=string_length)  :: units         = ''
      character(len=string_length)  :: long_name     = ''
      character(len=string_length)  :: standard_name = ''
      real(rk)                      :: fill_value    = default_fill_value
      real(rk)                      :: minimum       = default_minimum
      real(rk)                      :: maximum       = default_maximum
      integer,allocatable           :: dimensions(:)
      real(rk),pointer              :: data_0d       => null()
      real(rk),pointer,dimension(:) :: data_1d       => null()
      type (type_available_field),pointer :: next    => null()
   end type

   type type_used_field
      character(len=string_length)        :: output_name    = ''
      character(len=string_length)        :: source_name    = ''
      type (type_available_field),pointer :: source         => null()
      integer                             :: time_method = time_method_instantaneous
      real(rk)                            :: work_0d
      real(rk),allocatable,dimension(:)   :: work_1d
      integer                             :: ncid           = -1
      type (type_used_field),pointer      :: next           => null()
   end type

   type type_file
      character(len=max_path)        :: path          = ''
      integer                        :: time_unit     = time_unit_none
      integer                        :: time_step     = 0
      integer                        :: n             = 0  ! Number of model time steps included so far in next output
      integer                        :: itime         = 0  ! Next time index in NetCDF file
      integer                        :: ncid          = -1 ! NetCDF identifier for file
      integer                        :: time_id       = -1
      integer                        :: first_julian  = 0
      integer                        :: first_seconds = 0
      integer                        :: next_julian   = 0
      integer                        :: next_seconds  = 0
      type (type_used_field),pointer :: first_field   => null()
      type (type_file),pointer       :: next          => null()
   end type

   interface output_manager_register_field
      procedure register_field_0d
      procedure register_field_1d
   end interface

   type (type_available_field),pointer :: first_available_field
   type (type_file),           pointer :: first_file

contains

   subroutine output_manager_init(nlev)
      integer,intent(in) :: nlev
      nx = 1
      ny = 1
      nz = nlev
      nullify(first_available_field)
      nullify(first_file)
      call configure_from_yaml()
   end subroutine

   subroutine output_manager_clean()
      type (type_file),pointer :: file
      integer :: iret
      file => first_file
      do while (associated(file))
         iret = nf90_close(file%ncid); call check_err(iret)
         file => file%next
      end do
   end subroutine
   
   function add_file(path,time_unit,time_step) result(file)
      character(len=*), intent(in) :: path
      integer,          intent(in) :: time_unit,time_step

      type (type_file), pointer :: file

      allocate(file)
      file%path = path
      file%time_unit = time_unit
      file%time_step = time_step
      file%next => first_file
      first_file => file
   end function

   subroutine add_field(field, name, units, long_name, standard_name, fill_value, minimum, maximum, dimensions, used)
      type (type_available_field), target    :: field
      character(len=*), intent(in)           :: name, units, long_name
      character(len=*), intent(in), optional :: standard_name
      integer,optional, intent(in)           :: dimensions(:)
      real(rk),         intent(in), optional :: fill_value, minimum, maximum
      logical,          intent(out),optional :: used

      type (type_file),       pointer :: file
      type (type_used_field), pointer :: used_field

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

      ! Prepend to field list
      field%next => first_available_field
      first_available_field => field

      ! Look over used fields and determine whether they include the current available field.
      if (present(used)) used = .false.
      file => first_file
      do while (associated(file))
         used_field => file%first_field
         do while (associated(used_field))
            if (used_field%source_name==field%name) then
               used_field%source => field
               if (present(used)) used = .true.
            end if
            used_field => used_field%next
         end do
         file => file%next
      end do
   end subroutine add_field
   
   subroutine register_field_0d(name, units, long_name, standard_name, fill_value, minimum, maximum, data, used)
      character(len=*),          intent(in)  :: name, units, long_name
      character(len=*),optional, intent(in)  :: standard_name
      real(rk),        optional, intent(in)  :: fill_value, minimum, maximum
      real(rk),        optional, target      :: data
      logical,         optional, intent(out) :: used

      type (type_available_field), pointer :: field

      allocate(field)
      if (present(data)) field%data_0d => data
      call add_field(field, name, units, long_name, standard_name, fill_value, minimum, maximum, used=used)
   end subroutine register_field_0d
   
   subroutine register_field_1d(name, dimension, units, long_name, standard_name, fill_value, minimum, maximum, data, used)
      character(len=*),          intent(in)  :: name, units, long_name
      integer,                   intent(in)  :: dimension
      character(len=*),optional, intent(in)  :: standard_name
      real(rk),        optional, intent(in)  :: fill_value, minimum, maximum
      real(rk),        optional, target      :: data(:)
      logical,         optional, intent(out) :: used

      type (type_available_field), pointer :: field

      allocate(field)
      if (present(data)) then
         if (size(data)/=nz) call fatal_error('register_field_1d', &
            'length of data array provided for variable '//trim(name)//' does not match extents of the vertical domain.')
         field%data_1d => data
      end if
      call add_field(field, name, units, long_name, standard_name, fill_value, minimum, maximum, (/dimension/), used)
   end subroutine register_field_1d

   subroutine output_manager_save(julianday,secondsofday)
      integer,intent(in) :: julianday,secondsofday

      type (type_file),       pointer :: file
      type (type_used_field), pointer :: used_field
      integer                         :: yyyy,mm,dd

      file => first_file
      do while (associated(file))
         ! Determine whether output is required
         if ((julianday==file%next_julian.and.secondsofday>=file%next_seconds) .or. julianday>file%next_julian) then
            ! Output required
            if (file%ncid==-1) then
               call initialize_netcdf_output(file,julianday,secondsofday)
            else
               ! Perform temporal averaging where required.
               used_field => file%first_field
               do while (associated(used_field))
                  if (used_field%time_method==time_method_mean) then
                     ! This is a time-integrated field that needs to be incremented.
                     if (allocated(used_field%work_1d)) then
                        ! 1D field
                        used_field%work_1d = used_field%work_1d/file%n
                     else
                        ! 0D field
                        used_field%work_0d = used_field%work_0d/file%n
                     end if
                  end if
                  used_field => used_field%next
               end do
            end if

            ! Do NetCDF output
            call do_netcdf_output(file,julianday,secondsofday)

            ! Determine time (juian day, seconds of day) for next output.
            select case (file%time_unit)
               case (time_unit_second)
                  file%next_seconds = file%next_seconds + file%time_step
                  file%next_julian = file%next_julian + file%next_seconds/86400
                  file%next_seconds = mod(file%next_seconds,86400)
               case (time_unit_day)
                  file%next_julian = file%next_julian + file%time_step
               case (time_unit_month)
                  call calendar_date(julianday,yyyy,mm,dd)
                  mm = mm + file%time_step
                  yyyy = yyyy + (mm-1)/12
                  mm = mod(mm-1,12)+1
                  call julian_day(yyyy,mm,dd,file%next_julian)
               case (time_unit_year)
                  call calendar_date(julianday,yyyy,mm,dd)
                  yyyy = yyyy + file%time_step
                  call julian_day(yyyy,mm,dd,file%next_julian)
            end select

            ! Reset time step counter   
            file%n = 0

            ! Zero out time-step averaged fields (start of new time step)
            used_field => file%first_field
            do while (associated(used_field))
               if (used_field%time_method==time_method_mean) then
                  if (allocated(used_field%work_1d)) then
                     used_field%work_1d = 0.0_rk
                  else
                     used_field%work_0d = 0.0_rk
                  end if
               end if
               used_field => used_field%next
            end do
         end if

         ! Increment time-integrated fields
         used_field => file%first_field
         do while (associated(used_field))
            select case (used_field%time_method)
               case (time_method_mean,time_method_integrated)
                  ! This is a time-integrated field that needs to be incremented.
                  if (allocated(used_field%work_1d)) then
                     ! 1D field
                     used_field%work_1d = used_field%work_1d + used_field%source%data_1d
                  else
                     ! 0D field
                     used_field%work_0d = used_field%work_0d + used_field%source%data_0d
                  end if
            end select
            used_field => used_field%next
         end do
         file%n = file%n + 1
         file => file%next
      end do
   end subroutine output_manager_save

   subroutine initialize_netcdf_output(file,julianday,secondsofday)
      type (type_file),intent(inout) :: file
      integer,         intent(in)    :: julianday,secondsofday

      type (type_used_field), pointer :: used_field
      integer                         :: iret
      integer                         :: dims_ids(5), current_dim_ids(5)
      integer                         :: i
      integer                         :: yyyy,mm,dd
      character(len=19)               :: time_string

      ! First check whether all fields included in this file have been registered.
      used_field => file%first_field
      do while (associated(used_field))
         if (.not.associated(used_field%source)) call fatal_error('initialize_netcdf_output', &
            'File '//trim(file%path)//': requested field "'//trim(used_field%source_name)//'" has not been registered with output manager.')
         used_field => used_field%next
      end do
      
      ! Create NetCDF file
      iret = nf90_create(file%path,NF90_CLOBBER,file%ncid); call check_err(iret)

      ! Put in define mode
      !iret = nf90_redef(file%ncid); call check_err(iret)

      ! Create dimensions [TODO: only those used in the current file]
      iret = nf90_def_dim(file%ncid, 'lon',  1,  dims_ids(id_dim_lon )); call check_err(iret)
      iret = nf90_def_dim(file%ncid, 'lat',  1,  dims_ids(id_dim_lat )); call check_err(iret)
      iret = nf90_def_dim(file%ncid, 'z',    nz, dims_ids(id_dim_z   )); call check_err(iret)
      iret = nf90_def_dim(file%ncid, 'z1',   nz, dims_ids(id_dim_z1  )); call check_err(iret)
      iret = nf90_def_dim(file%ncid, 'time', NF90_UNLIMITED, dims_ids(id_dim_time)); call check_err(iret)

      ! Create coordinates
      file%first_julian = julianday
      file%first_seconds = secondsofday
      iret = nf90_def_var(file%ncid,'time',NF90_REAL,(/dims_ids(id_dim_time)/),file%time_id); call check_err(iret)
      call write_time_string(julianday,secondsofday,time_string)
      iret = nf90_put_att(file%ncid,file%time_id,'units','seconds since '//trim(time_string)); call check_err(iret)

      ! create variables
      used_field => file%first_field
      do while (associated(used_field))
         ! Map internal dimension indices to indices in NetCDF file.
         do i=1,size(used_field%source%dimensions)
            current_dim_ids(i) = dims_ids(used_field%source%dimensions(i))
         end do
         iret = nf90_def_var(file%ncid,used_field%output_name, NF90_REAL, current_dim_ids(1:size(used_field%source%dimensions)), used_field%ncid); call check_err(iret)
         iret = nf90_put_att(file%ncid,used_field%ncid,'units',trim(used_field%source%units)); call check_err(iret)
         iret = nf90_put_att(file%ncid,used_field%ncid,'long_name',trim(used_field%source%long_name)); call check_err(iret)
         if (used_field%source%standard_name/='') iret = nf90_put_att(file%ncid,used_field%ncid,'standard_name',trim(used_field%source%standard_name)); call check_err(iret)
         if (used_field%source%minimum/=default_minimum) iret = nf90_put_att(file%ncid,used_field%ncid,'valid_min',real(used_field%source%minimum,ncrk)); call check_err(iret)
         if (used_field%source%maximum/=default_maximum) iret = nf90_put_att(file%ncid,used_field%ncid,'valid_max',real(used_field%source%maximum,ncrk)); call check_err(iret)
         if (used_field%source%fill_value/=default_fill_value) iret = nf90_put_att(file%ncid,used_field%ncid,'_FillValue',real(used_field%source%fill_value,ncrk)); call check_err(iret)
         
         select case (used_field%time_method)
            case (time_method_instantaneous)
               iret = nf90_put_att(file%ncid,used_field%ncid,'cell_methods','time: point'); call check_err(iret)
            case (time_method_mean)
               ! Temporal mean: use initial value on first output.
               iret = nf90_put_att(file%ncid,used_field%ncid,'cell_methods','time: mean'); call check_err(iret)
               if (associated(used_field%source%data_1d)) then
                  allocate(used_field%work_1d(size(used_field%source%data_1d)))
                  used_field%work_1d = used_field%source%data_1d
               else
                  used_field%work_0d = used_field%source%data_0d
               end if
            case (time_method_integrated)
               ! Time integral: use zero at first output.
               iret = nf90_put_att(file%ncid,used_field%ncid,'cell_methods','time: sum'); call check_err(iret)
               if (associated(used_field%source%data_1d)) then
                  allocate(used_field%work_1d(size(used_field%source%data_1d)))
                  used_field%work_1d = 0.0_rk
               else
                  used_field%work_0d = 0.0_rk
               end if
         end select
         used_field => used_field%next
      end do

      ! Exit define mode
      iret = nf90_enddef(file%ncid); call check_err(iret)

      file%next_julian = julianday
      file%next_seconds = secondsofday
   end subroutine initialize_netcdf_output
   
   subroutine do_netcdf_output(file,julianday,secondsofday)
      type (type_file),intent(inout) :: file
      integer,         intent(in)    :: julianday,secondsofday

      type (type_used_field), pointer :: used_field
      integer                         :: iret
      real(ncrk)                      :: temp_time
      integer,dimension(4)            :: start,edges
      integer                         :: i

      ! Increment time index
      file%itime = file%itime + 1

      ! Store time coordinate
      temp_time = (julianday-file%first_julian)*real(86400,ncrk) + secondsofday-file%first_seconds
      start(1) = file%itime
      iret = nf90_put_var(file%ncid,file%time_id,temp_time,start); call check_err(iret)
      
      used_field => file%first_field
      do while (associated(used_field))
         ! Fill arrays with start index and count per dimension
         do i=1,size(used_field%source%dimensions)
            select case (used_field%source%dimensions(i))
               case (id_dim_lon)
                  start(i) = 1
                  edges(i) = nx
               case (id_dim_lat)
                  start(i) = 1
                  edges(i) = ny
               case (id_dim_z,id_dim_z1)
                  start(i) = 1
                  edges(i) = nz
               case (id_dim_time)
                  start(i) = file%itime
                  edges(i) = 1
            end select
         end do

         select case (used_field%time_method)
            case (time_method_instantaneous)
               ! Use source field directly
               if (associated(used_field%source%data_1d)) then
                  iret = nf90_put_var(file%ncid,used_field%ncid,used_field%source%data_1d,start,edges)
               else
                  iret = nf90_put_var(file%ncid,used_field%ncid,used_field%source%data_0d,start)
               end if
            case default
               ! Use work field
               if (allocated(used_field%work_1d)) then
                  iret = nf90_put_var(file%ncid,used_field%ncid,used_field%work_1d,start,edges)
               else
                  iret = nf90_put_var(file%ncid,used_field%ncid,used_field%work_0d,start)
               end if
         end select
         call check_err(iret)
         used_field => used_field%next
      end do
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

   subroutine check_err(iret)
      integer,intent(in) :: iret
      if (iret/=NF90_NOERR) call fatal_error('check_err',nf90_strerror(iret))
   end subroutine
   
   subroutine configure_from_yaml()
      character(len=yaml_error_length)   :: yaml_error
      class (type_node),         pointer :: node
      type (type_key_value_pair),pointer :: pair
      character(len=max_path)            :: file_path
      integer,parameter                  :: yaml_unit = 100

      ! Parse YAML.
      node => yaml_parse('output.yaml',yaml_unit,yaml_error)
      if (yaml_error/='') call fatal_error('configure_from_yaml',trim(yaml_error))

      ! Process root-level dictionary.
      select type (node)
         class is (type_dictionary)
            pair => node%first
            do while (associated(pair))
               if (pair%key=='') call fatal_error('configure_from_yaml','Empty file path specified.')
               select type (dict=>pair%value)
                  class is (type_dictionary)
                     call process_file(trim(pair%key),dict)
                  class default
                     call fatal_error('configure_from_yaml','Contents of '//trim(dict%path)//' must be a dictionary, not a single value.')
               end select
               pair => pair%next
            end do
         class default
            call fatal_error('configure_from_yaml','input.yaml must contain a dictionary with (variable name : information) pairs.')
      end select
   end subroutine configure_from_yaml

   subroutine process_file(path,mapping)
      character(len=*),       intent(in) :: path
      class (type_dictionary),intent(in) :: mapping

      type (type_error),  pointer :: config_error
      class (type_scalar),pointer :: scalar
      integer                     :: time_unit, time_step
      type (type_file),pointer    :: file
      class (type_dictionary),pointer :: variables
      type (type_key_value_pair),pointer :: pair

      ! Determine time unit
      scalar => mapping%get_scalar('time_unit',required=.true.,error=config_error)
      if (associated(config_error)) call fatal_error('process_file',config_error%message)
      select case (scalar%string)
         case ('day')
            time_unit = time_unit_day
         case ('month')
            time_unit = time_unit_month
         case ('year')
            time_unit = time_unit_year
         case default
            call fatal_error('process_file','Invalid value "'//trim(scalar%string)//'" specified for time_unit of file "'//trim(path)//'". Valid options are day, month, year.')
         end select

      ! Determine time step
      time_step = mapping%get_integer('time_step',error=config_error)
      if (associated(config_error)) call fatal_error('process_file',config_error%message)
      
      file => add_file(path,time_unit,time_step)

      ! Get dictionary with variables
      variables => mapping%get_dictionary('variables',required=.true.,error=config_error)
      if (associated(config_error)) call fatal_error('process_file',config_error%message)
      pair => variables%first
      do while (associated(pair))
         if (pair%key=='') call fatal_error('process_file','File "'//trim(path)//'": empty variable name specified.')
         select type (dict=>pair%value)
            class is (type_dictionary)
               call process_variable(file, trim(pair%key),dict)
            class is (type_null)
               call process_variable(file, trim(pair%key))
            class default
               call fatal_error('process_file','Contents of '//trim(dict%path)//' must be a dictionary, not a single value.')
         end select
         pair => pair%next
      end do
   end subroutine process_file

   subroutine process_variable(file,name,mapping)
      type (type_file),       intent(inout)       :: file
      character(len=*),       intent(in)          :: name
      class (type_dictionary),intent(in),optional :: mapping

      type (type_error),     pointer :: config_error
      type (type_used_field),pointer :: field

      allocate(field)
      field%output_name = name
      field%source_name = name
      if (present(mapping)) then
         ! Interpret variable-specific configuration information

         ! Name of source variable (may differ from NetCDF variable name)
         field%source_name = mapping%get_string('source',default=field%source_name,error=config_error)
         if (associated(config_error)) call fatal_error('process_variable',config_error%message)

         ! Time method
         field%time_method = mapping%get_integer('time_method',default=time_method_instantaneous,error=config_error)
         if (associated(config_error)) call fatal_error('process_variable',config_error%message)
      end if

      field%next => file%first_field
      file%first_field => field
   end subroutine process_variable
   
   subroutine fatal_error(location,error)
      character(len=*),intent(in) :: location,error
      
      FATAL trim(location)//': '//trim(error)
      stop 'output_manager::fatal_error'
   end subroutine
end module