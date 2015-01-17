module output_manager

   use time, only: calendar_date,julian_day
   use field_manager
   use output_manager_core
   use netcdf_output

   use fabm_config_types
   use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length

   implicit none

   public output_manager_init, output_manager_save, output_manager_clean

   private

   class (type_file),pointer :: first_file

contains

   subroutine output_manager_init()
      nullify(first_file)
      call configure_from_yaml()
   end subroutine

   subroutine output_manager_clean()
      class (type_file),pointer :: file
      file => first_file
      do while (associated(file))
         call file%finalize()
         file => file%next
      end do
   end subroutine

   subroutine output_manager_save(julianday,secondsofday)
      integer,intent(in) :: julianday,secondsofday

      class (type_file),      pointer :: file
      type (type_used_field), pointer :: used_field
      integer                         :: yyyy,mm,dd

      file => first_file
      do while (associated(file))
         ! Determine whether output is required
         if ((julianday==file%next_julian.and.secondsofday>=file%next_seconds) .or. julianday>file%next_julian) then
            ! Output required
            if (file%next_julian==-1) then

               ! First check whether all fields included in this file have been registered.
               used_field => file%first_field
               do while (associated(used_field))
                  select case (used_field%source%status)
                     case (0)
                        call output_manager_fatal_error('output_manager_save', 'File '//trim(file%path)//': &
                           requested field "'//trim(used_field%source%name)//'" has not been registered with field manager.')
                     case (1)
                        call output_manager_fatal_error('output_manager_save', 'File '//trim(file%path)//': &
                           data for requested field "'//trim(used_field%source%name)//'" have not been provided.')
                  end select
                  used_field => used_field%next
               end do

               ! Create output file
               call file%initialize(julianday,secondsofday)

               ! Store current time step so next time step can be computed correctly.
               file%next_julian = julianday
               file%next_seconds = secondsofday

               ! Initialize fields based on time integrals
               used_field => file%first_field
               do while (associated(used_field))
                  select case (used_field%time_method)
                     case (time_method_mean)
                        ! Temporal mean: use initial value on first output.
                        if (associated(used_field%source%data_1d)) then
                           allocate(used_field%work_1d(size(used_field%source%data_1d)))
                           used_field%work_1d = used_field%source%data_1d
                        else
                           used_field%work_0d = used_field%source%data_0d
                        end if
                     case (time_method_integrated)
                        ! Time integral: use zero at first output.
                        if (associated(used_field%source%data_1d)) then
                           allocate(used_field%work_1d(size(used_field%source%data_1d)))
                           used_field%work_1d = 0.0_rk
                        else
                           used_field%work_0d = 0.0_rk
                        end if
                  end select
                  used_field => used_field%next
               end do
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
            call file%save(julianday,secondsofday)

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
   
   subroutine configure_from_yaml()
      character(len=yaml_error_length)   :: yaml_error
      class (type_node),         pointer :: node
      type (type_key_value_pair),pointer :: pair
      character(len=max_path)            :: file_path
      integer,parameter                  :: yaml_unit = 100

      ! Parse YAML.
      node => yaml_parse('output.yaml',yaml_unit,yaml_error)
      if (yaml_error/='') call output_manager_fatal_error('configure_from_yaml',trim(yaml_error))

      ! Process root-level dictionary.
      select type (node)
         class is (type_dictionary)
            pair => node%first
            do while (associated(pair))
               if (pair%key=='') call output_manager_fatal_error('configure_from_yaml','Empty file path specified.')
               select type (dict=>pair%value)
                  class is (type_dictionary)
                     call process_file(trim(pair%key),dict)
                  class default
                     call output_manager_fatal_error('configure_from_yaml','Contents of '//trim(dict%path)//' must be a dictionary, not a single value.')
               end select
               pair => pair%next
            end do
         class default
            call output_manager_fatal_error('configure_from_yaml','input.yaml must contain a dictionary with (variable name : information) pairs.')
      end select
   end subroutine configure_from_yaml

   subroutine process_file(path,mapping)
      character(len=*),       intent(in) :: path
      class (type_dictionary),intent(in) :: mapping

      type (type_error),  pointer :: config_error
      class (type_scalar),pointer :: scalar
      integer                     :: time_unit, time_step
      class (type_file),pointer :: file
      class (type_dictionary),pointer :: variables
      type (type_key_value_pair),pointer :: pair

      ! Determine time unit
      scalar => mapping%get_scalar('time_unit',required=.true.,error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      select case (scalar%string)
         case ('day')
            time_unit = time_unit_day
         case ('month')
            time_unit = time_unit_month
         case ('year')
            time_unit = time_unit_year
         case default
            call output_manager_fatal_error('process_file','Invalid value "'//trim(scalar%string)//'" specified for time_unit of file "'//trim(path)//'". Valid options are day, month, year.')
         end select

      ! Determine time step
      time_step = mapping%get_integer('time_step',error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      
      allocate(type_netcdf_file::file)
      file%path = path
      file%time_unit = time_unit
      file%time_step = time_step
      file%next => first_file
      first_file => file

      ! Get dictionary with variables
      variables => mapping%get_dictionary('variables',required=.true.,error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      pair => variables%first
      do while (associated(pair))
         if (pair%key=='') call output_manager_fatal_error('process_file','File "'//trim(path)//'": empty variable name specified.')
         select type (dict=>pair%value)
            class is (type_dictionary)
               call process_variable(file, trim(pair%key),dict)
            class is (type_null)
               call process_variable(file, trim(pair%key))
            class default
               call output_manager_fatal_error('process_file','Contents of '//trim(dict%path)//' must be a dictionary, not a single value.')
         end select
         pair => pair%next
      end do
   end subroutine process_file

   subroutine process_variable(file,name,mapping)
      class (type_file),      intent(inout)       :: file
      character(len=*),       intent(in)          :: name
      class (type_dictionary),intent(in),optional :: mapping

      character(len=string_length)   :: source_name = ''
      type (type_error),     pointer :: config_error
      type (type_used_field),pointer :: field

      allocate(field)
      field%output_name = name
      if (present(mapping)) then
         ! Interpret variable-specific configuration information

         ! Name of source variable (may differ from NetCDF variable name)
         source_name = mapping%get_string('source',default=name,error=config_error)
         if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

         ! Time method
         field%time_method = mapping%get_integer('time_method',default=time_method_instantaneous,error=config_error)
         if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)
      else
         source_name = name
      end if
      field%source => field_manager_select_for_output(source_name)

      field%next => file%first_field
      file%first_field => field
   end subroutine process_variable

end module
