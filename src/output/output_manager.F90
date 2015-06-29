module output_manager

   use field_manager,type_base_node=>type_node
   use output_manager_core
   use netcdf_output

   use fabm_config_types
   use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length

   implicit none

   public output_manager_init, output_manager_save, output_manager_clean

   private

   class (type_file),pointer :: first_file

contains

   subroutine output_manager_init(field_manager)
      type (type_field_manager), target :: field_manager
      if (.not.associated(host)) call output_manager_fatal_error('output_manager_init','The host of an output manager must set the host pointer before calling output_manager_init')
      nullify(first_file)
      call configure_from_yaml(field_manager)
   end subroutine

   subroutine output_manager_clean()
      class (type_file),pointer :: file
      file => first_file
      do while (associated(file))
         call file%finalize()
         file => file%next
      end do
   end subroutine

   subroutine collect_from_categories(file)
      class (type_file), intent(inout) :: file
      class (type_output_category), pointer :: output_category
      class (type_output_item),pointer  :: output_item
      type (type_category_node) :: list
      class (type_base_node), pointer :: field_node, next_field_node

      output_category => file%first_category
      do while (associated(output_category))
         write (*,*) 'processing output category '//trim(output_category%name)
         list%first_child => null()
         call output_category%source%get_all_fields(list,output_category%output_level)
         field_node => list%first_child
         if (.not.associated(field_node)) call output_manager_fatal_error('collect_from_categories','No variables have been registered under output category "'//trim(output_category%name)//'".')
         do while (associated(field_node))
            select type (field_node)
            class is (type_field_node)
               write (*,*) '  adding '//trim(field_node%field%name)

               ! Create output field, set relevant properties, and prepend to output field list.
               output_item => file%create_field()
               output_item%time_method = output_category%time_method
               select type (output_item)
               class is (type_output_field)
                  output_item%source => field_node%field
                  output_item%output_name = trim(output_category%prefix)//trim(field_node%field%name)//trim(output_category%postfix)
                  output_item%next => file%first_field
                  file%first_field => output_item
               end select
            end select
            next_field_node => field_node%next_sibling
            deallocate(field_node)
            field_node => next_field_node
         end do
         output_category => output_category%next
      end do
   end subroutine collect_from_categories

   subroutine output_manager_save(julianday,secondsofday)
      integer,intent(in) :: julianday,secondsofday

      class (type_file),            pointer :: file
      class (type_output_field),    pointer :: output_field
      integer                               :: yyyy,mm,dd
      logical                               :: output_time_window

      file => first_file
      do while (associated(file))
#if 0
         ! If in output time-window
         output_time_window = (file%first_julian .eq. julianday .and. file%first_seconds .ge. secondsofday) .or. &
                              (file%first_julian .gt. julianday .and. julianday .lt. file%last_julian)      .or. &
                              (file%last_julian  .eq. julianday .and. file%last_seconds  .le. secondsofday)
#endif
         ! Determine whether output is required
         if ((julianday==file%next_julian.and.secondsofday>=file%next_seconds) .or. julianday>file%next_julian) then
            ! Output required
            if (file%next_julian==-1) then
               ! Add variables below selected categories to output
               call collect_from_categories(file)

               ! First check whether all fields included in this file have been registered.
               output_field => file%first_field
               do while (associated(output_field))
                  select case (output_field%source%status)
                     case (status_not_registered)
                        call output_manager_fatal_error('output_manager_save', 'File '//trim(file%path)//': &
                           requested field "'//trim(output_field%source%name)//'" has not been registered with field manager.')
                     case (status_registered_no_data)
                        call output_manager_fatal_error('output_manager_save', 'File '//trim(file%path)//': &
                           data for requested field "'//trim(output_field%source%name)//'" have not been provided.')
                  end select
                  output_field => output_field%next
               end do

               ! Create output file
               call file%initialize(julianday,secondsofday)

               ! Store current time step so next time step can be computed correctly.
               file%next_julian = julianday
               file%next_seconds = secondsofday

               ! Initialize fields based on time integrals
               output_field => file%first_field
               do while (associated(output_field))

                  ! Store instantaneous data by default.
                  output_field%data_0d => output_field%source%data_0d
                  output_field%data_1d => output_field%source%data_1d
                  output_field%data_2d => output_field%source%data_2d
                  output_field%data_3d => output_field%source%data_3d

                  if (output_field%time_method/=time_method_instantaneous) then
                     ! We are not storing the instantaneous value. Create a work array that will be stored instead.
                     if (associated(output_field%source%data_3d)) then
                        allocate(output_field%work_3d(size(output_field%source%data_3d,1),size(output_field%source%data_3d,2),size(output_field%source%data_3d,3)))
                        output_field%data_3d => output_field%work_3d
                     elseif (associated(output_field%source%data_2d)) then
                        allocate(output_field%work_2d(size(output_field%source%data_2d,1),size(output_field%source%data_2d,2)))
                        output_field%data_2d => output_field%work_2d
                     elseif (associated(output_field%source%data_1d)) then
                        allocate(output_field%work_1d(size(output_field%source%data_1d)))
                        output_field%data_1d => output_field%work_1d
                     else
                        output_field%data_0d => output_field%work_0d
                     end if
                  end if

                  select case (output_field%time_method)
                     case (time_method_mean)
                        ! Temporal mean: use initial value on first output.
                        if (associated(output_field%source%data_3d)) then
                           output_field%work_3d(:,:,:) = output_field%source%data_3d
                        elseif (associated(output_field%source%data_2d)) then
                           output_field%work_2d(:,:) = output_field%source%data_2d
                        elseif (associated(output_field%source%data_1d)) then
                           output_field%work_1d(:) = output_field%source%data_1d
                        else
                           output_field%work_0d = output_field%source%data_0d
                        end if
                     case (time_method_integrated)
                        ! Time integral: use zero at first output.
                        if (associated(output_field%source%data_3d)) then
                           output_field%work_3d(:,:,:) = 0.0_rk
                        elseif (associated(output_field%source%data_2d)) then
                           output_field%work_2d(:,:) = 0.0_rk
                        elseif (associated(output_field%source%data_1d)) then
                           output_field%work_1d(:) = 0.0_rk
                        else
                           output_field%work_0d = 0.0_rk
                        end if
                  end select
                  output_field => output_field%next
               end do
            else
               ! Perform temporal averaging where required.
               output_field => file%first_field
               do while (associated(output_field))
                  if (output_field%time_method==time_method_mean) then
                     ! This is a time-integrated field that needs to be incremented.
                     if (allocated(output_field%work_3d)) then
                        output_field%work_3d(:,:,:) = output_field%work_3d/file%n
                     elseif (allocated(output_field%work_2d)) then
                        output_field%work_2d(:,:) = output_field%work_2d/file%n
                     elseif (allocated(output_field%work_1d)) then
                        output_field%work_1d(:) = output_field%work_1d/file%n
                     else
                        output_field%work_0d = output_field%work_0d/file%n
                     end if
                  end if
                  output_field => output_field%next
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
               case (time_unit_hour)
                  file%next_seconds = file%next_seconds + file%time_step*3600
                  file%next_julian = file%next_julian + file%next_seconds/86400
                  file%next_seconds = mod(file%next_seconds,86400)
               case (time_unit_day)
                  file%next_julian = file%next_julian + file%time_step
               case (time_unit_month)
                  call host%calendar_date(julianday,yyyy,mm,dd)
                  mm = mm + file%time_step
                  yyyy = yyyy + (mm-1)/12
                  mm = mod(mm-1,12)+1
                  call host%julian_day(yyyy,mm,dd,file%next_julian)
               case (time_unit_year)
                  call host%calendar_date(julianday,yyyy,mm,dd)
                  yyyy = yyyy + file%time_step
                  call host%julian_day(yyyy,mm,dd,file%next_julian)
            end select

            ! Reset time step counter   
            file%n = 0

            ! Zero out time-step averaged fields (start of new time step)
            output_field => file%first_field
            do while (associated(output_field))
               if (output_field%time_method==time_method_mean) then
                  if (allocated(output_field%work_3d)) then
                     output_field%work_3d(:,:,:) = 0.0_rk
                  elseif (allocated(output_field%work_2d)) then
                     output_field%work_2d(:,:) = 0.0_rk
                  elseif (allocated(output_field%work_1d)) then
                     output_field%work_1d(:) = 0.0_rk
                  else
                     output_field%work_0d = 0.0_rk
                  end if
               end if
               output_field => output_field%next
            end do
         end if

         ! Increment time-integrated fields
         output_field => file%first_field
         do while (associated(output_field))
            select case (output_field%time_method)
               case (time_method_mean,time_method_integrated)
                  ! This is a time-integrated field that needs to be incremented.
                  if (allocated(output_field%work_3d)) then
                     output_field%work_3d(:,:,:) = output_field%work_3d + output_field%source%data_3d
                  elseif (allocated(output_field%work_2d)) then
                     output_field%work_2d(:,:) = output_field%work_2d + output_field%source%data_2d
                  elseif (allocated(output_field%work_1d)) then
                     output_field%work_1d(:) = output_field%work_1d + output_field%source%data_1d
                  else
                     output_field%work_0d = output_field%work_0d + output_field%source%data_0d
                  end if
            end select
            output_field => output_field%next
         end do
         file%n = file%n + 1
         file => file%next
      end do
   end subroutine output_manager_save
   
   subroutine configure_from_yaml(field_manager)
      type (type_field_manager), target :: field_manager
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
                     call process_file(field_manager,trim(pair%key),dict)
                  class default
                     call output_manager_fatal_error('configure_from_yaml','Contents of '//trim(dict%path)//' must be a dictionary, not a single value.')
               end select
               pair => pair%next
            end do
         class default
            call output_manager_fatal_error('configure_from_yaml','input.yaml must contain a dictionary with (variable name : information) pairs.')
      end select
   end subroutine configure_from_yaml

   subroutine process_file(field_manager,path,mapping)
      type (type_field_manager), target :: field_manager
      character(len=*),       intent(in) :: path
      class (type_dictionary),intent(in) :: mapping

      type (type_error),  pointer :: config_error
      class (type_scalar),pointer :: scalar
      integer                     :: time_unit, time_step
      class (type_file),pointer :: file

      ! Determine time unit
      scalar => mapping%get_scalar('time_unit',required=.true.,error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      select case (scalar%string)
         case ('second')
            time_unit = time_unit_second
         case ('hour')
            time_unit = time_unit_hour
         case ('day')
            time_unit = time_unit_day
         case ('month')
            time_unit = time_unit_month
         case ('year')
            time_unit = time_unit_year
         case default
            call output_manager_fatal_error('process_file','Invalid value "'//trim(scalar%string)//'" specified for time_unit of file "'//trim(path)//'". Valid options are second, day, month, year.')
      end select

      ! Determine time step
      time_step = mapping%get_integer('time_step',error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      
      allocate(type_netcdf_file::file)
      file%field_manager => field_manager
      file%path = path
      file%time_unit = time_unit
      file%time_step = time_step
      file%next => first_file
      first_file => file

      call process_group(file,mapping,time_method_instantaneous)
   end subroutine process_file

   recursive subroutine process_group(file,mapping,parent_time_method)
      class (type_file),      intent(inout) :: file
      class (type_dictionary),intent(in)    :: mapping
      integer,                intent(in)    :: parent_time_method

      class (type_list),pointer :: list
      type (type_list_item),pointer :: item
      type (type_error),  pointer :: config_error
      integer :: default_time_method

      default_time_method = mapping%get_integer('time_method',default=parent_time_method,error=config_error)

      ! Get list with groups [if any]
      list => mapping%get_list('groups',required=.false.,error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      if (associated(list)) then
         item => list%first
         do while (associated(item))
            select type (node=>item%node)
               class is (type_dictionary)
                  call process_group(file, node, default_time_method)
               class default
                  call output_manager_fatal_error('process_file','Elements below '//trim(list%path)//' must be dictionaries.')
            end select
            item => item%next
         end do
      end if

      ! Get list with variables
      list => mapping%get_list('variables',required=.true.,error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_file',config_error%message)
      item => list%first
      do while (associated(item))
         select type (node=>item%node)
            class is (type_dictionary)
               call process_variable(file, node)
            class default
               call output_manager_fatal_error('process_file','Elements below '//trim(list%path)//' must be dictionaries.')
         end select
         item => item%next
      end do
   end subroutine process_group
   
   subroutine process_variable(file,mapping)
      class (type_file),      intent(inout) :: file
      class (type_dictionary),intent(in)    :: mapping

      character(len=string_length) :: source_name
      type (type_error),        pointer :: config_error
      class (type_output_item),pointer  :: output_item
      integer                           :: n

      ! Name of source variable
      source_name = mapping%get_string('source',error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

      ! Determine whether to create an output field or an output category
      n = len_trim(source_name)
      if (source_name(n:n)=='*') then
         allocate(type_output_category::output_item)
      else
         output_item => file%create_field()
      end if

      ! Time method
      output_item%time_method = mapping%get_integer('time_method',default=time_method_instantaneous,error=config_error)
      if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

      select type (output_item)
      class is (type_output_field)
         ! Select this variable for output in the field manager.
         output_item%source => file%field_manager%select_for_output(source_name)

         ! Name of output variable (may differ from source name)
         output_item%output_name = mapping%get_string('name',default=source_name,error=config_error)
         if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

         output_item%next => file%first_field
         file%first_field => output_item
      class is (type_output_category)
         if (n==1) then
            output_item%name = ''
         else
            output_item%name = source_name(:n-2)
         end if

         ! Prefix for output name
         output_item%prefix = mapping%get_string('prefix',default='',error=config_error)
         if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

         ! Postfix for output name
         output_item%postfix = mapping%get_string('postfix',default='',error=config_error)
         if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

         ! Output level
         output_item%output_level = mapping%get_integer('output_level',default=output_level_default,error=config_error)
         if (associated(config_error)) call output_manager_fatal_error('process_variable',config_error%message)

         ! Select this category for output in the field manager.
         output_item%source => file%field_manager%select_category_for_output(output_item%name,output_item%output_level)

         output_item%next => file%first_category
         file%first_category => output_item
      end select

   end subroutine process_variable

end module
