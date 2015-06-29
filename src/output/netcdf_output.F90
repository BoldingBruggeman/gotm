module netcdf_output

   use field_manager
   use output_manager_core
   use netcdf

   implicit none

   public type_netcdf_file

   private

   REAL_4B,parameter :: dummy = 0.
   integer,parameter :: ncrk = kind(dummy)

   type,extends(type_file) :: type_netcdf_file
      integer :: itime         = 0  ! Next time index in NetCDF file
      integer :: ncid          = -1 ! NetCDF identifier for file
      integer :: time_id       = -1 ! Identifier of time dimension
      integer :: first_julian  = 0
      integer :: first_seconds = 0
      integer :: sync_interval = 1  ! Number of output time step between calls to nf90_sync (-1 to disbale syncing)
   contains
      procedure :: initialize
      procedure :: save
      procedure :: finalize
      procedure :: create_field
   end type

   type,extends(type_output_field) :: type_netcdf_field
      integer :: varid = -1
      integer,allocatable :: start(:)
      integer,allocatable :: edges(:)
      integer :: itimedim = -1
   end type

contains

   subroutine initialize(self,julianday,secondsofday)
      class (type_netcdf_file),intent(inout) :: self
      integer,                 intent(in)    :: julianday,secondsofday

      class (type_output_field), pointer :: output_field
      integer                            :: iret
      integer                            :: dims_ids(5), current_dim_ids(5)
      integer                            :: i
      integer                            :: yyyy,mm,dd
      character(len=19)                  :: time_string
      
      ! Create NetCDF file
      iret = nf90_create(trim(self%path)//'.nc',NF90_CLOBBER,self%ncid); call check_err(iret)

      ! Create dimensions [TODO: only those used in the current file]
      iret = nf90_def_dim(self%ncid, 'lon',  self%field_manager%dimension_length(id_dim_lon), dims_ids(id_dim_lon )); call check_err(iret)
      iret = nf90_def_dim(self%ncid, 'lat',  self%field_manager%dimension_length(id_dim_lat), dims_ids(id_dim_lat )); call check_err(iret)
      iret = nf90_def_dim(self%ncid, 'z',    self%field_manager%dimension_length(id_dim_z),   dims_ids(id_dim_z   )); call check_err(iret)
      iret = nf90_def_dim(self%ncid, 'z1',   self%field_manager%dimension_length(id_dim_z1),  dims_ids(id_dim_z1  )); call check_err(iret)
      iret = nf90_def_dim(self%ncid, 'time', NF90_UNLIMITED, dims_ids(id_dim_time)); call check_err(iret)

      ! Create coordinates
      self%first_julian = julianday
      self%first_seconds = secondsofday
      iret = nf90_def_var(self%ncid,'time',NF90_REAL,(/dims_ids(id_dim_time)/),self%time_id); call check_err(iret)
      call write_time_string(julianday,secondsofday,time_string)
      iret = nf90_put_att(self%ncid,self%time_id,'units','seconds since '//trim(time_string)); call check_err(iret)

      ! Create variables
      output_field => self%first_field
      do while (associated(output_field))
         select type (output_field)
         class is (type_netcdf_field)
            ! Map internal dimension indices to indices in NetCDF file.
            do i=1,size(output_field%source%dimensions)
               current_dim_ids(i) = dims_ids(output_field%source%dimensions(i))
            end do
            iret = nf90_def_var(self%ncid,output_field%output_name, NF90_REAL, current_dim_ids(1:size(output_field%source%dimensions)), output_field%varid); call check_err(iret)
            iret = nf90_put_att(self%ncid,output_field%varid,'units',trim(output_field%source%units)); call check_err(iret)
            iret = nf90_put_att(self%ncid,output_field%varid,'long_name',trim(output_field%source%long_name)); call check_err(iret)
            if (output_field%source%standard_name/='') iret = nf90_put_att(self%ncid,output_field%varid,'standard_name',trim(output_field%source%standard_name)); call check_err(iret)
            if (output_field%source%minimum/=default_minimum) iret = nf90_put_att(self%ncid,output_field%varid,'valid_min',real(output_field%source%minimum,ncrk)); call check_err(iret)
            if (output_field%source%maximum/=default_maximum) iret = nf90_put_att(self%ncid,output_field%varid,'valid_max',real(output_field%source%maximum,ncrk)); call check_err(iret)
            if (output_field%source%fill_value/=default_fill_value) iret = nf90_put_att(self%ncid,output_field%varid,'_FillValue',real(output_field%source%fill_value,ncrk)); call check_err(iret)
         
            select case (output_field%time_method)
               case (time_method_instantaneous)
                  iret = nf90_put_att(self%ncid,output_field%varid,'cell_methods','time: point'); call check_err(iret)
               case (time_method_mean)
                  iret = nf90_put_att(self%ncid,output_field%varid,'cell_methods','time: mean'); call check_err(iret)
               case (time_method_integrated)
                  iret = nf90_put_att(self%ncid,output_field%varid,'cell_methods','time: sum'); call check_err(iret)
            end select
   
            ! Fill arrays with start index and count per dimension
            allocate(output_field%start(size(output_field%source%dimensions)))
            allocate(output_field%edges(size(output_field%source%dimensions)))
            do i=1,size(output_field%source%dimensions)
               select case (output_field%source%dimensions(i))
                  case (id_dim_time)
                     output_field%start(i) = self%itime
                     output_field%edges(i) = 1
                     output_field%itimedim = i
                  case default
                     output_field%start(i) = 1
                     output_field%edges(i) = self%field_manager%dimension_length(output_field%source%dimensions(i))
               end select
            end do
         end select   
         output_field => output_field%next
      end do

      ! Exit define mode
      iret = nf90_enddef(self%ncid); call check_err(iret)

   end subroutine initialize

   function create_field(self) result(field)
      class (type_netcdf_file),intent(inout) :: self
      class (type_output_field), pointer :: field
      allocate(type_netcdf_field::field)
   end function create_field

   subroutine save(self,julianday,secondsofday)
      class (type_netcdf_file),intent(inout) :: self
      integer,                 intent(in)    :: julianday,secondsofday

      class (type_output_field), pointer :: output_field
      integer                            :: iret
      real(ncrk)                         :: temp_time

      ! Increment time index
      self%itime = self%itime + 1

      ! Store time coordinate
      temp_time = (julianday-self%first_julian)*real(86400,ncrk) + secondsofday-self%first_seconds
      iret = nf90_put_var(self%ncid,self%time_id,temp_time,(/self%itime/)); call check_err(iret)
      
      output_field => self%first_field
      do while (associated(output_field))
         select type (output_field)
         class is (type_netcdf_field)
            if (output_field%itimedim/=-1) output_field%start(output_field%itimedim) = self%itime
            if (associated(output_field%source%data_3d)) then
               iret = nf90_put_var(self%ncid,output_field%varid,output_field%data_3d,output_field%start,output_field%edges)
            elseif (associated(output_field%source%data_2d)) then
               iret = nf90_put_var(self%ncid,output_field%varid,output_field%data_2d,output_field%start,output_field%edges)
            elseif (associated(output_field%source%data_1d)) then
               iret = nf90_put_var(self%ncid,output_field%varid,output_field%data_1d,output_field%start,output_field%edges)
            else
               iret = nf90_put_var(self%ncid,output_field%varid,output_field%data_0d,output_field%start)
            end if
            call check_err(iret)
         end select
         output_field => output_field%next
      end do

      if (self%sync_interval>0 .and. mod(self%itime,self%sync_interval)==0) then
         iret = nf90_sync(self%ncid)
         call check_err(iret)
      end if
   end subroutine save

   subroutine finalize(file)
      class (type_netcdf_file),intent(inout) :: file
      integer :: iret
      iret = nf90_close(file%ncid); call check_err(iret)
   end subroutine finalize

   subroutine check_err(iret)
      integer,intent(in) :: iret
      if (iret/=NF90_NOERR) call output_manager_fatal_error('check_err',nf90_strerror(iret))
   end subroutine

end module netcdf_output