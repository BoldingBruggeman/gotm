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
   end type

contains

   subroutine initialize(self,julianday,secondsofday)
      class (type_netcdf_file),intent(inout) :: self
      integer,                 intent(in)    :: julianday,secondsofday

      type (type_used_field), pointer :: used_field
      integer                         :: iret
      integer                         :: dims_ids(5), current_dim_ids(5)
      integer                         :: i
      integer                         :: yyyy,mm,dd
      character(len=19)               :: time_string
      
      ! Create NetCDF file
      iret = nf90_create(self%path,NF90_CLOBBER,self%ncid); call check_err(iret)

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
      used_field => self%first_field
      do while (associated(used_field))
         ! Map internal dimension indices to indices in NetCDF file.
         do i=1,size(used_field%source%dimensions)
            current_dim_ids(i) = dims_ids(used_field%source%dimensions(i))
         end do
         iret = nf90_def_var(self%ncid,used_field%output_name, NF90_REAL, current_dim_ids(1:size(used_field%source%dimensions)), used_field%ncid); call check_err(iret)
         iret = nf90_put_att(self%ncid,used_field%ncid,'units',trim(used_field%source%units)); call check_err(iret)
         iret = nf90_put_att(self%ncid,used_field%ncid,'long_name',trim(used_field%source%long_name)); call check_err(iret)
         if (used_field%source%standard_name/='') iret = nf90_put_att(self%ncid,used_field%ncid,'standard_name',trim(used_field%source%standard_name)); call check_err(iret)
         if (used_field%source%minimum/=default_minimum) iret = nf90_put_att(self%ncid,used_field%ncid,'valid_min',real(used_field%source%minimum,ncrk)); call check_err(iret)
         if (used_field%source%maximum/=default_maximum) iret = nf90_put_att(self%ncid,used_field%ncid,'valid_max',real(used_field%source%maximum,ncrk)); call check_err(iret)
         if (used_field%source%fill_value/=default_fill_value) iret = nf90_put_att(self%ncid,used_field%ncid,'_FillValue',real(used_field%source%fill_value,ncrk)); call check_err(iret)
         
         select case (used_field%time_method)
            case (time_method_instantaneous)
               iret = nf90_put_att(self%ncid,used_field%ncid,'cell_methods','time: point'); call check_err(iret)
            case (time_method_mean)
               iret = nf90_put_att(self%ncid,used_field%ncid,'cell_methods','time: mean'); call check_err(iret)
            case (time_method_integrated)
               iret = nf90_put_att(self%ncid,used_field%ncid,'cell_methods','time: sum'); call check_err(iret)
         end select
         used_field => used_field%next
      end do

      ! Exit define mode
      iret = nf90_enddef(self%ncid); call check_err(iret)

   end subroutine initialize
   
   subroutine save(self,julianday,secondsofday)
      class (type_netcdf_file),intent(inout) :: self
      integer,                 intent(in)    :: julianday,secondsofday

      type (type_used_field), pointer :: used_field
      integer                         :: iret
      real(ncrk)                      :: temp_time
      integer,dimension(4)            :: start,edges
      integer                         :: i

      ! Increment time index
      self%itime = self%itime + 1

      ! Store time coordinate
      temp_time = (julianday-self%first_julian)*real(86400,ncrk) + secondsofday-self%first_seconds
      start(1) = self%itime
      iret = nf90_put_var(self%ncid,self%time_id,temp_time,start); call check_err(iret)
      
      used_field => self%first_field
      do while (associated(used_field))
         ! Fill arrays with start index and count per dimension
         do i=1,size(used_field%source%dimensions)
            select case (used_field%source%dimensions(i))
               case (id_dim_time)
                  start(i) = self%itime
                  edges(i) = 1
               case default
                  start(i) = 1
                  edges(i) = self%field_manager%dimension_length(used_field%source%dimensions(i))
            end select
         end do

         select case (used_field%time_method)
            case (time_method_instantaneous)
               ! Use source field directly
               if (associated(used_field%source%data_1d)) then
                  iret = nf90_put_var(self%ncid,used_field%ncid,used_field%source%data_1d,start,edges)
               else
                  iret = nf90_put_var(self%ncid,used_field%ncid,used_field%source%data_0d,start)
               end if
            case default
               ! Use work field
               if (allocated(used_field%work_1d)) then
                  iret = nf90_put_var(self%ncid,used_field%ncid,used_field%work_1d,start,edges)
               else
                  iret = nf90_put_var(self%ncid,used_field%ncid,used_field%work_0d,start)
               end if
         end select
         call check_err(iret)
         used_field => used_field%next
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