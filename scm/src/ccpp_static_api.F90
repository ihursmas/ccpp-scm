
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated API for the CCPP static build
!!
!
module ccpp_static_api

   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_tsinit_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_tsfinal_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_init_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_run_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_final_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_tsinit_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_tsfinal_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_init_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_run_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_final_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_radiation_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_radiation_tsinit_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_radiation_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_radiation_tsfinal_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_radiation_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_radiation_init_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_radiation_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_radiation_run_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_radiation_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_radiation_final_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_tsinit_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_tsfinal_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_init_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_run_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_final_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_tsinit_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_tsfinal_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_init_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_run_cap
   use ccpp_SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_cap, only: SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_final_cap
   use ccpp_types, only: one
   use scm_type_defs, only: physics
   use scm_physical_constants, only: con_pi
   use scm_physical_constants, only: con_g
   use scm_physical_constants, only: con_t0c
   use scm_physical_constants, only: con_hfus
   use scm_physical_constants, only: con_solr_2008
   use scm_physical_constants, only: con_solr_2002
   use scm_physical_constants, only: con_c
   use scm_physical_constants, only: con_plnk
   use scm_physical_constants, only: con_boltz
   use scm_physical_constants, only: con_rd
   use GFS_typedefs, only: LTP
   use scm_physical_constants, only: con_rerth
   use scm_physical_constants, only: con_p0
   use scm_physical_constants, only: con_omega
   use scm_physical_constants, only: con_cp
   use scm_physical_constants, only: con_rv
   use scm_physical_constants, only: con_fvirt
   use scm_physical_constants, only: con_ttp
   use scm_physical_constants, only: con_hvap
   use scm_physical_constants, only: con_eps
   use scm_physical_constants, only: con_thgni
   use scm_physical_constants, only: con_epsm1
   use scm_physical_constants, only: con_rog
   use scm_physical_constants, only: con_rocp
   use scm_physical_constants, only: con_tice
   use scm_physical_constants, only: con_sbc
   use scm_physical_constants, only: con_jcal
   use scm_physical_constants, only: con_rhw0
   use scm_physical_constants, only: rlapse
   use scm_physical_constants, only: rhowater
   use scm_physical_constants, only: karman
   use scm_physical_constants, only: con_epsq
   use scm_physical_constants, only: con_1ovg
   use scm_physical_constants, only: con_cliq
   use scm_physical_constants, only: con_cvap
   use scm_physical_constants, only: rainmin

   implicit none

   private
   public :: ccpp_physics_timestep_init,ccpp_physics_timestep_finalize,ccpp_physics_init,ccpp_physics_run,ccpp_physics_finalize

   contains

   subroutine ccpp_physics_timestep_init(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="SCM_GFS_v17_p8_ugwpv1_pumas") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_tsinit_cap(one=one,physics=physics,cdata=cdata,con_pi=con_pi)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_radiation_tsinit_cap()
            else if (trim(group_name)=="phys_ps") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_tsinit_cap(one=one,physics=physics,cdata=cdata,con_g=con_g,con_t0c=con_t0c,con_hfus=con_hfus)
            else if (trim(group_name)=="phys_ts") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_tsinit_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v17_p8_ugwpv1_pumas_tsinit_cap(one=one,physics=physics,cdata=cdata,con_pi=con_pi,con_g=con_g,con_t0c=con_t0c, &
                  con_hfus=con_hfus)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_timestep_init

   subroutine ccpp_physics_timestep_finalize(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="SCM_GFS_v17_p8_ugwpv1_pumas") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_tsfinal_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_radiation_tsfinal_cap()
            else if (trim(group_name)=="phys_ps") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_tsfinal_cap()
            else if (trim(group_name)=="phys_ts") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_tsfinal_cap()
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v17_p8_ugwpv1_pumas_tsfinal_cap()

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_timestep_finalize

   subroutine ccpp_physics_init(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="SCM_GFS_v17_p8_ugwpv1_pumas") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_init_cap(one=one,physics=physics,cdata=cdata,con_solr_2008=con_solr_2008,con_solr_2002=con_solr_2002, &
                  con_pi=con_pi,con_c=con_c,con_plnk=con_plnk,con_boltz=con_boltz,con_t0c=con_t0c, &
                  con_rd=con_rd,con_g=con_g,LTP=LTP)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_radiation_init_cap(one=one,physics=physics,cdata=cdata,con_pi=con_pi)
            else if (trim(group_name)=="phys_ps") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_init_cap(one=one,physics=physics,cdata=cdata,con_pi=con_pi,con_rerth=con_rerth,con_p0=con_p0, &
                  con_g=con_g,con_omega=con_omega,con_cp=con_cp,con_rd=con_rd,con_rv=con_rv, &
                  con_fvirt=con_fvirt)
            else if (trim(group_name)=="phys_ts") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_init_cap(one=one,physics=physics,cdata=cdata,con_g=con_g,con_rd=con_rd,con_rv=con_rv, &
                  con_cp=con_cp,con_ttp=con_ttp,con_hvap=con_hvap,con_hfus=con_hfus)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v17_p8_ugwpv1_pumas_init_cap(one=one,physics=physics,cdata=cdata,con_solr_2008=con_solr_2008,con_solr_2002=con_solr_2002, &
                  con_pi=con_pi,con_c=con_c,con_plnk=con_plnk,con_boltz=con_boltz,con_t0c=con_t0c, &
                  con_rd=con_rd,con_g=con_g,LTP=LTP,con_rerth=con_rerth,con_p0=con_p0,con_omega=con_omega, &
                  con_cp=con_cp,con_rv=con_rv,con_fvirt=con_fvirt,con_ttp=con_ttp,con_hvap=con_hvap, &
                  con_hfus=con_hfus)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_init

   subroutine ccpp_physics_run(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="SCM_GFS_v17_p8_ugwpv1_pumas") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_run_cap()
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_radiation_run_cap(one=one,physics=physics,cdata=cdata,LTP=LTP,con_eps=con_eps,con_pi=con_pi, &
                  con_rd=con_rd,con_g=con_g,con_ttp=con_ttp,con_thgni=con_thgni,con_epsm1=con_epsm1, &
                  con_fvirt=con_fvirt,con_rog=con_rog,con_rocp=con_rocp)
            else if (trim(group_name)=="phys_ps") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_run_cap(one=one,physics=physics,cdata=cdata,con_fvirt=con_fvirt,con_g=con_g,con_tice=con_tice, &
                  con_cp=con_cp,con_pi=con_pi,con_sbc=con_sbc,con_hvap=con_hvap,con_eps=con_eps, &
                  con_epsm1=con_epsm1,con_hfus=con_hfus,con_jcal=con_jcal,con_rd=con_rd,con_rhw0=con_rhw0, &
                  rlapse=rlapse,rhowater=rhowater,con_t0c=con_t0c,con_rocp=con_rocp,karman=karman, &
                  con_rv=con_rv,con_epsq=con_epsq)
            else if (trim(group_name)=="phys_ts") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_run_cap(one=one,physics=physics,cdata=cdata,con_1ovg=con_1ovg,con_fvirt=con_fvirt, &
                  con_cliq=con_cliq,con_cp=con_cp,con_cvap=con_cvap,con_eps=con_eps,con_epsm1=con_epsm1, &
                  con_g=con_g,con_hvap=con_hvap,con_rd=con_rd,con_rv=con_rv,con_t0c=con_t0c, &
                  con_pi=con_pi,rainmin=rainmin,rhowater=rhowater)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v17_p8_ugwpv1_pumas_run_cap(one=one,physics=physics,cdata=cdata,LTP=LTP,con_eps=con_eps,con_pi=con_pi, &
                  con_rd=con_rd,con_g=con_g,con_ttp=con_ttp,con_thgni=con_thgni,con_epsm1=con_epsm1, &
                  con_fvirt=con_fvirt,con_rog=con_rog,con_rocp=con_rocp,con_tice=con_tice, &
                  con_cp=con_cp,con_sbc=con_sbc,con_hvap=con_hvap,con_hfus=con_hfus,con_jcal=con_jcal, &
                  con_rhw0=con_rhw0,rlapse=rlapse,rhowater=rhowater,con_t0c=con_t0c,karman=karman, &
                  con_rv=con_rv,con_epsq=con_epsq,con_1ovg=con_1ovg,con_cliq=con_cliq,con_cvap=con_cvap, &
                  rainmin=rainmin)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_run

   subroutine ccpp_physics_finalize(cdata, suite_name, group_name, ierr)

      use ccpp_types, only : ccpp_t

      implicit none

      type(ccpp_t),               intent(inout) :: cdata
      character(len=*),           intent(in)    :: suite_name
      character(len=*), optional, intent(in)    :: group_name
      integer,                    intent(out)   :: ierr

      ierr = 0


      if (trim(suite_name)=="SCM_GFS_v17_p8_ugwpv1_pumas") then

         if (present(group_name)) then

            if (trim(group_name)=="time_vary") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_time_vary_final_cap(one=one,physics=physics,cdata=cdata)
            else if (trim(group_name)=="radiation") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_radiation_final_cap()
            else if (trim(group_name)=="phys_ps") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ps_final_cap(one=one,physics=physics,cdata=cdata)
            else if (trim(group_name)=="phys_ts") then
               ierr = SCM_GFS_v17_p8_ugwpv1_pumas_phys_ts_final_cap(one=one,physics=physics,cdata=cdata)
            else
               write(cdata%errmsg, '(*(a))') 'Group ' // trim(group_name) // ' not found'
               ierr = 1
            end if

         else

           ierr = SCM_GFS_v17_p8_ugwpv1_pumas_final_cap(one=one,physics=physics,cdata=cdata)

         end if

      else

         write(cdata%errmsg,'(*(a))') 'Invalid suite ' // trim(suite_name)
         ierr = 1

      end if

      cdata%errflg = ierr

   end subroutine ccpp_physics_finalize

end module ccpp_static_api
