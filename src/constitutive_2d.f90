!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : n_eqns , n_vars
  USE parameters_2d, ONLY : rheology_flag , rheology_model
  USE parameters_2d, ONLY : temperature_flag
  USE parameters_2d, ONLY : solid_transport_flag
  USE parameters_2d, ONLY : sed_vol_perc

  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  CHARACTER(LEN=20) :: phase1_name
  CHARACTER(LEN=20) :: phase2_name

  COMPLEX*16 :: h       !< height [m]
  COMPLEX*16 :: u       !< velocity (x direction)
  COMPLEX*16 :: v       !< velocity (y direction)
  COMPLEX*16 :: T       !< temperature
  COMPLEX*16 :: alphas  !< sediment volume fraction
  COMPLEX*16 :: rho_m   !< mixture density

  !> gravitational acceleration
  REAL*8 :: grav

  !> drag coefficients (Voellmy-Salm model)
  REAL*8 :: mu
  REAL*8 :: xi
  
  !> drag coefficients (plastic model)
  REAL*8 :: tau

  !> evironment temperature [K]
  REAL*8 :: T_env

  !> radiative coefficient
  REAL*8 :: rad_coeff

  !> friction coefficient
  REAL*8 :: frict_coeff

  !> fluid density [kg/m3]
  REAL*8 :: rho

  !> reference temperature [K]
  REAL*8 :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL*8 :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL*8 :: visc_par

  !> velocity boundary layer fraction of total thickness
  REAL*8 :: emme

  !> specific heat [J kg-1 K-1]
  REAL*8 :: c_p

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL*8 :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL*8 :: exp_area_fract

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL*8, PARAMETER :: SBconst = 5.67D-8

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL*8 :: emissivity

  !> thermal boundary layer fraction of total thickness
  REAL*8 :: enne

  !> temperature of lava-ground interface
  REAL*8 :: T_ground

  !> thermal conductivity [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL*8 :: thermal_conductivity

  !--- Lahars rheology model parameters

  !> 1st parameter for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha2

  !> 2nd parameter for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta2

  !> 1st parameter for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha1

  !> 2nd parameter for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta1

  !> Empirical resistance parameter
  REAL*8 :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL*8 :: n_td
 
  !> Specific weight of water
  REAL*8 :: rho_w

  !> Specific weight of sediments
  REAL*8 :: rho_s

  !> hindered settling 
  REAL*8 :: settling_vel

  !> erosion model coefficient
  REAL*8 :: erosion_coeff


CONTAINS

  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters_2d, ONLY : n_nh

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.
    implicit_flag(3) = .TRUE.

    IF ( solid_transport_flag ) THEN

       IF ( temperature_flag ) implicit_flag(5) = .TRUE.

    ELSE

       IF ( temperature_flag ) implicit_flag(4) = .TRUE.

    END IF

    n_nh = COUNT( implicit_flag )

  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h+B, u, v, xs , T \f$).
  !> \param[in]    r_qj     real conservative variables 
  !> \param[in]    c_qj     complex conservative variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE phys_var(r_qj,c_qj)
    
    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none
    
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)

    COMPLEX*16 :: qj(n_vars)

    IF ( present(c_qj) ) THEN

       qj = c_qj

    ELSE

       qj = DCMPLX(r_qj)

    END IF

    h = qj(1)

    IF ( solid_transport_flag ) THEN

       IF ( DBLE( h ) .GT. 1.D-25 ) THEN
          
          alphas = qj(4) / h
          
          IF ( temperature_flag ) T = qj(5) / h
          
       ELSE
          
          alphas = DSQRT(2.D0) * h * qj(4) / CDSQRT( h**4 + eps_sing )
          
          IF ( temperature_flag ) THEN
             
             T =  DSQRT(2.D0) * h * qj(5) / CDSQRT( h**4 + eps_sing )
             
          END IF
          
       END IF

       rho_m = ( DCMPLX(1.D0,0.D0) - alphas ) * rho_w + alphas * rho_s 

    ELSE

       IF ( temperature_flag ) THEN

          IF ( DBLE( h ) .GT. 1.D-25 ) THEN
          
             T = qj(4) / h
          
          ELSE
          
             T =  DSQRT(2.D0) * h * qj(4) / CDSQRT( h**4 + eps_sing )
             
          END IF
          
       END IF

       alphas = DCMPLX(0.D0,0.D0) 
       rho_m = DCMPLX(1.D0,0.D0) 

    END IF

    IF ( DBLE( h ) .GT. eps_sing ** 0.25D0 ) THEN

       u = qj(2) / ( rho_m * h )

       v = qj(3) / ( rho_m * h )

    ELSE

       u = DSQRT(2.D0) * h * ( qj(2) / rho_m ) / CDSQRT( h**4 + eps_sing )

       v = DSQRT(2.D0) * h * ( qj(3) / rho_m ) / CDSQRT( h**4 + eps_sing )

    END IF

    RETURN

  END SUBROUTINE phys_var

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine desingularize the velocities and then evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qj,vel_min,vel_max)
  
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)

    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    CALL phys_var(r_qj = qj)
    
    vel_min(1:n_eqns) = DBLE(u) - DSQRT( grav * DBLE(h) )
    vel_max(1:n_eqns) = DBLE(u) + DSQRT( grav * DBLE(h) )
        
  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine desingularize the velocities and then evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)
    
    CALL phys_var(r_qj = qj)
    
    vel_min(1:n_eqns) = DBLE(v) - DSQRT( grav * DBLE(h) )
    vel_max(1:n_eqns) = DBLE(v) + DSQRT( grav * DBLE(h) )
       
  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h+B \f$
  !> - qp(2) = \f$ u \f$
  !> - qp(3) = \f$ v \f$
  !> - qp(4) = \f$ xs \f$
  !> - qp(5) = \f$ T \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qp      physical variables  
  !> \date 15/08/2011
  !******************************************************************************
  
  SUBROUTINE qc_to_qp(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(r_qj = qc)
    
    ! It is important to have as reconstructed variable at the interfaces
    ! the total height h+B, instead of h. Only in this way, when the topography
    ! is not flat and the free surface is horizontal (steady condition), the  
    ! latter is kept flat by the reconstruction and the solution is stable.
    qp(1) = DBLE(h+B)
    qp(2) = DBLE(u)
    qp(3) = DBLE(v)

    IF ( solid_transport_flag ) THEN
       
       qp(4) = DBLE(alphas)
       
       IF ( temperature_flag ) qp(5) = DBLE(T)
       
    ELSE
       
       IF ( temperature_flag ) qp(4) = DBLE(T)
       
    END IF
    
    
  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h + B \f$
  !> - qp(2) = \f$ rho_m*u \f$
  !> - qp(3) = \f$ rho_m*v \f$
  !> - qp(4) = \f$ alphas \f$
  !> - qp(5) = \f$ T \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_h       !> height 
    REAL*8 :: r_u       !> velocity
    REAL*8 :: r_v       !> velocity
    REAL*8 :: r_T       !> temperature
    REAL*8 :: r_rho_m !> mixture density
    REAL*8 :: r_alphas  !> sediment volume fraction
 
    r_h = qp(1) - B
    r_u  = qp(2)
    r_v  = qp(3)


    qc(1) = r_h

    IF ( solid_transport_flag ) THEN
       
       r_alphas = qp(4)
       qc(4) = r_h * r_alphas 
       
       IF ( temperature_flag ) THEN
          
          r_T  = qp(5)
          qc(5) = r_h * r_T 
       
       END IF
       
       ! TO BE REMOVED
       ! r_alphas = sed_vol_perc/100.D0

       r_rho_m = r_alphas * rho_s + ( 1.D0 - r_alphas ) * rho_w

    ELSE

       IF ( temperature_flag ) THEN
          
          r_T  = qp(4)
          qc(4) = r_h * r_T 
       
       END IF

       r_alphas = 0.D0
       r_rho_m = 1.D0

    END IF

    qc(2) = r_h * r_rho_m * r_u
    qc(3) = r_h * r_rho_m * r_v


  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Conservative to alternative conservative variables
  !
  !> This subroutine evaluates from the conservative variables qc the array of 
  !> alternative conservative variables qc2:\n
  !> - qc2(1) = \f$ h+B \f$
  !> - qc2(2) = \f$ rho_m hu \f$
  !> - qc2(3) = \f$ rho_m hv \f$
  !> - qc2(4) = \f$ h alpha_s \f$
  !> - qc2(5) = \f$ hT \f$
  !> .
  !> The alternative conservative  variables are those used for the linear 
  !> reconstruction at the cell interfaces. Limiters are applied to the 
  !> reconstructed slopes.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qc2     alternative conservative variables  
  !> \date 02/04/2019
  !******************************************************************************

  SUBROUTINE qc_to_qc2(qc,B,qc2)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc2(n_vars)

    qc2(1) = qc(1) + B
    qc2(2:n_vars) = qc(2:n_vars)

  END SUBROUTINE qc_to_qc2

  !******************************************************************************
  !> \brief Reconstructed to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the array
  !> of alternative conser variables qc2, reconstructed at the interfaces:\n
  !> - qc2(1) = \f$ h + B \f$
  !> - qc2(2) = \f$ rho_m hu \f$
  !> - qc2(3) = \f$ rho_m hv \f$
  !> - qc2(4) = \f$ h \cdot alpha_s \f$
  !> - qc2(5) = \f$ h \cdot T \f$
  !> .
  !> \param[in]    qc2     alternative conservative variables  
  !> \param[out]   qc      conservative variables 
  !> \date 02/04/2019
  !******************************************************************************

  SUBROUTINE qc2_to_qc(qc2,B,qc)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc2(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)

    REAL*8 :: qc_temp(n_vars)

    REAL*8 :: r_h       !> height 
    REAL*8 :: r_u       !> velocity
    REAL*8 :: r_v       !> velocity
    REAL*8 :: r_alphas  !> sediment volume fraction
    REAL*8 :: r_T       !> temperature
    REAL*8 :: r_rho_m

    qc_temp(1) = qc2(1) - B
    qc_temp(2:n_vars) = qc2(2:n_vars)

    ! Desingularization
    CALL phys_var(r_qj = qc_temp)
          
    r_h = DBLE(h)
    r_u = DBLE(u)
    r_v = DBLE(v)
    r_alphas = DBLE(alphas)

    ! TO BE REMOVED
    ! r_alphas = sed_vol_perc/100.D0

    IF ( solid_transport_flag ) THEN

       r_rho_m = DBLE(rho_m)

    ELSE

       r_rho_m = 1.D0

    END IF

    qc(1) = r_h
    qc(2) = r_h * r_rho_m * r_u
    qc(3) = r_h * r_rho_m * r_v
    
    IF ( solid_transport_flag ) THEN
       
       qc(4) = r_h * r_alphas
       
       IF ( temperature_flag ) THEN
          
          r_T = DBLE(T)
          qc(5) = r_h * r_T 
          
       END IF
       
    ELSE
       
       IF ( temperature_flag ) THEN
          
          r_T = DBLE(T)
          qc(4) = r_h * r_T 
          
       END IF
       
    END IF
    
  END SUBROUTINE qc2_to_qc


  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qj, accordingly to the equations for the single temperature
  !> model introduced in Romenki et al. 2010.
  !> \date 01/06/2012
  !> \param[in]     c_qj     complex conservative variables 
  !> \param[in]     r_qj     real conservative variables 
  !> \param[out]    c_flux   complex analytical fluxes    
  !> \param[out]    r_flux   real analytical fluxes    
  !******************************************************************************
  
  SUBROUTINE eval_fluxes(c_qj,r_qj,c_flux,r_flux,dir)
    
    USE COMPLEXIFY
    IMPLICIT none

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_flux(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_flux(n_eqns)
    INTEGER, INTENT(IN) :: dir

    COMPLEX*16 :: qj(n_vars)
    COMPLEX*16 :: flux(n_eqns)
    COMPLEX*16 :: h_temp , u_temp , v_temp
    COMPLEX*16 :: alphas_temp , rho_m_temp

    INTEGER :: i 

    IF ( present(c_qj) .AND. present(c_flux) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_flux) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    h_temp = qj(1)

    pos_thick:IF ( DBLE(h_temp) .NE. 0.D0 ) THEN

       IF ( solid_transport_flag ) THEN
          
          alphas_temp = qj(4) / h_temp
          
          ! TO BE REMOVED
          ! alphas_temp = DCMPLX(sed_vol_perc/100.D0,0.D0)
          
          rho_m_temp = alphas_temp * rho_s + ( DCMPLX(1.D0,0.D0)                 &
               - alphas_temp ) * rho_w
          
       ELSE
          
          alphas_temp = DCMPLX( 0.D0 , 0.D0 )
          rho_m_temp = DCMPLX(1.D0,0.D0) 
          
       END IF
       
       IF ( dir .EQ. 1 ) THEN

          ! Velocity in x-direction
          u_temp = qj(2) / ( h_temp * rho_m_temp )

          ! Volumetric flux in x-direction: h*u
          flux(1) = qj(2) / rho_m_temp

          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = rho_m_temp * h_temp * u_temp**2 +                           &
               0.5D0 * rho_m_temp * grav * h_temp**2

          ! y-momentum flux in x-direction: rho*h*v*u
          flux(3) = u_temp * qj(3)

          IF ( solid_transport_flag ) THEN

             ! Volumetric flux of solid in x-direction: h * alphas * u
             flux(4) = u_temp * qj(4)

             ! Solid flux can't be larger than total flux
             IF ( flux(4) / flux(1) .GT. 1.D0 ) flux(4) = flux(1)
             
             ! Temperature flux in x-direction: U*T*h
             IF ( temperature_flag ) flux(5) = u_temp * qj(5)
             
          ELSE

             ! Temperature flux in x-direction: U*T*h
             IF ( temperature_flag ) flux(4) = u_temp * qj(4)
              
          END IF

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)

          v_temp = qj(3) / ( h_temp * rho_m_temp )
          
          flux(1) = qj(3) / rho_m_temp
          
          flux(2) = v_temp * qj(2)
          
          flux(3) = rho_m_temp * h_temp * v_temp**2 +                           &
               0.5D0 * rho_m_temp * grav * h_temp**2
          
          IF ( solid_transport_flag ) THEN

             flux(4) = v_temp * qj(4)

             ! Solid flux can't be larger than total flux
             IF ( flux(4) / flux(1) .GT. 1.D0 ) flux(4) = flux(1)

             ! Temperature flux in y-direction: V*T*h
             IF ( temperature_flag ) flux(5) = v_temp * qj(5)
             
          ELSE

             ! Temperature flux in y-direction: V*T*h
             IF ( temperature_flag ) flux(4) = v_temp * qj(4)

          END IF

       END IF

    ELSE
       
       flux(1:n_eqns) = 0.D0
       
    ENDIF pos_thick
 
    !WRITE(*,*) 'flux',DBLE(flux)
    !WRITE(*,*) 'qj',DBLE(qj)
    !WRITE(*,*) 'Bj',Bj
    !WRITE(*,*) 'rho_m_temp , h_temp', REAL(rho_m_temp) , REAL(h_temp)
    !READ(*,*)

    IF ( present(c_qj) .AND. present(c_flux) ) THEN

       c_flux = flux

    ELSEIF ( present(r_qj) .AND. present(r_flux) ) THEN

       r_flux = DBLE( flux )

    END IF

  END SUBROUTINE eval_fluxes

  !******************************************************************************
  !> \brief Non-Hyperbolic terms
  !
  !> This subroutine evaluates the non-hyperbolic terms (relaxation terms
  !> and forces) of the system of equations, both for real or complex 
  !> inputs. These terms are treated implicitely in the DIRK numerical
  !> scheme.
  !> \date 01/06/2012
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nonhyperbolic_terms( c_qj , c_nh_term_impl , r_qj ,           &
       r_nh_term_impl )

    USE COMPLEXIFY 

    IMPLICIT NONE

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)

    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX*16 :: qj(n_vars)

    COMPLEX*16 :: nh_term(n_eqns)

    COMPLEX*16 :: relaxation_term(n_eqns)

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i

    COMPLEX*16 :: mod_vel
    
    COMPLEX*16 :: gamma

    REAL*8 :: radiative_coeff

    COMPLEX*16 :: radiative_term

    REAL*8 :: convective_coeff

    COMPLEX*16 :: convective_term

    COMPLEX*16 :: conductive_coeff , conductive_term

    REAL*8 :: thermal_diffusivity

    REAL*8 :: h_threshold

    !--- Lahars rheology model variables
    
    !> Fluid viscosity
    COMPLEX*8 :: fluid_visc

    !> Total friction
    COMPLEX*8 :: s_f

    !> Viscous slope component of total Friction
    COMPLEX*8 :: s_v

    !> Turbulent dispersive slope component of total friction
    COMPLEX*8 :: s_td


    IF ( temperature_flag ) THEN

       IF ( ( thermal_conductivity .GT. 0.D0 ) .OR. ( emme .GT. 0.D0 ) ) THEN

          h_threshold = 1.D-10

       ELSE

          h_threshold = 0.D0

       END IF
       
    END IF
       

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the relaxation terms
    relaxation_term(1:n_eqns) = DCMPLX(0.D0,0.D0) 

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    IF (rheology_flag) THEN

       CALL phys_var(c_qj = qj)
    
       mod_vel = CDSQRT( u**2 + v**2 )
       
       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN
       
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
          
             ! IMPORTANT: grav3_surv is always negative 
             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  ( grav / xi ) * mod_vel ** 2
          
             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  ( grav / xi ) * mod_vel ** 2
          
          ENDIF
        
       ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN
       
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
          
             forces_term(2) = forces_term(2) - rho_m * tau * (u/mod_vel)
          
             forces_term(3) = forces_term(3) - rho_m * tau * (v/mod_vel)

          ENDIF

       ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( DBLE(h) .GT. h_threshold ) THEN
    
             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.D0 * nu_ref / h * CDEXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.D0 * nu_ref / h_threshold * CDEXP( - visc_par            &
                  * ( T - T_ref ) )
             
          END IF
          
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
          
             ! Last R.H.S. term in equation 2 from Costa & Macedonio, 2005
             forces_term(2) = forces_term(2) - rho_m * gamma * u
          
             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             forces_term(3) = forces_term(3) - rho_m * gamma * v

          ENDIF
          
       ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.D-20

          ! THIS IS ONLY USED WHEN SEDIMENT PERCENTAGE IS GIVEN AS INPUT AND
          ! IT IS CONSTANT WITHIN THE SIMULATION
          ! TO BE REMOVED
          ! alphas = DCMPLX(sed_vol_perc/100.D0,0.D0)

          ! Fluid viscosity
          fluid_visc = alpha1 * CDEXP( beta1 * alphas )

          IF ( h .GT. h_threshold ) THEN
             
             ! Viscous slope component
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * h**2 )
             
             ! Turbulent dispersive component
             s_td = rho_m * n_td**2 * mod_vel**2 / ( h**(4.D0/3.D0) )
          
          ELSE
             
             ! Viscous slope component
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * h_threshold**2 )
             
             ! Turbulent dispersive components
             s_td = rho_m * n_td**2 * (mod_vel**2) / ( h_threshold**(4.D0/3.D0) )
             
          END IF
          
          ! Total implicit friction slope
          s_f = s_v + s_td
         
          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - grav * h * ( u / mod_vel ) * s_f

             forces_term(3) = forces_term(3) - grav * h * ( v / mod_vel ) * s_f
  
          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          tau = 1.D-3 / ( 1.D0 + 10.D0 * h ) * mod_vel
          
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN
             
             forces_term(2) = forces_term(2) - rho_m * tau * ( u / mod_vel )
             forces_term(3) = forces_term(3) - rho_m * tau * ( v / mod_vel )

          END IF
          
       ENDIF
              
    ENDIF

    IF ( temperature_flag ) THEN

       CALL phys_var(c_qj = qj)

       IF ( DBLE(h) .GT. 0.d0 ) THEN

          ! Equation 8 from Costa & Macedonio, 2005
          radiative_coeff = emissivity * SBconst * exp_area_fract / ( rho * c_p )
       
       ELSE

          radiative_coeff = 0.D0

       END IF

       IF ( DBLE(T) .GT. T_env ) THEN

          ! First R.H.S. term in equation 4 from Costa & Macedonio, 2005
          radiative_term = - radiative_coeff * ( T**4 - T_env**4 )

       ELSE

          radiative_term = DCMPLX(0.D0,0.D0)

       END IF

       IF ( DBLE(h) .GT. 0.d0 ) THEN
       
          ! Equation 9 from Costa & Macedonio, 2005
          convective_coeff = atm_heat_transf_coeff * exp_area_fract             &
               / ( rho * c_p )

       ELSE

          convective_coeff = 0.D0

       END IF

       IF ( DBLE(T) .GT. T_env ) THEN

          ! Second R.H.S. term in equation 4 from Costa & Macedonio, 2005
          convective_term = - convective_coeff * ( T - T_env )

       ELSE

          convective_term =  DCMPLX(0.D0,0.D0)

       END IF

       IF ( DBLE(h) .GT. h_threshold ) THEN
    
          thermal_diffusivity = thermal_conductivity / ( rho * c_p ) 

          ! Equation 7 from Costa & Macedonio, 2005
          conductive_coeff = enne * thermal_diffusivity / h

       ELSE

          conductive_coeff =  DCMPLX(0.D0,0.D0)
          conductive_coeff = enne * thermal_diffusivity / DCMPLX(h_threshold,0.D0)

       END IF

       ! Third R.H.S. term in equation 4 from Costa & Macedonio, 2005
       IF ( DBLE(T) .GT. T_ground ) THEN

          conductive_term = - conductive_coeff * ( T - T_ground )

       ELSE

           conductive_term = DCMPLX(0.D0,0.D0)

        END IF

        IF ( solid_transport_flag ) THEN

           relaxation_term(5) = radiative_term + convective_term + conductive_term

        ELSE

           relaxation_term(4) = radiative_term + convective_term + conductive_term

        END IF

    END IF

    nh_term = relaxation_term + forces_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = DBLE( nh_term )

    END IF 

  END SUBROUTINE eval_nonhyperbolic_terms

  !******************************************************************************
  !> \brief Non-Hyperbolic semi-implicit terms
  !
  !> This subroutine evaluates the non-hyperbolic terms that are solved
  !> semi-implicitely by the solver. For example, any discontinuous term that
  !> appears in the friction terms.
  !> \date 20/01/2018
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nh_semi_impl_terms( grav3_surf , c_qj , c_nh_semi_impl_term , &
       r_qj , r_nh_semi_impl_term )

    USE COMPLEXIFY 


    IMPLICIT NONE

    REAL*8, INTENT(IN) :: grav3_surf

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_semi_impl_term(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_semi_impl_term(n_eqns)

    COMPLEX*16 :: qj(n_vars)

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i

    COMPLEX*16 :: mod_vel
    
    REAL*8 :: h_threshold

    !--- Lahars rheology model variables
    
    !> Yield strenght
    COMPLEX*8 :: tau_y

    !> Yield slope component of total friction
    COMPLEX*8 :: s_y


    IF ( present(c_qj) .AND. present(c_nh_semi_impl_term) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_semi_impl_term) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    IF (rheology_flag) THEN

       CALL phys_var(c_qj = qj)
    
       mod_vel = CDSQRT( u**2 + v**2 )

       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN

          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  mu * h * ( - grav * grav3_surf )
             
             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  mu * h * ( - grav * grav3_surf )
             
          END IF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN
          

       ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          
       ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.D-20

          ! THIS IS ONLY USED WHEN SEDIMENT PERCENTAGE IS GIVEN AS INPUT AND
          ! IT IS CONSTANT WITHIN THE SIMULATION
          ! alphas = DCMPLX( sed_vol_perc*1.D-2 , 0.D0 )

          ! Yield strength
          tau_y = alpha2 * CDEXP( beta2 * alphas )

          IF ( h .GT. h_threshold ) THEN
             
             ! Yield slope component
             s_y = tau_y / h
                       
          ELSE
             
             ! Yield slope component
             s_y = tau_y / h_threshold
                          
          END IF
  
          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - grav * h * ( u / mod_vel ) * s_y

             forces_term(3) = forces_term(3) - grav * h * ( v / mod_vel ) * s_y
  
          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          
       ENDIF
              
    ENDIF

    IF ( temperature_flag ) THEN


    END IF

    
    IF ( present(c_qj) .AND. present(c_nh_semi_impl_term) ) THEN

       c_nh_semi_impl_term = forces_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_semi_impl_term) ) THEN

       r_nh_semi_impl_term = DBLE( forces_term )

    END IF 

  END SUBROUTINE eval_nh_semi_impl_terms

  
  !******************************************************************************
  !> \brief Explicit Forces term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the left-hand side of the equations.
  !> \date 01/06/2012
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    expl_term          explicit term
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x , Bprimej_y , source_xy , qj ,          &
       expl_term )

    USE parameters_2d, ONLY : source_flag , vel_source , T_source
    
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: source_xy
    
    REAL*8, INTENT(IN) :: qj(n_eqns)                 !< conservative variables 
    REAL*8, INTENT(OUT) :: expl_term(n_eqns)         !< explicit forces 

    REAL*8 :: visc_heat_coeff

    expl_term(1:n_eqns) = 0.D0

    CALL phys_var(r_qj = qj)

    IF ( source_flag ) expl_term(1) = - source_xy * vel_source
    
    expl_term(2) = grav * DBLE(rho_m) * DBLE(h) * Bprimej_x
   
    expl_term(3) = grav * DBLE(rho_m) * DBLE(h) * Bprimej_y

    IF ( temperature_flag .AND. source_flag ) THEN

       IF ( solid_transport_flag ) THEN
    
          expl_term(5) = - source_xy * vel_source * T_source

       ELSE

          expl_term(4) = - source_xy * vel_source * T_source

       END IF

       IF ( rheology_model .EQ. 3 ) THEN
              
          IF ( DBLE(h) .GT. 0.D0 ) THEN
             
             ! Equation 10 from Costa & Macedonio, 2005
             visc_heat_coeff = emme * nu_ref / ( c_p * DBLE(h) ) 
             
          ELSE
             
             visc_heat_coeff = 0.D0
             
          END IF

          IF ( solid_transport_flag ) THEN
                          
             ! Viscous heating
             ! Last R.H.S. term in equation 4 from Costa & Macedonio, 2005
             expl_term(5) = expl_term(5) - visc_heat_coeff * ( DBLE(u)**2       &
                  + DBLE(v)**2 ) * DEXP( - visc_par * ( DBLE(T) - T_ref ) ) 
             
          ELSE

             ! Viscous heating
             ! Last R.H.S. term in equation 4 from Costa & Macedonio, 2005
             expl_term(4) = expl_term(4) - visc_heat_coeff * ( DBLE(u)**2       &
                  + DBLE(v)**2 ) * DEXP( - visc_par * ( DBLE(T) - T_ref ) ) 
             
          END IF

       END IF
          
    END IF
           
  END SUBROUTINE eval_expl_terms


  !******************************************************************************
  !> \brief Deposition term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 03/010/2018
  !> \param[in]     Bj                 topography
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    dep_term          explicit term
  !******************************************************************************

  SUBROUTINE eval_erosion_dep_term( qj , erosion_term , deposition_term )
    
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: qj(n_eqns)                 !< conservative variables 
    
    REAL*8, INTENT(OUT) :: erosion_term              !< erosion term
    REAL*8, INTENT(OUT) :: deposition_term           !< deposition term
    
    REAL*8 :: r_alphas


    !>  settling speed of the single grain
    !REAL*8 :: settling_single_vel

    !> Richardson-Zaki exponent [Richardson and Zaki, 1954]
    !REAL*8 :: RZ_exp

    !> settling particle Reynold number
    !REAL*8 :: Rey

    !> kinematic viscosity of the fluid
    !REAL*8 :: kin_visc

    !> grain diameter
    !REAL*8 :: diam

    !> dimensionless grain parameter
    !REAL*8 :: R

    REAL*8 :: mod_vel

    !diam = 1.D-3

    CALL phys_var(r_qj = qj)

    deposition_term = 0.D0
    erosion_term = 0.D0
    
    r_alphas = DBLE(alphas)

    !Rey = DSQRT( grav * diam**3 ) / kin_visc

    !RZ_exp = ( 4.7D0 + 0.41D0 * Rey**0.75D0 ) / ( 1.D0 + 0.175D0 * Rey**0.75D0 )

    !settling_single_vel = kin_visc / diam * 

    !settling_vel = ( 1.D0 - r_alphas )**(2.7D0 - 0.15D0*Rz_exp )

    deposition_term = r_alphas * settling_vel

    mod_vel = DSQRT( DBLE(u)**2 + DBLE(v)**2 )
  
    IF ( DBLE(h) .GT. 1.D-2) THEN
    
       !erosion_term = erosion_coeff * mod_vel * DBLE(h) * ( rho_s - DBLE(rho_m) )  &
       !     / ( rho_s - rho_w )
       
       erosion_term = erosion_coeff * mod_vel * DBLE(h) * ( 1.D0 - r_alphas )

    ELSE

       erosion_term = 0.D0

    END IF
    
    !WRITE(*,*) 'r_alphas, settling_vel',r_alphas, settling_vel
    !WRITE(*,*) 'eval_erosion_dep_term',erosion_term,deposition_term
    !READ(*,*)

    RETURN
  
  END SUBROUTINE eval_erosion_dep_term


  !******************************************************************************
  !> \brief Topography modification related term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 03/010/2018
  !> \param[in]     Bj                  topography
  !> \param[in]     qj                  conservative variables 
  !> \param[in]     deposition_avg_term averaged deposition terms 
  !> \param[in]     erosion_avg_term    averaged deposition terms 
  !> \param[out]    topo_term           explicit term
  !******************************************************************************

  SUBROUTINE eval_topo_term( qj , deposition_avg_term , erosion_avg_term ,      &
       eqns_term, topo_term )
    
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: qj(n_eqns)                  !< conservative variables 
    REAL*8, INTENT(IN) :: deposition_avg_term         !< deposition term
    REAL*8, INTENT(IN) :: erosion_avg_term            !< erosion term

    REAL*8, INTENT(OUT):: eqns_term(n_eqns)
    REAL*8, INTENT(OUT):: topo_term
    

    CALL phys_var(r_qj = qj)

    eqns_term(1:n_eqns) = 0.D0

    ! free surface (topography+flow) equation
    eqns_term(1) = erosion_avg_term - deposition_avg_term

    ! x-momenutm equation
    eqns_term(2) = - DBLE(u) * rho_s * deposition_avg_term

    ! y-momentum equation
    eqns_term(3) = - DBLE(v) * rho_s * deposition_avg_term

    ! solid phase thickness equation
    eqns_term(4) = erosion_avg_term - deposition_avg_term
  
    ! topography term
    topo_term = - erosion_avg_term + deposition_avg_term

  END SUBROUTINE eval_topo_term


END MODULE constitutive_2d

    
