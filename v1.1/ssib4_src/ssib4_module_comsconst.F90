MODULE ssib4_module_comsconst
!=======================================================================       
!                                                                               
!     ssib common block : 1-d version 
!                                                         12 august 2000       
!=======================================================================    
!        physical constants                                                     
!-----------------------------------------------------------------------       
  REAL, PARAMETER      ::    CPAIR    = 1004.6                  &
                            ,STEFAN   = 5.669 * 10E-9           &
                            ,GRAV     = 9.81                    &
                            ,VKC      = 0.4                     &
                            ,PIE      = 3.14159265              &
                            ,TIMCON   = PIE/86400.              &
                            ,CLAI     = 4.2 * 1000. * 0.2       &
                            ,CW       = 4.2 * 1000. * 1000.     &
                            ,TF       = 273.16                  &
                            ,GASR     = 287.05                  &
!-----------------------------------------------------------------------
!     N.B. : HLAT IS EXPRESSED IN J KG-1                                
!            SNOMEL IS EXPRESSED IN J M-1                               
!-----------------------------------------------------------------------
                            ,HLAT     = 2.52E6                  &
                            ,SNOMEL   = 370518.5 * 1000.
  INTEGER, PARAMETER   ::    ITRUNK   = 3

  real, parameter :: cpinv   = 1.0/CPAIR
  real, parameter :: hvapi   = 1.0/HLAT
  real, parameter :: AKAPPA  = GASR/CPAIR
  real :: psy
!
END MODULE
