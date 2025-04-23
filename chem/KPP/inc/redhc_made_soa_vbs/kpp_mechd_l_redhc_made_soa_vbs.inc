REAL(KIND=dp) :: roleho,roleo3,roleno3,rtolho,    &
                 rsoaalkho,rivocho,risoho,risoo3,risono3,           &
                 rc10h16ho,rc10h16o3,rc10h16no3,                  &
                 rsesqho,rsesqo3,rsesqno3,       &
                 ro2loss,rxo2no,rxo2xo2,rxo2ho2


REAL(KIND=dp) , DIMENSION(ldrog_soa2)  :: PRDROG

REAL(kind=dp) :: es, qvs, rh

REAL(KIND=dp) :: RATE_CVPOA
REAL(KIND=dp) :: RATE_CVASOA
REAL(KIND=dp) :: RATE_CVASOA0
REAL(KIND=dp) :: RATE_CVBSOA

INTEGER, SAVE :: icall

REAL(kind=dp), parameter :: rho_sulf = 1.7e3_dp 
REAL(kind=dp), parameter :: rho_orgc = 1.8e3_dp 
REAL(kind=dp), parameter :: rho_bc   = 1.0e3_dp 
REAL(kind=dp) :: aer_so4, aer_oc2, aer_bc2


RATE_CVPOA = REAL( config_flags%rate_cvpoa, KIND=dp )
RATE_CVASOA= REAL( config_flags%rate_cvasoa, KIND=dp )
RATE_CVASOA0= REAL( config_flags%rate_cvasoa0, KIND=dp )
RATE_CVBSOA= REAL( config_flags%rate_cvbsoa, KIND=dp )
