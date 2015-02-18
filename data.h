# define PI 3.141592653589793

/* NN defines the maximum number of trees in a quadrat */ 
# define NN 350

/* MAX_WAVE now defines the maximum number of optical property combination*/  
# define MAX_WAVE  3000 

// maximum of parametrs allowed ..

# define N_VZA 181  
# define N_SZA 5
# define N_PHI 2
# define N_OPTIC 10 
# define N_LAI 20  
# define N_D 20   
# define N_OMEGA 20
# define N_GAMMA 20
# define N_M2 20
# define N_HB 20 
# define N_HA 20
# define N_R 20
# define N_ANGLE 20
# define N_B 20
# define N_SHAPE 3 
# define N_QUADRAT 5











struct PARAMETER

{
	// Loop Parameters (to be read in a txt file)

	int NN_VZA; /* number of View zenith angle that will be read in input file */
	int NN_SZA;	/* number of solar zenith angle that will be read in input file */
	int NN_PHI;
	int NN_OPTIC; 
	int NN_LAI;  
	int NN_D;   
	int NN_OMEGA;
	int NN_GAMMA;
	int NN_M2;
	int NN_HB;
	int NN_HA;
	int NN_R;
	int NN_ANGLE;
	int NN_A;
	int NN_B; 
	int NN_SHAPE;
	int NN_QUADRAT;

	double OPTIC_NIRG[N_OPTIC];
	double OPTIC_NIRT[N_OPTIC]; 
	double OPTIC_NIRTT[N_OPTIC];
	double OPTIC_REDG[N_OPTIC]; 
	double OPTIC_REDT[N_OPTIC]; 
	double OPTIC_REDTT[N_OPTIC];


	double D_RANGE[N_D]; 
	double LAI_RANGE[N_LAI];
	double OMEGA_RANGE[N_OMEGA];
	double GAMMA_RANGE[N_GAMMA];  
	double M2_RANGE[N_M2];
	double HB_RANGE[N_HB];
	double HA_RANGE[N_HA];
	double R_RANGE[N_R];
	double ALPHA_G[N_ANGLE];
	double B_RANGE[N_B];
	double SHAPE_RANGE[N_SHAPE];
	double QUADRAT_RANGE[N_QUADRAT];


	double VZA_RANGE[N_VZA];
	double SZA_RANGE[N_SZA];
	double PHI_RANGE[N_PHI];



	/* LIBERTY PARAMETERS */
	double	m_D;
	double	m_XU;
	double	m_THICK;
	double	m_BASELINE;
	double	m_ELEMENT;
	double	m_C_FACTOR;
	double	m_L_FACTOR;
	double	m_P_FACTOR;
	double	m_W_FACTOR;

	char PIGMENT_FILE[120];
	char WATER_FILE[120];
	char ALBINO_FILE[120];
	char LIGCELL_FILE[120];
	char PROTEIN_FILE[120];
	int  LIBERTY_DEFAULT;

	double	A,C;                                        /* G = A* ZA + C */
	double	ALPHA;                                         /* apex angle */
	double	ALPHA_B,ALPHA_L;                    /* branch and leaf angle */
	double	ANGLE_RES;                      /* angle resolution of PLANE */
	double	B;                                          /* DOMAIN SIZE   */	
	double  Cp;                                  /* Used in Equ. 57 [1]  */
	double	D;                          /* Number of trees in the domain */	
	double  Band[10];
	double	DeltaLAI;                           /* Increase of LAI in MS */
	double	fB[10];      /* fraction of direct light (visible) 10 bands..*/ 
	double	Fr;              /* max Fo ...  1 == no overlaping at nadir */ 
	double	GAMMA_E;                            /* needle to shoot ratio */
	int		SPECTRAL; /* indicates if the model runs in 
						   0 : 4-scale mode with no multiple scattering (multispectral: 4 bands)
		                   1 : 4-scale mode with multiple scattering scheme (multispectral 4 bands)
						   2 : hyperspectral 
	                       3 : hyperspectral + LIBERTY mode 
						   4 : hyperspectral from file 
	                       5 : hyperspectral from file for background + LIBERTY for foliage */
	double	Ha;                                          /* stick height */
	double	Hb;                           /* cylinder or spheroid height */
 	double	Ll;                                /* Branch leaf area index */
	double	LAI;                                   /* LAI of the forest  */
	double	m2;                                /* Neyman grouping index  */	
	double	n;                                     /* number of quadrat  */
	double	OMEGA_E;                                   /* clumping index */
	double	PHI[1000];              /* contains the azimuthal differences */ 
	double	R;                                    /* radius of the crown */
	double	RATIO;                      /* foliage thickness/width ratio */
	double  Rb;               /* thickness of branches, used in ratio ...*/
	double	SZA;
	double	SZA2[1000];                              /* solar zenith angle */
	double  VZA[1000];
	double	VZA_MAX;                     /* maximum VZA in PLANE outpput */
	double  VZA_MIN;                     /* minimun VZA in PLANE outpput */
	double	Ws;                         /* typical foliage element width */
	double	G1,GZ1,T1,TZ1,G2,GZ2,T2,TZ2,G3,GZ3,T3,TZ3,G4,GZ4,T4,TZ4; 
					/* reflectivities */
		                      /* default: (1) RED
										  (2) NIR
							              (3) MIR (SWIR)   
	                                      (4) MIR (SWIR) */
 
	double	TT1,TT2,TT3,TT4;                              /* transmitivities */
	double  RG[MAX_WAVE];
	double  RT[MAX_WAVE];
	double  TT[MAX_WAVE];
	char    OUTPUT_FILE[128];
	char	COM_FILE[128];                                        /* input file name */
	char	COMMENTS[128];                            /* 128 characters comment line */
	char	ANGLE_FILE[128];     /* input file containing the angles to be simulated */
/*	char    OPTICAL_FILE[128];  Input spectra */
	char	GE_CHOICE[128];                                /* branch of G= A*ZA + C  */
	char	PLOT[128];                           /* Ploting style, PLANE, POLDER, 3D */
	char	SHAPE[128];                                /* SPHEROID of CONE+CYLINDER  */



	} ;   

 struct RESULT
	{
	double Medium; /* just for testing */
	double	Cs;                                           /* eq. 27 [2] */ 
	double	Cv;                                           /* eq. 28 [2] */
	double	DATA;                                /* POLDER measurements */ 
	double	DATA_BAND;                              /* POLDER data band */
	int		DIST;                                  /* tree distribution */
	double	E_r;                         /* mean distance between trees */
	int		Error[200];     /*                ERROR[1] = Px_tot problem */
	/*   (1) --> Neyman Px_tot 
	     (2) --> PG 
		 (3) --> PT
		 (4) --> FF0 
		 .... see README.htnl
    */
	double  BACKGROUND_REF[MAX_WAVE];
	double  FOLIAGE_REF[MAX_WAVE];
	double  FOLIAGE_TRANS[MAX_WAVE];
	double  Wave[MAX_WAVE];

	int Optic_index;
	double	Fs;                                           /* eq. 67 [1] */
	double	Ft;                                           /* eq. 52 [1] */
	double  Fn;               /* nadir maximum shadowed ground view ... */
	double  Fd;
	double	Fo;                                           /* eq. 31 [2] */
	double	Gv,Gs;                                         /* eq. 3 [2] */
	double  GFoliage;
	double	H;                       /* eq. 68 [1] and effective height */
	double	Hc;                                          /* cone height */
	double	Ls;                                          /* eq. 45 [1]  */
	double	Lb;                                           /* eq.  8 [2] */
	double	Lo;          /* LAI accumulated in sun direction */
	double	Lo_90;          /* LAI accumulated horizontaly in one crown */
	double	Lt;                                          /* eq. 42 [1]  */
	double	lambda_m;                  /* minimum gap size (eq. 52 [1]) */
	double  mu;
	double	OmegaT;                    /* clumping of trees, eq. 43 [1] */
	double  OmegaTotal;
	double	OmegaE_stand;                       /* Input clumping index */
	double	Ptrees[NN],Ptreev[NN];                        /* eq. 64 [1] */
	double	PT_Cold;                               /* Ptf in eq. 66 [1] */  
	double	PS;       /* prob. of seeing sunlit ground far from hotspot */  
	double	Pgap0;                       /* gap probability in a crown  */
	double	Pvg;                                        /* gap fraction */
	double  Pvg_mean;      /* mean gap based on LAI for view factors .. */
	double	Pig;                          /* prob. having sunlit ground */
	double  Pig_poisson; /* Same as Pig, but using Poisson distribution */
	double	Pti;                                             /*  eq. 39 */
	double	PVG;                                        /* gap fraction */
	double  PVG_NADIR;
	double	Ps,Pv,Psc,Pvc;    /* overlap probabilities eq 26. 27 .. [1] */ 
	double	PIG;                          /* prob. having sunlit ground */
	double	PT;           /*  prob./proportion of seeing sunlit foliage */
	double	PG;            /*  prob./proportion of seeing sunlit ground */
	double	PSG0_SUN;                  /* "sun" gap fraction no overlap */ 
	double	PSG0_VIEW;                       /* gap fraction no overlap */
	double	Px[NN];                  /* tree spatial distribution prob. */
	double	PgapV,PgapS;     /* gap inside crown, view or sun direction */
	double	PgapV_mean;     /* gap based on mean angle for view factors */
	double	phi;                                 /* azimuth diff. angle */
	double	PSG_HOT0;                        /*  gap fraction  at nadir */ 
	double	QQ1,QQ2,QQ1B,QQ2B;                    /* equ. 58 and 62 [1] */
	double  ro[MAX_WAVE];                 /* hyperspectral reflectance  */
	double	s;                                             /* mean path */  
	double	Sg,Sg_0,Sgc;         /* shadow from one crown on the ground */
	double	Ss,Sv;                                   /* see equ. 25 [1] */
	double	SZA_TMP;  /* used as a flag for not recomputing some values */
	double  SR;			                   /* Simple Ratio NIR / RED    */
	double	tic,tib,tac,tab;                   /* see eq. 15 and 39 [1] */
	double	V;                                          /* crown volume */
	double	vza;                                   /* view zenith angle */
	double	Vg,Vg_0,Vgc;             /* ground projected crown element  */
	double  Vg_0_mean; /* Vg_0 at mean angle based on LAI for view fact.*/
	double  Viewed_shadow;
	double	Wt;                                        /*   eq. 41 [1]  */  
	double	xi;                                          /* phase angle */
	double	ZT;      /* proportion/probability of seeing shaded foliage */
	double	ZG;      /* proportion/probability of seeing shaded ground  */
	
	double  TZ1,TZ2,TZ3,TZ4,GZ1,GZ2,GZ3,GZ4;

	} ;





