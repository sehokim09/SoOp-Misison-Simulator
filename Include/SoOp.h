
#ifndef __SoOp_H
#define __SoOp_H__
/*
** #ifdef __cplusplus
** namespace _42 {
** using namespace Kit;
** #endif
*/

#define Clight 299792458 // speed of light [m/s]
#define Boltzmann -228.6 // Boltzmann's constant [dBW/K/Hz]

#define MAXNUM_TX 8
#define MAXNUM_RX 6
#define ANTBODY 0 // Antenna body number

// antenna type
#define ANT_DIPOLE 0
#define ANT_PATCH 1

// polarizations
#define POL_H 0
#define POL_V 1
#define POL_X 2
#define POL_Y 3
#define POL_R 4
#define POL_L 5

// Tx Type
#define TX_MUOS 0
#define TX_ORBCOMM 1
#define TX_NOAA 2
#define TX_METEOR 3
#define TX_GPS 2
#define TX_GLONASS 3
#define TX_GALILEO 4
#define TX_BEIDOU 5

// Receiver
#define SC_RX1 0
#define SC_RX6 5

// MUOS
#define SC_MUOS1 16
#define SC_MUOS4 19
#define SC_MUOS5 9

// ORBCOMM
#define SC_FM4 20
#define SC_FM118 55

// NOAA
#define SC_NOAA1 42
#define SC_NOAA20 61

// METEOR
#define SC_METEOR_M1 62
#define SC_METEOR_M2_2 64

// GPS
#define SC_PRN1 56
#define SC_PRN32 86

// Glonass
#define SC_COSMOS2426 87
#define SC_COSMOS2534 111

// Galileo
#define SC_E1 112
#define SC_E36 137

// Beidou
#define SC_C1 138
#define SC_M22 185

#include <complex.h>
#include <stdint.h>

struct SurfaceDynType{
	char filename_snowCover[12][30]; // For monthly files
	char filename_soilTemp1[30];
	char filename_RZSM[30];
	uint8_t snowCover[3600][7200]; // 0.05 deg
	float soilTemp1[1624][3856]; // Soil Temperature Layer 1 [K] in 9km EASE grid2
	float RZSM[1624][3856]; // Root Zone Soil Moisture [m^3/m^3] in 9km EASE grid2
};

struct SurfaceSttType{
	char filename_DEM[40];
	char filename_landMask[30];
	char filename_landType[30];
	char filename_sand[30];
	char filename_clay[30];
	char filename_rho_b[30];
	int16_t height[10800][21600]; // ellipsoid height [m]
	uint8_t landMask[14616][34704]; // 1km EASE grid2
	uint8_t landType[3600][7200]; // 0.05 deg
	float sand[14616][34704]; // Mass Fractions of Sand (0<=S<=1) in 1km EASE grid2
	float clay[14616][34704]; // Mass Fractions of Clay (0<=C<=1) in 1km EASE grid2
	float rho_b[14616][34704]; // bulk density of soil sample in [g/cc] (grams/cubic centimeter)
};

struct SoilType{
	float m_v; // volumetric water content in g*cm^-3
	float sand; // sand textural component of a soil by weight (0<=S<=1)
	float clay; // clay textural component of a soil by weight (0<=C<=1)
	float T; // temperature in [C]
	double alpha; // empirically determined constant - from Peplinski 1995 p.804
	float rho_b; // bulk density of soil sample in [g/cc] (grams/cubic centimeter)
	double rho_s; // specific density of solid soil particles in [g/cm^3]
	double rho_w; // specific density of water in [g/cm^3]
	double rho_i; // specific density of ice in [g/cm^3]
	double e_s; // relative permittivity of solid soil, Equation 22 in Dobson 1985
	double e_i; // relative permittivity of inclusions (air, bound water, and free water) from L.Zhang 2003
	double complex e_c; // Complex dielectric constant
	double e_p; // Real part of dielectric constant
	double e_dp; // Imaginary part of dielectric constant
	double depthP; // penetration depth [m]
	double depthS; // skin depth [m]
};

struct ReflectType{
	double complex FVV; // Standard Fresnel coefficient
	double complex FHH;
	double complex FLR;
	double complex FRR;
	double gammaV; // reflectivity (Vertical polarization)
	double gammaH; // reflectivity (Horizontal polarization)
	double gammaLR; // reflectivity ( polarization)
	double gammaRR; // reflectivity ( polarization)
	double RcVV; // reflection coefficient
	double RcHH;
	double RcLR;
	double RcRR;
	double RxPwrR; // Received Signal Power from Reflected Signal [dBW]
	double SnrR; // SNR of Reflected Signal [dBW]
};

struct LocalSurfaceType{
	double theta;
	double phi;
	double range;
	double delay;
	double doppler;
	double incAngL2Tx;
	double incAngL2Rx;
	double scaAng;
	double scaVec[3];
	double lat; // Geodetic latitude [rad]
	double lon;	// Geodetic longitude [rad]
	double alt; // Altitude normal to ellipsoid [m]
	double dArea;
};
struct SpecularType{
	long Exists;
	long SpNum; // Label
	double RvB[3]; // Reflected signal vector in body frame
	double DvB[3]; // Direct signal vector in body frame
	double CEW[3][3]; // Direction Cosine Matrix from ECEF to ENU
	double CSE[3][3]; // Direction Cosine Matrix from ENU to Specular frame
	double CLS[3][3]; // Direction Cosine Matrix from Specular to local surface frame
	double CSW[3][3]; // Direction Cosine Matrix from ECEF to Specular frame
	double CLW[3][3]; // Direction Cosine Matrix from ECEF to local surface frame
	double PosN[3]; // Position of Specular Point expressed in ECI frame [m]
	double PosW[3]; // Position of Specular Point expressed in ECEF frame [m]
	double latGRS; // latitude of specular point expressed in GRS80 [rad]
	double lonGRS; // longitude of specular point expressed in GRS80 [rad]
	double altGRS; // altitude of specular point normal to GRS80 ellipsoid [m]
	double latSUR;
	double lonSUR;
	double altSUR;
	double PathLenD; // Direct Path Length [m]
	double PathLenGRS; // Specular Reflection Path Length on GRS80 [m]
	double PathLenSUR; // Specular Reflection Path Length on EARTH2014 Surface [m]
	double az; // Azimuth angle [rad]
	double elev; // Elevation angle [rad]
	double incAng; // Incidence angle [rad]
//	double slantRng; // Slant Range [m]
	double lookAngAzR; // Earth-view Antenna Azimuth Look Angle [rad]
	double lookAngElR; // Earth-view Antenna Elevation Look Angle [rad]
	double lookAngAzD; // Sky-view Antenna Azimuth Look Angle [rad]
	double lookAngElD; // Sky-view Antenna Elevation Look Angle [rad]
	double AntGainS; // Antenna Gain of Sky-view Antenna [dB]
	double AntGainE; // Antenna Gain of Earth-view Antenna [dB]
	uint8_t landMask; // Land mask (0: water, 1: land)
	uint8_t snowCover; // Snow mask
	uint8_t visibleTarget; // Visibility flag for Target
	double a; // 1st Fresnel Zone Semi-major Axis [m]
	double b; // 1st Fresnel Zone Semi-minor Axis [m]
	double area; // Area of 1st Fresnel zone area = pi*a*b [m^2]
	double K; // Gain of gradient function
	long cnt; // the number of iteration to find specular point
	double correction;
	double angErr;
	double gradient;
	double incAng2;
	double angBisNor;
	double RxPwrD; // Received Signal Power from Direct Signal [dBW]
	double SnrD; // SNR of Direct Signal [dBW]
	double EffReflCoef; // Effective Reflection Coefficient [dB]
	double VarSM; // Variance of Soil Moisture Estimate [dB]
	double VarERC; // Variance of Effective Reflection Coefficient [dB]
	double r;

	struct SoilType *Soil;
	struct ReflectType *Reflect;
	struct LocalSurfaceType **Local;
};

struct SoOpType{
	long NSp; // The Maximum Number of Specular Points
	long SpCnt; // The Number of Specular points
	long NValidSp; // Number of Valid Specular Points
	struct SpecularType *Sp;
	char SpriteFileName[40];
	char AntFileName[40]; // Antenna Pattern File Name
	long Nphi; // Number of rows (phi-angle Look-up) in Antenna Pattern File
	long Ntheta; // Number of columns (theta-angle Look-up) in Antenna Pattern File
    unsigned int SpriteTexTag;
	double Freq; // [Hz]
	double Wavelength; // [m]
	double Bandwidth; // [Hz]
	double CorIntTime; // Coherent Integration Time [sec]
	long AntTag; // Antenna Type - USER-DEFINED, DIPOLE, PATCH etc.
	double **AntPattern; // Antenna Pattern [dBi]
	double AntMaxGainS; // Antenna Max Gain of Sky-view Antenna [dB]
	double AntMaxGainE; // Antenna Max Gain of Earth-view Antenna [dB]
	double NoiseTemS; // Noise Temperature of Sky-view Antenna [K]
	double NoiseTemE; // Noise Temperature of Earth-view Antenna [K]
	double NoisePwrS; // Noise Power arriving at Sky-view Antenna [dBW]
	double NoisePwrE; // Noise Power arriving at Earth-view Antenna [dBW]
	double TotalRxPwrD; // Total Power of Received Direct Signals [dB]
	double TotalRxPwrR; // Total Power of Received Reflected Signals [dB]
	/* Local Surface Grid Input Parameters */
	double minSurfRes;
	long NthetaLocal;
	long NphiLocal;
};

struct RxType{
	long SC; // Assigned Spacecraft
    long NSoOp;   /* Number of SoOp Configurations */
	struct SoOpType *SoOp;
	double SwathWidth; // [m]
	double SwathArea; // [m^2]
	long Pol; // polarization
};

struct TxType{
	double Freq; // [Hz]
	double Wavelength; // [m]
	double EIRP; // [dBW]
	double Bandwidth; // [Hz]
	double NChannel; // the number of channels
	long Pol; // polarization
};

struct GsType{
	long GroundStation; // Assigned Ground Station
    long NSoOp;   /* Number of SoOp Configurations */
	struct SoOpType *SoOp;
	double SwathWidth; // [m]
	double SwathArea; // [m^2]
	double CorIntTime; // Coherent Integration Time [sec]
};

struct SpotBeamType{
	// Variables for global coverage
	double x[16];
	double y[16];
	double r0sq; // Reference radius square
	double r[17]; // actual radius of a spotbeam
	double el[4]; // Beam spacing elevation - origin + 3 layers of spot beams
	double az[3][5]; // Beam spacing azimuth - each layer has 5 azimuth angles
	double scaleFactor;
	long visibleSb[16];
	double latGRS[5][17];
	double lonGRS[5][17];
	double altGRS[5][17];
	// Variables for spotbeam on Target
	double halfBW; // Half beamwidth of spotbeam
	double halfBW_whole;
	double alt; // altitude of MUOS
	double r0Target; // radius of spotbeam
};

struct TargetAreaType{
	double radius;
	double lat;
	double lon;
	double alt; // Ellipsoid height [m]
	double PosW[3];
	double PosN[3];
	double vTarget2SpW[3]; // Vector from Target to Specular point in ECEF frame
	double mTarget2Sp; // Distance between Target and Specular point
	double CWT[3][3]; // Direction Cosine Matrix from ECEF to Target frame
};

struct PostProcType{
	uint16_t arc, TxID, site, visible, landMask;
	float m_v, Temp, sand, clay, rho_b;
	double t, lat, lon, alt, a, b, az, gammaLR, SnrR, SnrD, lookAngAzD, lookAngElD, lookAngAzR, lookAngElR;
};

/* Functions for Signals of Opportunity */
void FindSP_MPL(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX);
long FindSP_MPL_Spotbeam(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX);
long FindSP_MPL_Spotbeam_Gs(struct GroundStationType *GS, struct SpecularType *Sp, struct SCType *TX);
void FindSP_UD(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX);
long FindSP_UD_Spotbeam_Gs(struct GroundStationType *GS, struct SpecularType *Sp, struct SCType *TX);

long *FindVisibleTx(struct SCType *RX, long ITxStart, long ITxEnd, long *TxNum, long *TxCnt);
long *FindVisibleTx_Gs(struct GroundStationType *Gs, long ITxStart, long ITxEnd, long *TxNum, long *TxCnt);

uint8_t FindVisibleSpotbeam(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX);
long FindVisibleSpotbeam_Gs(struct GroundStationType *GS, struct SpecularType *Sp, struct SCType *TX);

void calcAntPattern(struct SoOpType *SoOp, struct SpecularType *Sp);
void CalcRxPwr(long RxNum, long SoOpNum, long SpNum, long TxNum, double Reflectivity);
void calcReflectivity(double complex e_c, struct ReflectType *R, double incAng);
void calcDielPeplinski(uint16_t rowSt, uint16_t colSt, uint16_t rowDy, uint16_t colDy, double Freq, struct SoilType *Soil);
void calcDielMironov(uint16_t rowSt, uint16_t colSt, uint16_t rowDy, uint16_t colDy, double f_Hz, struct SoilType *Soil);
//long *CntValidSpecularPt(long RxNum, long SoOpNum, long RoiNum);
//void CalcValidRxPwr(long RoiNum, double Reflectivity);
void FindSwathArea(long RxNum);
long FindClosestTx(long RxNum, long ITxStart, long ITxEnd);
long *FindBistaticGeo(long RxNum, long SoOpNum, long TxStart, long TxEnd);
void FindPathLength(long RxNum, long SoOpNum, long SpNum, long TxNum);

void LocalGeometry(double vSp2RxW[3], struct SCType *RX, struct SoOpType *SoOp, struct SpecularType *Sp, struct SCType *TX);
void ObsOnTarget(void);
void Obs_Calval_Snow_CONUS(void);
void Obs_CalVal_MUOS(void);
void ObsGlobeMUOS(void);
void ObsGlobeMOG(void);
void ObsCONUS_MUOS(void);

void ReportAttitude(void);
void ReportLLA_Rx(void);
void ReportPosN(void);
void ReportPosW(void);
void ReportOrbElements(void);

void TrackingSpotbeam(void);

void CmdReport(void);

void Reflectometry_MUOS(void);
void RepeatGroundTrack(struct OrbitType *O);

void UpdateOrbitCycle(void);
void UpdateSnowCover(void);
void LoadSnowCover(void);
void signalVec(long TxNum);
void LL2EzRC1km(double lat, double lon, uint16_t *row, uint16_t *col);
void LL2EzRC9km(double lat, double lon, uint16_t *row, uint16_t *col);

void InitFOVs(void);
void InitSpotbeams(void);
void InitSoOp(void);

void postProcessCalVal(uint8_t RxNum, uint8_t TxNum);

void RxFsw(struct SCType *S);

void GyroProcessing(struct AcType *AC);
void MagnetometerProcessing(struct AcType *AC);
void CssProcessing(struct AcType *AC);
void FssProcessing(struct AcType *AC);
void StarTrackerProcessing(struct AcType *AC);
void GpsProcessing(struct AcType *AC);
void AccelProcessing(struct AcType *AC);
void WheelProcessing(struct AcType *AC);
void MtbProcessing(struct AcType *AC);
void ThreeAxisAttitudeCommand(struct SCType *S);
/////////////////////////////////////////////
/*
** #ifdef __cplusplus
** }
** #endif
*/
#endif /* __SoOp_H__ */
