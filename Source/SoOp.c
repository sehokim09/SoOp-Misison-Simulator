#include "42.h"

#define EXTERN
#include "42GlutGui.h"
#undef EXTERN

/* Return X to the Y power.  */
double newpow(double x, double y){
    int sign=1;
    double r;
    if (x<0){
        sign = (-1);
        x *= sign;
    }
    r = pow(x, y);
    return (r*sign);
}

double Limit_max(double x, double max)
{
      return(x > max ? max : x);
}
void ECEFToWGS84_LAT(double p[3], double *glat)
{
      double a = 6378137.0;
      double f = 1.0/298.257223563;
      double b = a*(1.0-f);
      double e2 = f*(2.0-f);
      double ep2 = f*(2.0-f)/(1.0-f)/(1.0-f);
      double r,E2,F,G,C,S,P,Q,r0,V,Z0;

      double OneMinusE2,Z1,SpolyG,Qpoly;

      OneMinusE2 = 1.0-e2;

      r = sqrt(p[0]*p[0]+p[1]*p[1]);

      E2 = a*a-b*b;

      Z1 = b*p[2];

      F = 54.0*Z1*Z1;

      G = r*r+OneMinusE2*p[2]*p[2]-e2*E2;

      Z1 = e2*r/G;
      C = Z1*Z1*F/G;

      S = pow(1.0+C+sqrt(C*C+2.0*C),1.0/3.0);

      SpolyG = (S+1.0/S+1.0)*G;
      P = F/(3.0*SpolyG*SpolyG);

      Q = sqrt(1.0+2.0*e2*e2*P);

      Qpoly = 1.0+Q;
      r0 = -P*e2*r/Qpoly+sqrt(0.5*a*a*Qpoly/Q-P*OneMinusE2*p[2]*p[2]/(Q*Qpoly)-0.5*P*r*r);

      Z1 = r-e2*r0;
      Z1 *= Z1;

      V = sqrt(Z1+OneMinusE2*p[2]*p[2]);

      Z1 = b*b/a/V;
      Z0 = Z1*p[2];

      *glat = atan((p[2]+ep2*Z0)/r);
}
/**********************************************************************/
/* http://en.wikipedia.org/wiki/Geodetic_system#Geodetic_versus_geocentric_latitude */
void GRS80ToECEF(double glat, double glong, double alt, double p[3])
{
      double a = 6378137.0;
      double f = 1.0/298.257222101; // GRS80
//      double f = 1.0/298.257223563; // WGS84
      double e2 = f*(2.0-f);
      double X;
      double CosLat,SinLat,CosLng,SinLng;

      CosLat = cos(glat);
      SinLat = sin(glat);
      CosLng = cos(glong);
      SinLng = sin(glong);

      X = sqrt(1.0-e2*SinLat*SinLat);

      p[0] = (a/X+alt)*CosLat*CosLng;
      p[1] = (a/X+alt)*CosLat*SinLng;
      p[2] = (a/X*(1.0-e2)+alt)*SinLat;
}
/**********************************************************************/
/* J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates
   to geodetic coordinates," IEEE Transactions on Aerospace and
   Electronic Systems, vol. 30, pp. 957-961, 1994. */
void ECEFToGRS80(double p[3], double *glat, double *glong, double *alt)
{
      double a = 6378137.0;
      double f = 1.0/298.257222101; // GRS80
//      double f = 1.0/298.257223563; // WGS84
      double b = a*(1.0-f);
      double e2 = f*(2.0-f);
      double ep2 = f*(2.0-f)/(1.0-f)/(1.0-f);
      double r,E2,F,G,C,S,P,Q,r0,U,V,Z0;

      double OneMinusE2,Z1,SpolyG,Qpoly;

      OneMinusE2 = 1.0-e2;

      r = sqrt(p[0]*p[0]+p[1]*p[1]);

      E2 = a*a-b*b;

      Z1 = b*p[2];

      F = 54.0*Z1*Z1;

      G = r*r+OneMinusE2*p[2]*p[2]-e2*E2;

      Z1 = e2*r/G;
      C = Z1*Z1*F/G;

      S = pow(1.0+C+sqrt(C*C+2.0*C),1.0/3.0);

      SpolyG = (S+1.0/S+1.0)*G;
      P = F/(3.0*SpolyG*SpolyG);

      Q = sqrt(1.0+2.0*e2*e2*P);

      Qpoly = 1.0+Q;
      r0 = -P*e2*r/Qpoly+sqrt(0.5*a*a*Qpoly/Q-P*OneMinusE2*p[2]*p[2]/(Q*Qpoly)-0.5*P*r*r);

      Z1 = r-e2*r0;
      Z1 *= Z1;

      U = sqrt(Z1+p[2]*p[2]);

      V = sqrt(Z1+OneMinusE2*p[2]*p[2]);


      Z1 = b*b/a/V;
      Z0 = Z1*p[2];

      *alt = U*(1.0-Z1);
      *glat = atan((p[2]+ep2*Z0)/r);
      *glong = atan2(p[1],p[0]);
}

void ECEFToGRS80_DEM(double p[3], double *glat, double *glong, double *alt)
{
      double a = 6378137.0;
      double f = 1.0/298.257222101; // GRS80
      double b = a*(1.0-f);
      double e2 = f*(2.0-f);
      double ep2 = f*(2.0-f)/(1.0-f)/(1.0-f);
      double r,E2,F,G,C,S,P,Q,r0,U,V,Z0;

      double OneMinusE2,Z1,SpolyG,Qpoly;

      int row, col;

      OneMinusE2 = 1.0-e2;

      r = sqrt(p[0]*p[0]+p[1]*p[1]);

      E2 = a*a-b*b;

      Z1 = b*p[2];

      F = 54.0*Z1*Z1;

      G = r*r+OneMinusE2*p[2]*p[2]-e2*E2;

      Z1 = e2*r/G;
      C = Z1*Z1*F/G;

      S = pow(1.0+C+sqrt(C*C+2.0*C),1.0/3.0);

      SpolyG = (S+1.0/S+1.0)*G;
      P = F/(3.0*SpolyG*SpolyG);

      Q = sqrt(1.0+2.0*e2*e2*P);

      Qpoly = 1.0+Q;
      r0 = -P*e2*r/Qpoly+sqrt(0.5*a*a*Qpoly/Q-P*OneMinusE2*p[2]*p[2]/(Q*Qpoly)-0.5*P*r*r);

      Z1 = r-e2*r0;
      Z1 *= Z1;

      U = sqrt(Z1+p[2]*p[2]);

      V = sqrt(Z1+OneMinusE2*p[2]*p[2]);


      Z1 = b*b/a/V;
      Z0 = Z1*p[2];

      *glat = atan((p[2]+ep2*Z0)/r);
      *glong = atan2(p[1],p[0]);

      /* Integrating DEM */
      row = *glat*60+5400;
      col = *glong*60+10800;

      *alt = U*(1.0-Z1) + SurfaceStt->height[row][col] - 7137;
      printf("%lf, %lf\n", U*(1.0-Z1), *alt);
}

void FindSP_MPL(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX)
{
	double estCurr[3], corrEst[3], estNext[3], estGradi[3], corr[3];
	double PosWR[3], PosWT[3], ST[3], SR[3], SurfaceNorVec[3], BistaticVec[3];
	double MAG1, MAG2, Lat;
//	double temp[3];

//	double vecRx2SpN[3], vecNadir[3];
	// WGS84
    double a = 6378137.0;
    double f = 1.0/298.257223563;
    double b = a*(1.0-f);
    double r;
    // Tolerance
	double sigma = 0.1;
	double K;
	double correction = 20;
	long cnt = 0;

	K = Sp->K;
	//K = 400000.0;
	PosWR[0] = RX->PosRxAntW[0];
	PosWR[1] = RX->PosRxAntW[1];
	PosWR[2] = RX->PosRxAntW[2];

	PosWT[0] = TX->PosW[0];
	PosWT[1] = TX->PosW[1];
	PosWT[2] = TX->PosW[2];

	MAG1 = PosWR[0]*PosWR[0] + PosWR[1]*PosWR[1] + PosWR[2]*PosWR[2];
	MAG2 = PosWT[0]*PosWT[0] + PosWT[1]*PosWT[1] + PosWT[2]*PosWT[2];

	/* Unconstrained Weighted Initial Estimate */
	estCurr[0] = PosWR[0]/MAG1 + PosWT[0]/MAG2;
	estCurr[1] = PosWR[1]/MAG1 + PosWT[1]/MAG2;
	estCurr[2] = PosWR[2]/MAG1 + PosWT[2]/MAG2;

	/* Initial specular point guess, on the surface directly below R */
	MAG1 = sqrt(estCurr[0]*estCurr[0] + estCurr[1]*estCurr[1] + estCurr[2]*estCurr[2]);
    Lat = asin(estCurr[2]/MAG1); // Geocentric latitude
    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
//    r = a*sqrt((1-e*e)/(1-(e*e*cos(Lat)*cos(Lat))));

	estCurr[0] = r/MAG1*estCurr[0];
	estCurr[1] = r/MAG1*estCurr[1];
	estCurr[2] = r/MAG1*estCurr[2];

	while(correction > sigma && cnt<20000) {
		// Gain adaptation
		//if(correction<10) K=100000;

	    /* Step 1: Calculate gradient of path length function */
		ST[0] = estCurr[0]-PosWT[0];
		ST[1] = estCurr[1]-PosWT[1];
		ST[2] = estCurr[2]-PosWT[2];

		SR[0] = estCurr[0]-PosWR[0];
		SR[1] = estCurr[1]-PosWR[1];
		SR[2] = estCurr[2]-PosWR[2];

		MAG1 = sqrt(ST[0]*ST[0] + ST[1]*ST[1] + ST[2]*ST[2]);
		MAG2 = sqrt(SR[0]*SR[0] + SR[1]*SR[1] + SR[2]*SR[2]);

		estGradi[0] = ST[0]/MAG1 + SR[0]/MAG2;
		estGradi[1] = ST[1]/MAG1 + SR[1]/MAG2;
		estGradi[2] = ST[2]/MAG1 + SR[2]/MAG2;

		BistaticVec[0] = -estGradi[0];
		BistaticVec[1] = -estGradi[1];
		BistaticVec[2] = -estGradi[2];

		MAG1 = sqrt(BistaticVec[0]*BistaticVec[0] + BistaticVec[1]*BistaticVec[1] + BistaticVec[2]*BistaticVec[2]);
		BistaticVec[0] /= MAG1;
		BistaticVec[1] /= MAG1;
		BistaticVec[2] /= MAG1;

		SurfaceNorVec[0] = estCurr[0]*2.0/a/a;
		SurfaceNorVec[1] = estCurr[1]*2.0/a/a;
		SurfaceNorVec[2] = estCurr[2]*2.0/b/b;

	    /* Step 2 */
		// Tangential correction optimization
		//VxV(SurfaceNorVec,estGradi,temp);
		//VxV(temp, SurfaceNorVec, estGradi);
		//////////////////////////////////////////
		MAG2 = sqrt(SurfaceNorVec[0]*SurfaceNorVec[0] + SurfaceNorVec[1]*SurfaceNorVec[1] + SurfaceNorVec[2]*SurfaceNorVec[2]);
		SurfaceNorVec[0] /= MAG2;
		SurfaceNorVec[1] /= MAG2;
		SurfaceNorVec[2] /= MAG2;

		corrEst[0] = estCurr[0] - K*estGradi[0];
		corrEst[1] = estCurr[1] - K*estGradi[1];
		corrEst[2] = estCurr[2] - K*estGradi[2];

		/* Step 3 */
		MAG1 = sqrt(corrEst[0]*corrEst[0] + corrEst[1]*corrEst[1] + corrEst[2]*corrEst[2]);
	    Lat = asin(corrEst[2]/MAG1);
	    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));

		estNext[0] = r/MAG1*corrEst[0];
		estNext[1] = r/MAG1*corrEst[1];
		estNext[2] = r/MAG1*corrEst[2];

	    /* Step 4 */
		corr[0] = estNext[0]-estCurr[0];
		corr[1] = estNext[1]-estCurr[1];
		corr[2] = estNext[2]-estCurr[2];

		correction = sqrt(corr[0]*corr[0] + corr[1]*corr[1] + corr[2]*corr[2]);
//		printf("c:%le %le %le\n", estCurr[0], estCurr[1], estCurr[2]);
//		printf("n:%le %le %le\n", estNext[0], estNext[1], estNext[2]);
//		printf("n:%le %le %le\n", corr[0], corr[1], corr[2]);
		estCurr[0] = estNext[0];
		estCurr[1] = estNext[1];
		estCurr[2] = estNext[2];
		cnt++;
	}
	Sp->gradient = MAGV(estGradi);
	double RS[3], TS[3];
	//CopyUnitV(estCurr,uSpPosW);
	ST[0] = PosWT[0]-estCurr[0];
	ST[1] = PosWT[1]-estCurr[1];
	ST[2] = PosWT[2]-estCurr[2];

	SR[0] = PosWR[0]-estCurr[0];
	SR[1] = PosWR[1]-estCurr[1];
	SR[2] = PosWR[2]-estCurr[2];

	MAG1 = sqrt(ST[0]*ST[0] + ST[1]*ST[1] + ST[2]*ST[2]);
	MAG2 = sqrt(SR[0]*SR[0] + SR[1]*SR[1] + SR[2]*SR[2]);

	BistaticVec[0] = ST[0]/MAG1 + SR[0]/MAG2;
	BistaticVec[1] = ST[1]/MAG1 + SR[1]/MAG2;
	BistaticVec[2] = ST[2]/MAG1 + SR[2]/MAG2;
	UNITV(BistaticVec);

	SurfaceNorVec[0] = estCurr[0]*2.0/a/a;
	SurfaceNorVec[1] = estCurr[1]*2.0/a/a;
	SurfaceNorVec[2] = estCurr[2]*2.0/b/b;
	UNITV(SurfaceNorVec);

	RS[0] = PosWR[0] - estCurr[0];
	RS[1] = PosWR[1] - estCurr[1];
	RS[2] = PosWR[2] - estCurr[2];
	UNITV(RS);

	TS[0] = PosWT[0] - estCurr[0];
	TS[1] = PosWT[1] - estCurr[1];
	TS[2] = PosWT[2] - estCurr[2];
	UNITV(TS);

	Sp->incAng2 = acos(VoV(SurfaceNorVec, TS));
	Sp->angErr = Sp->incAng2 - acos(VoV(SurfaceNorVec, RS));
	Sp->angBisNor = acos(VoV(BistaticVec, SurfaceNorVec));
	Sp->correction = correction;
	Sp->cnt = cnt;
	//printf("%ld\n", cnt);
	Sp->PosW[0] = estCurr[0];
	Sp->PosW[1] = estCurr[1];
	Sp->PosW[2] = estCurr[2];

	/* Check visibility of specular point from receiver*/
//	vecRx2SpN[0] = Sp->PosN[0] - RX->PosN[0];
//	vecRx2SpN[1] = Sp->PosN[1] - RX->PosN[1];
//	vecRx2SpN[2] = Sp->PosN[2] - RX->PosN[2];

//	UNITV(vecRx2SpN);
//	vecNadir[0] = -RX->PosN[0];
//	vecNadir[1] = -RX->PosN[1];
//	vecNadir[2] = -RX->PosN[2];
//	UNITV(vecNadir);

//	Sp->lookAng = acos(VoV(vecRx2SpN,vecNadir));
	//printf("%lf\n", Sp->lookAng*R2D);

	//if(Sp->lookAng<FOV[RX->ID].Width/2)
	//	Sp->Exists = 1;
}

void FindSP_UD(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX)
{
	double estCurr[3], corrEst[3], estNext[3], estGradi[3], corr[3];
	double SurfaceNorVec[3], BistaticVec[3], PosWR[3], PosWT[3], TS[3], RS[3];
	double MAG1, MAG2, Lat;

	// Earth Ellipsoid
    double a = 6378137.0;
    double f = 1.0/298.257222101; // GRS80
//    double f = 1.0/298.257223563; // WGS84
    double b = a*(1.0-f);
    double r;
 //   double e = sqrt(2.0*f-f*f);
    // Tolerance
	double sigma = 0.1;
	double K;
	double correction = 20;
	long cnt = 0;

	K = Sp->K;

	PosWR[0] = RX->PosRxAntW[0];
	PosWR[1] = RX->PosRxAntW[1];
	PosWR[2] = RX->PosRxAntW[2];

	PosWT[0] = TX->PosW[0];
	PosWT[1] = TX->PosW[1];
	PosWT[2] = TX->PosW[2];

	MAG1 = PosWR[0]*PosWR[0] + PosWR[1]*PosWR[1] + PosWR[2]*PosWR[2];
	MAG2 = PosWT[0]*PosWT[0] + PosWT[1]*PosWT[1] + PosWT[2]*PosWT[2];

	/* Unconstrained Weighted Initial Estimate */
	estCurr[0] = PosWR[0]/MAG1 + PosWT[0]/MAG2;
	estCurr[1] = PosWR[1]/MAG1 + PosWT[1]/MAG2;
	estCurr[2] = PosWR[2]/MAG1 + PosWT[2]/MAG2;

	/* Constrain back to the surface */
	MAG1 = sqrt(estCurr[0]*estCurr[0] + estCurr[1]*estCurr[1] + estCurr[2]*estCurr[2]);
    Lat = asin(estCurr[2]/MAG1); // Geocentric latitude
    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
//    r = a*sqrt((1-e*e)/(1-(e*e*cos(Lat)*cos(Lat))));

	estCurr[0] = r/MAG1*estCurr[0];
	estCurr[1] = r/MAG1*estCurr[1];
	estCurr[2] = r/MAG1*estCurr[2];

	while(correction > sigma && cnt<20000) {
		// Gain adaptation
		if(correction<10){
		//	printf("Adjust gain (%lf)\n", correction);
			K=100000;
		}
		//if(cnt==19999) printf("Gain: %lf, cnt saturated (%lf)\n", K, correction);

	    /* Step 1: Calculate gradient of path length function */
		TS[0] = PosWT[0]-estCurr[0];
		TS[1] = PosWT[1]-estCurr[1];
		TS[2] = PosWT[2]-estCurr[2];

		RS[0] = PosWR[0]-estCurr[0];
		RS[1] = PosWR[1]-estCurr[1];
		RS[2] = PosWR[2]-estCurr[2];

		MAG1 = sqrt(TS[0]*TS[0] + TS[1]*TS[1] + TS[2]*TS[2]);
		MAG2 = sqrt(RS[0]*RS[0] + RS[1]*RS[1] + RS[2]*RS[2]);

		BistaticVec[0] = TS[0]/MAG1 + RS[0]/MAG2;
		BistaticVec[1] = TS[1]/MAG1 + RS[1]/MAG2;
		BistaticVec[2] = TS[2]/MAG1 + RS[2]/MAG2;

		MAG1 = sqrt(BistaticVec[0]*BistaticVec[0] + BistaticVec[1]*BistaticVec[1] + BistaticVec[2]*BistaticVec[2]);

		BistaticVec[0] /= MAG1;
		BistaticVec[1] /= MAG1;
		BistaticVec[2] /= MAG1;

		SurfaceNorVec[0] = estCurr[0]*2.0/a/a;
		SurfaceNorVec[1] = estCurr[1]*2.0/a/a;
		SurfaceNorVec[2] = estCurr[2]*2.0/b/b;

		MAG2 = sqrt(SurfaceNorVec[0]*SurfaceNorVec[0] + SurfaceNorVec[1]*SurfaceNorVec[1] + SurfaceNorVec[2]*SurfaceNorVec[2]);

		SurfaceNorVec[0] /= MAG2;
		SurfaceNorVec[1] /= MAG2;
		SurfaceNorVec[2] /= MAG2;

		estGradi[0] = BistaticVec[0] - SurfaceNorVec[0];
		estGradi[1] = BistaticVec[1] - SurfaceNorVec[1];
		estGradi[2] = BistaticVec[2] - SurfaceNorVec[2];

	    /* Step 2 */
		// Tangential correction optimization
		//VxV(SurfaceNorVec,estGradi,temp);
		//VxV(temp, SurfaceNorVec, estGradi);

		corrEst[0] = estCurr[0] + K*estGradi[0];
		corrEst[1] = estCurr[1] + K*estGradi[1];
		corrEst[2] = estCurr[2] + K*estGradi[2];

		/* Step 3 */
		MAG1 = sqrt(corrEst[0]*corrEst[0] + corrEst[1]*corrEst[1] + corrEst[2]*corrEst[2]);
	    Lat = asin(corrEst[2]/MAG1); // Geocentric latitude
	    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
	    //r = a*sqrt((1-e*e)/(1-(e*e*cos(Lat)*cos(Lat))));

		estNext[0] = r/MAG1*corrEst[0];
		estNext[1] = r/MAG1*corrEst[1];
		estNext[2] = r/MAG1*corrEst[2];

	    /* Step 4 */
		corr[0] = estNext[0]-estCurr[0];
		corr[1] = estNext[1]-estCurr[1];
		corr[2] = estNext[2]-estCurr[2];

		correction = sqrt(corr[0]*corr[0] + corr[1]*corr[1] + corr[2]*corr[2]);

		estCurr[0] = estNext[0];
		estCurr[1] = estNext[1];
		estCurr[2] = estNext[2];
		cnt++;
	}
	/*
	Sp->gradient = MAGV(estGradi);
//	double uSpPosW[3];
//	CopyUnitV(estCurr,uSpPosW);
	TS[0] = PosWT[0]-estCurr[0];
	TS[1] = PosWT[1]-estCurr[1];
	TS[2] = PosWT[2]-estCurr[2];
	UNITV(TS);

	RS[0] = PosWR[0]-estCurr[0];
	RS[1] = PosWR[1]-estCurr[1];
	RS[2] = PosWR[2]-estCurr[2];
	UNITV(RS);

	BistaticVec[0] = TS[0] + RS[0];
	BistaticVec[1] = TS[1] + RS[1];
	BistaticVec[2] = TS[2] + RS[2];
	UNITV(BistaticVec);

	SurfaceNorVec[0] = estCurr[0]*2.0/a/a;
	SurfaceNorVec[1] = estCurr[1]*2.0/a/a;
	SurfaceNorVec[2] = estCurr[2]*2.0/b/b;
	UNITV(SurfaceNorVec);

	Sp->incAng2 = acos(VoV(SurfaceNorVec, TS));
	Sp->angErr = Sp->incAng2 - acos(VoV(SurfaceNorVec, RS));
	Sp->angBisNor = acos(VoV(BistaticVec, SurfaceNorVec));
	Sp->correction = correction;
	Sp->cnt = cnt;
	//printf("%ld\n", cnt);
	*/

	Sp->PosW[0] = estCurr[0];
	Sp->PosW[1] = estCurr[1];
	Sp->PosW[2] = estCurr[2];

	//printf("%le\n",MAGV(estGradi));
	//printf("%lf\n", (tan(Sp->lat)-b/a*tan(Lat)));
}

long FindSP_UD_Spotbeam_Gs(struct GroundStationType *GS, struct SpecularType *Sp, struct SCType *TX)
{
	double estCurr[3], corrEst[3], estNext[3], estGradi[3], corr[3];
	double SurfaceNorVec[3], BistaticVec[3], PosWR[3], PosWT[3], TS[3], RS[3];
	double MAGR, MAGT, Lat;

	double vecRx2SpN[3], vecNadir[3];
	long nSpotBeam;
	// WGS84
    double a = 6378137.0;
    double f = 1.0/298.257223563;
    double b = a*(1.0-f);
    double r;
    // Tolerance
	double sigma = 0.1;
	double K;
	double MAG = 20;
	long cnt = 0;

	K = Sp->K;

	PosWR[0] = GS->PosW[0];
	PosWR[1] = GS->PosW[1];
	PosWR[2] = GS->PosW[2];

	PosWT[0] = TX->PosW[0];
	PosWT[1] = TX->PosW[1];
	PosWT[2] = TX->PosW[2];

	MAGR = PosWR[0]*PosWR[0] + PosWR[1]*PosWR[1] + PosWR[2]*PosWR[2];
	MAGT = PosWT[0]*PosWT[0] + PosWT[1]*PosWT[1] + PosWT[2]*PosWT[2];

	/* Unconstrained Weighted Initial Estimate */
	estCurr[0] = PosWR[0]/MAGR + PosWT[0]/MAGT;
	estCurr[1] = PosWR[1]/MAGR + PosWT[1]/MAGT;
	estCurr[2] = PosWR[2]/MAGR + PosWT[2]/MAGT;

	while(MAG > sigma && cnt<1000000) {
	//	if(MAG>10) K=500;
	//	else K=10;
	    /* Step 1: Calculate gradient of path length function */
		TS[0] = PosWT[0]-estCurr[0];
		TS[1] = PosWT[1]-estCurr[1];
		TS[2] = PosWT[2]-estCurr[2];

		RS[0] = PosWR[0]-estCurr[0];
		RS[1] = PosWR[1]-estCurr[1];
		RS[2] = PosWR[2]-estCurr[2];

		MAGT = sqrt(TS[0]*TS[0] + TS[1]*TS[1] + TS[2]*TS[2]);
		MAGR = sqrt(RS[0]*RS[0] + RS[1]*RS[1] + RS[2]*RS[2]);

		TS[0] /= MAGT;
		TS[1] /= MAGT;
		TS[2] /= MAGT;

		RS[0] /= MAGR;
		RS[1] /= MAGR;
		RS[2] /= MAGR;

		BistaticVec[0] = TS[0] + RS[0];
		BistaticVec[1] = TS[1] + RS[1];
		BistaticVec[2] = TS[2] + RS[2];

		MAG = sqrt(BistaticVec[0]*BistaticVec[0] + BistaticVec[1]*BistaticVec[1] + BistaticVec[2]*BistaticVec[2]);
		BistaticVec[0] /= MAG;
		BistaticVec[1] /= MAG;
		BistaticVec[2] /= MAG;

		SurfaceNorVec[0] = estCurr[0]*2/a/a;
		SurfaceNorVec[1] = estCurr[1]*2/a/a;
		SurfaceNorVec[2] = estCurr[2]*2/b/b;

		MAG = sqrt(SurfaceNorVec[0]*SurfaceNorVec[0] + SurfaceNorVec[1]*SurfaceNorVec[1] + SurfaceNorVec[2]*SurfaceNorVec[2]);
		SurfaceNorVec[0] /= MAG;
		SurfaceNorVec[1] /= MAG;
		SurfaceNorVec[2] /= MAG;

		estGradi[0] = BistaticVec[0] - SurfaceNorVec[0];
		estGradi[1] = BistaticVec[1] - SurfaceNorVec[1];
		estGradi[2] = BistaticVec[2] - SurfaceNorVec[2];

	    /* Step 2 */
		// Tangential correction optimization
	/*	VxV(estCurrNor,estGradi,temp);
		VxV(temp, estCurrNor, temp2);
		SxV(K, temp2, temp);

		corrEst[0] = estCurr[0] - temp[0];
		corrEst[1] = estCurr[1] - temp[1];
		corrEst[2] = estCurr[2] - temp[2];
	*/
		corrEst[0] = estCurr[0] + K*estGradi[0];
		corrEst[1] = estCurr[1] + K*estGradi[1];
		corrEst[2] = estCurr[2] + K*estGradi[2];

		/* Step 3 */
		MAG = sqrt(corrEst[0]*corrEst[0] + corrEst[1]*corrEst[1] + corrEst[2]*corrEst[2]);
	    Lat = asin(corrEst[2]/MAG);
	    r = a*b/sqrt(a*a*sin(Lat)*sin(Lat)+b*b*cos(Lat)*cos(Lat));
		corrEst[0] /= MAG;
		corrEst[1] /= MAG;
		corrEst[2] /= MAG;

		estNext[0] = r*corrEst[0];
		estNext[1] = r*corrEst[1];
		estNext[2] = r*corrEst[2];

	    /* Step 4 */
		corr[0] = estNext[0]-estCurr[0];
		corr[1] = estNext[1]-estCurr[1];
		corr[2] = estNext[2]-estCurr[2];

		MAG = sqrt(corr[0]*corr[0] + corr[1]*corr[1] + corr[2]*corr[2]);
//		printf("c:%le %le %le\n", estCurr[0], estCurr[1], estCurr[2]);
//		printf("n:%le %le %le\n", estNext[0], estNext[1], estNext[2]);
//		printf("n:%le %le %le\n", corr[0], corr[1], corr[2]);
		estCurr[0] = estNext[0];
		estCurr[1] = estNext[1];
		estCurr[2] = estNext[2];
		cnt++;
	}
	Sp->cnt = cnt;

	//printf("%ld\n", cnt);
	Sp->PosW[0] = estCurr[0];
	Sp->PosW[1] = estCurr[1];
	Sp->PosW[2] = estCurr[2];

	/* ECEF to ECI */
	MTxV(World[EARTH].CWN,estCurr,Sp->PosN);

	/* Check visibility of specular point from receiver*/
	vecRx2SpN[0] = Sp->PosN[0] - GS->PosN[0];
	vecRx2SpN[1] = Sp->PosN[1] - GS->PosN[1];
	vecRx2SpN[2] = Sp->PosN[2] - GS->PosN[2];

	UNITV(vecRx2SpN);
	vecNadir[0] = -GS->PosN[0];
	vecNadir[1] = -GS->PosN[1];
	vecNadir[2] = -GS->PosN[2];
	UNITV(vecNadir);

	//Sp->lookAng = acos(VoV(vecRx2SpN,vecNadir));
	//printf("%lf\n", Sp->lookAng*R2D);

	//if(Sp->lookAng<GS->FoV/2)
//	{
		nSpotBeam = FindVisibleSpotbeam_Gs(GS, Sp, TX);
		if(nSpotBeam>0)
			Sp->Exists = 1;
//	}

	/* ECEF to GRS80 */
	ECEFToGRS80(estCurr, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);
//	printf("%lf %lf %lf\n", Sp->lat*R2D, Sp->lon*R2D, Sp->alt*R2D);
	return nSpotBeam;
}

/**********************************************************************/
/* Calculate swath area of FOV. Assume spherical Earth, symmetric FOV */
void FindSwathArea(long RxNum)
{
	double alpha;
	double alpha1;
	double alpha2;
	double CBL[3][3];
	double vl[3] = {0.0, 0.0, 1.0};
	double vb[3];
	double theta, theta1, theta2, thetaMax;
	double halfFOV;
	double ratio;

	ratio = MAGV(SC[RxNum].PosN)/World[EARTH].rad;
	thetaMax = asin(1/ratio);
	halfFOV = FOV[RxNum].Width/2;

	/* DCM from Local to Body */
	MxMT(SC[RxNum].B[0].CN,SC[RxNum].CLN,CBL);
	/* +Z vector of Local frame is expressed in Body frame */
	MxV(CBL,vl,vb);
	/* Angle between +Z vector in Local and Body */
	theta = acos(VoV(vl,vb));

	/* Take into account FOV */
	if(theta > halfFOV){
		theta1 = theta - halfFOV;
		theta2 = Limit_max(theta + halfFOV,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 - alpha1;
	}
	else if(theta < halfFOV){
		theta1 = Limit_max(halfFOV - theta,thetaMax);
		theta2 = Limit_max(halfFOV + theta,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 + alpha1;
	}
	else{
		theta1 = theta;
		theta2 = Limit_max(halfFOV + theta,thetaMax);

		alpha1 = asin(sin(theta1)*ratio)-theta1;
		alpha2 = asin(sin(theta2)*ratio)-theta2;

		alpha = alpha2 + alpha1;
	}

	Rx[RxNum].SwathWidth = alpha/TwoPi*World[EARTH].rad;
	Rx[RxNum].SwathArea = TwoPi*World[EARTH].rad*World[EARTH].rad*(1-cos(alpha/2));
}

long FindClosestTx(long RxNum, long ITxStart, long ITxEnd)
{
	//double MaxToS = -1.0; /* Bogus */
	double ToS;
	double Rhat[3], RelPosN[3];
	long ITx, i;
	long TxNum = 0;
	double elevation;
	double MinElevation = -20*D2R;
	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(SC[RxNum].PosN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    for(ITx=ITxStart;ITx<ITxEnd+1;ITx++) {
       if (SC[ITx].Exists) {
          for(i=0;i<3;i++)
             RelPosN[i] = SC[ITx].PosN[i] - SC[RxNum].PosN[i];
          UNITV(RelPosN);
          ToS = VoV(RelPosN,Rhat);
          elevation = HalfPi - acos(ToS);
          if(elevation > MinElevation){
        	  MinElevation = elevation;
          //if (ToS > MaxToS) {
           //  MaxToS = ToS;
             TxNum = ITx;
          }
       }
    }
    //printf("ele:%f\n", MinElevation*R2D);
    return TxNum;
}

long *FindVisibleTxSoOp(struct SCType *RX, long *TxNum, long *TxCnt)
{
	double Rhat[3], RelPosN[3];
	long ITx, TxClosest[4];
	long cnt = 0, Flag = 0;
	double angLoS, angLoSmax = 1.221730476396031, angLoSmin = 1.221730476396031; // 70*D2R

	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(RX->PosRxAntN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    // MUOS
    for(ITx=SC_MUOS1;ITx<SC_MUOS4;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  TxClosest[cnt] = ITx;
        	  Flag = 1;
          }
    }
    if(Flag){
    	cnt++;
        Flag = 0;
    }

    // Orbcomm
    for(ITx=SC_FM4;ITx<SC_FM118;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  TxClosest[cnt] = ITx;
        	  Flag = 1;
          }
    }
    if(Flag){
    	cnt++;
        Flag = 0;
    }
/*
    // NOAA
    for(ITx=42;ITx<62;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  TxClosest[cnt] = ITx;
        	  Flag = 1;
          }
    }
    if(Flag){
    	cnt++;
        Flag = 0;
    }

    // METEOR
    for(ITx=62;ITx<65;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  TxClosest[cnt] = ITx;
        	  Flag = 1;
          }
    }
    if(Flag){
    	cnt++;
    }
*/
    *TxCnt = cnt;

    TxNum = (long *) calloc(cnt,sizeof(long));

	for(ITx=0; ITx<cnt; ITx++) {
		TxNum[ITx] = TxClosest[ITx];
	}
	return TxNum;
}


long *FindVisibleTxSG(struct SCType *RX, long *TxNum, long *TxCnt)
{
	static long bufferTx[30];
	static double bufferAngLoS[30];
	double Rhat[3], RelPosN[3];
	long ITx, TxClosest[20];
	long tempTx, cntTx = 0, cnt = 0, Flag = 0, i, j;
	double tempLoS, angLoS, angLoSmax = 1.221730476396031, angLoSmin = 1.221730476396031; // 70*D2R

	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(RX->PosRxAntN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    // MUOS
    for(ITx=SC_MUOS1;ITx<SC_MUOS4+1;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  TxClosest[cntTx] = ITx;
        	  Flag = 1;
          }
    }
    if(Flag){
    	cntTx++;
        Flag = 0;
    }

    // Orbcomm
    for(ITx=SC_FM4;ITx<SC_FM118+1;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  TxClosest[cntTx] = ITx;
        	  Flag = 1;
          }
    }
    if(Flag){
    	cntTx++;
        Flag = 0;
    }

    /* Aim at Tx closest to Rx's Zenith */
    // GPS
    for(ITx=SC_PRN1;ITx<SC_PRN32+1;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax) {
        	  bufferTx[cnt] = ITx;
        	  bufferAngLoS[cnt++] = angLoS;
          }
    }
    for(i=0; i<cnt; ++i) {
        for (j=i+1; j<cnt; ++j) {
            if (bufferAngLoS[i] > bufferAngLoS[j]) {
                tempLoS = bufferAngLoS[i];
                bufferAngLoS[i] = bufferAngLoS[j];
                bufferAngLoS[j] = tempLoS;

                tempTx = bufferTx[i];
                bufferTx[i] = bufferTx[j];
                bufferTx[j] = tempTx;
            }
        }
    }
    if(cnt>4){
    	for(i=0; i<4; i++)
    		TxClosest[cntTx++] = bufferTx[i];
    }
	else
		for(i=0; i<cnt; i++)
			TxClosest[cntTx++] = bufferTx[i];

    // Glonass
    cnt = 0;
    for(ITx=SC_COSMOS2426;ITx<SC_COSMOS2534+1;ITx++) {
		RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
		RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
		RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
		UNITV(RelPosN);
		angLoS = acos(VoV(RelPosN,Rhat));
		if(angLoS < angLoSmax) {
		  bufferTx[cnt] = ITx;
		  bufferAngLoS[cnt++] = angLoS;
		}
	}
	for(i=0; i<cnt; ++i) {
	  for (j=i+1; j<cnt; ++j) {
		  if (bufferAngLoS[i] > bufferAngLoS[j]) {
			  tempLoS = bufferAngLoS[i];
			  bufferAngLoS[i] = bufferAngLoS[j];
			  bufferAngLoS[j] = tempLoS;

			  tempTx = bufferTx[i];
			  bufferTx[i] = bufferTx[j];
			  bufferTx[j] = tempTx;
		  }
	  }
	}
	if(cnt>4){
	for(i=0; i<4; i++)
		TxClosest[cntTx++] = bufferTx[i];
	}
	else
		for(i=0; i<cnt; i++)
			TxClosest[cntTx++] = bufferTx[i];

    // Galileo
    cnt = 0;
    for(ITx=SC_E1;ITx<SC_E36+1;ITx++) {
		RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
		RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
		RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
		UNITV(RelPosN);
		angLoS = acos(VoV(RelPosN,Rhat));
		if(angLoS < angLoSmax) {
		  bufferTx[cnt] = ITx;
		  bufferAngLoS[cnt++] = angLoS;
		}
	}
	for(i=0; i<cnt; ++i) {
	  for (j=i+1; j<cnt; ++j) {
		  if (bufferAngLoS[i] > bufferAngLoS[j]) {
			  tempLoS = bufferAngLoS[i];
			  bufferAngLoS[i] = bufferAngLoS[j];
			  bufferAngLoS[j] = tempLoS;

			  tempTx = bufferTx[i];
			  bufferTx[i] = bufferTx[j];
			  bufferTx[j] = tempTx;
		  }
	  }
	}
	if(cnt>4){
	for(i=0; i<4; i++)
		TxClosest[cntTx++] = bufferTx[i];
	}
	else
		for(i=0; i<cnt; i++)
			TxClosest[cntTx++] = bufferTx[i];

    // Beidou
    cnt = 0;
    for(ITx=SC_C1;ITx<SC_M22+1;ITx++) {
		RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
		RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
		RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
		UNITV(RelPosN);
		angLoS = acos(VoV(RelPosN,Rhat));
		if(angLoS < angLoSmax) {
		  bufferTx[cnt] = ITx;
		  bufferAngLoS[cnt++] = angLoS;
		}
	}
	for(i=0; i<cnt; ++i) {
	  for (j=i+1; j<cnt; ++j) {
		  if (bufferAngLoS[i] > bufferAngLoS[j]) {
			  tempLoS = bufferAngLoS[i];
			  bufferAngLoS[i] = bufferAngLoS[j];
			  bufferAngLoS[j] = tempLoS;

			  tempTx = bufferTx[i];
			  bufferTx[i] = bufferTx[j];
			  bufferTx[j] = tempTx;
		  }
	  }
	}
	if(cnt>4){
	for(i=0; i<4; i++)
		TxClosest[cntTx++] = bufferTx[i];
	}
	else
		for(i=0; i<cnt; i++)
			TxClosest[cntTx++] = bufferTx[i];

	// Sum all
    *TxCnt = cntTx;

    TxNum = (long *) calloc(cntTx,sizeof(long));

	for(ITx=0; ITx<cntTx; ITx++) {
		TxNum[ITx] = TxClosest[ITx];
	}
	return TxNum;
}

void FindVisibleTxMUOS(struct SCType *RX, long *TxNum)
{
	double Rhat[3], RelPosN[3];
	long ITx;
	double angLoS, angLoSmax = 1.221730476396031, angLoSmin = 1.221730476396031; // 70*D2R

	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(RX->PosRxAntN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    // MUOS
    for(ITx=5;ITx<9;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax && angLoS < angLoSmin) {
        	  angLoSmin = angLoS;
        	  *TxNum = ITx;
          }
    }
}

long *FindVisibleTx(struct SCType *RX, long ITxStart, long ITxEnd, long *TxNum, long *TxCnt)
{
	static long buffer[150];
//	double ToS;
	double Rhat[3], RelPosN[3];
	long ITx;
	long cnt = 0;
	double angLoS, angLoSmax; //angHori, magRxPosN;

	angLoSmax = 70 * D2R;
	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(RX->PosRxAntN,Rhat);

	//angHori = Pi - asin(World[EARTH].rad/magRxPosN);
	//printf("%lf\n", angHori*R2D);
    /* Aim at Tx closest to Rx's Zenith */
    for(ITx=ITxStart;ITx<ITxEnd+1;ITx++) {
          RelPosN[0] = SC[ITx].PosN[0] - RX->PosRxAntN[0];
          RelPosN[1] = SC[ITx].PosN[1] - RX->PosRxAntN[1];
          RelPosN[2] = SC[ITx].PosN[2] - RX->PosRxAntN[2];
          UNITV(RelPosN);
          angLoS = acos(VoV(RelPosN,Rhat));
          if(angLoS < angLoSmax) {
        	  //if(LoS < FOV[RX->ID].Width/2) {
        	  //printf("LoS: %lf\n",LoS*R2D);
        	  buffer[cnt++] = ITx;
          }
    }

    *TxCnt = cnt;

    TxNum = (long *) calloc(cnt,sizeof(long));

	for(ITx=0; ITx<cnt; ITx++) {
		TxNum[ITx] = buffer[ITx];
	}
	return TxNum;
}

long *FindVisibleTx_Gs(struct GroundStationType *GS, long ITxStart, long ITxEnd, long *TxNum, long *TxCnt)
{
	static long buffer[25];
	double Rhat[3], RelPosN[3];
	long ITx, i;
	long cnt = 0;
	double LoS;

	/* SC[RxNum] should be Rx Sat */
    CopyUnitV(GS->PosN,Rhat);

    /* Aim at Tx closest to Rx's Zenith */
    for(ITx=ITxStart;ITx<ITxEnd+1;ITx++) {
          for(i=0;i<3;i++)
             RelPosN[i] = SC[ITx].PosN[i] - GS->PosN[i];
          UNITV(RelPosN);
          LoS = acos(VoV(RelPosN,Rhat));
          if(LoS < 70*D2R) {
        	  //printf("LoS: %lf\n",LoS*R2D);
        	  buffer[cnt++] = ITx;
          }
    }
    *TxCnt = cnt;

	TxNum = (long *) calloc(cnt,sizeof(long));

	for(ITx=0; ITx<cnt; ITx++) {
		TxNum[ITx] = buffer[ITx];
	}
	return TxNum;
}

uint8_t FindVisibleSpotbeamTarget(struct SCType *RX, struct SpecularType *Sp, struct TargetAreaType *TARGET)
{
	double vTarget2RxW[3], vTarget2RxT[3];
	double vTarget2SpW[3], vTarget2SpT[3];
	double radRx, radSp, x[3], y[3], z[3];

	CopyUnitV(TARGET->PosW, z);

	vTarget2RxW[0] = RX->PosRxAntW[0] - TARGET->PosW[0];
	vTarget2RxW[1] = RX->PosRxAntW[1] - TARGET->PosW[1];
	vTarget2RxW[2] = RX->PosRxAntW[2] - TARGET->PosW[2];

	vTarget2SpW[0] = Sp->PosW[0] - TARGET->PosW[0];
	vTarget2SpW[1] = Sp->PosW[1] - TARGET->PosW[1];
	vTarget2SpW[2] = Sp->PosW[2] - TARGET->PosW[2];

	// Direction Cosine Matrix from ECEF to Specular frame
	x[0] = vTarget2RxW[0] - VoV(vTarget2RxW, z)*z[0];
	x[1] = vTarget2RxW[1] - VoV(vTarget2RxW, z)*z[1];
	x[2] = vTarget2RxW[2] - VoV(vTarget2RxW, z)*z[2];
	UNITV(x);
	VxV(z,x,y);

	TARGET->CWT[0][0] = x[0];
	TARGET->CWT[0][1] = x[1];
	TARGET->CWT[0][2] = x[2];
	TARGET->CWT[1][0] = y[0];
	TARGET->CWT[1][1] = y[1];
	TARGET->CWT[1][2] = y[2];
	TARGET->CWT[2][0] = z[0];
	TARGET->CWT[2][1] = z[1];
	TARGET->CWT[2][2] = z[2];

	MxV(TARGET->CWT, vTarget2RxW, vTarget2RxT);
	if(vTarget2RxT[2]<0)
		return 1;
	MxV(TARGET->CWT, vTarget2SpW, vTarget2SpT);

	radRx = (SpotBeam->alt - vTarget2RxT[2]) * tan(SpotBeam->halfBW);
	radSp = (SpotBeam->alt - vTarget2SpT[2]) * tan(SpotBeam->halfBW);

	if((radRx - vTarget2RxT[0])>0 && (radSp - sqrt(vTarget2SpT[0]*vTarget2SpT[0]+vTarget2SpT[1]*vTarget2SpT[1]))>0)
	{
//		printf("Rx= %lf, radRx= %lf, Sp= %lf, radSp= %lf\n", vTarget2RxT[0]/1000, radSp/1000, sqrt(vTarget2SpT[0]*vTarget2SpT[0]+vTarget2SpT[1]*vTarget2SpT[1])/1000, radRx/1000);

		return 0; // visible: Inside spotbeam
	}
	else
		return 1; // invisible: Outside spotbeam
}

uint8_t FindVisibleSpotbeam(struct SCType *RX, struct SpecularType *Sp, struct SCType *TX)
{
	double vecTx2RxN[3], vecTx2RxB[3];
	double vecTx2SpN[3], vecTx2SpB[3];
	double radSqRx, radSqSp;
	uint8_t i, cnt;

	cnt = 0;

	for(i=0; i<3; i++){
		vecTx2RxN[i] = RX->PosRxAntN[i] - TX->PosN[i];
		vecTx2SpN[i] = Sp->PosN[i] - TX->PosN[i];
	}
	MxV(Orb[TX->RefOrb].CLN, vecTx2RxN, vecTx2RxB);
	MxV(Orb[TX->RefOrb].CLN, vecTx2SpN, vecTx2SpB);

	for(i=0; i<16; i++){
		radSqRx = (vecTx2RxB[0]-SpotBeam->x[i])*(vecTx2RxB[0]-SpotBeam->x[i])
			      +(vecTx2RxB[1]-SpotBeam->y[i])*(vecTx2RxB[1]-SpotBeam->y[i]);
		if(radSqRx < SpotBeam->r0sq)
		{
			radSqSp = (vecTx2SpB[0]-SpotBeam->x[i])*(vecTx2SpB[0]-SpotBeam->x[i])
				      +(vecTx2SpB[1]-SpotBeam->y[i])*(vecTx2SpB[1]-SpotBeam->y[i]);
			if(radSqSp < SpotBeam->r0sq)
			{
				SpotBeam->visibleSb[cnt++] = i;
			}
		}
	}
	//printf("%lu: ", cnt);
	//for(i=0; i<cnt; i++)
	//	printf("%lu ", SpotBeam->visibleSb[i]);
	//printf("\n");
	return cnt;
}

long FindVisibleSpotbeam_Gs(struct GroundStationType *GS, struct SpecularType *Sp, struct SCType *TX)
{
	double vecTx2RxN[3], vecTx2RxB[3];
	double vecTx2SpN[3], vecTx2SpB[3];
	double radSqRx, radSqSp;
	long i, cnt;

	cnt = 0;

	for(i=0; i<3; i++){
		vecTx2RxN[i] = GS->PosN[i] - TX->PosN[i];
		vecTx2SpN[i] = Sp->PosN[i] - TX->PosN[i];
	}
	MxV(Orb[TX->RefOrb].CLN, vecTx2RxN, vecTx2RxB);
	MxV(Orb[TX->RefOrb].CLN, vecTx2SpN, vecTx2SpB);

	for(i=0; i<16; i++){
		radSqRx = (vecTx2RxB[0]-SpotBeam->x[i])*(vecTx2RxB[0]-SpotBeam->x[i])
			      +(vecTx2RxB[1]-SpotBeam->y[i])*(vecTx2RxB[1]-SpotBeam->y[i]);
		if(radSqRx < SpotBeam->r0sq)
		{
			radSqSp = (vecTx2SpB[0]-SpotBeam->x[i])*(vecTx2SpB[0]-SpotBeam->x[i])
				      +(vecTx2SpB[1]-SpotBeam->y[i])*(vecTx2SpB[1]-SpotBeam->y[i]);
			if(radSqSp < SpotBeam->r0sq)
			{
				SpotBeam->visibleSb[cnt++] = i;
			}
		}
	}
	//printf("%lu: ", cnt);
	//for(i=0; i<cnt; i++)
	//	printf("%lu ", SpotBeam->visibleSb[i]);
	//printf("\n");
	return cnt;
}

void FindPathLength(long RxNum, long SoOpNum, long SpNum, long TxNum)
{
	long i;
	double VecTx2Rx[3], VecSp2Rx[3], VecSp2Tx[3];
	double MagVecSp2Rx, MagVecSp2Tx;

	if(Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].Exists) {
		for(i=0;i<3;i++) {
			VecTx2Rx[i] = SC[TxNum].PosN[i] - SC[RxNum].PosN[i];
			VecSp2Rx[i] = SC[RxNum].PosN[i] - Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PosN[i];
			VecSp2Tx[i] = SC[TxNum].PosN[i] - Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PosN[i];
		}
		MagVecSp2Rx = UNITV(VecSp2Rx);
		MagVecSp2Tx = UNITV(VecSp2Tx);
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenSUR = MAGV(VecTx2Rx);
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenSUR = MagVecSp2Rx + MagVecSp2Tx;
	}
	else {
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenSUR = 0.0;
		Rx[RxNum].SoOp[SoOpNum].Sp[SpNum].PathLenSUR = 0.0;
	}
}

void LocalGeometry(double vSp2RxW[3], struct SCType *RX, struct SoOpType *SoOp, struct SpecularType *Sp, struct SCType *TX)
{
	double vPS[3], vPW[3], vRxS[3], vTxS[3], vP2Rx[3], vP2Tx[3], velRxS[3], velTxS[3];
	double vTemp1[3], vTemp2[3];
	double mSp, mP2RxS, mP2TxS, mVelRxS, mVelTxS;
	long Itheta, Iphi;
	double dtheta, dphi, theta, phi, dAngMin;//, dArea;

	double x[3], y[3], z[3];

	// Specular frame
//	vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
//	vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
//	vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

	mSp = CopyUnitV(Sp->PosW, z);
	Sp->r = mSp;

	x[0] = vSp2RxW[0] - VoV(vSp2RxW, z)*z[0];
	x[1] = vSp2RxW[1] - VoV(vSp2RxW, z)*z[1];
	x[2] = vSp2RxW[2] - VoV(vSp2RxW, z)*z[2];
	UNITV(x);
	VxV(z,x,y);

	// Direction Cosine Matrix from ECEF to Specular frame
	Sp->CSW[0][0] = x[0];
	Sp->CSW[0][1] = x[1];
	Sp->CSW[0][2] = x[2];
	Sp->CSW[1][0] = y[0];
	Sp->CSW[1][1] = y[1];
	Sp->CSW[1][2] = y[2];
	Sp->CSW[2][0] = z[0];
	Sp->CSW[2][1] = z[1];
	Sp->CSW[2][2] = z[2];

	MxV(Sp->CSW, RX->PosRxAntW, vRxS);
	MxV(Sp->CSW, TX->PosW, vTxS);
	MxV(Sp->CSW, RX->VelRxAntW, velRxS);
	MxV(Sp->CSW, TX->VelW, velTxS);

	// Minimum angle step
	dAngMin = SoOp->minSurfRes/mSp;

	// delta theta
	dtheta = asin(Sp->a/mSp)/SoOp->NthetaLocal;
	if(dtheta < dAngMin)
		dtheta = dAngMin;

	// delta phi
	dphi = TwoPi/SoOp->NphiLocal;
	if(dphi < dAngMin)
		dphi = dAngMin;

	Sp->Local = (struct LocalSurfaceType **)calloc(SoOp->NthetaLocal, sizeof(struct LocalSurfaceType));
	for(Itheta=0; Itheta<SoOp->NthetaLocal; Itheta++){
    	Sp->Local[Itheta] = (struct LocalSurfaceType *)calloc(SoOp->NphiLocal, sizeof(struct LocalSurfaceType));
		for(Iphi=0; Iphi<SoOp->NphiLocal; Iphi++){
			theta = Itheta*dtheta;
			phi = Iphi*dphi;
			Sp->Local[Itheta][Iphi].theta = theta;
			Sp->Local[Itheta][Iphi].phi = phi;

			// Differential area over surface
			Sp->Local[Itheta][Iphi].dArea = mSp*mSp*sin(theta)*dtheta*dphi;

			/* Location of surface point in the specular frame */
			vPS[0] = mSp * sin(theta) * cos(phi);
			vPS[1] = mSp * sin(theta) * sin(phi);
			vPS[2] = mSp * cos(theta);

			MTxV(Sp->CSW, vPS, vPW);
			/* ECEF to GRS80 */
			ECEFToGRS80(vPW, &Sp->Local[Itheta][Iphi].lat, &Sp->Local[Itheta][Iphi].lon, &Sp->Local[Itheta][Iphi].alt);

			vP2Rx[0] = vRxS[0] - vPS[0];
			vP2Rx[1] = vRxS[1] - vPS[1];
			vP2Rx[2] = vRxS[2] - vPS[2];

			vP2Tx[0] = vTxS[0] - vPS[0];
			vP2Tx[1] = vTxS[1] - vPS[1];
			vP2Tx[2] = vTxS[2] - vPS[2];

			mP2RxS = UNITV(vP2Rx);
			mP2TxS = UNITV(vP2Tx);

			mVelRxS = VoV(vP2Rx,velRxS);
			mVelTxS = VoV(vP2Tx,velTxS);

			Sp->CLS[0][0] = sin(theta)*cos(phi);
			Sp->CLS[0][1] = sin(theta)*sin(phi);
			Sp->CLS[0][2] = cos(theta);
			Sp->CLS[1][0] = cos(theta)*cos(phi);
			Sp->CLS[1][1] = cos(theta)*sin(phi);
			Sp->CLS[1][2] = -sin(theta);
			Sp->CLS[2][0] = -sin(phi);
			Sp->CLS[2][1] = cos(phi);
			Sp->CLS[2][2] = 0;

			Sp->Local[Itheta][Iphi].range = mP2RxS+mP2TxS;
			Sp->Local[Itheta][Iphi].delay = Sp->Local[Itheta][Iphi].range/Clight - Sp->Local[0][0].delay;
			Sp->Local[Itheta][Iphi].doppler = -(mVelRxS + mVelTxS)/SoOp->Wavelength - Sp->Local[0][0].doppler;
			Sp->Local[Itheta][Iphi].incAngL2Rx = acos(VoV(Sp->CLS[0], vP2Rx));
			Sp->Local[Itheta][Iphi].incAngL2Tx = acos(VoV(Sp->CLS[0], vP2Tx));

			// Scattering Vector q
			vTemp1[0] = vP2Tx[0] + vP2Rx[0];
			vTemp1[1] = vP2Tx[1] + vP2Rx[1];
			vTemp1[2] = vP2Tx[2] + vP2Rx[2];

			MxV(Sp->CLS, vTemp1, vTemp2);
			Sp->Local[Itheta][Iphi].scaVec[0] = TwoPi/SoOp->Wavelength*vTemp2[0];
			Sp->Local[Itheta][Iphi].scaVec[1] = TwoPi/SoOp->Wavelength*vTemp2[1];
			Sp->Local[Itheta][Iphi].scaVec[2] = TwoPi/SoOp->Wavelength*vTemp2[2];

			// Scattering angle
			vP2Tx[0] = -vP2Tx[0];
			vP2Tx[1] = -vP2Tx[1];
			vP2Tx[2] = -vP2Tx[2];
			Sp->Local[Itheta][Iphi].scaAng = 0.5*acos(VoV(vP2Tx,vP2Rx));
		}
	}
}

void LocalGeometryDEM(double vSp2RxW[3], struct SCType *RX, struct SoOpType *SoOp, struct SpecularType *Sp, struct SCType *TX)
{
	double vPS[3], vPW[3], vRxS[3], vTxS[3], vP2Rx[3], vP2Tx[3], velRxS[3], velTxS[3];
	double vTemp1[3], vTemp2[3];
	double mSp, mP2RxS, mP2TxS, mVelRxS, mVelTxS;
	long Itheta, Iphi, thetaSp, phiSp;
	double dtheta, dphi, theta, phi;//, dAngMin;//, dArea;
	double x[3], y[3], z[3];

	uint16_t row, col;
	uint8_t updateFlag = 0;
	double rP, rE, sinLat2, localPathLen, minPathLen;
    double a = 6378137.0;
    double f = 1.0/298.257222101; // GRS80
//    double f = 1.0/298.257223563; // WGS84
    double e2 = f*(2.0-f);
	// Specular frame
//	vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
//	vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
//	vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

	mSp = CopyUnitV(Sp->PosW, z);
	Sp->r = mSp;

	// Direction Cosine Matrix from ECEF to Specular frame
	x[0] = vSp2RxW[0] - VoV(vSp2RxW, z)*z[0];
	x[1] = vSp2RxW[1] - VoV(vSp2RxW, z)*z[1];
	x[2] = vSp2RxW[2] - VoV(vSp2RxW, z)*z[2];
	UNITV(x);
	VxV(z,x,y);

	Sp->CSW[0][0] = x[0];
	Sp->CSW[0][1] = x[1];
	Sp->CSW[0][2] = x[2];
	Sp->CSW[1][0] = y[0];
	Sp->CSW[1][1] = y[1];
	Sp->CSW[1][2] = y[2];
	Sp->CSW[2][0] = z[0];
	Sp->CSW[2][1] = z[1];
	Sp->CSW[2][2] = z[2];

	MxV(Sp->CSW, RX->PosRxAntW, vRxS);
	MxV(Sp->CSW, TX->PosW, vTxS);
	//MxV(Sp->CWS, RX->VelRxAntW, velRxS);
	//MxV(Sp->CWS, TX->VelW, velTxS);

	// Minimum angle step
	//dAngMin = SoOp->minSurfRes/mSp;
	dtheta = SoOp->minSurfRes/mSp;
	// delta theta
	//dtheta = asin(Sp->a/mSp)/SoOp->NthetaLocal;
	//printf("%lf %lf", dAngMin, dtheta);
	//if(dtheta < dAngMin)
	//	dtheta = dAngMin;

	// delta phi
	dphi = TwoPi/SoOp->NphiLocal;
	//if(dphi < dAngMin)
	//	dphi = dAngMin;
	//printf(" %lf\n", dphi);
	/* Update Specular point integrating DEM */
	row = Sp->latGRS*R2Dx60+5399.9999;
	col = Sp->lonGRS*R2Dx60+10799.9999;
	//h_ellip = dem[row][col] + 6371000 - mSp;
	rP = SurfaceStt->height[row][col] + 6371000;
	Sp->altGRS = rP - mSp; // rE == mSp
	// Specular point on GRS80 Ellipsoid and add elevation
	minPathLen = sqrt(vRxS[0]*vRxS[0] + vRxS[1]*vRxS[1] + (vRxS[2]-rP)*(vRxS[2]-rP)) +
		         sqrt(vTxS[0]*vTxS[0] + vTxS[1]*vTxS[1] + (vTxS[2]-rP)*(vTxS[2]-rP));
	Sp->PathLenGRS = minPathLen;

	for(Itheta=0; Itheta<SoOp->NthetaLocal; Itheta++){
		for(Iphi=0; Iphi<SoOp->NphiLocal; Iphi++){
			theta = (Itheta+1)*dtheta;
			phi = Iphi*dphi;
			Sp->Local[Itheta][Iphi].theta = theta;
			Sp->Local[Itheta][Iphi].phi = phi;

			// Differential area over surface
			//Sp->Local[Itheta][Iphi].dArea = mSp*mSp*sin(theta)*dtheta*dphi;

			/* Location of surface point in the specular frame */
			vPS[0] = mSp * sin(theta) * cos(phi);
			vPS[1] = mSp * sin(theta) * sin(phi);
			vPS[2] = mSp * cos(theta);

			/* Specular frame to ECEF */
			MTxV(Sp->CSW, vPS, vPW);

			/* ECEF to GRS80 */
			ECEFToGRS80(vPW, &Sp->Local[Itheta][Iphi].lat, &Sp->Local[Itheta][Iphi].lon, &Sp->Local[Itheta][Iphi].alt);

			/* Integrate DEM */
			row = Sp->Local[Itheta][Iphi].lat*R2Dx60+5399.9999;
			col = Sp->Local[Itheta][Iphi].lon*R2Dx60+10799.9999;
			sinLat2 = sin(Sp->Local[Itheta][Iphi].lat)*sin(Sp->Local[Itheta][Iphi].lat);
			rE = a*sqrt((1-e2*(2-e2)*sinLat2) / (1-e2*sinLat2));
			Sp->Local[Itheta][Iphi].alt = SurfaceStt->height[row][col] + 6371000 - rE;

			/* GRS80 to ECEF */
	    	GRS80ToECEF(Sp->Local[Itheta][Iphi].lat, Sp->Local[Itheta][Iphi].lon, Sp->Local[Itheta][Iphi].alt, vPW);

			/* ECEF to Specular frame */
	    	MxV(Sp->CSW, vPW, vPS);

			vP2Rx[0] = vRxS[0] - vPS[0];
			vP2Rx[1] = vRxS[1] - vPS[1];
			vP2Rx[2] = vRxS[2] - vPS[2];

			vP2Tx[0] = vTxS[0] - vPS[0];
			vP2Tx[1] = vTxS[1] - vPS[1];
			vP2Tx[2] = vTxS[2] - vPS[2];

			mP2RxS = UNITV(vP2Rx);
			mP2TxS = UNITV(vP2Tx);

			/* Find minimum reflection path length */
			localPathLen = mP2RxS + mP2TxS;

			if(localPathLen < minPathLen){
				minPathLen = localPathLen;
				thetaSp = Itheta;
				phiSp = Iphi;
				updateFlag = 1;
			}
/*
			mVelRxS = VoV(vP2Rx,velRxS);
			mVelTxS = VoV(vP2Tx,velTxS);

			Sp->CLS[0][0] = sin(theta)*cos(phi);
			Sp->CLS[0][1] = sin(theta)*sin(phi);
			Sp->CLS[0][2] = cos(theta);
			Sp->CLS[1][0] = cos(theta)*cos(phi);
			Sp->CLS[1][1] = cos(theta)*sin(phi);
			Sp->CLS[1][2] = -sin(theta);
			Sp->CLS[2][0] = -sin(phi);
			Sp->CLS[2][1] = cos(phi);
			Sp->CLS[2][2] = 0;

			Sp->Local[Itheta][Iphi].range = localPathLen;
			Sp->Local[Itheta][Iphi].delay = Sp->Local[Itheta][Iphi].range/Clight - Sp->Local[0][0].delay;
			Sp->Local[Itheta][Iphi].doppler = -(mVelRxS + mVelTxS)/SoOp->Wavelength - Sp->Local[0][0].doppler;
			Sp->Local[Itheta][Iphi].incAngL2Rx = acos(VoV(Sp->CLS[0], vP2Rx));
			Sp->Local[Itheta][Iphi].incAngL2Tx = acos(VoV(Sp->CLS[0], vP2Tx));

			// Scattering Vector q
			vTemp1[0] = vP2Tx[0] + vP2Rx[0];
			vTemp1[1] = vP2Tx[1] + vP2Rx[1];
			vTemp1[2] = vP2Tx[2] + vP2Rx[2];

			MxV(Sp->CLS, vTemp1, vTemp2);
			Sp->Local[Itheta][Iphi].scaVec[0] = TwoPi/SoOp->Wavelength*vTemp2[0];
			Sp->Local[Itheta][Iphi].scaVec[1] = TwoPi/SoOp->Wavelength*vTemp2[1];
			Sp->Local[Itheta][Iphi].scaVec[2] = TwoPi/SoOp->Wavelength*vTemp2[2];

			// Scattering angle
			vP2Tx[0] = -vP2Tx[0];
			vP2Tx[1] = -vP2Tx[1];
			vP2Tx[2] = -vP2Tx[2];
			Sp->Local[Itheta][Iphi].scaAng = 0.5*acos(VoV(vP2Tx,vP2Rx));
*/
		}
	}
	/* Update new specular point integrating DEM */
	Sp->PathLenSUR = minPathLen;
	if(updateFlag){
		Sp->latSUR = Sp->Local[thetaSp][phiSp].lat;
		Sp->lonSUR = Sp->Local[thetaSp][phiSp].lon;
		Sp->altSUR = Sp->Local[thetaSp][phiSp].alt;
		//printf("(%ld, %ld) = %lf, %lf, %lf, snell: %lf\n",
		//		thetaSp, phiSp,
		//		(Sp->latGRS-Sp->latSUR)*R2D,
		//		(Sp->lonGRS-Sp->lonSUR)*R2D,
		//		(Sp->altGRS-Sp->altSUR)*R2D,
		//		(Sp->Local[thetaSp][phiSp].incAngL2Rx - Sp->Local[thetaSp][phiSp].incAngL2Tx)*R2D);
	}
	else{
		Sp->latSUR = Sp->latGRS;
		Sp->lonSUR = Sp->lonGRS;
		Sp->altSUR = Sp->altGRS;
		//printf("same\n");
	}
}

/*
long *CntValidSpecularPt(long RxNum, long SoOpNum, long RoiNum)
{
	long *ValidTxNum;
	long ISp, i;
	long j = 0;
	long cnt = 0;
	double theta;
	double unitR_Roi[3], unitR_Sp[3];

	ValidTxNum = calloc(Rx[RxNum].SoOp[SoOpNum].NSp,sizeof(long));

	// ECEF to ECI
	MTxV(World[EARTH].CWN,GroundStation[RoiNum].PosW,unitR_Roi);
	UNITV(unitR_Roi);

	for(ISp=0;ISp<Rx[RxNum].SoOp[SoOpNum].NSp;ISp++){
		if(Rx[RxNum].SoOp[SoOpNum].Sp[ISp].Exists){
			for(i=0;i<3;i++){
				unitR_Sp[i] = Rx[RxNum].SoOp[SoOpNum].Sp[ISp].PosN[i];
			}
			UNITV(unitR_Sp);
			theta = acos(VoV(unitR_Roi,unitR_Sp))*R2D;
			//printf("%f\n", theta);
//			if(theta<30){
//			if(theta<0.00078393){
			if(theta<0.0156){
				ValidTxNum[j] = Rx[RxNum].SoOp[SoOpNum].Sp[ISp].TxNum;
				cnt++;
				j++;
	            if(j == Rx[RxNum].SoOp[SoOpNum].NSp)
	            	break;
			}
		}
	}
	Rx[RxNum].SoOp[SoOpNum].NValidSp = cnt;
	return ValidTxNum;
}

void CalcValidRxPwr(long RoiNum, double Reflectivity)
{
	long IRx, ISoOp, ISp;
	long *ValidTxNum;

	for(IRx=0; IRx<NRx; IRx++){
		for(ISoOp=0; ISoOp<Rx[IRx].NSoOp; ISoOp++){
			ValidTxNum = CntValidSpecularPt(IRx, ISoOp, RoiNum);
			for(ISp=0; ISp<Rx[IRx].SoOp[ISoOp].NSp; ISp++){
				if(ValidTxNum[ISp]==0){
					// do nothing
				}
				else{
					FindPathLength(IRx, ISoOp, ISp, ValidTxNum[ISp]);
					CalcRxPwr(IRx, ISoOp, ISp, ValidTxNum[ISp], Reflectivity);
					Rx[IRx].SoOp[ISoOp].TotalRxPwrD += Rx[IRx].SoOp[ISoOp].Sp[ISp].RxPwrD;
					Rx[IRx].SoOp[ISoOp].TotalRxPwrR += Rx[IRx].SoOp[ISoOp].Sp[ISp].RxPwrR;
				}
			}
		}
	}
}
*/


void sph2car(double theta, double phi, double AR, double AT, double AP, double axis[3])
{
	double st, ct, sp, cp;

	st = sin(theta);
	ct = cos(theta);
	sp = sin(phi);
	cp = cos(phi);

	axis[0]=AR*st*cp + AT*ct*cp - AP*sp;
	axis[1]=AR*st*sp + AT*ct*sp + AP*cp;
	axis[2]=AR*ct - AT*st;
}

void tanUnitVec(double CB1[3][3], double CB2[3][3], double u[3], long pol1, long pol2,
				  double u1p1[3], double u1p2[3], double u2p1[3], double u2p2[3])
{
	double u1[3], u2[3];
	double th1, th2, ph1, ph2;
	double AR, AT, AP;
	double temp[3];
	double u1ph[3], u2ph[3]; // H-pol
	double u1th[3], u2th[3]; // V-pol
	double u1X[3], u2X[3]; // X-pol
	double u1Y[3], u2Y[3]; // Y-pol
	double u1R[3], u2R[3]; // R-pol
	double u1L[3], u2L[3]; // L-pol
	long i;

	MxV(CB1,u,u1);
	MxV(CB2,u,u2);

	th1 = acos(u1[2]);
	ph1 = atan2(u1[1], u1[0]);
	th2 = acos(-u2[2]);
	ph2 = atan2(-u2[1], -u2[0]);

	// unit vectors in phi direction [H-pol]
	AR = 0; AT = 0; AP = 1;
	// phi vector of frame 1 in reference frame
	sph2car(th1, ph1, AR, AT, AP, temp);
	MTxV(CB1, temp, u1ph);
	// phi vector of frame 2 in reference frame
	sph2car(th2, ph2, AR, AT, AP, temp);
	MTxV(CB2, temp, u2ph);

	// unit vectors in theta direction [V-pol]
	AR = 0; AT = 1; AP = 0;
	// theta vector of frame 1 in reference frame
	sph2car(th1, ph1, AR, AT, AP, temp);
	MTxV(CB1, temp, u1th);
	// theta vector of frame 2 in reference frame
	sph2car(th2, ph2, AR, AT, AP, temp);
	MTxV(CB2, temp, u2th);

	// unit vectors in Ludwig basis (reference: Y-pol)
	// Y-axis vector of frame 1 in reference frame
	AR = 0; AT = sin(ph1); AP = cos(ph1);
	sph2car(th1, ph1, AR, AT, AP, temp);
	MTxV(CB1, temp, u1Y);

	// Y-axis vector of frame 2 in reference frame
	AR = 0; AT = sin(ph2); AP = cos(ph2);
	sph2car(th2, ph2, AR, AT, AP, temp);
	MTxV(CB2, temp, u2Y);

	// unit vectors in Ludwig basis (cross: X-pol)
	// X-axis vector of frame 1 in reference frame
	AR = 0; AT = cos(ph1); AP = -sin(ph1);
	sph2car(th1, ph1, AR, AT, AP, temp);
	MTxV(CB1, temp, u1X);

	// X-axis vector of frame 2 in reference frame
	AR = 0; AT = cos(ph2); AP = -sin(ph2);
	sph2car(th2, ph2, AR, AT, AP, temp);
	MTxV(CB2, temp, u2X);

	// unit vectors in circular polarization (R-pol and L-pol)
	for(i=0; i<3; i++){
		// frame 1
		u1R[i] = 1/sqrt(2) * (u1X[i]-I*u1Y[i]); // port 1 RHCP
		u1L[i] = 1/sqrt(2) * (u1X[i]+I*u1Y[i]); // port 2 LHCP
		// frame 2
		u2R[i] = 1/sqrt(2) * (u2X[i]-I*u2Y[i]); // port 1 RHCP
		u2L[i] = 1/sqrt(2) * (u2X[i]+I*u2Y[i]); // port 2 LHCP
	}

	switch(pol1){
		case POL_H:
			for(i=0; i<3; i++){
				u1p1[i] = u1ph[i];
				u1p2[i] = u1th[i];
			}
			break;
		case POL_V:
			for(i=0; i<3; i++){
				u1p1[i] = u1th[i];
				u1p2[i] = u1ph[i];
			}
			break;
		case POL_X:
			for(i=0; i<3; i++){
				u1p1[i] = u1X[i];
				u1p2[i] = u1Y[i];
			}
			break;
		case POL_Y:
			for(i=0; i<3; i++){
				u1p1[i] = u1Y[i];
				u1p2[i] = u1X[i];
			}
			break;
		case POL_R:
			for(i=0; i<3; i++){
				u1p1[i] = u1R[i];
				u1p2[i] = u1L[i];
			}
			break;
		case POL_L:
			for(i=0; i<3; i++){
				u1p1[i] = u1L[i];
				u1p2[i] = u1R[i];
			}
			break;
	}

	switch(pol2){
		case POL_H:
			for(i=0; i<3; i++){
				u2p1[i] = u2ph[i];
				u2p2[i] = u2th[i];
			}
			break;
		case POL_V:
			for(i=0; i<3; i++){
				u2p1[i] = u2th[i];
				u2p2[i] = u2ph[i];
			}
			break;
		case POL_X:
			for(i=0; i<3; i++){
				u2p1[i] = u2X[i];
				u2p2[i] = u2Y[i];
			}
			break;
		case POL_Y:
			for(i=0; i<3; i++){
				u2p1[i] = u2Y[i];
				u2p2[i] = u2X[i];
			}
			break;
		case POL_R:
			for(i=0; i<3; i++){
				u2p1[i] = u2R[i];
				u2p2[i] = u2L[i];
			}
			break;
		case POL_L:
			for(i=0; i<3; i++){
				u2p1[i] = u2L[i];
				u2p2[i] = u2R[i];
			}
			break;
	}
}

void calcMueller(double u[2][2], double U[4][4])
{
	double u11, u12, u21, u22;

	// Get the elements of the 2 x 2 input matrix
	u11 = u[0][0]; u12 = u[0][1];
	u21 = u[1][0]; u22 = u[1][1];

	// Calculate the elements of 4 x 4 Mueller matrix
	U[0][0] = abs(u11)^2;
	U[0][1] = abs(u12)^2;
	U[0][2] = creal(u11 * conj(u12));
	U[0][3] = -cimag(u11 * conj(u12));

	U[1][0] = abs(u21)^2;
	U[1][1] = abs(u22)^2;
	U[1][2] = creal(u21*conj(u22));
	U[1][3] = -cimag(u21*conj(u22));

	U[2][0] = 2 * creal(u11 * conj(u21));
	U[2][1] = 2 * creal(u12 * conj(u22));
	U[2][2] = creal(u11 * conj(u22) + u12 * conj(u21));
	U[2][3] = -cimag(u11 * conj(u22) - u12 * conj(u21));

	U[3][0] = 2 * cimag(u11 * conj(u21));
	U[3][1] = 2 * cimag(u12 * conj(u22));
	U[3][2] = cimag(u11 * conj(u22) + u12 * conj(u21));
	U[3][3] = creal(u11 * conj(u22) - u12 * conj(u21));
}

void ObsGlobeMOG(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile[6][16];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	struct SoilType *Soil;
	double latT, lonT, altT, latR, lonR, altR;
	long *TxNum = NULL;

	uint8_t iRx, iTx, iSp;
	uint16_t rowSt, colSt, rowDy, colDy;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3], temp;
	double vNadir[3] = {0.0, 0.0, 1.0};

	static long *SpCntTotal;
	long SpCntMUOS, SpCntOrbcomm, SpCntGPS, SpCntGlonass, SpCntGalileo, SpCntBeidou;

    if (First) {
    	char filename[20], TxName[10];
    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(iTx=0; iTx<NTx; iTx++){
        	if(iTx == 0)
        		sprintf(TxName, "MUOS");
        	else if(iTx == 1)
        		sprintf(TxName, "Orbcomm");
        	else if(iTx == 2)
        		sprintf(TxName, "GPS");
        	else if(iTx == 3)
        		sprintf(TxName, "Glonass");
        	else if(iTx == 4)
        		sprintf(TxName, "Galileo");
        	else
        		sprintf(TxName, "Beidou");

			for(iRx=0; iRx<NRx; iRx++){
				sprintf(filename,"SoOp_%s_Rx%d.42", TxName, iRx);
				outFile[iTx][iRx] = FileOpen(InOutPath,filename,"w");
			}
        }
    }
    for(iRx=0; iRx<NRx; iRx++){
    	RX = &SC[Rx[iRx].SC];

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
    	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
    	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

    	//RX->VelRxAntN[0] = RX->VelN[0] + RX->B[ANTBODY].vn[0];
       	//RX->VelRxAntN[1] = RX->VelN[1] + RX->B[ANTBODY].vn[1];
       	//RX->VelRxAntN[2] = RX->VelN[2] + RX->B[ANTBODY].vn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
       	/* ECI to ECEF of Receiver's Velocity */
    	//MxV(World[EARTH].CWN,RX->VelRxAntN,RX->VelRxAntW);

    	/* ECEF to GRS80 */
    	//ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);

    	TxNum = FindVisibleTxSG(RX, TxNum, SpCntTotal); // max: 186
    	//TxNum = FindVisibleTxSoOp(RX, TxNum, SpCntTotal);

    	//SoOp
    	SpCntMUOS = 0;
    	SpCntOrbcomm = 0;
    	SpCntGPS = 0;
    	SpCntGlonass = 0;
    	SpCntGalileo = 0;
    	SpCntBeidou = 0;

    	for(iSp=0; iSp<*SpCntTotal; iSp++){
    		// GNSS
    		if(TxNum[iSp]>137){
    			SoOp = &Rx[iRx].SoOp[TX_BEIDOU];
				Sp = &SoOp->Sp[SpCntBeidou++];
				iTx = TX_BEIDOU;
    		}
    		else if(TxNum[iSp]>55 && TxNum[iSp]<87){
    			SoOp = &Rx[iRx].SoOp[TX_GPS];
    			Sp = &SoOp->Sp[SpCntGPS++];
    		    iTx = TX_GPS;
    		}
    		else if(TxNum[iSp]>86 && TxNum[iSp]<112){
				SoOp = &Rx[iRx].SoOp[TX_GLONASS];
				Sp = &SoOp->Sp[SpCntGlonass++];
				iTx = TX_GLONASS;
    		}
    		else if(TxNum[iSp]>111 && TxNum[iSp]<138){
    			SoOp = &Rx[iRx].SoOp[TX_GALILEO];
    			Sp = &SoOp->Sp[SpCntGalileo++];
    		    iTx = TX_GALILEO;
    		}
    		// SoOp
    		else if(TxNum[iSp]>15 && TxNum[iSp]<20){
    			SoOp = &Rx[iRx].SoOp[TX_MUOS];
    			Sp = &SoOp->Sp[SpCntMUOS++];
    		    iTx = TX_MUOS;
    		}
    		else if(TxNum[iSp]>19 && TxNum[iSp]<56){
    			SoOp = &Rx[iRx].SoOp[TX_ORBCOMM];
    			Sp = &SoOp->Sp[SpCntOrbcomm++];
    		    iTx = TX_ORBCOMM;
    		}

    		TX = &SC[TxNum[iSp]];


    		// ECI to ECEF of Transmitter's Position and Velocity
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
    		//MxV(World[EARTH].CWN,TX->VelN,TX->VelW);

    		/* ECEF to GRS80 */
    		//ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);

    		vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
    		vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
    		vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

    		/* ECEF to GRS80 */
    		ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

    		/* Update Specular Point Integrating DEM */
    		LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);

			// Coverage of CONUS
			if(Sp->latSUR<0.855211333477221 && Sp->latSUR>0.427164946867716
			&& Sp->lonSUR>-2.180570734210416 && Sp->lonSUR<-1.167593205124682){
//    		if(Sp->lonSUR>-2.180570734210416 && Sp->lonSUR<-1.167593205124682){
    			fprintf(outFile[iTx][iRx], "%d %lf %ld %lf %lf %lf\n",
						RX->Ncycle, AbsTime, TxNum[iSp], Sp->latSUR, Sp->lonSUR, Sp->altSUR);
			}

    		// GRS80 to ECEF
    		GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

    		// ECEF to ECI
    		MTxV(World[EARTH].CWN,Sp->PosW,Sp->PosN);
    		Sp->Exists = 1;
/*
			//signalVec(TxNum[ISp]);
			//Sp->Exists = 1;

			// Calculate Signal Angles of Arrival
			vRx2TxN[0] = TX->PosN[0] - RX->PosRxAntN[0]; // for Sky-view antenna look angle & slant range
			vRx2TxN[1] = TX->PosN[1] - RX->PosRxAntN[1]; // for Sky-view antenna look angle & slant range
			vRx2TxN[2] = TX->PosN[2] - RX->PosRxAntN[2]; // for Sky-view antenna look angle & slant range

			Sp->PathLenD = UNITV(vRx2TxN);
			// Sky-view (Direct signal) antenna look angle
//    		MxV(RX->B[ANTBODY].CN, vRx2TxN, vRx2TxA); // opposite direction vector of direct signal in Antenna body frame
			MxV(Orb[iRx].CLN, vRx2TxN, vRx2TxA); // opposite direction vector of reflected signal in Antenna body frame

			vRx2TxA[1] = -vRx2TxA[1]; // 180 degree rotation about x-axis
			vRx2TxA[2] = -vRx2TxA[2]; // 180 degree rotation about x-axis

			Sp->lookAngElD = acos(VoV(vRx2TxA,vNadir));
			Sp->lookAngAzD = atan2(vRx2TxA[1],vRx2TxA[0]);

    		// Check reflected signal angle of arrival
    		vRx2SpN[0] = Sp->PosN[0] - RX->PosRxAntN[0]; // for Earth-view antenna look angle
    		vRx2SpN[1] = Sp->PosN[1] - RX->PosRxAntN[1]; // for Earth-view antenna look angle
    		vRx2SpN[2] = Sp->PosN[2] - RX->PosRxAntN[2]; // for Earth-view antenna look angle
    		UNITV(vRx2SpN);
    		// Earth-view (Reflected signal) antenna look angle
//    		MxV(RX->B[ANTBODY].CN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame
    		MxV(Orb[iRx].CLN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame

    		Sp->lookAngElR = acos(VoV(vRx2SpA,vNadir));
   			Sp->lookAngAzR = atan2(vRx2SpA[1],vRx2SpA[0]);
			// Antenna Gain
			calcAntPattern(SoOp, Sp);

			// Calculate the First Fresnel Zone
			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0]; // for Fresnel zone
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1]; // for Fresnel zone
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2]; // for Fresnel zone
			// ENU to AE
			temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
			vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
			vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
			vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

			temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
			Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
			Sp->elev = atan2(vSp2RxE[2],temp);
			Sp->incAng = HalfPi-Sp->elev;
			// First Fresnel zone
			Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
			Sp->a = Sp->b / sin(Sp->elev);
			Sp->area = Pi* Sp->a * Sp->b;

			// Direct Signal
	//		nSpotBeam = FindVisibleSpotbeam(RX, Sp, TX);
			Sp->RxPwrD = Sp->AntGainS + Tx[TX_MUOS].EIRP + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926; //+ 10*log10(nSpotBeam);
			Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;// - 10*log10(SoOp->CorIntTime);

			// Reflected Signal
			Soil = Sp->Soil;
			R = Sp->Reflect;
			LL2EzRC1km(Sp->latSUR, Sp->lonSUR, &rowSt, &colSt);
			LL2EzRC9km(Sp->latSUR, Sp->lonSUR, &rowDy, &colDy);
			//printf("(%d,%d), (%d,%d)\n", rowSt, colSt, rowDy, colDy);

			//printf("(%d,%d) S:%f C:%f rho_b:%f, (%d,%d), RZSM:%f, T:%f\n", rowSt, colSt,
			//		SurfaceStt->sand[rowSt][colSt], SurfaceStt->clay[rowSt][colSt], SurfaceStt->rho_b[rowSt][colSt], rowDy, colDy,
	 		//		SurfaceDyn->RZSM[rowDy][colDy], SurfaceDyn->soilTemp1[rowDy][colDy]);
			if(rowSt>14615 || colSt>34703 || rowDy>1623 || colDy>3855 ||
			   isnan(SurfaceStt->sand[rowSt][colSt]) || isnan(SurfaceDyn->RZSM[rowDy][colDy]))
			{
				// do nothing
			}
			else{
				calcDielPeplinski(rowSt, colSt, rowDy, colDy, Tx[TX_MUOS].Freq, Soil);
				calcReflectivity(Soil->e_c, R, Sp->incAng);

				R->RxPwrR = Sp->AntGainE + Tx[TX_MUOS].EIRP + 10*log10(R->gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenSUR) - 21.984197280441926; // + 10*log10(nSpotBeam * R->gammaLR) ;
				R->SnrR = R->RxPwrR - SoOp->NoisePwrE;// - 10*log10(SoOp->CorIntTime);

				fprintf(outFile[iTx][iRx], "%lf %ld %lf %lf %lf %lf %lf %lf %lf\n",
						AbsTime, TxNum[iSp], Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->SnrD, R->SnrR, R->gammaLR, Sp->area);
			}
			*/
    	}
    	Rx[iRx].SoOp[TX_MUOS].SpCnt = SpCntMUOS;
    	Rx[iRx].SoOp[TX_ORBCOMM].SpCnt = SpCntOrbcomm;
    	Rx[iRx].SoOp[TX_GPS].SpCnt = SpCntGPS;
    	Rx[iRx].SoOp[TX_GLONASS].SpCnt = SpCntGlonass;
    	Rx[iRx].SoOp[TX_GALILEO].SpCnt = SpCntGalileo;
    	Rx[iRx].SoOp[TX_BEIDOU].SpCnt = SpCntBeidou;
    	free(TxNum);
    }
}

void ObsGlobeMUOS(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile1[5], *outFile2[5];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	struct SoilType *Soil;
	double latT, lonT, altT, latR, lonR, altR;
	long TxNum;

	uint8_t iRx;
	uint16_t rowSt, colSt, rowDy, colDy;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3], temp;
	double vNadir[3] = {0.0, 0.0, 1.0};

	static long *SpCntTotal;
	long SpCntMUOS;

    if (First) {
    	char filename[20];
    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(iRx=0; iRx<NRx; iRx++){
			sprintf(filename,"SoOpMUOS_Rx%d.42",iRx);
			outFile1[iRx] = FileOpen(InOutPath,filename,"w");
			// File for Positions
			//sprintf(filename,"PosMUOS.42");
			//outFile2 = FileOpen(InOutPath,filename,"w");
			// File for Angle of Arrival
			sprintf(filename,"AoAMUOS_Rx%d.42", iRx);
			outFile2[iRx] = FileOpen(InOutPath,filename,"w");
        }
    }
    for(iRx=0; iRx<NRx; iRx++){
    	RX = &SC[Rx[iRx].SC];
    	TxNum = 0;

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
    	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
    	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

    	RX->VelRxAntN[0] = RX->VelN[0] + RX->B[ANTBODY].vn[0];
       	RX->VelRxAntN[1] = RX->VelN[1] + RX->B[ANTBODY].vn[1];
       	RX->VelRxAntN[2] = RX->VelN[2] + RX->B[ANTBODY].vn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
       	/* ECI to ECEF of Receiver's Velocity */
    	MxV(World[EARTH].CWN,RX->VelRxAntN,RX->VelRxAntW);

    	/* ECEF to GRS80 */
    	ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);

    	//TxNum = FindVisibleTx(RX, 1, 5, TxNum, SpCntTotal); // max: 186
    	//TxNum = FindVisibleTxSoOp(RX, TxNum, SpCntTotal);
    	FindVisibleTxMUOS(RX, &TxNum);

    	//SoOp
    	SpCntMUOS = 0;

    	if(TxNum > 0)
    	{
    		// SoOp
    		SoOp = &Rx[iRx].SoOp[TX_MUOS];
    		Sp = &SoOp->Sp[SpCntMUOS++];
    		TX = &SC[TxNum];

    		// ECI to ECEF of Transmitter's Position and Velocity
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
    		MxV(World[EARTH].CWN,TX->VelN,TX->VelW);

    		/* ECEF to GRS80 */
    		ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);

    		vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
    		vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
    		vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

    		/* ECEF to GRS80 */
    		ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

    		/* Update Specular Point Integrating DEM */
    		LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);

    		/* GRS80 to ECEF */
    		GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

    		// ECEF to ECI
    		MTxV(World[EARTH].CWN,Sp->PosW,Sp->PosN);

			//signalVec(TxNum[ISp]);
			//Sp->Exists = 1;

			/* Calculate Signal Angles of Arrival */
			vRx2TxN[0] = TX->PosN[0] - RX->PosRxAntN[0]; // for Sky-view antenna look angle & slant range
			vRx2TxN[1] = TX->PosN[1] - RX->PosRxAntN[1]; // for Sky-view antenna look angle & slant range
			vRx2TxN[2] = TX->PosN[2] - RX->PosRxAntN[2]; // for Sky-view antenna look angle & slant range

			Sp->PathLenD = UNITV(vRx2TxN);
			// Sky-view (Direct signal) antenna look angle
//    		MxV(RX->B[ANTBODY].CN, vRx2TxN, vRx2TxA); // opposite direction vector of direct signal in Antenna body frame
			MxV(Orb[iRx].CLN, vRx2TxN, vRx2TxA); // opposite direction vector of reflected signal in Antenna body frame

			vRx2TxA[1] = -vRx2TxA[1]; // 180 degree rotation about x-axis
			vRx2TxA[2] = -vRx2TxA[2]; // 180 degree rotation about x-axis

			Sp->lookAngElD = acos(VoV(vRx2TxA,vNadir));
			Sp->lookAngAzD = atan2(vRx2TxA[1],vRx2TxA[0]);

    		// Check reflected signal angle of arrival
    		vRx2SpN[0] = Sp->PosN[0] - RX->PosRxAntN[0]; // for Earth-view antenna look angle
    		vRx2SpN[1] = Sp->PosN[1] - RX->PosRxAntN[1]; // for Earth-view antenna look angle
    		vRx2SpN[2] = Sp->PosN[2] - RX->PosRxAntN[2]; // for Earth-view antenna look angle
    		UNITV(vRx2SpN);
    		// Earth-view (Reflected signal) antenna look angle
//    		MxV(RX->B[ANTBODY].CN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame
    		MxV(Orb[iRx].CLN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame

    		Sp->lookAngElR = acos(VoV(vRx2SpA,vNadir));
   			Sp->lookAngAzR = atan2(vRx2SpA[1],vRx2SpA[0]);
			// Antenna Gain
			calcAntPattern(SoOp, Sp);

			/* Calculate the First Fresnel Zone */
			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0]; // for Fresnel zone
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1]; // for Fresnel zone
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2]; // for Fresnel zone
			// ENU to AE
			temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
			vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
			vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
			vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

			temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
			Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
			Sp->elev = atan2(vSp2RxE[2],temp);
			Sp->incAng = HalfPi-Sp->elev;
			// First Fresnel zone
			Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
			Sp->a = Sp->b / sin(Sp->elev);
			Sp->area = Pi* Sp->a * Sp->b;

			// Direct Signal
	//		nSpotBeam = FindVisibleSpotbeam(RX, Sp, TX);
			Sp->RxPwrD = Sp->AntGainS + Tx[TX_MUOS].EIRP + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926; //+ 10*log10(nSpotBeam);
			Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;// - 10*log10(SoOp->CorIntTime);

			// Reflected Signal
			Soil = Sp->Soil;
			R = Sp->Reflect;
			LL2EzRC1km(Sp->latSUR, Sp->lonSUR, &rowSt, &colSt);
			LL2EzRC9km(Sp->latSUR, Sp->lonSUR, &rowDy, &colDy);
	//		printf("(%d,%d) S:%f C:%f rho_b:%f, (%d,%d), RZSM:%f, T:%f\n", rowSt, colSt,
	//				SurfaceStt->sand[rowSt][colSt], SurfaceStt->clay[rowSt][colSt], SurfaceStt->rho_b[rowSt][colSt], rowDy, colDy,
	// 				SurfaceDyn->RZSM[rowDy][colDy], SurfaceDyn->soilTemp1[rowDy][colDy]);
			if(isnan(SurfaceStt->sand[rowSt][colSt]) || isnan(SurfaceStt->clay[rowSt][colSt]) || isnan(SurfaceStt->rho_b[rowSt][colSt])
			|| isnan(SurfaceDyn->RZSM[rowDy][colDy]) || isnan(SurfaceDyn->soilTemp1[rowDy][colDy]))
			{
				// do nothing
			}
			else{
				calcDielPeplinski(rowSt, colSt, rowDy, colDy, Tx[TX_MUOS].Freq, Soil);
				calcReflectivity(Soil->e_c, R, Sp->incAng);

				R->RxPwrR = Sp->AntGainE + Tx[TX_MUOS].EIRP + 10*log10(R->gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenSUR) - 21.984197280441926; // + 10*log10(nSpotBeam * R->gammaLR) ;
				R->SnrR = R->RxPwrR - SoOp->NoisePwrE;// - 10*log10(SoOp->CorIntTime);

				fprintf(outFile1[iRx], "%lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %f %f %f %f %f\n",
						AbsTime, TxNum, Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->SnrD, R->SnrR, R->gammaLR, R->RcLR,
						Sp->a, Sp->b, Sp->az, Sp->elev, Sp->area,
						Soil->m_v, Soil->T, Soil->sand, Soil->clay, Soil->rho_b);
			}

	//		fprintf(outFile2, "%lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	//				AbsTime, TxNum, latT, lonT, altT, latR, lonR, altR, Sp->latSUR, Sp->lonSUR, Sp->altSUR,
	//				TX->VelW[0], TX->VelW[1], TX->VelW[2], RX->VelRxAntW[0], RX->VelRxAntW[1], RX->VelRxAntW[2]);
			fprintf(outFile2[iRx], "%lf %ld %lf %lf %lf %lf\n",
					AbsTime, TxNum, Sp->lookAngAzD, Sp->lookAngElD, Sp->lookAngAzR, Sp->lookAngElR);
    	}
    	Rx[iRx].SoOp[TX_MUOS].SpCnt = SpCntMUOS;
    }
}

void ObsCONUS_MUOS(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile1[5], *outFile2[5];//, *outFile3[5];
    static uint8_t flagRecording[5] = {0};
    static uint8_t flagRecording2[5] = {0};

    static long recordingDoyLast[5] = {0}, recordingDoy[5] = {0}, recordingCnt[5] = {0};
    static long recordingDoyLast2[5] = {0}, recordingDoy2[5]= {0}, recordingCnt2[5] = {0}, snowArcID[5] = {0}, snowArcDuration[5] = {0};

    uint8_t validArc = 0;

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	struct SoilType *Soil;
	double latT, lonT, altT, latR, lonR, altR;
	long TxNum;

	uint8_t iRx, landMask;
	uint16_t rowSt, colSt, rowDy, colDy;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3], temp;
	double vNadir[3] = {0.0, 0.0, 1.0};

	static long *SpCntTotal;
	long SpCntMUOS;

    if (First) {
    	char filename[20];
    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(iRx=0; iRx<NRx; iRx++){
			sprintf(filename,"CONUS_MUOS_Case%d.42",iRx);
			outFile1[iRx] = FileOpen(InOutPath,filename,"w");
			sprintf(filename,"Snow_MUOS_Case%d.42",iRx);
			outFile2[iRx] = FileOpen(InOutPath,filename,"w");
			//sprintf(filename,"PosMUOS_Rx%d.42", iRx);
			//outFile3[iRx] = FileOpen(InOutPath,filename,"w");
        }
    }
    for(iRx=0; iRx<NRx; iRx++){
    	RX = &SC[Rx[iRx].SC];
    	TxNum = 0;

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
    	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
    	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);

    	FindVisibleTxMUOS(RX, &TxNum);

    	//SoOp
    	SpCntMUOS = 0;

    	if(TxNum > 0)
    	{
    		// SoOp
    		SoOp = &Rx[iRx].SoOp[TX_MUOS];
    		Sp = &SoOp->Sp[SpCntMUOS++];
    		TX = &SC[TxNum];

    		// ECI to ECEF of Transmitter's Position and Velocity
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);

    		FindSP_UD(RX, Sp, TX);

    		vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
    		vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
    		vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

    		/* ECEF to GRS80 */
    		ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

    		/* Update Specular Point Integrating DEM */
    		LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);

        	/* ECEF to GRS80 */
        	//ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);
    		/* ECEF to GRS80 */
    		//ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);
			/* GRS80 to ECEF */
			//GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

			// ECEF to ECI
			//MTxV(World[EARTH].CWN,Sp->PosW,Sp->PosN);

			// Coverage of CONUS
			if(Sp->latSUR<0.855211333477221 && Sp->latSUR>0.427164946867716
			&& Sp->lonSUR>-2.180570734210416 && Sp->lonSUR<-1.167593205124682){
				LL2EzRC9km(Sp->latSUR, Sp->lonSUR, &rowDy, &colDy);

				if(isnan(SurfaceDyn->RZSM[rowDy][colDy])){

				}
				else{
					if(flagRecording[iRx]==0 && recordingDoyLast[iRx] != doy){
						flagRecording[iRx] = 1;
						recordingDoy[iRx] = doy;
						recordingCnt[iRx] = 0;
					}
					if(flagRecording[iRx] && recordingDoy[iRx] == doy && recordingCnt[iRx] < 1200){
		    			fprintf(outFile1[iRx], "%lf %d %lf %lf %lf %ld %ld\n",
							AbsTime, TxNum, Sp->latSUR, Sp->lonSUR, Sp->altSUR, recordingDoy[iRx], recordingCnt[iRx]);
		    			recordingCnt[iRx]++;
					}
					else{
						flagRecording[iRx] = 0;
						recordingDoyLast[iRx] = recordingDoy[iRx];
					}
				}
			}

			// Snow cover
			rowDy = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
			colDy = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
			Sp->snowCover = SurfaceDyn->snowCover[rowDy][colDy];

			// Mean monthly snow 80% or more
			if(Sp->snowCover>79 && Sp->snowCover<101){
				if(flagRecording2[iRx]==0 && recordingDoyLast2[iRx] != doy){ // start recording
					flagRecording2[iRx] = 1;
					recordingDoy2[iRx] = doy;
					recordingCnt2[iRx] = 0;
				}
				if(flagRecording2[iRx] && recordingDoy2[iRx] == doy && recordingCnt2[iRx] < 1200){
					if(++snowArcDuration[iRx] > 29)
						validArc = 1;
					fprintf(outFile2[iRx], "%ld %lf %ld %lf %lf %lf %d %d %ld %ld %ld\n",
							snowArcID[iRx], AbsTime, TxNum, Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->snowCover, validArc, recordingDoy2[iRx], recordingCnt2[iRx], snowArcDuration[iRx]);
	    			recordingCnt2[iRx]++;
				}
				else{ // end of 1200 sec
					flagRecording2[iRx] = 0;
					recordingDoyLast2[iRx] = recordingDoy2[iRx];
					snowArcID[iRx]++;
					snowArcDuration[iRx] = 0;
				}
			}
			else{
				snowArcID[iRx]++;
				snowArcDuration[iRx] = 0;
			}
			//fprintf(outFile3[iRx], "%lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			//		AbsTime, TxNum, latT, lonT, altT, latR, lonR, altR, Sp->latSUR, Sp->lonSUR, Sp->altSUR);
    	}
//    	Rx[iRx].SoOp[TX_MUOS].SpCnt = SpCntMUOS;
    }
}

void Obs_CalVal_MUOS(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile1[5], *outFile2[5];//, *outFile3[5];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	struct SoilType *Soil;
	struct TargetAreaType *TARGET;
	double latR, lonR, altR;//, latT, lonT, altT;
	long TxNum;

	uint8_t iRx, iTarget, landMask;
	uint16_t rowSt, colSt, rowDy, colDy;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3], temp;
	double vNadir[3] = {0.0, 0.0, 1.0};

	static long *SpCntTotal;

    if (First) {
    	char filename[30];
    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(iRx=0; iRx<NRx; iRx++){
			sprintf(filename,"CalVal_MUOS_Case%d.42",iRx);
			outFile1[iRx] = FileOpen(InOutPath,filename,"w");
			sprintf(filename,"CONUS_MUOS_Case%d.42",iRx);
			outFile2[iRx] = FileOpen(InOutPath,filename,"w");
			//sprintf(filename,"POS_MUOS_Case%d.42",iRx);
			//outFile3[iRx] = FileOpen(InOutPath,filename,"w");
        }
    }
    for(iRx=0; iRx<NRx; iRx++){
    	RX = &SC[Rx[iRx].SC];
    	TxNum = 0;

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	//RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
    	//RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
    	//RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];
    	RX->PosRxAntN[0] = RX->PosN[0];
    	RX->PosRxAntN[1] = RX->PosN[1];
    	RX->PosRxAntN[2] = RX->PosN[2];

    	FindVisibleTxMUOS(RX, &TxNum);

    	if(TxNum > 0)
    	{
    		//printf("%lf %lf %lf\n", RX->B[ANTBODY].pn[0], RX->B[ANTBODY].pn[1], RX->B[ANTBODY].pn[2]);
    		// SoOp
    		SoOp = &Rx[iRx].SoOp[TX_MUOS];
    		Sp = &SoOp->Sp[0];
    		TX = &SC[TxNum];

           	/* ECI to ECEF of Receiver's Position */
        	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
        	/* ECEF to GRS80 */
        	ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);
    		// ECI to ECEF of Transmitter's Position
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
        	/* ECEF to GRS80 */
        	//ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);

    		vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
    		vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
    		vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

    		/* ECEF to GRS80 */
    		ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

    		/* Update Specular Point Integrating DEM */
    		LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);

			// Coverage of CONUS
			if(Sp->latSUR<0.855211333477221 && Sp->latSUR>0.427164946867716
			&& Sp->lonSUR>-2.180570734210416 && Sp->lonSUR<-1.167593205124682){
				fprintf(outFile2[iRx], "%d %lf %ld %lf %lf %lf\n",
						RX->Ncycle, AbsTime, TxNum, Sp->latSUR, Sp->lonSUR, Sp->altSUR);
			}

    		/* GRS80 to ECEF */
    		GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

	  		// ECEF to ECI
			MTxV(World[EARTH].CWN,Sp->PosW,Sp->PosN);

			/* Calculate Signal Angles of Arrival */
			vRx2TxN[0] = TX->PosN[0] - RX->PosRxAntN[0]; // for Sky-view antenna look angle & slant range
			vRx2TxN[1] = TX->PosN[1] - RX->PosRxAntN[1]; // for Sky-view antenna look angle & slant range
			vRx2TxN[2] = TX->PosN[2] - RX->PosRxAntN[2]; // for Sky-view antenna look angle & slant range

			Sp->PathLenD = UNITV(vRx2TxN);
			// Sky-view (Direct signal) antenna look angle
			//MxV(RX->B[ANTBODY].CN, vRx2TxN, vRx2TxA); // opposite direction vector of direct signal in Antenna body frame
			MxV(Orb[iRx].CLN, vRx2TxN, vRx2TxA); // opposite direction vector of reflected signal in Antenna body frame

			vRx2TxA[1] = -vRx2TxA[1]; // 180 degree rotation about x-axis
			vRx2TxA[2] = -vRx2TxA[2]; // 180 degree rotation about x-axis

			Sp->lookAngElD = acos(VoV(vRx2TxA,vNadir));
			Sp->lookAngAzD = atan2(vRx2TxA[1],vRx2TxA[0]);

    		// Check reflected signal angle of arrival
    		vRx2SpN[0] = Sp->PosN[0] - RX->PosRxAntN[0]; // for Earth-view antenna look angle
    		vRx2SpN[1] = Sp->PosN[1] - RX->PosRxAntN[1]; // for Earth-view antenna look angle
    		vRx2SpN[2] = Sp->PosN[2] - RX->PosRxAntN[2]; // for Earth-view antenna look angle
    		UNITV(vRx2SpN);
    		// Earth-view (Reflected signal) antenna look angle
    		//MxV(RX->B[ANTBODY].CN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame
    		MxV(Orb[iRx].CLN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame

    		Sp->lookAngElR = acos(VoV(vRx2SpA,vNadir));
   			Sp->lookAngAzR = atan2(vRx2SpA[1],vRx2SpA[0]);

			//if(Sp->lookAngElR*R2D>70)
			//	fprintf(outFile3[iRx], "%d %lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			//			RX->Ncycle, AbsTime, TxNum, latT, lonT, altT, latR, lonR, altR, Sp->latSUR, Sp->lonSUR, Sp->altSUR,
			//			Sp->lookAngAzD, Sp->lookAngElD, Sp->lookAngAzR, Sp->lookAngElR);
			//printf("R:(%lf, %lf)  ", Sp->lookAngElR*R2D, acos(VoV(vRx2SpN, test))*R2D);

   			// Antenna Gain
			calcAntPattern(SoOp, Sp);

			/* Calculate the First Fresnel Zone */
			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0]; // for Fresnel zone
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1]; // for Fresnel zone
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2]; // for Fresnel zone
			// ENU to AE
			temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
			vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
			vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
			vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

			temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
			Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
			Sp->elev = atan2(vSp2RxE[2],temp);
			Sp->incAng = HalfPi-Sp->elev;
			// First Fresnel zone
			Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
			Sp->a = Sp->b / sin(Sp->elev);
			Sp->area = Pi* Sp->a * Sp->b;

			// Direct Signal
			Sp->RxPwrD = Sp->AntGainS + Tx[TX_MUOS].EIRP + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926; //+ 10*log10(nSpotBeam);
			Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;

			// Reflected Signal
			Soil = Sp->Soil;
			R = Sp->Reflect;

			LL2EzRC1km(Sp->latSUR, Sp->lonSUR, &rowSt, &colSt);
			LL2EzRC9km(Sp->latSUR, Sp->lonSUR, &rowDy, &colDy);

			if(isnan(SurfaceStt->sand[rowSt][colSt]) || isnan(SurfaceStt->clay[rowSt][colSt]) || isnan(SurfaceStt->rho_b[rowSt][colSt])
			|| isnan(SurfaceDyn->RZSM[rowDy][colDy]) || isnan(SurfaceDyn->soilTemp1[rowDy][colDy]))
			{
				// do nothing
				Soil->m_v = NAN;
				Soil->T = NAN;
				Soil->sand = NAN;
				Soil->clay = NAN;
				Soil->rho_b = NAN;
				R->gammaLR = NAN;
				R->SnrR = NAN;
				landMask = 0;
			}
			else{
				calcDielPeplinski(rowSt, colSt, rowDy, colDy, Tx[TX_MUOS].Freq, Soil);
				calcReflectivity(Soil->e_c, R, Sp->incAng);

				R->RxPwrR = Sp->AntGainE + Tx[TX_MUOS].EIRP + 10*log10(R->gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenSUR) - 21.984197280441926; // + 10*log10(nSpotBeam * R->gammaLR) ;
				R->SnrR = R->RxPwrR - SoOp->NoisePwrE;
				landMask = 1;
			}

			// CalVal Overpass
			for(iTarget=0; iTarget<NTarget; iTarget++){
				TARGET = &Target[iTarget];
				if(FindVisibleSpotbeamTarget(RX, Sp, TARGET))
				{
					// do nothing
				}
				else{
					TARGET->vTarget2SpW[0] = Sp->PosW[0] - TARGET->PosW[0];
					TARGET->vTarget2SpW[1] = Sp->PosW[1] - TARGET->PosW[1];
					TARGET->vTarget2SpW[2] = Sp->PosW[2] - TARGET->PosW[2];

					TARGET->mTarget2Sp = MAGV(TARGET->vTarget2SpW);

					if(TARGET->mTarget2Sp > TARGET->radius)
						Sp->visibleTarget = 0;
					else
						Sp->visibleTarget = 1;

					// Land mask
/*					if(Sp->latSUR >= -1.483529864195180 && Sp->latSUR <= 1.483529864195180){
						// -85 < lat < 85
						//LL2EzRC1km(Sp->latSUR, Sp->lonSUR, &rowSt, &colSt);
						Sp->landMask = SurfaceStt->landMask[rowSt][colSt];
					}
					else{
						rowSt = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
						colSt = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
						Sp->landMask = SurfaceStt->landType[rowSt][colSt];
					}
*/
					fprintf(outFile1[iRx], "%d %lf %ld %d %d %lf %lf %lf %d %f %f %f %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
							RX->Ncycle, AbsTime, TxNum, iTarget, Sp->visibleTarget, Sp->latSUR, Sp->lonSUR, Sp->altSUR,
							landMask, Soil->m_v, Soil->T, Soil->sand, Soil->clay, Soil->rho_b,
							R->gammaLR, Sp->SnrD, R->SnrR,
							Sp->a, Sp->b, Sp->az,
							Sp->lookAngAzD, Sp->lookAngElD, Sp->lookAngAzR, Sp->lookAngElR);
				}
			}
    	}
    	//Rx[iRx].SoOp[TX_MUOS].SpCnt = SpCntMUOS;
    }
}

void Obs_Calval_Snow_CONUS(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile1[MAXNUM_TX][MAXNUM_RX];
    static FILE *outFile2[MAXNUM_TX][MAXNUM_RX];
    static FILE *outFile3[MAXNUM_TX][MAXNUM_RX];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	struct TargetAreaType *TARGET;
	double latT, lonT, altT, latR, lonR, altR;
	long ITx, IRx, ISp, iTarget;
	long *TxNum = NULL;
	uint8_t i;
	uint16_t row, col;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3];
	double vNadir[3] = {0.0, 0.0, 1.0};
	double temp;

	static long *SpCntTotal;
	long SpCntMUOS, SpCntOrbcomm, SpCntNOAA, SpCntMETEOR;
//	long SpCntGPS, SpCntGlonass, SpCntGalileo, SpCntBeidou;

    if (First) {
    	char filename[20], TxName[10];
    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(ITx=0; ITx<NTx; ITx++){
        	if(ITx == 0)
        		sprintf(TxName, "MUOS");
        	else if(ITx == 1)
        		sprintf(TxName, "Orbcomm");
        	else if(ITx == 2)
        		sprintf(TxName, "NOAA");
        	else if(ITx == 3)
        		sprintf(TxName, "METEOR");
        	else if(ITx == 4)
        		sprintf(TxName, "GPS");
        	else if(ITx == 5)
        		sprintf(TxName, "Glonass");
        	else if(ITx == 6)
        		sprintf(TxName, "Galileo");
        	else
        		sprintf(TxName, "Beidou");

            for(IRx=0; IRx<NRx; IRx++){
				// File for Cal/Val
				sprintf(filename,"%s_Rx%ld_calval.42", TxName, IRx);
				outFile1[ITx][IRx] = FileOpen(InOutPath,filename,"w");
				// File for Snow cover
				sprintf(filename,"%s_Rx%ld_SC.42", TxName, IRx);
				outFile2[ITx][IRx] = FileOpen(InOutPath,filename,"w");
				// File for Coverage of CONUS
				sprintf(filename,"%s_Rx%ld_CONUS.42", TxName, IRx);
				outFile3[ITx][IRx] = FileOpen(InOutPath,filename,"w");
            }
        }
    }
    for(IRx=0; IRx<NRx; IRx++){
    	RX = &SC[Rx[IRx].SC];

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
       	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
       	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

    	RX->VelRxAntN[0] = RX->VelN[0] + RX->B[ANTBODY].vn[0];
       	RX->VelRxAntN[1] = RX->VelN[1] + RX->B[ANTBODY].vn[1];
       	RX->VelRxAntN[2] = RX->VelN[2] + RX->B[ANTBODY].vn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
       	/* ECI to ECEF of Receiver's Velocity */
    	//MxV(World[EARTH].CWN,RX->VelRxAntN,RX->VelRxAntW);

    	/* ECEF to GRS80 */
    	ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);

    	//TxNum = FindVisibleTx(RX, 1, 64, TxNum, SpCntTotal); // max: 186
    	TxNum = FindVisibleTxSoOp(RX, TxNum, SpCntTotal);

    	//SoOp
    	SpCntMUOS = 0;
    	SpCntOrbcomm = 0;
    	SpCntNOAA = 0;
		SpCntMETEOR = 0;

    	for(ISp=0; ISp<*SpCntTotal; ISp++){
    		// GNSS
    		/*
    		if(TxNum[ISp]>123 && TxNum[ISp]<164){
    			SoOp = &Rx[IRx].SoOp[TX_BEIDOU];
				Sp = &SoOp->Sp[SpCntBeidou++];
				ITx = TX_BEIDOU;
    		}
    		else if(TxNum[ISp]>41 && TxNum[ISp]<73){
    			SoOp = &Rx[IRx].SoOp[TX_GPS];
    			Sp = &SoOp->Sp[SpCntGPS++];
    		    ITx = TX_GPS;
    		}
    		else if(TxNum[ISp]>72 && TxNum[ISp]<98){
				SoOp = &Rx[IRx].SoOp[TX_GLONASS];
				Sp = &SoOp->Sp[SpCntGlonass++];
				ITx = TX_GLONASS;
    		}
    		else if(TxNum[ISp]>97 && TxNum[ISp]<124){
    			SoOp = &Rx[IRx].SoOp[TX_GALILEO];
    			Sp = &SoOp->Sp[SpCntGalileo++];
    		    ITx = TX_GALILEO;
    		}*/
    		// SoOp
    		if(TxNum[ISp]<6){
    			SoOp = &Rx[IRx].SoOp[TX_MUOS];
    			Sp = &SoOp->Sp[SpCntMUOS++];
    		    ITx = TX_MUOS;
    		}
    		else if(TxNum[ISp]>5 && TxNum[ISp]<42){
    			SoOp = &Rx[IRx].SoOp[TX_ORBCOMM];
    			Sp = &SoOp->Sp[SpCntOrbcomm++];
    		    ITx = TX_ORBCOMM;
    		}
    		else if(TxNum[ISp]>41 && TxNum[ISp]<62){
    			SoOp = &Rx[IRx].SoOp[TX_NOAA];
    			Sp = &SoOp->Sp[SpCntNOAA++];
    		    ITx = TX_NOAA;
    		}
    		else{
    			SoOp = &Rx[IRx].SoOp[TX_METEOR];
    			Sp = &SoOp->Sp[SpCntMETEOR++];
    		    ITx = TX_METEOR;
    		}

    		TX = &SC[TxNum[ISp]];

    		// ECI to ECEF of Transmitter's Position and Velocity
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
    		//MxV(World[EARTH].CWN,TX->VelN,TX->VelW);

        	/* ECEF to GRS80 */
    		//ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);

			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0];
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1];
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2];

			/* ECEF to GRS80 */
			ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

			/* Update Specular Point Integrating DEM */
			LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);

			/* GRS80 to ECEF */
			GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

//			signalVec(TxNum[ISp]);
//			Sp->Exists = 1;
/*
			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0]; // for Fresnel zone
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1]; // for Fresnel zone
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2]; // for Fresnel zone

			// ENU to AE
			temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
			vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
			vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
			vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

			temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
			Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
			Sp->elev = atan2(vSp2RxE[2],temp);
//					Sp->incAng = HalfPi-Sp->elev;

			// First Fresnel zone
			Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
			Sp->a = Sp->b / sin(Sp->elev);
			Sp->area = Pi* Sp->a * Sp->b;

			// Coverage of CONUS
			if(Sp->latSUR<0.855211333477221 && Sp->latSUR>0.427164946867716
			&& Sp->lonSUR>-2.180570734210416 && Sp->lonSUR<-1.167593205124682){
				fprintf(outFile3[ITx][IRx][0], "%ld %lf %ld %lf %lf %lf %lf\n",
						RX->Ncycle, AbsTime, TxNum[ISp], Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->area);
			}

			// Snow cover
			row = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
			col = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
			Sp->snowCover = Surface->snowCover[row][col];
			// Mean monthly snow 80% or more
			if(Sp->snowCover>79 && Sp->snowCover<101)
				fprintf(outFile2[ITx][IRx][0], "%ld %lf %ld %lf %lf %lf %d %lf\n",
						RX->Ncycle, AbsTime, TxNum[ISp], Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->snowCover, Sp->area);
*/
			for(iTarget=0; iTarget<NTarget; iTarget++){
				TARGET = &Target[iTarget];
				if(FindVisibleSpotbeamTarget(RX, Sp, TARGET))
				{
					// do nothing
				}
				else{
					/* GRS80 to ECEF */
					//GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

					TARGET->vTarget2SpW[0] = Sp->PosW[0] - TARGET->PosW[0];
					TARGET->vTarget2SpW[1] = Sp->PosW[1] - TARGET->PosW[1];
					TARGET->vTarget2SpW[2] = Sp->PosW[2] - TARGET->PosW[2];

					TARGET->mTarget2Sp = MAGV(TARGET->vTarget2SpW);

					if(TARGET->mTarget2Sp > TARGET->radius)
						Sp->visibleTarget = 0;
					else
						Sp->visibleTarget = 1;

					/* ECEF to ECI */
//					MTxV(World[EARTH].CWN, Sp->PosW, Sp->PosN);

					// Land mask
					if(Sp->latSUR >= -1.483529864195180 && Sp->latSUR <= 1.483529864195180){
						// -85 < lat < 85
						LL2EzRC1km(Sp->latSUR, Sp->lonSUR, &row, &col);
						Sp->landMask = SurfaceStt->landMask[row][col];
					}
					else{
						row = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
						col = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
						Sp->landMask = SurfaceStt->landType[row][col];
					}

					fprintf(outFile1[ITx][IRx], "%d %lf %ld %ld %d %lf %lf %lf %lf %d\n",
							RX->Ncycle, AbsTime, TxNum[ISp], iTarget, Sp->visibleTarget, Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->area,
							Sp->landMask);
				}
			}
	    //	  printf("%lf %ld %ld %lf\n", AbsTime, Sp->landMask, Sp->visibleTarget, Orb[0].anom*R2D);
    	}
    	Rx[IRx].SoOp[TX_MUOS].SpCnt = SpCntMUOS;
    	Rx[IRx].SoOp[TX_ORBCOMM].SpCnt = SpCntOrbcomm;
    	Rx[IRx].SoOp[TX_NOAA].SpCnt = SpCntNOAA;
    	Rx[IRx].SoOp[TX_METEOR].SpCnt = SpCntMETEOR;

    	free(TxNum);
    }
}

void ObsOnTarget(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile[MAXNUM_TX][MAXNUM_RX];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	struct TargetAreaType *TARGET;
	double latT, lonT, altT, latR, lonR, altR;
	long ITx, IRx, ISp, iTarget;
	long *TxNum = NULL;
	uint8_t i;
	uint16_t row, col;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3];
	double vNadir[3] = {0.0, 0.0, 1.0};
	double temp;

	static long *SpCntTotal;
	long SpCntMUOS, SpCntORBCOMM, SpCntNOAA, SpCntMETEOR;
//	long SpCntGPS, SpCntGlonass, SpCntGalileo, SpCntBeidou;

    if (First) {
    	char filename[20], TxName[10];
    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(ITx=0; ITx<NTx; ITx++){
        	if(ITx == 0)
        		sprintf(TxName, "MUOS");
        	else if(ITx == 1)
        		sprintf(TxName, "ORBCOMM");
        	else if(ITx == 2)
        		sprintf(TxName, "NOAA");
        	else if(ITx == 3)
        		sprintf(TxName, "METEOR");
        	else if(ITx == 4)
        		sprintf(TxName, "GPS");
        	else if(ITx == 5)
        		sprintf(TxName, "Glonass");
        	else if(ITx == 6)
        		sprintf(TxName, "Galileo");
        	else
        		sprintf(TxName, "Beidou");

            for(IRx=0; IRx<NRx; IRx++){
				sprintf(filename,"%s_Rx%ld_calval.42", TxName, IRx);
				outFile[ITx][IRx] = FileOpen(InOutPath,filename,"w");
            }
        }
    }
    for(IRx=0; IRx<NRx; IRx++){
    	RX = &SC[Rx[IRx].SC];

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
       	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
       	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

    	RX->VelRxAntN[0] = RX->VelN[0] + RX->B[ANTBODY].vn[0];
       	RX->VelRxAntN[1] = RX->VelN[1] + RX->B[ANTBODY].vn[1];
       	RX->VelRxAntN[2] = RX->VelN[2] + RX->B[ANTBODY].vn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
       	/* ECI to ECEF of Receiver's Velocity */
    	//MxV(World[EARTH].CWN,RX->VelRxAntN,RX->VelRxAntW);

    	/* ECEF to GRS80 */
    	ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);

    	//TxNum = FindVisibleTx(RX, 1, 64, TxNum, SpCntTotal); // max: 186
    	TxNum = FindVisibleTxSoOp(RX, TxNum, SpCntTotal);

    	//SoOp
    	SpCntMUOS = 0;
    	SpCntORBCOMM = 0;
    	SpCntNOAA = 0;
		SpCntMETEOR = 0;
		// GNSS
//    	SpCntGPS = 0;
 //   	SpCntGlonass = 0;
  //  	SpCntGalileo = 0;
   // 	SpCntBeidou = 0;

    	for(ISp=0; ISp<*SpCntTotal; ISp++){
    		// GNSS
    		/*
    		if(TxNum[ISp]>123 && TxNum[ISp]<164){
    			SoOp = &Rx[IRx].SoOp[TX_BEIDOU];
				Sp = &SoOp->Sp[SpCntBeidou++];
				ITx = TX_BEIDOU;
    		}
    		else if(TxNum[ISp]>41 && TxNum[ISp]<73){
    			SoOp = &Rx[IRx].SoOp[TX_GPS];
    			Sp = &SoOp->Sp[SpCntGPS++];
    		    ITx = TX_GPS;
    		}
    		else if(TxNum[ISp]>72 && TxNum[ISp]<98){
				SoOp = &Rx[IRx].SoOp[TX_GLONASS];
				Sp = &SoOp->Sp[SpCntGlonass++];
				ITx = TX_GLONASS;
    		}
    		else if(TxNum[ISp]>97 && TxNum[ISp]<124){
    			SoOp = &Rx[IRx].SoOp[TX_GALILEO];
    			Sp = &SoOp->Sp[SpCntGalileo++];
    		    ITx = TX_GALILEO;
    		}*/
    		// SoOp
    		if(TxNum[ISp]<6){
    			SoOp = &Rx[IRx].SoOp[TX_MUOS];
    			Sp = &SoOp->Sp[SpCntMUOS++];
    		    ITx = TX_MUOS;
    		}
    		else if(TxNum[ISp]>5 && TxNum[ISp]<42){
    			SoOp = &Rx[IRx].SoOp[TX_ORBCOMM];
    			Sp = &SoOp->Sp[SpCntORBCOMM++];
    		    ITx = TX_ORBCOMM;
    		}
    		else if(TxNum[ISp]>41 && TxNum[ISp]<62){
    			SoOp = &Rx[IRx].SoOp[TX_NOAA];
    			Sp = &SoOp->Sp[SpCntNOAA++];
    		    ITx = TX_NOAA;
    		}
    		else{
    			SoOp = &Rx[IRx].SoOp[TX_METEOR];
    			Sp = &SoOp->Sp[SpCntMETEOR++];
    		    ITx = TX_METEOR;
    		}

    		TX = &SC[TxNum[ISp]];

    		// ECI to ECEF of Transmitter's Position and Velocity
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
    		//MxV(World[EARTH].CWN,TX->VelN,TX->VelW);

        	/* ECEF to GRS80 */
    		//ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);

			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0]; // for Fresnel zone
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1]; // for Fresnel zone
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2]; // for Fresnel zone

			/* ECEF to GRS80 */
			ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

			/* Update Specular Point Integrating DEM */
			LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);

			/* GRS80 to ECEF */
			GRS80ToECEF(Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->PosW);

			for(iTarget=0; iTarget<NTarget; iTarget++){
				TARGET = &Target[iTarget];
				if(FindVisibleSpotbeamTarget(RX, Sp, TARGET))
				{
					// do nothing
				}
				else{
					TARGET->vTarget2SpW[0] = Sp->PosW[0] - TARGET->PosW[0];
					TARGET->vTarget2SpW[1] = Sp->PosW[1] - TARGET->PosW[1];
					TARGET->vTarget2SpW[2] = Sp->PosW[2] - TARGET->PosW[2];

					TARGET->mTarget2Sp = MAGV(TARGET->vTarget2SpW);

					if(TARGET->mTarget2Sp > TARGET->radius)
						Sp->visibleTarget = 0;
					else
						Sp->visibleTarget = 1;

					/* ECEF to ECI */
//					MTxV(World[EARTH].CWN, Sp->PosW, Sp->PosN);

					for(i=0;i<3;i++) {
//						vRx2TxN[i] = TX->PosN[i] - RX->PosRxAntN[i]; // for Sky-view antenna look angle & slant range
//						vRx2SpN[i] = Sp->PosN[i] - RX->PosRxAntN[i]; // for Earth-view antenna look angle
						vSp2RxW[i] = RX->PosRxAntW[i] - Sp->PosW[i]; // for Fresnel zone
					}
/*					UNITV(vRx2SpN);

					Sp->PathLenD = UNITV(vRx2TxN);

					// Earth-view (Reflected signal) antenna look angle
					MxV(RX->B[ANTBODY].CN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame

					Sp->lookAngElR = acos(VoV(vRx2SpA,vNadir));
					Sp->lookAngAzR = atan2(vRx2SpA[1],vRx2SpA[0]);

					// Sky-view (Direct signal) antenna look angle
					MxV(RX->B[ANTBODY].CN, vRx2TxN, vRx2TxA); // opposite direction vector of direct signal in Antenna body frame

					vRx2TxA[1] = -vRx2TxA[1]; // 180 degree rotation about x-axis
					vRx2TxA[2] = -vRx2TxA[2]; // 180 degree rotation about x-axis

					Sp->lookAngElD = acos(VoV(vRx2TxA,vNadir));
					Sp->lookAngAzD = atan2(vRx2TxA[1],vRx2TxA[0]);

					// Antenna Gain
					calcAntPattern(SoOp, Sp);
*/
					// ENU to AE
					temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
					vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
					vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
					vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

					temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
					Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
					Sp->elev = atan2(vSp2RxE[2],temp);
//					Sp->incAng = HalfPi-Sp->elev;

					// First Fresnel zone
					Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
					Sp->a = Sp->b / sin(Sp->elev);
					Sp->area = Pi* Sp->a * Sp->b;
/*
					// Direct Signal
					Sp->RxPwrD = Sp->AntGainS + Tx[ITx].EIRP + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926;// -21.98419... = 10*log10(1/(4*pi)^2)
					Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;// - 10*log10(SoOp->CorIntTime);

					// Reflected Signal
					for(ISoil=0; ISoil<NSoil; ISoil++){
						R = &Sp->Reflect[ISoil];
						calcReflectivity(Soil[ITx][ISoil].e_c, R, Sp->incAng);

						R->RxPwrR = Sp->AntGainE + Tx[ITx].EIRP + 10*log10(R->gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenSUR) - 21.984197280441926;
						R->SnrR = R->RxPwrR - SoOp->NoisePwrE;// - 10*log10(SoOp->CorIntTime);

					//	fprintf(outFile[ITx][IRx][ISoil], "%lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					//			SimTime, iTarget, Sp->latSUR, Sp->lonSUR, Sp->a, Sp->b, Sp->az,
					//			Sp->AntGainS, Sp->AntGainE, Sp->PathLenD/Sp->PathLenSUR, Sp->SnrD, Sp->incAng,
					//			R->RcLR, R->gammaLR);
					}
*/
					if(Sp->latSUR >= -1.483529864195180 && Sp->latSUR <= 1.483529864195180){
						// -85 < lat < 85
						// Land mask
						LL2EzRC1km(Sp->latSUR, Sp->lonSUR, &row, &col);
						Sp->landMask = SurfaceStt->landMask[row][col];

						// Snow cover
						row = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
						col = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
						Sp->snowCover = SurfaceDyn->snowCover[row][col];
					}
					else{
						// Land mask
						row = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
						col = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
						Sp->landMask = SurfaceStt->landType[row][col];

						// Snow cover
						Sp->snowCover = SurfaceDyn->snowCover[row][col];
					}

					fprintf(outFile[ITx][IRx], "%d %lf %ld %ld %d %lf %lf %lf %lf %d %d\n",
							RX->Ncycle, AbsTime, TxNum[ISp], iTarget, Sp->visibleTarget, Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->area,
							Sp->landMask, Sp->snowCover);
				}
			}
    	}
    	free(TxNum);
    }
}

void ObsGlobe(void)
{
    static uint8_t First = TRUE;
    static FILE *outFile1[MAXNUM_TX][MAXNUM_RX];
    static FILE *outFile2[MAXNUM_TX][MAXNUM_RX];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	struct ReflectType *R;
	double latT, lonT, altT, latR, lonR, altR;
	long ITx, IRx, ISp;
	long *TxNum = NULL;
	uint8_t i, nSpotBeam;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3];
	double vNadir[3] = {0.0, 0.0, 1.0};
	double temp;

	static long *SpCntTotal;
	long SpCntMUOS, SpCntORBCOMM, SpCntGPS, SpCntGlonass, SpCntGalileo, SpCntBeidou, SpCntNOAA, SpCntMETEOR;
	double dist[3], dist2[3];
	uint16_t row, col;
	double rE, sinLat2;
    double a = 6378137.0;
//    double f = 1.0/298.257223563; // WGS84
    double f = 1.0/298.257222101; // GRS80
    double e2 = f*(2.0-f);
//	double RxPosN[3], RxPosW[3], SpPosW[3];
//	double TxPosN[3], TxPosW[3];
//	double vNorm[3];
//	double err1, theta_inc, theta_ref, theta_diff;

    if (First) {
    	char filename[20], TxName[10];

    	SpCntTotal = (long *) calloc(1, sizeof(long));
        First = FALSE;
        for(ITx=0; ITx<NTx; ITx++){
        	if(ITx == 0)
        		sprintf(TxName, "MUOS");
        	else if(ITx == 1)
        		sprintf(TxName, "ORBCOMM");
        	else if(ITx == 2)
        		sprintf(TxName, "GPS");
        	else if(ITx == 3)
        		sprintf(TxName, "Glonass");
        	else if(ITx == 4)
        		sprintf(TxName, "Galileo");
        	else if(ITx == 5)
        		sprintf(TxName, "Beidou");
        	else if(ITx == 6)
        		sprintf(TxName, "NOAA");
        	else
        		sprintf(TxName, "METEOR");
            for(IRx=0; IRx<NRx; IRx++){
				// File for Snow cover
				sprintf(filename,"%s_Rx%ld_SC.42", TxName, IRx);
				outFile1[ITx][IRx] = FileOpen(InOutPath,filename,"w");
				// File for Coverage of CONUS
				sprintf(filename,"%s_Rx%ld_CONUS.42", TxName, IRx);
				outFile2[ITx][IRx] = FileOpen(InOutPath,filename,"w");
            }
        }
    }
    for(IRx=0; IRx<NRx; IRx++){
    	RX = &SC[Rx[IRx].SC];

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
       	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
       	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

    	RX->VelRxAntN[0] = RX->VelN[0] + RX->B[ANTBODY].vn[0];
       	RX->VelRxAntN[1] = RX->VelN[1] + RX->B[ANTBODY].vn[1];
       	RX->VelRxAntN[2] = RX->VelN[2] + RX->B[ANTBODY].vn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
       	/* ECI to ECEF of Receiver's Velocity */
    	MxV(World[EARTH].CWN,RX->VelRxAntN,RX->VelRxAntW);

    	/* ECEF to GRS80 */
    	ECEFToGRS80(RX->PosRxAntW, &latR, &lonR, &altR);

    	TxNum = FindVisibleTx(RX, 1, 5, TxNum, SpCntTotal);

    	SpCntMUOS = 0;
    	SpCntORBCOMM = 0;
    	SpCntGPS = 0;
    	SpCntGlonass = 0;
    	SpCntGalileo = 0;
    	SpCntBeidou = 0;
    	SpCntNOAA = 0;
    	SpCntMETEOR = 0;

    	for(ISp=0; ISp<*SpCntTotal; ISp++){
    		if(TxNum[ISp]>123 && TxNum[ISp]<164){
    			SoOp = &Rx[IRx].SoOp[TX_BEIDOU];
				Sp = &SoOp->Sp[SpCntBeidou++];
				ITx = TX_BEIDOU;
    		}
    		else if(TxNum[ISp]>41 && TxNum[ISp]<73){
    			SoOp = &Rx[IRx].SoOp[TX_GPS];
    			Sp = &SoOp->Sp[SpCntGPS++];
    		    ITx = TX_GPS;
    		}
    		else if(TxNum[ISp]>72 && TxNum[ISp]<98){
				SoOp = &Rx[IRx].SoOp[TX_GLONASS];
				Sp = &SoOp->Sp[SpCntGlonass++];
				ITx = TX_GLONASS;
    		}
    		else if(TxNum[ISp]>97 && TxNum[ISp]<124){
    			SoOp = &Rx[IRx].SoOp[TX_GALILEO];
    			Sp = &SoOp->Sp[SpCntGalileo++];
    		    ITx = TX_GALILEO;
    		}
    		else if(TxNum[ISp]<6){
    			SoOp = &Rx[IRx].SoOp[TX_MUOS];
    			Sp = &SoOp->Sp[SpCntMUOS++];
    		    ITx = TX_MUOS;
    		}
    		else if(TxNum[ISp]>163 && TxNum[ISp]<184){
    			SoOp = &Rx[IRx].SoOp[TX_NOAA];
    			Sp = &SoOp->Sp[SpCntNOAA++];
    		    ITx = TX_NOAA;
    		}
    		else if(TxNum[ISp]>183){
    			SoOp = &Rx[IRx].SoOp[TX_METEOR];
    			Sp = &SoOp->Sp[SpCntMETEOR++];
    		    ITx = TX_METEOR;
    		}
    		else{
    			SoOp = &Rx[IRx].SoOp[TX_ORBCOMM];
    			Sp = &SoOp->Sp[SpCntORBCOMM++];
    		    ITx = TX_ORBCOMM;
    		}

    		Sp->SpNum = ISp;
    		TX = &SC[TxNum[ISp]];

    		//Sp->Exists = 0;

    		// ECI to ECEF of Transmitter's Position and Velocity
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
    		MxV(World[EARTH].CWN,TX->VelN,TX->VelW);

        	/* ECEF to GRS80 */
    		ECEFToGRS80(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);

    	//	Sp->Exists = 1;

    		//if(Sp->Exists){

			vSp2RxW[0] = RX->PosRxAntW[0] - Sp->PosW[0]; // for Fresnel zone
			vSp2RxW[1] = RX->PosRxAntW[1] - Sp->PosW[1]; // for Fresnel zone
			vSp2RxW[2] = RX->PosRxAntW[2] - Sp->PosW[2]; // for Fresnel zone

			/* ECEF to GRS80 */
			ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);

			/* Update Specular Point Integrating DEM */
			LocalGeometryDEM(vSp2RxW, RX, SoOp, Sp, TX);
			fprintf(outFile2[ITx][IRx], "%d %lf %ld %lf %lf %lf\n",
					RX->Ncycle, AbsTime, TxNum[ISp], Sp->latSUR, Sp->lonSUR, Sp->altSUR);
/*
			// Coverage of CONUS
			if(Sp->latSUR<0.855211333477221 && Sp->latSUR>0.427164946867716
			&& Sp->lonSUR>-2.180570734210416 && Sp->lonSUR<-1.167593205124682){
				fprintf(outFile2[ITx][IRx], "%d %lf %ld %lf %lf %lf\n",
						RX->Ncycle, AbsTime, TxNum[ISp], Sp->latSUR, Sp->lonSUR, Sp->altSUR);
			}
			*/
			/* ECEF to ECI */
			//MTxV(World[EARTH].CWN, Sp->PosW, Sp->PosN);
/*
			for(i=0;i<3;i++) {
				vRx2TxN[i] = TX->PosN[i] - RX->PosRxAntN[i]; // for Sky-view antenna look angle & slant range
				vRx2SpN[i] = Sp->PosN[i] - RX->PosRxAntN[i]; // for Earth-view antenna look angle
				vSp2RxW[i] = RX->PosRxAntW[i] - Sp->PosW[i]; // for Fresnel zone
			}
			UNITV(vRx2SpN);

			Sp->PathLenD = UNITV(vRx2TxN);

			// Earth-view (Reflected signal) antenna look angle
			MxV(RX->B[ANTBODY].CN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame

			Sp->lookAngElR = acos(VoV(vRx2SpA,vNadir));
			Sp->lookAngAzR = atan2(vRx2SpA[1],vRx2SpA[0]);

			// Sky-view (Direct signal) antenna look angle
			MxV(RX->B[ANTBODY].CN, vRx2TxN, vRx2TxA); // opposite direction vector of direct signal in Antenna body frame

			vRx2TxA[1] = -vRx2TxA[1]; // 180 degree rotation about x-axis
			vRx2TxA[2] = -vRx2TxA[2]; // 180 degree rotation about x-axis

			Sp->lookAngElD = acos(VoV(vRx2TxA,vNadir));
			Sp->lookAngAzD = atan2(vRx2TxA[1],vRx2TxA[0]);

			// Antenna Gain
			calcAntPattern(SoOp, Sp);

			// ENU to AE
			temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
			vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
			vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
			vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

			temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
			Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
			Sp->elev = atan2(vSp2RxE[2],temp);
			Sp->incAng = HalfPi-Sp->elev;

			// First Fresnel zone
			Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
			Sp->a = Sp->b / sin(Sp->elev);

			// Direct Signal
			if(ITx != TX_MUOS){
				Sp->RxPwrD = Sp->AntGainS + Tx[ITx].EIRP + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926;// -21.98419... = 10*log10(1/(4*pi)^2)

			}
			else{
				nSpotBeam = FindVisibleSpotbeam(RX, Sp, TX);
				Sp->RxPwrD = Sp->AntGainS + Tx[ITx].EIRP + 10*log10(nSpotBeam) + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926;
			}
			Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;// - 10*log10(SoOp->CorIntTime);

			// Reflected Signal
*//*		for(ISoil=0; ISoil<NSoil; ISoil++){
			R = &Sp->Reflect[ISoil];
			calcReflectivity(Soil[ITx][ISoil].e_c, R, Sp->incAng);

			if(ITx != TX_MUOS){
				R->RxPwrR = Sp->AntGainE + Tx[ITx].EIRP + 10*log10(R->gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenSUR) - 21.984197280441926;
			}
			else{
				R->RxPwrR = Sp->AntGainE + Tx[ITx].EIRP + 10*log10(nSpotBeam * R->gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenSUR) - 21.984197280441926;
			}
			R->SnrR = R->RxPwrR - SoOp->NoisePwrE;// - 10*log10(SoOp->CorIntTime);
*/
			//fprintf(outFile[ITx][IRx][ISoil], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			//		SimTime, Sp->latSUR, Sp->lonSUR, Sp->a, Sp->b, Sp->az,
			//		Sp->AntGainS, Sp->AntGainE, Sp->PathLenD/Sp->PathLenSUR, Sp->SnrD, Sp->incAng,
			//		R->RcLR, R->gammaLR);

/*			// DEM
			row = Sp->latSUR*R2Dx60+5399.9999;
			col = Sp->lonSUR*R2Dx60+10799.9999;
			sinLat2 = sin(Sp->latSUR)*sin(Sp->latSUR);
			rE = a*sqrt((1-e2*(2-e2)*sinLat2) / (1-e2*sinLat2));
*/
/*
			// Snow cover
			row = 1799.9999 - Sp->latSUR*R2Dx20;//(90-lat)*20;
			col = 3599.9999 + Sp->lonSUR*R2Dx20;//(lon+180)*20;
			Sp->snowCover = SurfaceDyn->snowCover[row][col];
			// Mean monthly snow 80% or more
			if(Sp->snowCover>79 && Sp->snowCover<101)
				fprintf(outFile1[ITx][IRx], "%d %lf %ld %lf %lf %lf %d\n",
						RX->Ncycle, AbsTime, TxNum[ISp], Sp->latSUR, Sp->lonSUR, Sp->altSUR, Sp->snowCover);
						*/
    	}
    	Rx[IRx].SoOp[TX_MUOS].SpCnt = SpCntMUOS;
    	Rx[IRx].SoOp[TX_ORBCOMM].SpCnt = SpCntORBCOMM;
    	Rx[IRx].SoOp[TX_GPS].SpCnt = SpCntGPS;
    	Rx[IRx].SoOp[TX_GLONASS].SpCnt = SpCntGlonass;
    	Rx[IRx].SoOp[TX_GALILEO].SpCnt = SpCntGalileo;
    	Rx[IRx].SoOp[TX_BEIDOU].SpCnt = SpCntBeidou;
    	Rx[IRx].SoOp[TX_NOAA].SpCnt = SpCntNOAA;
    	Rx[IRx].SoOp[TX_METEOR].SpCnt = SpCntMETEOR;

    	free(TxNum);
    }
}

void Reflectometry_MUOS(void)
{
    static char First = TRUE;
    static FILE *outFile[MAXNUM_RX];
    char filename[10];

    struct SCType *RX;
	struct SoOpType *SoOp;
	struct SpecularType *Sp;
	struct SCType *TX;
	double latT, lonT, altT, latR, lonR, altR;
	long IRx, ISp;
	long *TxNum = NULL;
	uint8_t i, nSpotBeam;

	double vRx2TxN[3], vRx2TxA[3], vRx2SpN[3], vRx2SpA[3], vSp2RxW[3], vSp2RxE[3], vSp2TxN[3];
	double vNadir[3] = {0.0, 0.0, 1.0};
	double mRx2Sp, mSp2Tx, temp;

	static long cnt = 0;

    if (First) {
        First = FALSE;
        for (IRx=0; IRx<NRx; IRx++){
            sprintf(filename,"MUOS%02ld.42",IRx);
            outFile[IRx] = FileOpen(InOutPath,filename,"w");
        }
    }
    for(IRx=0; IRx<NRx; IRx++){
    	SoOp = &Rx[IRx].SoOp[0];
    	RX = &SC[Rx[IRx].SC];

    	/* Position of Receiving Antenna wrt Rx's ref pt in ECI*/
    	RX->PosRxAntN[0] = RX->PosN[0] + RX->B[ANTBODY].pn[0];
       	RX->PosRxAntN[1] = RX->PosN[1] + RX->B[ANTBODY].pn[1];
       	RX->PosRxAntN[2] = RX->PosN[2] + RX->B[ANTBODY].pn[2];

    	RX->VelRxAntN[0] = RX->VelN[0] + RX->B[ANTBODY].vn[0];
       	RX->VelRxAntN[1] = RX->VelN[1] + RX->B[ANTBODY].vn[1];
       	RX->VelRxAntN[2] = RX->VelN[2] + RX->B[ANTBODY].vn[2];

       	/* ECI to ECEF of Receiver's Position */
    	MxV(World[EARTH].CWN,RX->PosRxAntN,RX->PosRxAntW);
       	/* ECI to ECEF of Receiver's Velocity */
    	MxV(World[EARTH].CWN,RX->VelRxAntN,RX->VelRxAntW);

    	/* ECEF to WGS84 */
    	ECEFToWGS84(RX->PosRxAntW, &latR, &lonR, &altR);

    //AbsTime>689775289
    //	if(AbsTime>690447223  || AbsTime<683726400.4) {
    //		printf("%lf %lf %lf\n", AbsTime, lat, lon);
    	//	fprintf(RxFile, "%lf %lf %lf\n", AbsTime, lat, lon);
    //	}
    	TxNum = FindVisibleTx(RX, 1, 4, TxNum, &SoOp->SpCnt);

    	for(ISp=0; ISp<SoOp->SpCnt; ISp++){
    		Sp = &SoOp->Sp[ISp];
    		Sp->SpNum = ISp;
    		TX = &SC[TxNum[ISp]];

    		Sp->Exists = 0;

    		// ECI to ECEF of Transmitter's Position
    		MxV(World[EARTH].CWN,TX->PosN,TX->PosW);
        	/* ECEF to WGS84 */
        	ECEFToWGS84(TX->PosW, &latT, &lonT, &altT);

    		FindSP_UD(RX, Sp, TX);
    		//FindSP_MPL(RX, Sp, TX);
    		nSpotBeam = FindVisibleSpotbeam(RX, Sp, TX);
    		//nSpotBeam = FindSP_MPL_Spotbeam(RX, Sp, TX);

    		/* ECEF to ECI */
    		MTxV(World[EARTH].CWN, Sp->PosW, Sp->PosN);

    		/* ECEF to GRS80 */
    		ECEFToGRS80(Sp->PosW, &Sp->latGRS, &Sp->lonGRS, &Sp->altGRS);


    		//printf("%lf %lf %le\n", Sp->lat*R2D, Sp->lon*R2D, Sp->alt);

    		if(Sp->Exists){
    			for(i=0;i<3;i++) {
    				vRx2TxN[i] = TX->PosN[i] - RX->PosRxAntN[i]; // for Sky-view antenna look angle & slant range
    				vRx2SpN[i] = Sp->PosN[i] - RX->PosRxAntN[i]; // for Earth-view antenna look angle
    				vSp2TxN[i] = TX->PosN[i] - Sp->PosN[i];
    				vSp2RxW[i] = RX->PosRxAntW[i] - Sp->PosW[i]; // for Fresnel zone
    			}
    			mRx2Sp = UNITV(vRx2SpN);
    			mSp2Tx = UNITV(vSp2TxN);

    			Sp->PathLenD = UNITV(vRx2TxN);
    			Sp->PathLenSUR = mRx2Sp + mSp2Tx;

    			// Earth-view (Reflected signal) antenna look angle
    			MxV(RX->B[ANTBODY].CN, vRx2SpN, vRx2SpA); // opposite direction vector of reflected signal in Antenna body frame

       			Sp->lookAngElR = acos(VoV(vRx2SpA,vNadir));
        		Sp->lookAngAzR = atan2(vRx2SpA[1],vRx2SpA[0]);

    			// Sky-view (Direct signal) antenna look angle
    			MxV(RX->B[ANTBODY].CN, vRx2TxN, vRx2TxA); // opposite direction vector of direct signal in Antenna body frame

        		vRx2TxA[1] = -vRx2TxA[1]; // 180 degree rotation about x-axis
    			vRx2TxA[2] = -vRx2TxA[2]; // 180 degree rotation about x-axis

    			Sp->lookAngElD = acos(VoV(vRx2TxA,vNadir));
    			Sp->lookAngAzD = atan2(vRx2TxA[1],vRx2TxA[0]);

    			// ECEF(dx,dy,dz From Tx to Sp) to ENU wrt (lat, lon) of Sp
    			//Sp->CWE[0][0] = -sin(Sp->lon);
    			//Sp->CWE[0][1] = cos(Sp->lon);
    			//Sp->CWE[0][2] = 0;
    			//Sp->CWE[1][0] = -cos(Sp->lon)*sin(Sp->lat);
    			//Sp->CWE[1][1] = -sin(Sp->lon)*sin(Sp->lat);
    			//Sp->CWE[1][2] = cos(Sp->lat);
    			//Sp->CWE[2][0] = cos(Sp->lon)*cos(Sp->lat);
    			//Sp->CWE[2][1] = sin(Sp->lon)*cos(Sp->lat);
    			//Sp->CWE[2][2] = sin(Sp->lat);

    			//MxV(Sp->CWE,vSp2RxW,vSp2RxE);
    			temp = cos(Sp->lonSUR)*vSp2RxW[0] + sin(Sp->lonSUR)*vSp2RxW[1];
    			vSp2RxE[0] = cos(Sp->lonSUR)*vSp2RxW[1] - sin(Sp->lonSUR)*vSp2RxW[0];
    			vSp2RxE[1] = cos(Sp->latSUR)*vSp2RxW[2] - sin(Sp->latSUR)*temp;
    			vSp2RxE[2] = cos(Sp->latSUR)*temp + sin(Sp->latSUR)*vSp2RxW[2];

    			// ENU to AE
    			temp = sqrt(vSp2RxE[0]*vSp2RxE[0] + vSp2RxE[1]*vSp2RxE[1]);
    			Sp->az = atan2(vSp2RxE[0],vSp2RxE[1]);
    			Sp->elev = atan2(vSp2RxE[2],temp);
    			Sp->incAng = HalfPi-Sp->elev;
//    			Sp->incAng = 0.5* acos((MagVecSp2Tx*MagVecSp2Tx + MagVecSp2Rx*MagVecSp2Rx - Sp->PathLenD*Sp->PathLenD)/(2*MagVecSp2Tx*MagVecSp2Rx));
//    			Sp->Reflect.gammaLR = 0.1;
    			//calcReflectivity(Soil[0][0].e_c, Sp);

    			Sp->b = sqrt(4*SoOp->Wavelength*altR*sin(Sp->elev)+SoOp->Wavelength*SoOp->Wavelength) / sin(Sp->elev)/2;
    			Sp->a = Sp->b / sin(Sp->elev);
    			//temp = sqrt(SoOp->Wavelength*MagVecSp2Rx*MagVecSp2Tx/Sp->PathLenR)/sin(Sp->elev); // Sgphaladi thesis

/*    			double ut1[3], ut2[3], ur1[3], ur2[3];
    			double u_tr[2][2], U_tr[4][4];

    			tanUnitVec(TX->B[0].CN, RX->B[ANTBODY].CN, VecTx2Rx, Tx[0].Pol, Rx[0].Pol, ut1, ut2, ur1, ur2);

    			u_tr[0][0] = VoV(ut1, conj(ur1));
    			u_tr[0][1] = VoV(ut2, conj(ur1));
    			u_tr[1][0] = VoV(ut1, conj(ur2));
    			u_tr[1][1] = VoV(ut2, conj(ur2));

    			calcMueller(u_tr, U_tr);*/

    			// -21.98419... = 10*log10(1/(4*pi)^2)

    			calcAntPattern(SoOp, Sp);

    			Sp->RxPwrD = Sp->AntGainS + Tx[0].EIRP + 10*log10(nSpotBeam) + 20*log10(SoOp->Wavelength/Sp->PathLenD) - 21.984197280441926;
//    			Sp->RxPwrR = Sp->AntGainE + Tx[0].EIRP + 10*log10(nSpotBeam * Sp->Reflect.gammaLR) + 20*log10(SoOp->Wavelength/Sp->PathLenR) - 21.984197280441926;

    			Sp->SnrD = Sp->RxPwrD - SoOp->NoisePwrS;// - 10*log10(SoOp->CorIntTime);
 //   			Sp->SnrR = Sp->RxPwrR - SoOp->NoisePwrE;// - 10*log10(SoOp->CorIntTime);

    			cnt++;
    	/*		if(First2 && cnt>1000){
    			    FILE *spFile;
    				long j;
    				First2 = FALSE;
    	            sprintf(filename,"LLA_sp.42");
    	            spFile = FileOpen(InOutPath,filename,"w");
        			LocalGeometry(vSp2RxW, RX, SoOp, Sp, TX);

					fprintf(spFile, "%lf %lf %lf\n", Sp->lat, Sp->lon, Sp->alt);
					ReportLLA();

        			for(i=0; i<SoOp->NthetaLocal; i++){
        				for(j=0; j<SoOp->NphiLocal; j++){
        	       			fprintf(outFile[IRx], "%lf %lf %lf %le %le %lf %lf %lf %ld %ld %lf %lf %lf\n",
        	       					Sp->r, Sp->Local[i][j].theta, Sp->Local[i][j].phi, Sp->Local[i][j].delay, Sp->Local[i][j].doppler,
									Sp->Local[i][j].lat, Sp->Local[i][j].lon, Sp->Local[i][j].alt,
									SoOp->NthetaLocal, SoOp->NphiLocal,
									Sp->a, Sp->b, Sp->az);
        				}
        			}
    			}*/
    			double RxPosW[3];
				double TxPosW[3];
				double v1[3];
				double err1;

				CopyUnitV(RX->PosRxAntW, RxPosW);
				CopyUnitV(TX->PosW, TxPosW);
				VxV(RxPosW, TxPosW, v1);
				err1 = VoV(v1, Sp->PosW)/1000;



    			//printf("%lf %lf\n", err1, err2);
    			fprintf(outFile[IRx], "%lf %lf %lf %lf\n", Sp->latSUR, Sp->lonSUR, err1, Sp->incAng);
    			//fprintf(outFile[IRx],"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
    			//		SimTime, latR, lonR, latT, lonT, Sp->lat, Sp->lon, Sp->SnrD, Sp->SnrR, Sp->Reflect.RcLR, Sp->incAng,
				//		SoOp->AntGainS, SoOp->AntGainE, Sp->lookAngAzD, Sp->lookAngElD, Sp->lookAngAzR, Sp->lookAngElR,
				//		Sp->a, Sp->b, Sp->az);
    		}
    	}
    	free(TxNum);
    }
}

void ReportLLA_Rx(void)
{
	static FILE *outFile[5];
	static char First = TRUE;
	double lat[5], lon[5], alt[5];
	uint8_t iRx;

	if (First) {
		First = FALSE;
    	char filename[20];
    	for(iRx=0; iRx<NRx; iRx++){
			sprintf(filename,"Rx%d_LLA.42",iRx);
			outFile[iRx] = FileOpen(InOutPath,filename,"w");
        }
	}

	for(iRx=0;iRx<NRx;iRx++){
		/* ECI to ECEF of Receiver's Position */
		MxV(World[EARTH].CWN,SC[iRx].PosN,SC[iRx].PosW);
		/* ECEF to WGS84 */
		ECEFToWGS84(SC[iRx].PosW, &lat[iRx], &lon[iRx], &alt[iRx]);

		fprintf(outFile[iRx], "%lf %lf %lf %lf\n", SimTime, lat[iRx]*R2D, lon[iRx]*R2D, alt[iRx]);
	}
}

void ReportLLA(void)
{
	static FILE *MUOS1file, *MUOS2file, *MUOS3file, *MUOS4file, *RxFile;
	static char First = TRUE;
	double lat[5], lon[5], alt[5];
	long i;

	if (First) {
		First = FALSE;
		RxFile = FileOpen(InOutPath, "LLA_Rx.42","w");
		MUOS1file = FileOpen(InOutPath,"LLA_MUOS1.42","w");
		MUOS2file = FileOpen(InOutPath,"LLA_MUOS2.42","w");
		MUOS3file = FileOpen(InOutPath,"LLA_MUOS3.42","w");
		MUOS4file = FileOpen(InOutPath,"LLA_MUOS4.42","w");
		fprintf(RxFile, "SimTime lat[deg] lon[deg] alt[m]\n");
		fprintf(MUOS1file, "SimTime lat[deg] lon[deg] alt[m]\n");
		fprintf(MUOS2file, "SimTime lat[deg] lon[deg] alt[m]\n");
		fprintf(MUOS3file, "SimTime lat[deg] lon[deg] alt[m]\n");
		fprintf(MUOS4file, "SimTime lat[deg] lon[deg] alt[m]\n");
	}

	if(OutFlag){
		for(i=0;i<5;i++){
	    	/* ECI to ECEF of Receiver's Position */
	    	MxV(World[EARTH].CWN,SC[i].PosN,SC[i].PosW);
	    	/* ECEF to WGS84 */
	    	ECEFToWGS84(SC[i].PosW, &lat[i], &lon[i], &alt[i]);
		}

		fprintf(RxFile, "%lf %lf %lf %lf\n", SimTime, lat[0]*R2D, lon[0]*R2D, alt[0]);
		fprintf(MUOS1file, "%lf %lf %lf %lf\n", SimTime, lat[1]*R2D, lon[1]*R2D, alt[1]);
		fprintf(MUOS2file, "%lf %lf %lf %lf\n", SimTime, lat[2]*R2D, lon[2]*R2D, alt[2]);
		fprintf(MUOS3file, "%lf %lf %lf %lf\n", SimTime, lat[3]*R2D, lon[3]*R2D, alt[3]);
		fprintf(MUOS4file, "%lf %lf %lf %lf\n", SimTime, lat[4]*R2D, lon[4]*R2D, alt[4]);
	}
}

void ReportPosN(void)
{
	static FILE *MUOS1file, *MUOS2file, *MUOS3file, *MUOS4file, *RxFile;
	static char First = TRUE;

	if (First) {
		First = FALSE;
		RxFile = FileOpen(InOutPath, "PosN_Rx.42","w");
		MUOS1file = FileOpen(InOutPath,"PosN_MUOS1.42","w");
		MUOS2file = FileOpen(InOutPath,"PosN_MUOS2.42","w");
		MUOS3file = FileOpen(InOutPath,"PosN_MUOS3.42","w");
		MUOS4file = FileOpen(InOutPath,"PosN_MUOS4.42","w");
	}

	if(OutFlag){
		fprintf(RxFile, "%lf %lf %lf %lf\n", SimTime, SC[0].PosN[0]/1000, SC[0].PosN[1]/1000, SC[0].PosN[2]/1000);
		fprintf(MUOS1file, "%lf %lf %lf %lf\n", SimTime, SC[1].PosN[0]/1000, SC[1].PosN[1]/1000, SC[1].PosN[2]/1000);
		fprintf(MUOS2file, "%lf %lf %lf %lf\n", SimTime, SC[2].PosN[0]/1000, SC[2].PosN[1]/1000, SC[2].PosN[2]/1000);
		fprintf(MUOS3file, "%lf %lf %lf %lf\n", SimTime, SC[3].PosN[0]/1000, SC[3].PosN[1]/1000, SC[3].PosN[2]/1000);
		fprintf(MUOS4file, "%lf %lf %lf %lf\n", SimTime, SC[4].PosN[0]/1000, SC[4].PosN[1]/1000, SC[4].PosN[2]/1000);
	}
}

void ReportPosW(void)
{
	static FILE *MUOS1file, *MUOS2file, *MUOS3file, *MUOS4file, *RxFile;
	static char First = TRUE;
	long i;

	if (First) {
		First = FALSE;
		RxFile = FileOpen(InOutPath, "PosW_Rx.42","w");
		MUOS1file = FileOpen(InOutPath,"PosW_MUOS1.42","w");
		MUOS2file = FileOpen(InOutPath,"PosW_MUOS2.42","w");
		MUOS3file = FileOpen(InOutPath,"PosW_MUOS3.42","w");
		MUOS4file = FileOpen(InOutPath,"PosW_MUOS4.42","w");
	}

	if(OutFlag){
		for(i=0;i<5;i++){
	    	/* ECI to ECEF of Receiver's Position */
	    	MxV(World[EARTH].CWN,SC[i].PosN,SC[i].PosW);
		}
		fprintf(RxFile, "%lf %lf %lf %lf\n", SimTime, SC[0].PosW[0]/1000, SC[0].PosW[1]/1000, SC[0].PosW[2]/1000);
		fprintf(MUOS1file, "%lf %lf %lf %lf\n", SimTime, SC[1].PosW[0]/1000, SC[1].PosW[1]/1000, SC[1].PosW[2]/1000);
		fprintf(MUOS2file, "%lf %lf %lf %lf\n", SimTime, SC[2].PosW[0]/1000, SC[2].PosW[1]/1000, SC[2].PosW[2]/1000);
		fprintf(MUOS3file, "%lf %lf %lf %lf\n", SimTime, SC[3].PosW[0]/1000, SC[3].PosW[1]/1000, SC[3].PosW[2]/1000);
		fprintf(MUOS4file, "%lf %lf %lf %lf\n", SimTime, SC[4].PosW[0]/1000, SC[4].PosW[1]/1000, SC[4].PosW[2]/1000);
	}
}

void ReportOrbElements(void)
{
	static FILE *outFile;
	long i;
	double rp, ra;

	outFile = FileOpen(InOutPath, "ORB_Rx.42","w");

	for(i=0;i<Norb;i++){
		rp = Orb[0].SMA*(1-Orb[i].ecc) - World[EARTH].rad;
		ra = Orb[0].SMA*(1+Orb[i].ecc) - World[EARTH].rad;
		fprintf(outFile, "%lf %lf %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", SimTime, AbsTime, i, (Orb[i].SMA-World[EARTH].rad)/1000, Orb[i].ecc, Orb[i].inc*R2D, Orb[i].MeanMotion, Orb[i].anom*R2D, Orb[i].ArgP*R2D, rp/1000, ra/1000, Orb[i].RAAN*R2D);
	}
	fclose(outFile);
}

void ReportAttitude(void)
{
	  static FILE *timefile,*AbsTimeFile;
	  static FILE *PosNfile,*VelNfile;
	  static FILE *qbnfile,*wbnfile;
	  static FILE *Hvbfile;
	  static FILE *RPYfile;
	  static FILE *Hwhlfile, *Twhlfile;
	  static FILE *Mmtbfile, *Tmtbfile;
	  static FILE *EclipseFile;
//	  static FILE *altFile;
//	  static long outFlagAttitude = 0;
//	  static long outFlagCnt = 0;

	  static char First = TRUE;
	  double CBL[3][3],Roll,Pitch,Yaw;
//	  double PosW[3],VelW[3],PosR[3],VelR[3];
	  double mtbTrq[3];

	  if (First) {
		 First = FALSE;
		 timefile = FileOpen(InOutPath,"time.42","w");
		 AbsTimeFile = FileOpen(InOutPath,"AbsTime.42","w");
		 PosNfile = FileOpen(InOutPath,"PosN.42","w");
		 VelNfile = FileOpen(InOutPath,"VelN.42","w");
		 qbnfile = FileOpen(InOutPath,"qbn.42","w");
		 wbnfile = FileOpen(InOutPath,"wbn.42","w");
		 Hvbfile = FileOpen(InOutPath,"Hvb.42","w");

		 RPYfile = FileOpen(InOutPath,"RPY.42","w");
		 Hwhlfile = FileOpen(InOutPath,"Hwhl.42","w");
		 Twhlfile = FileOpen(InOutPath,"Twhl.42","w");
		 Mmtbfile = FileOpen(InOutPath,"Mmtb.42","w");
		 Tmtbfile = FileOpen(InOutPath,"Tmtb.42","w");
		 EclipseFile = FileOpen(InOutPath, "EclipseFlag.42","w");
	  }

/*	  if(OutFlag) {
		  outFlagCnt++;
		  if(outFlagCnt == 10)
		  {
			  outFlagAttitude = 1;
			  outFlagCnt = 0;
		  }
	  }
*/
	  if (OutFlag) {
//		  outFlagAttitude = 0;
		 fprintf(timefile,"%lf\n",SimTime);
		 fprintf(AbsTimeFile,"%lf\n",AbsTime);
			fprintf(PosNfile,"%le %le %le\n",
			   SC[0].PosN[0],SC[0].PosN[1],SC[0].PosN[2]);
			fprintf(VelNfile,"%le %le %le\n",
			   SC[0].VelN[0],SC[0].VelN[1],SC[0].VelN[2]);

			fprintf(qbnfile,"%le %le %le %le\n",
			   SC[0].B[0].qn[0],SC[0].B[0].qn[1],SC[0].B[0].qn[2],SC[0].B[0].qn[3]);
			fprintf(wbnfile,"%le %le %le\n",
			   SC[0].B[0].wn[0],SC[0].B[0].wn[1],SC[0].B[0].wn[2]);
			fprintf(Hvbfile,"%18.12le %18.12le %18.12le\n",
			   SC[0].Hvb[0],SC[0].Hvb[1],SC[0].Hvb[2]);
			MxMT(SC[0].B[0].CN,SC[0].CLN,CBL);
			C2A(123,CBL,&Roll,&Pitch,&Yaw);
			fprintf(RPYfile,"%lf %lf %lf\n",Roll*R2D,Pitch*R2D,Yaw*R2D);
			fprintf(Hwhlfile,"%lf %lf %lf %lf\n",SC[0].Whl[0].H,SC[0].Whl[1].H,SC[0].Whl[2].H, SC[0].Whl[3].H);
			fprintf(Twhlfile,"%lf %lf %lf %lf\n",SC[0].Whl[0].Trq,SC[0].Whl[1].Trq,SC[0].Whl[2].Trq,SC[0].Whl[3].Trq);
			fprintf(Mmtbfile,"%lf %lf %lf\n",SC[0].MTB[0].M,SC[0].MTB[1].M,SC[0].MTB[2].M);
			mtbTrq[0] = MAGV(SC[0].MTB[0].Trq);
			mtbTrq[1] = MAGV(SC[0].MTB[1].Trq);
			mtbTrq[2] = MAGV(SC[0].MTB[2].Trq);
			fprintf(Tmtbfile,"%lf %lf %lf\n", mtbTrq[0], mtbTrq[1], mtbTrq[2]);

			fprintf(EclipseFile,"%ld\n", SC[0].Eclipse);

			CmdReport();
	  }
	  if (CleanUpFlag) {
		 fclose(timefile);
	  }

}

void projectEllipsoid(double  p[3],
                      const double w[3],
                      const double rad[3])
{
    double aa = rad[0] * rad[0];
    double bb = rad[1] * rad[1];
    double cc = rad[2] * rad[2];

    // we derive a lower limit for 'h' from  pX^2 + pY^2 + pZ^2 < max(radX,radY,radZ)^2
    double RR = fmax(aa, fmax(bb,cc));
    // 'hmin' is the minimum value that 'h' can have
    double hmin = sqrt( ( w[0]*w[0]*aa*aa + w[1]*w[1]*bb*bb + w[2]*w[2]*cc*cc ) / RR ) - RR;

    // we derive another lower limit for 'h' from  |pX| < radX
    hmin = fmax(hmin, ( fabs(w[0]) - rad[0] ) * rad[0]);

    // we derive another lower limit for 'h' from  |pY| < radY
    hmin = fmax(hmin, ( fabs(w[1]) - rad[1] ) * rad[1]);

    // we derive another lower limit for 'h' from  |pZ| < radZ
    hmin = fmax(hmin, ( fabs(w[2]) - rad[2] ) * rad[2]);

    if ( w[0]*w[0]/aa + w[1]*w[1]/bb + w[2]*w[2]/cc > 1  &&  hmin < 0 )
    {
        // if the point is outside, then 'h' should be positive:
        hmin = 0;
    }

    double h_old, h = hmin;
    //fprintf(stderr, "----- h %+f\n", h);

    /*
     Follow Newton's iteration to find the largest root.
     We start with h>0, and h should only increase
     */
    unsigned cnt = 0;
    do {
        double aah = aa + h;
        double bbh = bb + h;
        double cch = cc + h;

        double waX = w[0] / aah;
        double waY = w[1] / bbh;
        double waZ = w[2] / cch;

        double pXX = waX * waX * aa;
        double pYY = waY * waY * bb;
        double pZZ = waZ * waZ * cc;

        h_old = h;

        double   F = 1 - ( pXX         + pYY         + pZZ       );
        double  dF = 2 * ( pXX / aah   + pYY / bbh   + pZZ / cch );

        // Newton's method
        h -= F / dF;

        //fprintf(stderr, "  %i : h %+f  F %+e dh %+.20f\n", cnt, h_old, F, h-h_old);
        //fprintf(stderr, "       %+.10f   %+.10f   %+.10f   %+.10f\n", F, F/dF, ddF/dF, dddF/dF);

        if ( h < hmin )
        {
            h = 0.5 * ( h_old + hmin );
            continue;
        }

#if ( 1 )
        if ( cnt > 16 )
        {
            fprintf(stderr, "projectEllipsoid fails %u :  h %+f  F %.6e dh %.6e\n", cnt, h_old, F, h-h_old);
            //fprintf(stderr, "    pos  %+.10f     %+.10f       %+.10f\n", w[0], w[1], w[2]);
            //fprintf(stderr, "    F    %+.10f  dF %+.10f   ddF %+.10f\n", F, dF, ddF);
        }
#endif

        if ( ++cnt > 20 )
            break;

    } while ( h > h_old );

    // calculate the projection from h
    p[0] = w[0] * aa / ( aa + h );
    p[1] = w[1] * bb / ( bb + h );
    p[2] = w[2] * cc / ( cc + h );

#if ( 1 )
    // verify that projection is on ellipse
    double F = 1 - ( p[0]*p[0]/aa + p[1]*p[1]/bb + p[2]*p[2]/cc );
    fprintf(stderr, " %2i  >>> h %12.8f  F  %+e\n", cnt, h, F);
#endif
}

void TrackingSpotbeam(void)
{
	static FILE *outFile[5]; // # of MUOS
    static long First = TRUE;
	long iTx, iSb, iEl, iAz, iTh, NAz;
	double latT, lonT, altT;

	double a = 6378137.0;
    double f = 1.0/298.257222101; // GRS80
    double b = a*(1.0-f);

    double r, PosFpS[3], PosFpL[3], PosFpN[3], PosFpW[3], theta, A, B, C, az, el, CSL[3][3], CNL[3][3], z0;
    double vecConeS[3], vecConeL[3], vecConeN[3];

    if (First) {
    	char filename[20];
        First = FALSE;
		for(iTx=0; iTx<5; iTx++){
			sprintf(filename,"Spotbeam_%ld.42", iTx);
			outFile[iTx] = FileOpen(InOutPath,filename,"w");
		}
		//printf("Reference Radius of Spotbeam: %lf\n", sqrt(SpotBeam->r0sq));

		for(iTx=0; iTx<5; iTx++){
			iSb = 0;
			z0 = MAGV(SC[iTx+1].PosN);
			MT(Orb[iTx+1].CLN, CNL);

			// ECI to ECEF of Transmitter's Position
			MxV(World[EARTH].CWN, SC[iTx+1].PosN, SC[iTx+1].PosW);

	    	// ECEF to GRS80 of Transmitter's Position
	    	ECEFToGRS80(SC[iTx+1].PosW, &latT, &lonT, &altT);

	    	for(iEl=0; iEl<4; iEl++){
	    		if(iEl!=0){
	    			el = SpotBeam->el[iEl];
	    			NAz = 5;
	    		}
	    		else{
	    			el = 0;
	    			NAz = 1;
	    		}

	    		for(iAz=0; iAz<NAz; iAz++){
	    			if(iEl!=0)
	    				az = SpotBeam->az[iEl-1][iAz];
	    			else
	    				az = 0;

					// From LVLH to Spotbeam frame
					CSL[0][0] = cos(az)*cos(el); CSL[0][1] = sin(az)*cos(el); CSL[0][2] = -sin(el);
					CSL[1][0] = -sin(az);        CSL[1][1] = cos(az);         CSL[1][2] = 0;
					CSL[2][0] = cos(az)*sin(el); CSL[2][1] = sin(az)*sin(el); CSL[2][2] = cos(el);

					for(iTh=0; iTh<64; iTh++){
						theta = iTh*0.1;

						vecConeS[0] = cos(theta);
						vecConeS[1] = sin(theta);
						vecConeS[2] = 1/tan(SpotBeam->halfBW);

						MTxV(CSL, vecConeS, vecConeL);
						MTxV(Orb[iTx+1].CLN, vecConeL, vecConeN);

						A = (vecConeN[0]*vecConeN[0]+vecConeN[1]*vecConeN[1])/a/a+vecConeN[2]*vecConeN[2]/b/b;
						B = ((CNL[0][2]*vecConeN[0] + CNL[1][2]*vecConeN[1])/a/a + CNL[2][2]*vecConeN[2]/b/b)*z0;
						C = ((CNL[0][2]*CNL[0][2]+CNL[1][2]*CNL[1][2])/a/a + CNL[2][2]*CNL[2][2]/b/b)*z0*z0 - 1;

						r = (B - sqrt(B*B-A*C)) / A;

						PosFpS[0] = r * vecConeS[0];
						PosFpS[1] = r * vecConeS[1];
						PosFpS[2] = r * vecConeS[2];

						MTxV(CSL, PosFpS, PosFpL);

						PosFpL[2] -= z0;

						MxV(CNL, PosFpL, PosFpN);

						// ECI to ECEF of Transmitter's Position
						MxV(World[EARTH].CWN, PosFpN, PosFpW);

						/* ECEF to GRS80 */
						ECEFToGRS80(PosFpW, &SpotBeam->latGRS[iTx][iSb], &SpotBeam->lonGRS[iTx][iSb], &SpotBeam->altGRS[iTx][iSb]);

						SpotBeam->r[iSb] = r;

						fprintf(outFile[iTx],"%lf %ld %lf %lf %lf %lf %lf %lf %lf\n",
						   AbsTime, iSb, latT, lonT, altT,
						   SpotBeam->latGRS[iTx][iSb], SpotBeam->lonGRS[iTx][iSb], SpotBeam->altGRS[iTx][iSb], r);
					}
					iSb++;
	    		}
	    	}
	    	// Calculate entire beamwidth
			for(iTh=0; iTh<64; iTh++){
				theta = iTh*0.1;

				vecConeL[0] = cos(theta);
				vecConeL[1] = sin(theta);
				vecConeL[2] = 1/tan(SpotBeam->halfBW_whole);
				printf("%lf\n", SpotBeam->halfBW_whole*R2D);

				MTxV(Orb[iTx+1].CLN, vecConeL, vecConeN);

				A = (vecConeN[0]*vecConeN[0]+vecConeN[1]*vecConeN[1])/a/a+vecConeN[2]*vecConeN[2]/b/b;
				B = ((CNL[0][2]*vecConeN[0] + CNL[1][2]*vecConeN[1])/a/a + CNL[2][2]*vecConeN[2]/b/b)*z0;
				C = ((CNL[0][2]*CNL[0][2]+CNL[1][2]*CNL[1][2])/a/a + CNL[2][2]*CNL[2][2]/b/b)*z0*z0 - 1;

				r = (B - sqrt(B*B-A*C)) / A;

				PosFpL[0] = r * vecConeL[0];
				PosFpL[1] = r * vecConeL[1];
				PosFpL[2] = r * vecConeL[2] - z0;

				MxV(CNL, PosFpL, PosFpN);

				// ECI to ECEF of Transmitter's Position
				MxV(World[EARTH].CWN, PosFpN, PosFpW);

				/* ECEF to GRS80 */
				ECEFToGRS80(PosFpW, &SpotBeam->latGRS[iTx][iSb], &SpotBeam->lonGRS[iTx][iSb], &SpotBeam->altGRS[iTx][iSb]);

				SpotBeam->r[iSb] = r;

				fprintf(outFile[iTx],"%lf %ld %lf %lf %lf %lf %lf %lf %lf\n",
				   AbsTime, iSb, latT, lonT, altT,
				   SpotBeam->latGRS[iTx][iSb], SpotBeam->lonGRS[iTx][iSb], SpotBeam->altGRS[iTx][iSb], r);
			}
		}
    }
}

void CmdReport(void)
{
    static FILE *qbrFile, *therrFile, *werrFile, *WhlCmdFile, *MtqCmdFile;
    static long First = 1;

    if (First) {
       First = 0;
       qbrFile = FileOpen(InOutPath,"qbr.42","wt");
       therrFile = FileOpen(InOutPath,"therr.42","wt");
       werrFile = FileOpen(InOutPath,"werr.42","wt");
       WhlCmdFile = FileOpen(InOutPath, "WhlCmd.42", "wt");
       MtqCmdFile = FileOpen(InOutPath, "MtqCmd.42", "wt");
    }
    fprintf(qbrFile,"%le %le %le %le\n",
       SC[0].AC.qbr[0], SC[0].AC.qbr[1], SC[0].AC.qbr[2], SC[0].AC.qbr[3]);
    fprintf(werrFile,"%le %le %le \n",
       SC[0].AC.CfsCtrl.werr[0], SC[0].AC.CfsCtrl.werr[1], SC[0].AC.CfsCtrl.werr[2]);
    fprintf(therrFile,"%le %le %le \n",
       SC[0].AC.CfsCtrl.therr[0], SC[0].AC.CfsCtrl.therr[1], SC[0].AC.CfsCtrl.therr[2]);
    fprintf(WhlCmdFile,"%le %le %le \n",
       SC[0].AC.Tcmd[0], SC[0].AC.Tcmd[1], SC[0].AC.Tcmd[2]);
    fprintf(MtqCmdFile,"%le %le %le \n",
       SC[0].AC.Mcmd[0], SC[0].AC.Mcmd[1], SC[0].AC.Mcmd[2]);
}

void RepeatGroundTrack(struct OrbitType *O)
{
	double nPeriod, nDays = 0;
	double tolerance, dlon_track;
	double dlon, temp;
	double lon[2000], currLon;
	double wE = 7.2921151467E-5;
	uint32_t nOrbits = 0, i, j;
	uint8_t flagTrack = 1, flagRepeat = 0;

	tolerance = 0.1*D2R;
	dlon_track = 0.9*D2R;
	/* ECI to ECEF of Receiver's Position */
//	MxV(World[EARTH].CWN,S->PosN,PosWR);
	/* ECEF to WGS84 */
//	ECEFToWGS84(PosWR, &lat, &lon, &alt);

	lon[nOrbits] = 0;

	nPeriod = TwoPi / (O->MeanMotion + O->ArgPdot);
	dlon = nPeriod * (wE - O->RAANdot);

	printf("dlon: %lf rad, %lf deg\n", dlon, dlon*R2D);

	while(1)
	{
		currLon = lon[nOrbits++] + dlon;
		if(currLon >= TwoPi)
			currLon = currLon - TwoPi;
		if(currLon < 0)
			currLon = -currLon;

		lon[nOrbits] = currLon;

		for(i=0; i<nOrbits; i++){
			for(j=i+1; j<nOrbits; j++){
				if(lon[i]>lon[j]){
					temp = lon[i];
					lon[i] = lon[j];
					lon[j] = temp;
				}
			}
		}
		for(i=0; i<nOrbits; i++){
			if(lon[i+1]-lon[i]>dlon_track){
				flagTrack = 0;
			}
		}
		if(flagTrack)
			break;
		else
			flagTrack = 1;

		if(lon[nOrbits] <= tolerance){
			flagRepeat = 1;
			break;
		}
	}
	printf("0: %lf, 1: %lf\n", lon[0], lon[1]);
	nDays = nOrbits*nPeriod;
	if(flagTrack)
		printf("Distance btw ground tracks within 100km\n");
	if(flagRepeat)
		printf("SC passes over the same location.\n");
	printf("nDays: %lf days = %lf sec \n", nDays/86400, nDays);
	printf("nOrbits: %d\n", nOrbits);
	printf("nPeriod: %lf min\n", nPeriod/60);
	printf("tol: %lf\n", tolerance);
	printf("lon: %lf\n", lon[nOrbits]);
	printf("eccentricity: %lf\n", O->ecc);
	printf("Periapsis altitude: %lf\n", O->rmin-World[EARTH].rad);
	printf("Apoapsis altitude: %lf\n", (1+O->ecc)*O->SMA-World[EARTH].rad);
	printf("Inclination: %lf\n", O->inc*R2D);
	printf("RAAN: %lf\n", O->RAAN*R2D);
	printf("Argument of Periapsis: %lf\n", O->ArgP*R2D);
	printf("True Anomaly: %lf\n", O->anom*R2D);
//		AdvanceTime();
//		EnckeRK4(S);

//		O->RAAN = O->RAAN0 + O->RAANdot*(SimTime - 0.5/O->MeanMotion*sin(2.0*O->ArgP+2.0*O->anom));
//		O->ArgP = O->ArgP0 + O->ArgPdot*SimTime;

//		Eph2RV(O->MuPlusJ2,O->SLR,O->ecc,
//				O->inc,O->RAAN,O->ArgP,
//				AbsTime+DTSIM-O->tp,
//				O->PosN,O->VelN,&O->anom);

//		Ephemerides_Earth();
//		GravPertForce(S);

		/* ECI to ECEF of Receiver's Position */
//		MxV(World[EARTH].CWN,S->PosN,PosWR);
		/* ECEF to WGS84 */
//		ECEFToWGS84(PosWR, &lat, &lon, &alt);
}

void InitFOVs(void)
{
      FILE *infile;
      char junk[120],newline;
      char response[120],response1[120],response2[120];
      double Ang1,Ang2,Ang3;
      long Seq;
      long i;

      infile = FileOpen(InOutPath,"Inp_FOV.txt","r");
      fscanf(infile,"%[^\n] %[\n]",junk,&newline);
      fscanf(infile,"%[^\n] %[\n]",junk,&newline);
      fscanf(infile,"%ld %[^\n] %[\n]",&Nfov,junk,&newline);
      FOV = (struct FovType *) calloc(Nfov,sizeof(struct FovType));
      for(i=0;i<Nfov;i++) {
         fscanf(infile,"%[^\n] %[\n]",junk,&newline);
         fscanf(infile,"\"%[^\"]\" %[^\n] %[\n]",
            FOV[i].Label,junk,&newline);
         fscanf(infile,"%ld %lf %[^\n] %[\n]",
            &FOV[i].Nv,&FOV[i].Length,junk,&newline);
         fscanf(infile,"%lf %lf %[^\n] %[\n]",
            &FOV[i].Width,&FOV[i].Height,junk,&newline);
         if (FOV[i].Width >= 180.0) {
            printf("FOV[%ld] Width >= 180 deg.  This is not allowed.  Bailing out.\n",i);
            exit(1);
         }
         if (FOV[i].Height >= 180.0) {
            printf("FOV[%ld] Width >= 180 deg.  This is not allowed.  Bailing out.\n",i);
            exit(1);
         }
         FOV[i].Width *= D2R;
         FOV[i].Height *= D2R;
         fscanf(infile,"%f %f %f %f %[^\n] %[\n]",
            &FOV[i].Color[0],&FOV[i].Color[1],&FOV[i].Color[2],&FOV[i].Color[3],
            junk,&newline);
         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
         FOV[i].Type = DecodeString(response);
         fscanf(infile,"%s %s %[^\n] %[\n]",
            response1,response2,junk,&newline);
         FOV[i].NearExists = DecodeString(response1);
         FOV[i].FarExists = DecodeString(response2);
         fscanf(infile,"%ld %ld %[^\n] %[\n]",
            &FOV[i].SC,&FOV[i].Body,junk,&newline);
         if (FOV[i].SC >= Nsc) {
            printf("FOV[%ld].SC is out of range.\n",i);
            exit(1);
         }
         if (SC[FOV[i].SC].Exists && FOV[i].Body >= SC[FOV[i].SC].Nb) {
            printf("FOV[%ld].Body is out of range.\n",i);
            exit(1);
         }
         FOV[i].RefOrb = SC[FOV[i].SC].RefOrb;
         if (!SC[FOV[i].SC].Exists) {
            FOV[i].NearExists = FALSE;
            FOV[i].FarExists = FALSE;
         }

         fscanf(infile,"%lf %lf %lf %[^\n] %[\n]",
            &FOV[i].pb[0],&FOV[i].pb[1],&FOV[i].pb[2],junk,&newline);
         fscanf(infile,"%lf %lf %lf %ld %[^\n] %[\n]",
            &Ang1,&Ang2,&Ang3,&Seq,junk,&newline);
            A2C(Seq,Ang1*D2R,Ang2*D2R,Ang3*D2R,FOV[i].CB);
      }
      fclose(infile);
}

void InitSpotbeams(void)
{
	double r_beams, r0, rS1, rS2, rS3; // Beam spacing
	int i;
	double theta, temp;
	double a = 6378137.0;
    double f = 1.0/298.257222101; // GRS80
    double b = a*(1.0-f);
	r_beams = 1+4*cos(36*D2R);
//	r_beams = 3+2*sin(18*D2R);
//	r_beams = sqrt(1+8*cos(36*D2R)*cos(36*D2R)+8*sin(18*D2R)*sin(18*D2R)+8*sin(18*D2R));
	r0 = b/r_beams;
	rS1 = 2*r0*cos(36*D2R);
	rS2 = (2+2*sin(18*D2R))*r0;
	rS3 = 2*rS1;

	SpotBeam->el[0] = 0;
	SpotBeam->el[1] = atan2(rS1, SpotBeam->alt);
	SpotBeam->el[2] = atan2(rS2, SpotBeam->alt);
	SpotBeam->el[3] = atan2(rS3, SpotBeam->alt);

	for(i=0; i<5; i++){
		SpotBeam->az[0][i] = (-18+72*i)*D2R;
		SpotBeam->az[1][i] = (18+72*i)*D2R;
		SpotBeam->az[2][i] = (-18+72*i)*D2R;
	}

	SpotBeam->r0sq = r0*r0*SpotBeam->scaleFactor*SpotBeam->scaleFactor;
	SpotBeam->x[0] = 0; SpotBeam->y[0] = 0;
	SpotBeam->x[1] = rS1; SpotBeam->y[1] = 0;
	SpotBeam->x[2] = rS1*cos(72*D2R); SpotBeam->y[2] = rS1*sin(72*D2R);
	SpotBeam->x[3] = rS1*cos(144*D2R); SpotBeam->y[3] = rS1*sin(144*D2R);
	SpotBeam->x[4] = rS1*cos(216*D2R); SpotBeam->y[4] = rS1*sin(216*D2R);
	SpotBeam->x[5] = rS1*cos(288*D2R); SpotBeam->y[5] = rS1*sin(288*D2R);
	SpotBeam->x[6] = rS3; SpotBeam->y[6] = 0;
	SpotBeam->x[7] = rS2*cos(36*D2R); SpotBeam->y[7] = rS2*sin(36*D2R);
	SpotBeam->x[8] = rS3*cos(72*D2R); SpotBeam->y[8] = rS3*sin(72*D2R);
	SpotBeam->x[9] = rS2*cos(108*D2R); SpotBeam->y[9] = rS2*sin(108*D2R);
	SpotBeam->x[10] = rS3*cos(144*D2R); SpotBeam->y[10] = rS3*sin(144*D2R);
	SpotBeam->x[11] = rS2*cos(180*D2R); SpotBeam->y[11] = rS2*sin(180*D2R);
	SpotBeam->x[12] = rS3*cos(216*D2R); SpotBeam->y[12] = rS3*sin(216*D2R);
	SpotBeam->x[13] = rS2*cos(252*D2R); SpotBeam->y[13] = rS2*sin(252*D2R);
	SpotBeam->x[14] = rS3*cos(288*D2R); SpotBeam->y[14] = rS3*sin(288*D2R);
	SpotBeam->x[15] = rS2*cos(324*D2R); SpotBeam->y[15] = rS2*sin(324*D2R);

	theta = 18*D2R;

	for(i=0; i<16; i++){
		temp = cos(theta)*SpotBeam->x[i] + sin(theta)*SpotBeam->y[i];
		SpotBeam->y[i] = cos(theta)*SpotBeam->y[i] - sin(theta)*SpotBeam->x[i];
		SpotBeam->x[i] = temp;
	}

	SpotBeam->r0Target = SpotBeam->alt * tan(SpotBeam->halfBW);
	SpotBeam->halfBW_whole = asin(b/(SpotBeam->alt + b));
}

void antPattern_Dipole(struct SoOpType *SoOp)
{
    int i, j;
    double phi;

    SoOp->AntPattern = (double **)calloc(180, sizeof(double));
    for(i=0; i<180; i++)
    	SoOp->AntPattern[i] = (double *)calloc(360, sizeof(double));

    for(i=0; i<180; i++){
    	for(j=0; j<360; j++){
    		phi = D2R*j;
    		SoOp->AntPattern[i][j] = SoOp->AntMaxGainE + 10*log10(cos(HalfPi*cos(phi))*cos(HalfPi*cos(phi))/sin(phi)/sin(phi));
    	}
    }
}

void antPattern_Patch(struct SoOpType *SoOp)
{
    int i, j;
    double phi;

    SoOp->AntPattern = (double **)calloc(180, sizeof(double));
    for(i=0; i<180; i++)
    	SoOp->AntPattern[i] = (double *)calloc(360, sizeof(double));

    for(i=0; i<180; i++){
    	for(j=0; j<360; j++){
    		phi = D2R*j;
    		SoOp->AntPattern[i][j] = SoOp->AntMaxGainE + 10*log10(3282.81*((sin(phi*1.7724)*sin(phi*1.7724))/(phi*1.7724)/(phi*1.7724)));
    	}
    }
}

void antPattern_User(struct SoOpType *SoOp)
{
	FILE *infile;
    const char s[2] = ",";
    char *token;
    int i, j, k;
    char line[3500];

    SoOp->AntPattern = (double **)calloc(SoOp->Ntheta, sizeof(double));
    for(i=0; i<SoOp->Ntheta; i++)
    	SoOp->AntPattern[i] = (double *)calloc(SoOp->Nphi, sizeof(double));

    k = 0;
    j = 0;

    infile = FileOpen(InOutPath,SoOp->AntFileName,"r");
	while(fgets(line, sizeof line, infile) != NULL)
	{
		token = strtok(line, s);
		while(token!=NULL)
		{
			SoOp->AntPattern[j][k] = atof(token);
		    //printf("(th,ph)=(%d,%d) g: %lf \n", j, k, SoOp->AntPattern[j][k]);
			k++;
			token = strtok(NULL, s);
		}
		j++;
		k=0;
	}
	fclose(infile);
}

void calcAntPattern(struct SoOpType *SoOp, struct SpecularType *Sp)
{
	long thetaInt, phiInt;

	// Earth-view Antenna Gain
	thetaInt = round(Sp->lookAngElR*R2D);
	phiInt = round(Sp->lookAngAzR*R2D);

	if(phiInt<0)
		phiInt = phiInt + 360;
	if(phiInt>359)
		phiInt = 0;
	if(thetaInt>179)
		thetaInt = 0;
	//printf("R:(%ld,%ld)      ", thetaInt, phiInt);
	Sp->AntGainE = SoOp->AntPattern[thetaInt][phiInt];

	// Sky-view Antenna Gain
	thetaInt = round(Sp->lookAngElD*R2D);
	phiInt = round(Sp->lookAngAzD*R2D);

	if(phiInt<0)
		phiInt = phiInt + 360;
	if(phiInt>359)
		phiInt = 0;
	if(thetaInt>179)
		thetaInt = 0;
	//printf("D:(%ld,%ld)\n", thetaInt, phiInt);
	Sp->AntGainS = SoOp->AntPattern[thetaInt][phiInt];
}

/**********************************************************************/
/* Calculate complex dielectric constant using Dolsan model           */
/* from Equation 14 in Hallikinen 1985                                */
/* e_c = complex dielectric constant                                  */
/* incAng = incidence angle (rad)                                     */
/* F = Fresnel coefficients                                           */
/* Gamma = Reflectivity                                               */
/* This calculation is based on Equation 4.132 in Ulaby Book          */
/* Assumption: Relative permeability is approximated as 1 in this calculation */
void calcReflectivity(double complex e_c, struct ReflectType *R, double incAng)
{
	/* Calculate standard Fresnel coefficients */
	// From Equation 4.132a in Ulaby Book
	R->FHH = (cos(incAng)-sqrt(e_c-sin(incAng)*sin(incAng)))/
		  (cos(incAng)+sqrt(e_c-sin(incAng)*sin(incAng)));
	// From Equation 4.132b in Ulaby Book
	R->FVV = (e_c*cos(incAng)-sqrt(e_c-sin(incAng)*sin(incAng)))/
		  (e_c*cos(incAng)+sqrt(e_c-sin(incAng)*sin(incAng)));
	R->FLR = 0.5*(R->FVV-R->FHH);
	R->FRR = 0.5*(R->FVV+R->FHH);
		//printf("FVV: %f + i%f\n", creal(R->FVV), cimag(R->FVV));
		//printf("FHH: %f + i%f\n", creal(R->FHH), cimag(R->FHH));

	// From Equation 4.132a in Ulaby Book
	//Gamma_h = (abs((cos(theta_rad)-sqrt(e_c-(sin(theta_rad)).^2))./...
	//    (cos(theta_rad)+sqrt(e_c-(sin(theta_rad)).^2)))).^2;
	// From Equation 4.132b in Ulaby Book
	//Gamma_v = (abs((e_c.*cos(theta_rad)-sqrt(e_c-(sin(theta_rad)).^2))./...
	//    (e_c.*cos(theta_rad)+sqrt(e_c-(sin(theta_rad)).^2)))).^2;
	//From Zhang Fusion 2010
	//Gamma = (Gamma_v + Gamma_h)./2;

	R->gammaV = R->FVV * conj(R->FVV);
	R->gammaH = R->FHH * conj(R->FHH);
	R->gammaLR = R->FLR * conj(R->FLR);
	R->gammaRR = R->FRR * conj(R->FRR);

	R->RcHH = cabs(R->FHH);
	R->RcVV = cabs(R->FVV);
	R->RcLR = cabs(R->FLR);
	R->RcRR = cabs(R->FRR);
}

/**********************************************************************************/
/*   INPUTS:
 *   f_Hz:       Frequency (Hertz)
 *   VSM:        Volumetric soil moisture (cm3/cm3) [0,1]
 *   clay_ratio: Mass fraction of clay content in soil
 *
 *   Implemented from the following paper:
 *   V. L. Mironov and S. V Fomin, �Temperature dependable microwave
 *   dielectric model for moist soils,� PIERS Proceedings, March, pp. 23-27, 2009.*/
void calcDielMironov(uint16_t rowSt, uint16_t colSt, uint16_t rowDy, uint16_t colDy, double f_Hz, struct SoilType *Soil)
{
	float C = SurfaceStt->clay[rowSt][colSt];
	float Temp = SurfaceDyn->soilTemp1[rowDy][colDy] - 273.15; // Kelvin to Celsius
	float VSM = SurfaceDyn->RZSM[rowDy][colDy];

	// Mironov's regression expressions based on Curtis, Dobson, and Hallikainen datasets
	double nd, kd, mvt, eps0b, taub, sigb, sigu, eps0u, tauu;
	nd = 1.634 - 0.539e-2 * C + 0.2748e-4 * C * C;     // Eqn 17
	kd = 0.03952 - 0.04038e-2 * C ;                    // Eqn 18
	mvt = 0.02863 + 0.30673e-2 * C ;                   // Eqn 19
	eps0b = 79.8 - 85.4e-2 * C + 32.7e-4 * C * C;      // Eqn 20
	taub = 1.062e-11 + 3.450e-12 * 1e-2 * C ;          // Eqn 21
	sigb = 0.3112 + 0.467e-2 * C ;                     // Eqn 22
	sigu = 0.3631 + 1.217e-2 * C ;                     // Eqn 23
	eps0u = 100 ;                                      // Eqn 24
	tauu = 8.5e-12 ;                                   // Eqn 25

	// Debye relaxation equations for water as a function of frequency
	double eps0 = 8.854e-12; // Vacuum permittivity
	double epsinf = 4.9;
	double epsb_real, epsb_imag, epsu_real, epsu_imag;
	// Section IV - Epn 16
	epsb_real = epsinf + ( (eps0b - epsinf) / (1 + (TwoPi * f_Hz * taub)*(TwoPi * f_Hz * taub)) );

	epsb_imag = (eps0b - epsinf) / (1 + (TwoPi * f_Hz * taub)*(TwoPi * f_Hz * taub))
	    		* (TwoPi * f_Hz * taub) + sigb / (TwoPi * eps0 * f_Hz);

	epsu_real = epsinf + ( (eps0u - epsinf) / (1 + (TwoPi * f_Hz * tauu)*(TwoPi * f_Hz * tauu)) );

	epsu_imag = (eps0u - epsinf) / (1 + (TwoPi * f_Hz * tauu)*(TwoPi * f_Hz * tauu) )
	    		* (TwoPi * f_Hz * tauu) + sigu / (TwoPi * eps0 * f_Hz);

	// Refractive indices - Eqn 14
	double nb, kb, nu, ku;
	nb = 1/sqrt(2) * sqrt( sqrt(epsb_real*epsb_real + epsb_imag*epsb_imag) + epsb_real );
	kb = 1/sqrt(2) * sqrt( sqrt(epsb_real*epsb_real + epsb_imag*epsb_imag) - epsb_real );
	nu = 1/sqrt(2) * sqrt( sqrt(epsu_real*epsu_real + epsu_imag*epsu_imag) + epsu_real );
	ku = 1/sqrt(2) * sqrt( sqrt(epsu_real*epsu_real + epsu_imag*epsu_imag) - epsu_real );

	// n(*) are refractive indices, k(*) are normalized attenuation coefficients
	// m: moist soil
	// d: dry soil
	// b: bound soil water (BSW)
	// u: unbound (free) soil water (FSW)

	double nm, km, er_r_real, er_r_imag, tmp;
	if(VSM>mvt){
		nm = nd + (nb - 1) * mvt + (nu - 1) * (VSM - mvt);   // Eqn 12
		km = kd + kb * mvt + ku * (VSM - mvt);               // Eqn 13
	}
	else{
		nm = nd + (nb - 1) * VSM;         //  Eqn 12
		km = kd + kb * VSM;               // Eqn 13
	}
	Soil->e_p = nm*nm - km*km;          // Eqn 11
	Soil->e_dp = 2 * nm * km;           // Eqn 11

	// Combine the dielectric constant (complex number)
	Soil->e_c = Soil->e_p - I*Soil->e_dp;
	Soil->m_v = VSM;
	Soil->clay = C;
}

/**********************************************************************/
/* Calculate complex dielectric constant using a semi-empirical model */
/* described in Dobson 1985 for frequency range 1.4-18 GHz with a     */
/* correction for frequency range 0.3 - 1.3 GHz given in Peplinski 1995 */
// Input: EASE-grid2.0 row and column for static and dynamic data, Signal Frequency in Hertz
void calcDielPeplinski(uint16_t rowSt, uint16_t colSt, uint16_t rowDy, uint16_t colDy, double Freq, struct SoilType *Soil)
{
	double alpha = 0.65; // empirically determined constant - from Peplinski 1995 p.804
//	double rho_b = 1.2; // bulk density of soil sample in [g/cc] (grams/cubic centimeter)
	double rho_s = 2.66; // specific density of solid soil particles in [g/cm^3]
//	double rho_w = 0.997; // specific density of water in [g/cm^3]
//	double rho_i = 0.917; // specific density of ice in [g/cm^3]
	double e_s; // relative permittivity of solid soil, Equation 22 in Dobson 1985
	double e_i = 3.15; // relative permittivity of inclusions (air, bound water, and free water) from L.Zhang 2003

	double e_p_fw; // real part of the relative dielectric constant of free water
	double e_dp_fw; // imaginary part of the relative dielectric constant of free water
	double e_winf = 4.9; // high-frequency (or optical) limit of e_p_fw from Equation E.16 in Ulaby Vol. III, Appendix E.2
	double e_0 = 8.85418782E-12; // permittivity of free space in [m^-3 kg^-1 s^4 A^2]
	double e_w0; //static dielectric constant of pure water
	double tau_w; // relaxation time of pure water [sec]

	double sigma_eff; // Effective conductivity (dependent on frequency)

	double mvu; // Unfrozen volumetric moisture content (100v/v)
	double mvi; // Ice volumetric content (100v/v)

	double beta_p, beta_dp; // Soil dependent constants
	double TK; // Soil temperature in [K]

	float sand = SurfaceStt->sand[rowSt][colSt];
	float clay = SurfaceStt->clay[rowSt][colSt];
	float rho_b = SurfaceStt->rho_b[rowSt][colSt];
	float Temp = SurfaceDyn->soilTemp1[rowDy][colDy] - 273.15; // Kelvin to Celsius
	float RZSM = SurfaceDyn->RZSM[rowDy][colDy];

	e_s = (1.01 + 0.44*rho_s)*(1.01 + 0.44*rho_s) - 0.062;

	// Equation 4 in Peplinski 1995
	beta_p = 1.2748 - 0.519*sand - 0.152*clay;
	// Equation 5 in Peplinski 1995
	beta_dp = 1.33797 - 0.603*sand - 0.166*clay;

	/* Relative dielectric constant of free water, given by Debye-type dispersion equation (Ulaby Vol.III, Appendix E.2) */
	// Equation E.19 in Ulaby Vol.III, Appendix E.2
	e_w0 = 88.045 - 0.4147*Temp + 6.295E-4*Temp*Temp + 1.075E-5*Temp*Temp*Temp;
	// Equation E.17 in Ulaby Vol.III, Appnedix E.2
	tau_w = (1.1109E-10 - 3.824E-12*Temp + 6.938E-14*Temp*Temp - 5.096E-16*Temp*Temp*Temp)/(TwoPi);

	// Real part
	// Equation 6 in Peplinski 1995
	e_p_fw = e_winf + (e_w0 - e_winf) / ( 1 + (TwoPi*Freq*tau_w)*(TwoPi*Freq*tau_w));

	if(Freq>1.4E9)
		// Equation 8 in Peplinski 1995
		sigma_eff = -1.645 + 1.939*rho_b - 2.25622*sand + 1.549*clay;
	else
		// Equation 10 in Peplinski 1995
		sigma_eff = 0.0467 + 0.2204*rho_b - 0.4111*sand + 0.6614*clay;


	// Frozen soil - from Zhang 2003. Note that paper used T in K, here we change to deg C  // need to revise...
	TK = Temp + 273.15;
	//mvu = A * pow(abs(TK-273.2), -B) * rho_b / rho_w;// eq (2) from Zhang 2003
	//mvi = (Soil->m_v - mvu) * (rho_w/rho_i); // eq (3) from Zhang - believed to be corrected.
	mvu = RZSM;
	mvi = 0;

	// Imaginary part
	// Equation 7 in Peplinski 1995
	e_dp_fw = (TwoPi*Freq*tau_w*(e_w0 - e_winf))/(1+(TwoPi*Freq*tau_w)*(TwoPi*Freq*tau_w))
			+ sigma_eff/(TwoPi*e_0*Freq) * (rho_s - rho_b)/(rho_s*RZSM);

	/* Complex dielectric constant */
	// Real part (Equation 2 in Peplinski 1995)
	//epsilon_m_p = (1 + rho_b./rho_s.*(epsilon_s.^alpha-1) + m_v.^beta_p .* epsilon_p_fw.^alpha - m_v).^(1/alpha);
	// updated (6) from Zhang 2003
	Soil->e_p = pow(1 + rho_b /rho_s*(pow(e_s,alpha)-1) + pow(mvu,beta_p) * pow(e_p_fw,alpha)
			- RZSM + mvi*pow(e_i,alpha) - mvi,1/alpha);

	// Imaginary part (Equation 3 in Peplinski 1995)
	// epsilon_m_dp = (m_v.^beta_dp .* epsilon_dp_fw.^alpha).^(1/alpha);
	// updated from Zhang 2003 (my interpretation ....) also note that this gives a small imag component (~E-16).
	Soil->e_dp = newpow(pow(mvu,beta_dp)*newpow(e_dp_fw,alpha), 1/alpha);

	// Correction to the real part if f < 1.4 GHz (Equation 9 in Peplinski 1995)
	if (Freq < 1.4E9)
		Soil->e_p = 1.15 * Soil->e_p - 0.68;
	Soil->e_c = Soil->e_p - I*Soil->e_dp;

	Soil->m_v = RZSM;
	Soil->T = Temp;
	Soil->sand = sand;
	Soil->clay = clay;
	Soil->rho_b = rho_b;
}

/************************************************************************************/
/* Calculate penetration and skin depth based on frequency and relative complex     */
/* dielectric constant.                                                             */
/* This calculation is based on Eq.11.43-11.45 (p.847) in Ulaby Book, Vol.II, 1986. */
void calcDepth(struct SoilType Sm, struct TxType Transmitter)
{
	double attenuation; // field attenuation coefficient

	// Equation 11.44
	attenuation = TwoPi/Transmitter.Wavelength * cabs(cimag(csqrt(Sm.e_c)));
	// attenuation = pi./wavelength.*(imag(e_c))./(sqrt((real(e_c))));

	Sm.depthS = 1/attenuation;

	// Equation 11.43
	Sm.depthP = cabs(1/(2*attenuation));
	// Sm.depthP = Transmitter.Wavelength * sqrt(Sm.e_p) / (TwoPi*Sm.e_dp);

	//double e, del;
	//e = sqrt(Sm.e_p*Sm.e_p + Sm.e_dp*Sm.e_dp);
	//del = atan2(Sm.e_dp,Sm.e_p);
	//Sm.depthP = cabs(Transmitter.Wavelength/(2*TwoPi*sqrt(e)*sin(del/2)));
}

/**************************************************************************************/
/* Estimate the reflectivity and penetration depth for a single semi-infinite medium. */
/* Soil dielectric constant from (Peplinski 1995) will be used.                       */
void coherent_error_model(void)
{

}
/*
void surfaceGeometry(double r, double theta, double phi)
{
	double r_s[3]; // A local surface frame

	r_s[0] = r * sin(theta) * cos(phi);
	r_s[1] = r * sin(theta) * sin(phi);
	r_s[2] = r * cos(theta);

	double C_spec_surf[3][3]; // local rotation matrix between the specular frame and the surface frame

	C_spec_surf[0][0] = sin(theta) * cos(phi);
	C_spec_surf[0][1] = sin(theta) * sin(phi);
	C_spec_surf[0][2] = cos(theta);

	C_spec_surf[1][0] = cos(theta) * cos(phi);
	C_spec_surf[1][1] = cos(theta) * sin(phi);
	C_spec_surf[1][2] = -sin(theta);

	C_spec_surf[2][0] = -sin(phi);
	C_spec_surf[2][1] = cos(phi);
	C_spec_surf[2][2] = 0;

	double dA;

	dA = r*r*sin(theta); // differential area over the surface

	double r_surf_rx[3], r_surf_tx[3], temp[3];

}
*/

void LoadSoilData(void)
{
	FILE * pFile;
	long lSize;
	size_t result;

	/*.. Root Zone Soil Moisture ..*/
	pFile = FileOpen(AncillaryPath, SurfaceDyn->filename_RZSM, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread (SurfaceDyn->RZSM, sizeof(float), lSize/sizeof(float), pFile);
	if (result != lSize/sizeof(float)) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	printf("RZSM loaded!\n");

	/*.. Soil Temperature ..*/
	pFile = FileOpen(AncillaryPath, SurfaceDyn->filename_soilTemp1, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread (SurfaceDyn->soilTemp1, sizeof(float), lSize/sizeof(float), pFile);
	if (result != lSize/sizeof(float)) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	printf("Soil Temperature Layer 1 loaded!\n");
}

void LoadSnowCover(void)
{
	FILE * pFile;
	long lSize;
	size_t result;
	double lat, lon;
	uint16_t row, col;

	/*.. Snow Cover ..*/
	pFile = FileOpen(AncillaryPath, SurfaceDyn->filename_snowCover[Month-1], "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);

	// copy the file into the buffer:
	result = fread (SurfaceDyn->snowCover, 1, lSize, pFile);
	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

	// terminate
	fclose (pFile);

	// Sanity Check
	lat = 70.525;
	lon = 56.475;
	row = 1799.9999 - lat*20;//(90-lat)*20;
	col = 3599.9999 + lon*20;//(lon+180)*20;

//	if(Surface->snow[row][col]!=6){
//		printf("Error loading snow mask\n");
//		exit(1);
//	}
//	else
		//printf("(%d,%d) = %d\n", row, col, LandMask[row][col]);
		printf("Snow mask %ld loaded!\n", Month);
}

void UpdateOrbitCycle(void)
{
	SC[0].Ncycle = floor(SimTime/Orb[0].Period)+1;
	SC[1].Ncycle = floor(SimTime/Orb[1].Period)+1;
	SC[2].Ncycle = floor(SimTime/Orb[2].Period)+1;
	SC[3].Ncycle = floor(SimTime/Orb[3].Period)+1;
	SC[4].Ncycle = floor(SimTime/Orb[4].Period)+1;
	SC[5].Ncycle = floor(SimTime/Orb[5].Period)+1;
	SC[6].Ncycle = floor(SimTime/Orb[6].Period)+1;
	SC[7].Ncycle = floor(SimTime/Orb[7].Period)+1;
	SC[8].Ncycle = floor(SimTime/Orb[8].Period)+1;
	SC[9].Ncycle = floor(SimTime/Orb[9].Period)+1;
	SC[10].Ncycle = floor(SimTime/Orb[10].Period)+1;
	SC[11].Ncycle = floor(SimTime/Orb[11].Period)+1;
	SC[12].Ncycle = floor(SimTime/Orb[12].Period)+1;
	SC[13].Ncycle = floor(SimTime/Orb[13].Period)+1;
	SC[14].Ncycle = floor(SimTime/Orb[14].Period)+1;
	SC[15].Ncycle = floor(SimTime/Orb[15].Period)+1;
}

void UpdateSnowCover(void)
{
	if(Day==1 && Hour==0 && Minute==0 && Second==0){
		LoadSnowCover();
	}
}

/* Convert WGS84 geodetic coordinates to row, column index for 1km resolution EASE grid 2.0
   Input: lat, lon in radian */
void LL2EzRC1km(double lat, double lon, uint16_t *row, uint16_t *col)
{
	double e = 0.081819190842621; // eccentricity sqrt(2f-f^2)
	double c = 1000.90; // interdistance pixel
	double nc = 34704; // Number of columns
	double nr = 14616; // Number of lines
	double s0 = 7307.5; //(nr-1)/2;
	double r0 = 17351.5; //(nc-1)/2;
	double q, x, y;
	double aXk0 = 5528256.639292833; // a*k0
	double e2 = 0.006694379990141; // e*e
	double oneMinuse2 = 0.993305620009859; // 1-e2
	double oneOver2e = 6.111035746634866; // 1/(2*e)
	double aOver2k0 = 3679336.384424151; // a/(2*k0)
	double sinLat = sin(lat);
//	double a = 6378137;
//	double phi1 = 30*D2R;
//	double k0 = cos(phi1) / sqrt(1 - (e*e*sin(phi1)*sin(phi1)));

//	q = (1-e*e) * ( (sin(lat)/(1-e*e*sin(lat)*sin(lat))) - (1/(2*e))*log((1-e*sin(lat))/(1+e*sin(lat))) );
//	x= a*k0 * lon;
//	y= q * a/(2*k0);

	q = oneMinuse2 * ( (sinLat/(1-e2*sinLat*sinLat)) - oneOver2e*log((1-e*sinLat)/(1+e*sinLat)) );
	x= aXk0 * lon;
	y= q * aOver2k0;

	*col = round(r0 + (x/c)); // as Brodzik et al formula => r [0 1387]
	*row = round(s0 - (y/c)); // as Brodzik et al. formula => s [0 583] +1 included in nl
}

/* Convert WGS84 geodetic coordinates to row, column index for 9km resolution EASE grid 2.0
   Input: lat, lon in radian */
void LL2EzRC9km(double lat, double lon, uint16_t *row, uint16_t *col)
{
	double e = 0.081819190842621; // eccentricity sqrt(2f-f^2)
	double c = 9008.05; // interdistance pixel
	double nc = 3856; // Number of columns
	double nr = 1624; // Number of lines
	double s0 = 811.5; //(nr-1)/2;
	double r0 = 1927.5; //(nc-1)/2;
	double q, x, y;
	double aXk0 = 5528256.639292833; // a*k0
	double e2 = 0.006694379990141; // e*e
	double oneMinuse2 = 0.993305620009859; // 1-e2
	double oneOver2e = 6.111035746634866; // 1/(2*e)
	double aOver2k0 = 3679336.384424151; // a/(2*k0)
	double sinLat = sin(lat);
//	double a = 6378137;
//	double phi1 = 30*D2R;
//	double k0 = cos(phi1) / sqrt(1 - (e*e*sin(phi1)*sin(phi1)));

//	q = (1-e*e) * ( (sin(lat)/(1-e*e*sin(lat)*sin(lat))) - (1/(2*e))*log((1-e*sin(lat))/(1+e*sin(lat))) );
//	x= a*k0 * lon;
//	y= q * a/(2*k0);

	q = oneMinuse2 * ( (sinLat/(1-e2*sinLat*sinLat)) - oneOver2e*log((1-e*sinLat)/(1+e*sinLat)) );
	x= aXk0 * lon;
	y= q * aOver2k0;

	*col = round(r0 + (x/c)); // as Brodzik et al formula => r [0 1387]
	*row = round(s0 - (y/c)); // as Brodzik et al. formula => s [0 583] +1 included in nl
}

void signalVec(long TxNum)
{
	struct SpecularType *Sp;
	struct SCType *RX, *TX;
	double DvN[3], RvN[3];

	RX = &SC[0];
	TX = &SC[TxNum];
	Sp = &Rx[0].SoOp[TX_MUOS].Sp[0];

	/* ECEF to ECI */
	MTxV(World[EARTH].CWN, Sp->PosW, Sp->PosN);

	// Direct signal vector
	DvN[0] = TX->PosN[0] - RX->PosN[0];
	DvN[1] = TX->PosN[1] - RX->PosN[1];
	DvN[2] = TX->PosN[2] - RX->PosN[2];
	UNITV(DvN);
	MxV(RX->B[0].CN,DvN,Sp->DvB);

	// Reflected signal vector
	RvN[0] = Sp->PosN[0] - RX->PosN[0];
	RvN[1] = Sp->PosN[1] - RX->PosN[1];
	RvN[2] = Sp->PosN[2] - RX->PosN[2];
	UNITV(RvN);
	MxV(RX->B[0].CN,RvN,Sp->RvB);
}

/* Initialize static data of surface */
void InitStaticSurface(void)
{
	FILE * pFile;
	long lSize;
	size_t result;
	double lat, lon;
	uint16_t row, col;

	/* Initialize DEM EARTH2014 1 arc-min resolution */
	pFile = FileOpen(AncillaryPath, SurfaceStt->filename_DEM, "r" );
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread (SurfaceStt->height,2,lSize/2,pFile);
	if (result != lSize/2) {fputs ("DEM file reading error",stderr); exit (3);}
	fclose (pFile);
	// Sanity Check
	lat = 35.1917;
	lon = 172.3417;
	row = 5399.9999 + lat*60;//(90+lat)*60-0.001;
	col = 10799.9999 + lon*60;//(lon+180)*60-0.001;
	if(SurfaceStt->height[row][col]!=65) {
		fputs ("DEM Reading error",stderr); exit (1);
	}
	else
		printf("DEM loaded!\n");

	/* Land Mask - 1km EASE grid 2 */
	pFile = FileOpen(AncillaryPath, SurfaceStt->filename_landMask, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread (SurfaceStt->landMask, 1, lSize, pFile);
	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	// Sanity Check
	lat = 65.415679931640620*D2R;
	lon = 22.432573318481445*D2R;
	LL2EzRC1km(lat, lon, &row, &col);
	if(SurfaceStt->landMask[row][col]!=1){
		printf("Error loading land mask\n");
		exit(1);
	}
	else
		printf("Land mask loaded!\n");

	/*.. Land Type ..*/
	pFile = FileOpen(AncillaryPath, SurfaceStt->filename_landType, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread (SurfaceStt->landType, 1, lSize, pFile);
	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	// Sanity Check
	lat = 11.4499;
	lon = 36.7001;
	row = 1799.9999 - lat*20;//(90-lat)*20;
	col = 3599.9999 + lon*20;//(lon+180)*20;
	if(SurfaceStt->landType[row][col]!=9){
		printf("Error loading land mask\n");
		exit(1);
	}
	else
		printf("Land type loaded!\n");

	/*.. Soil - Sand ..*/
	pFile = FileOpen(AncillaryPath, SurfaceStt->filename_sand, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread(SurfaceStt->sand, sizeof(float), lSize/sizeof(float), pFile);
	if (result != lSize/sizeof(float)) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	printf("Mass fractions of sand loaded!\n");

	/*.. Soil - Clay ..*/
	pFile = FileOpen(AncillaryPath, SurfaceStt->filename_clay, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread(SurfaceStt->clay, sizeof(float), lSize/sizeof(float), pFile);
	if (result != lSize/sizeof(float)) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	printf("Mass fractions of clay loaded!\n");

	/*.. Soil - Bulk Density ..*/
	pFile = FileOpen(AncillaryPath, SurfaceStt->filename_rho_b, "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	// copy the file into the buffer:
	result = fread(SurfaceStt->rho_b, sizeof(float), lSize/sizeof(float), pFile);
	if (result != lSize/sizeof(float)) {fputs ("Reading error",stderr); exit (3);}
	fclose (pFile);
	printf("Bulk density loaded!\n");
}

void InitAntGain(struct SoOpType *SoOp)
{
	switch(SoOp->AntTag){
	case USER_DEFINED:
		antPattern_User(SoOp);
		break;
	case ANT_DIPOLE:
		break;
	case ANT_PATCH:
		break;
	}
//	printf("Ant Gain: %lf\n", SoOp->AntPattern[0][0]);
}

void InitSoOp(void)
{
    FILE *infile;
    char junk[120],newline, response[120];
	long IRx, ITx, IGs, ISoOp, ISp, ITarget, Itheta;
	double K, temp;
	double rE, sinLat2;
    double a = 6378137.0;
    double f = 1.0/298.257222101; // GRS80
    double e2 = f*(2.0-f);
	int row, col;

	char SpriteFileName[40], AntFileName[40];
	double Freq, Bandwidth, CorIntTime, AntMaxGainS, AntMaxGainE, NoiseTemS, NoiseTemE, minSurfRes;
	long NSoOp, NSp, AntTag, Nphi, Ntheta, NthetaLocal, NphiLocal;

/* .. Read from file Inp_SoOp.txt */
    infile=FileOpen(InOutPath,"Inp_SoOp.txt","r");

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    /* .. Dynamic Data .. */
    SurfaceDyn = (struct SurfaceDynType *) calloc(1,sizeof(struct SurfaceDynType));

    // Soil
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_RZSM,junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_soilTemp1,junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    // Snow
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[0],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[1],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[2],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[3],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[4],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[5],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[6],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[7],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[8],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[9],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[10],junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceDyn->filename_snowCover[11],junk,&newline);

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    /* .. Static Data .. */
    // Digital Elevation Model
	SurfaceStt = (struct SurfaceSttType *) calloc(1,sizeof(struct SurfaceSttType));
	fscanf(infile,"%s %[^\n] %[\n]",SurfaceStt->filename_DEM,junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

	// Land Surface Model
	fscanf(infile,"%s %[^\n] %[\n]",SurfaceStt->filename_landMask,junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceStt->filename_landType,junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    // Soil Texture
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceStt->filename_sand,junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceStt->filename_clay,junk,&newline);
    fscanf(infile,"%s %[^\n] %[\n]",SurfaceStt->filename_rho_b,junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    /* .. Target .. */
	fscanf(infile,"%ld %[^\n] %[\n]",&NTarget,junk,&newline);
	Target = (struct TargetAreaType *) calloc(NTarget, sizeof(struct TargetAreaType));
	for (ITarget=0; ITarget<NTarget; ITarget++){
		fscanf(infile,"%lf %lf %lf %[^\n] %[\n]",&Target[ITarget].radius, &Target[ITarget].lat, &Target[ITarget].lon, junk,&newline);
		row = Target[ITarget].lat*60+5399.9999;
		col = Target[ITarget].lon*60+10799.9999;
		Target[ITarget].lat *= D2R;
		Target[ITarget].lon *= D2R;
		sinLat2 = sin(Target[ITarget].lat)*sin(Target[ITarget].lat);
		rE = a*sqrt((1-e2*(2-e2)*sinLat2) / (1-e2*sinLat2));
		Target[ITarget].alt = SurfaceStt->height[row][col] + 6371000 - rE;

		// GRS80 to ECEF
		Target[ITarget].radius *= 1000;
    	GRS80ToECEF(Target[ITarget].lat, Target[ITarget].lon, Target[ITarget].alt, Target[ITarget].PosW);
    	/* ECEF to ECI */
    	MTxV(World[EARTH].CWN, Target[ITarget].PosW, Target[ITarget].PosN);
	}

	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    /* .. Number of Transmitters */
    fscanf(infile,"%ld %[^\n] %[\n]",&NTx,junk,&newline);
    printf("NTx is %ld\n", NTx);
    Tx = (struct TxType *) calloc(NTx,sizeof(struct TxType));
	if (Tx == NULL) {
	 printf("SoOp Tx calloc returned null pointer.  Bailing out!\n");
	 exit(1);
	}
	/* .. Initialization for Transmitters */
	for (ITx=0;ITx<NTx;ITx++){
	    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
		/* .. Tx Frequency .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&temp,junk,&newline); // read in [MHz]
		Tx[ITx].Freq = temp * 1000000; // [MHz] to [Hz]
		/* .. EIRP .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&Tx[ITx].EIRP,junk,&newline);

		Tx[ITx].Wavelength = Clight / Tx[ITx].Freq;
	}

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    /* .. MUOS spot beam radius .. */
    SpotBeam = (struct SpotBeamType *) calloc(1,sizeof(struct SpotBeamType));
	fscanf(infile,"%lf %[^\n] %[\n]",&SpotBeam->scaleFactor,junk,&newline);
	fscanf(infile,"%lf %lf %[^\n] %[\n]",&SpotBeam->halfBW, &SpotBeam->alt, junk,&newline);
	SpotBeam->halfBW *= D2R/2;
	SpotBeam->alt *= 1000;
	InitSpotbeams();

	fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);

    /* .. Number of Receivers */
    fscanf(infile,"%ld %[^\n] %[\n]",&NRx,junk,&newline);
    Rx = (struct RxType *) calloc(NRx,sizeof(struct RxType));
	if (Rx == NULL) {
	 printf("SoOp Rx calloc returned null pointer.  Bailing out!\n");
	 exit(1);
	}
	printf("NRx is %ld\n", NRx);
	/* .. Number of SoOp Configurations .. */
	fscanf(infile,"%ld %[^\n] %[\n]",&NSoOp,junk,&newline);

	/* .. Initialization for Receivers */
	for (IRx=0;IRx<NRx;IRx++){
		/* .. Assign Spacecraft .. */
		Rx[IRx].SC = IRx;
		Rx[IRx].NSoOp = NSoOp;
		Rx[IRx].SoOp = (struct SoOpType *) calloc(NSoOp,sizeof(struct SoOpType));
		if (Rx[IRx].SoOp == NULL) {
		   printf("Rx[%ld].SoOp calloc returned null pointer.  Bailing out!\n", IRx);
		   exit(1);
		}
	}

	for (ISoOp=0;ISoOp<NSoOp;ISoOp++){
		fscanf(infile,"%[^\n] %[\n]",junk,&newline);
		/* .. Number of Specular Points */
		fscanf(infile,"%ld %[^\n] %[\n]",&NSp,junk,&newline);
		/* .. Gain of Gradient Function .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&K,junk,&newline);
		/* .. Sprite File Name for Specular Point .. */
		fscanf(infile,"%s %[^\n] %[\n]",SpriteFileName,junk,&newline);
		/* .. Rx Frequency .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&temp,junk,&newline); // read in [MHz]
		Freq = temp * 1000000; // [MHz] to [Hz]
		/* .. Bandwidth of Channel .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&temp,junk,&newline); // read in [kHz]
		Bandwidth = temp * 1000; // [kHz] to [Hz]
		/* .. Coherent Integration Time .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&CorIntTime,junk,&newline);
		/* .. Antenna Type .. */
		fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline);
		AntTag = DecodeString(response);
		if(AntTag == USER_DEFINED){
			/* .. File Name for Antenna Pattern .. */
			fscanf(infile,"%s %[^\n] %[\n]",AntFileName,junk,&newline);
			/* .. Number of rows (phi-angle Look-up) in Antenna Pattern File .. */
			fscanf(infile,"%ld %[^\n] %[\n]",&Nphi,junk,&newline);
			/* .. Number of columns (theta-angle Look-up) in Antenna Pattern File .. */
			fscanf(infile,"%ld %[^\n] %[\n]",&Ntheta,junk,&newline);
			fscanf(infile,"%[^\n] %[\n]",junk,&newline);
			fscanf(infile,"%[^\n] %[\n]",junk,&newline);
		}
		else{
			fscanf(infile,"%[^\n] %[\n]",junk,&newline);
			fscanf(infile,"%[^\n] %[\n]",junk,&newline);
			fscanf(infile,"%[^\n] %[\n]",junk,&newline);
			/* .. Sky-view Antenna Gain for Direct Signal .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&AntMaxGainS,junk,&newline);
			/* .. Earth-view Antenna Gain for Reflected Signal .. */
			fscanf(infile,"%lf %[^\n] %[\n]",&AntMaxGainE,junk,&newline);
		}
		/* .. Noise Temperature of Sky-view Channel .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&NoiseTemS,junk,&newline);
		/* .. Noise Temperature of Earth-view Channel .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&NoiseTemE,junk,&newline);
		/* .. Minimum Local Surface Resolution .. */
		fscanf(infile,"%lf %[^\n] %[\n]",&minSurfRes,junk,&newline);
		/* .. Number of Grid Points in theta direction .. */
		fscanf(infile,"%ld %[^\n] %[\n]",&NthetaLocal,junk,&newline);
		/* .. Number of Grid Points in phi direction .. */
		fscanf(infile,"%ld %[^\n] %[\n]",&NphiLocal,junk,&newline);

		for (IRx=0;IRx<NRx;IRx++){
			Rx[IRx].SoOp[ISoOp].NSp = NSp;
			Rx[IRx].SoOp[ISoOp].Sp = (struct SpecularType *) calloc(NSp,sizeof(struct SpecularType));
			if (Rx[IRx].SoOp[ISoOp].Sp == NULL) {
			   printf("Rx[%ld].SoOp[%ld].Sp calloc returned null pointer.  Bailing out!\n", IRx, ISoOp);
			   exit(1);
			}
			for(ISp=0; ISp<Rx[IRx].SoOp[ISoOp].NSp; ISp++){
				Rx[IRx].SoOp[ISoOp].Sp[ISp].Soil = (struct SoilType *)calloc(1, sizeof(struct SoilType));
				if (Rx[IRx].SoOp[ISoOp].Sp[ISp].Soil == NULL) {
				   printf("Rx[%ld].SoOp[%ld].Sp[%ld]->Soil calloc returned null pointer.  Bailing out!\n", IRx, ISoOp, ISp);
				   exit(1);
				}
				Rx[IRx].SoOp[ISoOp].Sp[ISp].Reflect = (struct ReflectType *)calloc(1, sizeof(struct ReflectType));
				if (Rx[IRx].SoOp[ISoOp].Sp[ISp].Reflect == NULL) {
				   printf("Rx[%ld].SoOp[%ld].Sp[%ld]->Reflect calloc returned null pointer.  Bailing out!\n", IRx, ISoOp, ISp);
				   exit(1);
				}
			}
			/* .. Sprite File Name for Specular Point .. */
			strcpy(Rx[IRx].SoOp[ISoOp].SpriteFileName, SpriteFileName);
			/* .. Rx Frequency .. */
			Rx[IRx].SoOp[ISoOp].Freq = Freq;
			/* .. Bandwidth of Channel .. */
			Rx[IRx].SoOp[ISoOp].Bandwidth = Bandwidth;
			/* .. Coherent Integration Time .. */
			Rx[IRx].SoOp[ISoOp].CorIntTime = CorIntTime;
			/* .. Antenna Type .. */
			Rx[IRx].SoOp[ISoOp].AntTag = AntTag;

			if(AntTag == USER_DEFINED){
				/* .. File Name for Antenna Pattern .. */
				strcpy(Rx[IRx].SoOp[ISoOp].AntFileName, AntFileName);
				/* .. Number of rows (phi-angle Look-up) in Antenna Pattern File .. */
				Rx[IRx].SoOp[ISoOp].Nphi = Nphi;
				/* .. Number of columns (theta-angle Look-up) in Antenna Pattern File .. */
				Rx[IRx].SoOp[ISoOp].Ntheta = Ntheta;
			}
			else{
				/* .. Sky-view Antenna Gain for Direct Signal .. */
				Rx[IRx].SoOp[ISoOp].AntMaxGainS = AntMaxGainS;
				/* .. Earth-view Antenna Gain for Reflected Signal .. */
				Rx[IRx].SoOp[ISoOp].AntMaxGainE = AntMaxGainE;
			}
			/* .. Noise Temperature of Sky-view Channel .. */
			Rx[IRx].SoOp[ISoOp].NoiseTemS = NoiseTemS;
			/* .. Noise Temperature of Earth-view Channel .. */
			Rx[IRx].SoOp[ISoOp].NoiseTemE = NoiseTemE;
			/* .. Minimum Local Surface Resolution .. */
			Rx[IRx].SoOp[ISoOp].minSurfRes = minSurfRes;
			/* .. Number of Grid Points in theta direction .. */
			Rx[IRx].SoOp[ISoOp].NthetaLocal = NthetaLocal;
			/* .. Number of Grid Points in phi direction .. */
			Rx[IRx].SoOp[ISoOp].NphiLocal = NphiLocal;

			for(ISp=0; ISp<NSp; ISp++){
				Rx[IRx].SoOp[ISoOp].Sp[ISp].K = K;
				Rx[IRx].SoOp[ISoOp].Sp[ISp].Local = (struct LocalSurfaceType **)calloc(NthetaLocal, sizeof(struct LocalSurfaceType));
				for(Itheta=0; Itheta<NthetaLocal; Itheta++)
					Rx[IRx].SoOp[ISoOp].Sp[ISp].Local[Itheta] = (struct LocalSurfaceType *)calloc(NphiLocal, sizeof(struct LocalSurfaceType));
			}
			Rx[IRx].SoOp[ISoOp].Wavelength = Clight / Freq;
//				Rx[IRx].SoOp[ISoOp].NoisePwrS = 10*log10(Rx[IRx].SoOp[ISoOp].NoiseTemS) + 10*log10(Rx[IRx].SoOp[ISoOp].Bandwidth) - 228.6;
//				Rx[IRx].SoOp[ISoOp].NoisePwrE = 10*log10(Rx[IRx].SoOp[ISoOp].NoiseTemE) + 10*log10(Rx[IRx].SoOp[ISoOp].Bandwidth) - 228.6;
			Rx[IRx].SoOp[ISoOp].NoisePwrS = 10*log10(NoiseTemS) - 10*log10(CorIntTime) - 228.6;
			Rx[IRx].SoOp[ISoOp].NoisePwrE = 10*log10(NoiseTemE) - 10*log10(CorIntTime) - 228.6;
			InitAntGain(&Rx[IRx].SoOp[ISoOp]);
		}
	}

    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    fscanf(infile,"%[^\n] %[\n]",junk,&newline);
    /* .. Number of Ground Station */
	fscanf(infile,"%ld %[^\n] %[\n]",&NGs,junk,&newline);
	printf("NGs is %ld\n", NGs);
	Gs = (struct GsType *) calloc(NGs,sizeof(struct GsType));
	if (Gs == NULL) {
	 printf("SoOp Rx calloc returned null pointer.  Bailing out!\n");
	 exit(1);
	}
	if (NGs == 0) {
	   for(IGs=0;IGs<33;IGs++) fscanf(infile,"%[^\n] %[\n]",junk,&newline);
	}
	else {
	/* .. Initialization for Ground Station */
		for (IGs=0;IGs<NGs;IGs++){
			fscanf(infile,"%[^\n] %[\n]",junk,&newline);
			/* .. Assign Ground Station .. */
			fscanf(infile,"%ld %[^\n] %[\n]",&Gs[IGs].GroundStation,junk,&newline);
			/* .. Number of SoOp Configurations .. */
			fscanf(infile,"%ld %[^\n] %[\n]",&Gs[IGs].NSoOp,junk,&newline);
			Gs[IGs].SoOp = (struct SoOpType *) calloc(Gs[IGs].NSoOp,sizeof(struct SoOpType));
			if (Gs[IGs].SoOp == NULL) {
			   printf("Gs[%ld].SoOp calloc returned null pointer.  Bailing out!\n", IGs);
			   exit(1);
			}
			for (ISoOp=0;ISoOp<Gs[IGs].NSoOp;ISoOp++){
				fscanf(infile,"%[^\n] %[\n]",junk,&newline);
				/* .. Number of Specular Points */
				fscanf(infile,"%ld %[^\n] %[\n]",&Gs[IGs].SoOp[ISoOp].NSp,junk,&newline);
				Gs[IGs].SoOp[ISoOp].Sp = (struct SpecularType *) calloc(Gs[IGs].SoOp[ISoOp].NSp,sizeof(struct SpecularType));
				if (Gs[IGs].SoOp[ISoOp].Sp == NULL) {
				   printf("Gs[%ld].SoOp[%ld].Sp calloc returned null pointer.  Bailing out!\n", IGs, ISoOp);
				   exit(1);
				}
				/* .. Gain of Gradient Function .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&K,junk,&newline);
				/* .. Sprite File Name for Specular Point .. */
				fscanf(infile,"%s %[^\n] %[\n]",Gs[IGs].SoOp[ISoOp].SpriteFileName,junk,&newline);
				/* .. Rx Frequency .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&temp,junk,&newline); // read in [MHz]
				Gs[IGs].SoOp[ISoOp].Freq = temp * 1000000; // [MHz] to [Hz]
				/* .. Bandwidth of Channel .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&temp,junk,&newline); // read in [kHz]
				Gs[IGs].SoOp[ISoOp].Bandwidth = temp * 1000; // [kHz] to [Hz]
				/* .. Coherent Integration Time .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&Rx[IRx].SoOp[ISoOp].CorIntTime,junk,&newline);
				/* .. Sky-view Antenna Gain for Direct Signal .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&Gs[IGs].SoOp[ISoOp].AntMaxGainS,junk,&newline);
				/* .. Earth-view Antenna Gain for Reflected Signal .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&Gs[IGs].SoOp[ISoOp].AntMaxGainE,junk,&newline);
				/* .. Noise Temperature of Sky-view Channel .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&Gs[IGs].SoOp[ISoOp].NoiseTemS,junk,&newline);
				/* .. Noise Temperature of Earth-view Channel .. */
				fscanf(infile,"%lf %[^\n] %[\n]",&Gs[IGs].SoOp[ISoOp].NoiseTemE,junk,&newline);

				for(ISp=0; ISp<Gs[IGs].SoOp[ISoOp].NSp; ISp++)
					Gs[IGs].SoOp[ISoOp].Sp[ISp].K = K;
				Gs[IGs].SoOp[ISoOp].Wavelength = Clight / Gs[IGs].SoOp[ISoOp].Freq;
				Gs[IGs].SoOp[ISoOp].NoisePwrS = 10*log10(Gs[IGs].SoOp[ISoOp].NoiseTemS) + 10*log10(Gs[IGs].SoOp[ISoOp].Bandwidth) - 228.6;
				Gs[IGs].SoOp[ISoOp].NoisePwrE = 10*log10(Gs[IGs].SoOp[ISoOp].NoiseTemE) + 10*log10(Gs[IGs].SoOp[ISoOp].Bandwidth) - 228.6;
			}
		}
	}
	fclose(infile);

	/* .. Initialize Surface Static Data .. */
	InitStaticSurface();
    /* .. Initialize Soil Dynamic Data .. */
	LoadSoilData();
}

void RxFsw(struct SCType *S)
{
	struct AcType *AC;
	struct AcCfsCtrlType *C;
	struct AcJointType *G;
    double L1[3],L2[3],L3[3];
    double HxB[3];
//    double Hb[3];
	long i, j;

	AC = &S->AC;
	C = &AC->CfsCtrl;
	G = &AC->G[0];

	if (C->Init) {
		C->Init = 0;
		for(i=0;i<3;i++) FindPDGains(AC->MOI[i][i],0.1,0.7,&C->Kr[i],&C->Kp[i]);
		C->Kunl = 1.0E6;
		FindPDGains(100.0,0.2,1.0,&G->AngRateGain[0],&G->AngGain[0]);
		G->MaxAngRate[0] = 1.0*D2R;
		G->MaxTrq[0] = 10.0;
	}

	/* .. Sensor Processing */
	GyroProcessing(AC);
	MagnetometerProcessing(AC);
	CssProcessing(AC);
	FssProcessing(AC);
	StarTrackerProcessing(AC);
	GpsProcessing(AC);

	/* .. Commanded Attitude */
	      if (AC->GPS[0].Valid) {
	         CopyUnitV(AC->PosN,L3);
	         VxV(AC->PosN,AC->VelN,L2);
	         UNITV(L2);
	         UNITV(L3);
	         for(i=0;i<3;i++) {
	            L2[i] = -L2[i];
	            L3[i] = -L3[i];
	         }
	         VxV(L2,L3,L1);
	         UNITV(L1);
	         for(i=0;i<3;i++) {
	            AC->CLN[0][i] = L1[i];
	            AC->CLN[1][i] = L2[i];
	            AC->CLN[2][i] = L3[i];
	         }
	         C2Q(AC->CLN,AC->qln);
	         AC->wln[1] = -MAGV(AC->VelN)/MAGV(AC->PosN);
	      }
	      else {
	         for(i=0;i<3;i++) {
	            for(j=0;j<3;j++) {
	               AC->CLN[i][j] = 0.0;
	            }
	            AC->CLN[i][i] = 1.0;
	            AC->qln[i] = 0.0;
	            AC->wln[i] = 0.0;
	         }
	         AC->qln[3] = 1.0;
	      }

	/* .. Attitude Control */
	      if (AC->StValid) {
	         QxQT(AC->qbn,AC->qln,AC->qbr);
	         RECTIFYQ(AC->qbr);
	      }
	      else {
	         for(i=0;i<3;i++) AC->qbr[i] = 0.0;
	         AC->qbr[3] = 1.0;
	      }
	      for(i=0;i<3;i++) {
	         C->therr[i] = Limit(2.0*AC->qbr[i],-0.1,0.1);
	         C->werr[i] = AC->wbn[i] - AC->wln[i];
	         AC->Tcmd[i] = Limit(-C->Kr[i]*C->werr[i] - C->Kp[i]*C->therr[i],-0.1,0.1);
	      }
	/* .. Momentum Management */
	      for(i=0;i<3;i++) {
	         AC->Hvb[i] = AC->MOI[i][i]*AC->wbn[i];
	         for(j=0;j<AC->Nwhl;j++) AC->Hvb[i] += AC->Whl[j].Axis[i]*AC->Whl[j].H;
	      }
	      VxV(AC->Hvb,AC->bvb,HxB);
	      for(i=0;i<3;i++) AC->Mcmd[i] = C->Kunl*HxB[i];

/* .. Commanded Attitude */
	/* Find qrn, wrn and joint angle commands */
//    ThreeAxisAttitudeCommand(S);

    /* .. Attitude Control */
    /* Form Error Signals */
/*    QxQT(AC->qbn,AC->Cmd.qrn,AC->qbr);
    RECTIFYQ(AC->qbr);

	for(i=0;i<3;i++) {
		C->therr[i] = Limit(2.0*AC->qbr[i],-0.05,0.05);
		C->werr[i] = AC->wbn[i] - AC->Cmd.wrn[i];
		AC->Tcmd[i] = Limit(-C->Kr[i]*C->werr[i] - C->Kp[i]*C->therr[i],-0.1,0.1);
	}
*/
	/* .. Momentum Management */
//	for(i=0;i<3;i++) Hb[i] = AC->MOI[i][i]*AC->wbn[i] + AC->Whl[i].H;
//	VxV(Hb,AC->bvb,HxB);
//	for(i=0;i<3;i++) AC->Mcmd[i] = C->Kunl*HxB[i];

	/* .. Solar Array Steering */
//	G->Cmd.Ang[0] = atan2(AC->svb[0],AC->svb[2]);

	/* .. Actuator Processing */
	WheelProcessing(AC);
	MtbProcessing(AC);
}

void postProcessCalVal(uint8_t RxNum, uint8_t TxNum)
{
    FILE *infile, *outfile;
	uint16_t arc[2000] = {0}, iArc = 0, NArc;
	uint32_t iLine = 0, NLine = 0;

	uint8_t flagValid = 0, flagWater =0 , flagVisible = 0;

	char filename[30], chr;
	char TxName[4][10] = {"MUOS", "Orbcomm", "NOAA", "METEOR"};

	sprintf(filename, "CalVal_%s_Case%d.42", TxName[TxNum], RxNum);
	infile=FileOpen(InOutPath,filename,"r");
	chr = getc(infile);
	while(chr!=EOF){
		if(chr=='\n')
			NLine++;
		chr = getc(infile);
	}
	printf("# lines of %s = %d\n", filename, NLine);
	rewind(infile);

    Data = (struct PostProcType *) calloc(NLine,sizeof(struct PostProcType));

	fscanf(infile, "%hu %lf %hu %hu %hu %lf %lf %lf %hu %f %f %f %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&Data[iLine].arc, &Data[iLine].t, &Data[iLine].TxID, &Data[iLine].site, &Data[iLine].visible,
			&Data[iLine].lat, &Data[iLine].lon, &Data[iLine].alt,
			&Data[iLine].landMask, &Data[iLine].m_v, &Data[iLine].Temp, &Data[iLine].sand, &Data[iLine].clay, &Data[iLine].rho_b,
			&Data[iLine].gammaLR, &Data[iLine].SnrD, &Data[iLine].SnrR,
			&Data[iLine].a, &Data[iLine].b, &Data[iLine].az,
			&Data[iLine].lookAngAzD, &Data[iLine].lookAngElD, &Data[iLine].lookAngAzR, &Data[iLine].lookAngElR);
	for(iLine = 1; iLine<NLine; iLine++){
		fscanf(infile, "%hu %lf %hu %hu %hu %lf %lf %lf %hu %f %f %f %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
				&Data[iLine].arc, &Data[iLine].t, &Data[iLine].TxID, &Data[iLine].site, &Data[iLine].visible,
				&Data[iLine].lat, &Data[iLine].lon, &Data[iLine].alt,
				&Data[iLine].landMask, &Data[iLine].m_v, &Data[iLine].Temp, &Data[iLine].sand, &Data[iLine].clay, &Data[iLine].rho_b,
				&Data[iLine].gammaLR, &Data[iLine].SnrD, &Data[iLine].SnrR,
				&Data[iLine].a, &Data[iLine].b, &Data[iLine].az,
				&Data[iLine].lookAngAzD, &Data[iLine].lookAngElD, &Data[iLine].lookAngAzR, &Data[iLine].lookAngElR);

		if(Data[iLine].arc > Data[iLine-1].arc){
			flagValid = 0;
			flagVisible = 0;
			flagWater = 0;
		}

		if(Data[iLine].visible)
			flagVisible = 1;

		if(Data[iLine].landMask==0)
			flagWater = 1;

		if(flagVisible && flagWater && flagValid==0){
			arc[iArc++] = Data[iLine].arc;
			flagValid = 1;
			flagVisible = 0;
			flagWater = 0;
		}
	}
	NArc = iArc;
	iArc = 0;

	printf("# valid arcs of %s = %d\n", filename, NArc);

	strcat(filename,".arc");
	outfile = FileOpen(InOutPath,filename,"w");

	for(iLine = 1; iLine<NLine; iLine++){
		if(Data[iLine].arc > arc[iArc] && iArc < NArc){
			iArc++;
		}
		if(Data[iLine].arc == arc[iArc]){
			fprintf(outfile, "%hu %lf %hu %hu %hu %lf %lf %lf %hu %f %f %f %f %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					Data[iLine].arc, Data[iLine].t, Data[iLine].TxID, Data[iLine].site, Data[iLine].visible,
					Data[iLine].lat, Data[iLine].lon, Data[iLine].alt,
					Data[iLine].landMask, Data[iLine].m_v, Data[iLine].Temp, Data[iLine].sand, Data[iLine].clay, Data[iLine].rho_b,
					Data[iLine].gammaLR, Data[iLine].SnrD, Data[iLine].SnrR,
					Data[iLine].a, Data[iLine].b, Data[iLine].az,
					Data[iLine].lookAngAzD, Data[iLine].lookAngElD, Data[iLine].lookAngAzR, Data[iLine].lookAngElR);
		}
	}
	printf("%s saved!\n", filename);
}
