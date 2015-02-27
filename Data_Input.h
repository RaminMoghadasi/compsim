#ifndef DATA_INPUT_H
#define DATA_INPUT_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#include "Globals.h"
#include "Init_Calcs.h"
#include "Do_Cycle.h"


extern int Nx, Ny, Nz;		//Reservoir Dimension
extern FType *gridDim;
extern FType ***porosity;		//Porosity Distribution
extern FType ****perm;		//Permeability tensor distribution
extern FType refP, cpor, dcpor;		//Compressibility
extern int PNc, UNc, Nc;		//Number of Components
extern unsigned char PR, SRK;		//EOS Type, boolean
extern FType **fluidProp;
extern int Nswt, Nsgt;
extern FType **swt, **sgt;
extern int initCond;		//Initial Condition State
extern FType ****sat, ****bsat;	//Saturation
extern FType ****P;		//Pressure
extern FType *****comp;		//Phases Composition
extern int refL;
extern int wellNO;		//Number of wells
extern int **welli;
extern FType **wellf;
extern FType ****IFTran;
extern FType resTemp;
extern FType *****blockFProps;
extern FType *****trans;
extern FType ****Wtran;
extern FType ****relPerm;
extern FType watRo, watMu;
//extern FType **jac;
extern FType *ans;
extern FType ****preProp;
extern FType *****dbMu;
extern FType *****dE;
extern char *****transS;
extern FType ****dSatProps;
extern char ***phaseStat;
extern FType *Unk;
extern FType WOCHeight;
extern FType totalTime;
extern FType ***blockH;
extern FType *****bcomp;
extern FType ****bP;
extern FType *CSRjac;
extern int *CSRrow, *CSRcol;
extern int ***pJHolder;
extern FType *preCon;
extern int *preConRow;
extern int *preConIndex;
extern FType **bic;
extern FType *TStepMarker;
extern FType **STcomp;
extern char ***bphaseStat;
extern int TSMCtrl;
extern FType ****dRelPerm;
extern FType ****dWtran;
//extern FType **fullPre;
extern int *preConCSRrow;
extern int *preConCSRcol;
extern FType ***tor;
extern FType *****diffusion;
extern FType ***Bift;



int File_Search(FILE *, char *);
int Read_Word(FILE *, char *);
void Allocation(void);
void ECLStar(char *, int *, FType *);
void RestoreRST(FType *);

void Data_Input(FILE *fp) {
	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
	register int i, j, k, n;
	FType tempL;
	int eclint;
	
	//Dimension
	if (!File_Search(fp, "GRID")) TerM("No GRID keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect GRID keyword format in the input file!");
	Nx=atoi(str);
	if (!Read_Word(fp, str)) TerM("Incorrect GRID keyword format in the input file!");
	Ny=atoi(str);
	if (!Read_Word(fp, str)) TerM("Incorrect GRID keyword format in the input file!");
	Nz=atoi(str);

	//Number of Components
	if (!File_Search(fp, "NC")) TerM("No NC keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect NC keyword format in the input file!");
	PNc=atoi(str);
	if (!Read_Word(fp, str)) TerM("Incorrect NC keyword format in the input file!");
	UNc=atoi(str);
	Nc=PNc+UNc;


	//Saturation Tables
	if (!File_Search(fp, "SWT")) TerM("No SWT keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect SWT keyword format in the input file!");
	Nswt=atoi(str);

	if (!File_Search(fp, "SGT")) TerM("No SGT keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect SGT keyword format in the input file!");
	Nsgt=atoi(str);

	//Wells
	if (!File_Search(fp, "WELLS")) TerM("No WELLS keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect WELLS keyword format in the input file!");
	wellNO=atoi(str);
	

	Allocation();
	

	//Reservoir Temperature
	if (!File_Search(fp, "RESTEMP")) TerM("No RESTEMP keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect RESTEMP keyword format in the input file!");
	resTemp=atof(str)+273.15;


	//Water Properties
	if (!File_Search(fp, "WATERPROPS")) TerM("No WATERPROPS keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect WATERPROPS keyword format in the input file!");
	watRo=atof(str);
	if (!Read_Word(fp, str)) TerM("Incorrect WATERPROPS keyword format in the input file!");
	watMu=atof(str);


	//////////block sizes (finite difference)//////////////
	if (!File_Search(fp, "DI")) TerM("No DI keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect DI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Nx; i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect DI keyword format in the input file!");
		gridDim[i]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect DI keyword format in the input file!");
		tempL=atof(str1);
		for (i=0; i<Nx; i++) gridDim[i]=tempL;
	}
	else {
		TerM("Incorrect DI keyword format in the input file!");
	}

	if (!File_Search(fp, "DJ")) TerM("No DJ keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect DJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=Nx; i<(Nx+Ny); i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect DJ keyword format in the input file!");
		gridDim[i]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect DJ keyword format in the input file!");
		tempL=atof(str1);
		for (i=Nx; i<(Nx+Ny); i++) gridDim[i]=tempL;
	}
	else {
		TerM("Incorrect DJ keyword format in the input file!");
	}

	if (!File_Search(fp, "DK")) TerM("No DK keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect DK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=(Nx+Ny); i<(Nx+Ny+Nz); i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect DK keyword format in the input file!");
		gridDim[i]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect DK keyword format in the input file!");
		tempL=atof(str1);
		for (i=(Nx+Ny); i<(Nx+Ny+Nz); i++) gridDim[i]=tempL;
	}
	else {
		TerM("Incorrect DK keyword format in the input file!");
	}

	CalcBlockHeight();


	//Tortuosity
	if (!File_Search(fp, "TOR")) {
		puts("All tortuosity values reset to unity!");
		for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) tor[i][j][k]=1;
	}
	else {
		if (!Read_Word(fp, str)) TerM("Incorrect TOR keyword format in the input file!");
		if (!strcmp(str, "VAR")) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) { 
			if (!Read_Word(fp, str1)) TerM("Incorrect TOR keyword format in the input file!");
			tor[i][j][k]=atof(str1);		
		}
		else if (!strcmp(str, "CON")){
			if (!Read_Word(fp, str1)) TerM("Incorrect TOR keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) tor[i][j][k]=tempL;
		}
		else if (!strcmp(str, "IVAR")){
			for (i=0; i<Nx; i++) {
				if (!Read_Word(fp, str1)) TerM("Incorrect TOR keyword format in the input file!");
				tempL=atof(str1);
				for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) tor[i][j][k]=tempL;
			}
		}
		else if (!strcmp(str, "JVAR")){
			for (j=0; i<Ny; j++) {
				if (!Read_Word(fp, str1)) TerM("Incorrect TOR keyword format in the input file!");
				tempL=atof(str1);
				for (k=0; k<Nz; k++) for (i=0; i<Nx; i++) tor[i][j][k]=tempL;
			}
		}
		else if (!strcmp(str, "KVAR")){
			for (k=0; k<Nz; k++) {
				if (!Read_Word(fp, str1)) TerM("Incorrect TOR keyword format in the input file!");
				tempL=atof(str1);
				for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) tor[i][j][k]=tempL;
			}
		}
		else {
			TerM("Incorrect TOR keyword format in the input file!");
		}
	}
	

	//POROSITY
	if (!File_Search(fp, "POR")) TerM("No POR keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect POR keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect POR keyword format in the input file!");
		porosity[i][j][k]=atof(str1);
		if (!porosity[i][j][k]) porosity[i][j][k]=1e-5;
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect POR keyword format in the input file!");
		tempL=atof(str1);
		for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) porosity[i][j][k]=tempL;
	}
	else if (!strcmp(str, "IVAR")){
		for (i=0; i<Nx; i++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect POR keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) porosity[i][j][k]=tempL;
		}
	}
	else if (!strcmp(str, "JVAR")){
		for (j=0; i<Ny; j++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect POR keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (i=0; i<Nx; i++) porosity[i][j][k]=tempL;
		}
	}
	else if (!strcmp(str, "KVAR")){
		for (k=0; k<Nz; k++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect POR keyword format in the input file!");
			tempL=atof(str1);
			for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) porosity[i][j][k]=tempL;
		}
	}
	else {
		TerM("Incorrect POR keyword format in the input file!");
	}

	//PERMEABILITY
	if (!File_Search(fp, "PERMI")) TerM("No PERMI keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect PERMI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect PERMI keyword format in the input file!");
		perm[i][j][k][0]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect PERMI keyword format in the input file!");
		tempL=atof(str1);
		for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) perm[i][j][k][0]=tempL;
	}
	else if (!strcmp(str, "IVAR")){
		for (i=0; i<Nx; i++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMI keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) perm[i][j][k][0]=tempL;
		}
	}
	else if (!strcmp(str, "JVAR")){
		for (j=0; i<Ny; j++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMI keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (i=0; i<Nx; i++) perm[i][j][k][0]=tempL;
		}
	}
	else if (!strcmp(str, "KVAR")){
		for (k=0; k<Nz; k++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMI keyword format in the input file!");
			tempL=atof(str1);
			for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) perm[i][j][k][0]=tempL;
		}
	}
	else {
		TerM("Incorrect PERMI keyword format in the input file!");
	}

	if (!File_Search(fp, "PERMJ")) TerM("No PERMJ keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect PERMJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect PERMI keyword format in the input file!");
		perm[i][j][k][1]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect PERJ keyword format in the input file!");
		tempL=atof(str1);
		for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) perm[i][j][k][1]=tempL;
	}
	else if (!strcmp(str, "IVAR")){
		for (i=0; i<Nx; i++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMJ keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) perm[i][j][k][1]=tempL;
		}
	}
	else if (!strcmp(str, "JVAR")){
		for (j=0; i<Ny; j++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMJ keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (i=0; i<Nx; i++) perm[i][j][k][1]=tempL;
		}
	}
	else if (!strcmp(str, "KVAR")){
		for (k=0; k<Nz; k++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMJ keyword format in the input file!");
			tempL=atof(str1);
			for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) perm[i][j][k][1]=tempL;
		}
	}
	else {
		TerM("Incorrect PERMJ keyword format in the input file!");
	}

	if (!File_Search(fp, "PERMK")) TerM("No PERMK keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect PERMK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) { 
		if (!Read_Word(fp, str1)) TerM("Incorrect PERMK keyword format in the input file!");
		perm[i][j][k][2]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("Incorrect PERMK keyword format in the input file!");
		tempL=atof(str1);
		for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) perm[i][j][k][2]=tempL;
	}
	else if (!strcmp(str, "IVAR")){
		for (i=0; i<Nx; i++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMK keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) perm[i][j][k][2]=tempL;
		}
	}
	else if (!strcmp(str, "JVAR")){
		for (j=0; i<Ny; j++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMK keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (i=0; i<Nx; i++) perm[i][j][k][2]=tempL;
		}
	}
	else if (!strcmp(str, "KVAR")){
		for (k=0; k<Nz; k++) {
			if (!Read_Word(fp, str1)) TerM("Incorrect PERMK keyword format in the input file!");
			tempL=atof(str1);
			for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) perm[i][j][k][2]=tempL;
		}
	}
	else {
		TerM("Incorrect PERMK keyword format in the input file!");
	}

	//Compressibility
	if (!File_Search(fp, "REFP")) {
		puts("Warning: No REFP keyword in the input file!");
		refP=0;
	}
	else {
		if (!Read_Word(fp, str)) TerM("Incorrect REFP keyword format in the input file!");
		refP=atof(str);
	}
	if (!File_Search(fp, "CPOR")) {
		puts("Warning: No CPOR keyword in the input file!");
		cpor=0;
	}
	else {
		if (!Read_Word(fp, str)) TerM("Incorrect CPOR keyword format in the input file!");
		cpor=atof(str);
	}
	if (!File_Search(fp, "DCPOR")) {
		puts("Warning: No DCPOR keyword in the input file!");
		dcpor=0;
	}
	else {
		if (!Read_Word(fp, str)) TerM("Incorrect DCPOR keyword format in the input file!");
		dcpor=atof(str);
	}

	//EOS Type
	if (!File_Search(fp, "MODEL")) TerM("No EOS MODEL keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect EOS MODEL keyword format in the input file!");
	if (!strcmp(str, "PR")) {
		PR=1;
		SRK=0;
	}
	else if (!strcmp(str, "SRK")) {
		PR=0;
		SRK=1;
	}
	else{
		TerM("Incorrect EOS keyword format in the input file!");
	}

	if (UNc) {
		if (!File_Search(fp, "SG")) puts("Warning: No SG keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect SG keyword format in the input file!");
			fluidProp[i][SG]=atof(str);
		}

		if (!File_Search(fp, "TB")) puts("Warning: No TB keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect TB keyword format in the input file!");
			fluidProp[i][TB]=atof(str);
		}

		if (!File_Search(fp, "MW")) puts("Warning: No MW keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect MW keyword format in the input file!");
			fluidProp[i][MW]=atof(str);
		}

		if (!File_Search(fp, "AC")) puts("Warning: No AC keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect AC keyword format in the input file!");
			fluidProp[i][AC]=atof(str);
		}

		if (!File_Search(fp, "PCRIT")) puts("Warning: No PCRIT keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect PCRIT keyword format in the input file!");
			fluidProp[i][PCRIT]=atof(str);
		}

		if (!File_Search(fp, "TCRIT")) puts("Warning: No TCRIT keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect TCRIT keyword format in the input file!");
			fluidProp[i][TCRIT]=atof(str);
		}

		if (!File_Search(fp, "VCRIT")) puts("Warning: No VCRIT keyword in the input file!");
		else for (i=PNc; i<Nc; i++) {
			if (!Read_Word(fp, str)) TerM("Incorrect VCRIT keyword format in the input file!");
			fluidProp[i][VCRIT]=atof(str);
		}

	}
	CPlus_Props();

	if (!File_Search(fp, "COMPNAME")) TerM("No COMPNAME keyword in the input file!");
	for (i=0; i<PNc; i++) {
		if (!Read_Word(fp, str)) TerM("Incorrect COMPNAME keyword format in the input file!");
		if (!strcmp(str, "C1")) {
			fluidProp[i][MW]=16.043;
			fluidProp[i][TCRIT]=190.56;
			fluidProp[i][PCRIT]=4599000;
			fluidProp[i][VCRIT]=98.6;
			fluidProp[i][AC]=0.008;	//edited
			fluidProp[i][PARACHOR]=74.05;

		}
		else if (!strcmp(str, "C2")) {
			fluidProp[i][MW]=30.070;
			fluidProp[i][TCRIT]=305.32;
			fluidProp[i][PCRIT]=4872000;
			fluidProp[i][VCRIT]=145.5;
			fluidProp[i][AC]=0.0995;
			fluidProp[i][PARACHOR]=112.91;
		}
		else if (!strcmp(str, "C3")) {
			fluidProp[i][MW]=44.096;
			fluidProp[i][TCRIT]=369.83;
			fluidProp[i][PCRIT]=4248000;
			fluidProp[i][VCRIT]=200;
			fluidProp[i][AC]=0.1523;
			fluidProp[i][PARACHOR]=154.03;
		}
		else if (!strcmp(str, "iC4")) {
			fluidProp[i][MW]=58.123;
			fluidProp[i][TCRIT]=408.14;
			fluidProp[i][PCRIT]=3648000;
			fluidProp[i][VCRIT]=262.7;
			fluidProp[i][AC]=0.1770;
			fluidProp[i][PARACHOR]=185.32;
		}
		else if (!strcmp(str, "nC4")) {
			fluidProp[i][MW]=58.123;
			fluidProp[i][TCRIT]=425.12;
			fluidProp[i][PCRIT]=3796000;
			fluidProp[i][VCRIT]=255;
			fluidProp[i][AC]=0.2002;
			fluidProp[i][PARACHOR]=193.90;
		}
		else if (!strcmp(str, "iC5")) {
			fluidProp[i][MW]=72.150;
			fluidProp[i][TCRIT]=460.43;
			fluidProp[i][PCRIT]=3381000;
			fluidProp[i][VCRIT]=305.8;
			fluidProp[i][AC]=0.2275;
			fluidProp[i][PARACHOR]=229.37;
		}
		else if (!strcmp(str, "nC5")) {
			fluidProp[i][MW]=72.150;
			fluidProp[i][TCRIT]=469.7;
			fluidProp[i][PCRIT]=3370000;
			fluidProp[i][VCRIT]=313;
			fluidProp[i][AC]=0.2515;
			fluidProp[i][PARACHOR]=236.00;
		}
		else if (!strcmp(str, "nC6")) {
			fluidProp[i][MW]=86.177;
			fluidProp[i][TCRIT]=507.6;
			fluidProp[i][PCRIT]=3025000;
			fluidProp[i][VCRIT]=371;
			fluidProp[i][AC]=0.3013;
			fluidProp[i][PARACHOR]=276.71;
		}
		else if (!strcmp(str, "nC7")) {
			fluidProp[i][MW]=100.204;
			fluidProp[i][TCRIT]=540.2;
			fluidProp[i][PCRIT]=2740000;
			fluidProp[i][VCRIT]=428;
			fluidProp[i][AC]=0.3495;
			fluidProp[i][PARACHOR]=318.44;
		}
		else if (!strcmp(str, "nC8")) {
			fluidProp[i][MW]=114.231;
			fluidProp[i][TCRIT]=568.7;
			fluidProp[i][PCRIT]=2490000;
			fluidProp[i][VCRIT]=486;
			fluidProp[i][AC]=0.3996;
			fluidProp[i][PARACHOR]=359.33;
		}
		else if (!strcmp(str, "nC9")) {
			fluidProp[i][MW]=128.258;
			fluidProp[i][TCRIT]=594.6;
			fluidProp[i][PCRIT]=2290000;
			fluidProp[i][VCRIT]=544;
			fluidProp[i][AC]=0.4435;
			fluidProp[i][PARACHOR]=399.57;
		}
		else if (!strcmp(str, "nC10")) {
			fluidProp[i][MW]=142.285;
			fluidProp[i][TCRIT]=617.7;
			fluidProp[i][PCRIT]=2110000;
			fluidProp[i][VCRIT]=600;
			fluidProp[i][AC]=0.4923;
			fluidProp[i][PARACHOR]=440.69;
		}
		else if (!strcmp(str, "nC11")) {
			fluidProp[i][MW]=156.312;
			fluidProp[i][TCRIT]=639;
			fluidProp[i][PCRIT]=1949000;
			fluidProp[i][VCRIT]=659;
			fluidProp[i][AC]=0.5303;
			fluidProp[i][PARACHOR]=482.00;
		}
		else if (!strcmp(str, "nC12")) {
			fluidProp[i][MW]=170.338;
			fluidProp[i][TCRIT]=658;
			fluidProp[i][PCRIT]=1820000;
			fluidProp[i][VCRIT]=716;
			fluidProp[i][AC]=0.5764;
			fluidProp[i][PARACHOR]=522.26;
		}
		else if (!strcmp(str, "nC13")) {
			fluidProp[i][MW]=184.365;
			fluidProp[i][TCRIT]=675;
			fluidProp[i][PCRIT]=1680000;
			fluidProp[i][VCRIT]=775;
			fluidProp[i][AC]=0.6174;
			fluidProp[i][PARACHOR]=536.77;
		}
		else if (!strcmp(str, "nC14")) {
			fluidProp[i][MW]=198.392;
			fluidProp[i][TCRIT]=693;
			fluidProp[i][PCRIT]=1570000;
			fluidProp[i][VCRIT]=830;
			fluidProp[i][AC]=0.6430;
			fluidProp[i][PARACHOR]=606.05;
		}
		else if (!strcmp(str, "nC15")) {
			fluidProp[i][MW]=212.419;
			fluidProp[i][TCRIT]=708;
			fluidProp[i][PCRIT]=1480000;
			fluidProp[i][VCRIT]=889;
			fluidProp[i][AC]=0.6863;
			fluidProp[i][PARACHOR]=647.43;
		}
		else if (!strcmp(str, "nC16")) {
			fluidProp[i][MW]=226.446;
			fluidProp[i][TCRIT]=723;
			fluidProp[i][PCRIT]=1400000;
			fluidProp[i][VCRIT]=944;
			fluidProp[i][AC]=0.7174;
			fluidProp[i][PARACHOR]=688.50;
		}
		else if (!strcmp(str, "nC17")) {
			fluidProp[i][MW]=240.473;
			fluidProp[i][TCRIT]=736;
			fluidProp[i][PCRIT]=1340000;
			fluidProp[i][VCRIT]=1000;
			fluidProp[i][AC]=0.7697;
			fluidProp[i][PARACHOR]=730.05;
		}
		else if (!strcmp(str, "nC18")) {
			fluidProp[i][MW]=254.5;
			fluidProp[i][TCRIT]=747;
			fluidProp[i][PCRIT]=1270000;
			fluidProp[i][VCRIT]=1060;
			fluidProp[i][AC]=0.8114;
			fluidProp[i][PARACHOR]=771.95;
		}
		else if (!strcmp(str, "nC19")) {
			fluidProp[i][MW]=268.527;
			fluidProp[i][TCRIT]=758;
			fluidProp[i][PCRIT]=1210000;
			fluidProp[i][VCRIT]=1120;
			fluidProp[i][AC]=0.8522;
			fluidProp[i][PARACHOR]=813.85;
		}
		else if (!strcmp(str, "nC20")) {
			fluidProp[i][MW]=282.553;
			fluidProp[i][TCRIT]=768;
			fluidProp[i][PCRIT]=1160000;
			fluidProp[i][VCRIT]=1170;
			fluidProp[i][AC]=0.9069;
			fluidProp[i][PARACHOR]=853.67;
		}
		else if (!strcmp(str, "nC21")) {
			fluidProp[i][MW]=296.580;
			fluidProp[i][TCRIT]=781.7;
			fluidProp[i][PCRIT]=1147000;
			fluidProp[i][VCRIT]=1198;
			fluidProp[i][AC]=0.9220;
			fluidProp[i][PARACHOR]=897.64;
		}
		else if (!strcmp(str, "nC22")) {
			fluidProp[i][MW]=310.610;
			fluidProp[i][TCRIT]=791.8;
			fluidProp[i][PCRIT]=1101000;
			fluidProp[i][VCRIT]=1253;
			fluidProp[i][AC]=0.9550;
			fluidProp[i][PARACHOR]=939.55;
		}
		else if (!strcmp(str, "nC23")) {
			fluidProp[i][MW]=324.630;
			fluidProp[i][TCRIT]=801.3;
			fluidProp[i][PCRIT]=1059000;
			fluidProp[i][VCRIT]=1307;
			fluidProp[i][AC]=0.9890;
			fluidProp[i][PARACHOR]=981.43;
		}
		else if (!strcmp(str, "nC24")) {
			fluidProp[i][MW]=338.680;
			fluidProp[i][TCRIT]=810.4;
			fluidProp[i][PCRIT]=1019000;
			fluidProp[i][VCRIT]=1362;
			fluidProp[i][AC]=1.0190;
			fluidProp[i][PARACHOR]=1023.40;
		}
		else if (!strcmp(str, "CO2")) {
			fluidProp[i][MW]=44.010;
			fluidProp[i][TCRIT]=304.19;
			fluidProp[i][PCRIT]=7382000;
			fluidProp[i][VCRIT]=94;
			fluidProp[i][AC]=0.2276;
			fluidProp[i][PARACHOR]=82.00;
		}
		else if (!strcmp(str, "O2")) {
			fluidProp[i][MW]=31.999;
			fluidProp[i][TCRIT]=154.58;
			fluidProp[i][PCRIT]=5043000;
			fluidProp[i][VCRIT]=73.4;
			fluidProp[i][AC]=0.0218;
		}
		else if (!strcmp(str, "N2")) {
			fluidProp[i][MW]=28.014;
			fluidProp[i][TCRIT]=126.1;
			fluidProp[i][PCRIT]=3394000;
			fluidProp[i][VCRIT]=90.1;
			fluidProp[i][AC]=0.0403;
			fluidProp[i][PARACHOR]=61.12;
		}
		else if (!strcmp(str, "H2S")) {
			fluidProp[i][MW]=34.082;
			fluidProp[i][TCRIT]=373.53;
			fluidProp[i][PCRIT]=8963000;
			fluidProp[i][VCRIT]=98.5;
			fluidProp[i][AC]=0.0827;
			fluidProp[i][PARACHOR]=85.50;
		}
		else if (!strcmp(str, "SO2")) {
			fluidProp[i][MW]=64.065;
			fluidProp[i][TCRIT]=430.75;
			fluidProp[i][PCRIT]=7884000;
			fluidProp[i][VCRIT]=122;
			fluidProp[i][AC]=0.2451;
		}
		else if (!strcmp(str, "H2")) {
			fluidProp[i][MW]=2.016;
			fluidProp[i][TCRIT]=33.18;
			fluidProp[i][PCRIT]=1313000;
			fluidProp[i][VCRIT]=64.2;
			fluidProp[i][AC]=0.2150;
		}
		else if (!strcmp(str, "H2O")) {
			fluidProp[i][MW]=18.015;
			fluidProp[i][TCRIT]=647.13;
			fluidProp[i][PCRIT]=22055000;
			fluidProp[i][VCRIT]=56;
			fluidProp[i][AC]=0.3449;
		}
		else {
			TerM("Unknown Component!");
		}
	}
	EOS_Init();
	
	File_Search(fp, "SWT");
	Read_Word(fp, str);
	for (i=0; i<Nswt; i++) for (j=0; j<SAT_TABLE; j++) {
		if (!Read_Word(fp, str)) TerM("Incorrect SWT keyword format in the input file!");
		tempL=atof(str);
		if (tempL) swt[i][j]=tempL;
		else swt[i][j]=RELPERM0;
		//tempL=atof(str);
	}


	File_Search(fp, "SGT");
	Read_Word(fp, str);
	for (i=0; i<Nsgt; i++) for (j=0; j<SAT_TABLE; j++) {
		if (!Read_Word(fp, str)) TerM("Incorrect SGT keyword format in the input file!");
		tempL=atof(str);
		if (tempL) sgt[i][j]=tempL;
		else sgt[i][j]=RELPERM0;
		//sgt[i][j]=atof(str);
	}

	for (i=0; i<Nc; i++)
		for (j=0; j<Nc; j++) bic[i][j]=1;
	if (!File_Search(fp, "BIC")) puts("Warning: No BIC keyword in the input file!");
	else {
		for (i=1; i<Nc; i++)
			for (j=0; j<i; j++) {
				if (!Read_Word(fp, str)) TerM("Incorrect BIC keyword format in the input file!");
				tempL=1-atof(str);
				bic[i][j]=tempL;
				bic[j][i]=tempL;
			}
	}

	if (!File_Search(fp, "INITCOND")) TerM("No INITCOND keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect INITCOND keyword format in the input file!");
	initCond=atoi(str);
	switch (initCond) {
		case 0:		//all
			if (!File_Search(fp, "IPRESS")) TerM("No IPRESS keyword in the input file!");		//initial pressure
			eclint=0;
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) {
				if (!eclint) {
					if (!Read_Word(fp, str)) TerM("Incorrect IPRESS keyword format in the input file!");
					ECLStar(str, &eclint, &tempL);
				}
				P[i+1][j+1][k+1][1]=tempL;
				eclint--;
			}

			if (!File_Search(fp, "IWS")) TerM("No IWS keyword in the input file!");		//Initial water saturation
			eclint=0;
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) {
				if (!eclint) {
					if (!Read_Word(fp, str)) TerM("Incorrect IWS keyword format in the input file!");
					ECLStar(str, &eclint, &tempL);
				}
				sat[i][j][k][0]=tempL;
				eclint--;
			}

			if (!File_Search(fp, "IGC")) TerM("No IGC keyword in the input file!");		//Initial global composition
			eclint=0;
			for (n=0; n<Nc; n++) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) {
				if (!eclint) {
					if (!Read_Word(fp, str)) TerM("Incorrect IGC keyword format in the input file!");
					ECLStar(str, &eclint, &tempL);
				}
				comp[i][j][k][n][2]=tempL;
				eclint--;
			}
			AllFlash();
		break;

		case 1:		//same as 0 but depth variation only
			if (!File_Search(fp, "IPRESS")) TerM("No IPRESS keyword in the input file!");		//initial pressure
			for (k=0; k<Nz; k++) {
				if (!Read_Word(fp, str)) TerM("Incorrect IPRESS keyword format in the input file!");
				tempL=atof(str);
				for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) P[i+1][j+1][k+1][1]=tempL;
			}

			if (!File_Search(fp, "IWS")) TerM("No IWS keyword in the input file!");		//Initial water saturation
			for (k=0; k<Nz; k++) {
				if (!Read_Word(fp, str)) TerM("Incorrect IWS keyword format in the input file!");
				tempL=atof(str);				
				for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) sat[i][j][k][0]=atof(str);
			}

			if (!File_Search(fp, "IGC")) TerM("No IGC keyword in the input file!");		//Initial global composition
			for (k=0; k<Nz; k++) 
				for (n=0; n<Nc; n++) {
					if (!Read_Word(fp, str)) TerM("Incorrect IGC keyword format in the input file!");
					tempL=atof(str);
					for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) comp[i][j][k][n][2]=tempL;
			}
			
			AllFlash();
		break;

		case 2:		//compositional grading
			if (!File_Search(fp, "REFLAYER")) TerM("No REFLAYER keyword in the input file!");
			if (!Read_Word(fp, str)) TerM("Incorrect REFLAYER keyword format in the input file!");
			refL=atoi(str);

			if (!File_Search(fp, "REFPRES")) TerM("No REFPRES keyword in the input file!");
			if (!Read_Word(fp, str)) TerM("Incorrect REFPRES keyword format in the input file!");
			tempL=atof(str);
			for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) P[i+1][j+1][refL+1][1]=tempL;

			if (!File_Search(fp, "WOCHEIGHT")) TerM("No WOCHEIGHT keyword in the input file!");
			if (!Read_Word(fp, str)) TerM("Incorrect WOCHEIGHT keyword format in the input file!");
			WOCHeight=atof(str);
			
			if (!File_Search(fp, "REFCOMP")) TerM("No REFCOMP keyword in the input file!");
			for(n=0; n<Nc; n++) {
				if (!Read_Word(fp, str)) TerM("Incorrect REFCOMP keyword format in the input file!");
				tempL=atof(str);
				for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) comp[i][j][refL][n][2]=tempL;
			}
			//GCE_NR();  //compistional grading module ->not exist yet
		break;
	}
	
	
		
	
	/*
	Well Type
	Well X block
	Well Y block
	Well Z starting block		Well Z ending block

	Enter all well rates in positive form, the sign will be determined automatically by the well type
	*/
	if (!File_Search(fp, "WELLS")) puts("Warning: No WELL keyword in the input file!");
	Read_Word(fp, str);
	n=0;
	for(i=0; i<wellNO; i++){
		for(j=0; j<WELL_I; j++){
			if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
			welli[i][j]=atoi(str)-1;
		}
		(welli[i][0])++;
		switch (welli[i][0]) {
			case 0:		//Constant Pwf- Production
				for(j=0; j<5; j++){		//Pwf, Skin, Rw, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>2) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}
				}				
				break;

			case 1:		//Constant Pwf- Water injection
				for(j=0; j<5; j++){		//Pwf, Skin, Rw, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>2) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}
				}				

				break;

			case 2:		//Constant Pwf- Gas injection
				for(j=0; j<(Nc+5); j++){		//Pwf, Skin, Rw, Nc Composition, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>(2+Nc)) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}
				}			

				break;

			case 3:		//Constant Production Oil Rate
				for(j=0; j<3; j++){		//Qo, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>0) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}
				}			

				break;

			case 4:		//Constant Production Total Rate
				for(j=0; j<3; j++){		//QT, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>0) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}		
				}							

				break;

			case 5:		//Constant Water Injection Rate
				for(j=0; j<3; j++){		//Qw, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>0) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}
				}
				break;

			case 6:		//Constant Gas Injection Rate
				for(j=0; j<(Nc+3); j++){		//Qg, Nc Components, T_Start, T_End
					if (!Read_Word(fp, str)) TerM("Incorrect WELL keyword format in the input file!");
					wellf[i][j]=atof(str);

					if (j>Nc) {
						wellf[i][j]*=TIMECHFACT;
						TStepMarker[n]=wellf[i][j];
						n++;
					}
				}				

				break;

			default:
				TerM("Incorrect WELL keyword format in the input file!");
		}		
	}

	
	//Time step and total time
	if (!File_Search(fp, "DT")) TerM("No DT keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect DT keyword format in the input file!");
	Dt=atof(str)*TIMECHFACT;		//Hour In
	if (!File_Search(fp, "TOTALT")) TerM("No TOTALT keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("Incorrect TOTALT keyword format in the input file!");
	totalTime=atof(str)*TIMECHFACT;		//Hour In
	TStepMarker[2*wellNO]=totalTime;	
}

int File_Search(FILE *fp, char *rSeek) {
	register int i;
	char str[MAX_STRING_LENGTH];

	rewind(fp);
	*str='\0';
	do {
		i=Read_Word(fp, str);
		if (!strcmp(str, rSeek)) return -1;
	} while (i);

	return 0;		//Nothing found
}

int Read_Word(FILE *fp, char *rWord) {
	register unsigned char ch;
	register int i=0;

	*rWord='\0';
	do {
		ch=getc(fp);
		if (feof(fp)) {
			return 0;		//Nothing has been read
		}
	} while ((ch<33) || (ch>126));
	
	while ((ch>32) && (ch<127)) {
		*(rWord+i)=ch;
		i++;
		ch=getc(fp);
		if (feof(fp)) {
			*(rWord+i)='\0';
			return 1;		//Read, but end of file also encountered
		}		
	}
	
	*(rWord+i)='\0';
	return -1;		//Correct execution
}



void Allocation(void) {
	register int i, j, k, n;

	bCSRSize=Nx*Ny*Nz*(16*Nc*Nc+31*Nc+17)-2*(2*Nc*Nc+4*Nc+2)*(Nx*Ny+Ny*Nz+Nz*Nx);
	
	if ((gridDim=(FType *) malloc((Nx+Ny+Nz)*sizeof(FType)))==NULL) TerM("Can not allocate memory for gridDim");

	if ((porosity=(FType ***) malloc(Nx*sizeof(FType **)))==NULL) TerM("Can not allocate memory for porosity");
	for (i=0;i<Nx;i++) {
		if ((porosity[i]=(FType **) malloc(Ny*sizeof(FType *)))==NULL) TerM("Can not allocate memory for porosity");
		for (j=0;j<Ny;j++) if ((porosity[i][j]=(FType *) malloc(Nz*sizeof(FType)))==NULL) TerM("Can not allocate memory for porosity");
	}

	if ((Bift=(FType ***) malloc(Nx*sizeof(FType **)))==NULL) TerM("Can not allocate memory for Bift");
	for (i=0;i<Nx;i++) {
		if ((Bift[i]=(FType **) malloc(Ny*sizeof(FType *)))==NULL) TerM("Can not allocate memory for Bift");
		for (j=0;j<Ny;j++) if ((Bift[i][j]=(FType *) malloc(Nz*sizeof(FType)))==NULL) TerM("Can not allocate memory for Bift");
	}

	if ((tor=(FType ***) malloc(Nx*sizeof(FType **)))==NULL) TerM("Can not allocate memory for tortuosity");
	for (i=0;i<Nx;i++) {
		if ((tor[i]=(FType **) malloc(Ny*sizeof(FType *)))==NULL) TerM("Can not allocate memory for tortuosity");
		for (j=0;j<Ny;j++) if ((tor[i][j]=(FType *) malloc(Nz*sizeof(FType)))==NULL) TerM("Can not allocate memory for tortuosity");
	}

	if ((perm=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for perm");
	for (i=0;i<Nx;i++) {
		if ((perm[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for perm");
		for (j=0;j<Ny;j++) {
			if ((perm[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for perm");
			for (k=0;k<Nz;k++) if ((perm[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for perm");
		}
	}

	if ((fluidProp=(FType **) malloc(Nc*sizeof(FType *)))==NULL) TerM("Can not allocate memory for components");
	for (i=0;i<Nc;i++) {
		if ((fluidProp[i]=(FType *) malloc(COMP_PROPS*sizeof(FType)))==NULL) TerM("Can not allocate memory for components");
	}


	if ((swt=(FType **) malloc(Nswt*sizeof(FType *)))==NULL) TerM("Can not allocate memory for liquid saturation table");
	for (i=0;i<Nswt;i++) {
		if ((swt[i]=(FType *) malloc(SAT_TABLE*sizeof(FType)))==NULL) TerM("Can not allocate memory for liquid saturation table");
	}

	if ((sgt=(FType **) malloc(Nsgt*sizeof(FType *)))==NULL) TerM("Can not allocate memory for gas saturation table");
	for (i=0;i<Nsgt;i++) {
		if ((sgt[i]=(FType *) malloc(SAT_TABLE*sizeof(FType)))==NULL) TerM("Can not allocate memory for gas saturation table");
	}

	if ((comp=(FType *****) malloc(Nx*sizeof(FType ****)))==NULL) TerM("Can not allocate memory for composition");
	for (i=0;i<Nx;i++) {
		if ((comp[i]=(FType ****) malloc(Ny*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for composition");
		for (j=0;j<Ny;j++) {
			if ((comp[i][j]=(FType ***) malloc(Nz*sizeof(FType **)))==NULL) TerM("Can not allocate memory for composition");
			for (k=0;k<Nz;k++) {
				if ((comp[i][j][k]=(FType **) malloc(Nc*sizeof(FType *)))==NULL) TerM("Can not allocate memory for composition");
				for (n=0; n<Nc; n++) if ((comp[i][j][k][n]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for composition");
			}
		}
	}

	if ((diffusion=(FType *****) malloc(Nx*sizeof(FType ****)))==NULL) TerM("Can not allocate memory for diffusion");
	for (i=0;i<Nx;i++) {
		if ((diffusion[i]=(FType ****) malloc(Ny*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for diffusion");
		for (j=0;j<Ny;j++) {
			if ((diffusion[i][j]=(FType ***) malloc(Nz*sizeof(FType **)))==NULL) TerM("Can not allocate memory for diffusion");
			for (k=0;k<Nz;k++) {
				if ((diffusion[i][j][k]=(FType **) malloc(Nc*sizeof(FType *)))==NULL) TerM("Can not allocate memory for diffusion");
				for (n=0; n<Nc; n++) if ((diffusion[i][j][k][n]=(FType *) malloc(2*sizeof(FType)))==NULL) TerM("Can not allocate memory for diffusion");
			}
		}
	}

	if ((P=(FType ****) malloc((Nx+2)*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for press");
	for (i=0;i<(Nx+2);i++) {
		if ((P[i]=(FType ***) malloc((Ny+2)*sizeof(FType **)))==NULL) TerM("Can not allocate memory for press");
		for (j=0;j<(Ny+2);j++) {
			if ((P[i][j]=(FType **) malloc((Nz+2)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for press");
			for (k=0;k<(Nz+2);k++) if ((P[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for press");
		}
	}


	if ((sat=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for saturation");
	for (i=0;i<Nx;i++) {
		if ((sat[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for saturation");
		for (j=0;j<Ny;j++) {
			if ((sat[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for saturation");
			for (k=0;k<Nz;k++) if ((sat[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for saturation");
		}
	}

	if ((bsat=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for saturation backup");
	for (i=0;i<Nx;i++) {
		if ((bsat[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for saturation backup");
		for (j=0;j<Ny;j++) {
			if ((bsat[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for saturation backup");
			for (k=0;k<Nz;k++) if ((bsat[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for saturation backup");
		}
	}

	if ((bcomp=(FType *****) malloc(Nx*sizeof(FType ****)))==NULL) TerM("Can not allocate memory for composition");
	for (i=0;i<Nx;i++) {
		if ((bcomp[i]=(FType ****) malloc(Ny*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for composition");
		for (j=0;j<Ny;j++) {
			if ((bcomp[i][j]=(FType ***) malloc(Nz*sizeof(FType **)))==NULL) TerM("Can not allocate memory for composition");
			for (k=0;k<Nz;k++) {
				if ((bcomp[i][j][k]=(FType **) malloc(Nc*sizeof(FType *)))==NULL) TerM("Can not allocate memory for composition");
				for (n=0; n<Nc; n++) if ((bcomp[i][j][k][n]=(FType *) malloc(2*sizeof(FType)))==NULL) TerM("Can not allocate memory for composition");
			}
		}
	}

	if ((bP=(FType ****) malloc((Nx+2)*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for bpress");
	for (i=0;i<(Nx+2);i++) {
		if ((bP[i]=(FType ***) malloc((Ny+2)*sizeof(FType **)))==NULL) TerM("Can not allocate memory for bpress");
		for (j=0;j<(Ny+2);j++) {
			if ((bP[i][j]=(FType **) malloc((Nz+2)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for bpress");
			for (k=0;k<(Nz+2);k++) if ((bP[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for bpress");
		}
	}
	if (wellNO) {
		if ((wellf=(FType **) malloc(wellNO*sizeof(FType *)))==NULL) TerM("Can not allocate memory for wells data");
		for (i=0;i<wellNO;i++) {
			if ((wellf[i]=(FType *) malloc((Nc+5)*sizeof(FType)))==NULL) TerM("Can not allocate memory for wells data");
		}
		if ((welli=(int **) malloc(wellNO*sizeof(int *)))==NULL) TerM("Can not allocate memory for wells data");
		for (i=0;i<wellNO;i++) {
			if ((welli[i]=(int *) malloc(WELL_I*sizeof(int)))==NULL) TerM("Can not allocate memory for wells data");
		}
			
	}


	//////////////////////////////////////////////////////
	/////////////////////NON_INPUT////////////////////////
	//////////////////////////////////////////////////////
	if ((IFTran=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for IFTran");
	for (i=0;i<Nx;i++) {
		if ((IFTran[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for IFTran");
		for (j=0;j<Ny;j++) {
			if ((IFTran[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for IFTran");
			for (k=0;k<Nz;k++) if ((IFTran[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for IFTran");
		}
	}

	if ((blockFProps=(FType *****) malloc(Nx*sizeof(FType ****)))==NULL) TerM("Can not allocate memory for blockFProps");
	for (i=0;i<Nx;i++) {
		if ((blockFProps[i]=(FType ****) malloc(Ny*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for blockFProps");
		for (j=0;j<Ny;j++) {
			if ((blockFProps[i][j]=(FType ***) malloc(Nz*sizeof(FType **)))==NULL) TerM("Can not allocate memory for blockFProps");
			for (k=0;k<Nz;k++) {
				if ((blockFProps[i][j][k]=(FType **) malloc(BLOCK_F_PROPS*sizeof(FType *)))==NULL) TerM("Can not allocate memory for blockFProps");
				for (n=0; n<BLOCK_F_PROPS; n++) if ((blockFProps[i][j][k][n]=(FType *) malloc(2*sizeof(FType)))==NULL) TerM("Can not allocate memory for blockFProps");
			}
		}
	}

	if ((trans=(FType *****) malloc((Nx+1)*sizeof(FType ****)))==NULL) TerM("Can not allocate memory for trans");
	for (i=0;i<(Nx+1);i++) {
		if ((trans[i]=(FType ****) malloc((Ny+1)*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for trans");
		for (j=0;j<(Ny+1);j++) {
			if ((trans[i][j]=(FType ***) malloc((Nz+1)*sizeof(FType **)))==NULL) TerM("Can not allocate memory for trans");
			for (k=0;k<(Nz+1);k++) {
				if ((trans[i][j][k]=(FType **) malloc(12*sizeof(FType*)))==NULL) TerM("Can not allocate memory for trans");
				for (n=0; n<12; n++) if ((trans[i][j][k][n]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for trans");
			}
		}
	}

	if ((transS=(char *****) malloc((Nx+1)*sizeof(char ****)))==NULL) TerM("Can not allocate memory for transS");
	for (i=0;i<(Nx+1);i++) {
		if ((transS[i]=(char ****) malloc((Ny+1)*sizeof(char ***)))==NULL) TerM("Can not allocate memory for transS");
		for (j=0;j<(Ny+1);j++) {
			if ((transS[i][j]=(char ***) malloc((Nz+1)*sizeof(char **)))==NULL) TerM("Can not allocate memory for transS");
			for (k=0;k<(Nz+1);k++) {
				if ((transS[i][j][k]=(char **) malloc(3*sizeof(char*)))==NULL) TerM("Can not allocate memory for transS");
				for (n=0; n<3; n++) if ((transS[i][j][k][n]=(char *) malloc(3*sizeof(char)))==NULL) TerM("Can not allocate memory for transS");
			}
		}
	}

	if ((Wtran=(FType ****) malloc((Nx+1)*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for Wtran");
	for (i=0;i<(Nx+1);i++) {
		if ((Wtran[i]=(FType ***) malloc((Ny+1)*sizeof(FType **)))==NULL) TerM("Can not allocate memory for Wtran");
		for (j=0;j<(Ny+1);j++) {
			if ((Wtran[i][j]=(FType **) malloc((Nz+1)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for Wtran");
			for (k=0;k<(Nz+1);k++) if ((Wtran[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for Wtran");
		}
	}

	if ((dWtran=(FType ****) malloc((Nx+1)*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for dWtran");
	for (i=0;i<(Nx+1);i++) {
		if ((dWtran[i]=(FType ***) malloc((Ny+1)*sizeof(FType **)))==NULL) TerM("Can not allocate memory for dWtran");
		for (j=0;j<(Ny+1);j++) {
			if ((dWtran[i][j]=(FType **) malloc((Nz+1)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for dWtran");
			for (k=0;k<(Nz+1);k++) if ((dWtran[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for dWtran");
		}
	}


	if ((relPerm=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for relPerm");
	for (i=0;i<Nx;i++) {
		if ((relPerm[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for relPerm");
		for (j=0;j<Ny;j++) {
			if ((relPerm[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for relPerm");
			for (k=0;k<Nz;k++) if ((relPerm[i][j][k]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for relPerm");
		}
	}

	if ((dRelPerm=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for dRelPerm");
	for (i=0;i<Nx;i++) {
		if ((dRelPerm[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for dRelPerm");
		for (j=0;j<Ny;j++) {
			if ((dRelPerm[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for dRelPerm");
			for (k=0;k<Nz;k++) if ((dRelPerm[i][j][k]=(FType *) malloc(5*sizeof(FType)))==NULL) TerM("Can not allocate memory for dRelPerm");
		}
	}

	/*if ((preCon=(FType **) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for Jacobian matrix");
	for (i=0;i<(Nx*Ny*Nz*(2*Nc+4));i++) {
		if ((preCon[i]=(FType *) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Jacobian matrix");
	}*/

	if ((CSRjac=(FType *) malloc(bCSRSize*sizeof(FType)))==NULL) TerM("Can not allocate memory for Sparse Jacobian matrix");
	if ((CSRrow=(int *) malloc((Nx*Ny*Nz*(2*Nc+4)+1)*sizeof(int)))==NULL) TerM("Can not allocate memory for Sparse Jacobian Row Index matrix");
	if ((CSRcol=(int *) malloc(bCSRSize*sizeof(int)))==NULL) TerM("Can not allocate memory for Sparse Jacobian Column Index matrix");

	if ((preCon=(FType *) malloc(Nx*Ny*Nz*((2*Nc+4)*(2*Nc+3)+MRMAXNONZERO)*sizeof(FType)))==NULL) TerM("Can not allocate memory for preconditioner matrix");
	if ((preConCSRrow=(int *) malloc((Nx*Ny*Nz*(2*Nc+4)+1)*sizeof(int)))==NULL) TerM("Can not allocate memory for preConCSRrow matrix");
	if ((preConCSRcol=(int *) malloc(Nx*Ny*Nz*((2*Nc+4)*(2*Nc+3)+MRMAXNONZERO)*sizeof(int)))==NULL) TerM("Can not allocate memory for preConCSRrow matrix");


	if ((preConRow=(int *) malloc((Nx*Ny*Nz+1)*sizeof(int)))==NULL) TerM("Can not allocate memory for preConRow matrix");
	if ((preConIndex=(int *) malloc((Nx*Ny*Nz+1)*sizeof(int)))==NULL) TerM("Can not allocate memory for preConIndex matrix");

	
	if ((ans=(FType *) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType)))==NULL) TerM("Can not allocate memory for answer matrix");
	//if ((Unk=(FType *) malloc(Nx*Ny*Nz*(2*Nc+4)*sizeof(FType)))==NULL) TerM("Can not allocate memory for unknown matrix");

	if ((pJHolder=(int ***) malloc(Nx*sizeof(int **)))==NULL) TerM("Can not allocate memory for pJHolder");
	for (i=0;i<Nx;i++) {
		if ((pJHolder[i]=(int **) malloc(Ny*sizeof(int *)))==NULL) TerM("Can not allocate memory for pJHolder");
		for (j=0;j<Ny;j++) if ((pJHolder[i][j]=(int *) malloc(Nz*sizeof(FType)))==NULL) TerM("Can not allocate memory for pJHolder");
	}


	if ((preProp=(FType ****) malloc(Nx*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for preProp");
	for (i=0;i<Nx;i++) {
		if ((preProp[i]=(FType ***) malloc(Ny*sizeof(FType **)))==NULL) TerM("Can not allocate memory for preProp");
		for (j=0;j<Ny;j++) {
			if ((preProp[i][j]=(FType **) malloc(Nz*sizeof(FType *)))==NULL) TerM("Can not allocate memory for preProp");
			for (k=0;k<Nz;k++) if ((preProp[i][j][k]=(FType *) malloc((Nc+1)*sizeof(FType)))==NULL) TerM("Can not allocate memory for preProp");
		}
	}

	if ((dE=(FType *****) malloc(Nx*sizeof(FType ****)))==NULL) TerM("Can not allocate memory for dE");
	for (i=0;i<Nx;i++) {
		if ((dE[i]=(FType ****) malloc(Ny*sizeof(FType ***)))==NULL) TerM("Can not allocate memory for dE");
		for (j=0;j<Ny;j++) {
			if ((dE[i][j]=(FType ***) malloc(Nz*sizeof(FType **)))==NULL) TerM("Can not allocate memory for dE");
			for (k=0;k<Nz;k++) {
				if ((dE[i][j][k]=(FType **) malloc((Nc+1)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for dE");
				for (n=0; n<(Nc+1); n++) if ((dE[i][j][k][n]=(FType *) malloc(2*sizeof(FType)))==NULL) TerM("Can not allocate memory for dE");
			}
		}
	}

	if ((satJac=(FType **) malloc((Nc+1)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for Flash Jacobian matrix");
	for (i=0;i<(Nc+1);i++) {
		if ((satJac[i]=(FType *) malloc((Nc+1)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Flash Jacobian matrix");
	}

	if ((satAns=(FType *) malloc((Nc+1)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Flash Answer matrix");
	if ((Xm=(FType *) malloc((Nc+1)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Flash Composition matrix");
	if ((Xms=(FType *) malloc((Nc+1)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Flash Composition save matrix");

	if ((phaseStat=(char ***) malloc(Nx*sizeof(char **)))==NULL) TerM("Can not allocate memory for phaseStat");
	for (i=0;i<Nx;i++) {
		if ((phaseStat[i]=(char **) malloc(Ny*sizeof(char *)))==NULL) TerM("Can not allocate memory for phaseStat");
		for (j=0;j<Ny;j++) if ((phaseStat[i][j]=(char *) malloc(Nz*sizeof(char)))==NULL) TerM("Can not allocate memory for phaseStat");
	}

	if ((bphaseStat=(char ***) malloc(Nx*sizeof(char **)))==NULL) TerM("Can not allocate memory for bphaseStat");
	for (i=0;i<Nx;i++) {
		if ((bphaseStat[i]=(char **) malloc(Ny*sizeof(char *)))==NULL) TerM("Can not allocate memory for bphaseStat");
		for (j=0;j<Ny;j++) if ((bphaseStat[i][j]=(char *) malloc(Nz*sizeof(char)))==NULL) TerM("Can not allocate memory for bphaseStat");
	}


	if ((Ki=(FType *) malloc(Nc*sizeof(FType)))==NULL) TerM("Can not allocate memory for Initial Flash K matrix");
	
	if ((blockH=(FType ***) malloc((Nx+2)*sizeof(FType **)))==NULL) TerM("Can not allocate memory for blockH");
	for (i=0;i<(Nx+2);i++) {
		if ((blockH[i]=(FType **) malloc((Ny+2)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for blockH");
		for (j=0;j<(Ny+2);j++) if ((blockH[i][j]=(FType *) malloc((Nz+2)*sizeof(FType)))==NULL) TerM("Can not allocate memory for blockH");
	}

	if ((bic=(FType **) malloc(Nc*sizeof(FType *)))==NULL) TerM("Can not allocate memory for binary interaction coefficient matrix");
	for (i=0;i<Nc;i++) {
		if ((bic[i]=(FType *) malloc(Nc*sizeof(FType)))==NULL) TerM("Can not allocate memory for binary interaction coefficient matrix");
	}

	if ((TStepMarker=(FType *) malloc((2*wellNO+1)*sizeof(FType)))==NULL) TerM("Can not allocate memory for Time Step Marker");

	if ((STcomp=(FType **) malloc(Nc*sizeof(FType *)))==NULL) TerM("Can not allocate memory for STcomp");
	for (i=0;i<Nc;i++) {
		if ((STcomp[i]=(FType *) malloc(3*sizeof(FType)))==NULL) TerM("Can not allocate memory for STcomp");
	}

	/*k=(2*Nc+4)*Nx*Ny*Nz;
	if ((fullPre=(FType **) malloc(k*sizeof(FType *)))==NULL) TerM("Can not allocate memory for fullPre matrix");
	for (i=0;i<k;i++) {
		if ((fullPre[i]=(FType *) malloc(k*sizeof(FType)))==NULL) TerM("Can not allocate memory for fullPre matrix");
	}*/
}

void ECLStar(char *str, int *times, FType *mvalue) {
	register int i, j;
	char Sind=0;
	char str1[MAX_STRING_LENGTH], str2[MAX_STRING_LENGTH];
	
	i=0;
	j=0;
	while (str[i]!='\0') {
		if (str[i]=='*') {
			str1[i]='\0';
			Sind=-1;
		}
		else {
			if (Sind) {
				str2[j]=str[i];
				j++;
			}
			else {
				str1[i]=str[i];
			}
		}
		i++;
	}

	if (Sind) {
		str2[j]='\0';
		*times=atoi(str1);
		*mvalue=atof(str2);
	}
	else {
		str1[i]='\0';
		*times=1;
		*mvalue=atof(str1);
	}
}

void RestoreRST(FType *simTime) {
	register int Ix, Iy, Iz, i;
	FILE *fp;

	if ((fp=fopen("test.rst","rb"))!=NULL) {

		for (Ix=0; Ix<(Nx+2); Ix++)
			for (Iy=0; Iy<(Ny+2); Iy++)
				for (Iz=0; Iz<(Nz+2); Iz++) 
					for (i=0; i<3; i++) {
						fread(&P[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
						//if (ferror(fp)) puts("Error in File in!");
				}

				for (Ix=0; Ix<Nx; Ix++)
					for (Iy=0; Iy<Ny; Iy++)
						for (Iz=0; Iz<Nz; Iz++) {
							for (i=0; i<3; i++) {
								fread(&sat[Ix][Iy][Iz][i], sizeof(FType), 1, fp);
								//if (ferror(fp)) puts("Error in File in!");
							}

							for (i=0; i<Nc; i++) {
								fread(&comp[Ix][Iy][Iz][i][0], sizeof(FType), 1, fp);
								//if (ferror(fp)) puts("Error in File in!");
								fread(&comp[Ix][Iy][Iz][i][1], sizeof(FType), 1, fp);
								//if (ferror(fp)) puts("Error in File in!");
							}
							fread(&phaseStat[Ix][Iy][Iz], sizeof(char), 1, fp);
							//if (ferror(fp)) puts("Error in File in!");
						}
		fread(simTime, sizeof(FType), 1, fp);
		//if (ferror(fp)) puts("Error in File in!");
		fread(&Dt, sizeof(FType), 1, fp);
		//if (ferror(fp)) puts("Error in File in!");
		fread(&TSMCtrl, sizeof(int), 1, fp);
		//if (ferror(fp)) puts("Error in File in!");

		fread(&SumQoProduced, sizeof(FType), 1, fp);
		fread(&SumQgInjected, sizeof(FType), 1, fp);

		fclose(fp);

		
	}
}
#endif