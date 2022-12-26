#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h> 
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <chrono>

#define IN_FILE_1 "DEM_seitei.txt"
#define IN_FILE_2 "MPS_seitei.txt"
#define IN_FILE_3 "initialwall.txt"
#define OUT_FILE "initial.txt"


int main() {
	FILE* in1;
	FILE* in2;
	FILE* in3;
	FILE* out;

	int nPwall = 0;
	int nPsolid = 0;
	int nPfluid = 0;
	int nP = 0;

	if (fopen_s(&out, OUT_FILE, "a") != 0) {
		printf("initial.txtが開けません\n");
	}
	else {
		printf_s("initial.txt opend!\n");
	}

	if (fopen_s(&in1, IN_FILE_1, "r") != 0) {
		printf("DEM_seitei.txtが開けません\n");
	}
	else {
		fscanf_s(in1, "%d", &nPsolid);
	}

	if (fopen_s(&in2, IN_FILE_2, "r") != 0) {
		printf("MPS_seitei.txtが開けません\n");
	}
	else {
		fscanf_s(in2, "%d", &nPfluid);
	}

	if (fopen_s(&in3, IN_FILE_3, "r") != 0) {
		printf("initialwall.txtが開けません\n");
	}
	else {
		fscanf_s(in3, "%d", &nPwall);
	}

	nP = nPsolid + nPfluid + nPwall;
	fprintf_s(out, "%d\n", nP);
	printf("固体粒子：%d個\n流体粒子：%d個\n壁粒子：%d個\n総粒子%：%d個\n", nPsolid, nPfluid, nPwall, nP);


	for (int j = 0; j < nPsolid; j++) {
		int a[1];
		float b[11];
		int c[1];
		float g[1];
		fscanf_s(in1, " %d %d %f %f %f %f %f %f %f", &a[0], &c[0], &b[0], &b[1], &b[2], &b[8], &b[9], &b[10], &g[0]);
		fprintf_s(out, " %d %d %f %f %f %f %f %f %f\n", a[0], c[0], b[0], b[1], b[2], b[8], b[9], b[10], g[0]);
	}
	fclose(in1);

	for (int j = 0; j < nPfluid; j++) {
		int a[1];
		float b[11];
		int c[1];
		float g[1];
		fscanf_s(in2, " %d %d %f %f %f %f %f %f %f", &a[0], &c[0], &b[0], &b[1], &b[2], &b[8], &b[9], &b[10], &g[0]);
		fprintf_s(out, " %d %d %f %f %f %f %f %f %f\n", a[0], c[0], b[0], b[1], b[2], b[8], b[9], b[10], g[0]);
	}
	fclose(in2);

	for (int j = 0; j < nPwall; j++) {
		int a[1];
		float b[11];
		int c[1];
		float g[1];
		fscanf_s(in3, " %d %d %f %f %f %f %f %f %f", &a[0], &c[0], &b[0], &b[1], &b[2], &b[8], &b[9], &b[10], &g[0]);
		fprintf_s(out, " %d %d %f %f %f %f %f %f %f\n", a[0], c[0], b[0], b[1], b[2], b[8], b[9], b[10], g[0]);
	}
	fclose(in3);

	fclose(out);
	printf_s("Finished!!!\n");
	
	return 0;

}