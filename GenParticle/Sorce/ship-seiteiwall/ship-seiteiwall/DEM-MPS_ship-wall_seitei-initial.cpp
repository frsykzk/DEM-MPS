#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include<iostream>
#include<iomanip>
using namespace std;


int main(void) {
	int Nall = 0;//合計粒子数
	int NumSLD = 0;//固体粒子数
	int NumFLD = 0;//流体粒子数
	///////////////////////////////////////////固体読み込み////////////////////////////////////////////////////////
	{
		int p;
		char filename2[10000];
		char filenamea[10000];


		{
			sprintf_s(filenamea, "DEM_ship.csv");
			ifstream file(filenamea);

			vector<vector<string>> values;
			string str;

			if (file.fail()) {
				cerr << "failed." << endl;
				exit(0);
			}

			while (getline(file, str)) {
				//コメント箇所は除く
				if ((p = str.find("//")) != str.npos) continue;
				vector<string> inner;

				//コンマがあるかを探し、そこまでをvaluesに格納
				while ((p = str.find(",")) != str.npos) {
					inner.push_back(str.substr(0, p));

					//strの中身は", "の2文字を飛ばす
					str = str.substr(p + 1);
				}
				inner.push_back(str);
				values.push_back(inner);
			}

			sprintf_s(filename2, "DEM_test.txt");

			ifstream file2(filename2);

			std::ofstream ofs(filename2);

			for (unsigned int i = 1; i < values.size(); ++i) {
				NumSLD++;
				Nall += 1;
				for (unsigned int j = 7; j < values[i].size(); ++j) { //consid 有り：j=6　無し：j=5
					ofs << values[i][j] << "	";

				}
				ofs << std::endl;
			}
			fprintf(stderr, "solid particle =\t%d now \n", NumSLD); fflush(stderr);
			fprintf(stderr, "ALL1=%d \n", Nall); fflush(stderr);
		}
	}
	///////////////////////////////////////////固体読み込み////////////////////////////////////////////////////////


	///////////////////////////////////////////流体読み込み////////////////////////////////////////////////////////
	{
		int p;
		char filename2[10000];
		char filenamea[10000];


		{
			sprintf_s(filenamea, "MPS_ship.csv");
			ifstream file(filenamea);

			vector<vector<string>> values;
			string str;

			if (file.fail()) {
				cerr << "failed." << endl;
				exit(0);
			}

			while (getline(file, str)) {
				//コメント箇所は除く
				if ((p = str.find("//")) != str.npos) continue;
				vector<string> inner;

				//コンマがあるかを探し、そこまでをvaluesに格納
				while ((p = str.find(",")) != str.npos) {
					inner.push_back(str.substr(0, p));

					//strの中身は", "の2文字を飛ばす
					str = str.substr(p + 1);
				}
				inner.push_back(str);
				values.push_back(inner);
			}

			sprintf_s(filename2, "MPS_test.txt");

			ifstream file2(filename2);

			std::ofstream ofs(filename2);

			for (unsigned int i = 1; i < values.size(); ++i) {
				NumFLD++;
				Nall += 1;
				for (unsigned int j = 7; j < values[i].size(); ++j) { //consid 有り：j=6　無し：j=5
					ofs << values[i][j] << "	";

				}
				ofs << std::endl;
			}
			fprintf(stderr, "water particle =\t%d now \n", NumFLD); fflush(stderr);
			fprintf(stderr, "ALL2=%d \n", Nall); fflush(stderr);
		}
	}
	///////////////////////////////////////////流体読み込み////////////////////////////////////////////////////////


	///////////////////////////////////////////壁粒子数読み込み////////////////////////////////////////////////////////
	{
		char filename4[10000];
		sprintf_s(filename4, "initialwall.txt");
		ifstream IN4(filename4);
		int ppp;//壁粒子数
		IN4 >> ppp;
		Nall += ppp;
		fprintf(stderr, "Now +wall =\t%d \n", Nall); fflush(stderr);
	}

	fprintf(stderr, "all particle =\t%d \n", Nall); fflush(stderr);
	FILE* output2;
	fopen_s(&output2, "initial.txt", "a");
	fprintf(output2, "%d\n", Nall);
	fclose(output2);
	///////////////////////////////////////////壁粒子数読み込み////////////////////////////////////////////////////////

	
	///////////////////////////////////////////固体出力////////////////////////////////////////////////////////
	{
		char filename3[10000];

		float* rx;
		float* ry;
		float* rz;
		rx = new float[NumSLD];
		ry = new float[NumSLD];
		rz = new float[NumSLD];

		int ii;

		sprintf_s(filename3, "DEM_test.txt");
		ifstream IN(filename3);

		for (ii = 0; ii < NumSLD; ii++)
		{
			IN >> rx[ii] >> ry[ii] >> rz[ii];
		}

		fprintf(stderr, "myrankDEM %d =\t%d\n", 0, NumSLD); fflush(stderr);

		FILE* output2;
		fopen_s(&output2, "initial.txt", "a");
		for (ii = 0; ii < NumSLD; ii++)
		{
			fprintf(output2, "%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", ii, 0, rx[ii], ry[ii], rz[ii], 0.0, 0.0, 0.0, 0.0);//固体Typ=0
		}
		fclose(output2);

		delete(rx); delete(ry); delete(rz);

	}
	///////////////////////////////////////////固体出力////////////////////////////////////////////////////////


	///////////////////////////////////////////流体出力////////////////////////////////////////////////////////
	{
		char filename3[10000];

		float* rx;
		float* ry;
		float* rz;
		rx = new float[NumFLD];
		ry = new float[NumFLD];
		rz = new float[NumFLD];

		int ii;

		sprintf_s(filename3, "MPS_test.txt");
		ifstream IN(filename3);

		for (ii = 0; ii < NumFLD; ii++)
		{
			IN >> rx[ii] >> ry[ii] >> rz[ii];
		}

		fprintf(stderr, "myrankMPS %d =\t%d\n", 0, NumFLD); fflush(stderr);

		FILE* output2;
		fopen_s(&output2, "initial.txt", "a");
		for (ii = 0; ii < NumFLD; ii++)
		{
			fprintf(output2, "%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", ii, 1, rx[ii], ry[ii], rz[ii], 0.0, 0.0, 0.0, 0.0);//流体Typ=1
		}
		fclose(output2);

		delete(rx); delete(ry); delete(rz);

	}
	///////////////////////////////////////////流体出力////////////////////////////////////////////////////////
	

	///////////////////////////////////////////壁出力////////////////////////////////////////////////////////
	{

		char filename4[10000];
		sprintf_s(filename4, "initialwall.txt");
		ifstream IN4(filename4);
		int ppp;
		int null;
		IN4 >> ppp;
		float* rx = new float[int(ppp)];
		float* ry = new float[int(ppp)];
		float* rz = new float[int(ppp)];
		float* n = new float[int(ppp)];
		float* nx = new float[int(ppp)];
		float* ny = new float[int(ppp)];
		float* nz = new float[int(ppp)];
		fprintf(stderr, "wall =\t%d \n", ppp); fflush(stderr);
		for (int ii = 0; ii < ppp; ii++)
		{
			IN4 >> null >> null >> rx[ii] >> ry[ii] >> rz[ii] >> nx[ii] >> ny[ii] >> nz[ii] >> n[ii];
		}

		{
			FILE* output2;
			fopen_s(&output2, "initial.txt", "a");
			for (int ii = 0; ii < ppp; ii++)
			{
				fprintf(output2, "%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", ii, 2, rx[ii], ry[ii], rz[ii], nx[ii], ny[ii], nz[ii], n[ii]); //壁Typ = 2
			}
			fclose(output2);
		}

		delete(rx); delete(ry); delete(rz);
		delete(nx); delete(ny); delete(nz); delete(n);

	}
	///////////////////////////////////////////壁出力////////////////////////////////////////////////////////

	return 0;
}
