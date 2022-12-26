#pragma once
#pragma warning

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
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
//#include <device_atomic_functions.hpp>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include <chrono>
#include <random>


#define THREADS 512

//////////////////////////////出力フラグ///////////////////////////////////////

#define Multi_flg 0 //固液混相流なら１, 単相解析なら0
#define DEM_flg 1 //MPS単体で計算すると0 
#define MPS_flg 1 //DEM単体で計算するとき0
#define move_WLL 0 //動く壁の時１　壁とobjを出力するフラグ

//////////////////////////////出力フラグ///////////////////////////////////////

#define IN_FILE_0 "./input/physical.txt"
#define IN_FILE_1 "./input/DEM_physical.txt"
#define IN_FILE_2 "./input/MPS_physical.txt"
#define IN_FILE_3 "./input/initial.txt"

#define Dns_Num 6
#define GST -1//ゴースト粒子
#define SLD 0//流体粒子
#define FLD 1//流体粒子
#define WLL 2//壁粒子
#define OBJ 3//動壁粒子
#define OBJ2 4//動壁粒子
#define MRR 5//ミラー粒子

#define Dns_SLD 2500.0f	//各粒子の密度
#define Dns_FLD 1000.0f
#define Dns_WLL 1000.0f
#define Dns_OBJ 1000.0f
#define Dns_OBJ2 1000.0f

#define pi 3.141592f

#define WEI(dst, re) ((re/dst) + (dst/re) - 2.0f)
#define WEI_grad(dst, re) ((re/dst) - (dst/re))

#define NCP 20//固体粒子最大接触数(十分な任意の大きさ)		最密充填で12になるはず
#define NumMRR 10//1つの流体粒子から生成されるミラー粒子の最大値(十分な任意の大きさ)

#define Surface 0//壁面
#define Edge 1//壁角


typedef float real;
//typedef double real;

typedef struct {
	real x, y, z;
}treal3;

typedef struct {
	char x, y, z;
}tchar3;

typedef struct {
	real max, min;
}treal_m2;

typedef struct {
	real* x;
	real* y;
	real* z;
}areal3;

typedef struct {
	int* x;
	int* y;
	int* z;
}aint3;

typedef struct {
	real* x;
	real* y;
	real* z;
}host_vec3;


//CUDA error check
#define CHECK(call)																							\
{																														\
	const cudaError_t error = call;																		\
	if (error != cudaSuccess) {																				\
		printf_s("Error: %s:%d, ", __FILE__, __LINE__);										\
		printf_s("code:%d, reason: %s\n", error, cudaGetErrorString(error));			\
		exit(1);																										\
	}																													\
}


/*メモ
・for(;;)の中ではTyp判定にcontinue処理使わないこと
・DEM_objを参考に、力も算出・出力できるようにしたい
・setparaで重力G.yを使うときは絶対値にすること(読み込んだデータは負の数になってるから)

*/