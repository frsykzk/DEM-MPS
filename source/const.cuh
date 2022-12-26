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

//////////////////////////////�o�̓t���O///////////////////////////////////////

#define Multi_flg 0 //�ŉt�������Ȃ�P, �P����͂Ȃ�0
#define DEM_flg 1 //MPS�P�̂Ōv�Z�����0 
#define MPS_flg 1 //DEM�P�̂Ōv�Z����Ƃ�0
#define move_WLL 0 //�����ǂ̎��P�@�ǂ�obj���o�͂���t���O

//////////////////////////////�o�̓t���O///////////////////////////////////////

#define IN_FILE_0 "./input/physical.txt"
#define IN_FILE_1 "./input/DEM_physical.txt"
#define IN_FILE_2 "./input/MPS_physical.txt"
#define IN_FILE_3 "./input/initial.txt"

#define Dns_Num 6
#define GST -1//�S�[�X�g���q
#define SLD 0//���̗��q
#define FLD 1//���̗��q
#define WLL 2//�Ǘ��q
#define OBJ 3//���Ǘ��q
#define OBJ2 4//���Ǘ��q
#define MRR 5//�~���[���q

#define Dns_SLD 2500.0f	//�e���q�̖��x
#define Dns_FLD 1000.0f
#define Dns_WLL 1000.0f
#define Dns_OBJ 1000.0f
#define Dns_OBJ2 1000.0f

#define pi 3.141592f

#define WEI(dst, re) ((re/dst) + (dst/re) - 2.0f)
#define WEI_grad(dst, re) ((re/dst) - (dst/re))

#define NCP 20//�ő̗��q�ő�ڐG��(�\���ȔC�ӂ̑傫��)		�Ŗ��[�U��12�ɂȂ�͂�
#define NumMRR 10//1�̗��̗��q���琶�������~���[���q�̍ő�l(�\���ȔC�ӂ̑傫��)

#define Surface 0//�ǖ�
#define Edge 1//�Ǌp


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


/*����
�Efor(;;)�̒��ł�Typ�����continue�����g��Ȃ�����
�EDEM_obj���Q�l�ɁA�͂��Z�o�E�o�͂ł���悤�ɂ�����
�Esetpara�ŏd��G.y���g���Ƃ��͐�Βl�ɂ��邱��(�ǂݍ��񂾃f�[�^�͕��̐��ɂȂ��Ă邩��)

*/