#pragma once

#include "const.cuh"

class DEMPS {
public:
	///////////���ʊ֐�////////////////
	void RdDat();
	void WrtDat();//���q�f�[�^�o��
	void WrtDatWLL();//�ǃf�[�^�o��
	void WrtDat2();//seitei�o��
	void Output(char, int);//�o��
	void Output2(char, int);//seitei�o��
	void AlcBkt();
	void SetPara();
	void MkBkt();
	void ClcDEMPS();
	void checkPCL();//device�֐��ł����������H
	void DevicetoHost();
	void memory_free();
	///////////���ʊ֐�////////////////


	///////////DEM�֐�////////////////
	void ColForce();//�ڐG�͌v�Z
	void update();//�ő����q�X�V
	///////////DEM�֐�////////////////


	///////////MPS�֐�////////////////
	void VscTrm();
	void UpPcl1();
	void ChkCol();
	void MkPrs();
	void PrsGrdTrm();
	void UpPcl2();
	void Surface_Edge();
	void ResetMRR();
	void GenMRR_nonslip();
	void MkBkt_MRR();
	///////////MPS�֐�////////////////


	///////////���ݍ�p�֐�////////////////
	void SLD_FLD();//���ݍ�p�͌v�Z
	void SF_update();//���q�X�V
	///////////���ݍ�p�֐�////////////////


	///////////physical.txt////////////////
	treal3 MINc;//�v�Z�͈�min
	treal3 MAXc;//�v�Z�͈�max
	treal3 G;//�d��
	real FIN_TIM;//�v�Z�I���^�C��
	real output_time;//�t�@�C���o�͎��ԊԊu
	real WLL_PCL_DST;//�Ǘ��q���a
	///////////DEM_physical.txt////////////////


	///////////DEM_physical.txt////////////////
	real DEMdt;//���ԍ��݊Ԋu		���肷�鐔�l�ɐݒ�@MPS�̃^�C���X�e�b�v�Ɣ�ׂď����������̗p
	real k;//�΂˒萔
	real mu;//���C�W��
	real SLD_PCL_DST;//�ő̗��q���a
	///////////DEM_physical.txt////////////////


	///////////DEM_physical.txt////////////////
	real Ma;//�}�b�n��
	real CRT_NUM;//�N�[������
	real COL_RAT;//���̏Փ˂̔����W��
	real DST_LMT_RAT;//���E�ڋߋ���
	real KNM_VSC;//���̂̓��S���W��
	real FLD_PCL_DST;//���̗��q���a
	///////////DEM_physical.txt////////////////





	///////////���ʕϐ�////////////////
	areal3 Pos;//���W
	areal3 d_Pos;
	areal3 Vel;//���x
	areal3 d_Vel;
	char* Typ;//���q�^�C�v�@
	char* d_Typ;

	real dt;//�^�C���X�e�b�v
	int iF = 0;//�t�@�C���ԍ�
	real TIM = 0.0f;//���݂̌v�Z����
	real outtime = 0.0f;//�o�͊Ԋu
	char outout_filename[256];
	FILE* fp;
	real OPT_FQC = 1.0f;//�t�@�C���o�͉�

	real DB, DB2, DBinv;//���q�T���̕����@�̃u���b�N�̕�
	int nBx, nBy, nBz, nBxy, nBxyz;//�o�b�P�g�p
	int* d_bfst, * d_blst, * d_nxt;//�����N���X�g�p

	int nP;//�������q��
	int nPWLL;
	int nPSLD;
	int nPFLD;
	int nPOBJ;
	int nPOBJ2;

	real PCL_DST;//�o�P�b�g��鎞�Ɏg���@��ԑ傫�����q�a

	real* D;//���q���a�z��
	real* d_D;

	real Vol_SLD;//���q�̐ρ@
	real Vol_FLD;
	///////////���ʕϐ�////////////////


	///////////DEM�ϐ�////////////////
	areal3 Ftotal;//�ڐG��
	areal3 d_Ftotal;
	areal3 Omega;//�p���x
	areal3 d_Omega;
	areal3 Torque;//�g���N
	areal3 d_Torque;

	areal3 ep;
	areal3 d_ep;
	int* pair;
	int* d_pair;

	real m;//����1���q����

	real kn;//�΂˒萔
	real kt;

	real eta_n;//�S�������W��
	real eta_t;

	real I;//���̊������[�����g

	real ulmax;
	real usmax;
	///////////DEM�ϐ�////////////////


	///////////MPS�ϐ�////////////////
	areal3 Acc;//�����x
	areal3 d_Acc;
	real* Prs;//����
	real* d_Prs;
	real* pav;//���͕���
	real* d_pav;
	real* n;//���q�����x
	real* d_n;
	real* Dns;//���q���x
	real* d_Dns;

	areal3 PosM;//���W
	areal3 d_PosM;
	areal3 VelM;//���x
	areal3 d_VelM;
	real* PrsM;//����
	real* d_PrsM;
	char* TypM;//�~���[�^�C�v�@�~���[�@�S�[�X�g
	char* d_TypM;

	int* d_FromWLL;//�ǂ̕Ǘ��q(�ԍ�)���琶�����ꂽ�~���[����ۑ�

	char* WLLSE;
	char* d_WLLSE;//�ǃ^�C�v�@surface�@edge
	areal3 WLLVec;//�ǂ̖@���x�N�g��(������)
	areal3 d_WLLVec;

	int* d_bfstM, * d_blstM, * d_nxtM;//�~���[�����N���X�g�p

	real COL;

	real r;//�e�����a
	real r2;//�e�����a��2��
	real rp;//���͗p
	real rp2;
	real rlim;//�ڋߋ֎~���a
	real rlim2;//�ڋߋ֎~���a��2��

	real n0;//�������q���x
	real lmd;//���v���V�A���W��
	real n0_grad;
	real lmd_grad;

	real FLDumax;//���E���x
	real SND;//����
	real Prs_coef;//���̈��͂����߂�Ƃ��̌W��
	real Vsc_coef;//�S�����̌W��
	real Pmax;

	///////////MPS�ϐ�////////////////
};
