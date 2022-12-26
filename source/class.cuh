#pragma once

#include "const.cuh"

class DEMPS {
public:
	///////////共通関数////////////////
	void RdDat();
	void WrtDat();//粒子データ出力
	void WrtDatWLL();//壁データ出力
	void WrtDat2();//seitei出力
	void Output(char, int);//出力
	void Output2(char, int);//seitei出力
	void AlcBkt();
	void SetPara();
	void MkBkt();
	void ClcDEMPS();
	void checkPCL();//device関数でもいいかも？
	void DevicetoHost();
	void memory_free();
	///////////共通関数////////////////


	///////////DEM関数////////////////
	void ColForce();//接触力計算
	void update();//固相粒子更新
	///////////DEM関数////////////////


	///////////MPS関数////////////////
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
	///////////MPS関数////////////////


	///////////相互作用関数////////////////
	void SLD_FLD();//相互作用力計算
	void SF_update();//粒子更新
	///////////相互作用関数////////////////


	///////////physical.txt////////////////
	treal3 MINc;//計算範囲min
	treal3 MAXc;//計算範囲max
	treal3 G;//重力
	real FIN_TIM;//計算終了タイム
	real output_time;//ファイル出力時間間隔
	real WLL_PCL_DST;//壁粒子直径
	///////////DEM_physical.txt////////////////


	///////////DEM_physical.txt////////////////
	real DEMdt;//時間刻み間隔		安定する数値に設定　MPSのタイムステップと比べて小さい方を採用
	real k;//ばね定数
	real mu;//摩擦係数
	real SLD_PCL_DST;//固体粒子直径
	///////////DEM_physical.txt////////////////


	///////////DEM_physical.txt////////////////
	real Ma;//マッハ数
	real CRT_NUM;//クーラン数
	real COL_RAT;//剛体衝突の反発係数
	real DST_LMT_RAT;//限界接近距離
	real KNM_VSC;//流体の動粘性係数
	real FLD_PCL_DST;//流体粒子直径
	///////////DEM_physical.txt////////////////





	///////////共通変数////////////////
	areal3 Pos;//座標
	areal3 d_Pos;
	areal3 Vel;//速度
	areal3 d_Vel;
	char* Typ;//粒子タイプ　
	char* d_Typ;

	real dt;//タイムステップ
	int iF = 0;//ファイル番号
	real TIM = 0.0f;//現在の計算時間
	real outtime = 0.0f;//出力間隔
	char outout_filename[256];
	FILE* fp;
	real OPT_FQC = 1.0f;//ファイル出力回数

	real DB, DB2, DBinv;//粒子探索の分割法のブロックの幅
	int nBx, nBy, nBz, nBxy, nBxyz;//バッケト用
	int* d_bfst, * d_blst, * d_nxt;//リンクリスト用

	int nP;//初期粒子数
	int nPWLL;
	int nPSLD;
	int nPFLD;
	int nPOBJ;
	int nPOBJ2;

	real PCL_DST;//バケット作る時に使う　一番大きい粒子径

	real* D;//粒子直径配列
	real* d_D;

	real Vol_SLD;//粒子体積　
	real Vol_FLD;
	///////////共通変数////////////////


	///////////DEM変数////////////////
	areal3 Ftotal;//接触力
	areal3 d_Ftotal;
	areal3 Omega;//角速度
	areal3 d_Omega;
	areal3 Torque;//トルク
	areal3 d_Torque;

	areal3 ep;
	areal3 d_ep;
	int* pair;
	int* d_pair;

	real m;//粉体1粒子質量

	real kn;//ばね定数
	real kt;

	real eta_n;//粘性減衰係数
	real eta_t;

	real I;//粉体慣性モーメント

	real ulmax;
	real usmax;
	///////////DEM変数////////////////


	///////////MPS変数////////////////
	areal3 Acc;//加速度
	areal3 d_Acc;
	real* Prs;//圧力
	real* d_Prs;
	real* pav;//圧力平均
	real* d_pav;
	real* n;//粒子数密度
	real* d_n;
	real* Dns;//粒子密度
	real* d_Dns;

	areal3 PosM;//座標
	areal3 d_PosM;
	areal3 VelM;//速度
	areal3 d_VelM;
	real* PrsM;//圧力
	real* d_PrsM;
	char* TypM;//ミラータイプ　ミラー　ゴースト
	char* d_TypM;

	int* d_FromWLL;//どの壁粒子(番号)から生成されたミラーかを保存

	char* WLLSE;
	char* d_WLLSE;//壁タイプ　surface　edge
	areal3 WLLVec;//壁の法線ベクトル(内向き)
	areal3 d_WLLVec;

	int* d_bfstM, * d_blstM, * d_nxtM;//ミラーリンクリスト用

	real COL;

	real r;//影響半径
	real r2;//影響半径の2乗
	real rp;//圧力用
	real rp2;
	real rlim;//接近禁止半径
	real rlim2;//接近禁止半径の2乗

	real n0;//初期粒子密度
	real lmd;//ラプラシアン係数
	real n0_grad;
	real lmd_grad;

	real FLDumax;//限界速度
	real SND;//音速
	real Prs_coef;//仮の圧力を求めるときの係数
	real Vsc_coef;//粘性項の係数
	real Pmax;

	///////////MPS変数////////////////
};
