#include "class.cuh"

void DEMPS::RdDat() {
	//////////////////////physical.txt///////////////////////////////
	FILE* in0;

	if (fopen_s(&in0, IN_FILE_0, "r") != 0) {
		printf("physical.txtが開けません\n");
	}
	else {
		real scan[50];
		fscanf_s(in0, "%f %f %f", &scan[0], &scan[1], &scan[2]);//最小範囲
		fscanf_s(in0, "%f %f %f", &scan[3], &scan[4], &scan[5]);//最大範囲
		fscanf_s(in0, "%f %f %f", &scan[6], &scan[7], &scan[8]);//重力加速度
		fscanf_s(in0, "%f", &scan[9]);//計算終了時間
		fscanf_s(in0, "%f", &scan[10]);//ファイル出力間隔
		fscanf_s(in0, "%f", &scan[11]);//壁粒子直径
		fclose(in0);

		WLL_PCL_DST = real(scan[11]);
		MINc.x = real(scan[0] - 3.1f * WLL_PCL_DST);
		MINc.y = real(scan[1] - 3.1f * WLL_PCL_DST);
		MINc.z = real(scan[2] - 3.1f * WLL_PCL_DST);
		MAXc.x = real(scan[3] + 3.1f * WLL_PCL_DST);
		MAXc.y = real(scan[4] + 3.1f * WLL_PCL_DST);
		MAXc.z = real(scan[5] + 3.1f * WLL_PCL_DST);
		G.x = real(scan[6]);
		G.y = real(scan[7]);
		G.z = real(scan[8]);
		FIN_TIM = real(scan[9]);
		output_time = real(scan[10]);

		printf_s("MINc.x = %f  MINc.y = %f  MINc.z = %f\n", MINc.x, MINc.y, MINc.z);
		printf_s("MAXc.x = %f  MAXc.y = %f  MAXc.z = %f\n", MAXc.x, MAXc.y, MAXc.z);
		printf_s("G.x = %f  G.y = %f  G.z = %f\n", G.x, G.y, G.z);
		printf_s("FIN_TIM = %f\n", FIN_TIM);
		printf_s("output_time = %f\n", output_time);
		printf_s("WLL_PCL_DST = %f\n\n", WLL_PCL_DST);

	}
	//////////////////////MPS_physical.txt///////////////////////////////


	//////////////////////DEM_physical.txt///////////////////////////////
	FILE* in1;

	if (fopen_s(&in1, IN_FILE_1, "r") != 0) {
		printf("DEM_physical.txtが開けません\n");
	}
	else {
		real scan[50];
		fscanf_s(in1, "%f", &scan[0]);//時間刻み幅
		fscanf_s(in1, "%f", &scan[1]);//ばね定数
		fscanf_s(in1, "%f", &scan[2]);//摩擦係数
		fscanf_s(in1, "%f", &scan[3]);//固体粒子直径
		fclose(in1);

		DEMdt = real(scan[0]);
		k = real(scan[1]);
		mu = real(scan[2]);
		SLD_PCL_DST = real(scan[3]);

		printf("DEMdt = %f\n", DEMdt);
		printf("k = %f\n", k);
		printf("mu = %f\n", mu);
		printf("SLD_PCL_DST = %f\n\n", SLD_PCL_DST);
	}
	//////////////////////DEM_physical.txt///////////////////////////////


	//////////////////////MPS_physical.txt///////////////////////////////
	FILE* in2;

	if (fopen_s(&in2, IN_FILE_2, "r") != 0) {
		printf("MPS_physical.txtが開けません\n");
	}
	else {
		real scan[50];
		fscanf_s(in2, "%f", &scan[0]);//マッハ数
		fscanf_s(in2, "%f", &scan[1]);//クーラン数
		fscanf_s(in2, "%f", &scan[2]);//反発係数
		fscanf_s(in2, "%f", &scan[3]);//限界接近距離
		fscanf_s(in2, "%f", &scan[4]);//動粘性係数
		fscanf_s(in2, "%f", &scan[5]);//流体粒子直径
		fclose(in2);

		Ma = real(scan[0]);
		CRT_NUM = real(scan[1]);
		COL_RAT = real(scan[2]);
		DST_LMT_RAT = real(scan[3]);
		KNM_VSC = real(scan[4]);
		FLD_PCL_DST = real(scan[5]);

		printf("Ma = %f\n", Ma);
		printf("CRT_NUM = %f\n", CRT_NUM);
		printf("COL_RAT= %f\n", COL_RAT);
		printf("DST_LMT_RAT = %f\n", DST_LMT_RAT);
		printf("KNM_VSC  = %f\n", KNM_VSC);
		printf("FLD_PCL_DST = %f\n\n", FLD_PCL_DST);
	}
	//////////////////////MPS_physical.txt///////////////////////////////



	//////////////////////initial.txt///////////////////////////////
	FILE* in3;
	if (fopen_s(&in3, IN_FILE_3, "r") != 0) {
		printf_s("initial.txtが開けません\n");
	}
	else {
		fscanf_s(in3, "%d", &nP);//総粒子数取得
		std::cout << "総粒子数(ファイル先頭値) nP = " << nP << std::endl;

		//粒子
		Pos.x = (real*)malloc(sizeof(real) * (nP));
		Pos.y = (real*)malloc(sizeof(real) * (nP));
		Pos.z = (real*)malloc(sizeof(real) * (nP));

		Vel.x = (real*)malloc(sizeof(real) * (nP));
		Vel.y = (real*)malloc(sizeof(real) * (nP));
		Vel.z = (real*)malloc(sizeof(real) * (nP));

		Omega.x = (real*)malloc(sizeof(real) * (nP));
		Omega.y = (real*)malloc(sizeof(real) * (nP));
		Omega.z = (real*)malloc(sizeof(real) * (nP));

		Ftotal.x = (real*)malloc(sizeof(real) * (nP));
		Ftotal.y = (real*)malloc(sizeof(real) * (nP));
		Ftotal.z = (real*)malloc(sizeof(real) * (nP));

		Torque.x = (real*)malloc(sizeof(real) * (nP));
		Torque.y = (real*)malloc(sizeof(real) * (nP));
		Torque.z = (real*)malloc(sizeof(real) * (nP));

		Acc.x = (real*)malloc(sizeof(real) * (nP));
		Acc.y = (real*)malloc(sizeof(real) * (nP));
		Acc.z = (real*)malloc(sizeof(real) * (nP));

		Prs = (real*)malloc(sizeof(real) * (nP));
		pav = (real*)malloc(sizeof(real) * (nP));

		D = (real*)malloc(sizeof(real) * (nP));
		Typ = (char*)malloc(sizeof(char) * (nP));
		Dns = (real*)malloc(sizeof(real) * (Dns_Num));
		//粒子


		//壁
		WLLVec.x = (real*)malloc(sizeof(real) * (nP));
		WLLVec.y = (real*)malloc(sizeof(real) * (nP));
		WLLVec.z = (real*)malloc(sizeof(real) * (nP));
		WLLSE = (char*)malloc(sizeof(char) * (nP));
		//壁

		//ミラー
		PosM.x = (real*)malloc(sizeof(real) * (nP * NumMRR));
		PosM.y = (real*)malloc(sizeof(real) * (nP * NumMRR));
		PosM.z = (real*)malloc(sizeof(real) * (nP * NumMRR));
		VelM.x = (real*)malloc(sizeof(real) * (nP * NumMRR));
		VelM.y = (real*)malloc(sizeof(real) * (nP * NumMRR));
		VelM.z = (real*)malloc(sizeof(real) * (nP * NumMRR));
		PrsM = (real*)malloc(sizeof(real) * (nP * NumMRR));
		TypM = (char*)malloc(sizeof(char) * (nP * NumMRR));
		//ミラー


		int nPsolid = 0;
		int nPfluid = 0;
		int nPwall = 0;
		int nPobj = 0;
		int nPobj2 = 0;
		int nPtmp = 0;
		for (int i = 0; i < nP; i++) {
			int a[1];
			float b[11];
			int c[1];
			float g[1];
			fscanf_s(in3, " %d %d %f %f %f %f %f %f %f", &a[0], &c[0], &b[0], &b[1], &b[2], &b[8], &b[9], &b[10], &g[0]);
			const treal3 pos = { b[0], b[1], b[2] };
			if (pos.x<MAXc.x && pos.x>MINc.x && pos.y<MAXc.y && pos.y>MINc.y && pos.z<MAXc.z && pos.z>MINc.z) {
				Typ[nPtmp] = char(c[0]);
				Pos.x[nPtmp] = real(b[0]); Pos.y[nPtmp] = real(b[1]); Pos.z[nPtmp] = real(b[2]);
				Vel.x[nPtmp] = Vel.y[nPtmp] = Vel.z[nPtmp] = 0.0f;
				Acc.x[nPtmp] = Acc.y[nPtmp] = Acc.z[nPtmp] = 0.0f;
				Prs[nPtmp] = 0.0f;
				pav[nPtmp] = 0.0f;

				if (Typ[nPtmp] == SLD) { nPsolid += 1; WLLVec.x[nPtmp] = 0.0f; WLLVec.y[nPtmp] = 0.0f; WLLVec.z[nPtmp] = 0.0f; D[nPtmp] = SLD_PCL_DST; }
				else if (Typ[nPtmp] == FLD) { nPfluid += 1; WLLVec.x[nPtmp] = 0.0f; WLLVec.y[nPtmp] = 0.0f; WLLVec.z[nPtmp] = 0.0f;  D[nPtmp] = FLD_PCL_DST; }
				else if (Typ[nPtmp] == WLL) { nPwall += 1; WLLVec.x[nPtmp] = real(b[8]); WLLVec.y[nPtmp] = real(b[9]); WLLVec.z[nPtmp] = real(b[10]); D[nPtmp] = WLL_PCL_DST; }
				else if (Typ[nPtmp] == OBJ) { nPobj += 1; WLLVec.x[nPtmp] = real(b[8]); WLLVec.y[nPtmp] = real(b[9]); WLLVec.z[nPtmp] = real(b[10]);  D[nPtmp] = WLL_PCL_DST; }
				else if (Typ[nPtmp] == OBJ2) { nPobj2 += 1; WLLVec.x[nPtmp] = real(b[8]); WLLVec.y[nPtmp] = real(b[9]); WLLVec.z[nPtmp] = real(b[10]);  D[nPtmp] = WLL_PCL_DST; }
				nPtmp += 1;
			}
		}
		nP = nPtmp;
		nPWLL = nPwall;
		nPSLD = nPsolid;
		nPFLD = nPfluid;
		nPOBJ = nPobj;
		nPOBJ2 = nPobj2;

		ep.x = (real*)malloc(sizeof(real) * (NCP * nPSLD));
		ep.y = (real*)malloc(sizeof(real) * (NCP * nPSLD));
		ep.z = (real*)malloc(sizeof(real) * (NCP * nPSLD));
		pair = (int*)malloc(sizeof(int) * (NCP * nPSLD));

		std::cout << "総固体粒子数 nPSLD = " << nPSLD << std::endl;
		std::cout << "総流体粒子数 nPFLD = " << nPFLD << std::endl;
		std::cout << "総壁粒子数 nPWLL = " << nPWLL << std::endl;
		std::cout << "総動壁粒子数 nPOBJ = " << nPOBJ << std::endl;
		std::cout << "総動壁弐粒子数 nPOBJ2 = " << nPOBJ2 << std::endl;
		std::cout << "総粒子数 nP = " << nP << std::endl;
	}
	fclose(in3);
	//////////////////////initial.txt///////////////////////////////
}


__global__ void d_initialize_int_array(const int n, int* i_array, const int a) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= n) { return; }
	i_array[i] = a;
}

__global__ void d_initialize_real_array(const int n, real* i_array, real a) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= n) { return; }
	i_array[i] = a;
}


void DEMPS::Output(const char typ, const int outputflg) {
	if (outputflg == 0)//固体粒子出力
	{
		sprintf(outout_filename, "./output/outputSLD%05d.csv", iF);
		printf("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,V.x,V.y,V.z,Vel,F.x,F.y,F.z,Force\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", Pos.x[i], Pos.y[i], Pos.z[i],
						Vel.x[i], Vel.y[i], Vel.z[i], sqrt(Vel.x[i] * Vel.x[i] + Vel.y[i] * Vel.y[i] + Vel.z[i] * Vel.z[i]),
						Ftotal.x[i], Ftotal.y[i], Ftotal.z[i], sqrt(Ftotal.x[i] * Ftotal.x[i] + Ftotal.y[i] * Ftotal.y[i] + Ftotal.z[i] * Ftotal.z[i]));
				}
			}
		}
		fclose(fp);
	}

	else if (outputflg == 1)//流体粒子出力
	{
		sprintf_s(outout_filename, "./output/outputFLD%05d.csv", iF);
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Vel.x,Vel.y,Vel.z,Vel,Prs,pav\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f\n", Pos.x[i], Pos.y[i], Pos.z[i], Vel.x[i], Vel.y[i], Vel.z[i], sqrt(Vel.x[i] * Vel.x[i] + Vel.y[i] * Vel.y[i] + Vel.z[i] * Vel.z[i]), Prs[i], pav[i] / OPT_FQC);
				}
			}
		}
		fclose(fp);
	}

	else if (outputflg == 2)//壁粒子出力　将来的には壁に働く力(固体衝突力・水圧力)を出力したい
	{
		sprintf(outout_filename, "./output/outputWLL%05d.csv", iF);
		printf("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Prs,pav,SE\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%d\n", Pos.x[i], Pos.y[i], Pos.z[i], Prs[i], pav[i] / OPT_FQC, WLLSE[i]);
				}
			}
		}
		fclose(fp);
	}

	else if (outputflg == 3)//動壁粒子出力　将来的には壁に働く力(固体衝突力・水圧力)を出力したい
	{
		sprintf(outout_filename, "./output/outputOBJ%05d.csv", iF);
		printf("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Prs,pav,SE\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%d\n", Pos.x[i], Pos.y[i], Pos.z[i], Prs[i], pav[i] / OPT_FQC, WLLSE[i]);
				}
			}
		}
		fclose(fp);
	}

	else if (outputflg == 4)//動壁粒子弐出力　将来的には壁に働く力(固体衝突力・水圧力)を出力したい
	{
		sprintf(outout_filename, "./output/outputOBJ2%05d.csv", iF);
		printf("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "Pos.x,Pos.y,Pos.z,Prs,pav,SE\n");
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%f,%f,%f,%f,%f,%d\n", Pos.x[i], Pos.y[i], Pos.z[i], Prs[i], pav[i] / OPT_FQC, WLLSE[i]);
				}
			}
		}
		fclose(fp);
	}

	else if (outputflg == 5)//ミラー粒子出力
	{
		sprintf_s(outout_filename, "./output/outputMRR%05d.csv", iF);
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf_s("%sが開けません\n", outout_filename);
		}
		else {
			fprintf_s(fp, "PosM.x,PosM.y,PosM.z,VelM,PrsM\n");
			for (int i = 0; i < nP; i++) {
				for (int k = 0; k < NumMRR; k++) {
					int kiNM = k + i * NumMRR;
					if (TypM[kiNM] == typ) {
						fprintf_s(fp, "%f,%f,%f,%f,%f,\n", PosM.x[kiNM], PosM.y[kiNM], PosM.z[kiNM], sqrt(VelM.x[kiNM] * VelM.x[kiNM] + VelM.y[kiNM] * VelM.y[kiNM] + VelM.z[kiNM] * VelM.z[kiNM]), PrsM[kiNM]);
					}
				}
			}
		}
		fclose(fp);
	}

}


void DEMPS::Output2(const char typ, const int outputflg) {//seitei.txt出力
	if (outputflg == 0)//SLD_seitei
	{
		sprintf(outout_filename, "./output/DEM_seitei.txt");
		printf("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf("%sが開けません\n", outout_filename);
		}
		else {
			int sldnp = 0;
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) { sldnp += 1; }//固体粒子カウント
			}
			fprintf(fp, "%d\n", sldnp);
			int SLDnP = 0;
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf(fp, "%d %d %f %f %f %f %f %f %f\n", SLDnP, Typ[i], Pos.x[i], Pos.y[i], Pos.z[i], 0.0f, 0.0f, 0.0f, 0.0f);
					SLDnP += 1;
				}
			}
		}
		fclose(fp);
	}

	else if (outputflg == 1)//FLD_seitei
	{
		sprintf_s(outout_filename, "./output/MPS_seitei.txt");
		printf_s("Filename = %s\n", outout_filename);

		if (fopen_s(&fp, outout_filename, "w") != 0) {
			printf_s("%sが開けません\n", outout_filename);
		}
		else {
			int fldnp = 0;
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) { fldnp += 1; }//流体粒子カウント
			}
			fprintf_s(fp, "%d\n", fldnp);
			int FLDnP = 0;
			for (int i = 0; i < nP; i++) {
				if (Typ[i] == typ) {
					fprintf_s(fp, "%d %d %f %f %f %f %f %f %f\n", FLDnP, Typ[i], Pos.x[i], Pos.y[i], Pos.z[i], 0.0f, 0.0f, 0.0f, 0.0f);
					FLDnP += 1;
				}
			}
		}
		fclose(fp);
	}

}


void DEMPS::WrtDat(void) {
#if DEM_flg
	Output(SLD, 0);//固体粒子出力
#endif

#if MPS_flg
	Output(FLD, 1);//水粒子出力
	Output(MRR, 5);//生成ミラー確認
#endif

#if move_WLL
	Output(WLL, 2);//水粒子出力
	Output(OBJ, 3);//水粒子出力
	Output(OBJ2, 4);//水粒子出力
#endif

	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	((d_initialize_real_array << <blocks_nP, threads >> > (nP, d_pav, 0.0f)));//平均圧力初期化
	//DEMの力も時間平均とりたい
	//CHECK(cudaDeviceSynchronize());
	printf_s("WrtDat finished!\n\n");
}


void DEMPS::WrtDatWLL(void) {
	Output(WLL, 2);
	printf("WrtDatWLL finished!\n\n");
}


void DEMPS::WrtDat2(void) {//seitei.txt
	Output2(SLD, 0);
	Output2(FLD, 1);
	printf_s("WrtDat2 finished!\n\n");
}


void DEMPS::AlcBkt() {

	PCL_DST = WLL_PCL_DST;
	if (SLD_PCL_DST > PCL_DST) { PCL_DST = SLD_PCL_DST; }
	if (FLD_PCL_DST > PCL_DST) { PCL_DST = FLD_PCL_DST; }//でかい方に合わせてバケット作る

	r = FLD_PCL_DST * 3.1f;
	r2 = r * r;
	rp = FLD_PCL_DST * 2.1f;//圧力用
	rp2 = rp * rp;

	DB = PCL_DST * 3.1f;
	DB2 = DB * DB;
	DBinv = 1.0f / DB;

	nBx = (int)((MAXc.x - MINc.x) * DBinv) + 3;
	nBy = (int)((MAXc.y - MINc.y) * DBinv) + 3;
	nBz = (int)((MAXc.z - MINc.z) * DBinv) + 3;

	nBxy = nBx * nBy;
	nBxyz = nBx * nBy * nBz;
	printf_s("nBx:%d  nBy:%d  nBz:%d  nBxy:%d  nBxyz:%d\n", nBx, nBy, nBz, nBxy, nBxyz);

	(cudaMalloc((void**)&d_bfst, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_blst, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_nxt, sizeof(int) * (nP)));
	(cudaMalloc((void**)&d_bfstM, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_blstM, sizeof(int) * nBxyz));
	(cudaMalloc((void**)&d_nxtM, sizeof(int) * (nP * NumMRR)));
}


void DEMPS::SetPara() {
	///////////////////////////////////////////////////////共通///////////////////////////////////////////////////////////////////

	//音速とタイムステップ
	real max_height = MINc.y;
	real min_height = MAXc.y;
	for (int i = 0; i < nP; i++) {//初期化
		if (Typ[i] == FLD) {
			if (max_height < Pos.y[i]) { max_height = Pos.y[i]; }
			else if (min_height > Pos.y[i]) { min_height = Pos.y[i]; }
		}
	}
	real lmm = max_height - min_height;
	ulmax = sqrt(2.0f * abs(G.y) * lmm);

	real max_height2 = MINc.y;
	real min_height2 = MAXc.y;
	for (int i = 0; i < nP; i++) {//初期化
		if (Typ[i] == SLD) {
			if (max_height2 < Pos.y[i]) { max_height2 = Pos.y[i]; }
			else if (min_height2 > Pos.y[i]) { min_height2 = Pos.y[i]; }
		}
	}
	usmax = sqrt(2.0f * abs(G.y) * (max_height2 - min_height2));


#pragma omp parallel for
	for (int i = 0; i < Dns_Num; i++) {//初期化
		Dns[i] = 0.0f;
	}
	Dns[SLD] = Dns_SLD;
	Dns[FLD] = Dns_FLD;
	Dns[WLL] = Dns_WLL;
	Dns[OBJ] = Dns_OBJ;
	Dns[OBJ2] = Dns_OBJ2;
	Dns[MRR] = Dns_FLD;

#pragma omp parallel for
	for (int i = 0; i < nP; i++) {
		Ftotal.x[i] = Ftotal.y[i] = Ftotal.z[i] = 0.0f;
		Omega.x[i] = Omega.y[i] = Omega.z[i] = 0.0f;
		Torque.x[i] = Torque.y[i] = Torque.z[i] = 0.0f;

		if (Typ[i] == SLD) { D[i] = SLD_PCL_DST; }
		else if (Typ[i] == FLD) { D[i] = FLD_PCL_DST; }
		else if (Typ[i] == WLL) { D[i] = WLL_PCL_DST; }
		else if (Typ[i] == OBJ) { D[i] = WLL_PCL_DST; }
		else if (Typ[i] == OBJ2) { D[i] = WLL_PCL_DST; }
	}

	Vol_SLD = pi * SLD_PCL_DST * SLD_PCL_DST * SLD_PCL_DST / 6.0f;
	Vol_FLD = pi * FLD_PCL_DST * FLD_PCL_DST * FLD_PCL_DST / 6.0f;

	///////////////////////////////////////////////////////共通///////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////DEM///////////////////////////////////////////////////////////////////

	m = Dns[SLD] * pi * SLD_PCL_DST * SLD_PCL_DST * SLD_PCL_DST / 6.0f;

	I = 0.1f * m * SLD_PCL_DST * SLD_PCL_DST;

#if 1
	real alpha = 0.27;

	kn = k;
	kt = kn / (2.0f * (1.0f + alpha));

	eta_n = 2.0f * sqrt(m * kn);
	eta_t = eta_n / sqrt(2.0f * (1.0f + alpha));
#else
	kn = 1000;
	kt = kn * 0.25f;

	eta_n = sqrt(2.0f * m * kn);
	eta_t = sqrt(2.0f * m * kt);

#endif


#pragma omp parallel for
	for (int i = 0; i < nPSLD; i++) { for (int k = 0; k < NCP; k++) { ep.x[k + i * NCP] = 0.0f; } }
#pragma omp parallel for
	for (int i = 0; i < nPSLD; i++) { for (int k = 0; k < NCP; k++) { ep.y[k + i * NCP] = 0.0f; } }
#pragma omp parallel for
	for (int i = 0; i < nPSLD; i++) { for (int k = 0; k < NCP; k++) { ep.z[k + i * NCP] = 0.0f; } }
#pragma omp parallel for
	for (int i = 0; i < nPSLD; i++) { for (int k = 0; k < NCP; k++) { pair[k + i * NCP] = -2; } }

	printf("m:%.10f\nkn:%f  kt:%f\neta_n:%f  eta_t:%f\n\n", m, kn, kt, eta_n, eta_t);

	///////////////////////////////////////////////////////DEM///////////////////////////////////////////////////////////////////	



	///////////////////////////////////////////////////////MPS///////////////////////////////////////////////////////////////////

	//初期粒子密度
	real tn0 = 0.0f;
	real tn0_grad = 0.0f;
	real tlmd = 0.0f;
	real tlmd_grad = 0.0f;

	for (int ix = -10; ix < 10; ix++) {
		for (int iy = -10; iy < 10; iy++) {
			for (int iz = -10; iz < 10; iz++) {
				real x = real(FLD_PCL_DST) * real(ix);
				real y = real(FLD_PCL_DST) * real(iy);
				real z = real(FLD_PCL_DST) * real(iz);
				real dist2 = x * x + y * y + z * z;
				if (dist2 == 0.0f) { continue; }
				real dist = sqrt(dist2);
				if (dist2 < rp2) {
					tn0_grad += WEI_grad(dist, rp);
					tlmd_grad += dist2 * WEI_grad(dist, rp);
				}
				if (dist2 < r2) {
					tn0 += WEI(dist, r);
					tlmd += dist2 * WEI(dist, r);
				}
			}
		}
	}
	n0 = tn0;
	lmd = tlmd / tn0;
	n0_grad = tn0_grad;
	lmd_grad = tlmd_grad / tn0_grad;
	printf_s("n0:%f\nn0_grad:%f\n", n0, n0_grad);

	//接近禁止距離
	rlim = FLD_PCL_DST * DST_LMT_RAT;
	rlim2 = rlim * rlim;

	COL = real(COL_RAT + 1.0f);//反発係数設定

#pragma omp parallel for
	for (int i = 0; i < nP; i++) {//初期化
		for (int k = 0; k < NumMRR; k++) {
			int kiNM = k + i * NumMRR;
			PosM.x[kiNM] = PosM.y[kiNM] = PosM.z[kiNM] = 0.0f;
			VelM.x[kiNM] = VelM.y[kiNM] = VelM.z[kiNM] = 0.0f;
			PrsM[kiNM] = 0.0f;
			TypM[kiNM] = GST;
		}
		WLLSE[i] = Surface;
	}



	SND = ulmax / Ma;
	Prs_coef = SND * SND / n0_grad;
	Vsc_coef = 2.0f * 3.0f * KNM_VSC / n0 / lmd;
	Pmax = 1000.0f * abs(G.y) * lmm;//9.80665f 

	///////////////////////////////////////////////////////MPS///////////////////////////////////////////////////////////////////


	//dt = CRT_NUM * PCL_DST / SND;
	dt = 0.00001f;//固定
	//dt = real(PCL_DST / (ulmax + SND));//計算間隔(ダムブレイクの最大速度を基準にする)
	if (DEMdt < dt) { dt = DEMdt; }

	printf_s("タイムステップ:dt=%f\n最大流体速度:ulmax=%f\n最大固体速度:usmax=%f\n音速=%f\nVsc_coef=%f\nPrs_coef=%f\nPmax=%f\n", dt, ulmax, usmax, SND, Vsc_coef, Prs_coef, Pmax);


	//配列確保とVRAMへの転送
	(cudaMalloc((void**)&d_Typ, sizeof(char) * nP));
	(cudaMalloc((void**)&d_Pos.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Pos.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Pos.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Vel.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Vel.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Vel.z, sizeof(real) * nP));

	(cudaMemcpy(d_Typ, Typ, sizeof(char) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Pos.x, Pos.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Pos.y, Pos.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Pos.z, Pos.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Vel.x, Vel.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Vel.y, Vel.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Vel.z, Vel.z, sizeof(real) * nP, cudaMemcpyHostToDevice));

	///////////////////////////////////////////////////////DEM///////////////////////////////////////////////////////////////////
	(cudaMalloc((void**)&d_D, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Ftotal.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Ftotal.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Ftotal.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Omega.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Omega.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Omega.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Torque.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Torque.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Torque.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_ep.x, sizeof(real) * (NCP * nPSLD)));
	(cudaMalloc((void**)&d_ep.y, sizeof(real) * (NCP * nPSLD)));
	(cudaMalloc((void**)&d_ep.z, sizeof(real) * (NCP * nPSLD)));
	(cudaMalloc((void**)&d_pair, sizeof(int) * (NCP * nPSLD)));


	(cudaMemcpy(d_D, D, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Ftotal.x, Ftotal.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Ftotal.y, Ftotal.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Ftotal.z, Ftotal.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Omega.x, Omega.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Omega.y, Omega.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Omega.z, Omega.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Torque.x, Torque.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Torque.y, Torque.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Torque.z, Torque.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_ep.x, ep.x, sizeof(real) * (NCP * nPSLD), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_ep.y, ep.y, sizeof(real) * (NCP * nPSLD), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_ep.z, ep.z, sizeof(real) * (NCP * nPSLD), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_pair, pair, sizeof(int) * (NCP * nPSLD), cudaMemcpyHostToDevice));
	///////////////////////////////////////////////////////DEM///////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////MPS///////////////////////////////////////////////////////////////////
	(cudaMalloc((void**)&d_Acc.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Acc.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Acc.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Prs, sizeof(real) * nP));
	(cudaMalloc((void**)&d_pav, sizeof(real) * nP));
	(cudaMalloc((void**)&d_Dns, sizeof(real) * Dns_Num));

	(cudaMalloc((void**)&d_WLLVec.x, sizeof(real) * nP));
	(cudaMalloc((void**)&d_WLLVec.y, sizeof(real) * nP));
	(cudaMalloc((void**)&d_WLLVec.z, sizeof(real) * nP));
	(cudaMalloc((void**)&d_WLLSE, sizeof(char) * nP));

	(cudaMalloc((void**)&d_FromWLL, sizeof(int) * (nP * NumMRR)));//デバイスのみ(転送なし)

	(cudaMalloc((void**)&d_TypM, sizeof(char) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PosM.x, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PosM.y, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PosM.z, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_VelM.x, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_VelM.y, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_VelM.z, sizeof(real) * (nP * NumMRR)));
	(cudaMalloc((void**)&d_PrsM, sizeof(real) * (nP * NumMRR)));


	(cudaMemcpy(d_Acc.x, Acc.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Acc.y, Acc.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Acc.z, Acc.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Prs, Prs, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_pav, pav, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_Dns, Dns, sizeof(real) * Dns_Num, cudaMemcpyHostToDevice));

	(cudaMemcpy(d_WLLVec.x, WLLVec.x, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_WLLVec.y, WLLVec.y, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_WLLVec.z, WLLVec.z, sizeof(real) * nP, cudaMemcpyHostToDevice));
	(cudaMemcpy(d_WLLSE, WLLSE, sizeof(char) * nP, cudaMemcpyHostToDevice));

	(cudaMemcpy(d_TypM, TypM, sizeof(char) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PosM.x, PosM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PosM.y, PosM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PosM.z, PosM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_VelM.x, VelM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_VelM.y, VelM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_VelM.z, VelM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	(cudaMemcpy(d_PrsM, PrsM, sizeof(real) * (nP * NumMRR), cudaMemcpyHostToDevice));
	///////////////////////////////////////////////////////MPS///////////////////////////////////////////////////////////////////

}


__global__ void d_MkBkt(const int nP, const  int nBx, const  int nBxy, const  real DBinv, int* d_bfst, int* d_blst, int* d_nxt, const char* d_Typ, areal3 d_Pos, const treal3 MINc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	if (d_Typ[i] == GST) { return; }

	int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
	int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
	int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
	int ib = iz * nBxy + iy * nBx + ix;
	const int j = atomicExch(&d_blst[ib], i);
	if (j == -1) { d_bfst[ib] = i; }
	else { d_nxt[j] = i; }

}


void DEMPS::MkBkt() {//粒子をバケットに収納
	//printf_s("MkBkt start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_bfst, -1)));
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_blst, -1)));
	((d_initialize_int_array << <blocks_nP, threads >> > (nP, d_nxt, -1)));
	//CHECK(cudaDeviceSynchronize());

	d_MkBkt << <blocks_nP, threads >> > (nP, nBx, nBxy, DBinv, d_bfst, d_blst, d_nxt, d_Typ, d_Pos, MINc);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("MkBkt finished!\n\n");
}


#if MPS_flg
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////MPS_Function/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void d_Surface_Edge(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_WLLVec, char* d_WLLSE, const real WLL_PCL_DST,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)//壁粒子surface edge判定
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	real r2 = (1.0f * WLL_PCL_DST) * (1.0f * WLL_PCL_DST);//隣接壁粒子判定距離の2乗
	char Typ = d_Typ[i];

	if (Typ == WLL) {
		treal3 pos;
		pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
		treal3 vec;
		vec.x = d_WLLVec.x[i];	vec.y = d_WLLVec.y[i];	vec.z = d_WLLVec.z[i];

		int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
		int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
		int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							if (d_Typ[j] == WLL) {
								treal3 p;//i,jの距離の成分
								p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];
								real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
								if (dist2 < r2) {
									if ((d_WLLVec.x[j] != vec.x) || (d_WLLVec.y[j] != vec.y) || (d_WLLVec.z[j] != vec.z)) {//近くの粒子と法線ベクトル違うならedge
										d_WLLSE[i] = Edge;
									}
								}
							}
						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
	}

	else if (Typ == OBJ) {
		treal3 pos;
		pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
		treal3 vec;
		vec.x = d_WLLVec.x[i];	vec.y = d_WLLVec.y[i];	vec.z = d_WLLVec.z[i];

		int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
		int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
		int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							if (d_Typ[j] == OBJ) {
								treal3 p;//i,jの距離の成分
								p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];
								real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
								if (dist2 < r2) {
									if ((d_WLLVec.x[j] != vec.x) || (d_WLLVec.y[j] != vec.y) || (d_WLLVec.z[j] != vec.z)) {//近くの粒子と法線ベクトル違うならedge
										d_WLLSE[i] = Edge;
									}
								}
							}
						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
	}

	else if (Typ == OBJ2) {
		treal3 pos;
		pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
		treal3 vec;
		vec.x = d_WLLVec.x[i];	vec.y = d_WLLVec.y[i];	vec.z = d_WLLVec.z[i];

		int ix = (int)((d_Pos.x[i] - MINc.x) * DBinv) + 1;
		int iy = (int)((d_Pos.y[i] - MINc.y) * DBinv) + 1;
		int iz = (int)((d_Pos.z[i] - MINc.z) * DBinv) + 1;
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							if (d_Typ[j] == OBJ2) {
								treal3 p;//i,jの距離の成分
								p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];
								real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
								if (dist2 < r2) {
									if ((d_WLLVec.x[j] != vec.x) || (d_WLLVec.y[j] != vec.y) || (d_WLLVec.z[j] != vec.z)) {//近くの粒子と法線ベクトル違うならedge
										d_WLLSE[i] = Edge;
									}
								}
							}
						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
	}

}


void DEMPS::Surface_Edge() {//surface-edge判定
	//printf_s("Surface_Edge  start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_Surface_Edge << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_WLLVec, d_WLLSE, WLL_PCL_DST, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	//CHECK(cudaDeviceSynchronize());

	cudaMemcpy(WLLSE, d_WLLSE, sizeof(char) * nP, cudaMemcpyDeviceToHost);

	//printf_s("Surface_Edge finished!\n\n");
}


__global__ void d_ResetMRR(const int nP_NumMRR, char* d_TypM, areal3 d_PosM, areal3 d_VelM, real* d_PrsM, int* d_FromWLL)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP_NumMRR) { return; }
	d_TypM[i] = GST;
	d_PosM.x[i] = d_PosM.y[i] = d_PosM.z[i] = 0.0f;
	d_VelM.x[i] = d_VelM.y[i] = d_VelM.z[i] = 0.0f;
	d_PrsM[i] = 0.0f;
	d_FromWLL[i] = -1;
}


void DEMPS::ResetMRR() {//ミラー粒子削除
	//printf_s("ResetMRR start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP * NumMRR);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP_NumMRR(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_ResetMRR << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), d_TypM, d_PosM, d_VelM, d_PrsM, d_FromWLL);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("ResetMRR finished!\n\n");
}


__global__ void d_GenMRR_nonslip(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, real* d_Prs, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, areal3 d_VelM, real* d_PrsM, int* d_FromWLL, char* d_WLLSE, const real r2,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	char Typ = d_Typ[i];
	if (Typ != FLD) { return; }//流体粒子から面対称ミラー生成

	real r2_2 = 1.0f * r2;
	treal3 pos;
	pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
	int WLLexist = 0;//近傍に壁粒子が存在したら１
	int iNM = i * NumMRR;
	int Edgeexist = 0;
	int Edge_unique[NumMRR];
	for (int k = 0; k < NumMRR; k++) { Edge_unique[k] = -1; }

	//最近傍壁粒子探索
	int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
	int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
	int iz = (int)((pos.z - MINc.z) * DBinv) + 1;
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						char typ = d_Typ[j];
						if ((typ == WLL) || (typ == OBJ) || (typ == OBJ2)) {
							treal3 p;
							p.x = pos.x - d_Pos.x[j];	p.y = pos.y - d_Pos.y[j];	p.z = pos.z - d_Pos.z[j];//壁粒子→流体粒子のベクトル
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
							if (dist2 <= r2_2) {//影響半径内に壁有り
								WLLexist = 1;//壁有りフラグ
								treal3 wallvec;
								wallvec.x = d_WLLVec.x[j];	wallvec.y = d_WLLVec.y[j];	wallvec.z = d_WLLVec.z[j];

								int unique = 0;//面のミラー生成
								for (int k = 0; k < NumMRR; k++) {//unique判定(iからk番目に生成されたミラーと法線ベクトルが重複していないか判定)
									int jNum = d_FromWLL[k + iNM];
									if (jNum != -1) {
										if ((d_WLLVec.x[jNum] == wallvec.x) && (d_WLLVec.y[jNum] == wallvec.y) && (d_WLLVec.z[jNum] == wallvec.z)) {
											unique = 1;//同じ法線ベクトルが登録済み													
										}
									}
								}
								if (unique == 0) {//jがuniqueであれば登録
									for (int k = 0; k < NumMRR; k++) {
										int kiNM = k + iNM;
										if (d_FromWLL[kiNM] == -1) {
											d_FromWLL[kiNM] = j;
											break;
										}
									}
								}

								if (d_WLLSE[j] == Edge) {//角のミラー生成
									Edgeexist = 1;
									int E_unique = 0;
									for (int k = 0; k < NumMRR; k++) {
										int q = Edge_unique[k];
										if (q != -1) {
											if ((d_WLLVec.x[q] == wallvec.x) && (d_WLLVec.y[q] == wallvec.y) && (d_WLLVec.z[q] == wallvec.z)) {
												E_unique = 1;//同じ法線ベクトルが登録済み														
											}
										}
									}
									if (E_unique == 0) {
										for (int k = 0; k < NumMRR; k++) {
											if (Edge_unique[k] == -1) {
												Edge_unique[k] = j;
												break;
											}
										}
									}
								}

							}
						}
					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

	//ミラー粒子生成
	if (WLLexist == 1) {//近くに壁がある

		for (int k = 0; k < NumMRR; k++) {
			int kiNM = k + iNM;
			int FromNum = d_FromWLL[kiNM];//壁粒子番号jをFromNumとしてレジスタに登録
			if (FromNum == -1) { continue; }//壁粒子j(FromNum)からミラー生成	
			treal3 posw;
			posw.x = d_Pos.x[FromNum];	posw.y = d_Pos.y[FromNum];	posw.z = d_Pos.z[FromNum];
			treal3 WLLvec;
			WLLvec.x = d_WLLVec.x[FromNum];	WLLvec.y = d_WLLVec.y[FromNum];	WLLvec.z = d_WLLVec.z[FromNum];
			treal3 PWvec;//最近傍壁粒子と流体粒子iの相対座標ベクトル
			PWvec.x = pos.x - posw.x;		PWvec.y = pos.y - posw.y;		PWvec.z = pos.z - posw.z;

			real PW_WLL = PWvec.x * WLLvec.x + PWvec.y * WLLvec.y + PWvec.z * WLLvec.z;
			if (PW_WLL < 0) { continue; }//内外判定　内側にはミラー生成しない(流体が壁の外にある)

			real distance = PW_WLL / sqrt(WLLvec.x * WLLvec.x + WLLvec.y * WLLvec.y + WLLvec.z * WLLvec.z);//距離のスカラー値 法線ベクトルなので分母は１になってるはず
			d_TypM[kiNM] = MRR;
			d_PosM.x[kiNM] = posw.x + PWvec.x - 2.0f * distance * WLLvec.x;//法線ベクトルを用いて対称位置に生成
			d_PosM.y[kiNM] = posw.y + PWvec.y - 2.0f * distance * WLLvec.y;
			d_PosM.z[kiNM] = posw.z + PWvec.z - 2.0f * distance * WLLvec.z;
			d_VelM.x[kiNM] = -d_Vel.x[i];//non-slip
			d_VelM.y[kiNM] = -d_Vel.y[i];
			d_VelM.z[kiNM] = -d_Vel.z[i];
			d_PrsM[kiNM] = d_Prs[i];
			/*d_VelM.x[kiNM] = 0.0f;
			d_VelM.y[kiNM] = 0.0f;
			d_VelM.z[kiNM] = 0.0f;*/
			//d_PrsM[kiNM] = 0.0f;
		}//ここまで面ミラー生成

#if 0 //Edge生成フラグ
		if (Edgeexist == 1) {
			int NumEdge = 0;
			for (int k = 0; k < NumMRR; k++) {
				if (Edge_unique[k] != -1) { NumEdge++; }
			}

			if (NumEdge >= 2) {//Edge粒子が影響半径内に2種類以上あったら合成ベクトルから角のミラー生成
				treal3 posw;
				treal3 synvec;
				treal3 PWvec;
				real distance;
				for (int k = 0; k < NumMRR; k++) {//合成ベクトル計算+最後のedgeを生成元に設定	
					int FromNum = Edge_unique[k];
					if (FromNum == -1) { continue; }
					posw.x = d_Pos.x[FromNum];	posw.y = d_Pos.y[FromNum];	posw.z = d_Pos.z[FromNum];
					synvec.x += d_WLLVec.x[FromNum];		synvec.y += d_WLLVec.y[FromNum];		synvec.z += d_WLLVec.z[FromNum];
					PWvec.x = pos.x - posw.x;		PWvec.y = pos.y - posw.y;		PWvec.z = pos.z - posw.z;
				}

				real abs_synvec = sqrt(synvec.x * synvec.x + synvec.y * synvec.y + synvec.z * synvec.z);
				synvec.x /= abs_synvec;	synvec.y /= abs_synvec;	synvec.z /= abs_synvec;//合成ベクトルの大きさを１にする

				real PW_WLL = PWvec.x * synvec.x + PWvec.y * synvec.y + PWvec.z * synvec.z;
				if (PW_WLL < 0) { return; }//内外判定　内側にはミラー生成しない(流体が壁の外にある)

				distance = PW_WLL / sqrt(synvec.x * synvec.x + synvec.y * synvec.y + synvec.z * synvec.z);//距離のスカラー値 法線ベクトルなので分母は１になってるはず

				int MRR_space = -1;
				for (int k = 0; k < NumMRR; k++) {
					int kiNM = k + iNM;
					if (d_FromWLL[kiNM] == -1) { MRR_space = kiNM; break; }
				}

				d_TypM[MRR_space] = MRR;
				d_PosM.x[MRR_space] = posw.x + PWvec.x - 2.0f * distance * synvec.x;//法線ベクトルを用いて対称位置に生成
				d_PosM.y[MRR_space] = posw.y + PWvec.y - 2.0f * distance * synvec.y;
				d_PosM.z[MRR_space] = posw.z + PWvec.z - 2.0f * distance * synvec.z;
				/*d_VelM.x[MRR_space] = -d_Vel.x[i];
				d_VelM.y[MRR_space] = -d_Vel.y[i];
				d_VelM.z[MRR_space] = -d_Vel.z[i];
				d_PrsM[MRR_space] = d_Prs[i];*/
				d_VelM.x[MRR_space] = 0.0f;
				d_VelM.y[MRR_space] = 0.0f;
				d_VelM.z[MRR_space] = 0.0f;
				d_PrsM[MRR_space] = d_Prs[i];
			}
		}
#endif
	}

}


void DEMPS::GenMRR_nonslip() {//ミラー全成分生成
	//printf_s("GenMRR_nonslip start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP * NumMRR);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nPMRR(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_GenMRR_nonslip << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Prs, d_WLLVec, d_TypM, d_PosM, d_VelM, d_PrsM, d_FromWLL, d_WLLSE, r2, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("GenMRR_nonslip finished!\n\n");
}


__global__ void d_MkBkt_MRR(const int nP_NumMRR, const  int nBx, const  int nBxy, const  real DBinv, int* d_bfstM, int* d_blstM, int* d_nxtM, const char* d_TypM, const  areal3 d_PosM, const treal3 MINc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP_NumMRR) { return; }
	if (d_TypM[i] == GST) { return; }
	int ix = (int)((d_PosM.x[i] - MINc.x) * DBinv) + 1;
	int iy = (int)((d_PosM.y[i] - MINc.y) * DBinv) + 1;
	int iz = (int)((d_PosM.z[i] - MINc.z) * DBinv) + 1;
	int ib = iz * nBxy + iy * nBx + ix;
	const int j = atomicExch(&d_blstM[ib], i);
	if (j == -1) { d_bfstM[ib] = i; }
	else { d_nxtM[j] = i; }

}

//壁の内側にできたミラーをはじく
__global__ void d_MRRinout(const int nP_NumMRR, const  int nBx, const  int nBxy, const  real DBinv, int* d_bfst, int* d_blst, int* d_nxt, int* d_bfstM, int* d_blstM, int* d_nxtM, char* d_TypM, const  areal3 d_PosM, const char* d_Typ, const areal3 d_Pos, const areal3 d_WLLVec, const real r, const treal3 MINc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP_NumMRR) { return; }

	treal3 posM;
	posM.x = d_PosM.x[i];	posM.y = d_PosM.y[i];	posM.z = d_PosM.z[i];
	int WLLexist = 0;
	real walldist = 10.0f * r;
	treal3 posw = { 0.0f };
	treal3 wallvec = { 0.0f };

	int ix = (int)((posM.x - MINc.x) * DBinv) + 1;
	int iy = (int)((posM.y - MINc.y) * DBinv) + 1;
	int iz = (int)((posM.z - MINc.z) * DBinv) + 1;

	//最近傍壁粒子探索
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始

					char typ = d_Typ[j];
					if ((typ == WLL) || (typ == OBJ) || (typ == OBJ2)) {
						treal3 posj;
						posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];
						treal3 p;
						p.x = posM.x - posj.x;	p.y = posM.y - posj.y;	p.z = posM.z - posj.z;
						real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
						real dist = sqrt(dist2);
						if (dist <= r) {//影響半径内に壁有り
							if (dist < walldist) {
								WLLexist = 1;//壁有りフラグ jにする
								walldist = dist;
								wallvec.x = d_WLLVec.x[j];	wallvec.y = d_WLLVec.y[j];	wallvec.z = d_WLLVec.z[j];
								posw.x = posj.x;	posw.y = posj.y;	posw.z = posj.z;
							}
						}
					}

					j = d_nxt[j];
					if (j == -1) break;
				}
			}
		}
	}

	real inout = (posM.x - posw.x) * wallvec.x + (posM.y - posw.y) * wallvec.y + (posM.z - posw.z) * wallvec.z;//壁粒子→ミラー粒子のベクトル * 壁法線ベクトル
	if ((WLLexist == 1) && (inout >= 0)) { d_TypM[i] = GST; }//内外判定に引っかかったらGSTに

}

void DEMPS::MkBkt_MRR() {//粒子をバケットに収納
	//printf_s("MkBkt_MRR start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP * NumMRR);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP_NumMRR(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_bfstM, -1)));
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_blstM, -1)));
	((d_initialize_int_array << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), d_nxtM, -1)));
	d_MkBkt_MRR << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), nBx, nBxy, DBinv, d_bfstM, d_blstM, d_nxtM, d_TypM, d_PosM, MINc);

	d_MRRinout << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), nBx, nBxy, DBinv, d_bfst, d_blst, d_nxt, d_bfstM, d_blstM, d_nxtM, d_TypM, d_PosM, d_Typ, d_Pos, d_WLLVec, r, MINc);

	/*((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_bfstM, -1)));//inoutでGSTになったやつはバケットに入れない　＝＞　呼び出し側でGSTならcontinueの処理入れた
	((d_initialize_int_array << <blocks_nBxyz, threads >> > (nBxyz, d_blstM, -1)));
	((d_initialize_int_array << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), d_nxtM, -1)));
	d_MkBkt_MRR << <blocks_nP_NumMRR, threads >> > ((nP * NumMRR), nBx, nBxy, DBinv, d_bfstM, d_blstM, d_nxtM, d_TypM, d_PosM, MINc);*/
	//CHECK(cudaDeviceSynchronize());

	//printf_s("MkBkt_MRR finished!\n\n");
}


__global__ void d_VscTrm(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, areal3 d_VelM, const real r, const real PCL_DST, const real n0, const real KNM_VSC, const real Vsc_coef, treal3 G,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	if (d_Typ[i] != FLD) { return; }

	treal3 pos;		treal3 vel;
	pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
	vel.x = d_Vel.x[i];		vel.y = d_Vel.y[i];		vel.z = d_Vel.z[i];
	treal3 Acc;
	Acc.x = Acc.y = Acc.z = 0.0f;//加速度の一時計算
	real invn0 = 1.0f / n0;

	int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
	int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
	int iz = (int)((pos.z - MINc.z) * DBinv) + 1;

	//流体ループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						if (d_Typ[j] == FLD) {
							treal3 posj;
							posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
							real dist2;// i,j距離の2乗
							dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
							real dist;//粒子の距離(絶対値)
							dist = sqrt(dist2);
							if (dist <= r) {
								/*real w = WEI(dist, r);//i,jの重み
								Acc.x += (d_Vel.x[j] - vel.x) * w;
								Acc.y += (d_Vel.y[j] - vel.y) * w;
								Acc.z += (d_Vel.z[j] - vel.z) * w;*/
								Acc.x += 2.0f * (d_Vel.x[j] - vel.x) * PCL_DST / dist / dist / dist;
								Acc.y += 2.0f * (d_Vel.y[j] - vel.y) * PCL_DST / dist / dist / dist;
								Acc.z += 2.0f * (d_Vel.z[j] - vel.z) * PCL_DST / dist / dist / dist;
							}

						}
					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}
	//ミラーループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfstM[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					//内外判定
					treal3 posMj;
					posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];


					treal3 p;//i,jの距離の成分
					p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
					real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
					real dist = sqrt(dist2);//粒子の距離(絶対値)
					if (dist <= r) {
						/*real w = WEI(dist, r);//i,jの重み
						Acc.x += (d_VelM.x[j] - vel.x) * w;
						Acc.y += (d_VelM.y[j] - vel.y) * w;
						Acc.z += (d_VelM.z[j] - vel.z) * w;*/
						Acc.x += 2.0f * (d_VelM.x[j] - vel.x) * PCL_DST / dist / dist / dist;
						Acc.y += 2.0f * (d_VelM.y[j] - vel.y) * PCL_DST / dist / dist / dist;
						Acc.z += 2.0f * (d_VelM.z[j] - vel.z) * PCL_DST / dist / dist / dist;
					}

					j = d_nxtM[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

	real n0KNM = invn0 * KNM_VSC;
	d_Acc.x[i] = n0KNM * Acc.x + G.x;
	d_Acc.y[i] = n0KNM * Acc.y + G.y;
	d_Acc.z[i] = n0KNM * Acc.z + G.z;

}


void DEMPS::VscTrm() {//粘性項・外力項の計算
	//printf_s("VscTrm start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_VscTrm << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_WLLVec, d_TypM, d_PosM, d_VelM, r, PCL_DST, n0, KNM_VSC, Vsc_coef, G, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("VscTrm finished!\n\n");
}


__device__ void d_ChkPcl(const int i, char* d_Typ, areal3 d_Pos, areal3 d_Vel, const  areal3 d_Acc, real* d_Prs, const  treal3 MINc, const  treal3 MAXc, const real PCL_DST, const real ulmax)//&要る？
{
	//if (d_Typ[i] != FLD) { return; }//呼び出し側でこの処理してる　壁動かすときは変更するかも

	if (d_Pos.x[i] < (MAXc.x - 3.1f * PCL_DST) && d_Pos.x[i] > (MINc.x + 3.1f * PCL_DST) &&
		d_Pos.y[i] < (MAXc.y - 3.1f * PCL_DST) && d_Pos.y[i] > (MINc.y + 3.1f * PCL_DST) &&
		d_Pos.z[i] < (MAXc.z - 3.1f * PCL_DST) && d_Pos.z[i] > (MINc.z + 3.1f * PCL_DST)) {//最大速度に制限

		treal3 Utmp;
		Utmp.x = d_Vel.x[i];	Utmp.y = d_Vel.y[i];	Utmp.z = d_Vel.z[i];
		real U = Utmp.x * Utmp.x + Utmp.y * Utmp.y + Utmp.z * Utmp.z;
		U = sqrt(U);
		if (U > ulmax) {
			Utmp.x *= ulmax / U;	Utmp.y *= ulmax / U;	Utmp.z *= ulmax / U;
			d_Vel.x[i] = Utmp.x;		d_Vel.y[i] = Utmp.y;		d_Vel.z[i] = Utmp.z;
		}

	}
	else {//発散している粒子は削除
		d_Typ[i] = GST;
		d_Pos.x[i] = d_Pos.y[i] = d_Pos.z[i] = 0.0f;
		d_Vel.x[i] = d_Vel.y[i] = d_Vel.z[i] = 0.0f;
		d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;
		d_Prs[i] = 0.0f;
	}

}


__global__ void d_UpPcl1(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, real* d_Prs, const  treal3 MINc, const  treal3 MAXc, const  real dt, const real  PCL_DST, const real ulmax)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	if (d_Typ[i] != FLD) { return; }

	d_Vel.x[i] += d_Acc.x[i] * dt;	d_Vel.y[i] += d_Acc.y[i] * dt;	d_Vel.z[i] += d_Acc.z[i] * dt;
	d_Pos.x[i] += d_Vel.x[i] * dt;	d_Pos.y[i] += d_Vel.y[i] * dt;	d_Pos.z[i] += d_Vel.z[i] * dt;
	d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;

	d_ChkPcl(i, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, MINc, MAXc, PCL_DST, ulmax);

}


void DEMPS::UpPcl1() {//仮の粒子移動
	//printf_s("UpPcl1 start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_UpPcl1 << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, MINc, MAXc, dt, PCL_DST, ulmax);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("UpPcl1 finished!\n\n");
}


__global__ void d_ChkCol(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, areal3 d_VelM, real* d_Dns, const real PCL_DST, const real r, const real rlim2, const real COL,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	char Typ = d_Typ[i];
	if (Typ != FLD) { return; }

	treal3 pos;		treal3 vec;	treal3 vec2;
	pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
	vec.x = d_Vel.x[i];		vec.y = d_Vel.y[i];		vec.z = d_Vel.z[i];
	vec2.x = d_Vel.x[i];		vec2.y = d_Vel.y[i];		vec2.z = d_Vel.z[i];
	real mi = d_Dns[Typ];
	real rlim_wall = 0.45f * PCL_DST;
	real rlim_wall2 = rlim_wall * rlim_wall;

	int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
	int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
	int iz = (int)((pos.z - MINc.z) * DBinv) + 1;

	//流体ループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						char typ = d_Typ[j];
						if (typ == FLD) {
							//内外判定
							treal3 posj;
							posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							if (dist2 < rlim2) {
								real fDT = (vec.x - d_Vel.x[j]) * p.x + (vec.y - d_Vel.y[j]) * p.y + (vec.z - d_Vel.z[j]) * p.z;
								if (fDT > 0.0f) {
									real mj = d_Dns[typ];
									fDT *= COL * mj / (mi + mj) / dist2;
									vec2.x -= p.x * fDT;	vec2.y -= p.y * fDT;	vec2.z -= p.z * fDT;
								}
							}

						}
						else if ((typ == WLL) || (typ == OBJ) || (typ == OBJ2)) {

							treal3 posj;
							posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];
							treal3 p;//i,jの距離の成分
							p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							if (dist2 < rlim_wall2) {
								real fDT = (vec.x - d_Vel.x[j]) * p.x + (vec.y - d_Vel.y[j]) * p.y + (vec.z - d_Vel.z[j]) * p.z;
								if (fDT > 0.0f) {
									real mj = d_Dns[typ];
									fDT *= COL * mj / (mi + mj) / dist2;
									vec2.x -= p.x * fDT;	vec2.y -= p.y * fDT;	vec2.z -= p.z * fDT;
								}
							}
						}

					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}
	//ミラーループ
	real mj = d_Dns[MRR];
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfstM[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					//内外判定
					treal3 posMj;
					posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];

					treal3 p;//i,jの距離の成分
					p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
					real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
					if (dist2 < rlim2) {
						real fDT = (vec.x - d_VelM.x[j]) * p.x + (vec.y - d_VelM.y[j]) * p.y + (vec.z - d_VelM.z[j]) * p.z;
						if (fDT > 0.0f) {
							fDT *= COL * mj / (mi + mj) / dist2;
							vec2.x -= p.x * fDT;	vec2.y -= p.y * fDT;	vec2.z -= p.z * fDT;
						}
					}

					j = d_nxtM[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}
	d_Acc.x[i] = vec2.x;	d_Acc.y[i] = vec2.y;	d_Acc.z[i] = vec2.z;

}


void __global__ d_UpChkCol(const int nP, char* d_Typ, areal3 d_Vel, areal3 d_Acc) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	if (d_Typ[i] != FLD) { return; }
	d_Vel.x[i] = d_Acc.x[i];		d_Vel.y[i] = d_Acc.y[i];		d_Vel.z[i] = d_Acc.z[i];
}

void DEMPS::ChkCol() {//仮の粒子移動
	//printf_s("ChkCol start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_ChkCol << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_WLLVec, d_TypM, d_PosM, d_VelM, d_Dns, PCL_DST, r, rlim2, COL, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	d_UpChkCol << <blocks_nP, threads >> > (nP, d_Typ, d_Vel, d_Acc);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("ChkCol finished!\n\n");
}


__global__ void d_MkPrs(const int nP, char* d_Typ, areal3 d_Pos, real* d_Dns, real* d_Prs, areal3 d_WLLVec, char* d_TypM, areal3 d_PosM, const real rp, const real n0_grad, const real Pmax, const real Prs_coef,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	char Typ = d_Typ[i];
	if (Typ != FLD) { return; }

	treal3 pos;
	pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
	real ni = 0.0f;

	int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
	int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
	int iz = (int)((pos.z - MINc.z) * DBinv) + 1;



	//流体ループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						if (d_Typ[j] == FLD) {
							//内外判定
							treal3 posj;
							posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							real dist = sqrt(dist2);;
							if (dist < rp) {
								ni += WEI_grad(dist, rp);
							}

						}
					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

	//ミラーループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfstM[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					//内外判定
					treal3 posMj;
					posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];

					treal3 p;//i,jの距離の成分
					p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
					real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
					real dist = sqrt(dist2);
					if (dist < rp) {
						ni += WEI_grad(dist, rp);
					}

					j = d_nxtM[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

	real mi = d_Dns[Typ];
	real pressure = (ni > n0_grad) * (ni - n0_grad) * Prs_coef * mi;

	if (pressure > Pmax)
	{
		pressure = Pmax;//最大値抑制
	}
	else if (pressure < 0.0f)
	{
		pressure = 0.0f;//負圧を0に
	}

	d_Prs[i] = pressure;

}


void DEMPS::MkPrs() {//仮の粒子移動
	//printf_s("MkPrs start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_MkPrs << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Dns, d_Prs, d_WLLVec, d_TypM, d_PosM, rp, n0_grad, Pmax, Prs_coef, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("MkPrs finished!\n\n");
}


#if 0
__global__ void d_PrsGrdTrm(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Acc, real* d_Prs, areal3 d_WLLVec, real* d_Dns, char* d_TypM, areal3 d_PosM, real* d_PrsM, const real rp, const real rp2, const real n0_grad,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nP) {
		char Typ = d_Typ[i];
		if (Typ == FLD) {
			real invn0 = 1.0f / n0_grad;
			real invro = -3.0f / d_Dns[Typ];//-DIM/ρ係数
			treal3 pos;
			pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
			real Prs_min = d_Prs[i];
			/*real invLeft[9] = {0.0f};//圧力勾配項左ブロックのインバースする前の3*3テンソル
			real Left[9] = { 0.0f };//圧力勾配左ブロックの3*3テンソル
			real Right[3] = { 0.0f };//圧力勾配右ブロックの圧力かけてるベクトル*/

			real invLeft0 = 0.0f;
			real invLeft1 = 0.0f;
			real invLeft2 = 0.0f;
			real invLeft3 = 0.0f;
			real invLeft4 = 0.0f;
			real invLeft5 = 0.0f;
			real invLeft6 = 0.0f;
			real invLeft7 = 0.0f;
			real invLeft8 = 0.0f;

			real Left0 = 0.0f;
			real Left1 = 0.0f;
			real Left2 = 0.0f;
			real Left3 = 0.0f;
			real Left4 = 0.0f;
			real Left5 = 0.0f;
			real Left6 = 0.0f;
			real Left7 = 0.0f;
			real Left8 = 0.0f;

			real Right0 = 0.0f;
			real Right1 = 0.0f;
			real Right2 = 0.0f;

			int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
			int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
			int iz = (int)((pos.z - MINc.z) * DBinv) + 1;


			//近傍最小圧力抽出
			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
									if (dist2 < rp2) {
										real prs = d_Prs[j];
										if (prs < Prs_min) {
											Prs_min = prs;
										}
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}


			//圧力項計算  入部・仲座モデル
			//流体ループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfst[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								if (d_Typ[j] == FLD) {
									//内外判定
									treal3 posj;
									posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

									treal3 p;//i,jの距離の成分
									p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
									real dist2;// i,j距離の2乗
									dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
									if (dist2 < rp2) {
										real dist = sqrt(dist2);
										real w = WEI_grad(dist, rp) / dist2;
										real Prsj_min = d_Prs[j] - Prs_min;
										//real Prsj_min = d_Prs[j] + d_Prs[i];
										invLeft0 += w * p.x * p.x;
										invLeft1 += w * p.x * p.y;
										invLeft2 += w * p.x * p.z;
										invLeft3 += w * p.y * p.x;
										invLeft4 += w * p.y * p.y;
										invLeft5 += w * p.y * p.z;
										invLeft6 += w * p.z * p.x;
										invLeft7 += w * p.z * p.y;
										invLeft8 += w * p.z * p.z;
										Right0 += w * Prsj_min * p.x;
										Right1 += w * Prsj_min * p.y;
										Right2 += w * Prsj_min * p.z;
									}

								}
							}
							j = d_nxt[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}
			//ミラーループ
			for (int jz = iz - 1; jz <= iz + 1; jz++) {
				for (int jy = iy - 1; jy <= iy + 1; jy++) {
					for (int jx = ix - 1; jx <= ix + 1; jx++) {
						int jb = jz * nBxy + jy * nBx + jx;
						int j = d_bfstM[jb];
						if (j == -1) continue;
						for (;;) {//粒子iの近傍粒子jのループ開始
							if (j != i) {
								//内外判定
								treal3 posMj;
								posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];


								treal3 p;//i,jの距離の成分
								p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
								real dist2;// i,j距離の2乗
								dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
								if (dist2 < rp2) {
									real dist = sqrt(dist2);
									real w = WEI_grad(dist, rp) / dist2;
									real PrsMj_min = d_PrsM[j] - Prs_min;
									//real PrsMj_min = d_PrsM[j] + d_Prs[i];
									invLeft0 += w * p.x * p.x;
									invLeft1 += w * p.x * p.y;
									invLeft2 += w * p.x * p.z;
									invLeft3 += w * p.y * p.x;
									invLeft4 += w * p.y * p.y;
									invLeft5 += w * p.y * p.z;
									invLeft6 += w * p.z * p.x;
									invLeft7 += w * p.z * p.y;
									invLeft8 += w * p.z * p.z;
									Right0 += w * PrsMj_min * p.x;
									Right1 += w * PrsMj_min * p.y;
									Right2 += w * PrsMj_min * p.z;
								}

							}
							j = d_nxtM[j];
							if (j == -1) break;
						}//粒子iの近傍粒子jのループ終了
					}
				}
			}


			invLeft0 = invn0 * invLeft0;//n0で割る
			invLeft1 = invn0 * invLeft1;
			invLeft2 = invn0 * invLeft2;
			invLeft3 = invn0 * invLeft3;
			invLeft4 = invn0 * invLeft4;
			invLeft5 = invn0 * invLeft5;
			invLeft6 = invn0 * invLeft6;
			invLeft7 = invn0 * invLeft7;
			invLeft8 = invn0 * invLeft8;

			//if ((invLeft0 == 0.0f) && (invLeft1 == 0.0f) && (invLeft2 == 0.0f) && (invLeft3 == 0.0f) && (invLeft4 == 0.0f) && (invLeft5 == 0.0f) && (invLeft6 == 0.0f) && (invLeft7 == 0.0f) && (invLeft8 == 0.0f)) { return; }//ディターミナント計算で0割しないように

			real DET_Left = 1.0f / ((invLeft0 * invLeft4 * invLeft8 + invLeft1 * invLeft5 * invLeft6 + invLeft2 * invLeft3 * invLeft7) - (invLeft2 * invLeft4 * invLeft6 + invLeft0 * invLeft5 * invLeft7 + invLeft1 * invLeft3 * invLeft8));//ここからインバース計算
			if (DET_Left == 0.0f) { return; }//ディターミナント計算で0割しないように
			Left0 = DET_Left * (invLeft4 * invLeft8 - invLeft5 * invLeft7);
			Left1 = DET_Left * -(invLeft1 * invLeft8 - invLeft2 * invLeft7);
			Left2 = DET_Left * (invLeft1 * invLeft5 - invLeft2 * invLeft4);
			Left3 = DET_Left * -(invLeft3 * invLeft8 - invLeft5 * invLeft6);
			Left4 = DET_Left * (invLeft0 * invLeft8 - invLeft2 * invLeft6);
			Left5 = DET_Left * -(invLeft0 * invLeft5 - invLeft2 * invLeft3);
			Left6 = DET_Left * (invLeft3 * invLeft7 - invLeft4 * invLeft6);
			Left7 = DET_Left * -(invLeft0 * invLeft7 - invLeft1 * invLeft6);
			Left8 = DET_Left * (invLeft0 * invLeft4 - invLeft1 * invLeft3);

			Right0 = invn0 * Right0;//n0で割る
			Right1 = invn0 * Right1;
			Right2 = invn0 * Right2;



			d_Acc.x[i] = invro * (Left0 * Right0 + Left1 * Right1 + Left2 * Right2);
			d_Acc.y[i] = invro * (Left3 * Right0 + Left4 * Right1 + Left5 * Right2);
			d_Acc.z[i] = invro * (Left6 * Right0 + Left7 * Right1 + Left8 * Right2);


		}
	}
}

#else
__global__ void d_PrsGrdTrm(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Acc, real* d_Prs, areal3 d_WLLVec, real* d_Dns, char* d_TypM, areal3 d_PosM, real* d_PrsM, const real rp, const real rp2, const real n0_grad,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_nxt, const int* d_bfstM, const int* d_nxtM)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	char Typ = d_Typ[i];
	if (Typ != FLD) { return; }

	//real invn0 = 1.0f / n0_grad;
	real invro = -3.0f / d_Dns[Typ];//-DIM/ρ係数
	treal3 pos;
	pos.x = d_Pos.x[i];	pos.y = d_Pos.y[i];	pos.z = d_Pos.z[i];
	real Prs_min = d_Prs[i];

	treal3 Acc = { 0.0f };
	real A3 = 3.0f / n0_grad;//Dimension / n0

	int ix = (int)((pos.x - MINc.x) * DBinv) + 1;
	int iy = (int)((pos.y - MINc.y) * DBinv) + 1;
	int iz = (int)((pos.z - MINc.z) * DBinv) + 1;



	//近傍最小圧力抽出
	//流体ループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						if (d_Typ[j] == FLD) {
							//内外判定
							treal3 posj;
							posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
							real dist2 = p.x * p.x + p.y * p.y + p.z * p.z;// i,j距離の2乗
							if (dist2 < rp2) {
								real prs = d_Prs[j];
								if (prs < Prs_min) {
									Prs_min = prs;
								}
							}

						}
					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

	//流体ループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						if (d_Typ[j] == FLD) {

							treal3 posj;
							posj.x = d_Pos.x[j];	posj.y = d_Pos.y[j];	posj.z = d_Pos.z[j];

							treal3 p;//i,jの距離の成分
							p.x = posj.x - pos.x;	p.y = posj.y - pos.y;	p.z = posj.z - pos.z;
							real dist2;// i,j距離の2乗
							dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
							if (dist2 < rp2) {
								real dist = sqrt(dist2);
								real w = WEI_grad(dist, rp) / dist2;
								real Prsj_min = d_Prs[j] - Prs_min;
								Acc.x += Prsj_min * w * p.x;
								Acc.y += Prsj_min * w * p.y;
								Acc.z += Prsj_min * w * p.z;
								/*real Prs = d_Prs[j] +d_Prs[i];
								Acc.x += Prs * w * p.x;
								Acc.y += Prs * w * p.y;
								Acc.z += Prs * w * p.z;*/
							}

						}
					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}
	//ミラーループ
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfstM[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						//内外判定
						treal3 posMj;
						posMj.x = d_PosM.x[j];	posMj.y = d_PosM.y[j];	posMj.z = d_PosM.z[j];

						treal3 p;//i,jの距離の成分
						p.x = posMj.x - pos.x;	p.y = posMj.y - pos.y;	p.z = posMj.z - pos.z;
						real dist2;// i,j距離の2乗
						dist2 = p.x * p.x + p.y * p.y + p.z * p.z;
						if (dist2 < rp2) {
							real dist = sqrt(dist2);
							real w = WEI_grad(dist, rp) / dist2;
							real PrsMj_min = d_PrsM[j] - Prs_min;
							Acc.x += PrsMj_min * w * p.x;
							Acc.y += PrsMj_min * w * p.y;
							Acc.z += PrsMj_min * w * p.z;
							/*real Prs = d_PrsM[j] + d_Prs[i];
							Acc.x += Prs * w * p.x;
							Acc.y += Prs * w * p.y;
							Acc.z += Prs * w * p.z;*/
						}

					}
					j = d_nxtM[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

	d_Acc.x[i] = invro * Acc.x * A3;
	d_Acc.y[i] = invro * Acc.y * A3;
	d_Acc.z[i] = invro * Acc.z * A3;

}
#endif


void DEMPS::PrsGrdTrm() {//仮の粒子移動

	//printf_s("PrsGrdTrm start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_PrsGrdTrm << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Acc, d_Prs, d_WLLVec, d_Dns, d_TypM, d_PosM, d_PrsM, rp, rp2, n0_grad, MINc, DBinv, nBx, nBxy, d_bfst, d_nxt, d_bfstM, d_nxtM);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("PrsGrdTrm finished!\n\n");
}


__global__ void d_UpPcl2(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, real* d_Prs, real* d_pav, const  treal3 MINc, const  treal3 MAXc, const  real dt, const real  PCL_DST, const real ulmax)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	if (d_Typ[i] != FLD) { return; }

	treal3 acc;
	acc.x = d_Acc.x[i];	acc.y = d_Acc.y[i];	acc.z = d_Acc.z[i];
	d_Vel.x[i] += acc.x * dt;	d_Vel.y[i] += acc.y * dt;	d_Vel.z[i] += acc.z * dt;
	d_Pos.x[i] += acc.x * dt * dt;	d_Pos.y[i] += acc.y * dt * dt;	d_Pos.z[i] += acc.z * dt * dt;
	d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;

	d_pav[i] += d_Prs[i];//時間平均圧力加算

	d_ChkPcl(i, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, MINc, MAXc, PCL_DST, ulmax);

}

void DEMPS::UpPcl2() {//圧力修正分の粒子移動
	//printf_s("UpPcl2 start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_UpPcl2 << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc, d_Prs, d_pav, MINc, MAXc, dt, PCL_DST, ulmax);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("UpPcl2 finished!\n\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////MPS_Function/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif


#if DEM_flg
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////DEM_Function/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__device__ int sign(real a, real b) {//aの絶対値にbの符号をつける
	int c = (b >= 0.0f) - (b <= 0.0f);
	return abs(a) * c;
}

__device__ real PairNumber(int* pair, int i, int j) {//粒子iの接触履歴照合
	int contact;
	int count = 0;

	for (int p = 0; p < NCP; p++) {//接触履歴あり　場所特定
		if (pair[p + i * NCP] == j) {
			contact = p;
			break;
		}
		else { count += 1; }
	}

	if (count == NCP) {//接触履歴なし　新規登録
		for (int q = 0; q < NCP; q++) {
			if (pair[q + i * NCP] == -2) {
				pair[q + i * NCP] = j;
				contact = q;
				break;
			}
		}
	}
	return contact;
}


__global__ void d_ColForce(const int nP, const int nPSLD, const int nPWLL, const real* d_D, const char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Ftotal, areal3 d_Omega, areal3 d_Torque,
	areal3 d_ep, int* d_pair, const real m, const real eta_n, const real eta_t, real kn, real kt, const real mu, const real dt,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	if (d_Typ[i] != SLD) { return; }

	real Di = d_D[i];
	treal3 Posi;	Posi.x = d_Pos.x[i];	Posi.y = d_Pos.y[i];	Posi.z = d_Pos.z[i];
	treal3 Veli;		Veli.x = d_Vel.x[i];		Veli.y = d_Vel.y[i];		Veli.z = d_Vel.z[i];
	treal3 Omegai;	Omegai.x = d_Omega.x[i];	Omegai.y = d_Omega.y[i];	Omegai.z = d_Omega.z[i];

	int ix = (int)((Posi.x - MINc.x) * DBinv) + 1;
	int iy = (int)((Posi.y - MINc.y) * DBinv) + 1;
	int iz = (int)((Posi.z - MINc.z) * DBinv) + 1;
	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						if (d_Typ[j] != FLD) {
							treal3 r_delta_i;
							treal3 r_delta_j;
							treal3 Xi;

							treal3 dist;
							real L;
							real L_2;
							treal3 n;

							treal3 Fcol;
							real Tr = 0.0f;
							treal3 dp;

							real Dj = d_D[j];
							treal3 Posj;	Posj.x = d_Pos.x[j];	Posj.y = d_Pos.y[j];	Posj.z = d_Pos.z[j];
							treal3 Velj;		Velj.x = d_Vel.x[j];		Velj.y = d_Vel.y[j];		Velj.z = d_Vel.z[j];
							treal3 Omegaj;	Omegaj.x = d_Omega.x[j];	Omegaj.y = d_Omega.y[j];	Omegaj.z = d_Omega.z[j];

							dist.x = Posj.x - Posi.x;   dist.y = Posj.y - Posi.y;   dist.z = Posj.z - Posi.z;

							L = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);

							n.x = dist.x / L;//l   
							n.y = dist.y / L;//m
							n.z = dist.z / L;//n

							L_2 = sqrt(n.x * n.x + n.y * n.y);//方向余弦の中身

							if (L - 0.5f * (Di + Dj) < 0.0f) {//接触判定

								int couple = PairNumber(d_pair, i, j) + i * NCP; //粒子i,jのpair配列番号 (i - nPWLL)から変更　前は壁粒子が先頭だったため
								treal3 ep;	ep.x = d_ep.x[couple];		ep.y = d_ep.y[couple];		ep.z = d_ep.z[couple];

								if (L_2 == 0) {//ジンバルロック対策

									r_delta_i.x = -Omegai.z * dt;
									r_delta_i.y = Omegai.y * dt;
									r_delta_i.z = Omegai.x * dt;
									r_delta_j.x = -Omegaj.z * dt;
									r_delta_j.y = Omegaj.y * dt;
									r_delta_j.z = Omegaj.x * dt;//角変位

									Xi.x = -(Veli.z - Velj.z) * dt;
									Xi.y = (Veli.y - Velj.y) * dt + (r_delta_i.z * Di + r_delta_j.z * Dj) * 0.5f;
									Xi.z = (Veli.x - Velj.x) * dt - (r_delta_i.y * Di + r_delta_j.y * Dj) * 0.5f;

									//////////////////////ローカルx///////////////////////////////
									ep.x += kn * Xi.x;
									dp.x = eta_n * Xi.x / dt;
									if (ep.x < 0.0f) { ep.x = dp.x = 0.0f; }
									d_ep.x[couple] = ep.x;//バネ更新
									Fcol.x = ep.x + dp.x;
									//////////////////////ローカルx///////////////////////////////

									//////////////////////ローカルy///////////////////////////////
									ep.y += kt * Xi.y;
									dp.y = eta_t * Xi.y / dt;
									if (ep.x < 0.0f) { ep.y = dp.y = 0.0f; }
									if (abs(ep.y) > mu * ep.x) { ep.y = mu * sign(ep.x, ep.y);		dp.y = 0.0f; }
									d_ep.y[couple] = ep.y;
									Fcol.y = ep.y + dp.y;
									//////////////////////ローカルy///////////////////////////////

									//////////////////////ローカルz///////////////////////////////
									ep.z += kt * Xi.z;
									dp.z = eta_t * Xi.z / dt;
									if (ep.x < 0.0f) { ep.z = dp.z = 0.0f; }
									if (abs(ep.z) > mu * ep.x) { ep.z = mu * sign(ep.x, ep.z);		dp.z = 0.0f; }
									d_ep.z[couple] = ep.z;
									Fcol.z = ep.z + dp.z;
									//////////////////////ローカルz///////////////////////////////


									//////////////////////ワールド///////////////////////////////
									d_Ftotal.x[i] -= Fcol.z;
									d_Ftotal.y[i] -= Fcol.y;
									d_Ftotal.z[i] -= -Fcol.x;

									d_Torque.x[i] -= Fcol.y * Di * 0.5f;
									d_Torque.y[i] -= -Fcol.z * Di * 0.5f;
									d_Torque.z[i] -= -Tr;
									//////////////////////ワールド///////////////////////////////

								}
								else {//通常

									r_delta_i.x = (n.x * Omegai.x + n.y * Omegai.y + n.z * Omegai.z) * dt;
									r_delta_i.y = (-n.y * Omegai.x / L_2 + n.x * Omegai.y / L_2) * dt;
									r_delta_i.z = (-n.x * n.z * Omegai.x / L_2 - n.y * n.z * Omegai.y / L_2 + L_2 * Omegai.z) * dt;
									r_delta_j.x = (n.x * Omegaj.x + n.y * Omegaj.y + n.z * Omegaj.z) * dt;
									r_delta_j.y = (-n.y * Omegaj.x / L_2 + n.x * Omegaj.y / L_2) * dt;
									r_delta_j.z = (-n.x * n.z * Omegaj.x / L_2 - n.y * n.z * Omegaj.y / L_2 + L_2 * Omegaj.z) * dt;//角変位

									Xi.x = (n.x * (Veli.x - Velj.x) + n.y * (Veli.y - Velj.y) + n.z * (Veli.z - Velj.z)) * dt;
									Xi.y = (-n.y * (Veli.x - Velj.x) / L_2 + n.x * (Veli.y - Velj.y) / L_2) * dt + (r_delta_i.z * Di + r_delta_j.z * Dj) * 0.5f;
									Xi.z = (-n.x * n.z * (Veli.x - Velj.x) / L_2 - n.y * n.z * (Veli.y - Velj.y) / L_2 + (Veli.z - Velj.z) * L_2) * dt - (r_delta_i.y * Di + r_delta_j.y * Dj) * 0.5f;

									//////////////////////ローカルx///////////////////////////////
									ep.x += kn * Xi.x;
									dp.x = eta_n * Xi.x / dt;
									if (ep.x < 0.0f) { ep.x = dp.x = 0.0f; }
									d_ep.x[couple] = ep.x;//バネ更新
									Fcol.x = ep.x + dp.x;
									//////////////////////ローカルx///////////////////////////////

									//////////////////////ローカルy///////////////////////////////
									ep.y += kt * Xi.y;
									dp.y = eta_t * Xi.y / dt;
									if (ep.x < 0.0f) { ep.y = dp.y = 0.0f; }
									if (abs(ep.y) > mu * ep.x) { ep.y = mu * sign(ep.x, ep.y);		dp.y = 0.0f; }
									d_ep.y[couple] = ep.y;
									Fcol.y = ep.y + dp.y;
									//////////////////////ローカルy///////////////////////////////

									//////////////////////ローカルz///////////////////////////////
									ep.z += kt * Xi.z;
									dp.z = eta_t * Xi.z / dt;
									if (ep.x < 0.0f) { ep.z = dp.z = 0.0f; }
									if (abs(ep.z) > mu * ep.x) { ep.z = mu * sign(ep.x, ep.z);		dp.z = 0.0f; }
									d_ep.z[couple] = ep.z;
									Fcol.z = ep.z + dp.z;
									//////////////////////ローカルz///////////////////////////////


									//////////////////////ワールド///////////////////////////////
									d_Ftotal.x[i] -= n.x * Fcol.x - n.y * Fcol.y / L_2 - n.x * n.z * Fcol.z / L_2;
									d_Ftotal.y[i] -= n.y * Fcol.x + n.x * Fcol.y / L_2 - n.y * n.z * Fcol.z / L_2;
									d_Ftotal.z[i] -= n.z * Fcol.x + Fcol.z * L_2;

									d_Torque.x[i] -= n.x * Tr - (-n.y * Fcol.z / L_2 + n.x * n.z * Fcol.y / L_2) * Di * 0.5f;
									d_Torque.y[i] -= n.y * Tr - (n.x * Fcol.z / L_2 + n.y * n.z * Fcol.y / L_2) * Di * 0.5f;
									d_Torque.z[i] -= n.z * Tr + Fcol.y * L_2 * Di * 0.5f;
									//////////////////////ワールド///////////////////////////////

								}
							}
							else {//接触していない時,epを0,pairを-1にしておく。
								for (int k = 0; k < NCP; k++) {
									int kiNPWLLNCP = k + i * NCP;
									if (d_pair[kiNPWLLNCP] == j) {
										d_ep.x[kiNPWLLNCP] = 0.0f;
										d_ep.y[kiNPWLLNCP] = 0.0f;
										d_ep.z[kiNPWLLNCP] = 0.0f;
										d_pair[kiNPWLLNCP] = -2;
										break;
									}
								}
							}
						}
					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

}


void DEMPS::ColForce() {
	//printf_s("ColForce start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = (nP);	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_ColForce << <blocks_nP, threads >> > (nP, nPSLD, nPWLL, d_D, d_Typ, d_Pos, d_Vel, d_Ftotal, d_Omega, d_Torque, d_ep, d_pair, m, eta_n, eta_t, kn, kt, mu, dt, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("ColForce finished!\n\n");
}


__global__ void d_update(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Ftotal, areal3 d_Omega, areal3 d_Torque, const real usmax,
	const real m, const real dt, const treal3 MINc, const treal3 MAXc, const treal3 G, const real PCL_DST, const real I)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i >= nP) { return; }

	if (d_Typ[i] != SLD) { return; }//粉体速度更新

		/*treal3 Vtmp;
		Vtmp.x = d_Vel.x[i];		Vtmp.y = d_Vel.y[i];		Vtmp.z = d_Vel.z[i];*/

	d_Vel.x[i] += d_Ftotal.x[i] * dt / m;
	d_Vel.y[i] += d_Ftotal.y[i] * dt / m + G.y * dt;
	d_Vel.z[i] += d_Ftotal.z[i] * dt / m;

	treal3 Utmp;
	Utmp.x = d_Vel.x[i];	Utmp.y = d_Vel.y[i];	Utmp.z = d_Vel.z[i];
	real U = Utmp.x * Utmp.x + Utmp.y * Utmp.y + Utmp.z * Utmp.z;
	U = sqrt(U);
	if (U > usmax) {
		Utmp.x *= usmax / U;	Utmp.y *= usmax / U;	Utmp.z *= usmax / U;
		d_Vel.x[i] = Utmp.x;		d_Vel.y[i] = Utmp.y;		d_Vel.z[i] = Utmp.z;
	}

	/*d_Pos.x[i] += 0.5f * (Vtmp.x + d_Vel.x[i]) * dt;//前ステップと現在の速度の平均値分だけ移動させる
	d_Pos.y[i] += 0.5f * (Vtmp.y + d_Vel.y[i]) * dt;
	d_Pos.z[i] += 0.5f * (Vtmp.z + d_Vel.z[i]) * dt;*/

	d_Pos.x[i] += d_Vel.x[i] * dt;
	d_Pos.y[i] += d_Vel.y[i] * dt;
	d_Pos.z[i] += d_Vel.z[i] * dt;

	if (d_Pos.x[i] > MAXc.x - 3.0f * PCL_DST) { d_Typ[i] = GST; }
	else if (d_Pos.y[i] > MAXc.y - 3.0f * PCL_DST) { d_Typ[i] = GST; }
	else if (d_Pos.z[i] > MAXc.z - 3.0f * PCL_DST) { d_Typ[i] = GST; }
	else if (d_Pos.x[i] < MINc.x + 3.0f * PCL_DST) { d_Typ[i] = GST; }
	else if (d_Pos.y[i] < MINc.y + 3.0f * PCL_DST) { d_Typ[i] = GST; }
	else if (d_Pos.z[i] < MINc.z + 3.0f * PCL_DST) { d_Typ[i] = GST; }

	d_Omega.x[i] += d_Torque.x[i] * dt / I;
	d_Omega.y[i] += d_Torque.y[i] * dt / I;
	d_Omega.z[i] += d_Torque.z[i] * dt / I;

	d_Ftotal.x[i] = 0.0f;
	d_Ftotal.y[i] = 0.0f;
	d_Ftotal.z[i] = 0.0f;

	d_Torque.x[i] = 0.0f;
	d_Torque.y[i] = 0.0f;
	d_Torque.z[i] = 0.0f;

}


void DEMPS::update() {
	//printf_s("update start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_update << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Ftotal, d_Omega, d_Torque, usmax, m, dt, MINc, MAXc, G, PCL_DST, I);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("update finished!\n\n");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////DEM_Function/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif


#if Multi_flg
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Multi_Function/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
//座標返還
__global__ void d_SLD_FLD(const int nP, real* d_D, real* d_Dns, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Ftotal, areal3 d_Omega, areal3 d_Torque, const real dt,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i >= nP) { return; }
	char Typ = d_Typ[i];
	if (Typ != SLD) { return; }

	treal3 Posi;	Posi.x = d_Pos.x[i];	Posi.y = d_Pos.y[i];	Posi.z = d_Pos.z[i];
	treal3 Veli;		Veli.x = d_Vel.x[i];		Veli.y = d_Vel.y[i];		Veli.z = d_Vel.z[i];
	treal3 Omegai;	Omegai.x = d_Omega.x[i];	Omegai.y = d_Omega.y[i];	Omegai.z = d_Omega.z[i];
	real Di = d_D[i];
	real dnsi = d_Dns[Typ];


	int ix = (int)((Posi.x - MINc.x) * DBinv) + 1;
	int iy = (int)((Posi.y - MINc.y) * DBinv) + 1;
	int iz = (int)((Posi.z - MINc.z) * DBinv) + 1;

	for (int jz = iz - 1; jz <= iz + 1; jz++) {
		for (int jy = iy - 1; jy <= iy + 1; jy++) {
			for (int jx = ix - 1; jx <= ix + 1; jx++) {
				int jb = jz * nBxy + jy * nBx + jx;
				int j = d_bfst[jb];
				if (j == -1) continue;
				for (;;) {//粒子iの近傍粒子jのループ開始
					if (j != i) {
						char typ = d_Typ[j];
						if (typ == FLD) {
							treal3 vs = { 0.0f };//ローカル固体速度
							treal3 vl = { 0.0f };//ローカル流体速度

							treal3 r_delta_i;//ローカル回転速度 Tgl * Omega

							treal3 dist;
							real L;
							real L_2;
							treal3 n;

							treal3 Posj;	Posj.x = d_Pos.x[j];	Posj.y = d_Pos.y[j];	Posj.z = d_Pos.z[j];
							treal3 Velj;		Velj.x = d_Vel.x[j];		Velj.y = d_Vel.y[j];		Velj.z = d_Vel.z[j];
							real Dj = d_D[j];
							real dnsj = d_Dns[typ];

							dist.x = Posj.x - Posi.x;   dist.y = Posj.y - Posi.y;   dist.z = Posj.z - Posi.z;

							L = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);

							n.x = dist.x / L;//l   
							n.y = dist.y / L;//m
							n.z = dist.z / L;//n

							L_2 = sqrt(n.x * n.x + n.y * n.y);//方向余弦の中身

							treal3 Fpint = { 0.0f };//液体=>固体の力
							real area = 0.5f * (Di + Dj);
							if (L < area) {//接触判定
								real w = WEI(L, area) / wtotal;
								if (L_2 == 0) {//ジンバルロック対策
									r_delta_i.x = -Omegai.z;
									r_delta_i.y = Omegai.y;
									r_delta_i.z = Omegai.x;

									vs.x = -Veli.z;
									vs.y = Veli.y + r_delta_i.z * L;
									vs.z = Veli.x - r_delta_i.y * L;
									vl.x = -Velj.z;
									vl.y = Velj.y;
									vl.z = Velj.x;

									Fpint.x = w * (dnsi * vs.x - dnsj * vl.x) * dt;
									Fpint.y = w * (dnsi * vs.y - dnsj * vl.y) * dt;
									Fpint.z = w * (dnsi * vs.z - dnsj * vl.z) * dt;

									d_Ftotal.x[i] -= Fpint.z;
									d_Ftotal.y[i] -= Fpint.y;
									d_Ftotal.z[i] -= -Fpint.x;

									d_Torque.x[i] -= Fpint.y * L;
									d_Torque.y[i] -= -Fpint.z * L;

									//atomicAdd(&d_Ftotal.x[j], Fpint.z);
									//atomicAdd(&d_Ftotal.y[j], Fpint.y);
									//atomicAdd(&d_Ftotal.z[j], -Fpint.x);
									d_Ftotal.x[j] += Fpint.z;
									d_Ftotal.y[j] += Fpint.y;
									d_Ftotal.z[j] += -Fpint.x;

								}
								else {//通常
									r_delta_i.x = (n.x * Omegai.x + n.y * Omegai.y + n.z * Omegai.z);
									r_delta_i.y = (-n.y * Omegai.x / L_2 + n.x * Omegai.y / L_2);
									r_delta_i.z = (-n.x * n.z * Omegai.x / L_2 - n.y * n.z * Omegai.y / L_2 + L_2 * Omegai.z);

									vs.x = (n.x * Veli.x + n.y * Veli.y + n.z * Veli.z);
									vs.y = (-n.y * Veli.x / L_2 + n.x * Veli.y / L_2) + r_delta_i.z * L;
									vs.z = (-n.x * n.z * Veli.x / L_2 - n.y * n.z * Veli.y / L_2 + Veli.z * L_2) - r_delta_i.y * L;
									vl.x = (n.x * Velj.x + n.y * Velj.y + n.z * Velj.z);
									vl.y = (-n.y * Velj.x / L_2 + n.x * Velj.y / L_2);
									vl.z = (-n.x * n.z * Velj.x / L_2 - n.y * n.z * Velj.y / L_2 + Velj.z * L_2);

									Fpint.x = w * (dnsi * vs.x - dnsj * vl.x) * dt;
									Fpint.y = w * (dnsi * vs.y - dnsj * vl.y) * dt;
									Fpint.z = w * (dnsi * vs.z - dnsj * vl.z) * dt;

									d_Ftotal.x[i] -= n.x * Fpint.x - n.y * Fpint.y / L_2 - n.x * n.z * Fpint.z / L_2;
									d_Ftotal.y[i] -= n.y * Fpint.x + n.x * Fpint.y / L_2 - n.y * n.z * Fpint.z / L_2;
									d_Ftotal.z[i] -= n.z * Fpint.x + Fpint.z * L_2;

									d_Torque.x[i] -= - (-n.y * Fpint.z / L_2 + n.x * n.z * Fpint.y / L_2) * L;
									d_Torque.y[i] -= - (n.x * Fpint.z / L_2 + n.y * n.z * Fpint.y / L_2) * L;
									d_Torque.z[i] -= Fpint.y * L_2 * L;

									//atomicAdd(&d_Ftotal.x[j], n.x * Fpint.x - n.y * Fpint.y / L_2 - n.x * n.z * Fpint.z / L_2);
									//atomicAdd(&d_Ftotal.y[j], n.y * Fpint.x + n.x * Fpint.y / L_2 - n.y * n.z * Fpint.z / L_2);
									//atomicAdd(&d_Ftotal.z[j], n.z * Fpint.x + Fpint.z * L_2);
									d_Ftotal.x[j] += n.x * Fpint.x - n.y * Fpint.y / L_2 - n.x * n.z * Fpint.z / L_2;
									d_Ftotal.y[j] += n.y * Fpint.x + n.x * Fpint.y / L_2 - n.y * n.z * Fpint.z / L_2;
									d_Ftotal.z[j] += -n.z * Fpint.x + Fpint.z * L_2;

								}
							}

						}

					}
					j = d_nxt[j];
					if (j == -1) break;
				}//粒子iの近傍粒子jのループ終了
			}
		}
	}

}


//回転込みの運動量交換　制作ちゅう
__global__ void d_SLD_FLD(const int nP, real* d_D, real* d_Dns, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, areal3 d_Ftotal, areal3 d_Omega, areal3 d_Torque, const real dt, const real Vol_SLD, const real Vol_FLD,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	char Typ = d_Typ[i];

	if (Typ == SLD) {
		treal3 Posi;	Posi.x = d_Pos.x[i];	Posi.y = d_Pos.y[i];	Posi.z = d_Pos.z[i];
		treal3 Veli;		Veli.x = d_Vel.x[i];		Veli.y = d_Vel.y[i];		Veli.z = d_Vel.z[i];
		treal3 vec2;		vec2.x = d_Vel.x[i];		vec2.y = d_Vel.y[i];		vec2.z = d_Vel.z[i];
		real Di = d_D[i];
		//real dnsi = d_Dns[Typ];
		real mi = d_Dns[Typ] * Vol_SLD;

		int ix = (int)((Posi.x - MINc.x) * DBinv) + 1;
		int iy = (int)((Posi.y - MINc.y) * DBinv) + 1;
		int iz = (int)((Posi.z - MINc.z) * DBinv) + 1;

		//流体ループ
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							char typ = d_Typ[j];
							if (typ == FLD) {
								treal3 Posj;	Posj.x = d_Pos.x[j];	Posj.y = d_Pos.y[j];	Posj.z = d_Pos.z[j];
								treal3 dist;	 dist.x = Posj.x - Posi.x;   dist.y = Posj.y - Posi.y;   dist.z = Posj.z - Posi.z;
								real L = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);
								real Dj = d_D[j];
								real area = 0.5f * (Di + Dj);
								if (L < area) {
									treal3 Velj;		Velj.x = d_Vel.x[j];		Velj.y = d_Vel.y[j];		Velj.z = d_Vel.z[j];
									treal3 vec;
									vec.x = Veli.x;
									vec.y = Veli.y;
									vec.z = Veli.z;
									//real dnsj = d_Dns[typ];
									real mj = d_Dns[typ] * Vol_FLD;
									real fDT = (vec.x - Velj.x) * dist.x + (vec.y - Velj.y) * dist.y + (vec.z - Velj.z) * dist.z;
									if (fDT > 0.0f) {//反発だけの場合の条件？
										fDT *= mj / (mi + mj) / L;
										vec2.x -= dist.x * fDT;	vec2.y -= dist.y * fDT;	vec2.z -= dist.z * fDT;
									}
								}

							}

						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
		d_Acc.x[i] = vec2.x;	d_Acc.y[i] = vec2.y;	d_Acc.z[i] = vec2.z;
	}

	else if (Typ == FLD) {
		treal3 Posi;	Posi.x = d_Pos.x[i];	Posi.y = d_Pos.y[i];	Posi.z = d_Pos.z[i];
		treal3 Veli;		Veli.x = d_Vel.x[i];		Veli.y = d_Vel.y[i];		Veli.z = d_Vel.z[i];
		treal3 vec2;		vec2.x = d_Vel.x[i];		vec2.y = d_Vel.y[i];		vec2.z = d_Vel.z[i];
		real Di = d_D[i];
		//real dnsi = d_Dns[Typ];
		real mi = d_Dns[Typ] * Vol_FLD;

		int ix = (int)((Posi.x - MINc.x) * DBinv) + 1;
		int iy = (int)((Posi.y - MINc.y) * DBinv) + 1;
		int iz = (int)((Posi.z - MINc.z) * DBinv) + 1;

		//流体ループ
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							char typ = d_Typ[j];
							if (typ == SLD) {
								treal3 Posj;	Posj.x = d_Pos.x[j];	Posj.y = d_Pos.y[j];	Posj.z = d_Pos.z[j];
								treal3 dist;	 dist.x = Posj.x - Posi.x;   dist.y = Posj.y - Posi.y;   dist.z = Posj.z - Posi.z;
								real L = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);
								real Dj = d_D[j];
								real area = 0.5f * (Di + Dj);
								if (L < area) {
									treal3 Velj;
									Velj.x = d_Vel.x[j];
									Velj.y = d_Vel.y[j]; 
									Velj.z = d_Vel.z[j]; 
									//real dnsj = d_Dns[typ];
									real mj = d_Dns[typ] * Vol_SLD;
									real fDT = (Veli.x - Velj.x) * dist.x + (Veli.y - Velj.y) * dist.y + (Veli.z - Velj.z) * dist.z;
									if (fDT > 0.0f) {//反発だけの場合の条件？
										fDT *= mj / (mi + mj) / L;
										vec2.x -= dist.x * fDT;	vec2.y -= dist.y * fDT;	vec2.z -= dist.z * fDT;
									}
								}

							}

						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
		d_Acc.x[i] = vec2.x;	d_Acc.y[i] = vec2.y;	d_Acc.z[i] = vec2.z;
	}

}


void DEMPS::SLD_FLD() {
	//printf_s("SLD_FLD start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = (nP);	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_SLD_FLD << <blocks_nP, threads >> > (nP, d_D, d_Dns, d_Typ, d_Pos, d_Vel, d_Acc, d_Ftotal, d_Omega, d_Torque, dt, Vol_SLD, Vol_FLD, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	CHECK(cudaDeviceSynchronize());

	//printf_s("SLD_FLD finished!\n\n");
}

#else
//chkcolと同じ感じで計算
__global__ void d_SLD_FLD(const int nP, real* d_D, real* d_Dns, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc, areal3 d_Ftotal, areal3 d_Omega, areal3 d_Torque, const real dt, const real Vol_SLD, const real Vol_FLD,
	const treal3 MINc, const real DBinv, const int nBx, const int nBxy, const int* d_bfst, const int* d_blst, const int* d_nxt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	char Typ = d_Typ[i];

	if (Typ == SLD) {
		treal3 Posi;	Posi.x = d_Pos.x[i];	Posi.y = d_Pos.y[i];	Posi.z = d_Pos.z[i];
		treal3 Veli;		Veli.x = d_Vel.x[i];		Veli.y = d_Vel.y[i];		Veli.z = d_Vel.z[i];
		treal3 vec2;		vec2.x = d_Vel.x[i];		vec2.y = d_Vel.y[i];		vec2.z = d_Vel.z[i];
		//treal3 Omegai;	Omegai.x = d_Omega.x[i];	Omegai.y = d_Omega.y[i];	Omegai.z = d_Omega.z[i];
		real Di = d_D[i];
		//real dnsi = d_Dns[Typ];
		real mi = d_Dns[Typ] * Vol_SLD;

		int ix = (int)((Posi.x - MINc.x) * DBinv) + 1;
		int iy = (int)((Posi.y - MINc.y) * DBinv) + 1;
		int iz = (int)((Posi.z - MINc.z) * DBinv) + 1;

		//流体ループ
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							char typ = d_Typ[j];
							if (typ == FLD) {
								treal3 Posj;	Posj.x = d_Pos.x[j];	Posj.y = d_Pos.y[j];	Posj.z = d_Pos.z[j];
								treal3 dist;	 dist.x = Posj.x - Posi.x;   dist.y = Posj.y - Posi.y;   dist.z = Posj.z - Posi.z;
								real L = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);
								real Dj = d_D[j];
								real area = 0.5f * (Di + Dj);
								if (L < area) {
									treal3 Velj;		Velj.x = d_Vel.x[j];		Velj.y = d_Vel.y[j];		Velj.z = d_Vel.z[j];
									treal3 vec;		
									vec.x = Veli.x;// +(dist.y * Omegai.z - dist.z * Omegai.y);
									vec.y = Veli.y;// +(dist.z * Omegai.x - dist.x * Omegai.z);
									vec.z = Veli.z;// +(dist.x * Omegai.y - dist.y * Omegai.x);
									//real dnsj = d_Dns[typ];
									real mj = d_Dns[typ] * Vol_FLD;
									real fDT = (vec.x - Velj.x) * dist.x + (vec.y - Velj.y) * dist.y + (vec.z - Velj.z) * dist.z;
									if (fDT > 0.0f) {//反発だけの場合の条件？
										fDT *= mj / (mi + mj) / L;
										vec2.x -= dist.x * fDT;	vec2.y -= dist.y * fDT;	vec2.z -= dist.z * fDT;
									}
								}

							}

						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
		d_Acc.x[i] = vec2.x;	d_Acc.y[i] = vec2.y;	d_Acc.z[i] = vec2.z;
	}

	else if (Typ == FLD) {
		treal3 Posi;	Posi.x = d_Pos.x[i];	Posi.y = d_Pos.y[i];	Posi.z = d_Pos.z[i];
		treal3 Veli;		Veli.x = d_Vel.x[i];		Veli.y = d_Vel.y[i];		Veli.z = d_Vel.z[i];
		treal3 vec2;		vec2.x = d_Vel.x[i];		vec2.y = d_Vel.y[i];		vec2.z = d_Vel.z[i];
		real Di = d_D[i];
		//real dnsi = d_Dns[Typ];
		real mi= d_Dns[Typ] * Vol_FLD;

		int ix = (int)((Posi.x - MINc.x) * DBinv) + 1;
		int iy = (int)((Posi.y - MINc.y) * DBinv) + 1;
		int iz = (int)((Posi.z - MINc.z) * DBinv) + 1;

		//流体ループ
		for (int jz = iz - 1; jz <= iz + 1; jz++) {
			for (int jy = iy - 1; jy <= iy + 1; jy++) {
				for (int jx = ix - 1; jx <= ix + 1; jx++) {
					int jb = jz * nBxy + jy * nBx + jx;
					int j = d_bfst[jb];
					if (j == -1) continue;
					for (;;) {//粒子iの近傍粒子jのループ開始
						if (j != i) {
							char typ = d_Typ[j];
							if (typ == SLD) {
								treal3 Posj;	Posj.x = d_Pos.x[j];	Posj.y = d_Pos.y[j];	Posj.z = d_Pos.z[j];
								treal3 dist;	 dist.x = Posj.x - Posi.x;   dist.y = Posj.y - Posi.y;   dist.z = Posj.z - Posi.z;
								real L = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z);
								real Dj = d_D[j];
								real area = 0.5f * (Di + Dj);
								if (L < area) {
									//treal3 Omegaj;	Omegaj.x = d_Omega.x[j];	Omegaj.y = d_Omega.y[j];	Omegaj.z = d_Omega.z[j];
									treal3 Velj;
									Velj.x = d_Vel.x[j]; //- (dist.y * Omegaj.z - dist.z * Omegaj.y);//distの符号が逆になってる？から
									Velj.y = d_Vel.y[j]; //-(dist.z * Omegaj.x - dist.x * Omegaj.z);
									Velj.z = d_Vel.z[j]; //-(dist.x * Omegaj.y - dist.y * Omegaj.x);
									//real dnsj = d_Dns[typ];
									real mj = d_Dns[typ] * Vol_SLD;
									real fDT = (Veli.x - Velj.x) * dist.x + (Veli.y - Velj.y) * dist.y + (Veli.z - Velj.z) * dist.z;
									if (fDT > 0.0f) {//反発だけの場合の条件？
									fDT *= mj / (mi + mj) / L;
									vec2.x -= dist.x * fDT;	vec2.y -= dist.y * fDT;	vec2.z -= dist.z * fDT;
									}
								}

							}

						}
						j = d_nxt[j];
						if (j == -1) break;
					}//粒子iの近傍粒子jのループ終了
				}
			}
		}
		d_Acc.x[i] = vec2.x;	d_Acc.y[i] = vec2.y;	d_Acc.z[i] = vec2.z;
	}

}
#endif

void DEMPS::SLD_FLD() {
	//printf_s("SLD_FLD start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = (nP);	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_SLD_FLD << <blocks_nP, threads >> > (nP, d_D, d_Dns, d_Typ, d_Pos, d_Vel, d_Acc, d_Ftotal, d_Omega, d_Torque, dt, Vol_SLD, Vol_FLD, MINc, DBinv, nBx, nBxy, d_bfst, d_blst, d_nxt);
	CHECK(cudaDeviceSynchronize());

	//printf_s("SLD_FLD finished!\n\n");
}

#if 0
__global__ void d_SF_update(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Ftotal, areal3 d_Omega, areal3 d_Torque, const real m, const real dt, const real I)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }

	char Typ = d_Typ[i];
	if (Typ == SLD) {
		treal3 acc;
		acc.x = d_Ftotal.x[i] / Dns_SLD;
		acc.y = d_Ftotal.y[i] / Dns_SLD;
		acc.z = d_Ftotal.z[i] / Dns_SLD;//質量 or 密度 ?

		d_Vel.x[i] += acc.x * dt;
		d_Vel.y[i] += acc.y * dt;
		d_Vel.z[i] += acc.z * dt;

		d_Pos.x[i] += acc.x * dt * dt;
		d_Pos.y[i] += acc.y * dt * dt;
		d_Pos.z[i] += acc.z * dt * dt;

		d_Omega.x[i] += d_Torque.x[i] * dt / I;
		d_Omega.y[i] += d_Torque.y[i] * dt / I;
		d_Omega.z[i] += d_Torque.z[i] * dt / I;

	}
	else if (Typ == FLD) {
		treal3 acc;
		acc.x = d_Ftotal.x[i] / Dns_FLD;
		acc.y = d_Ftotal.y[i] / Dns_FLD;
		acc.z = d_Ftotal.z[i] / Dns_FLD;

		d_Vel.x[i] += acc.x * dt;
		d_Vel.y[i] += acc.y * dt;
		d_Vel.z[i] += acc.z * dt;

		d_Pos.x[i] += acc.x * dt * dt;
		d_Pos.y[i] += acc.y * dt * dt;
		d_Pos.z[i] += acc.z * dt * dt;
	}

	d_Ftotal.x[i] = 0.0f;
	d_Ftotal.y[i] = 0.0f;
	d_Ftotal.z[i] = 0.0f;

	d_Torque.x[i] = 0.0f;
	d_Torque.y[i] = 0.0f;
	d_Torque.z[i] = 0.0f;

}


void DEMPS::SF_update() {
	//printf_s("update start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_SF_update << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Ftotal, d_Omega, d_Torque, m, dt, I);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("update finished!\n\n");
}
#else
__global__ void d_SF_update(const int nP, char* d_Typ, areal3 d_Pos, areal3 d_Vel, areal3 d_Acc)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nP) { return; }
	char Typ = d_Typ[i];
	if ((Typ == SLD) || (Typ == FLD)) {
		d_Vel.x[i] = d_Acc.x[i];		d_Vel.y[i] = d_Acc.y[i];		d_Vel.z[i] = d_Acc.z[i];
		d_Acc.x[i] = d_Acc.y[i] = d_Acc.z[i] = 0.0f;
	}
	
}


void DEMPS::SF_update() {
	//printf_s("update start!\n");
	////////////////cudaスレッド設定/////////////////////
	dim3 threads(THREADS, 1, 1);
	int TOTAL_THREADS = nBxyz;	int BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nBxyz(BLOCKS, 1, 1);
	TOTAL_THREADS = (nP);	BLOCKS = TOTAL_THREADS / THREADS + 1;
	dim3 blocks_nP(BLOCKS, 1, 1);
	//////////////////////////////////////////////
	d_SF_update << <blocks_nP, threads >> > (nP, d_Typ, d_Pos, d_Vel, d_Acc);
	//CHECK(cudaDeviceSynchronize());

	//printf_s("update finished!\n\n");
}
#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Multi_Function/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif


void DEMPS::DevicetoHost() {
	//printf_s("DevicetoHost start!\n");
	CHECK(cudaMemcpy(Typ, d_Typ, sizeof(char) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Pos.x, d_Pos.x, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Pos.y, d_Pos.y, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Pos.z, d_Pos.z, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Vel.x, d_Vel.x, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Vel.y, d_Vel.y, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Vel.z, d_Vel.z, sizeof(real) * nP, cudaMemcpyDeviceToHost));

#if 0 //updateで初期化してしまってる　drillのやつ参考に要変更
	CHECK(cudaMemcpy(Ftotal.x, d_Ftotal.x, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Ftotal.y, d_Ftotal.y, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(Ftotal.z, d_Ftotal.z, sizeof(real) * nP, cudaMemcpyDeviceToHost));
#endif


#if MPS_flg
	CHECK(cudaMemcpy(Prs, d_Prs, sizeof(real) * nP, cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(pav, d_pav, sizeof(real) * nP, cudaMemcpyDeviceToHost));

	CHECK(cudaMemcpy(TypM, d_TypM, sizeof(char) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PosM.x, d_PosM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PosM.y, d_PosM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PosM.z, d_PosM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(VelM.x, d_VelM.x, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(VelM.y, d_VelM.y, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(VelM.z, d_VelM.z, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
	CHECK(cudaMemcpy(PrsM, d_PrsM, sizeof(real) * (nP * NumMRR), cudaMemcpyDeviceToHost));
#endif


	//printf_s("DevicetoHost finished!\n\n");
}





void DEMPS::ClcDEMPS() {
	iF = 0;
	TIM = 0.0f;
	outtime = 0.0f;
	OPT_FQC = 1.0f;

#if MPS_flg
	MkBkt();
	Surface_Edge();//壁動かすときは毎ステップ実行できるように移動させること！
#endif

	WrtDatWLL();
	WrtDat();
	printf_s("\n");

	iF++;
	OPT_FQC = 0.0f;


	while (1) {

		//////////////////////出力///////////////////////////
		if (outtime >= output_time) {
			outtime -= output_time;

			DevicetoHost();

			WrtDat();
			printf_s("Time = %f\n\n", TIM);
			iF++;
			if (TIM >= FIN_TIM) {
				WrtDat2();
				break;
			}
			OPT_FQC = 0.0f;
		}
		//////////////////////出力///////////////////////////

#if MPS_flg
		//////////////////////MPS///////////////////////////
		MkBkt();

		ResetMRR();
		GenMRR_nonslip();
		MkBkt_MRR();

		VscTrm();

		UpPcl1();

		MkBkt();

		ChkCol();

		MkBkt();//いらん？

		MkPrs();

		ResetMRR();
		GenMRR_nonslip();
		MkBkt_MRR();

		PrsGrdTrm();

		UpPcl2();
		//////////////////////MPS///////////////////////////
#endif


#if Multi_flg
		//////////////////////相互作用///////////////////////////

		SLD_FLD();

		SF_update();

		//////////////////////相互作用///////////////////////////
#endif


#if DEM_flg
		//////////////////////DEM///////////////////////////
		MkBkt();

		ColForce();

		update();
		//////////////////////DEM///////////////////////////
#endif


		outtime += dt;
		TIM += dt;
		OPT_FQC += 1.0f;
		//printf_s("time=%f\n", TIM);
	}


}


void DEMPS::memory_free() {
	free(Pos.x); free(Pos.y); free(Pos.z);
	free(Vel.x); free(Vel.y); free(Vel.z);
	free(Omega.x); free(Omega.y); free(Omega.z);
	free(Ftotal.x); free(Ftotal.y); free(Ftotal.z);
	free(Torque.x); free(Torque.y); free(Torque.z);
	free(Typ);
	free(ep.x); free(ep.y); free(ep.z);
	free(D);
	free(pair);

	cudaFree(d_Pos.x); cudaFree(d_Pos.y); cudaFree(d_Pos.z);
	cudaFree(d_Vel.x); cudaFree(d_Vel.y); cudaFree(d_Vel.z);
	cudaFree(d_Omega.x); cudaFree(d_Omega.y); cudaFree(d_Omega.z);
	cudaFree(d_Ftotal.x); cudaFree(d_Ftotal.y); cudaFree(d_Ftotal.z);
	cudaFree(d_Torque.x); cudaFree(d_Torque.y); cudaFree(d_Torque.z);
	cudaFree(d_Typ);
	cudaFree(d_ep.x); cudaFree(d_ep.y); cudaFree(d_ep.z);
	cudaFree(d_D);
	cudaFree(d_pair);


	free(Acc.x);	free(Acc.y);	free(Acc.z);
	free(Prs);
	free(pav);
	free(TypM);
	free(PosM.x);		free(PosM.y);		free(PosM.z);
	free(VelM.x);		free(VelM.y);		free(VelM.z);
	free(PrsM);
	free(Dns);
	free(WLLSE);
	free(WLLVec.x);		free(WLLVec.y);		free(WLLVec.z);

	cudaFree(d_Acc.x);	cudaFree(d_Acc.y);	cudaFree(d_Acc.z);
	cudaFree(d_Prs);
	cudaFree(d_pav);
	cudaFree(d_TypM);
	cudaFree(d_PosM.x);		cudaFree(d_PosM.y);		cudaFree(d_PosM.z);
	cudaFree(d_VelM.x);		cudaFree(d_VelM.y);		cudaFree(d_VelM.z);
	cudaFree(d_PrsM);
	cudaFree(d_Dns);
	cudaFree(d_WLLSE);
	cudaFree(d_WLLVec.x);		cudaFree(d_WLLVec.y);		cudaFree(d_WLLVec.z);
}