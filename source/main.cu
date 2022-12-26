#include "class.cuh"

int main(int argc, char** argv) {
	std::chrono::system_clock::time_point  start, end; // 型は auto で可
	start = std::chrono::system_clock::now(); // 計測開始時間

	if (Multi_flg == 1) { printf_s("混相流解析を実行します.\n\n"); }
	else{ 
		if (DEM_flg == 1) { printf_s("DEMを実行します.\n\n"); }
		if (MPS_flg == 1) { printf_s("MPSを実行します.\n\n"); }
	}
	
	printf_s("使用可能な最大スレッド数：%d\n", omp_get_max_threads());
#pragma omp parallel
	for (int i = 0; i < omp_get_max_threads(); i++) {
		printf_s("Hello!! CPU Thread %d\n", i);
	}
	printf_s("\n");

	DEMPS obj;

	obj.RdDat();
	printf_s("RdDat finished!\n\n");

	obj.AlcBkt();
	printf_s("AlkBkt finished!\n\n");

	obj.SetPara();
	printf_s("SetPara finished!\n\n");

	obj.ClcDEMPS();

	printf_s("Coaculation finished!!!\n\n");

	end = std::chrono::system_clock::now();  // 計測終了時間
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //処理に要した時間を秒に変換
	std::cout << elapsed << " second" << std::endl;

	int t = elapsed;
	int h = t / 3600;   t %= 3600;
	int m = t / 60;     t %= 60;
	int s = t;
	std::cout << h << "h " << m << "m " << s << "s " << std::endl;

	obj.memory_free();

	return 0;
}
