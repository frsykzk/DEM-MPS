#include "class.cuh"

int main(int argc, char** argv) {
	std::chrono::system_clock::time_point  start, end; // �^�� auto �ŉ�
	start = std::chrono::system_clock::now(); // �v���J�n����

	if (Multi_flg == 1) { printf_s("��������͂����s���܂�.\n\n"); }
	else{ 
		if (DEM_flg == 1) { printf_s("DEM�����s���܂�.\n\n"); }
		if (MPS_flg == 1) { printf_s("MPS�����s���܂�.\n\n"); }
	}
	
	printf_s("�g�p�\�ȍő�X���b�h���F%d\n", omp_get_max_threads());
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

	end = std::chrono::system_clock::now();  // �v���I������
	double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); //�����ɗv�������Ԃ�b�ɕϊ�
	std::cout << elapsed << " second" << std::endl;

	int t = elapsed;
	int h = t / 3600;   t %= 3600;
	int m = t / 60;     t %= 60;
	int s = t;
	std::cout << h << "h " << m << "m " << s << "s " << std::endl;

	obj.memory_free();

	return 0;
}
