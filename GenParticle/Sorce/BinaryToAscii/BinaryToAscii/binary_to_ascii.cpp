#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>

using namespace std;

int main() {

	char comment[80];
	unsigned int num_triangles;
	unsigned short reserved;

	float vec_x, vec_y, vec_z;
	float pos_x1, pos_x2, pos_x3;
	float pos_y1, pos_y2, pos_y3;
	float pos_z1, pos_z2, pos_z3;

	ifstream IN("initial.stl", ios::in | ios::binary);
	//ifstream IN("tank_full.stl", ios::in | ios::binary);
	//ifstream IN("tank_up.stl", ios::in | ios::binary);
	//ifstream IN("tank_down.stl", ios::in | ios::binary);


	if (!IN) {
		cout << "ファイル読み込み失敗" << endl;
		return -1;
	}

	ofstream OUT("initial.csv");

	OUT << " i , f , x , y , z , nx , ny , nz , nc  " << endl;

	IN.read(comment, 80);

	IN.read((char *)&num_triangles, sizeof(unsigned int));


	cout << "number_of_triangles :" << num_triangles << endl;
	FILE *output3;
	output3 = fopen("initial.txt", "a");
	//output3 = fopen("initialwall.txt", "a");

	fprintf(output3, "%d\n", num_triangles);
	for (int i = 0; i < num_triangles; i++)
	{
		IN.read((char *)&vec_x, sizeof(float));
		IN.read((char *)&vec_y, sizeof(float));
		IN.read((char *)&vec_z, sizeof(float));
		IN.read((char *)&pos_x1, sizeof(float));
		IN.read((char *)&pos_y1, sizeof(float));
		IN.read((char *)&pos_z1, sizeof(float));
		IN.read((char *)&pos_x2, sizeof(float));
		IN.read((char *)&pos_y2, sizeof(float));
		IN.read((char *)&pos_z2, sizeof(float));
		IN.read((char *)&pos_x3, sizeof(float));
		IN.read((char *)&pos_y3, sizeof(float));
		IN.read((char *)&pos_z3, sizeof(float));
		IN.read((char *)&reserved, sizeof(unsigned short));

		float g_x = (pos_x1 + pos_x2 + pos_x3) / 3.0;//重心座標
		float g_y = (pos_y1 + pos_y2 + pos_y3) / 3.0;
		float g_z = (pos_z1 + pos_z2 + pos_z3) / 3.0;

		float vec_g_x = vec_x;//法線ベクトル
		float vec_g_y = vec_y;
		float vec_g_z = vec_z;

		float distx1, distx2, distx3;
		float disty1, disty2, disty3;
		float distz1, distz2, distz3;
		float dist1, dist2, dist3;

		distx1 = abs(pos_x1 - pos_x2);
		disty1 = abs(pos_y1 - pos_y2);
		distz1 = abs(pos_z1 - pos_z2);

		distx2 = abs(pos_x2 - pos_x3);
		disty2 = abs(pos_y2 - pos_y3);
		distz2 = abs(pos_z2 - pos_z3);

		distx3 = abs(pos_x3 - pos_x1);
		disty3 = abs(pos_y3 - pos_y1);
		distz3 = abs(pos_z3 - pos_z1);

		dist1 = sqrt(distx1*distx1 + disty1*disty1 + distz1*distz1);
		dist2 = sqrt(distx2*distx2 + disty2*disty2 + distz2*distz2);
		dist3 = sqrt(distx3*distx3 + disty3*disty3 + distz3*distz3);

		float s = (dist1 + dist2 + dist3) / 2.0;
		float T = sqrt(s*(s - dist1)*(s - dist2)*(s - dist3));//三角形メッシュ面積（ヘロンの公式）
		float r = T / s;//内接円半径

		float max_dist = 0;

		if (max_dist < dist1)
			max_dist = dist1;

		if (max_dist < dist2)
			max_dist = dist2;

		if (max_dist < dist3)
			max_dist = dist3;

		float aspect_ratio = max_dist / (2 * sqrt(3)*r);//三角形メッシュ縦横比
		//g_x /= 10.0f;  g_y /= 10.0f; g_z /= 10.0f;
		
		//if ((T > 0.0000000002) && (aspect_ratio < 30000.0))
		{
	//		OUT << i << "," << 1 << "," << g_x << "," << g_y << "," << g_z << "," << vec_g_x << "," << vec_g_y << "," << vec_g_z << "," << 0.5 << endl;
		
	//	fprintf(output3, "%d	%d	%f	%f	%f	%f	%f	%f	%f	%f	\n", i, 1, g_x, g_y, g_z, vec_g_x, vec_g_y, vec_g_z, 0.5, T);
		fprintf(output3, "%d	%d	%10.4e	%10.4e	%10.4e	%10.4e	%10.4e	%10.4e	%10.4e\n", i, 2, g_x, g_y, g_z, -vec_g_x, -vec_g_y, -vec_g_z, T);
		}
 	}


	OUT.close();
	IN.close();

	system("pause");
	return 0;
}
