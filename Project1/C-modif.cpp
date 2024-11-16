#pragma comment(lib,"graphics.lib")
#include "graphics.h"
#include <iostream>
#include <iomanip>         // ���������ʽ
#include <vector>			// ʹ�� vector ����
#include <algorithm>		// ����
#include <random>			// ���������
#include <fstream>		    // �ļ�����

using namespace std;

////////////////////////FillRectangleWithShapes//////////////////////////
//����һ�����࣬�����ɴ���չ���������
class FillRectangleWithShapes   
{
public:
	double region_width, region_height;		// ���ο���
	double total_area;									// �������
	double region[8];										// ��������
	int num_sieve;											// ����ȼ���Ŀ
	double** grading;                  // �ۼ������������ά����
	double filled_area_ratio;       // �Ѿ�������
	double needed_area;           // Ŀ��������
	char projname[50];               // ��Ŀ����

public:
	FillRectangleWithShapes(double _region_width, double _region_height,
		int _num_sieve, double _grading[][2], double filled_area_ratio,
		const char* _projname = "Noname");                   //���湹�캯��

	FillRectangleWithShapes(const char* filename);  //��ȡ�����ļ��������

	void setname(const char* s) { strcpy_s(projname, s); } //����Ŀ��������

	~FillRectangleWithShapes();                         //�����ٹ�������ռ�
};

FillRectangleWithShapes::FillRectangleWithShapes(const char * projname)
{
	char fn[54]; //�ļ�����

	strcpy_s(this->projname, 50, projname);  //��Ŀ����

	//Ĭ�϶����ļ�������Ŀ����ͬ������׺Ϊin����projname.in
	strcpy_s(fn, 50, this->projname);
	strcat_s(fn, 50, ".in");

	ifstream fin(fn);
	fin >> this->region_width >> this->region_height;  //��ȡ���
	
	fin >> this->num_sieve; //��ȡ����ȼ�
	
	int i = 0;
	for (auto each : { 0.0, 0.0, this->region_width, 0.0, this->region_width, 
		this->region_height, 0.0, this->region_height })
		this->region[i++] = each;

	this->grading = new double* [num_sieve];  //��̬�������飬��Ҫ�ͷţ�

	for (int i = 0; i < this->num_sieve; i++)
	{
		this->grading[i] = new double[2];
		fin >> this->grading[i][0] >> this->grading[i][1];
		//cout << this->grading[i][0] <<"    "<< this->grading[i][1] << endl;
	}

	fin >> this->filled_area_ratio;  //��ȡ�����
	this->total_area = region_width * region_height * this->filled_area_ratio;
	//cout << this->filled_area_ratio << endl;
	fin.close();
}

FillRectangleWithShapes::FillRectangleWithShapes(double _region_width, 
	double _region_height, int _num_sieve, double _grading[][2], 
	double filled_area_ratio, const char* _projname)
{
	this->region_width = _region_width;
	this->region_height = _region_height;
	strcpy_s(this->projname, _projname);
	this->filled_area_ratio = filled_area_ratio;
	this->total_area = region_width * region_height * this->filled_area_ratio;

	int i = 0;
	for (auto each : { 0.0, 0.0, this->region_width, 0.0, this->region_width, 
		this->region_height, 0.0, this->region_height })
		this->region[i++] = each;
	this->num_sieve = _num_sieve;
	this->grading = new double* [num_sieve];  //��̬�������飬��Ҫ�ͷţ�
	for (i = 0; i < num_sieve; i++)
	{
		this->grading[i] = new double[2];
		this->grading[i][0] = _grading[i][0];
		this->grading[i][1] = _grading[i][1];
	}
}

FillRectangleWithShapes::~FillRectangleWithShapes()
{
	//�����������Դ�������ʱ����Ŀռ���л���
	//���ɺ�ϰ��
	for (int i = 0; i < this->num_sieve; i++)
	{
		delete this->grading[i];
	}
	delete this->grading;
}

////////////////////////FillRectangleWithCircles/////////////////////
class FillRectangleWithCircles :public FillRectangleWithShapes
{
private:
	vector<double> radius, x_pos, y_pos;        //�洢Բ����Ϣ
	void generate_circle_radius();  //�������Բ�İ뾶
	void generate_circle_position(); //ȷ��Բ��λ��
	bool is_circle_in_region(double xi, double yi, double  r);
	bool is_overlap_two_circles(double x1, double y1, double r1,
		double x2, double y2, double r2, double tor=1.0);
	
	//�б�ĳһ��Բ���������Բ�Ƿ����ص�
	bool is_overlap_circle_circles(double x1, double y1, double r1, int istop);
	
public:
	FillRectangleWithCircles(double _region_width, double _region_height, 
		int _num_sieve, double _grading[][2], double filled_area_ratio, 
		const char* _projname = "Noname") :\
		FillRectangleWithShapes(_region_width, _region_height, _num_sieve,
			_grading, filled_area_ratio, _projname)
	{
		; //�̳��˸���Ĺ��캯��
	}
	//FillRectangleWithShapes(const char *filename) 
	FillRectangleWithCircles(const char* _projname) :
		FillRectangleWithShapes( _projname)
	{
		; //�̳��˸���Ĺ��캯��2
	}

	void out2file(); //������ļ�
	void out2graph(int wind_w = 1920, int win_h = 1080); //�������Ļ
	void generator(); //����������еĹ���
};

bool FillRectangleWithCircles::is_overlap_circle_circles(double x1, double y1, 
	double r1, int istop)
{   //�ص�false  ���ص���true
	// istop ��ʾ��0����istop��Բ���뵱ǰԲ�Ƚ�
	if (istop <= 0) return true;
	bool isoverlap = true; // 
	for (int i = 0; i <= istop; i++)
	{
		if (is_overlap_two_circles(x1, y1, r1, 
			this->x_pos[i], this->y_pos[i], this->radius[i]))
		{
			return false;
		}
	}
	return isoverlap;
}

bool FillRectangleWithCircles::is_overlap_two_circles(double x1, double y1, 
	double r1, double x2, double y2, double r2, double tor)
{
	//�ж�����Բ�Ƿ��ཻ
	//tor��һ�����ڵ���1��ϵ��������1��ʾ����Բ֮����һ����϶
	//tor=1��ʾ��϶����Ϊ0
	double dist = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));
	if (dist < 1.1 * (r1 + r2)) 
		return true;
	else
		return false;
}

bool FillRectangleWithCircles::is_circle_in_region(double xi, double yi, double  ri)
{
	// �ж�Բ�Ƿ���ȫ�� region ��
	// ������򷵻�true�����򷵻�false
	if ((xi - ri < region[0]) || (xi + ri > region[2]) || 
		(yi - ri < region[1]) || (yi + ri > region[5]))
		return false;
	else
		return true;
}

void FillRectangleWithCircles::generate_circle_radius()
{
	const double EPS = 1e-3;
	const double PI = acos(-1);
	vector<double> needed_area;	// ��Ҫ��Ŀ����������
								// ÿһ���ȼ�ֱ������Ҫ�������
	vector<double> occupied_area;	// �Ѿ��������
								// ÿ���ȼ�Բ�Ѿ�ռ�е����
	vector<double> tmp;      
	double min_rad, max_rad, r;

	//���ݼ��䣬��ʼ��ÿ���ȼ���������Ҫ����������������
	for (int i = 0; i < num_sieve - 1; i++)
	{
		needed_area.push_back((grading[i + 1][1] - grading[i][1]) * total_area);
		occupied_area.push_back(0.0);
	}

	random_device rnd; // ����ϵͳ��CPUӲ��ȡ�������
	mt19937 rng(rnd()); // ��������Ӷ��������������
	uniform_real_distribution<double> uni(0.0, 1.0); // ������ȷֲ�

	for (int i = num_sieve - 2; i >= 0; i--)
	{
		// ���չ��ϼ���ѭ��

		// ��ǰ������С�ĳߴ磨ֱ����
		min_rad = grading[i][0];
		// ��ǰ�������ĳߴ磨ֱ����
		max_rad = grading[i + 1][0];

		// ѭ���жϵ�ǰ��������Ƿ����������0.1%
		// ȷ�� occupied_area �� needed_area һ����
		// ��������Բ������������ֵ��ܶ�
		while (abs(occupied_area[i] - needed_area[i]) / needed_area[i] > EPS)
		{
			occupied_area[i] = 0;
			tmp.clear();  // ��ʱ����ÿ���뾶����������ɵİ뾶
			while (occupied_area[i] < needed_area[i])
			{
				// uni(rng)����һ��0-1���������
				r = uni(rng) * (max_rad - min_rad) + min_rad;
				//cout << r << endl;
				tmp.push_back(r);
				occupied_area[i] = occupied_area[i] + r * r * PI;
			}
		}
		 
		for (auto iter = tmp.begin(); iter != tmp.end(); iter++)
			radius.push_back(*iter);
	}
	sort(radius.begin(), radius.end(), greater<double>()); // �Ӵ�С����
	//this->radius = radius;
}

// region:   width (x ����), height (y ����)
void FillRectangleWithCircles::generate_circle_position()
{
	double x, y;                          // ��������
	double min_x = region[0];    // ������߽�
	double max_x = region[2];   // �����ұ߽�
	double min_y = region[1];    // �����±߽�
	double max_y = region[5];   // �����ϱ߽�

	random_device rnd; // ����ϵͳ��CPUӲ��ȡ�������
	mt19937 rng(rnd()); // ��������Ӷ�������������������Զ�����������������
	uniform_real_distribution<double> uni(0.0, 1.0); // ������ȷֲ�

	for (size_t i = 0; i < radius.size(); i++)
	{
		cout << "\r��ǰ����Բ��λ�ã���ɣ�" << setw(6) << fixed <<
			setprecision(2) << (i + 1.0) / radius.size() * 100 << "% \t";

		do
		{   //�������λ��
			x = uni(rng) * (max_x - min_x) + min_x;
			y = uni(rng) * (max_y - min_y) + min_y;
		} while (!(is_circle_in_region(x, y, radius[i]) &&  
			is_overlap_circle_circles(x, y, radius[i], i - 1)));
		// ��������������������

		x_pos.push_back(x);
		y_pos.push_back(y);
	}
	cout << endl;
}

void FillRectangleWithCircles::out2file() //������ļ�
{
	char filename[54];
	strcpy_s(filename, 50, this->projname);
	strcat_s(filename, 50, ".txt");
	ofstream file(filename);
	if (file.is_open())
	{
		for (int i = 0; i < 8; i++)
			file << region[i] << '\t';
		file << endl;
		for (size_t i = 0; i < radius.size(); i++)
		{
			file << 1 << '\t' << x_pos[i] << '\t' << y_pos[i] << '\t' 
				<< radius[i] << endl;
		}
	}
	file.close();
}

void FillRectangleWithCircles::out2graph(int win_w, int win_h) //������ļ�
{
	// ����Ϊ��ͼģ��
	// ��ʼ��ͼ�δ���
	initwindow(win_w, win_h, this->projname);
	
	// ���ñ���Ϊ��ɫ
	setbkcolor(WHITE);
	cleardevice(); // ����豸��Ӧ�ñ�����ɫ

	// ���û�ͼ��ɫΪ��ɫ
	setcolor(RED);
	setfillstyle(SOLID_FILL, RED);

	double scale = 0; //ͼ��������ʾ����
	double maxw = region[2] - region[0];// ���
	double maxh = region[7] - region[1];//�߶�
	if (maxw / maxh > 1.0 * win_w / win_h)//ͼ�����С
		scale = maxw / win_w;
	else
		scale = maxh / win_h;

	int intregion[8 + 2];  //��ʾ���������飬��Ҫ�ص�ԭ��
	for (size_t i = 0; i < 8 + 2; i++)
	{
		if (i < 8)
			if (i % 2 == 0) //���������
				intregion[i] = int((region[i] - region[0]) / scale); 
			else               //����������
				intregion[i] = win_h - int((region[i] - region[1]) / scale);  
		else
			intregion[i] = intregion[i - 8];
	}

	// ���Ʋ���䱳�������
	drawpoly(4, intregion);
	floodfill(int(intregion[0] / 2.0 + intregion[2] / 2.0),
		int(intregion[1] / 2.0 + intregion[7] / 2.0), WHITE); // �����ɫ

	// ����ͬ��Բ��
	for (size_t i = 0; i < radius.size(); i++)
	{
		int xi = int((x_pos[i] - region[0]) / scale);
		int yi = win_h - int((y_pos[i] - region[1]) / scale);
		int ri = int(radius[i] / scale);
		//floodfill(xi, yi, WHITE); // ��䱳��ɫ
		circle(xi, yi, ri);       // Բ�����ĺͰ뾶
		floodfill(xi, yi, RED);   // ���Ϊ��ɫ
	}

	// ��ͼ
	char fn[54];
	strcpy_s(fn, 50, this->projname);
	strcat_s(fn, 50, ".bmp");

	writeimagefile(fn, 0, 0, getmaxx(), getmaxy()); //��֧��bmp��ʽ���

	// �ȴ��û��������˳�
	getch();
	closegraph();
}

void FillRectangleWithCircles::generator() //������ļ�
{
	this->generate_circle_radius();
	cout << "��Ŀ��" << this->projname 
		<< "�����й���Բ�İ뾶������ɣ�" << endl;

	this->generate_circle_position();

	cout << "��Ŀ��" << this->projname
		<< "�����й���Բ��λ��������ɣ�" << endl;

	this->out2file();
	cout << "��Ŀ��" << this->projname 
		<< "�����������ļ���" << endl;

	this->out2graph(1320, 1080);
	cout << "��Ŀ��" << this->projname 
		<< "����Ѿ���ͼ��" << endl;

	cout << "��Ŀ��" << this->projname 
		<< "�Ѿ�ȫ����ϣ�" << endl;

	cout << "=========THE END==========" << endl;
}


////////////////////////////////FillRectangleWithZtriangle////////
class FillRectangleWithZtriangle :public FillRectangleWithShapes
{
private:
	vector<double> radius, x_pos, y_pos, angle;        //�洢�����ε���Ϣ
	void generate_Ztriangle_radius();  //�������Բ�İ뾶
	void generate_Ztriangle_position(); //ȷ��Բ��λ��
	bool is_Ztriangle_in_region(double xi, double yi, double  r);
	bool is_overlap_two_Ztriangle(double x1, double y1, double r1,
		double x2, double y2, double r2, double tor = 1.0);

	//�б�ĳһ��Բ���������Բ�Ƿ����ص�
	bool is_overlap_Ztriangle_Ztriangles(double x1, double y1, double r1, int istop);

public:
	FillRectangleWithZtriangle(double _region_width, double _region_height,
		int _num_sieve, double _grading[][2], double filled_area_ratio,
		const char* _projname = "Noname") :\
		FillRectangleWithShapes(_region_width, _region_height, _num_sieve,
			_grading, filled_area_ratio, _projname)
	{
		; //�̳��˸���Ĺ��캯��
	}
	//FillRectangleWithShapes(const char *filename) 
	FillRectangleWithZtriangle(const char* _projname) :
		FillRectangleWithShapes(_projname)
	{
		; //�̳��˸���Ĺ��캯��2
	}

	void out2file(); //������ļ�
	void out2graph(int wind_w = 1920, int win_h = 1080); //�������Ļ
	void generator(); //����������еĹ���
};


void FillRectangleWithZtriangle::generate_Ztriangle_radius()  //�������Բ�İ뾶
{

}

int main()
{
	double region_w = 20, region_h = 20;
	const int num_sieve = 4;
	double grading[num_sieve][2] = { {0.2, 0},{0.5, 0.3},{1.0, 0.8},{1.5, 1.0} };
	double areafilledratio = 0.55; //������ռ��

	FillRectangleWithShapes fws01(region_w, region_h, num_sieve, grading, 
		areafilledratio); //test ���࣬���������塣

	FillRectangleWithCircles fwc01(region_w, region_h, num_sieve, grading, 
	areafilledratio, "Test01");

	FillRectangleWithCircles fwc02("Test02");

	fwc01.generator();

	//fwc02.generator();

	//FillRectangleWithZtriangle fwzt(region_w, region_h, num_sieve, grading,areafilledratio);
	//fwzt.generator();

	return 0;
}