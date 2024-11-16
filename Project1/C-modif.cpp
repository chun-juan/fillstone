#pragma comment(lib,"graphics.lib")
#include "graphics.h"
#include <iostream>
#include <iomanip>         // 控制输出格式
#include <vector>			// 使用 vector 容器
#include <algorithm>		// 排序
#include <random>			// 生成随机数
#include <fstream>		    // 文件操作

using namespace std;

////////////////////////FillRectangleWithShapes//////////////////////////
//定义一个基类，将来由此拓展出多个子类
class FillRectangleWithShapes   
{
public:
	double region_width, region_height;		// 矩形宽、高
	double total_area;									// 矩形面积
	double region[8];										// 矩形坐标
	int num_sieve;											// 级配等级数目
	double** grading;                  // 累计粒径面积，二维数组
	double filled_area_ratio;       // 已经填充面积
	double needed_area;           // 目标填充面积
	char projname[50];               // 项目名称

public:
	FillRectangleWithShapes(double _region_width, double _region_height,
		int _num_sieve, double _grading[][2], double filled_area_ratio,
		const char* _projname = "Noname");                   //常规构造函数

	FillRectangleWithShapes(const char* filename);  //读取控制文件构造对象

	void setname(const char* s) { strcpy_s(projname, s); } //给项目设置名字

	~FillRectangleWithShapes();                         //处理开辟过的数组空间
};

FillRectangleWithShapes::FillRectangleWithShapes(const char * projname)
{
	char fn[54]; //文件名字

	strcpy_s(this->projname, 50, projname);  //项目名称

	//默认读入文件名与项目名称同名，后缀为in，如projname.in
	strcpy_s(fn, 50, this->projname);
	strcat_s(fn, 50, ".in");

	ifstream fin(fn);
	fin >> this->region_width >> this->region_height;  //读取宽高
	
	fin >> this->num_sieve; //读取级配等级
	
	int i = 0;
	for (auto each : { 0.0, 0.0, this->region_width, 0.0, this->region_width, 
		this->region_height, 0.0, this->region_height })
		this->region[i++] = each;

	this->grading = new double* [num_sieve];  //动态开辟数组，需要释放！

	for (int i = 0; i < this->num_sieve; i++)
	{
		this->grading[i] = new double[2];
		fin >> this->grading[i][0] >> this->grading[i][1];
		//cout << this->grading[i][0] <<"    "<< this->grading[i][1] << endl;
	}

	fin >> this->filled_area_ratio;  //读取填充率
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
	this->grading = new double* [num_sieve];  //动态开辟数组，需要释放！
	for (i = 0; i < num_sieve; i++)
	{
		this->grading[i] = new double[2];
		this->grading[i][0] = _grading[i][0];
		this->grading[i][1] = _grading[i][1];
	}
}

FillRectangleWithShapes::~FillRectangleWithShapes()
{
	//析构函数，对创建对象时申请的空间进行回收
	//养成好习惯
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
	vector<double> radius, x_pos, y_pos;        //存储圆的信息
	void generate_circle_radius();  //随机生成圆的半径
	void generate_circle_position(); //确定圆的位置
	bool is_circle_in_region(double xi, double yi, double  r);
	bool is_overlap_two_circles(double x1, double y1, double r1,
		double x2, double y2, double r2, double tor=1.0);
	
	//判别某一个圆与既有所有圆是否有重叠
	bool is_overlap_circle_circles(double x1, double y1, double r1, int istop);
	
public:
	FillRectangleWithCircles(double _region_width, double _region_height, 
		int _num_sieve, double _grading[][2], double filled_area_ratio, 
		const char* _projname = "Noname") :\
		FillRectangleWithShapes(_region_width, _region_height, _num_sieve,
			_grading, filled_area_ratio, _projname)
	{
		; //继承了父类的构造函数
	}
	//FillRectangleWithShapes(const char *filename) 
	FillRectangleWithCircles(const char* _projname) :
		FillRectangleWithShapes( _projname)
	{
		; //继承了父类的构造函数2
	}

	void out2file(); //输出到文件
	void out2graph(int wind_w = 1920, int win_h = 1080); //输出到屏幕
	void generator(); //打包集成所有的功能
};

bool FillRectangleWithCircles::is_overlap_circle_circles(double x1, double y1, 
	double r1, int istop)
{   //重叠false  不重叠是true
	// istop 表示从0到第istop个圆，与当前圆比较
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
	//判断两个圆是否相交
	//tor是一个大于等于1的系数，大于1表示两个圆之间有一个间隙
	//tor=1表示间隙可以为0
	double dist = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));
	if (dist < 1.1 * (r1 + r2)) 
		return true;
	else
		return false;
}

bool FillRectangleWithCircles::is_circle_in_region(double xi, double yi, double  ri)
{
	// 判断圆是否完全在 region 中
	// 如果在则返回true，否则返回false
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
	vector<double> needed_area;	// 需要的目标填充面积；
								// 每一个等级直径，需要的总面积
	vector<double> occupied_area;	// 已经填充的面积
								// 每个等级圆已经占有的面积
	vector<double> tmp;      
	double min_rad, max_rad, r;

	//根据级配，初始化每个等级骨料所需要的面积、已填充的面积
	for (int i = 0; i < num_sieve - 1; i++)
	{
		needed_area.push_back((grading[i + 1][1] - grading[i][1]) * total_area);
		occupied_area.push_back(0.0);
	}

	random_device rnd; // 操作系统，CPU硬件取随机种子
	mt19937 rng(rnd()); // 用随机种子定义随机数生成器
	uniform_real_distribution<double> uni(0.0, 1.0); // 定义均匀分布

	for (int i = num_sieve - 2; i >= 0; i--)
	{
		// 按照骨料级配循环

		// 当前级配最小的尺寸（直径）
		min_rad = grading[i][0];
		// 当前级配最大的尺寸（直径）
		max_rad = grading[i + 1][0];

		// 循环判断当前级配骨料是否填满，误差0.1%
		// 确保 occupied_area 与 needed_area 一样，
		// 否则最后的圆总面积会比期望值大很多
		while (abs(occupied_area[i] - needed_area[i]) / needed_area[i] > EPS)
		{
			occupied_area[i] = 0;
			tmp.clear();  // 临时保存每个半径区间随机生成的半径
			while (occupied_area[i] < needed_area[i])
			{
				// uni(rng)产生一个0-1的随机数；
				r = uni(rng) * (max_rad - min_rad) + min_rad;
				//cout << r << endl;
				tmp.push_back(r);
				occupied_area[i] = occupied_area[i] + r * r * PI;
			}
		}
		 
		for (auto iter = tmp.begin(); iter != tmp.end(); iter++)
			radius.push_back(*iter);
	}
	sort(radius.begin(), radius.end(), greater<double>()); // 从大到小排序
	//this->radius = radius;
}

// region:   width (x 方向), height (y 方向)
void FillRectangleWithCircles::generate_circle_position()
{
	double x, y;                          // 生成坐标
	double min_x = region[0];    // 区域左边界
	double max_x = region[2];   // 区域右边界
	double min_y = region[1];    // 区域下边界
	double max_y = region[5];   // 区域上边界

	random_device rnd; // 操作系统，CPU硬件取随机种子
	mt19937 rng(rnd()); // 用随机种子定义随机数生成器（可以定义其他的生成器）
	uniform_real_distribution<double> uni(0.0, 1.0); // 定义均匀分布

	for (size_t i = 0; i < radius.size(); i++)
	{
		cout << "\r当前生成圆的位置，完成：" << setw(6) << fixed <<
			setprecision(2) << (i + 1.0) / radius.size() * 100 << "% \t";

		do
		{   //产生随机位置
			x = uni(rng) * (max_x - min_x) + min_x;
			y = uni(rng) * (max_y - min_y) + min_y;
		} while (!(is_circle_in_region(x, y, radius[i]) &&  
			is_overlap_circle_circles(x, y, radius[i], i - 1)));
		// 不满足条件，反复迭代

		x_pos.push_back(x);
		y_pos.push_back(y);
	}
	cout << endl;
}

void FillRectangleWithCircles::out2file() //输出到文件
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

void FillRectangleWithCircles::out2graph(int win_w, int win_h) //输出到文件
{
	// 以下为绘图模块
	// 初始化图形窗口
	initwindow(win_w, win_h, this->projname);
	
	// 设置背景为白色
	setbkcolor(WHITE);
	cleardevice(); // 清空设备以应用背景颜色

	// 设置绘图颜色为红色
	setcolor(RED);
	setfillstyle(SOLID_FILL, RED);

	double scale = 0; //图像缩放显示比例
	double maxw = region[2] - region[0];// 宽度
	double maxh = region[7] - region[1];//高度
	if (maxw / maxh > 1.0 * win_w / win_h)//图像宽大高小
		scale = maxw / win_w;
	else
		scale = maxh / win_h;

	int intregion[8 + 2];  //显示的轮廓数组，需要回到原点
	for (size_t i = 0; i < 8 + 2; i++)
	{
		if (i < 8)
			if (i % 2 == 0) //处理横坐标
				intregion[i] = int((region[i] - region[0]) / scale); 
			else               //处理纵坐标
				intregion[i] = win_h - int((region[i] - region[1]) / scale);  
		else
			intregion[i] = intregion[i - 8];
	}

	// 绘制并填充背景多边形
	drawpoly(4, intregion);
	floodfill(int(intregion[0] / 2.0 + intregion[2] / 2.0),
		int(intregion[1] / 2.0 + intregion[7] / 2.0), WHITE); // 填充颜色

	// 画不同的圆块
	for (size_t i = 0; i < radius.size(); i++)
	{
		int xi = int((x_pos[i] - region[0]) / scale);
		int yi = win_h - int((y_pos[i] - region[1]) / scale);
		int ri = int(radius[i] / scale);
		//floodfill(xi, yi, WHITE); // 填充背景色
		circle(xi, yi, ri);       // 圆的中心和半径
		floodfill(xi, yi, RED);   // 填充为红色
	}

	// 出图
	char fn[54];
	strcpy_s(fn, 50, this->projname);
	strcat_s(fn, 50, ".bmp");

	writeimagefile(fn, 0, 0, getmaxx(), getmaxy()); //仅支持bmp格式输出

	// 等待用户按键后退出
	getch();
	closegraph();
}

void FillRectangleWithCircles::generator() //输出到文件
{
	this->generate_circle_radius();
	cout << "项目：" << this->projname 
		<< "的所有骨料圆的半径生成完成！" << endl;

	this->generate_circle_position();

	cout << "项目：" << this->projname
		<< "的所有骨料圆的位置生成完成！" << endl;

	this->out2file();
	cout << "项目：" << this->projname 
		<< "结果已输出到文件！" << endl;

	this->out2graph(1320, 1080);
	cout << "项目：" << this->projname 
		<< "结果已经出图！" << endl;

	cout << "项目：" << this->projname 
		<< "已经全部完毕！" << endl;

	cout << "=========THE END==========" << endl;
}


////////////////////////////////FillRectangleWithZtriangle////////
class FillRectangleWithZtriangle :public FillRectangleWithShapes
{
private:
	vector<double> radius, x_pos, y_pos, angle;        //存储三角形的信息
	void generate_Ztriangle_radius();  //随机生成圆的半径
	void generate_Ztriangle_position(); //确定圆的位置
	bool is_Ztriangle_in_region(double xi, double yi, double  r);
	bool is_overlap_two_Ztriangle(double x1, double y1, double r1,
		double x2, double y2, double r2, double tor = 1.0);

	//判别某一个圆与既有所有圆是否有重叠
	bool is_overlap_Ztriangle_Ztriangles(double x1, double y1, double r1, int istop);

public:
	FillRectangleWithZtriangle(double _region_width, double _region_height,
		int _num_sieve, double _grading[][2], double filled_area_ratio,
		const char* _projname = "Noname") :\
		FillRectangleWithShapes(_region_width, _region_height, _num_sieve,
			_grading, filled_area_ratio, _projname)
	{
		; //继承了父类的构造函数
	}
	//FillRectangleWithShapes(const char *filename) 
	FillRectangleWithZtriangle(const char* _projname) :
		FillRectangleWithShapes(_projname)
	{
		; //继承了父类的构造函数2
	}

	void out2file(); //输出到文件
	void out2graph(int wind_w = 1920, int win_h = 1080); //输出到屏幕
	void generator(); //打包集成所有的功能
};


void FillRectangleWithZtriangle::generate_Ztriangle_radius()  //随机生成圆的半径
{

}

int main()
{
	double region_w = 20, region_h = 20;
	const int num_sieve = 4;
	double grading[num_sieve][2] = { {0.2, 0},{0.5, 0.3},{1.0, 0.8},{1.5, 1.0} };
	double areafilledratio = 0.55; //填充面积占比

	FillRectangleWithShapes fws01(region_w, region_h, num_sieve, grading, 
		areafilledratio); //test 基类，无其他意义。

	FillRectangleWithCircles fwc01(region_w, region_h, num_sieve, grading, 
	areafilledratio, "Test01");

	FillRectangleWithCircles fwc02("Test02");

	fwc01.generator();

	//fwc02.generator();

	//FillRectangleWithZtriangle fwzt(region_w, region_h, num_sieve, grading,areafilledratio);
	//fwzt.generator();

	return 0;
}