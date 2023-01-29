#ifndef MinecraftParticleMiDiPlayerLines
#define MinecraftParticleMiDiPlayerLines namespace MLinesParticle{
#define MinecraftParticleMiDiPlayerLinesEnd }
#pragma warning(disable : 4996)
#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <windows.h>
#include <vector>
#include <cmath>
#include <cstdio>
using namespace std;

MinecraftParticleMiDiPlayerLines
class BaseLines
{
public:
	using fp64 = double;
	using fp32 = float;
	using int64 = long long;
	BaseLines() = default;
	BaseLines(const string& pty) : ParticleType(pty) {}
	BaseLines(string&& pty) : ParticleType(pty) {}

	void open(FILE* file) { OutPutFile = file; }

	void ChangeBaseLambda(fp64 NewLambda)
	{
		BaseLambda = NewLambda;
	}

	void InitButterFlyStats(int64 amount)
	{
		vectorButterFly = vector<vector<double>>(amount, vector<double>(3, 0.0));
		initTmp = vector<bool>(amount, false);
	}

	void LoadMidiFromTxt() {
		if (!TXTFILE)
			return;
		int TEMP_NOTE[2];
		for (int i = 0; fscanf_s(TXTFILE, "%d %d", TEMP_NOTE, TEMP_NOTE + 1) != EOF;) {
			if (MIDI_INFO_PARTICLE.empty()) {
				vector<int64> TEM_NOTE;
				TEM_NOTE.push_back(*TEMP_NOTE);
				TEM_NOTE.push_back(*(TEMP_NOTE + 1));
				MIDI_INFO_PARTICLE.push_back(TEM_NOTE);
			}
			else {
				if (*TEMP_NOTE == MIDI_INFO_PARTICLE[i][0]) {
					MIDI_INFO_PARTICLE[i].push_back(*(TEMP_NOTE + 1));
				}
				else {
					++i;
					vector<int64> TEM_NOTE;
					TEM_NOTE.push_back(*TEMP_NOTE);
					TEM_NOTE.push_back(*(TEMP_NOTE + 1));
					MIDI_INFO_PARTICLE.push_back(TEM_NOTE);
				}
			}
		}
	}

	//模式A（每一次分支出的线条不聚合，且每一次分支的起点为读入的第一个音符的坐标）
	void AMode() const {//A
		for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
			for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
				if (i == 0) {
					StraightLine(*StartPos, MIDI_INFO_PARTICLE[i][0], YPos, YPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j]);//将此函数名与下一个else中函数名改为Particle.h中的任何一个绘制线条分组中非连续相切圆的函数可改变线条类型
				}
				else {
					StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][1], MIDI_INFO_PARTICLE[i][j]);
				}
			}
		}
	}

	//模式B（每一次分支出的线条不聚合，且每一次分支的起点为读入次序相同的音符，若本行音符数目大于上一行，则多出部分从读入的最后一个起始）
	void BMode() const {//B
		for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
			for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
				if (i == 0) {
					StraightLine(*StartPos, MIDI_INFO_PARTICLE[i][0], YPos, YPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j]);//将此函数名与下一个else中函数名改为Particle.h中的任何一个绘制线条分组中非连续相切圆的函数可改变线条类型
				}
				else {
					if (i > (int)MIDI_INFO_PARTICLE[i - 1].size()) {
						StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][(int)MIDI_INFO_PARTICLE[i - 1].size() - 1], MIDI_INFO_PARTICLE[i][j]);
					}
					else {
						StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][j]);
					}

				}
			}
		}
	}

	//模式C（每一次分支出的线条会聚合，聚合点读入的第一个音符的坐标）
	void CMode() const {//C
		for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
			for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
				if (i == 0) {
					StraightLine(*StartPos, MIDI_INFO_PARTICLE[i][0], YPos, YPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j]);//将此函数名与下一个else中函数名改为Particle.h中的任何一个绘制线条分组中非连续相切圆的函数可改变线条类型
				}
				else {
					StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][1], MIDI_INFO_PARTICLE[i][j]);
				}
			}
			if ((int)MIDI_INFO_PARTICLE[i - 1].size() > (int)MIDI_INFO_PARTICLE[i].size()) {
				for (int j = (int)MIDI_INFO_PARTICLE[i].size() - 1; j < (int)MIDI_INFO_PARTICLE[i - 1].size(); j++) {
					StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][1]);
				}
			}
		}
	}

	//模式D（每一次分支出的线条会聚合，聚合点规则类似于B模式）
	void DMode() const {//D
		for (int i = 0; i < (int)MIDI_INFO_PARTICLE.size(); i++) {
			for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); j++) {
				if (i == 0) {
					StraightLine(*StartPos, MIDI_INFO_PARTICLE[i][0], YPos, YPos, *(StartPos + 1), MIDI_INFO_PARTICLE[i][j]);//将此函数名与下一个else中函数名改为Particle.h中的任何一个绘制线条分组中非连续相切圆的函数可改变线条类型
				}
				else {
					if (i > (int)MIDI_INFO_PARTICLE[i - 1].size()) {
						StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][(int)MIDI_INFO_PARTICLE[i - 1].size() - 1], MIDI_INFO_PARTICLE[i][j]);
					}
					else {
						StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][j]);
					}
				}
			}
			if ((int)MIDI_INFO_PARTICLE[i - 1].size() > (int)MIDI_INFO_PARTICLE[i].size()) {
				for (int j = (int)MIDI_INFO_PARTICLE[i].size() - 1; j < (int)MIDI_INFO_PARTICLE[i - 1].size(); j++) {
					StraightLine(MIDI_INFO_PARTICLE[i - 1][0], MIDI_INFO_PARTICLE[i][0], YPos, YPos, MIDI_INFO_PARTICLE[i - 1][j], MIDI_INFO_PARTICLE[i][(int)MIDI_INFO_PARTICLE[i].size() - 1]);
				}
			}
		}
	}

	//海星线
	void StarLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) {
		const fp64 deltaT = abs(x2 - x1);
		const fp64 S = sqrt(square((x2 - x1)) + square((z2 - z1)));
		const fp64 theta = atan((z2 - z1) / deltaT);
		const fp64 H = (S / 6);
		const int64 N = (628 / (int64)deltaT) + 1;
		if (!Startype) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 3.14 \"s1,dis,cr,cg,cb=%.10lf+t,(-2-2^(sin(5*t)))*%.10lf,sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.005 %lld 20\n", int64(x1), (x1 + x2) / 2.0, y1 + 0.5, (z2 + z1) / 2.0, theta, H, N);
			Startype = true;
		}
		else {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 3.14 \"s1,dis,cr,cg,cb=%.10lf-t,(-2-2^(sin(5*t)))*%.10lf,sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.005 %lld 20\n", int64(x1), (x1 + x2) / 2.0, y1 + 0.5, (z2 + z1) / 2.0, theta, H, N);
			Startype = false;
		}
	}

	//直线
	void StraightLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 S = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 vy = ((y2 - y1) / deltaT);//参数方程（Y轴）
		const auto x = (int64)x1;//起始时间
		const auto NumBC = (fp64)((int64)(S / deltaT)) * 10.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10f %.10f %.10f 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, %.10f*t+0.5, %.10f*t+0.5, sin(t/7)/4+0.75, cos(t/5)/4+0.75, sin(t/3)/4+0.75\" %.10f %lld 25\n", x, x1, y1, z1, deltaT2, vy, vz, NumBC_, (int64)NumBC);
	}

	//笛卡儿叶形线
	void FoliumOfDescartes(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const fp64 S = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const fp64 theta = atan((z2 - z1) / deltaT);
		const fp64 H = S / 71.0;
		const fp64 A = S * 2.0 * sqrt(2.0) / 71.0;
		const int64 N = 30280 / (16 * (int64)deltaT) + 1;
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -0.7287 2.2997 \"s1,s2,dis,cr,cg,cb=%.10lf,t+PI/4,(3*%.10lf*sin(t)*cos(t))/(sin(t)*sin(t)*sin(t)+cos(t)*cos(t)*cos(t)),sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.0016 %lld 20\n", int64(x1), (x1 + x2) / 2.0 + 0.5, y1 + H + 0.5, (z2 + z1) / 2.0 + 0.5, theta, A, N);
	}

	//正弦曲线
	void SinStraightLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 S = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 vy = ((y2 - y1) / deltaT);//参数方程（Y轴）
		const auto x = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(S / deltaT)) * 15.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		const fp64 Sin_A_ = S / 10;//正弦曲线的振幅
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5,%.10lf*sin(t/%.10lf)+%.10lf*t+0.5, %.10lf*t+0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, Sin_A_, deltaT / (6 * pi), vy, vz, NumBC_, (int64)NumBC);
	}

	//抛物线
	void ParabolaLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 s = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const auto x = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		double H;//抛物线最高点
		const double deltaY = y2 - y1;//Y轴相对坐标
		const double deltaX = x2 - x1;//X轴相对坐标
		if (y2 > y1) 
			H = 0.2 * s + deltaY;//高度的定义
		else 
			H = 0.2 * s;//高度的定义
		double a;//y=bx(x-a)
		double b;//y=bx(x-a)
		if (!equal(deltaY, 0.0)) {
			a = ((4 * H * deltaX) + sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			if (a / 2 < 0 || a / 2 > x2)
				a = ((4 * H * deltaX) - sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			b = (-4 * H) / (a * a);
		}
		else{
			a = deltaX;
			b = (-4 * H) / (a * a);
		}//求y轴参数方程
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, (%.10lf*t*(t-%.10lf))+0.5, %.10lf*t+0.5, 1, 1, 1\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, b, a, vz, NumBC_, (int64)NumBC);
	}

	//螺线
	void SpiralLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		//求平均速度（直线移动速度）
		const fp64 Lambda = BaseLambda;
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 vx = ((x2 - x1) / deltaT);//参数方程（X为参数）
		const fp64 vy = ((y2 - y1) / deltaT);//参数方程（Y轴）
		const fp64 s = sqrt(((x2 - x1) * (x2 - x1)) + ((z2 - z1) * (z2 - z1)));//路程
		const fp64 h = y2 - y1;//高度
		const fp64 theta = atan((x2 - x1) / (z2 - z1));//水平偏转角
		const fp64 phi = atan(h / s);//竖直偏转角
		const fp64 omega = 2 * deltaT * pi / (Lambda * deltaT);//转速
		const auto x_1 = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		fp64 r_r = 3;
		if (s > 40)
			r_r = s / 8.0;
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = 0+(%.10lf * t - (%.10lf*sin(t*%.10lf)^2 * ((cos(%.10lf) * cos(%.10lf * t)) + (sin(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5), %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * (cos(%.10lf) * sin(%.10lf * t))) + 0.5, %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * ((sin(%.10lf) * cos(%.10lf * t)) - (cos(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x_1, x1, y1, z1, deltaT2, vx, r_r, pi / deltaT, theta, omega, theta, phi, omega, vy, r_r, pi / deltaT, phi, omega, vz, r_r, pi / deltaT, theta, omega, theta, phi, omega, NumBC_, (int64)NumBC);

	}

	//螺线(带直线）
	void SpiralLineStraight(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		//求平均速度（直线移动速度）
		const fp64 Lambda = BaseLambda;
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 vx = ((x2 - x1) / deltaT);//参数方程（X为参数）
		const fp64 vy = ((y2 - y1) / deltaT);//参数方程（Y轴）
		const fp64 s = sqrt(((x2 - x1) * (x2 - x1)) + ((z2 - z1) * (z2 - z1)));//路程
		const fp64 h = y2 - y1;//高度
		const fp64 theta = atan((x2 - x1) / (z2 - z1));//水平偏转角
		const fp64 phi = atan(h / s);//竖直偏转角
		const fp64 omega = 2 * deltaT * pi / (Lambda * deltaT);//转速
		const auto x_1 = (int64)x1;//起始时间
		const auto x = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		fp64 r_r = 3;
		if (s > 40)
			r_r = s / 8.0;
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = 0+(%.10lf * t - (%.10lf*sin(t*%.10lf)^2 * ((cos(%.10lf) * cos(%.10lf * t)) + (sin(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5), %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * (cos(%.10lf) * sin(%.10lf * t))) + 0.5, %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * ((sin(%.10lf) * cos(%.10lf * t)) - (cos(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x_1, x1, y1, z1, deltaT2, vx, r_r, pi / deltaT, theta, omega, theta, phi, omega, vy, r_r, pi / deltaT, phi, omega, vz, r_r, pi / deltaT, theta, omega, theta, phi, omega, NumBC_, (int64)NumBC);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, %.10lf*t+0.5, %.10lf*t+0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, vy, vz, NumBC_, (int64)NumBC);
	}

	//螺线组
	void SpiralLines(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		//求平均速度（直线移动速度）
		const fp64 Lambda = BaseLambda;
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 vx = ((x2 - x1) / deltaT);//参数方程（X为参数）
		const fp64 vy = ((y2 - y1) / deltaT);//参数方程（Y轴）
		const fp64 s = sqrt(((x2 - x1) * (x2 - x1)) + ((z2 - z1) * (z2 - z1)));//路程
		const fp64 h = y2 - y1;//高度
		const fp64 theta = atan((x2 - x1) / (z2 - z1));//水平偏转角
		const fp64 phi = atan(h / s);//竖直偏转角
		const fp64 omega = 4 * deltaT * pi / (Lambda * s);//转速
		const auto x_1 = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		fp64 r_r = 4;
		if (s > 40)
			r_r = s / 8.0;
		//螺线最大半径
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = 0+(%.10lf * t - (%.10lf*sin(t*%.10lf)^2 * ((cos(%.10lf) * cos(%.10lf * t)) + (sin(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5), %.10lf * t + (%.10lf*sin(t*%.10lf) * (cos(%.10lf) * sin(%.10lf * t))) + 0.5, %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * ((sin(%.10lf) * cos(%.10lf * t)) - (cos(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x_1, x1, y1, z1, deltaT2, vx, r_r, pi / deltaT, theta, omega, theta, phi, omega, vy, r_r, pi / deltaT, phi, omega, vz, r_r, pi / deltaT, theta, omega, theta, phi, omega, NumBC_, (int64)NumBC);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = 0+(%.10lf * t - (%.10lf*sin(t*%.10lf)^2 * ((cos(%.10lf) * cos(%.10lf * t + (3.1415926535))) + (sin(%.10lf) * sin(%.10lf) * sin(%.10lf * t+ (3.1415926535))))) + 0.5), %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * (cos(%.10lf) * sin(%.10lf * t+ (3.1415926535)))) + 0.5, %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * ((sin(%.10lf) * cos(%.10lf * t+ (3.1415926535))) - (cos(%.10lf) * sin(%.10lf) * sin(%.10lf * t+ (3.1415926535))))) + 0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x_1, x1, y1, z1, deltaT2, vx, r_r, pi / deltaT, theta, omega, theta, phi, omega, vy, r_r, pi / deltaT, phi, omega, vz, r_r, pi / deltaT, theta, omega, theta, phi, omega, NumBC_, (int64)NumBC);

	}

	//螺线组(带中心直线）
	void SpiralLinesWithStraightLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		//求平均速度（直线移动速度）
		const fp64 Lambda = BaseLambda;
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 vx = ((x2 - x1) / deltaT);//参数方程（X为参数）
		const fp64 vy = ((y2 - y1) / deltaT);//参数方程（Y轴）
		const fp64 s = sqrt(((x2 - x1) * (x2 - x1)) + ((z2 - z1) * (z2 - z1)));//路程
		const fp64 h = y2 - y1;//高度
		const fp64 theta = atan((x2 - x1) / (z2 - z1));//水平偏转角
		const fp64 phi = atan(h / s);//竖直偏转角
		const fp64 omega = 4 * deltaT * pi / (Lambda * deltaT);//转速
		const auto x_1 = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		fp64 r_r = 3;
		if (s > 40) {
			r_r = s / 8.0;
		}//螺线最大半径
		if (r_r > 5) {
			r_r = 5;
		}
		const auto x = (int64)x1;//起始时间
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = 0+(%.10lf * t - (%.10lf*sin(t*%.10lf)^2 * ((cos(%.10lf) * cos(%.10lf * t)) + (sin(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5), %.10lf * t + (%.10lf*sin(t*%.10lf) * (cos(%.10lf) * sin(%.10lf * t))) + 0.5, %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * ((sin(%.10lf) * cos(%.10lf * t)) - (cos(%.10lf) * sin(%.10lf) * sin(%.10lf * t)))) + 0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x_1, x1, y1, z1, deltaT2, vx, r_r, pi / deltaT, theta, omega, theta, phi, omega, vy, r_r, pi / deltaT, phi, omega, vz, r_r, pi / deltaT, theta, omega, theta, phi, omega, NumBC_, (int64)NumBC);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = 0+(%.10lf * t - (%.10lf*sin(t*%.10lf)^2 * ((cos(%.10lf) * cos(%.10lf * t + (3.1415926535))) + (sin(%.10lf) * sin(%.10lf) * sin(%.10lf * t+ (3.1415926535))))) + 0.5), %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * (cos(%.10lf) * sin(%.10lf * t+ (3.1415926535)))) + 0.5, %.10lf * t + (%.10lf*sin(t*%.10lf)^2 * ((sin(%.10lf) * cos(%.10lf * t+ (3.1415926535))) - (cos(%.10lf) * sin(%.10lf) * sin(%.10lf * t+ (3.1415926535))))) + 0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x_1, x1, y1, z1, deltaT2, vx, r_r, pi / deltaT, theta, omega, theta, phi, omega, vy, r_r, pi / deltaT, phi, omega, vz, r_r, pi / deltaT, theta, omega, theta, phi, omega, NumBC_, (int64)NumBC);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, %.10lf*t+0.5, %.10lf*t+0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, vy, vz, NumBC_, (int64)NumBC);
	}

	//抛物线（下方向）
	void ReverseParabolaLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 s = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const auto x = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		double H;//抛物线最高点
		const double deltaY = y2 - y1;//Y轴相对坐标
		const double deltaX = x2 - x1;//X轴相对坐标
		if (y2 < y1)
			H = -0.2 * s + deltaY;
		else
			H = -0.2 * s;
		double a;//y=bx(x-a)
		double b;//y=bx(x-a)
		if (!equal(deltaY, 0.0)) {
			a = ((4 * H * deltaX) + sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			if (a / 2 < 0 || a / 2 > x2) {
				a = ((4 * H * deltaX) - sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			}
			b = (-4 * H) / (a * a);
		}
		else{
			a = deltaX;
			b = (-4 * H) / (a * a);
		}//求y轴参数方程
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, (%.10lf*t*(t-%.10lf))+0.5, %.10lf*t+0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, b, a, vz, NumBC_, (int64)NumBC);
	}

	//（有螺旋圈阵）的抛物线
	void SpiralRoundParabolaLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 s = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const auto x = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		double H;//抛物线最高点
		const double deltaY = y2 - y1;//Y轴相对坐标
		const double deltaX = x2 - x1;//X轴相对坐标
		if (y2 > y1) 
			H = 0.2 * s + deltaY;//高度的定义
		else 
			H = 0.2 * s;//高度的定义
		double a;//y=bx(x-a)
		double b;//y=bx(x-a)
		if (!equal(deltaY, 0.0)) {
			a = ((4 * H * deltaX) + sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			if (a / 2 < 0 || a / 2 > x2) {
				a = ((4 * H * deltaX) - sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			}
			b = (-4 * H) / (a * a);
		}
		else{
			a = deltaX;
			b = (-4 * H) / (a * a);
		}//求y轴参数方程
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, (%.10lf*t*(t-%.10lf))+0.5, %.10lf*t+0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, b, a, vz, NumBC_, (int64)NumBC);
	}

	//（有螺旋圈阵）的抛物线（反向）
	void ReverseSpiralRoundParabolaLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		const fp64 deltaT = abs(x2 - x1);//求时间间隔（Tick）
		const auto deltaT2 = (int64)abs(x2 - x1);//转换为整形
		const fp64 vz = ((z2 - z1) / deltaT);//参数方程（Z轴）
		const fp64 s = sqrt(square((x2 - x1)) + square((z2 - z1)));//求位移
		const auto x = (int64)x1;//起始时间
		const fp64 NumBC = (fp64)((int64)(s / deltaT)) * 20.0; //每一刻计算的次数
		const fp64 NumBC_ = 1.0000000 / NumBC;//每次计算的参数递增量
		double H;//抛物线最高点
		const double deltaY = y2 - y1;//Y轴相对坐标
		const double deltaX = x2 - x1;//X轴相对坐标
		if (y2 < y1) 
			H = -0.2 * s + deltaY;
		else 
			H = -0.2 * s;
		double a;//y=bx(x-a)
		double b;//y=bx(x-a)
		if (!equal(deltaY, 0.0)) {
			a = ((4 * H * deltaX) + sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			if (a / 2 < 0 || a / 2 > x2) {
				a = ((4 * H * deltaX) - sqrt((16 * H * H * deltaX * deltaX) - 16 * H * deltaY * deltaX * deltaX)) / (2 * deltaY);
			}
			b = (-4 * H) / (a * a);
		}
		else{
			a = deltaX;
			b = (-4 * H) / (a * a);
		}//求y轴参数方程
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %lld \"x, y, z, cr, cg, cb = t+0.5, (%.10lf*t*(t-%.10lf))+0.5, %.10lf*t+0.5, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" %.10lf %lld 25\n", x, x1, y1, z1, deltaT2, b, a, vz, NumBC_, (int64)NumBC);
	}
	
	void Butter_Fly_To_Next_Point(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2, int64 track) {
		fp64 x = x1 + 0.5, y = y1 + 0.5, z = z1 + 0.5;
		x2 += 0.5;
		x1 += 0.5;
		y1 += 0.5;
		y2 += 0.5;
		z1 += 0.5;
		z2 += 0.5;
		const int64 a = (int64)round(x2 - x1) * 4;
		fp64 wing = 0.6;
		fp64 wingAddValue = 0.0;
		vector<double> vectorButterFlyTarget = { (x2 - x1) / (fp64)a, (y2 - y1) / (fp64)a,(z2 - z1) / (fp64)a };
		if (!initTmp[track]) {
			vectorButterFly[track] = { (x2 - x1) / (fp64)a, (y2 - y1) / (fp64)a,(z2 - z1) / (fp64)a };
			for (int64 i = 0; i < a; i++) {
				if (i < (a / 2) - 1) {
					wing = 0.75 * cos((fp64)i * pi / 10.0);
				}
				else if (i == (a / 2) - 1) {
					wingAddValue = (0.75 - wing) / ((fp64)a - (((fp64)a / 2) - 1));
					wing += wingAddValue;
				}
				else {
					wing += wingAddValue;
				}
				const fp64 s2 = atan(vectorButterFly[track][1] / vectorButterFly[track][2]);
				const fp64 s1 = atan(vectorButterFly[track][2] / vectorButterFly[track][0]);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 6.29 \"s1,dis,cr,cg,cb,s2=%.10lf+t,(E^cos(t)-2*cos(4*t)+(sin(t/12))^5)*%.10lf,0.25+0.25*sin(t),0.3+0.3*cos(t),0.6+0.1*sin(t),%.10lf+abs(%.10lf*sin(t))\" 0.01 1\n", (int64)(x1 * 4) + i, x, y + (log((x2 - x1) / 1.5) - log((x2 - x1) / 1.5) * cos((pi * 2.0 * (fp64)i / (fp64)a))) / 2.0, z, s1, ButterFlySize, s2, wing);
				x += vectorButterFly[track][0];
				y += vectorButterFly[track][1];
				z += vectorButterFly[track][2];
			}
			initTmp[track] = true;
		}
		else {
			const int64 j = a / 4;
			for (int64 i = 0; i < a; i++) {
				if (i < (a / 2) - 1) {
					wing = 0.75 * cos((fp64)i * pi / 10.0);
				}
				else if (i == (a / 2) - 1) {
					wingAddValue = (0.75 - wing) / ((fp64)a - (((fp64)a / 2) - 1));
					wing += wingAddValue;
				}
				else {
					wing += wingAddValue;
				}
				if (i < j) {
					vectorButterFlyTarget = { (x2 - x) / (fp64)a, (y2 - y) / (fp64)a,(z2 - z) / (fp64)a };
					vectorButterFly[track][0] += (vectorButterFlyTarget[0] - vectorButterFly[track][0]) / (fp64)((fp64)j - (fp64)i);
					vectorButterFly[track][1] += (vectorButterFlyTarget[1] - vectorButterFly[track][1]) / (fp64)((fp64)j - (fp64)i);
					vectorButterFly[track][2] += (vectorButterFlyTarget[2] - vectorButterFly[track][2]) / (fp64)((fp64)j - (fp64)i);
				}
				else if (i == j) {
					vectorButterFly[track] = { (x2 - x) / ((fp64)a - (fp64)j), (y2 - y) / ((fp64)a - (fp64)j),(z2 - z) / ((fp64)a - (fp64)j) };
				}
				const fp64 s1 = atan(vectorButterFly[track][2] / vectorButterFly[track][0]);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 6.29 \"s1,dis,cr,cg,cb,s2=%.10lf+t,(E^cos(t)-2*cos(4*t)+(sin(t/12))^5)*%.10lf,0.25+0.25*sin(t),0.3+0.3*cos(t),0.6+0.1*sin(t),%.10lf+abs(%.10lf*sin(t))\" 0.01 1\n", (int64)(x1 * 4) + i, x, y + (log((x2 - x1) / 1.5) - log((x2 - x1) / 1.5) * cos((pi * 2.0 * (fp64)i / (fp64)a))) / 2.0, z, s1, ButterFlySize, 0.0, wing);
				x += vectorButterFly[track][0];
				y += vectorButterFly[track][1];
				z += vectorButterFly[track][2];
			}
		}

	}

	void Butter_Fly_To_Next_Point_Sim(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2, int64 track) {
		fp64 x = x1 + 0.5, y = y1 + 0.5, z = z1 + 0.5;
		x2 += 0.5;
		x1 += 0.5;
		y1 += 0.5;
		y2 += 0.5;
		z1 += 0.5;
		z2 += 0.5;
		const int64 a = (int64)round(x2 - x1) * 4;
		fp64 wing = 0.6;
		fp64 wingAddValue = 0.0;
		vector<double> vectorButterFlyTarget = { (x2 - x1) / (fp64)a, (y2 - y1) / (fp64)a,(z2 - z1) / (fp64)a };
		if (!initTmp[track]) {
			vectorButterFly[track] = { (x2 - x1) / (fp64)a, (y2 - y1) / (fp64)a,(z2 - z1) / (fp64)a };
			for (int64 i = 0; i < a; i++) {
				if (i < (a / 2) - 1) {
					wing = 0.75 * cos((fp64)i * pi / 40.0);
				}
				else if (i == (a / 2) - 1) {
					wingAddValue = (0.75 - wing) / ((fp64)a - (((fp64)a / 2) - 1));
					wing += wingAddValue;
				}
				else {
					wing += wingAddValue;
				}
				const fp64 s2 = atan(vectorButterFly[track][1] / vectorButterFly[track][2]);
				const fp64 s1 = atan(vectorButterFly[track][2] / vectorButterFly[track][0]);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 6.29 \"s1,dis,cr,cg,cb,s2=%.10lf+t,(E^cos(t)-2*cos(4*t)+(sin(t/12))^5)*%.10lf,0.25+0.25*sin(t),0.3+0.3*cos(t),0.6+0.1*sin(t),%.10lf+abs(%.10lf*sin(t))\" 0.01 1\n", (int64)x1 + i, x, y + (log((x2 - x1) / 1.5) - log((x2 - x1) / 1.5) * cos((pi * 2.0 * (fp64)i / (fp64)a))) / 2.0, z, s1, ButterFlySize, s2, wing);
				x += vectorButterFly[track][0];
				y += vectorButterFly[track][1];
				z += vectorButterFly[track][2];
			}
			initTmp[track] = true;
		}
		else {
			const int64 j = a / 4;
			for (int64 i = 0; i < a; i++) {
				if (i < (a / 2) - 1) {
					wing = 0.75 * cos((fp64)i * pi / 10.0);
				}
				else if (i == (a / 2) - 1) {
					wingAddValue = (0.75 - wing) / ((fp64)a - (((fp64)a / 2) - 1));
					wing += wingAddValue;
				}
				else {
					wing += wingAddValue;
				}
				if (i < j) {
					vectorButterFlyTarget = { (x2 - x) / (fp64)a, (y2 - y) / (fp64)a,(z2 - z) / (fp64)a };
					vectorButterFly[track][0] += (vectorButterFlyTarget[0] - vectorButterFly[track][0]) / (fp64)((fp64)j - (fp64)i);
					vectorButterFly[track][1] += (vectorButterFlyTarget[1] - vectorButterFly[track][1]) / (fp64)((fp64)j - (fp64)i);
					vectorButterFly[track][2] += (vectorButterFlyTarget[2] - vectorButterFly[track][2]) / (fp64)((fp64)j - (fp64)i);
				}
				else if (i == j) {
					vectorButterFly[track] = { (x2 - x) / ((fp64)a - (fp64)j), (y2 - y) / ((fp64)a - (fp64)j),(z2 - z) / ((fp64)a - (fp64)j) };
				}
				const fp64 s1 = atan(vectorButterFly[track][2] / vectorButterFly[track][0]);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 6.29 \"s1,dis,cr,cg,cb,s2=%.10lf+t,(E^cos(t)-2*cos(4*t)+(sin(t/12))^5)*%.10lf,0.25+0.25*sin(t),0.3+0.3*cos(t),0.6+0.1*sin(t),%.10lf+abs(%.10lf*sin(t))\" 0.01 1\n", (int64)x1 + i, x, y + (log((x2 - x1) / 1.5) - log((x2 - x1) / 1.5) * cos((pi * 2.0 * (fp64)i / (fp64)a))) / 2.0, z, s1, ButterFlySize, 0.0, wing);
				x += vectorButterFly[track][0];
				y += vectorButterFly[track][1];
				z += vectorButterFly[track][2];
			}
		}
	}

	//上下抛物线中心螺线链接
	void SpecialParabolaLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		ParabolaLine(x1, x2, y1, y2, z1, z2);
		ReverseParabolaLine(x1, x2, y1 - 20, y2 - 20, z1, z2);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"x, y, z, cr, cg, cb = cos(t * (8 * PI / 10)), t, sin(t * (8 * PI / 10)), sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 - 20 + 0.5, z2 + 0.5);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"x, y, z, cr, cg, cb = cos(t * (8 * PI / 10)), 0-t, sin(t * (8 * PI / 10)), sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 + 0.5, z2 + 0.5);
	}

	//上下抛物线中心螺线链接（有螺旋圈阵）
	void special_parabola_line_round(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		SpiralRoundParabolaLine(x1, x2, y1, y2, z1, z2);
		ReverseSpiralRoundParabolaLine(x1, x2, y1 - 20, y2 - 20, z1, z2);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"x, y, z, cr, cg, cb = cos(t * (8 * PI / 10)), t, sin(t * (8 * PI / 10)), sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 - 20 + 0.5, z2 + 0.5);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"x, y, z, cr, cg, cb = cos(t * (8 * PI / 10)), 0-t, sin(t * (8 * PI / 10)), sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 + 0.5, z2 + 0.5);
	}

	//上下直线中心螺线链接
	void SpecialStraightLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		StraightLine(x1, x2, y1, y2, z1, z2);
		StraightLine(x1, x2, y1 - 20, y2 - 20, z1, z2);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"x, y, z, cr, cg, cb = cos(t * (8 * PI / 10)), t, sin(t * (8 * PI / 10)), sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 - 20 + 0.5, z2 + 0.5);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"x, y, z, cr, cg, cb = cos(t * (8 * PI / 10)), 0-t, sin(t * (8 * PI / 10)), sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 + 0.5, z2 + 0.5);
	}

	//上下抛物线中心直线链接
	void special_parabola_line_s(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		ParabolaLine(x1, x2, y1, y2, z1, z2);
		ReverseParabolaLine(x1, x2, y1 - 20, y2 - 20, z1, z2);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"y, cr, cg, cb = t, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 - 20 + 0.5, z2 + 0.5);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"y, cr, cg, cb = 0-t, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 + 0.5, z2 + 0.5);
	}

	//上下抛物线中心直线链接（有螺旋圈阵）
	void SpecialParabolaLinesRound(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		SpiralRoundParabolaLine(x1, x2, y1, y2, z1, z2);
		ReverseSpiralRoundParabolaLine(x1, x2, y1 - 20, y2 - 20, z1, z2);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"y, cr, cg, cb = t, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 - 20 + 0.5, z2 + 0.5);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"y, cr, cg, cb = 0-t, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 + 0.5, z2 + 0.5);
	}

	//上下直线中心直线链接
	void SpecialStraightLines(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2) const {
		StraightLine(x1, x2, y1, y2, z1, z2);
		StraightLine(x1, x2, y1 - 20, y2 - 20, z1, z2);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"y, cr, cg, cb = t, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 - 20 + 0.5, z2 + 0.5);
		fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 10 \"y, cr, cg, cb = 0-t, sin(t/7)/2+0.5, cos(t/5)/2+0.5, sin(t/3)/2+0.5\" 0.1 20 25\n", int64(x2), x2 + 0.5, y2 + 0.5, z2 + 0.5);
	}

	int64 StartPos[2] = { 0,0 };
	fp64 ButterFlySize = 0.5;
	fp64 YPos = 100.0;
private:
	static bool equal(fp64 a, fp64 b, fp64 n = 5)
	{
		return ((a - b) < pow(10, -n) && (a - b) > -pow(10, -n));
	}
	static double square(double input) {
		return input * input;
	}
	string ParticleType = "minecraft:end_rod";
	bool Startype = false;
	FILE* OutPutFile = nullptr;
	fp64 BaseLambda = 20.0;
	const double pi = 3.1415926535897932384626433832795;
	vector<vector<double>> vectorButterFly;
	vector<bool> initTmp;
	FILE* TXTFILE = nullptr;
	vector<vector<int64>> MIDI_INFO_PARTICLE;
};

class CoiledTangentRound
{
public:
	using fp64 = double;
	using fp32 = float;
	using int64 = long long;
	using uint64 = unsigned long long;
	CoiledTangentRound() = default;
	CoiledTangentRound(const string& pty) : ParticleType(pty) {}
	CoiledTangentRound(string&& pty) : ParticleType(pty) {}

	void open(FILE* file)
	{
		OutPutFile = file;
		base.open(file);
	}

	void ChangeBaseLambda(fp64 NewLambda)
	{
		BaseLambda = NewLambda;
	}

	void InitTrack(int64 TrackAmount)
	{
		if (TrackAmount < 1)
			return;
		LastCenter = vector<Point>(TrackAmount, Point());
		Init = vector<int64>(TrackAmount, 1);
		Init[0] = 0;
		Mode = vector<int64>(TrackAmount, 1);
	}

	void ChangeLambda(fp64 newLambda) { base.ChangeBaseLambda(newLambda); }

	void LoadMidiFromTxt() {
		if (!TXTFILE)
			return;
		int TEMP_NOTE[2];
		for (int i = 0; fscanf_s(TXTFILE, "%d %d", TEMP_NOTE, TEMP_NOTE + 1) != EOF;) {
			if (MIDI_INFO_PARTICLE.empty()) {
				vector<int64> TEM_NOTE;
				TEM_NOTE.push_back(*TEMP_NOTE);
				TEM_NOTE.push_back(*(TEMP_NOTE + 1));
				MIDI_INFO_PARTICLE.push_back(TEM_NOTE);
			}
			else {
				if (*TEMP_NOTE == MIDI_INFO_PARTICLE[i][0]) {
					MIDI_INFO_PARTICLE[i].push_back(*(TEMP_NOTE + 1));
				}
				else {
					++i;
					vector<int64> TEM_NOTE;
					TEM_NOTE.push_back(*TEMP_NOTE);
					TEM_NOTE.push_back(*(TEMP_NOTE + 1));
					MIDI_INFO_PARTICLE.push_back(TEM_NOTE);
				}
			}
		}
	}

	void Round() {//ROUND
		for (size_t i = 1; i < MIDI_INFO_PARTICLE.size(); ++i) {
			if ((int)MIDI_INFO_PARTICLE[i].size() > (int)MIDI_INFO_PARTICLE[i - 1].size()) {
				InitReset(MIDI_INFO_PARTICLE[i - 1].size());
			}
			for (int j = 1; j < (int)MIDI_INFO_PARTICLE[i].size(); ++j) {
				if (j < (int)MIDI_INFO_PARTICLE[i - 1].size()) {
					RoundLine(static_cast<fp64>(MIDI_INFO_PARTICLE[i - 1][0]), static_cast<fp64>(MIDI_INFO_PARTICLE[i][0]), YPos, YPos, static_cast<fp64>(MIDI_INFO_PARTICLE[i - 1][j]), static_cast<fp64>(MIDI_INFO_PARTICLE[i][j]), j);
				}
				else {
					RoundLine(static_cast<fp64>(MIDI_INFO_PARTICLE[i - 1][0]), static_cast<fp64>(MIDI_INFO_PARTICLE[i][0]), YPos, YPos, static_cast<fp64>(MIDI_INFO_PARTICLE[i - 1][MIDI_INFO_PARTICLE[i - 1].size() - 1]), static_cast<fp64>(MIDI_INFO_PARTICLE[i][j]), j);
				}
			}
		}
	}
	fp64 YPos = 100.0;
private:
	struct Note
	{
		int64 time = 0;
		int64 pitch = 0;
		Note(int64 t, int64 p) :time(t), pitch(p) {}
	};
	static bool equal(fp64 a, fp64 b, fp64 n = 5)
	{
		return ((a - b) < pow(10, -n) && (a - b) > -pow(10, -n));
	}
	static fp64 square(fp64 input) {
		return input * input;
	}
	struct Line
	{
		fp64 k;
		fp64 b;
		bool Test;
		Line() :k(0.0), b(0.0), Test(false) {}
	};
	struct Point
	{
		fp64 z;
		fp64 x;
		Point() :z(0.0), x(0.0) {}
	};
	static void GetLine(fp64 z0, fp64 z1, fp64 x0, fp64 x1, Line& out_Line)
	{
		const fp64 m = z1 - z0;
		if (equal(m, 0.0))
		{
			out_Line.k = 1000000;
			out_Line.b = x0 - (1000000 * z0);
			out_Line.Test = false;
		}
		else
		{
			out_Line.k = (x1 - x0) / (z1 - z0);
			out_Line.b = x0 - (out_Line.k * z0);
			out_Line.Test = true;
		}
	}
	static void GetMidLine(const Line& Line_in, Line& Line_out, fp64 z0, fp64 z1, fp64 x0, fp64 x1)
	{
		if (!Line_in.Test)
		{
			Line_out.k = 0;
			Line_out.b = (x1 + x0) * 0.5;
			Line_out.Test = true;
		}
		else if (!equal(Line_in.k, 0.0))
		{
			Line_out.k = -1 / Line_in.k;
			Line_out.b = ((x1 + x0) * 0.5) - (Line_out.k * ((z1 + z0) * 0.5));
			Line_out.Test = true;
		}
		else
		{
			Line_out.k = 1000000;
			Line_out.b = ((x1 + x0) * 0.5) - (Line_out.k * ((z1 + z0) * 0.5));
			Line_out.Test = false;
		}
	}
	static void GetCross(const Line& L1, const Line& L2, Point& Cross, fp64 Z)
	{
		if ((!L1.Test) && equal(L2.k, 0.0))
		{
			Cross.z = Z;
			Cross.x = L2.b;
		}
		else if (!L1.Test)
		{
			Cross.z = Z;
			Cross.x = (L2.k * Z) + L2.b;
		}
		else
		{
			Cross.z = ((L2.b - L1.b) / (L1.k - L2.k));
			Cross.x = (L1.k * Cross.z) + L1.b;
		}
	}
	void RoundLine(fp64 x1, fp64 x2, fp64 y1, fp64 y2, fp64 z1, fp64 z2, int64 track)
	{
		fp64 R;
		Point Center;
		Line L1;
		Line L2;
		Line L3;
		if (Init[track] == 0)
		{
			GetLine(z1, z2, x1, x2, L1);
			GetMidLine(L1, L2, z1, z2, x1, x2);
			L3.k = 0;
			L3.b = x1;
			if (L2.k >= -0.000001 && L2.k <= 0.000001)
			{
				base.ParabolaLine(x1, x2, y1, y2, z1, z2);
				return;
			}
			GetCross(L3, L2, Center, LastCenter[track].z);
			LastCenter[track] = Center;
			Init[track] = 1;
			if (z2 > z1)
				Mode[track] = 1;
			else
				Mode[track] = -1;
			const fp64 R1 = sqrt(((Center.x - x1) * (Center.x - x1)) + (((Center.z - z1) * (Center.z - z1))));
			const fp64 R2 = sqrt(((Center.x - x2) * (Center.x - x2)) + (((Center.z - z2) * (Center.z - z2))));
			R = (R1 + R2) * 0.5;
			if (R > 200) {
				base.ParabolaLine(x1, x2, y1, y2, z1, z2);
				return;
			}
		}
		else
		{
			GetLine(z1, z2, x1, x2, L1);
			GetMidLine(L1, L2, z1, z2, x1, x2);
			GetLine(LastCenter[track].z, z1, LastCenter[track].x, x1, L3);
			const fp64 test_line = L3.k - L2.k;
			if (test_line <= 0.000001 && test_line >= 0.000001)
			{
				base.ParabolaLine(x1, x2, y1, y2, z1, z2);
				Init[track] = 0;
				return;
			}
			GetCross(L3, L2, Center, LastCenter[track].z);
			R = sqrt(((Center.x - x1) * (Center.x - x1)) + (((Center.z - z1) * (Center.z - z1))));
			if (Center.z >= LastCenter[track].z)
			{
				if (z1 <= Center.z && z1 >= LastCenter[track].z)
					Mode[track] = -1 * Mode[track];
			}
			else
			{
				if (z1 >= Center.z && z1 <= LastCenter[track].z)
					Mode[track] = -1 * Mode[track];
			}
			LastCenter[track] = Center;
			if (R > 200) {
				base.ParabolaLine(x1, x2, y1, y2, z1, z2);
				return;
			}
		}
		const fp64 C = 2 * pi * R;
		const fp64 deltaT = abs(x2 - x1) / SPEED;
		fp64 Lambda = (fp64)((int64)((C * 20 / deltaT)));
		if (Lambda < 40)
			Lambda = 40;
		const fp64 S = sqrt(((x2 - x1) * (x2 - x1)) + ((z2 - z1) * (z2 - z1)));
		const fp64 Theta = acos(S / (2 * R));
		const fp64 Phi = pi - (2 * Theta);
		const fp64 TTT = ((Center.x - L1.b) / L1.k);
		fp64 Omega = 0.0;
		if (Mode[track] == 1)
		{
			if (TTT <= Center.z)
				Omega = Phi / (deltaT * Lambda);
			else
				Omega = ((2 * pi) - Phi) / (deltaT * Lambda);
		}
		if (Mode[track] == -1)
		{
			if (TTT >= Center.z)
				Omega = Phi / (deltaT * Lambda);
			else
				Omega = ((2 * pi) - Phi) / (deltaT * Lambda);
		}
		fp64 Omega_0;
		if (Init[track] != 0) {
			Omega_0 = atan(L3.k);
			if (z1 < Center.z)
				Omega_0 = Omega_0 + pi;
		}
		else {
			if (z1 < Center.z)
				Omega_0 = pi;
			else
				Omega_0 = 0;
		}
		const auto deltaT2 = (int64)x1;//转换为整形
		if(KH)
		{
			if (Mode[track] == -1) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %.10lf \"x, y, z, cr, cg, cb = %.10lf * sin(%.10lf + %.10lf * t)+0.5, %.10lf, %.10lf * cos(%.10lf + %.10lf * t)+0.5, cos(2*PI*t/(%.10lf)) / 4 + 0.75, cos(2*PI*t/(%.10lf)) / 4 + 0.75, sin(2*PI*t/(%.10lf)) / 4 + 0.75\" 1 %lld 25 \"(vx, vz) = ((random(), random()) - 0.5) * t / 5\" %.10lf\n", deltaT2, Center.x, y1 + 0.5, Center.z, deltaT * Lambda, R, Omega_0, Omega, 0.0, R, Omega_0, Omega, deltaT * Lambda, deltaT * Lambda, deltaT * Lambda, (int64)Lambda, 0.5);
			}
			else {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %.10lf \"x, y, z, cr, cg, cb = %.10lf * sin(%.10lf - %.10lf * t)+0.5, %.10lf, %.10lf * cos(%.10lf - %.10lf * t)+0.5, cos(2*PI*t/(%.10lf)) / 4 + 0.75, cos(2*PI*t/(%.10lf)) / 4 + 0.75, sin(2*PI*t/(%.10lf)) / 4 + 0.75\" 1 %lld 25 \"(vx, vz) = ((random(), random()) - 0.5) * t / 5\" %.10lf\n", deltaT2, Center.x, y1 + 0.5, Center.z, deltaT * Lambda, R, Omega_0, Omega, 0.0, R, Omega_0, Omega, deltaT * Lambda, deltaT * Lambda, deltaT * Lambda, (int64)Lambda, 0.5);
			}
		}
		else
		{
			if (Mode[track] == -1) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %.10lf \"x, y, z, cr, cg, cb = %.10lf * sin(%.10lf + %.10lf * t)+0.5, %.10lf, %.10lf * cos(%.10lf + %.10lf * t)+0.5, cos(2*PI*t/(%.10lf)) / 4 + 0.75, cos(2*PI*t/(%.10lf)) / 4 + 0.75, sin(2*PI*t/(%.10lf)) / 4 + 0.75\" 1 %lld 25\n", deltaT2, Center.x, y1 + 0.5, Center.z, deltaT * Lambda, R, Omega_0, Omega, y1 + 0.5, R, Omega_0, Omega, deltaT * Lambda, deltaT * Lambda, deltaT * Lambda, (int64)Lambda);
			}
			else {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %.10lf \"x, y, z, cr, cg, cb = %.10lf * sin(%.10lf - %.10lf * t)+0.5, %.10lf, %.10lf * cos(%.10lf - %.10lf * t)+0.5, cos(2*PI*t/(%.10lf)) / 4 + 0.75, cos(2*PI*t/(%.10lf)) / 4 + 0.75, sin(2*PI*t/(%.10lf)) / 4 + 0.75\" 1 %lld 25\n", deltaT2, Center.x, y1 + 0.5, Center.z, deltaT * Lambda, R, Omega_0, Omega, y1 + 0.5, R, Omega_0, Omega, deltaT * Lambda, deltaT * Lambda, deltaT * Lambda, (int64)Lambda);
			}
		}
		Init[track] = 1;
	}
	void InitReset(uint64 TruckAmount) {
		for (uint64 i = TruckAmount - 1; i < Init.size(); ++i) {
			LastCenter[i] = LastCenter[i - 1];
			Mode[i] = Mode[i - 1];
		}
	}
	string ParticleType = "minecraft:end_rod";
	FILE* OutPutFile = nullptr;
	fp64 BaseLambda = 20.0;
	const double pi = 3.1415926535897932384626433832795;
	vector<Point> LastCenter;
	vector<int64> Init;
	vector<int64> Mode;
	BaseLines base;
	fp64 SPEED = 1;
	bool KH = false;
	FILE* TXTFILE = nullptr;
	vector<vector<int64>> MIDI_INFO_PARTICLE;
};
MinecraftParticleMiDiPlayerLinesEnd

namespace MEndPointParticle
{
	class EParticle
	{
	public:
		using int64 = long long;
		using fp64 = double;
		void print_round_move(long double ex, long double ey, long double ez) {
			int64 x_3 = (int64)ex;
			long double COR;//RGB颜色通道
			long double COG;//RGB颜色通道
			long double COB;//RGB颜色通道

			COR = 0.3;//随机RGB颜色通道
			COG = 0.7;//随机RGB颜色通道
			COB = 0.4;//随机RGB颜色通道

			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex tickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf 1 0 0 0 0 2 \"s1,dis=t*10,t*1.5\" 0.01 20 25\n", x_3, ex + 0.5, ey + 0.5, ez + 0.5, COR, COG, COB);

		}

		void print_round(long double ex, long double ey, long double ez) {


			long double Omega1 = 2 * pi / 30;
			long double theta = 0.0;

			for (int64 A = 0; A <= 30; A++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex tickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0.75686274 0.8235294 0.9411764 1 0 0 0 0 3 \"s1,dis=%.10lf,t\" 0.1 1 1\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5, theta);
				theta += Omega1;
			}
			//particleex tickpolarparameter minecraft:end_rod 0 8 0 1 1 1 1 0 0 0 0 2 "s1,dis=0,t" 0.1 1 1
		}

		void summon_firework(long double ex, long double ey, long double ez) {

			fprintf(OutPutFile, "execute if score @p Timer matches %lld run summon minecraft:firework_rocket %.10lf %.10lf %.10lf {FireworksItem:{tag:{Fireworks:{Flight:0,Explosions:[{Type:0,Colors:[I;50431],FadeColors:[I;10239]},{Type:0,Colors:[I;65526],FadeColors:[I;16777215]},{Type:0,Colors:[I;12976383],FadeColors:[I;16711880]}]}},id:\"minecraft:firework_rocket\",Count:1},Life:0,LifeTime:1}\n", (int64)ex, ex, ey, ez);



			///summon minecraft:firework_rocket %.10lf %.10lf %.10lf {FireworksItem:{tag:{Fireworks:{Flight:0,Explosions:[{Type:0,Colors:[I;50431],FadeColors:[I;10239]},{Type:0,Colors:[I;65526],FadeColors:[I;16777215]},{Type:0,Colors:[I;12976383],FadeColors:[I;16711880]}]}},id:"minecraft:firework_rocket",Count:1},Life:0,LifeTime:1}
		}

		void F_Rose(long double ex, long double ey, long double ez) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -0.7853981634 0 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 0.7853981634 \"s1,dis,cr,cg,cb=(0.7853981634-t),5*cos(2*(0.7853981634-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0.7853981634 1.5707963268 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 1.5707963268 2.3561944902 \"s1,dis,cr,cg,cb=(2.3561944902+1.5707963268-t),5*cos(2*(2.3561944902+1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 2.3561944902 3.1415926535 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -3.1415926535 -2.3561944902 \"s1,dis,cr,cg,cb=(-2.3561944902-3.1415926535-t),5*cos(2*(-2.3561944902-3.1415926535-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -2.3561944902 -1.5707963268 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -1.5707963268 -0.7853981634 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),5*cos(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -3.1415926535 3.1415926535 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),5*cos(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.01,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void F_D_Rose(long double ex, long double ey, long double ez) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -0.7853981634 0 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 0.7853981634 \"s1,dis,cr,cg,cb=(0.7853981634-t),5*cos(2*(0.7853981634-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0.7853981634 1.5707963268 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 1.5707963268 2.3561944902 \"s1,dis,cr,cg,cb=(2.3561944902+1.5707963268-t),5*cos(2*(2.3561944902+1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 2.3561944902 3.1415926535 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -3.1415926535 -2.3561944902 \"s1,dis,cr,cg,cb=(-2.3561944902-3.1415926535-t),5*cos(2*(-2.3561944902-3.1415926535-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -2.3561944902 -1.5707963268 \"s1,dis,cr,cg,cb=t,5*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -1.5707963268 -0.7853981634 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),5*cos(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -0.7853981634 0 \"s1,dis,cr,cg,cb=t,3*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 0.7853981634 \"s1,dis,cr,cg,cb=(0.7853981634-t),3*cos(2*(0.7853981634-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0.7853981634 1.5707963268 \"s1,dis,cr,cg,cb=t,3*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 1.5707963268 2.3561944902 \"s1,dis,cr,cg,cb=(2.3561944902+1.5707963268-t),3*cos(2*(2.3561944902+1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 2.3561944902 3.1415926535 \"s1,dis,cr,cg,cb=t,3*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -3.1415926535 -2.3561944902 \"s1,dis,cr,cg,cb=(-2.3561944902-3.1415926535-t),3*cos(2*(-2.3561944902-3.1415926535-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -2.3561944902 -1.5707963268 \"s1,dis,cr,cg,cb=t,3*cos(2*t),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -1.5707963268 -0.7853981634 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),3*cos(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -3.1415926535 3.1415926535 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),5*cos(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.02,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -3.1415926535 3.1415926535 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),3*cos(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.06,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void L_Rose(long double ex, long double ey, long double ez) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -0.7853981634 0 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 0.7853981634 \"s1,dis,cr,cg,cb=(0.7853981634-t+0.7853981634),4*sin(2*(0.7853981634-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0.7853981634 1.5707963268 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 1.5707963268 2.3561944902 \"s1,dis,cr,cg,cb=(2.3561944902+1.5707963268-t+0.7853981634),4*sin(2*(2.3561944902+1.5707963268-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 2.3561944902 3.1415926535 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -3.1415926535 -2.3561944902 \"s1,dis,cr,cg,cb=(-2.3561944902-3.1415926535-t+0.7853981634),4*sin(2*(-2.3561944902-3.1415926535-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -2.3561944902 -1.5707963268 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -1.5707963268 -0.7853981634 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t+0.7853981634),4*sin(2*(-0.7853981634-1.5707963268-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -3.1415926535 3.1415926535 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),4*sin(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void L_D_Rose(long double ex, long double ey, long double ez) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -0.7853981634 0 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 0.7853981634 \"s1,dis,cr,cg,cb=(0.7853981634-t+0.7853981634),4*sin(2*(0.7853981634-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0.7853981634 1.5707963268 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 1.5707963268 2.3561944902 \"s1,dis,cr,cg,cb=(2.3561944902+1.5707963268-t+0.7853981634),4*sin(2*(2.3561944902+1.5707963268-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 2.3561944902 3.1415926535 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -3.1415926535 -2.3561944902 \"s1,dis,cr,cg,cb=(-2.3561944902-3.1415926535-t+0.7853981634),4*sin(2*(-2.3561944902-3.1415926535-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -2.3561944902 -1.5707963268 \"s1,dis,cr,cg,cb=t+0.7853981634,4*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -1.5707963268 -0.7853981634 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t+0.7853981634),4*sin(2*(-0.7853981634-1.5707963268-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -0.7853981634 0 \"s1,dis,cr,cg,cb=t+0.7853981634,6*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 0.7853981634 \"s1,dis,cr,cg,cb=(0.7853981634-t+0.7853981634),6*sin(2*(0.7853981634-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0.7853981634 1.5707963268 \"s1,dis,cr,cg,cb=t+0.7853981634,6*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 1.5707963268 2.3561944902 \"s1,dis,cr,cg,cb=(2.3561944902+1.5707963268-t+0.7853981634),6*sin(2*(2.3561944902+1.5707963268-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 2.3561944902 3.1415926535 \"s1,dis,cr,cg,cb=t+0.7853981634,6*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -3.1415926535 -2.3561944902 \"s1,dis,cr,cg,cb=(-2.3561944902-3.1415926535-t+0.7853981634),6*sin(2*(-2.3561944902-3.1415926535-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -2.3561944902 -1.5707963268 \"s1,dis,cr,cg,cb=t+0.7853981634,6*sin(2*(t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 -1.5707963268 -0.7853981634 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t+0.7853981634),6*sin(2*(-0.7853981634-1.5707963268-t+0.7853981634)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 5 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -3.1415926535 3.1415926535 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),4*sin(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 -3.1415926535 3.1415926535 \"s1,dis,cr,cg,cb=(-0.7853981634-1.5707963268-t),6*sin(2*(-0.7853981634-1.5707963268-t)),sin(t*1.5)/4+0.5,sin(t*2)/4+0.5,sin(t*2.5)/4+0.5\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.08,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void L_F_Rose(long double ex, long double ey, long double ez) {
			L_Rose(ex, ey, ez);
			F_Rose(ex, ey, ez);
		}

		void L_F_D_Rose(long double ex, long double ey, long double ez) {
			L_D_Rose(ex, ey, ez);
			F_D_Rose(ex, ey, ez);
		}

		void Butter_Fly_Spread(long double ex, long double ey, long double ez, int64 Temround, long double Butterflysize) {
			for (int64 i = 0; i < Butterflysize * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0 0 0 %.10lf \"s1,dis,cr,cg,cb=t,(E^cos(t)-2*cos(4*t)+(sin(t/12))^5)*%.10lf,0.25+0.25*sin(t),0.3+0.3*cos(t),0.6+0.1*sin(t)\" 0.02 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, pi * 2 * (fp64)Temround, (fp64)i / 10.0);
			}
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 %.10lf \"s1,dis,cr,cg,cb=t,(E^cos(t)-2*cos(4*t)+(sin(t/12))^5)*%.10lf,0.25+0.25*sin(t),0.3+0.3*cos(t),0.6+0.1*sin(t)\" 0.02 25\n", (int64)ex + (int64)Butterflysize * 10, ex + 0.5, ey + 0.5, ez + 0.5, pi * 2 * (fp64)Temround, (fp64)Butterflysize * 10.0 / 10.0);
		}

		void Flower_3d(long double ex, long double ey, long double ez) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 6.28 \"s1,dis,s2,cr,cg,cb=t,5*sin(4*t),0.5,0.75+0.25*sin(t),0.2+0.1*cos(t),0.2+0.2*cos(t)\" 0.01 45 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 6.28 \"s1,dis,s2,cr,cg,cb=t,5*sin(4*t),0.5,0.75+0.25*sin(t),0.2+0.1*cos(t),0.2+0.2*cos(t)\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void Flower_2d(long double ex, long double ey, long double ez) {
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 6.28 \"s1,dis,s2,cr,cg,cb=t,5*sin(4*t),0,0.75+0.25*sin(t),0.2+0.1*cos(t),0.2+0.2*cos(t)\" 0.01 45 15\n", (int64)ex, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbapolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 6.28 \"s1,dis,s2,cr,cg,cb=t,5*sin(4*t),0,0.75+0.25*sin(t),0.2+0.1*cos(t),0.2+0.2*cos(t)\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 15, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void Flower_3d_And_2d(long double ex, long double ey, long double ez) {
			Flower_3d(ex, ey, ez);
			Flower_2d(ex, ey, ez);
		}

		void SeaStar(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0 1 0 0.1 0 0 6.28 \"s1, dis = t, (-2 - 2 ^ (sin(5 * t)))*%.10lf\" 0.01 10\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 20);
			}
		}

		void SeaStar_Flash(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0 1 0 0 0 0 6.28 \"s1, dis = t, (-2 - 2 ^ (sin(5 * t)))*%.10lf\" 0.01 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 20);
			}
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0 1 0 0 0 0 6.28 \"s1, dis = t, (-2 - 2 ^ (sin(5 * t)))\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 20, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void SmallFlower(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0 1 0 0.1 0 0 6.28 \"s1, dis = t, (5 - 2 ^ (sin(5 * t)))*%.10lf\" 0.01 10\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 20);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0 1 0 0.1 0 0 6.28 \"s1, dis = t, %.10lf\" 0.01 10\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 20);
			}
		}

		void SmallFlower_Flash(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t, (5 - 2 ^ (sin(5 * t)))*%.10lf\" 0.01 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 20);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t, %.10lf\" 0.01 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 20);
			}
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t, (5 - 2 ^ (sin(5 * t)))\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 19, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t, 1\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 20, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void D_SmallFlower_Flash(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t, (5 - 2 ^ (sin(5 * t)))*%.10lf\" 0.02 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 40);
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t+0.3141592, (5 - 2 ^ (sin(5 * t)))*%.10lf\" 0.04 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 120);
			}
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t, (5 - 2 ^ (sin(5 * t)))/2\" 0.02 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 19, ex + 0.5, ey + 0.5, ez + 0.5);
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0.75 1 1 0 0 0 0 6.28 \"s1, dis = t+0.3141592, (5 - 2 ^ (sin(5 * t)))/6\" 0.04 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 20, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void SeaStar_to_F(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0 1 0 0.1 0 0 6.28 \"s1, dis = t, (-%.10lf - 2 ^ (sin(5 * t)))*0.25\" 0.01 10\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i);
			}
		}

		void SeaStar_to_F_Flash(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0.2 1 0 0 0 0 6.28 \"s1, dis = t, (-%.10lf - 2 ^ (sin(5 * t)))*0.25\" 0.01 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i);
			}
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.8 0.2 1 0 0 0 0 6.28 \"s1, dis = t, (-20 - 2 ^ (sin(5 * t)))*0.25\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 20, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void Eight_Flower__Flash(long double ex, long double ey, long double ez) {
			long double R = 2.0;
			for (int64 i = 0; i <= (int64)R * 10; i++) {
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.6 0 1 0 0 0 0 6.28 \"s1,dis=t,(2*abs(sin(2*t))+abs(sin(4*t)))*%.10lf\" 0.01 1\n", (int64)ex + i, ex + 0.5, ey + 0.5, ez + 0.5, (fp64)i / 10);
			}
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex polarparameter minecraft:end_rod %.10lf %.10lf %.10lf 1 0.6 0 1 0 0 0 0 6.28 \"s1,dis=t,(2*abs(sin(2*t))+abs(sin(4*t)))*2\" 0.01 0 \"a=0.05;(vx,,vy,,vz)=(-sin(a),0,-cos(a),,0,1,0,,cos(a),0,-sin(a))*(x*2*sin(a),,0.04,,z*2*sin(a))\"\n", (int64)ex + 20, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void Eight_flower_small_flower_star(long double ex, long double ey, long double ez) {
			Eight_Flower__Flash(ex, ey, ez);
			SeaStar_to_F_Flash(ex, ey, ez);
			SmallFlower_Flash(ex, ey, ez);
			SeaStar_Flash(ex, ey, ez);
		}

		void brust(long double ex, long double ey, long double ez) {
			int64 x_3 = (int64)ex;
			fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex rgbatickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf 0 0 0 0 2 \"s1,dis,cr,cg,cb=0,0,random(),random(),random()\" 0.01 20 25 \"(vx, vy, vz) = ((random(), random(), random()) - 0.5) * t / 20\"\n", x_3, ex + 0.5, ey + 0.5, ez + 0.5);
		}

		void brust_2(long double ex, long double ey, long double ez) {
			int64 x_3 = (int64)ex;
			long double COR;//RGB颜色通道
			long double COG;//RGB颜色通道
			long double COB;//RGB颜色通道
			int64 t;
			double S1, S2;


			COR = 0.5;//随机RGB颜色通道
			COG = 0.7;//随机RGB颜色通道
			COB = 0;//随机RGB颜色通道
			for (int64 i = 1; i < 20; i++) {
				t = (rand() % 2) + 5;
				S1 = (double)(rand() % 6280) / 6280.0;
				S2 = ((double)(rand() % 2355) / 2355) + 0.3925;
				fprintf(OutPutFile, "execute if score @p Timer matches %lld run particleex tickpolarparameter minecraft:end_rod %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf 1 0 0 0 0 %lld \"s1,s2,dis=%.5f,%.5f-(t/10),t\" 0.1 2 25 \"(vx, vz) = ((random(), random()) - 0.5) * t / 30\" 0.5\n", x_3, ex + 0.5, ey + 0.5, ez + 0.5, COR, COG, COB, t, S1, S2);
			}
		}
		const double pi = 3.1415926535;
		FILE* OutPutFile = nullptr;
	};
}

#endif