#pragma once
#ifndef ELASTICROD_H
#define ELASTICROD_H

#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <cassert>
#include <dlib/matrix.h>
#include <dlib/optimization.h>
#include <gl/glut.h>

using namespace Eigen;
using namespace std;

typedef dlib::matrix<double, 0, 1> col_vector;

typedef struct _bishop_frame {
	Vector3d t;
	Vector3d u;
	Vector3d v;
}BishopFrame;

typedef struct _material_frame {
	Vector3d t;
	Vector3d m1;
	Vector3d m2;
}MaterialFrame;

//假设一条ElasticRod上有n个点
class ElasticRod
{
public:
	ElasticRod(std::vector<Vector3d> ppos_list,
			   std::vector<bool> pis_clamped,
			   Vector3d pu0,
		       Matrix<double,2,2> pbending_stiffness,
			   double pbeta);
	~ElasticRod();

public:
	std::vector<Vector3d>		pos_list;//控制点的位置序列							0~n-1  n个点
	std::vector<Vector3d>		edge_list;//centerline边序列,同时也是标架t序列		0~n-2  n-1条边
	std::vector<double>			elength_list;//边长序列								0~n-2  n-1条边  ei表示pi和pi+1
	
	std::vector<bool>			is_clamped;//是否固定								0~n-1  n个点
	std::vector<MaterialFrame> 	mframe_list;//material_frame序列					0~n-2  n-1个标架
	std::vector<BishopFrame>	bframe_list;//bishop_frame序列						0~n-2  n-1个标架
	std::vector<double>			voronoi_list;//Voronoi域序列，每个控制点都有其对应的voronoi域  每个点都有

	col_vector					theta_list;//旋转角theta序列                        0~n-2 n-1个theta
	std::vector<Vector3d>		kb_list;										  //1~n-2 n-2个kb (index 0~n-3)
	std::vector<Vector3d>       force_list;  

	Vector3d					u0;
	Matrix<double, 2, 2>		bending_stiffness;
	double						tot_energy;
	double						beta;

public:
	std::vector<Matrix<double, 3, 3>> minusGKb;
	std::vector<Matrix<double, 3, 3>> eqGKb;
	std::vector<Matrix<double, 3, 3>> plusGKb;

	std::vector<Vector3d> minusGH;
	std::vector<Vector3d> eqGH;
	std::vector<Vector3d> plusGH;

public:
	void ComputeEdgeList();			//计算center line的切向场以及边长序列
	void ComputeVoronoiList();		//Voronoi域以控制点为单位
	void ComputeBishopFrame();		 
	void ParallelTransport();
	void ComputeKbList();
	void SetupMaterialFrame();
	void ComputeEnergy();
	void ComputeGradientKb();
	void ComputeGradientHolonomy();
	void ComputeForce();
	void MinimizeTheta();

	Vector3d Rotate(Vector3d u, double angle, Vector3d axis);
	Vector2d ComputeW(Vector3d kb, MaterialFrame mframe);
	double ComputedEdqj(unsigned int j);
	Matrix<double, 2, 3> GetGradientW(unsigned int i,unsigned int j,unsigned int k);
	Vector3d GetGradientHolonomy(unsigned int i, unsigned int j);

};


struct Evaluate
{
	Evaluate(ElasticRod* _rod) { rod = _rod;}
	ElasticRod* rod;

	double operator() (const col_vector theta)const
	{
		int j = 0;
		for (int i = 0; i < theta.size(); i++)
		{
			if (rod->is_clamped[i]) continue;
			else  rod->theta_list(i) = theta(j++); //有问题，写反了!!!s
		}//bframe_list is only dependent on pos_list,so it will not change before x opt.
		//cout << "Iteration ";

		//for (int i = 0; i < rod->theta_list.size(); i++)
		//{

//			cout << rod->theta_list(i) << " ";
	//	}
		//cout << endl;

		
		rod->SetupMaterialFrame();
		rod->ComputeEnergy();
		return rod->tot_energy;
	}
};

struct EvaluateGradient
{
	EvaluateGradient(ElasticRod* _rod) { rod = _rod; }
	ElasticRod* rod;

	col_vector operator() (const col_vector theta)const
	{
		int j = 0;
		for (int i = 0; i < theta.size(); i++)
		{
			if (rod->is_clamped[i]) continue;
			else rod->theta_list(i) = theta(j++);
		}
		rod->SetupMaterialFrame();
	
		col_vector derivate;
		derivate.set_size(theta.size());
		for (int i = 0; i < theta.size(); i++)
		{
			derivate(i) = rod->ComputedEdqj(i);
		}
		
		cout << "Gradient ";

		for (int i = 0; i < derivate.size(); i++)
		{

			cout << derivate(i) << " ";
		}
		cout << endl;

		return derivate;
	}
};






#endif