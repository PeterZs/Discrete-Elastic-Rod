#include "ElasticRod.h"
#include "render.h"

int main(int argc,char* argv[])
{
	std::vector<Vector3d> pos_list;
	std::vector<bool> is_clamped;
	
	Vector3d u0(1,0,0);
	Matrix<double, 2, 2> bending_stiffness;
	double beta = 1.0;

	bending_stiffness(0, 0) = 1.0; bending_stiffness(0, 1) = 0;
	bending_stiffness(1, 0) = 0; bending_stiffness(1, 1) = 1.5;

	is_clamped.push_back(true);
	is_clamped.push_back(false);
	is_clamped.push_back(false);
	is_clamped.push_back(false);
	is_clamped.push_back(false);
	is_clamped.push_back(false);
	is_clamped.push_back(false);

	pos_list.push_back(Vector3d(0, 0, 0));
	pos_list.push_back(Vector3d(0, 1, 0.5));
	pos_list.push_back(Vector3d(0, 2, 0));
	pos_list.push_back(Vector3d(0, 3, 0.5));
	pos_list.push_back(Vector3d(0, 4, 0));
	pos_list.push_back(Vector3d(0, 5, 0.5));
	pos_list.push_back(Vector3d(0, 6, 0));

	ElasticRod elastic_rod(pos_list,is_clamped,u0,bending_stiffness,beta);	
		
	elastic_rod.MinimizeTheta();
	
	glutInit(&argc, (char**)argv);
	DrawElasticRod(elastic_rod);
	return 0;
}