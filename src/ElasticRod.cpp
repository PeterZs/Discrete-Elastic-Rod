#include "ElasticRod.h"

using namespace std;

ElasticRod::ElasticRod(vector<Vector3d>ppos_list,			//���캯��
					   vector<bool> pis_clamped,
					   Vector3d pu0,
					   Matrix<double,2,2> pbending_stiffness,
					   double pbeta)
{
	pos_list.assign(ppos_list.begin(), ppos_list.end());
	is_clamped.assign(pis_clamped.begin(), pis_clamped.end());
	u0 = pu0;
	bending_stiffness = pbending_stiffness;
	pbeta = beta;

	ComputeEdgeList(); //���ݶ������
	theta_list.set_size(edge_list.size());
	for (int i = 0; i < edge_list.size(); i++) theta_list(i) = 0;
	theta_list(0) = 0.3;
	ComputeVoronoiList(); //���ݶ�����l_i
	ComputeKbList(); //���ݶ�����kb_list
	ComputeBishopFrame(); //������Ť��bishop frame����
	mframe_list.resize(bframe_list.size());
	ComputeGradientKb();//����Ҫtheta�������Kb�ݶ�
	ComputeGradientHolonomy();//����Ҫtheta�������holonomy���ݶ�
	//���ϵ�����������ȡ����pos_list�ģ������pos opt֮ǰ��Ϊ��ֵ��
}

ElasticRod::~ElasticRod(){}

void ElasticRod::ComputeEdgeList()							//�ɸ����ĵ㹹��߳�,ei��ʾ���Ƶ�iΪ��� 0~n-2 n-1����
{
	for (unsigned int i=0;i<pos_list.size()-1;i++)
	{
		Vector3d temp;
		temp = pos_list[i + 1] - pos_list[i];
		edge_list.push_back(temp);
		double temp_length = temp.norm();
		elength_list.push_back(temp_length);
	}
	return;
}

void ElasticRod::ComputeVoronoiList()						//����voronoi�������l_i������D l_i=e_{i-1}+e_{i}
{
	double l0 = elength_list[0];
	voronoi_list.push_back(l0);
	for (unsigned int i = 1; i < pos_list.size()-1; i++)
	{
		double voronoi_length;
		voronoi_length = elength_list[i - 1] + elength_list[i];
		voronoi_list.push_back(voronoi_length);
	}
	int n = pos_list.size();
	double ln = elength_list[n - 2];
	voronoi_list.push_back(ln);
	return;
}

void ElasticRod::ComputeKbList()							//1~n-2 n-2��kb
{
	Vector3d dumb_head;
	kb_list.push_back(dumb_head);//Ϊ��ʹ�±걣��һ��,��Ϊ��һ����û��kb
	for (unsigned int i = 1; i < edge_list.size(); i++)
	{
		Vector3d kb = 2*edge_list[i - 1].cross(edge_list[i]);
		kb = kb / (edge_list[i - 1].dot(edge_list[i]) + elength_list[i - 1] * elength_list[i]);
		kb_list.push_back(kb);
	}
	return;
}

Vector3d ElasticRod::Rotate(Vector3d u,double angle, Vector3d axis)
{
	AngleAxisd rotation_vector(angle, axis);
	u = rotation_vector * u;
	return u;
}

void ElasticRod::ComputeBishopFrame() //���� bishop frame ����
{
	BishopFrame temp0;
	temp0.t = edge_list[0];//����
	temp0.u = u0;//��ֵ
	temp0.v = temp0.t.cross(temp0.u);
	bframe_list.push_back(temp0);

	for (unsigned int i = 1; i < edge_list.size(); i++)  //n-2  ÿ��������һ�����
	{
		Vector3d axis = edge_list[i - 1].cross(edge_list[i]);
		axis.normalize();
		double angle = acos(edge_list[i - 1].dot(edge_list[i]) / (elength_list[i - 1] * elength_list[i]));

		BishopFrame temp;
		temp.u = Rotate(bframe_list[i - 1].u, angle, axis);
		temp.t = edge_list[i];
		temp.v = temp.t.cross(temp.u);
		bframe_list.push_back(temp);
	}
	return;
}

void ElasticRod::ParallelTransport()  //Ĭ��bishop frame[0]��֪   ����ģ����ڸ���bishop frame���У������ò���
{
	for (unsigned int i = 1; i < edge_list.size(); i++)  //n-2
	{
		Vector3d axis = edge_list[i - 1].cross(edge_list[i]);
		axis.normalize();
		double angle = acos(edge_list[i - 1].dot(edge_list[i]) / (elength_list[i - 1] * elength_list[i]));
		
		bframe_list[i].u = Rotate(bframe_list[i - 1].u, angle, axis);
		bframe_list[i].t = edge_list[i];
		bframe_list[i].v = bframe_list[i].t.cross(bframe_list[i].u);
	}
	return;
}

void ElasticRod::SetupMaterialFrame()//��Ҫtheta_list
{
	for (unsigned int i = 0; i < bframe_list.size(); i++)
	{
		mframe_list[i].t = bframe_list[i].t;
		mframe_list[i].m1 = bframe_list[i].u*cos(theta_list(i)) + bframe_list[i].v*sin(theta_list(i));
		mframe_list[i].m2 = -bframe_list[i].u*sin(theta_list(i)) + bframe_list[i].v*cos(theta_list(i));
	}
	return;
}

Vector2d ElasticRod::ComputeW(Vector3d kb, MaterialFrame mframe)//����kb��material frame�е�����
{
	Vector2d W;
	W[0] = kb.dot(mframe.m2);
	W[1] = -kb.dot(mframe.m1);
	return W;
}

void ElasticRod::ComputeEnergy() //����������  ��������Ϊmaterial frame��theta list�Ѿ������ú�
{
	double energy = 0;
	Vector2d W;
	for (unsigned int i = 1; i < edge_list.size(); i++) //�ӵڶ����߿�ʼ����Ϊ�и�i-1
	{
		//bending energy
		//Vector2d Ĭ����2*1��������
		double term = 0;

		W = ComputeW(kb_list[i], mframe_list[i - 1]);
		term += (W.transpose()*bending_stiffness*W);

		W = ComputeW(kb_list[i], mframe_list[i]);
		term += (W.transpose()*bending_stiffness*W);

		energy += term / 2.0 / voronoi_list[i];

		//twisting energy
		double mi = (theta_list(i) - theta_list(i - 1));
		energy += beta * mi*mi / voronoi_list[i];
	}
	tot_energy = energy;
	return;
}


double ElasticRod::ComputedEdqj(unsigned int j)		//����������theta�ǵ�ƫ��
{
	double dEdqj = 0;
	Matrix<double, 2, 2> J;
	J(0, 0) = 0.0; J(0, 1) = -1.0; J(1, 0) = 1.0; J(1, 1) = 0.0;
	if (j > 0)
	{
		double term = 0;
		Vector2d W = ComputeW(kb_list[j], mframe_list[j]);
		term = (W.transpose()*J*bending_stiffness*W);
		term /= voronoi_list[j];
		term += 2 * beta*(theta_list(j) - theta_list(j - 1)) / voronoi_list[j];
		dEdqj += term;
	}
	if (j < theta_list.size() - 1)
	{
		double term = 0;
		Vector2d W = ComputeW(kb_list[j + 1], mframe_list[j]);
		term = W.transpose()*J*bending_stiffness*W;
		term /= voronoi_list[j + 1];
		term -= 2 * beta*(theta_list(j + 1) - theta_list(j)) / voronoi_list[j + 1];
		dEdqj += term;
	}
	return dEdqj;
}


void ElasticRod::ComputeGradientKb() //����kb�Ŀռ��ݶ�
{
	vector<Matrix<double, 3, 3>> edgeMatrix;//get the skew-symmetric matrix
	for (unsigned int i = 0; i < edge_list.size(); i++)//[e_i] ��˾���
	{
		Matrix<double, 3, 3> temp;
		temp(0, 0) = 0;                temp(0, 1) = -edge_list[i][2]; temp(0, 2) = edge_list[i][1];
		temp(1, 0) = edge_list[i][2];  temp(1, 1) = 0;                temp(1, 2) = -edge_list[i][0];
		temp(2, 0) = -edge_list[i][1]; temp(2, 1) = edge_list[i][0];  temp(2, 2) = 0;
		
		edgeMatrix.push_back(temp);
	}
	minusGKb.resize(edge_list.size());
	eqGKb.resize(edge_list.size());
	plusGKb.resize(edge_list.size());

	for (unsigned int i = 1; i < edge_list.size(); i++)
	{
		double denominator = (elength_list[i - 1] * elength_list[i]) + edge_list[i - 1].dot(edge_list[i]);
		assert(denominator != 0);
		minusGKb[i] = (2 * edgeMatrix[i] + kb_list[i] * edge_list[i].transpose()) / denominator;
		plusGKb[i] = (2 * edgeMatrix[i - 1] - kb_list[i] * edge_list[i - 1].transpose()) / denominator;
		eqGKb[i] = -(minusGKb[i] + plusGKb[i]);
	}
	return;
}

void ElasticRod::ComputeGradientHolonomy()//holonomy���ݶ��±��1��ʼ����Ϊ��kb�ҹ�
{
	minusGH.resize(edge_list.size());
	eqGH.resize(edge_list.size());
	plusGH.resize(edge_list.size());

	for (unsigned int i = 1; i < edge_list.size(); i++)
	{
		minusGH[i] = kb_list[i] / 2.0 / elength_list[i - 1];
		plusGH[i] = -kb_list[i] / 2.0 / elength_list[i];
		eqGH[i] = -(minusGH[i] + plusGH[i]);
	}
	return;
}


Vector3d ElasticRod::GetGradientHolonomy(unsigned int i, unsigned int j)//ע��˴�jΪ�ϱ꣬��Ӧ�����з���(10)
{
	Vector3d GH(0, 0, 0);
	if (j >= i - 1&&i-1>=0&&i-1<plusGH.size())
	{
		GH += plusGH[i - 1];
	}

	if (j >= i&&i<eqGH.size()&&i>=0)
	{
		GH += eqGH[i];
	}
	if (j >= i + 1&&i+1<minusGH.size()&&i+1>=0)
	{
		GH += minusGH[i + 1];
	}
	return GH;
}

Matrix<double, 2, 3> ElasticRod::GetGradientW(unsigned int i, unsigned int j, unsigned int k) // \nabla_iw_k^j
{
	assert(j == (k - 1) || j == k);
	
	Matrix<double, 2, 3> GW;
	GW=MatrixXd::Zero(2, 3);
	
	Matrix<double, 2, 2> J;
	J(0, 0) = 0.0; J(0, 1) = -1.0; J(1, 0) = 1.0; J(1, 1) = 0.0;

	if (k >= i - 1 && k <= i + 1)//k���±괦�����⡣
	{
		GW(0, 0) = mframe_list[j].m2[0];
		GW(0, 0) = mframe_list[j].m2[1];
		GW(0, 0) = mframe_list[j].m2[2];

		GW(1, 0) = -mframe_list[j].m1[0];
		GW(1, 1) = -mframe_list[j].m1[1];
		GW(1, 2) = -mframe_list[j].m1[2];

		if (k == i - 1) GW = GW * plusGKb[k];
		else if (k == i) GW = GW * eqGKb[k];
		else if (k == i + 1)GW = GW * minusGKb[k];
	}
	Vector2d W = ComputeW(kb_list[k], mframe_list[j]);
	Vector3d GH = GetGradientHolonomy(i, j);
	GW = GW - J * W*GH.transpose();
	return GW;
}

void ElasticRod::ComputeForce()  //According to the formula in $7.1 General case 
{
	Matrix<double, 2, 2> J;
	J(0, 0) = 0.0; J(0, 1) = -1.0; J(1, 0) = 1.0; J(1, 1) = 0.0;

	for (unsigned int i = 0; i < pos_list.size(); i++)
	{
		Vector3d force = Vector3d::Zero();

		if (is_clamped[i] == true)
		{
			force_list.push_back(force);
			continue;
		}

		for (unsigned int k = std::max((int)i - 1, 1); k < edge_list.size(); k++)//�����±�ֵ��˼��һ��
		{
			Vector2d W;
			Matrix<double, 2, 3> GW;
			Vector3d term;

			W = ComputeW(kb_list[k], mframe_list[k - 1]);
			GW = GetGradientW(i, k - 1, k);
			term = GW.transpose()*bending_stiffness*W;

			W = ComputeW(kb_list[k], mframe_list[k]);
			GW = GetGradientW(i, k, k);
			term += GW.transpose()*bending_stiffness*W;

			term = term / voronoi_list[k];

			force -= term;
		}
		bool clamped = false;
		for (unsigned int i = 0; i < is_clamped.size(); i++)
		{
			if (is_clamped[i])
			{
				clamped = true;
				break;
			}
		}

		if (clamped)//��Ҫ�ж�
		{
			unsigned int n = edge_list.size();
			Vector3d GH = GetGradientHolonomy(i, n);
			double dEdqn = ComputedEdqj(n-1);
			force += dEdqn * GH;
		}

		//���ɼ��϶Ե��������������

		force_list.push_back(force);
	}

}

void ElasticRod::MinimizeTheta()
{
	Evaluate _evaluate(this);
	EvaluateGradient _evaluate_gradient(this);

	int clamped_count = 0;
	for (unsigned int i = 0; i < theta_list.size(); i++)
	{
		if (is_clamped[i]||is_clamped[i+1]) clamped_count++;
	}

	col_vector theta(theta_list.size()-clamped_count);
	
	cout <<"opt theta size:" << theta.size() << endl;
	int j = 0;
	for (unsigned int i = 0; i < theta_list.size(); i++)
	{
		if (is_clamped[i]||is_clamped[i+1]) continue;
		theta(j++) = theta_list(i);
	}

	double min_energy = dlib::find_min(dlib::bfgs_search_strategy(),

									   dlib::objective_delta_stop_strategy(1e-16f, 10000),
									   _evaluate,
									   _evaluate_gradient,
									   theta,
									   0.0);

	std::cout << "Theta List: " <<endl<< theta_list << endl;
	std::cout << "Min Energy: "<<min_energy<<endl;
	
	j = 0;
	for (unsigned int i = 0; i < theta_list.size(); i++)
	{
		if (is_clamped[i]||is_clamped[i+1]) continue;
		theta_list(i) = theta(j++);
	}
	SetupMaterialFrame();
	ComputeForce();
}


