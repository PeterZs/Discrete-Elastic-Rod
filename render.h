#pragma once
#ifndef _RENDER_H_
#define _RENDER_H_

#include "interactive.h"

vector<Vector3d> draw_pos_list;
vector<Vector3d> draw_force_list;
vector<BishopFrame> draw_bframe_list;
vector<MaterialFrame> draw_mframe_list;

void DrawElasticRod()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(move_x, move_y, move_z, 0, 0, 0, 0, 0, 1);

	glPushMatrix();
	glRotatef(angle_x*radius, 0, 0, 1);
	glRotatef(angle_y * radius, 0, 1, 0);

	//axis
	glColor3d(0, 0, 0);
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(5, 0, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 5, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 5);
	glEnd();

	//rod
	glLineWidth(3.0f);
	glColor3d(1, 0, 0);
	glBegin(GL_LINE_STRIP);
	for (unsigned int i = 0; i < draw_pos_list.size(); i++)
	{
		glVertex3d(draw_pos_list[i](0), draw_pos_list[i](1), draw_pos_list[i](2));
	}
	glEnd();

	//force
	glLineWidth(2.0f);
	glColor3d(0, 1, 0);
	glBegin(GL_LINES);
	for (unsigned int i = 0; i < draw_pos_list.size(); i++)
	{
		glVertex3d(draw_pos_list[i](0), draw_pos_list[i](1), draw_pos_list[i](2));
		draw_force_list[i].normalize();
		draw_force_list[i] /= 2;
		glVertex3d(draw_pos_list[i](0) + draw_force_list[i](0), draw_pos_list[i](1) + draw_force_list[i](1), draw_pos_list[i](2) + draw_force_list[i](2));
	}
	glEnd();

	//bframe
	glColor3d(0, 0, 1);
	glBegin(GL_LINES);
	for (unsigned int i = 0; i < draw_bframe_list.size(); i++)
	{
		Vector3d origin = (draw_pos_list[i] + draw_pos_list[i + 1]) / 2.0;
		draw_bframe_list[i].u.normalize();
		draw_bframe_list[i].v.normalize();
		draw_bframe_list[i].t.normalize();

		Vector3d u = origin + draw_bframe_list[i].u/2.0;
		Vector3d v = origin + draw_bframe_list[i].v/2.0;
		Vector3d t = origin + draw_bframe_list[i].t/ 2.0;

		glVertex3d(origin(0), origin(1), origin(2));
		glVertex3d(u(0), u(1), u(2));
		glVertex3d(origin(0), origin(1), origin(2));
		glVertex3d(v(0), v(1), v(2));
		glVertex3d(origin(0), origin(1), origin(2));
		glVertex3d(t(0), t(1), t(2));
	}
	glEnd();

	//mframe
	/*
	glColor3d(1, 1, 0);
	glBegin(GL_LINES);
	for (unsigned int i = 0; i < draw_bframe_list.size(); i++)
	{
		Vector3d origin = (draw_pos_list[i] + draw_pos_list[i + 1]) / 2.0;
		draw_mframe_list[i].m1.normalize();
		draw_mframe_list[i].m2.normalize();
		draw_mframe_list[i].t.normalize();

		Vector3d m1 = origin + draw_mframe_list[i].m1 / 2.0;
		Vector3d m2 = origin + draw_mframe_list[i].m2 / 2.0;
		Vector3d t = origin + draw_bframe_list[i].t / 2.0;

		glVertex3d(origin(0), origin(1), origin(2));
		glVertex3d(m1(0), m1(1), m1(2));
		glVertex3d(origin(0), origin(1), origin(2));
		glVertex3d(m2(0), m2(1), m2(2));
		glVertex3d(origin(0), origin(1), origin(2));
		glVertex3d(t(0), t(1), t(2));
	}
	glEnd();
	*/

	glPopMatrix();
	glutSwapBuffers();
}

void init()
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("DER");
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_FLAT);
}

void reshape(int width, int height)
{
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0f, (GLfloat)width / (GLfloat)height, 1.0f, 100.f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void DrawElasticRod(ElasticRod& elastic_rod)
{
	draw_pos_list.assign(elastic_rod.pos_list.begin(), elastic_rod.pos_list.end());
	draw_force_list.assign(elastic_rod.force_list.begin(), elastic_rod.force_list.end());
	draw_bframe_list.assign(elastic_rod.bframe_list.begin(), elastic_rod.bframe_list.end());
	draw_mframe_list.assign(elastic_rod.mframe_list.begin(), elastic_rod.mframe_list.end());

	cout <<endl<< "Force List:" << endl;
	for (unsigned int i = 0; i < draw_force_list.size(); i++)
	{
		cout << draw_force_list[i].norm() << endl;
	}
	cout << endl << "dedq list:" << endl;
	for (int i = 0; i < elastic_rod.theta_list.size(); i++)
	{
		cout << elastic_rod.ComputedEdqj(i) << endl;
	}


	init();
	glutDisplayFunc(DrawElasticRod);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(mykeyboard);
	glutMouseFunc(mymouse);
	glutMotionFunc(rotate);
	glutMainLoop();
}


#endif