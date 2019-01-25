#pragma once
#ifndef INTERACTIVE_H
#define INTERACTIVE_H

#define PI 3.1415926535
double move_x = 8.0, move_y = 8.0, move_z = 8.0;//the position of the camera
int oldposx = 0, oldposy = 0;
double radius = 0.3;
double angle_x = 0, angle_y = 0;

void mymouse(int button, int state, int x, int y);
void rotate(int x, int y);
void mykeyboard(unsigned char key, int x, int y);

void mymouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		oldposx = x; oldposy = y;
	}
}

void rotate(int x, int y)
{
	int delta_x = x - oldposx;
	int delta_y = y - oldposy;
	angle_x += (double)delta_x;
	angle_y += (double)delta_y;
	if (angle_x >= 1080) angle_x -= 1080;
	if (angle_y >= 1080) angle_y -= 1080;
	oldposx = x;
	oldposy = y;
	glutPostRedisplay();
}

void mykeyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 's':move_x += 0.5; break;
	case 'w':move_x -= 0.5; break;
	case 'd':move_y += 0.5; break;
	case 'a':move_y -= 0.5; break;
	case 'q':move_z -= 0.5; break;
	case 'e':move_z += 0.5; break;
	}
	glutPostRedisplay();
}

#endif // !INTERACTIVE_H