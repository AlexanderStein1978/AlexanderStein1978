#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QPixmap>
#include <QThread>
#include <QMutex>

struct Particle
{
	Particle *next, *prev;
	int xp, yp, zp;
	double X, Y, Z, vX, vY, vZ;
};

class Picture : public QWidget
{
	public:
		Picture(QWidget *parent = 0);
		~Picture();
		QPixmap *getPixmap();
		
	protected:
		void paintEvent(QPaintEvent *event);
		
	private:
		QPixmap *Map;
};

class Calculation : public QThread
{
	Q_OBJECT
	
	public:
		Calculation(QObject* parent = 0);
		~Calculation();
		void initialize();
		void run();
		void setStepSize(double h);
		double getStepSize();
		double getEnergy();
		double setEnergy(double E);
		void getSize(int &width, int &height);
		void getScales(double &ScF, int &MaxZ);

		QMutex mutex;
		
	public slots:
		void stop();
		
	signals:
		void PictureChanged(double *XP, double *YP, double *ZP, int N);
		
	private:
		void geta(double *tx, double *ty, double *tz, double *ax, double *ay, double *az);
		
		double Energy, *Pot, *dPdR, Rm, RM, MaxX, MaxY, MaxZ, ScF, PS, U, T, E, h, Re;
		double *XP, *YP, *ZP;
		int N, XS, YS, ZS, GridSizeDiv, nx, ny, nz;
		Particle *P, ****G, **D;
		bool Run;
};

class Window : public QWidget
{
	Q_OBJECT
	
	public:
		Window();
		~Window();
		
	public slots:
		void run();
		void restart();
		void draw(double *XP, double *YP, double *ZP, int N);
		
	protected:
		void closeEvent(QCloseEvent *event);
		
	private:
		QLineEdit *StepE, *EnE;
		QPushButton *Start, *Restart, *End;
		Picture *Pict;
		Calculation Calc;
};

#endif // WINDOW_H
