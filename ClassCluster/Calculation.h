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
    double X, Y, Z, vX, vY, vZ, lX, lY, lZ, lvX, lvY, lvZ, aaX, aaY, aaZ, lU;
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
		void move();
		void setSpeed(double S);
		double getSpeed();
		
		inline bool getMove()
		{
			return Move;
		}
		
		QMutex mutex;
		
	public slots:
		void stop();
		void rotate();
		
	signals:
		void PictureChanged(double *XP, double *YP, double *ZP, int N);
		
	private:
		void geta(double *tx, double *ty, double *tz, double *ax, double *ay, double *az);
        double getE(const Particle* const P, const double X, const double Y, const double Z, const int mx, const int my, const int mz, const bool lastPos) const;
        void checkE(const Particle* const P, const double tx, const double ty, const double tz, double& bx, double& by, double& bz, double& curMinU,
                      const int mx, const int my, const int mz) const;
        void walkDownhil(const double targetU, const Particle* const currentParticle, double& rx, double &ry, double& rz, const int mx, const int my, const int mz) const;
        void correctLocalE();

		double Energy, *Pot, *dPdR, Rm, RM, MaxX, MaxY, MaxZ, ScF, PS, U, T, E, h, Re;
		double *XP, *YP, *ZP, *RepF, *RepP, **MAR, Speed, YMid;
		int N, XS, YS, ZS, GridSizeDiv, nx, ny, nz, **MG, *MD, MXS, MZS, PXS, PYS, PZS, NPot;
		Particle *P, ****G, **D;
		bool Run, rotated, *Fixed, Move;
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
		void move();
		void speedChanged();
		void draw(double *XP, double *YP, double *ZP, int N);
		
	protected:
		void closeEvent(QCloseEvent *event);
		
	private:
		QLineEdit *StepE, *EnE, *Speed;
		QPushButton *Start, *Restart, *End, *Rotate, *Move;
		Picture *Pict;
		Calculation Calc;
};

#endif // WINDOW_H
