#ifndef CALCULATION_H
#define CALCULATION_H

#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QPixmap>
#include <QThread>
#include <QMutex>

class Particle;


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
        void stop();
        void rotate();
        void triggerSnapShot();
        Particle* getParticles(int &N);
		
		inline bool getMove()
		{
			return Move;
		}

        inline void setPotRangeScale(const double newScale)
        {
            potRangeScale = PS * newScale;
        }
		
		QMutex mutex;
		
	signals:
		void PictureChanged(double *XP, double *YP, double *ZP, int N);
        void WriteSnapShot(Particle* P, int N);
		
	private:
        enum Positions{temporaryPos, lastPos, currentPos, particles};

        void getU(const Particle* const P1, const Particle* const P2, double &U, const double* const tx, const double* const ty,
                  const double* const tz, Positions pos, double* ax = NULL, double* ay = NULL, double* az = NULL) const;
		void geta(double *tx, double *ty, double *tz, double *ax, double *ay, double *az);
        double getE(const Particle* const P, const double X, const double Y, const double Z, const bool lastPos) const;
        void checkE(const Particle* const P, const double tx, const double ty, const double tz, double& bx, double& by, double& bz, double& curMinU) const;
        void walkDownhil(const double targetU, const Particle* const currentParticle, double& rx, double &ry, double& rz) const;
        void correctLocalE();
        void initializeParticle(Particle &cP, const int x, const int z, const double X, const double Y, const double Z,
                                const double XF, const double YF, const double ZF) const;
        void WriteSnapshot();

        static void updateDelta(double& tuUpdate, double& delta, const double newValue);

        const double PS;
        double Energy, *Pot, *dPdR, Rm, RM, MaxX, MaxY, MaxZ, ScF, U, T, E, h, Re;
        double *XP, *YP, *ZP, *RepF, *RepP, **MAR, Speed, YMid, potRangeScale;
		int N, XS, YS, ZS, GridSizeDiv, nx, ny, nz, **MG, *MD, MXS, MZS, PXS, PYS, PZS, NPot;
		Particle *P, ****G, **D;
        bool Run, rotated, *Fixed, Move, writeSnapShot;
};

#endif // CALCULATION_H
