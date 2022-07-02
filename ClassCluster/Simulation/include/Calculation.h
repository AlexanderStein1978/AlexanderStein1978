#ifndef CALCULATION_H
#define CALCULATION_H

#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QPixmap>
#include <QThread>
#include <QMutex>

class Particle;
class Potential;
class CalculationTest;
class WatchPoint;
class Vector;
class PotentialDefinerInputData;
class QFile;
class QTextStream;

struct PotStruct;


class Calculation : public QThread
{
	Q_OBJECT
	
	public:
        enum PotRole{ClosestTwo, NextTwo, SecondOrder, Remaining, NumPot};

        Calculation(PotStruct* PotSs = nullptr, QObject* parent = 0);
		~Calculation();
		void initialize();
		void run();
		void setStepSize(double h);
		double getStepSize();
        double getPotentialEnergy() const;
        double getKineticEnergy() const;
        double setKineticEnergy(const double newT);
		void getSize(int &width, int &height);
		void getScales(double &ScF, int &MaxZ);
		void move();
		void setSpeed(double S);
		double getSpeed();
        void stop();
        void rotate();
        void triggerSnapShot();
        Particle* getParticles(int &N);
        bool setPotential(const PotRole role, PotStruct &Pot);
        void setLayerDistance(double newDistance);
        bool arePotentialsOK();
        void CalcEndpointsOfEnergyDefinitionAxis(const int particeIndex, const Vector& direction, Vector& end1, Vector& end2) const;
        void GetAxisEnergies(PotentialDefinerInputData& data);
        
        void setParticleWatchPoint(WatchPoint* point)
        {
            ParticleWatchPoint = point;
        }
        
        inline void setParticleWatchIndex(const int indexToWatch)
        {
            watchParticle = indexToWatch;
        }
		
		inline bool getMove()
		{
			return Move;
		}

        inline void setPotRangeScale(const double newScale)
        {
            potRangeScale = PS / newScale;
        }

        inline double getRe() const
        {
            return Re;
        }

        inline void setWatchParticle(const int index)
        {
            watchParticle = index;
        }

        inline static int getNumACalcsPerIt()
        {
            return 4;
        }

        inline int getNumParticles() const
        {
            return N;
        }

        inline int getNumXDimParticles() const
        {
            return PXS;
        }
		
		QMutex mutex;
		
	signals:
        void PictureChanged(Vector* Pos, int N);
        void EnergiesChanged(double kineticEnergy, double totalEnergy);
        void WriteSnapShot(Particle* P, int N);
		
	private:

        const double Error_Double;

        friend class CalculationTest;

        struct MARStruct
        {
            double R;
            double index;
        };

        enum Positions{temporaryPos, lastPos, currentPos, particles};
        enum Result{Success, Error};

        void rk4(Vector *t0, Vector *dvt, Vector *a, Vector *dt, Vector *dm, Vector *dvm, const double lh);
        Result getU(const Particle* const P1, const Particle* const P2, double &U, const Vector* const t0, Positions pos,
                  Vector *a) const;
        Result geta(Vector *t0, Vector *a);
        bool wasStepOK() const;
        double getE(const Particle* const P, const Vector& R, const bool lastPos) const;
        void correctLocalE();
        void initializeParticle(Particle &cP, const int x, const int z, const Vector& R, const Vector& Fact) const;
        void WriteSnapshot();
        static void updateDelta(double& tuUpdate, double& delta, const double newValue);
        void calcMAR();
        bool updateBindings();
        static double dist(const Particle *const P1, const Particle *const P2);
        static bool isNotBound(const Particle *const P1, const Particle *const P2);
        void checkPotentials();
        void checkPotential(const PotRole role);
        void updateBlock(int particleIndex);

        bool potentialOK[NumPot];
        int N, XS, YS, ZS, GridSizeDiv, nx, ny, nz, **MG, *MD, MXS, MZS, PXS, PYS, PZS, NPot, watchParticle, particleWatchStep;
        const double PS;
        double Energy, **Pot, **dPdR, Rm, RM, MaxX, MaxY, MaxZ, ScF, U, T, E, h, Re;
        double Speed, YMid, potRangeScale;
        Vector *Pos;
		Particle *P, ****G, **D;
        WatchPoint* ParticleWatchPoint;
        MARStruct **MAR;
        bool Run, rotated, *Fixed, Move, writeSnapShot;
        QFile* DebugLogFile;
        QTextStream* DebugLog;
};

#endif // CALCULATION_H
