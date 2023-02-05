#ifndef CALCULATION_H
#define CALCULATION_H

#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QPixmap>
#include <QThread>
#include <QMutex>

#include "particle.h"


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
        enum PotRole{ClosestTwo, NextTwo, SecondOrder, Remaining, Angular, NumPot};

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
		void getScales(double& rScF, int& rMaxX, int&rMaxY, int& rMaxZ);
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
        int TranslateParticleIndex(int index) const;
        
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

        inline bool IsRotated() const
        {
            return rotated;
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
        Result getU(Particle* const P1, Particle* const P2, double &U, const Vector* const t0, Positions pos, Vector *a, const bool collectCandidates) const;
        Result geta(Vector *t0, Vector *a, const bool collectCandidates);
        void addCandidate(Particle* const currPart, Particle *const candidate, const double dist) const;
        void removeBinding(Particle *const part, const int index) const;
        bool bindToRadical(Particle *const CP, Particle *const CanP, const double lastDist, const bool force) const;
        bool wasStepOK() const;
        double getE(Particle * const P, const Vector& R, const bool lastPos, const bool collectCandidates) const;
        void correctLocalE();
        void initializeParticle(Particle &cP, const int x, const int z, const Vector& R, const Vector& Fact) const;
        void WriteSnapshot();
        static void updateDelta(double& tuUpdate, double& delta, const double newValue);
        void calcMAR();
        bool updateBindings();
        bool UpdateBindings();
        bool goRight(Particle *const CP, Particle *const CanP, const int CPCanIndex, const int CanPBIndex);
        bool goLeft(Particle *const CP, Particle *const CanP, const int CPCanIndex, const int CanPBIndex);
        static double dist(const Particle *const P1, const Particle *const P2);
        static bool isNotBound(const Particle *const P1, const Particle *const P2);
        void checkPotentials();
        void checkPotential(const PotRole role);
        void updateBlock(int particleIndex);
        int* createRandomParticleOrder();
        void verifyNoBindingDoubled() const;
        bool isBindingDoubled(const int particleIndex) const;
        void updateBindingPairs();
        static int getPartnerBindingIndex(const Particle* const P, const int index);
        void applyMove(const double lh);
        static Particle* getClosestBound(const Particle *const P1, const Vector& Pos);
        static int getBindingIndexAtBound(const Particle *const P1, const int index);
        void getGridAtPos(const Vector& Pos, int& x, int& y, int& z) const;
        void doParticleLayerSwitch(Particle *const cP, const Vector& dist);

        bool potentialOK[NumPot];
        int N, XS, YS, ZS, GridSizeDiv, nx, ny, nz, **MG, *MD, MXS, MZS, PXS, PYS, PZS, NPot, watchParticle, particleWatchStep;
        const double PS;
        double Energy, **Pot, **dPdR, Rm, RM, MaxX, MaxY, MaxZ, ScF, U, T, E, h, Re;
        double Speed, YMid, potRangeScale, mLayerDistance;
        Vector *Pos;
		Particle *P, ****G, **D, **FixedWallPos;
        WatchPoint* ParticleWatchPoint;
        MARStruct **MAR;
        bool Run, rotated, Move, writeSnapShot, mRotationChanged;
        QFile* DebugLogFile;
        QTextStream* DebugLog;
};

#endif // CALCULATION_H
