#ifndef PARTICLE_H
#define PARTICLE_H


struct Particle
{
    enum boundParticles{NBound = 4};
    Particle *next, *prev, *bound[NBound];
    int xp, yp, zp, NB;
    double X, Y, Z, vX, vY, vZ, lX, lY, lZ, lvX, lvY, lvZ, aaX, aaY, aaZ, E, deltaE, T, deltaT, U, deltaU, corrX, corrY, corrZ;
};

#endif // PARTICLE_H
