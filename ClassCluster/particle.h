#ifndef PARTICLE_H
#define PARTICLE_H


struct Particle
{
    Particle *next, *prev;
    int xp, yp, zp;
    double X, Y, Z, vX, vY, vZ, lX, lY, lZ, lvX, lvY, lvZ, aaX, aaY, aaZ;
};

#endif // PARTICLE_H
