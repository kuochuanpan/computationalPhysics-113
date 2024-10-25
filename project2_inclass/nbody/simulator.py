import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from particles import Particles
from numba import jit, njit, prange, set_num_threads

"""
The N-Body Simulator class is responsible for simulating the motion of N bodies



"""

class NBodySimulator:

    def __init__(self, particles: Particles):
        
        self.particles = particles
        self.setup()

        return

    def setup(self, G=1,
                    rsoft=0.01,
                    method="RK4",
                    io_freq=10,
                    io_header="nbody",
                    io_screen=True,
                    visualization=False):
        """
        Customize the simulation enviroments.

        :param G: the graivtational constant
        :param rsoft: float, a soften length
        :param meothd: string, the numerical scheme
                       support "Euler", "RK2", and "RK4"

        :param io_freq: int, the frequency to outupt data.
                        io_freq <=0 for no output. 
        :param io_header: the output header
        :param io_screen: print message on screen or not.
        :param visualization: on the fly visualization or not. 
        """ 
        self.G = G
        self.rsoft = rsoft
        self.method = method
        self.io_freq = io_freq
        self.io_header = io_header
        self.io_screen = io_screen
        self.visualization = visualization
        self.nparticles = self.particles.nparticles # optional
        self.time = self.particles.time             # optional

        if method.lower() == "euler":
            self._advance_particles = self._advance_particles_Euler
        elif method.lower() == "rk2":
            self._advance_particles = self._advance_particles_RK2
        elif method.lower() == "rk4":
            self._advance_particles = self._advance_particles_RK4
        else:
            raise ValueError(f"Method not supported. {method}")

        return

    def evolve(self, dt:float, tmax:float):
        """
        Start to evolve the system

        :param dt: float, the time step
        :param tmax: float, the total time to evolve
        
        """
        particles = self.particles
        time  = particles.time
        nstep = 0

        while time < tmax:

            # 1. print some info on the screen
            if self.io_screen:
                print(f"Step: {nstep} | Time: {time:.3f}")
            # 2. check io
            if nstep % self.io_freq == 0:

                # recalculate the accelerations
                acc = self._calculate_acceleration(particles.nparticles, 
                                                   particles.masses, 
                                                   particles.positions)
                self.particles.accelerations = acc
                
                str_nstep = str(nstep).zfill(5)
                particles.output(f"{self.io_header}_{str_nstep}.txt")

            # 3. advance the particles
            self._advance_particles(dt, particles)

            # 4. update the time  
            nstep += 1 
            time += dt
            particles.time = time









        print("Simulation is done!")
        return

    def _calculate_acceleration(self, nparticles, masses, positions):
        """
        Calculate the acceleration of the particles
        """
        rsoft = self.rsoft
        accelerations = np.zeros_like(positions)
        
        for i in range(nparticles):
            for j in range(nparticles):
                if j>i:
                    rij = positions[j] - positions[i]
                    r = np.sqrt(np.sum(rij**2))
                    force = self.G * masses[j,0] * masses[i,0] / (r**2 + rsoft**2)**1.5
                    accelerations[i,:] = accelerations[i,:] + (force * rij)/masses[i,0]
                    accelerations[j,:] = accelerations[j,:] - (force * rij)/masses[j,0]

        return accelerations
        
    def _advance_particles_Euler(self, dt, particles):

        nparticles = particles.nparticles
        masses = particles.masses
        pos = particles.positions
        vel = particles.velocities
        acc = self._calculate_acceleration(nparticles,masses, pos)

        # update the positions
        pos = pos + vel*dt
        vel = vel + acc*dt

        particles.positions = pos
        particles.velocities = vel

        # update the accelerations (only for data IO)
        #acc = self._calculate_acceleration(nparticles,masses, pos)
        #particles.accelerations = acc

        return particles

    def _advance_particles_RK2(self, dt, particles):

        # TODO





        return particles

    def _advance_particles_RK4(self, dt, particles):
        
        #TODO








        return particles



if __name__ == "__main__":
    
    particles = Particles(N=2)
    particles.positions = np.array([[0,0,0],[1,0,0]])
    simulator = NBodySimulator(particles)
    simulator.setup(G=1, 
                    rsoft=0.01, 
                    method="Euler", 
                    io_freq=10, 
                    io_header="nbody", 
                    io_screen=True, 
                    visualization=False)
    simulator.evolve(dt=0.01, tmax=1.0)
    print("Done!")
