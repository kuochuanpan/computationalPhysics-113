import numpy as np
import matplotlib.pyplot as plt


class Particles:
    """
    Particle class to store particle properties
    """
    def __init__(self, N):
        self.nparticles = N
        self.time = 0
        self._tags = np.arange(N)
        self._masses = np.ones((N,1))
        self._positions = np.zeros((N,3))
        self._velocities = np.zeros((N,3))
        self._accelerations = np.zeros((N,3))        
        return

    @property
    def tags(self):
        return self._tags
    
    @tags.setter
    def tags(self, new_tags):
        # check if the length of the new_tags same as the self.nparticles and
        # check if the size of the new_tags is a one dimensional array
        if len(new_tags) != self.nparticles or len(new_tags.shape) != 1:
            raise ValueError("Number of tags must match number of particles")
        self._tags = new_tags
        return
    
    @property
    def masses(self):
        return self._masses
    
    @masses.setter
    def masses(self, new_masses):
        # check if the shape of the new_masses is (nparticles, 1)
        if new_masses.shape != (self.nparticles, 1):
            raise ValueError("The shape of the masses must be (nparticles, 1)")
        
        # additional check
        # new features

        self._masses = new_masses
        return
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self, new_positions):
        # check if the shape of the new_positions is (nparticles, 3)
        # TODO: 

        self._positions = new_positions
        return
    
    @property
    def velocities(self):
        return self._velocities
    
    @velocities.setter
    def velocities(self, new_velocities):
        # check if the shape of the new_velocities is (nparticles, 3)
        # TODO:
        
        self._velocities = new_velocities
        return
    
    @property
    def accelerations(self):
        return self._accelerations
    
    @accelerations.setter
    def accelerations(self, new_accelerations):
        # check if the shape of the new_accelerations is (nparticles, 3)
        # TODO:

        self._accelerations = new_accelerations
        return
    
    def output(self, filename):
        # output the particle properties to a file
        np.savetxt(filename, 
                   np.hstack((self.tags.reshape(-1,1), 
                              self.masses, 
                              self.positions, 
                              self.velocities, 
                              self.accelerations)), 
                              header="tag mass x y z vx vy vz ax ay az")
        return

    def draw(self, dimension=2,save=False): 
        # draw the particles
        if dimension == 2:
            plt.figure(1, figsize=(6,6))
            plt.scatter(self.positions[:,0], self.positions[:,1], 
                        s=10*self.masses, c='k', alpha=0.5)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('equal')
            
            if save:
                plt.savefig("particles_00000.png") # TODO
            else:
                    plt.show()

            # clear the fig to avoid overlap
            plt.clf()
            plt.cla()
            return
        elif dimension == 3:
            # TODO: 3D plot
            pass
        else:
            raise ValueError("Dimension must be 2 or 3")

        return

if __name__=="__main__":

    pts = Particles(N=5)
    print(pts.nparticles)
    print(pts.masses)
    print(pts.positions)
    print(pts.velocities)
    print(pts.accelerations)
    print(pts.tags)
    pts.output("test.txt")
    pts.draw(dimension=2,save=True)
    print("All tests passed!")