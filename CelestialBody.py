import numpy as np
import numpy.linalg as linAlg


class CelestialBody(object):
    '''
    Below is the constructor of the class, taking as parameters the name, mass, initial velocity,
    initial position and time-step, and thus creating a CelestialBody object, to be later used
    in the Solar System simulation. In addition, in the constructor, a NumPy array of 3 rows and
    2 columns is initialised to zeros, to be later used to store the three values of acceleration
    needed to implement the Beeman algorithm (i.e. a(t-dt), a(t) and a(t+dt) ).
    Also, a list called periodTimes is used to store the times when a CelestialBody object crosses
    the axis from the 4th quadrant back to the 1st one, indicating a full orbital period.
    '''

    def __init__(self, name, mass, velocity, position, dt):
        self.name = name
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.dt = dt
        self.accelerations = np.zeros((3, 2))
        self.periodTimes = [0]  # time t=0 is the first value

    '''
    Below is a method called gravLaw, which uses Newton's Law of Universal Gravitation ( Fg = G*M*m / r^2 ), 
    to calculate the force on a CelestialBody object, and thus by dividing over its mass, 
    calculate and return its acceleration. The method takes as arguments the mass and position of 
    the other body, and uses these values in the calculations. 
    '''

    def gravLaw(self, otherMass, otherPosition):
        relativePos = otherPosition - self.position
        Gconstant = (6.67408e-11)
        # multiplied by the direction vector
        numerator = Gconstant * self.mass * otherMass * relativePos
        # modulus of relative position but cubed, since the numerator is multiplied by the direction vector
        denominator = (linAlg.norm(relativePos)) ** 3
        force = (numerator / denominator)
        return (force / self.mass)

    '''
    Below are two distinct methods called updatePosition and updateVelocity, which use the three-step
    Beeman scheme equations to predict the position at the next time step by combining the current acceleration 
    with the acceleration from the previous time step. This new position can then be used to calculate the new 
    acceleration which, in turn, predicts the new velocity. 
    accelerations[0] corresponds to a(t-dt), accelerations[1] corresponds to a(t) and
    accelerations[2] corresponds to a(t+dt).
    '''

    def updatePos(self):
        self.position = self.position + self.velocity * self.dt + (1 / 6) * (
                    4 * self.accelerations[1] - self.accelerations[0]) * ((self.dt) ** 2)

    def updateVel(self):
        self.velocity = self.velocity + (self.dt / 6) * (
                    2 * self.accelerations[2] + 5 * self.accelerations[1] - self.accelerations[0])
