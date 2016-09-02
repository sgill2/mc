from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import simtk.unit as unit
import numpy as np
def zero_masses( system, firstres, lastres):
    for index in range(firstres, lastres):
        system.setParticleMass(index, 0*daltons)

def beta(temperature):
    kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
    kT = kB * temperature
    beta = 1.0 / kT
    return beta



def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        forcegroups[force] = i
    return forcegroups

def getEnergyDecomposition(context, forcegroups):
    energies = {}
    for f, i in forcegroups.items():
        energies[f] = context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
    return energies

class SmartDarting(object):
    def __init__(self, temperature, ligList):
        self.dartboard = []
        self.dart_size = 0.2*unit.nanometers
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
        kT = kB * temperature
        beta = 1.0 / kT
        self.beta = beta
        self.firstres = ligList[0]
        self.lastres = ligList[-1]
        self.total_mass = 0
        self.mass_list = None
        self.acceptance = 0


    def get_particle_masses(self, system, firstres=None, lastres=None):
        if firstres == None:
            firstres = self.firstres
        if lastres == None:
            lastres = self.lastres
        mass_list = []
        total_mass = 0*unit.dalton
        for index in range(firstres, lastres):
            mass = system.getParticleMass(index)
            total_mass = total_mass + mass
            mass_list.append([mass])
        total_mass = np.sum(mass_list)
        mass_list = np.asarray(mass_list)
        mass_list.reshape((-1,1))
        total_mass = np.array(total_mass)
        total_mass = np.sum(mass_list)
        temp_list = np.zeros(((lastres-firstres), 1))
        for index in range(lastres-firstres):
            mass_list[index] = (np.sum(mass_list[index])).value_in_unit(unit.daltons)
        mass_list =  mass_list*unit.daltons
        self.total_mass = total_mass
        self.mass_list = mass_list
        return total_mass, mass_list


    def calculate_com(self, pos_state, total_mass=None, mass_list=None, firstres=None, lastres=None):
        """
        This controls the ability to run a ncmc simulation with MD
        Arguments
        ---------
        total_mass: simtk.unit.quantity.Quantity in units daltons, contains the total masses of the particles for COM calculation
        mass_list:  nx1 np.array in units daltons, contains the masses of the particles for COM calculation
        pos_state:  nx3 np. array in units.nanometers, returned from state.getPositions
        firstres:   int, first residue of ligand
        lastres:    int, last residue of ligand
    
        Returns
        -------
        rotation : nx3 np.array in units.nm
            positions of ligand after random rotation
        """
        if firstres == None:
            firstres = self.firstres
        if lastres == None:
            lastres = self.lastres
        if total_mass == None:
            total_mass = self.total_mass
        if mass_list == None:
            mass_list = self.mass_list
        #choose ligand indicies
        copy_orig = copy.deepcopy(pos_state)
        lig_coord = copy_orig[firstres:lastres]
        lig_coord = lig_coord.value_in_unit(unit.nanometers)*unit.nanometers
        copy_coord = copy.deepcopy(lig_coord)
        #mass corrected coordinates (to find COM)
        mass_corrected = mass_list / total_mass * copy_coord
        sum_coord = mass_corrected.sum(axis=0).value_in_unit(unit.nanometers)
        com_coord = [0.0, 0.0, 0.0]*unit.nanometers
        #units are funky, so do this step to get them to behave right
        for index in range(3):
            com_coord[index] = sum_coord[index]*unit.nanometers
        #remove COM from ligand coordinates to then perform rotation
        return com_coord    



    def add_dart(self, dart):
        self.dartboard.append(dart)

    def calc_from_center(self, com):
        distList = []
        diffList = []
        indexList = []
        for dart in self.dartboard:
#            diff = dart - com
            diff = com - dart

            print 'diff, dart, com', diff, dart, com
            dist = np.sqrt(np.sum((diff)*(diff)))*unit.nanometers
#            print 'dist', dist
            distList.append(dist)
            diffList.append(diff)
        selected = []
        for index, entry in enumerate(distList):
#            print distList
#            print type(distList)
            if entry <= self.dart_size:
                selected.append(entry)
                diff = diffList[index]
                indexList.append(index)
#            if entry._value > 2.5:
#                print('bugged')
#                exit()
        if len(selected) == 1:
            return selected[0], diffList[indexList[0]]
        elif len(selected) == 0:
            return None, diff
        elif len(selected) >= 2:
            print ('sphere size overlap, check darts')
            exit()

    def redart(self, changevec):
        dartindex = np.random.randint(len(self.dartboard))
        dvector = self.dartboard[dartindex]
#        chboard = dvector + diff   #EDIT!!!!!!! should be dvector + diff
        #chboard is the new dart location () moved by the changevector
        chboard = dvector + changevec   
#        chboard = dvector
        print 'chboard', chboard
        return chboard

    def dartmove(self, context, firstres=None, lastres=None):
        if firstres == None:
            firstres = self.firstres
        if lastres == None:
            lastres = self.lastres

        stateinfo = context.getState(True, True, False, True, True, False)
        oldDartPos = stateinfo.getPositions(asNumpy=True)
        oldDartPE = stateinfo.getPotentialEnergy()
        center = self.calculate_com(oldDartPos)
        selectedboard, changevec = dboard.calc_from_center(com=center)
        print('changevec', changevec)
        if selectedboard != None:
        #notes
        #comMove is where the com ends up after accounting from where it was from the original dart center
        #basically where it's final displacement location
            newDartPos = copy.deepcopy(oldDartPos)
            comMove = dboard.redart(changevec)
            print('comMove', comMove)
            print('center', center)
            vecMove = comMove - center
            print('vecMove', vecMove)
#            print('vecMove*2', vecMove/2)
#            print('tile', np.reshape(np.tile(vecMove, (lastres-firstres)),(-1,3) ) )
            tiled_vector = np.reshape(np.tile(vecMove, (lastres-firstres)), (-1,3))
#            print('tiledvec', tiled_vector)
#            print('selectedres', newDartPos[firstres:lastres])
#            print('add', newDartPos[firstres:lastres] + tiled_vector)
#            print newDartPos._value 
            add = newDartPos[firstres:lastres]._value + tiled_vector
#            newDartPos[firstres:lastres] = (newDartPos[firstres:lastres] + tiled_vector)
#            print('newDartpos', newDartPos)
            newDartPos = newDartPos._value
            newDartPos[firstres:lastres] = add
            print 'worked'

    #        print('dartmove', dartmove)
    #        print('changevec', changevec)
            #print dartmove
            print newDartPos
            #newDartPos[firstres:lastres] = dartmove
            context.setPositions(newDartPos)
            newDartInfo = context.getState(True, True, False, True, True, False)
            newDartPE = newDartInfo.getPotentialEnergy()
            logaccept = -1.0*(newDartPE - oldDartPE) * self.beta
            randnum = math.log(np.random.random())
            print('logaccept', logaccept, randnum)
            print('old/newPE', oldDartPE, newDartPE)
            if logaccept >= randnum:
                print('move accepted!')
                self.acceptance = self.acceptance+1
            else:
                print('rejected')
                context.setPositions(oldDartPos)










	
    

pdb = PDBFile('circle.pdb')
forcefield = ForceField('circle.xml')
system = forcefield.createSystem(pdb.topology,
         constraints=HBonds)
print 'removing', system.getForce(0)
system.removeForce(1)
system.removeForce(0)
numParticles = system.getNumParticles()
###custom nonbonded
pairwiseForce = CustomNonbondedForce("q/(r^2) + 4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2); q = q1*q2")
print('energy', pairwiseForce.getEnergyFunction())
rangeparticles = range(numParticles)
pairwiseForce.addInteractionGroup(rangeparticles[-3:], rangeparticles[:-3])
print 'Force num particles', pairwiseForce.getNumParticles()
pairwiseForce.addPerParticleParameter("sigma")
pairwiseForce.addPerParticleParameter("epsilon")
pairwiseForce.addPerParticleParameter("q")
for i in rangeparticles:
    pairwiseForce.addParticle()
    pairwiseForce.setParticleParameters(i,[0.324999852378,0.71128, 0])
#    print pairwiseForce.getParticleParameters(i)
#    print pairwiseForce.getPerParticleParameterName(i)
for i in rangeparticles[-3:]:
    pairwiseForce.setParticleParameters(i,[0.32,0.7, -1])
for i in rangeparticles[-2:]:
    pairwiseForce.setParticleParameters(i,[0.15,0.50, 0.75])
#pairwiseForce.setParticleParameters(140,[0.324999852378,0.71128, 1])
for i in [125, 128, 131, 134, 126, 129, 132, 135]:
    pairwiseForce.setParticleParameters(i,[0.324999852378,0.71128, 1])

pairwiseForce.setParticleParameters(132,[0.324999852378,0.71128, 1])
pairwiseForce.setParticleParameters(122,[0.324999852378,0.71128, 1])



system.addForce(pairwiseForce)
###

###harmonic
if 1:
    harmonic = HarmonicBondForce()
    harmonic.addBond(rangeparticles[-3], rangeparticles[-2], 0.30, 10000)
    harmonic.addBond(rangeparticles[-3], rangeparticles[-1], 0.30, 10000)
    system.addForce(harmonic)


##angle
if 1:
    angleForce = CustomAngleForce("0.5*k*(theta-theta0)^2")
    angleForce.addPerAngleParameter("k")
    angleForce.addPerAngleParameter("theta0")
    angleForce.addAngle(rangeparticles[-1], rangeparticles[-3], rangeparticles[-2], [1000, 1.5707963268])
    system.addForce(angleForce)




print('pairwise force')
integrator = LangevinIntegrator(100*kelvin, 1/picosecond, 0.002*picoseconds)
zero_masses(system, 0, numParticles-3)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(100*kelvin)
#simulation.minimizeEnergy()
numParticles = simulation.system.getNumParticles()
print(numParticles)
#zero_masses(simulation.system, 0, numParticles)
print simulation.system.getParticleMass(0)
stateinfo = simulation.context.getState(True, True, False, True, True, False)
#print stateinfo.getPositions()
simulation.reporters.append(DCDReporter('output.dcd', 10))
simulation.reporters.append(StateDataReporter('info.csv', 10, step=True,
        potentialEnergy=True, temperature=False))

dboard = SmartDarting(temperature=300*unit.kelvin, ligList=[rangeparticles[-3], rangeparticles[-1]+1])
dboard.get_particle_masses(system=simulation.system, firstres=None, lastres=None)            
dboard.dart_size = 0.20*unit.nanometers
dboard.add_dart((np.array([3.37, 2.64, 3.43]))*unit.nanometers)
dboard.add_dart((np.array([3.40, 3.34, 2.61]))*unit.nanometers)

print('running')
fgrps=forcegroupify(system)
print getEnergyDecomposition(simulation.context, fgrps) 
counter=0
for n in range(50):
    simulation.step(2000)
#    stateinfo = simulation.context.getState(True, True, True, True, True, True)
#    forces = stateinfo.getForces()
#print forces
#    print forces[-1]
#    stateinfo = simulation.context.getState(True, True, False, True, True, False)
#    oldDartPos = stateinfo.getPositions(asNumpy=True)
#    oldDartPE = stateinfo.getPotentialEnergy()
#   center = oldDartPos[-1]
#    print center
#    print type(center)
#    print dboard.dartboard[0]
#    print type(dboard.dartboard[0])
#    selectedboard, changevec = dboard.calc_from_center(com=center)
#    if selectedboard != None:
#        counter = counter+1
#        newDartPos = oldDartPos[:]
#        dartmove = dboard.redart(changevec)
#        print('dartmove', dartmove)
#        print('changevec', changevec)
#        newDartPos[-1] = dartmove
#        simulation.context.setPositions(newDartPos)
#        newDartInfo = simulation.context.getState(True, True, False, True, True, False)
#        newDartPE = newDartInfo.getPotentialEnergy()
#        print('old/newPE', oldDartPE, newDartPE)
    dboard.dartmove(simulation.context)
    simulation.step(500)
#        if counter == 800:
#            exit()
print dboard.acceptance
        
    
print getEnergyDecomposition(simulation.context, fgrps)

