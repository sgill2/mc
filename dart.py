from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
def zero_masses( system, firstres, lastres):
    for index in range(firstres, lastres):
        system.setParticleMass(index, 0*daltons)

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
    def __init__(self):
        self.dartboard = []
        self.dart_size = 0.2*nanometers

    def add_dart(self, dart):
        self.dartboard.append(dart)

    def calc_from_center(self, com):
        distList = []
        diffList = []
        for dart in self.dartboard:
            diff = dart - com
            print 'diff', diff
            dist = np.sqrt(np.sum((diff)*(diff)))*nanometers
            print 'dist', dist
            distList.append(dist)
            diffList.append(diff)
        selected = []
        for index, entry in enumerate(distList):
            print distList
            print type(distList)
            if entry <= self.dart_size:
                selected.append(entry)
                diff = diffList[index]
        if len(selected) >= 2:
            print ('sphere size overlap, check darts')
            exit()
        elif len(selected) == 1:
            return selected[0], diff
        elif len(selected) == 0:
            return None, diff

    def redart(self, changevec):
        dartindex = np.random.randint(len(self.dartboard))
        dvector = self.dartboard[dartindex]
#        chboard = dvector + diff   #EDIT!!!!!!! should be dvector + diff
        chboard = dvector + changevec   

        return chboard

            

dboard = SmartDarting()            
dboard.dart_size = 0.20*nanometers
dboard.add_dart((np.array([6.25, 5.99, 6.5]))*nanometers)
dboard.add_dart((np.array([5.75, 5.99, 6.5]))*nanometers)
dboard.add_dart((np.array([6.18, 6.33, 5.98]))*nanometers)
dboard.add_dart((np.array([5.82, 6.33, 5.98]))*nanometers)
dboard.add_dart((np.array([6.18, 5.66, 5.99]))*nanometers)
dboard.add_dart((np.array([5.82, 5.66, 5.99]))*nanometers)
dboard.add_dart((np.array([6.20, 6.00, 5.54]))*nanometers)
dboard.add_dart((np.array([5.80, 6.00, 5.54]))*nanometers)



print 'darts', dboard.dartboard


	
    

pdb = PDBFile('circle.pdb')
forcefield = ForceField('circle.xml')
system = forcefield.createSystem(pdb.topology,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
numParticles = system.getNumParticles()

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.02*picoseconds)
zero_masses(system, 0, numParticles-1)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)
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
print('running')
fgrps=forcegroupify(system)
print getEnergyDecomposition(simulation.context, fgrps) 
counter=0
for n in range(2000):
    simulation.step(50)
#    stateinfo = simulation.context.getState(True, True, True, True, True, True)
#    forces = stateinfo.getForces()
#print forces
#    print forces[-1]
    stateinfo = simulation.context.getState(True, True, False, True, True, False)
    oldDartPos = stateinfo.getPositions(asNumpy=True)
    oldDartPE = stateinfo.getPotentialEnergy()
    center = oldDartPos[-1]
    print center
    print type(center)
    print dboard.dartboard[0]
    print type(dboard.dartboard[0])
    selectedboard, changevec = dboard.calc_from_center(com=center)
    if selectedboard != None:
        counter = counter+1
        newDartPos = oldDartPos[:]
        dartmove = dboard.redart(changevec)
        print('dartmove', dartmove)
        print('changevec', changevec)
        newDartPos[-1] = dartmove
        simulation.context.setPositions(newDartPos)
        newDartInfo = simulation.context.getState(True, True, False, True, True, False)
        newDartPE = newDartInfo.getPotentialEnergy()
        print('old/newPE', oldDartPE, newDartPE)
        simulation.step(50)
        if counter == 10:
            simulation.step(100)
            exit()

        
    
print getEnergyDecomposition(simulation.context, fgrps)

