from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

#notes 2500 particle seems confined

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

pdb = PDBFile('circle.pdb')
forcefield = ForceField('circle.xml')
system = forcefield.createSystem(pdb.topology,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
numParticles = system.getNumParticles()

integrator = LangevinIntegrator(5000*kelvin, 1/picosecond, 0.02*picoseconds)
zero_masses(system, 0, numParticles-1)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(2500*kelvin)
#simulation.minimizeEnergy()
numParticles = simulation.system.getNumParticles()
print(numParticles)
#zero_masses(simulation.system, 0, numParticles)
print simulation.system.getParticleMass(0)
stateinfo = simulation.context.getState(True, True, True, True, True, True)
#print stateinfo.getPositions()
simulation.reporters.append(DCDReporter('output.dcd', 10))
simulation.reporters.append(StateDataReporter('info.csv', 10, step=True,
        potentialEnergy=True, temperature=False))
print('running')
fgrps=forcegroupify(system)
print getEnergyDecomposition(simulation.context, fgrps) 
for n in range(50000):
    simulation.step(50)
#    stateinfo = simulation.context.getState(True, True, True, True, True, True)
#    forces = stateinfo.getForces()
#print forces
#    print forces[-1]
print getEnergyDecomposition(simulation.context, fgrps)

