class PdbMaker(object):
    def __init__(self, numAtoms):
        self.atom = 'ATOM'
        self.serial = None
        self.atomname = None
#        self.alternate = 'A'
        self.resname = None
        self.chain = None
        self.seqnum = None
        self.xyz = None
        self.occupancy = '1.00'
        self.tempFactor = '0.00'
        self.elementSymbol = None
        self.numAtoms = numAtoms
        self.placement = {'atom' : [['ATOM' for x in range(numAtoms)], 1, 4 ], 
            'serial' : [[str(x) for x in range(1,numAtoms+1)], 7, 11], 
            'atomname' : [['N' for x in range(numAtoms)], 13, 16],  
            'resname' : [['TMP' for x in range(numAtoms)], 18, 20],  
            'seqnum' : [[str(x) for x in range(0, numAtoms)], 23, 26], 
            'x' : [None, 31, 38], 'y' : [None, 39, 46], 'z' : [None, 47, 54], 
            'occupancy' : [['1.00' for x in range(numAtoms)], 55, 60], 
            'tempFactor' : [['0.00' for x in range(numAtoms)], 61, 66], 
            'elementSymbol' : [['N' for x in range(numAtoms)], 77, 78]     }
        self.spacing = {'atom' : [0], 'serial' : [0], 
            'atomname' : [0],  'resname' : [0],  
            'seqnum' : [0], 'x' : [0], 'y' : [0], 'z' : [0],'occupancy' : [0],
            'tempFactor' : [0], 'elementSymbol' : [0]     }



    def calc_spacing(self):
        for dictKey in self.placement.keys():
            print dictKey
            tempSpace = []

            for index in range(self.numAtoms):
#                print self.placement[dictKey][2]
#                print self.placement[dictKey][1]
#                print self.placement[dictKey][0][index]
#                print len(self.placement[dictKey][index][0])
                space =  (self.placement[dictKey][2] - self.placement[dictKey][1] + 1) - len(self.placement[dictKey][0][index])
                tempSpace.append(space)
            print tempSpace
            self.spacing[dictKey] = tempSpace[:]

    def makePDB(self, output_file='output.pdb'):
        with open(output_file, 'w') as f:
            for line in range(self.numAtoms):
                pdbstr = (self.placement['atom'][0][line] + '  ' + ' '*self.spacing['serial'][line] + self.placement['serial'][0][line] + 
                    ' ' + ' '*self.spacing['atomname'][line] + self.placement['atomname'][0][line] + ' ' + 
                    ' '*self.spacing['resname'][line] + self.placement['resname'][0][line] + ' ' + 'A' + 
                    ' '*self.spacing['seqnum'][line] + self.placement['seqnum'][0][line] + '    ' +
                    ' '*self.spacing['x'][line] + self.placement['x'][0][line] +
                    ' '*self.spacing['y'][line] + self.placement['y'][0][line] +
                    ' '*self.spacing['z'][line] + self.placement['z'][0][line] +
                    ' '*self.spacing['occupancy'][line] + self.placement['occupancy'][0][line] +
                    ' '*self.spacing['tempFactor'][line] + self.placement['tempFactor'][0][line] + 10*' ' + 
                    self.placement['elementSymbol'][0][line]) + '\n'
                f.write(pdbstr)






