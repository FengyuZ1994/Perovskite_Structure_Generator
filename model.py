
class Atom:
  def __init__(self, name,x,y,z,element):
    self.name=name
    self.x = x
    self.y = y
    self.z = z
    self.element=element

class Structure:
    def __init__(self,name,atomList,acell,bcell,ccell):
        self.name=name
        self.atomList=atomList
        self.acell=acell
        self.bcell=bcell
        self.ccell=ccell

class Ligand:
    def __init__(self,atomList):
        self.atomList=atomList
