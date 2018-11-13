   def __init__(self, U, MassConversion, save=False):

      self.U = U
      self.MassConversion = MassConversion

      # galaxy or cluster catalog name
      #self.name = "cmasssouth"
      #self.nameLong = "CMASS South"
      
      # path to the input catalog
      #self.pathInCatalog = "./input/?.txt"
      
      # Output path
      self.pathOut = "./output/catalog/"+self.name
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      # catalog path
      self.pathOutCatalog = self.pathOut + "/catalog.txt"
      
      # Figures path
      self.pathFig = "./figures/catalog/"+self.name
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      
      if save:
         self.readInputCatalog()
         self.addHaloMass()
         self.writeCatalog()
      
      self.loadCatalog()




   def copyCatalog(self, cat):
      
      newCat = Catalog(cat.U, cat.MassConversion, )
      self.U = cat.U
      self.MassConversion = cat.MassConversion

      self.name = cat.name
      self.nameLong = cat.nameLong

#      self.pathOut = cat.pathOut
#      self.pathOutCatalog = cat.pathOutCatalog
#
#      self.pathFig = cat.pathFig
