require 'rbbt-util'
require 'rbbt/resource'

module TFCLass
  extend Resource
  self.subdir = 'share/databases/TFCLass'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  TFCLass.claim TFCLass.tf_tg, :proc do 
    url="http://tfclass.bioinf.med.uni-goettingen.de/suplementary/TFClass_ontologies.zip"
  end
end

iif TFCLass.tf_tg.produce.find if __FILE__ == $0

